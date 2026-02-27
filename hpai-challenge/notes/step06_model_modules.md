# Modularising the DAGs (Workflow Step 6)

Following Abbott et al., this step decomposes the refined DAGs (step 05) into sub-models that can be independently developed, validated, and diagnosed before combining. The principle: start with simple sub-models, validate each in isolation, then add complexity incrementally. This makes it possible to diagnose poor fit, misspecification, or convergence issues before they become entangled.

---

## Decomposition rationale

The full model (step 05) has 7 estimated parameters, three transmission pathways, and a shared observation model. Fitting everything simultaneously makes it hard to tell whether problems originate in the process model, the observation model, or the data. Modularisation addresses this by building up the model in progressive fitting stages, each adding one mechanism.

The key insight from step 05 is that the three transmission pathways — spillover, local transmission, and movement — enter additively into the force of infection. This means we can "switch off" pathways by fixing their parameters to zero, producing nested sub-models that share the same observation model and likelihood structure.

Each module is a **complete, fittable model** (process + observation), not just a code component. This means every module can go through the full validation cycle: prior predictive checks, synthetic data recovery, fit to real data, posterior predictive checks.

---

## Module 1: Spillover only

### Purpose

Test whether spillover from wild birds alone can explain the observed case pattern. This is the simplest possible model: no farm-to-farm transmission, no movement. If spillover alone produces reasonable case counts and temporal dynamics, it establishes a baseline. If it fails (e.g. cannot explain spatial clustering or late-outbreak cases far from HRZ), that motivates adding local transmission.

### Parameters

**Estimated (5)**:

| Parameter | Symbol | Prior | Role |
|---|---|---|---|
| Spillover onset day | $t_0$ | $\text{Normal}(15, 5)$ truncated to $[1, 44]$ | Day spillover begins ($t = 1$ is 1 Dec 2025) |
| HRZ spillover rate | $\phi_{\text{hrz}}$ | $\text{LogNormal}(\log(10^{-3}), 1.0)$ | Daily per-farm spillover in HRZ at onset |
| Non-HRZ spillover rate | $\phi_{\text{non}}$ | $\text{LogNormal}(\log(10^{-4}), 1.0)$ | Daily per-farm spillover outside HRZ |
| Spillover decay rate | $\delta$ | $\text{Exponential}(\text{rate} = 50)$ (mean $= 0.02\ \text{day}^{-1}$) | Post-onset decline in spillover |
| Duck susceptibility | $\beta_{\text{duck}}$ | $\text{Beta}(2, 8)$ | Relative susceptibility (chicken = 1) |

**Fixed at zero**:

| Parameter | Value | Effect |
|---|---|---|
| $\beta$ | 0 | No local transmission |
| $p_{\text{mov}}$ | 0 | No movement transmission |

**Other fixed parameters**: as step 05 ($r = 1.0$/day, $\tau_{\min} = 1$ day, $\varepsilon = 0.5$, $\sigma_{\text{test}} = 0.9$).

### Process equations

For each susceptible farm $j$ at time $t$:

$$\psi(t) = \begin{cases} 0 & \text{if } t < t_0 \\ \exp(-\delta \cdot (t - t_0)) & \text{if } t \geq t_0 \end{cases}$$

$$\lambda_j(t) = \beta_{\text{species}}[j] \times \phi_j \times \psi(t)$$

where $\phi_j = \phi_{\text{hrz}}$ if $j \in \text{HRZ}$, else $\phi_{\text{non}}$; and $\beta_{\text{species}}[j] = \beta_{\text{duck}}$ for duck farms, 1 for chicken farms. Zone biosecurity reduction applies as in step 05.

### Observation model

Shared across all modules (described in full below in §Observation model).

### Data inputs

- `population.csv`: farm locations, species
- `activity.csv`: time-varying susceptibility
- `hrz_32626.geojson`: HRZ indicator
- `cases.csv`: confirmation times (fit target)
- `prev_culls.csv`: removal times (defines at-risk periods)

### Expected behaviour and validation criteria

- **Prior predictive checks**: simulated epidemics should produce O(10–1000) cases across plausible parameter draws, not 0 or all farms.
- **Synthetic data recovery**: generate data from known spillover parameters with $\beta = 0$; confirm posteriors recover true values.
- **Fit to real data**: expect $\phi_{\text{hrz}} > \phi_{\text{non}}$ (HRZ farms at higher risk). Expect the model to capture the early temporal pattern but struggle to explain spatial clustering of cases near previously infected farms.
- **Posterior predictive checks**: compare observed vs predicted (a) epidemic curve, (b) proportion of HRZ vs non-HRZ cases over time, (c) spatial distribution of cases. Spatial clustering beyond what spillover can explain motivates Module 2.
- **Diagnostic to trigger Module 2**: if posterior predictive p-value for nearest-neighbour distance between consecutive cases is < 0.05, or if there is systematic under-prediction of late-outbreak cases, local transmission is needed.

---

## Module 2: Spillover + local transmission

### Purpose

Add farm-to-farm spatial transmission. This is the core scientific model — most outbreak dynamics should be captured here. The question is whether the spatial kernel and transmission rate can explain the observed spatial-temporal clustering of cases.

### Parameters

**Estimated (7)** — the full parameter set from step 05:

| Parameter | Symbol | Prior | Role |
|---|---|---|---|
| Spillover onset day | $t_0$ | $\text{Normal}(15, 5)$ truncated to $[1, 44]$ | Day spillover begins |
| HRZ spillover rate | $\phi_{\text{hrz}}$ | $\text{LogNormal}(\log(10^{-3}), 1.0)$ | Daily per-farm spillover in HRZ at onset |
| Non-HRZ spillover rate | $\phi_{\text{non}}$ | $\text{LogNormal}(\log(10^{-4}), 1.0)$ | Daily per-farm spillover outside HRZ |
| Spillover decay rate | $\delta$ | $\text{Exponential}(\text{rate} = 50)$ (mean $= 0.02\ \text{day}^{-1}$) | Post-onset decline in spillover |
| Spatial transmission rate | $\beta$ | $\text{LogNormal}(\log(10^{-4}), 1.5)$ | Farm-to-farm transmission intensity |
| Spatial kernel scale | $\alpha$ | $\text{LogNormal}(\log(3500), 0.5)$ | Characteristic distance (metres) |
| Duck susceptibility | $\beta_{\text{duck}}$ | $\text{Beta}(2, 8)$ | Relative susceptibility (chicken = 1) |

**Fixed at zero**:

| Parameter | Value | Effect |
|---|---|---|
| $p_{\text{mov}}$ | 0 | No movement transmission |

**Other fixed parameters**: as step 05.

### Process equations

For each susceptible farm $j$ at time $t$:

$$\lambda_j(t) = \beta_{\text{species}}[j] \times \left[\phi_j \times \psi(t) + \beta \sum_{i \in \mathcal{I}(t)} w(t - T_i^I) \times \exp(-d_{ij} / \alpha)\right]$$

where $\mathcal{I}(t)$ is the set of farms infected and not yet removed at time $t$, and $w(\tau)$ is the within-farm infectiousness profile from step 02:

$$w(\tau) = \begin{cases} 0 & \text{if } \tau < \tau_{\min} \\ 1 - \exp(-r \cdot (\tau - \tau_{\min})) & \text{if } \tau \geq \tau_{\min} \end{cases}$$

Zone biosecurity and culling apply as in step 05.

### Data inputs

Same as Module 1, plus:
- Pairwise distance matrix between farms (precomputed from `population.csv`)

### Expected behaviour and validation criteria

- **Prior predictive checks**: simulated epidemics should show spatial clustering (cases near previous cases), not just spatially uniform spillover.
- **Synthetic data recovery**: generate data from known $(\beta, \alpha)$ values; confirm posteriors recover both. Monitor $|\text{corr}(\beta, \alpha)|$ — if > 0.8, the $\beta_0$ reparameterisation (step 05, §4) is needed.
- **Fit to real data**: expect $\beta > 0$ (local transmission present). Expect $\phi_{\text{hrz}}$ and $\phi_{\text{non}}$ to shift relative to Module 1 as some cases previously attributed to spillover are reassigned to local transmission.
- **Posterior predictive checks**: compare observed vs predicted (a) epidemic curve, (b) spatial distribution, (c) nearest-neighbour distances, (d) proportion of HRZ/non-HRZ cases. The spatial checks that failed in Module 1 should improve.
- **Key diagnostics**:
  - $\beta$–$\alpha$ posterior correlation: if $|\text{corr}| > 0.8$, switch to $\beta_0$ reparameterisation.
  - $\beta_{\text{duck}}$ posterior vs prior: if 95% CrI covers > 80% of prior range, fall back to scenario analysis.
  - Spatial residuals: any systematic patterns (e.g. specific regions with excess unexplained cases) may indicate missing mechanisms.
- **Diagnostic to trigger Module 3**: if cases connected by recorded movements are systematically under-predicted, or if long-range transmission events cannot be explained by the exponential kernel, movement transmission is needed.

---

## Module 3: Full model

### Purpose

Add movement transmission for mechanistic completeness. This is the full model from step 05. Movement transmission is included with fixed $p_{\text{mov}}$ rather than estimated, so the parameter count remains 7. The purpose of this module is to check whether explicitly accounting for movements changes the parameter estimates or predictions meaningfully, and to provide the complete model for counterfactual analysis (research questions Q4, Q5).

### Parameters

**Estimated (7)**: identical to Module 2.

**Fixed (no longer zero)**:

| Parameter | Value | Role |
|---|---|---|
| $p_{\text{mov}}$ | 0.01 | Per-movement transmission probability |
| $\sigma_{\text{test}}$ | 0.9 | Pre-shipment testing sensitivity (HRZ only) |

**Other fixed parameters**: as step 05.

### Process equations

The full force of infection from step 05:

$$\lambda_j(t) = \beta_{\text{species}}[j] \times \left[\phi_j \times \psi(t) + \beta \sum_{i \in \mathcal{I}(t)} w(t - T_i^I) \times \exp(-d_{ij} / \alpha) + \sum_{i \in \mathcal{I}(t)} M_{i \to j}(t) \times p_{\text{eff}}(i,t) \times w(t - T_i^I)\right]$$

where $M_{i \to j}(t)$ is the number of recorded movements from farm $i$ to farm $j$ on day $t$, and:

$$p_{\text{eff}}(i,t) = \begin{cases} 0 & \text{if } i \text{ in regulated zone at } t \\ p_{\text{mov}} \times (1 - \sigma_{\text{test}}) & \text{if } i \in \text{HRZ} \\ p_{\text{mov}} & \text{otherwise} \end{cases}$$

Zone biosecurity and culling apply as in step 05.

### Data inputs

Same as Module 2, plus:
- `movement.csv`: source, destination, date (broiler_1 → broiler_2 only)
- Regulated zone status by farm and date (derived from cases + zone rules)

### Expected behaviour and validation criteria

- **Synthetic data recovery**: generate data with known movement contribution; confirm that including movements does not distort spillover/local parameter estimates.
- **Fit to real data**: compare posteriors to Module 2. If $\beta$ and $\alpha$ shift substantially, the movement pathway was absorbing some transmission previously attributed to the spatial kernel. If posteriors are similar, movements contribute little to the overall fit (expected, given $p_{\text{mov}}$ is small and fixed).
- **Posterior predictive checks**: as Module 2, plus check whether the 6 pre-shipment-detected cases are consistent with the model's predicted movement hazard.
- **Sensitivity analysis**: vary $p_{\text{mov}}$ over $[0.001, 0.05]$; report how $\beta$ and $\alpha$ posteriors change.
- **Model for Q4/Q5**: this module provides the baseline for counterfactual simulations (modified culling, modified delays).

---

## Observation model (shared)

All three modules use the same observation model, as specified in steps 04 and 05.

### Compound delay

Infection-to-confirmation delay is fixed at a literature-informed value $d = \mu_E + \mu_{ID} + (D \to C)$, combining the latent/amplification period, detection delay, and suspicion-to-confirmation delay. Default: $d = 3.5 + 5.0 + 2.0 = 10.5$ days (rounded to nearest integer per case for back-calculation).

### Likelihood

For each observed case $j$ with confirmation time $T_j^C$:

- Back-calculate infection day: $T_j^I = T_j^C - d$
- Case contribution: $\ell_j = \log\left(1 - \exp(-\lambda_j(T_j^I))\right)$ (infection probability on back-calculated day)

For each non-case farm $j$ (never confirmed, not preventively culled before study end):

- Survival contribution: $\ell_j = -\sum_{t=1}^{T_j^{\text{end}}} \lambda_j(t)$ (survived the cumulative hazard)

Total log-likelihood: $\mathcal{L} = \sum_j \ell_j$.

### Inference approach

Deterministic back-calculation with fixed delays as the initial approach (step 05, §2). Data augmentation MCMC (sampling latent $T_j^I$ jointly with parameters) is the target approach once the base model structure is validated.

### Delay sensitivity analysis

Required for all modules. Grid over plausible delay values:
- $\mu_E \in \{2, 3, 4, 5\}$ days
- $\mu_{ID} \in \{3, 4, 5, 6, 7\}$ days

Report how posteriors for key parameters ($\beta$, $\phi_{\text{hrz}}$, $\phi_{\text{non}}$) vary across the grid.

---

## Module dependencies and development order

### Nesting structure

The modules are strictly nested:

$$\text{Module 1} \subset \text{Module 2} \subset \text{Module 3}$$

Module 1 is Module 2 with $\beta = 0$. Module 2 is Module 3 with $p_{\text{mov}} = 0$. This means:

- Shared code: the observation model, likelihood structure, spillover component, removal/culling logic, and species modifier are identical across all modules.
- Progression: moving from Module $n$ to Module $n+1$ adds one transmission pathway without modifying existing components.
- Comparison: because modules are nested, posterior predictive performance can be compared directly to assess whether each additional pathway improves fit.

### Development order

1. **Data preprocessing** — load and join data sources, build distance matrix, construct movement network, impute missing preventive cull dates, assign HRZ membership. Shared across all modules.

2. **Module 1 (spillover only)** — implement spillover hazard + observation model + likelihood. Validate the observation model and likelihood machinery on a simple process before adding transmission.

3. **Module 2 (spillover + local)** — add the spatial transmission kernel. This is where most development effort goes, as the kernel computation is the most expensive part of the likelihood.

4. **Module 3 (full model)** — add movement hazard lookup. Relatively straightforward given the Module 2 infrastructure.

### Validation at each stage

Each module goes through four validation steps before proceeding:

1. **Prior predictive checks** — simulate from the prior; verify that simulated epidemics span a plausible range of outcomes (not all trivial or all catastrophic).
2. **Synthetic data recovery** — generate data from known parameters; fit the model; verify posteriors concentrate around true values with appropriate coverage.
3. **Fit to real data** — run inference on observed case data; check convergence diagnostics ($\hat{R}$, ESS, trace plots).
4. **Posterior predictive checks** — simulate from the posterior; compare simulated data to observed data on key summary statistics.

Only proceed to the next module once the current module passes all four checks.

---

## Candidate DAG variants

Step 05 identified two candidate DAGs: DAG A (exponential kernel, base) and DAG B (Cauchy kernel, fat-tailed). The modular structure applies identically to both — the kernel function is the only difference, and it affects only Modules 2 and 3.

Development proceeds with DAG A (exponential) as the default. DAG B is fitted as a variant of Module 2 once DAG A's Module 2 is validated. Comparison via posterior predictive checks (particularly spatial clustering statistics) determines whether the heavier tail improves fit.

---

## Link to research questions

| Question | Module required |
|----------|----------------|
| Q1 (descriptive epidemiology) | None (pre-modelling) |
| Q2 (forecasting) | Module 3 (full model) |
| Q3 (duck susceptibility) | Module 2 or 3 ($\beta_{\text{duck}}$ estimated in both) |
| Q4 (counterfactual: modified culling) | Module 3 (modify intervention rules) |
| Q5 (counterfactual: modified delays) | Module 3 (modify $\delta_{\text{reactive}}$) |

---
