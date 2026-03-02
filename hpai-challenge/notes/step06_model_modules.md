# Modularising the DAGs (Workflow Step 6)

This step decomposes the refined DAGs (step 05) into independent modules that can be developed, validated, and — where possible — fitted separately. When something goes wrong during fitting, modular structure makes it possible to isolate the problem to a specific component.

The question driving the decomposition is: what is separable? A module should have its own parameters and data inputs, and be something you can validate (at minimum through prior predictive checks) or fit on its own. Components that lack estimated parameters or that only make sense within a composed model are shared infrastructure.

---

## Modules

### Module A: Spillover

**Purpose**: external introduction of infection from wild birds, with temporal profile and zone-dependent rates.

**Parameters (estimated)**: $t_0$, $\phi_{\text{hrz}}$, $\phi_{\text{non}}$, $\delta$.

| Parameter | Symbol | Prior | Role |
|---|---|---|---|
| Spillover onset day | $t_0$ | $\text{Normal}(15, 5)$ truncated to $[1, 44]$ | Day spillover begins ($t = 1$ is 1 Dec 2025) |
| HRZ spillover rate | $\phi_{\text{hrz}}$ | $\text{LogNormal}(\log(10^{-3}), 1.0)$ | Daily per-farm spillover in HRZ at onset |
| Non-HRZ spillover rate | $\phi_{\text{non}}$ | $\text{LogNormal}(\log(10^{-4}), 1.0)$ | Daily per-farm spillover outside HRZ |
| Spillover decay rate | $\delta$ | $\text{Exponential}(\text{rate} = 50)$ (mean $= 0.02\ \text{day}^{-1}$) | Post-onset decline in spillover |

**Inputs**: farm HRZ membership (`hrz_32626.geojson`), time $t$.

**Output**: $\text{hazard}_{\text{spillover},j}(t)$ for each susceptible farm $j$.

**Equations**:

$$\psi(t) = \begin{cases} 0 & \text{if } t < t_0 \\ \exp(-\delta \cdot (t - t_0)) & \text{if } t \geq t_0 \end{cases}$$

$$\text{hazard}_{\text{spillover},j}(t) = \phi_j \times \psi(t)$$

where $\phi_j = \phi_{\text{hrz}}$ if $j \in \text{HRZ}$, else $\phi_{\text{non}}$.

**Why this is a module**: spillover has no dependence on the infection state of other farms. It is a complete standalone model — given parameter values and a set of farms with HRZ labels, it produces hazard values without needing any other component. Composed with the observation module, it forms a fittable model on its own.

**Standalone validation**:
- Unit tests: $\psi(t) = 0$ for $t < t_0$; correct HRZ assignment; exponential decay after onset.
- Prior predictive checks: simulated spillover-only epidemics should produce O(10–1000) cases, with early cases concentrated in HRZ.
- Synthetic data recovery: generate spillover-only data; confirm posteriors recover true values.
- Fit to real data: estimate $t_0$, $\phi_{\text{hrz}}$, $\phi_{\text{non}}$, $\delta$ (plus $\beta_{\text{duck}}$ from shared components). Expect to capture early dynamics and HRZ excess but fail on spatial clustering near previously infected farms.

### Module B: Transmission

**Purpose**: farm-to-farm spread, combining distance-dependent local transmission and movement-based transmission. These are grouped into a single module because they answer the same question (given infected farms, which susceptible farms get infected next?), share the same infectiousness profile $w(\tau)$, and are partially confounded — for nearby farm pairs that also have movement links, infections cannot be cleanly attributed to one pathway. The joint posterior of $(\beta, \alpha, p_{\text{mov}})$ must be examined together.

**Parameters (estimated)**: $\beta$, $\alpha$, $p_{\text{mov}}$.

| Parameter | Symbol | Prior | Role |
|---|---|---|---|
| Spatial transmission rate | $\beta$ | $\text{LogNormal}(\log(10^{-4}), 1.5)$ | Farm-to-farm transmission intensity |
| Spatial kernel scale | $\alpha$ | $\text{LogNormal}(\log(3500), 0.5)$ | Characteristic distance (metres) |
| Per-movement transmission probability | $p_{\text{mov}}$ | $\text{Beta}(2, 20)$ (mean $\approx 0.09$) | Probability that a movement from an infected farm infects the destination |

Step 05 fixed $p_{\text{mov}} = 0.01$ on identifiability grounds (confounded with the spatial kernel for nearby farm pairs). In a Bayesian framework, identifiability is a result rather than a precondition — the posterior either updates from the prior or it does not. Estimating $p_{\text{mov}}$ and checking whether the data are informative is preferable to fixing it a priori. Long-range movements (where the spatial kernel contributes negligibly) and the 6 pre-shipment-detected cases provide some identifying variation. If the posterior does not meaningfully update from the prior, that is itself a reportable finding.

**Fixed parameters**: $r = 1.0$/day (within-farm growth rate), $\tau_{\min} = 1$ day (hard latent period), $\sigma_{\text{test}} = 0.9$ (pre-shipment testing sensitivity).

**Inputs**: pairwise distances $d_{ij}$ (from `population.csv`), movement records from `movement.csv` (source, destination, date), infection times $T_i^I$ and removal times $T_i^R$, HRZ status of source farms, regulated zone status by farm and date.

**Output**: $\text{hazard}_{\text{transmission},j}(t)$ for each susceptible farm $j$.

#### Sub-component: local transmission

$$\text{hazard}_{\text{local},j}(t) = \beta \sum_{i \in \mathcal{I}(t)} w(t - T_i^I) \times K(d_{ij})$$

where $\mathcal{I}(t)$ is the set of farms infected and not yet removed at time $t$, and:

$$w(\tau) = \begin{cases} 0 & \text{if } \tau < \tau_{\min} \\ 1 - \exp(-r \cdot (\tau - \tau_{\min})) & \text{if } \tau \geq \tau_{\min} \end{cases}$$

Default kernel: $K(d) = \exp(-d / \alpha)$ (exponential, DAG A). Alternative: $K(d) = 1/(1 + (d/\alpha)^2)$ (Cauchy, DAG B). The module interface is the same for both — only $K(d)$ changes.

#### Sub-component: movement transmission

$$\text{hazard}_{\text{movement},j}(t) = \sum_{i \in \mathcal{I}(t)} M_{i \to j}(t) \times p_{\text{eff}}(i,t) \times w(t - T_i^I)$$

where $M_{i \to j}(t)$ is the number of recorded movements from farm $i$ to farm $j$ on day $t$, and:

$$p_{\text{eff}}(i,t) = \begin{cases} 0 & \text{if } i \text{ in regulated zone at } t \\ p_{\text{mov}} \times (1 - \sigma_{\text{test}}) & \text{if } i \in \text{HRZ} \\ p_{\text{mov}} & \text{otherwise} \end{cases}$$

Conditions evaluated in order. **Regulated zones** (3 km protection + 10 km surveillance around confirmed cases, lasting 28 days) are distinct from the **HRZ** (static wild bird spillover zone from `hrz_32626.geojson`). A farm can be in both; the regulated zone condition takes precedence.

#### Combined output

$$\text{hazard}_{\text{transmission},j}(t) = \text{hazard}_{\text{local},j}(t) + \text{hazard}_{\text{movement},j}(t)$$

**Why this is a module**: transmission needs infection times but nothing about spillover rates, HRZ membership (for the local sub-component), or species. It cannot generate infections from nothing — it needs seed infections — but given back-calculated infection times it is conditionally fittable: given the infection history up to day $t$, what is the probability of each new infection on day $t+1$?

**Conditional validation**:
- Unit tests: kernel decays with distance; infectiousness ramps up correctly; removed farms excluded; movements filtered by regulated zone status; pre-shipment testing modifier applied to HRZ sources.
- Synthetic data recovery: generate data from known $(\beta, \alpha, p_{\text{mov}})$ with spillover providing seed infections; confirm posteriors recover all three.
- Conditional fit: fix infection times at back-calculated values; estimate $(\beta, \alpha, p_{\text{mov}})$ from the spatial-temporal pattern.
- Identifiability diagnostics: examine the joint posterior of $(\beta, \alpha, p_{\text{mov}})$. Monitor $|\text{corr}(\beta, \alpha)|$ — if > 0.8, activate the $\beta_0$ reparameterisation (step 05, §4). Check whether $p_{\text{mov}}$ posterior updates from its prior; if not, report the degree of non-identifiability.
- Compare kernel variants: fit exponential (DAG A) and Cauchy (DAG B); compare via posterior predictive checks on spatial clustering.
- Check whether the 6 pre-shipment-detected cases are consistent with the posterior predictive distribution of movement-pathway infections.

### Module C: Observation

**Purpose**: map latent infection times to observed confirmation times, and define the likelihood that connects process modules to data.

**Fixed parameters**: compound delay $d = \mu_E + \mu_{ID} + (D \to C)$, combining latent/amplification period, detection delay, and suspicion-to-confirmation delay. Default: $d = 3.5 + 5.0 + 2.0 = 10.5$ days (rounded to nearest integer per case).

**Inputs**: infection times $T_j^I$ (from process modules or back-calculation), confirmation times $T_j^C$ (from `cases.csv`), at-risk periods (from removal component).

**Output**: log-likelihood $\mathcal{L}$.

**Likelihood structure**:

For each observed case $j$ with confirmation time $T_j^C$:
- Back-calculate infection day: $T_j^I = T_j^C - d$
- Case contribution: $\ell_j = \log\left(1 - \exp(-\lambda_j(T_j^I))\right)$

For each non-case farm $j$ (never confirmed, not preventively culled before study end):
- Survival contribution: $\ell_j = -\sum_{t=1}^{T_j^{\text{end}}} \lambda_j(t)$

Total: $\mathcal{L} = \sum_j \ell_j$.

**Inference approach**: deterministic back-calculation with fixed delays initially. Data augmentation MCMC (sampling latent $T_j^I$ jointly with parameters) is the target once the model structure is validated.

**Why this is a module**: the observation model can be validated independently of any process model. The compound delay assumptions can be checked against observed $T^S \to T^C$ intervals. Prior predictive checks on back-calculated infection time distributions require no process model.

**Standalone validation**:
- Check that back-calculated infection times are consistent with $T^S \to T^C$ data (median ~3 days from suspicion to confirmation; the remaining delay should be plausible for latent + recognition periods).
- Delay sensitivity analysis: grid over $\mu_E \in \{2, 3, 4, 5\}$ and $\mu_{ID} \in \{3, 4, 5, 6, 7\}$. Report how back-calculated infection times shift and whether any configurations produce implausible results (e.g. infection after confirmation).

---

## Shared components

These components lack estimated parameters of their own and serve as infrastructure that modules plug into.

### Species susceptibility

Multiplicative modifier on total hazard: $\beta_{\text{species}}[j] = 1$ for chicken farms, $\beta_{\text{duck}}$ for duck farms, applied uniformly across all pathways.

**Parameter (estimated)**: $\beta_{\text{duck}}$ with prior $\text{Beta}(2, 8)$. Mean 0.2 reflects the observed ~4× higher attack rate in chickens vs ducks; this conflates susceptibility, infectiousness, and detectability (step 02).

**Diagnostic**: if the posterior 95% CrI covers > 80% of the prior range, the parameter is not meaningfully informed by the data — fall back to scenario analysis with $\beta_{\text{duck}} \in \{0.5, 1.0, 1.5\}$. The 80% threshold is a pragmatic heuristic; any value in 70–90% serves the same purpose.

Species susceptibility is estimable in any composed model that includes at least one process module, but is not independently fittable — it requires a hazard model underneath.

### Removal

Tracks when farms leave the at-risk population. No estimated parameters; delays imputed from data.

**Reactive culling**: $T_j^R = T_j^C + \delta_{\text{reactive}}$, where $\delta_{\text{reactive}}$ is from `cases.csv` (median ~2 days).

**Preventive culling** (from 1 Jan): farms within 1 km of a confirmed case are culled at $T_j^R = T_{\text{trigger}} + \delta_{\text{prev}}$, with $\delta_{\text{prev}}$ estimated from 12 complete records; 40/52 missing dates imputed.

**Zone biosecurity**: farms inside an active surveillance zone (10 km around a confirmed case, lasting 28 days) have their hazard reduced by $(1 - \varepsilon)$ where $\varepsilon = 0.5$.

---

## Composition

The full force of infection is the sum of process module outputs, modified by species susceptibility and zone biosecurity:

$$\lambda_j(t) = \beta_{\text{species}}[j] \times \left[\text{hazard}_{\text{spillover},j}(t) + \text{hazard}_{\text{transmission},j}(t)\right]$$

If farm $j$ is inside an active surveillance zone: $\lambda_j(t) \mathrel{*}= (1 - \varepsilon)$.

The removal component determines each farm's at-risk period. The observation module provides the likelihood. Either process module can be composed independently with the observation module — omitting a module is equivalent to its hazard contribution being zero.

### Fitting strategy

Fitting proceeds by composing progressively richer models:

1. **A + C** (spillover + observation) — 5 estimated parameters ($t_0$, $\phi_{\text{hrz}}$, $\phi_{\text{non}}$, $\delta$, $\beta_{\text{duck}}$). Validates the observation model and likelihood machinery on a simple process. Expected to capture early dynamics but fail on spatial clustering.

2. **A + B + C** (+ transmission) — 8 estimated parameters (full set). Core model. Diagnostics: $\beta$–$\alpha$ correlation, $p_{\text{mov}}$ informativeness, spatial residuals, $\beta_{\text{duck}}$ informativeness.

At each stage, four validation checks before proceeding:
1. **Prior predictive checks** — simulated epidemics span a plausible range.
2. **Synthetic data recovery** — posteriors recover known parameter values.
3. **Fit to real data** — convergence diagnostics ($\hat{R}$, ESS, trace plots).
4. **Posterior predictive checks** — simulated data match observed summaries.

### Development order

1. **Data preprocessing** — shared across all compositions. Load and join data sources, build distance matrix, construct movement network, impute missing preventive cull dates, assign HRZ membership.
2. **Module C (observation)** — validate delay assumptions against $T^S \to T^C$ data; implement likelihood structure.
3. **Module A (spillover)** — implement and validate standalone; compose with C for first fit.
4. **Module B (transmission)** — implement local sub-component first (kernel computation is the expensive part), then movement sub-component; compose with A + C for full fit.

---

## Link to research questions

| Question | Composition required |
|----------|---------------------|
| Q1 (descriptive epidemiology) | None (pre-modelling) |
| Q2 (forecasting) | A + B + C (full) |
| Q3 (duck susceptibility) | A + B + C ($\beta_{\text{duck}}$ estimated) |
| Q4 (counterfactual: modified culling) | A + B + C (modify removal rules) |
| Q5 (counterfactual: modified delays) | A + B + C (modify $\delta_{\text{reactive}}$) |

---
