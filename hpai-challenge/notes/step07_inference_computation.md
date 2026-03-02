# Inference and Computation (Workflow Step 7)

This step specifies how to estimate the parameters defined in steps 05–06: the likelihood formulation, the treatment of latent infection times, the sampling algorithm, and the software. The model has 8 continuous estimated parameters (step 06, §Fitting strategy) and an analytical force-of-infection survival likelihood. The goal is to choose an inference approach that is efficient for this structure and identifies a clear pathway for increasing complexity if needed.

---

## Likelihood formulation

The likelihood follows from Module C (step 06). Each farm contributes one term depending on whether it was a confirmed case or survived uninfected.

**Case contribution** — for each confirmed farm $j$ with back-calculated infection day $T_j^I = T_j^C - d$ (where $d = 10.5$ days, step 06):

$$\ell_j = -\sum_{t=1}^{T_j^I-1}\lambda_j(t)\;+\;\log\left(1 - \exp\left(-\lambda_j(T_j^I)\right)\right)$$

The first term is survival from day 1 to day $T_j^I - 1$ (the farm was not infected on any earlier day); the second is infection on day $T_j^I$.

**Non-case survival** — for each farm $j$ that was never confirmed and not preventively culled before the end of the observation period:

$$\ell_j = -\sum_{t=1}^{T_j^{\text{end}}} \lambda_j(t)$$

where $T_j^{\text{end}}$ is the last day farm $j$ is in the susceptible pool (end of study, or date of preventive culling).

**Total log-likelihood**:

$$\mathcal{L} = \sum_j \ell_j$$

### Properties relevant to inference

| Property | Status | Implication |
|---|---|---|
| Analytical | Yes — no simulation or numerical integration | Fast evaluation; exact gradients |
| Differentiable | Yes — w.r.t. all 8 continuous parameters | Compatible with gradient-based samplers (NUTS) |
| Sparse | Only farms near infectious cases get non-zero spatial hazard | Precomputation of neighbour lists reduces cost |
| Modular | Spillover and transmission hazards enter additively | Modules can be included or excluded by zeroing their contribution |

The likelihood evaluation cost scales as $O(N_{\text{cases}} \times N_{\text{infectious neighbours}})$ for the spatial component, plus $O(N_{\text{farms}})$ for the spillover survival terms. The spillover survival can be aggregated into bulk categories (step 06, §Shared components: 2 HRZ classes × 2 species × active/inactive ≈ 8 bins), reducing the non-case sum from ~9000 farms to ~8 category counts.

---

## Treatment of latent infection times

The compound delay $d = \mu_E + \mu_{ID} + (D \to C)$ determines the back-calculated infection day for each case. Because $d$ enters through integer rounding ($T_j^I = T_j^C - \text{round}(d)$), it has zero gradient with respect to continuous parameters — it cannot be estimated by gradient-based samplers.

This is the most consequential methodological choice in the inference design. Three approaches are available.

### Option A: deterministic back-calculation (recommended first)

Fix $d = 10.5$ days (step 05, §Observation DAG Refinements) and back-calculate all infection times deterministically: $T_j^I = T_j^C - 11$ (after rounding).

**Advantages**:
- No additional latent variables — the model has only 8 continuous parameters.
- Fully compatible with NUTS; no mixed discrete/continuous sampling.
- Fast: single likelihood evaluation per MCMC step.

**Limitation**: ignores uncertainty in the delay. If key conclusions are sensitive to the assumed delay, this approach is insufficient.

**Mitigation**: grid sensitivity analysis over $\mu_E \in \{2, 3, 4, 5\}$ and $\mu_{ID} \in \{3, 4, 5, 6, 7\}$ (20 combinations, step 05 §Delay sensitivity analysis). Each grid point is a separate MCMC run with fixed delays. If posteriors for the transmission and spillover parameters shift materially across the grid, delay uncertainty must be propagated formally.

### Option B: data augmentation MCMC

Sample 103 discrete latent infection days $T_j^I$ jointly with the 8 continuous parameters. Each $T_j^I$ is drawn from a truncated set of plausible days (e.g. $T_j^C - [7, 14]$).

**Advantages**:
- Properly propagates delay uncertainty into parameter estimates.
- The "correct" Bayesian treatment.

**Limitations**:
- Adds 103 discrete latent variables, each requiring a Metropolis–Hastings update (no gradient information).
- Requires a Gibbs sampler composing NUTS (for continuous parameters) with MH (for discrete infection times): `Gibbs(NUTS(:continuous_params), MH(:infection_times))` in Turing.jl.
- Potential mixing problems: infection times are correlated with transmission parameters (later infection → higher required hazard), and discrete variables can trap the sampler in local modes.

### Option C: analytical marginalisation

For each case, sum the likelihood over plausible infection days:

$$\ell_j = \log \sum_{d_j \in \mathcal{D}} \pi(d_j) \times \left(1 - \exp(-\lambda_j(T_j^C - d_j))\right) \times \prod_{t < T_j^C - d_j} \exp(-\lambda_j(t))$$

where $\mathcal{D}$ is the set of plausible delays (e.g. 7–14 days) and $\pi(d_j)$ is a prior on the delay.

**Advantages**:
- Eliminates discrete latent variables entirely — the model remains continuous-only and NUTS-compatible.
- Properly accounts for delay uncertainty (within the specified support).

**Limitations**:
- Increases likelihood evaluation cost by a factor of $|\mathcal{D}|$ per case (~8× for $d \in [7, 14]$).
- The `logsumexp` operation slightly complicates gradient computation but remains differentiable.

### Recommended pathway

Start with **Option A** (deterministic back-calculation). This is sufficient for initial model development and validation. The grid sensitivity analysis will reveal whether conclusions depend on delay assumptions. If they do, transition to **Option C** (analytical marginalisation) — it retains NUTS compatibility while propagating delay uncertainty. Reserve **Option B** (data augmentation) for cases where the delay prior itself needs to be estimated from the data, which requires additional information not currently available.

| Approach | Parameters | Latent variables | Sampler | When to use |
|---|---|---|---|---|
| A: Fixed delays | 8 continuous | 0 | NUTS | Initial development; grid sensitivity |
| C: Marginalise | 8 continuous | 0 | NUTS | If sensitivity analysis reveals delay dependence |
| B: Data augmentation | 8 continuous | 103 discrete | Gibbs(NUTS, MH) | If delay prior itself is estimable |

---

## MCMC algorithm

### Sampler: NUTS

The No-U-Turn Sampler is the default choice for low-dimensional continuous posterior distributions with analytical gradients. With 8 parameters, the model sits comfortably within the regime where NUTS is most efficient.

**Automatic differentiation**: ForwardDiff.jl (forward-mode). Forward-mode AD has cost proportional to the number of input parameters; for $\leq 10$ parameters it is faster than reverse-mode. Set chunk size equal to the number of estimated parameters (5 for the A+C composition, 8 for A+B+C) so the full gradient is computed in a single forward pass.

### Sampling configuration

| Setting | Value | Rationale |
|---|---|---|
| Chains | 4 | Standard for convergence diagnostics; parallelisable |
| Warmup | 500 iterations | Sufficient for NUTS adaptation given low dimension |
| Sampling | 1000 iterations | 4000 total post-warmup draws across chains |
| Max tree depth | 10 | Default; monitor for saturation |
| Target acceptance | 0.8 | Turing.jl default for NUTS |

### Initialisation

Use a coarse grid search over the prior to find an approximate MAP estimate, then initialise all chains near this point with small random perturbations. This reduces warmup time and avoids chains starting in low-probability regions where the sampler wastes iterations adapting.

For the A+C composition (5 parameters), a grid of $\sim 5^5 = 3125$ evaluations is feasible (likelihood evaluation $< 2$ ms each, total $\sim 6$ seconds). For the full A+B+C model (8 parameters), a structured grid is too large; instead use $\sim 1000$ random draws from the prior and pick the best.

### Convergence criteria

| Diagnostic | Threshold | Action if failed |
|---|---|---|
| $\hat{R}$ (split) | $< 1.01$ | Increase sampling; check parameterisation |
| Bulk ESS | $> 400$ per parameter | Increase sampling; thin less aggressively |
| Tail ESS | $> 400$ per parameter | Increase sampling; check for heavy-tailed posteriors |
| Divergent transitions | 0 | Increase `adapt_delta`; check model specification |
| Tree depth saturation | $< 5\%$ at max depth | Increase max tree depth |

All diagnostics computed via MCMCDiagnosticTools.jl (part of the Turing ecosystem).

### Parameter-specific diagnostics

These diagnostics are motivated by the identifiability analysis in step 05 (§3–4) and should be checked after every fit.

| Diagnostic | Trigger | Response |
|---|---|---|
| $\lvert\text{corr}(\beta, \alpha)\rvert > 0.8$ | $\beta$–$\alpha$ ridge | Activate $\beta_0$ reparameterisation (step 05, §4) |
| $\beta_{\text{duck}}$ posterior 95% CrI covers $> 80\%$ of prior range | Not informed by data | Fall back to scenario analysis: $\beta_{\text{duck}} \in \{0.5, 1.0, 1.5\}$ |
| $p_{\text{mov}}$ posterior $\approx$ prior | Not identified | Report non-identifiability; consider fixing $p_{\text{mov}}$ |
| $t_0$–$\delta$ posterior correlation $> 0.7$ | Onset-decay trade-off | Tighten $t_0$ prior based on first-case timing |

---

## Software

### Recommended: Turing.jl

Turing.jl is the primary framework for this model. The key capabilities that make it suitable:

- **Custom likelihood** via `@addlogprob!` — the FOI survival likelihood does not fit standard distribution families and must be coded directly.
- **Composable samplers** — if the model later requires mixed discrete/continuous inference (Option B above), Turing supports `Gibbs(NUTS(...), MH(...))` without changing the model specification.
- **Julia ecosystem** — spatial preprocessing (distances, zone membership), data handling, and likelihood computation all in one language, avoiding inter-language overhead.
- **ForwardDiff integration** — automatic differentiation is built in; no manual gradient derivation required.

### Alternatives considered

| Framework | Strengths | Limitation for this model |
|---|---|---|
| Stan | Mature NUTS implementation; excellent diagnostics | No native discrete latent variable support; would require marginalisation (Option C) from the start. Julia-to-Stan bridge adds friction. |
| Custom Julia MCMC | Maximum flexibility for likelihood structure | Sacrifices tested NUTS implementation, adaptation, and diagnostic tools. Significant development effort for little gain given Turing exists. |

### Dependencies

The core inference stack:

```julia
using Turing              # Model specification, NUTS sampler
using Distributions       # Prior distributions
using MCMCDiagnosticTools # R-hat, ESS, convergence checks
using MCMCChains          # Chain storage, summary statistics
```

Preprocessing and data handling (shared with other steps):

```julia
using CSV, DataFrames     # Data loading
using GeoJSON             # HRZ boundary
using NearestNeighbors    # Efficient distance queries for spatial kernel
```

---

## Computational strategy

### Precomputation (done once before MCMC)

These quantities are constant across all MCMC iterations and should be computed and cached during data preparation.

| Quantity | Description | Cost |
|---|---|---|
| Pairwise distances | $d_{ij}$ for all farm pairs within a cutoff radius | $O(N^2)$ but sparse; use KD-tree for radius query |
| Neighbour lists | For each farm, the set of farms within the kernel's effective range | From distance matrix; stored as flat arrays for cache efficiency |
| HRZ membership | Static wild-bird spillover zone; binary flag per farm from `hrz_32626.geojson` | $O(N)$ point-in-polygon |
| Activity status | Per-farm, per-day indicator of bird presence | From `activity.csv`; $O(N \times T)$ |
| Movement network | Source–destination–date tuples, filtered by zone restrictions | From `movement.csv`; prefilter by date range |
| Regulated zone status | Dynamic movement-restriction zones (3 km protection + 10 km surveillance) around confirmed cases, 28-day duration. Distinct from the static HRZ; a farm can be in both, but regulated-zone status takes precedence in $p_{\text{eff}}$ calculations (step 06, Module B) | Derived from confirmed case dates and locations; recomputed as new cases are confirmed |
| Bulk spillover bins | Farm counts by (HRZ status $\times$ species $\times$ active), per day | Enables $O(1)$ non-case survival computation per bin instead of $O(N)$ per farm |

### Per-iteration optimisations

- **Sparse spatial hazard**: only compute $K(d_{ij})$ for farm pairs where $i$ is currently infectious. With $\sim 10$ farms infectious on any given day and $\sim 50$ neighbours each within the kernel's effective range, this is $\sim 500$ kernel evaluations per time step rather than $N^2$.
- **Bulk non-case survival**: the spillover survival contribution for non-case farms in the same (HRZ, species, activity) bin is identical. Aggregate into bin counts and compute once per bin.
- **Incremental hazard updates**: when the sampler proposes new parameter values, recompute only the hazard components that changed. For A+C fits, Module B is absent and spatial hazard computation is skipped entirely.

### Parallelisation

Run 4 chains in parallel using Julia's multi-threading or `Distributed` workers. Each chain is independent; no inter-chain communication until convergence diagnostics. On a 4-core machine, wall-clock time is approximately $1/4$ of single-chain time.

---

## Sensitivity analysis plan

Sensitivity analysis addresses the assumptions that are fixed rather than estimated. Each analysis reruns the full MCMC and compares posteriors for the core parameters ($\phi_{\text{hrz}}$, $\beta$, $\alpha$) and derived quantities (e.g. proportion of cases attributable to spillover vs transmission).

### Delay sensitivity (highest priority)

Grid over $\mu_E \in \{2, 3, 4, 5\}$ and $\mu_{ID} \in \{3, 4, 5, 6, 7\}$ (20 combinations). The $D \to C$ delay is fixed at 2 days from data. Report:
- Shift in posterior medians for $\beta$, $\alpha$, $\phi_{\text{hrz}}$, $\phi_{\text{non}}$.
- Whether $\beta_{\text{duck}}$ informativeness changes.
- Whether the $\beta$–$\alpha$ correlation exceeds the 0.8 reparameterisation trigger.

If posteriors are stable across the grid, the fixed-delay assumption (Option A) is defensible and no transition to Option C is needed.

### Prior sensitivity

For each estimated parameter, multiply the prior standard deviation by 2 and refit. Report changes in posterior medians and 95% CrIs. Parameters whose posteriors shift substantially are prior-sensitive and warrant more careful prior elicitation or additional data.

### Kernel comparison

Fit both DAG A (exponential) and DAG B (Cauchy) as defined in step 05 (§Candidate DAGs). Compare via:
- Posterior predictive checks on spatial clustering of cases.
- WAIC or LOO-CV if posterior predictive discrimination is unclear.

### Detection stationarity

Refit with two delay regimes (step 05, §Stationarity of detection delays):
- **Stationary**: constant $d = 10.5$ days throughout.
- **Non-stationary**: $d = 10.5$ days before 1 Jan 2026, $d = 8.5$ days after (reflecting heightened surveillance and faster recognition).

Report whether key conclusions ($\beta_{\text{duck}}$, spillover/transmission partitioning) are robust to this change.

### Zone effectiveness

Refit with $\varepsilon \in \{0.3, 0.5, 0.7\}$ (step 05, §Fixed parameters). This parameter reduces hazard for farms inside active regulated zones (dynamic 3 km protection + 10 km surveillance zones around confirmed cases; distinct from the static HRZ) and is otherwise unidentifiable from the data.

---

## Summary of decisions

| Consideration | Decision | Rationale |
|---|---|---|
| Likelihood | FOI survival (analytical) | Differentiable, exact, sparse; no simulation needed |
| Latent infection times | Option A (fixed delays) initially | 8 continuous params only; NUTS-compatible; grid sensitivity to validate |
| Transition pathway | A → C → B | Increase complexity only if sensitivity analysis demands it |
| Sampler | NUTS via Turing.jl | Low dimension, differentiable likelihood, proven in iteration 1 |
| AD mode | ForwardDiff (forward-mode) | Optimal for $\leq 10$ parameters; chunk size = param count |
| Chains | 4 × 1000 samples, 500 warmup | Standard; sufficient for $\hat{R}$ and ESS diagnostics |
| Initialisation | Grid/random search for approximate MAP | Reduces warmup; avoids poor starting regions |
| Software | Turing.jl | Custom likelihood, composable samplers, Julia ecosystem |
| Sensitivity | Delays (20-point grid), priors (2× width), kernels (exp vs Cauchy), detection stationarity, zone effectiveness | Covers all major fixed assumptions |

---
