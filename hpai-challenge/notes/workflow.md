# HPAI Modelling Workflow — following Abbott et al.

## Stage 0: Data source characterisation

### Available data sources

**1. Confirmed case notifications** (`cases.csv`)
- Metadata: Routine passive surveillance + pre-shipment testing + contact tracing
- Scope: All confirmed HPAI cases on Jolly Island; farm-level
- Resolution: Daily (date_suspicious, date_confirmed, cull_start); farm-level spatial
- Quality: Complete for confirmed cases; detection method recorded; 96/103 passive, 6 preshipment, 1 contact tracing
- Utility: Directly informs transmission dynamics, spatial spread, species differences
- Practical: 103 cases over 24 days (22 Dec – 13 Jan)

**2. Mortality ledgers** (3 example farms in narrative PDF)
- Metadata: Daily bird mortality counts from 3 infected farms
- Scope: Only 3 farms — illustrative, not systematic
- Resolution: Daily; per-farm
- Quality: High for these 3 farms; consistent exponential growth pattern
- Utility: Informs within-farm infectiousness growth rate (r ≈ 1.0/day)
- Practical: Very limited sample; used to fix r, not for inference

**3. Population register** (`population.csv`)
- Metadata: Complete farm census
- Scope: All 9,160 poultry farms; farm locations, species, production type, capacity
- Resolution: Farm-level; UTM coordinates; species (chicken/duck), production type
- Quality: Complete census, high quality
- Utility: Defines the susceptible population, spatial structure, species composition
- Practical: Static (no temporal dimension)

**4. Activity data** (`activity.csv`)
- Metadata: Which farms have birds present on which dates
- Scope: All farms, full period
- Resolution: Daily; per-farm
- Quality: Complete
- Utility: Defines which farms are at risk on each day
- Practical: Essential for correct susceptible-days calculation

**5. Movement data** (`movement.csv`)
- Metadata: Recorded movements between farms
- Scope: All recorded movements during the period
- Resolution: Daily; source-destination pairs
- Quality: Complete for recorded movements; may not capture all contacts
- Utility: Potential transmission pathway; limited signal (only 6 preshipment detections)
- Practical: 7,187 movements over the period

**6. Preventive culling records** (`prev_culls.csv`)
- Metadata: Farms culled preventively (not confirmed cases)
- Scope: 52 entries, 12 with recorded cull dates
- Resolution: Farm-level; daily
- Quality: Incomplete (many missing cull dates)
- Utility: Removes farms from susceptible pool; important for correct survival calculation
- Practical: Small numbers

**7. High-risk zone polygon** (`hrz_32626.geojson`)
- Metadata: Spatial boundary of the HRZ
- Scope: Island-wide
- Resolution: Polygon boundary
- Quality: Definitive (regulatory boundary)
- Utility: Classifies farms for spillover rate and testing regime
- Practical: 1,962 of 9,160 farms in HRZ

---

## Stage 1: Research question & target estimands

### Research questions (from the challenge)

**Q1** (descriptive): General description of the epidemic — epi curve, spatial distribution, species breakdown. *Done — output already generated.*

**Q2** (prediction): 4-week temporal and spatial prediction of new cases by production type.
- Target estimand: Daily incidence of new confirmed cases by date, species, and production type, with uncertainty.
- Spatial: Per-farm or per-district infection probability.

**Q3** (mechanistic): Relative contribution of galliform vs palmiped farms to spread.
- Target estimand: Species-specific susceptibility ratio β_duck/β_chicken.
- Secondary: Proportion of total transmission attributable to each species.

**Q4** (counterfactual): Should preventive culling focus on chicken farms only?
- Target estimand: Difference in epidemic trajectory (cases, removals) between baseline and chicken-only preventive culling.

**Q5** (counterfactual): Would dropping preventive culling and speeding up reactive culling be better?
- Target estimand: Difference in epidemic trajectory between baseline and faster-reactive scenario.

### Scope and limitations
- Farm-level model (not individual bird)
- Single outbreak on synthetic island
- 44-day observation period, 28-day prediction horizon
- No serological or genomic data available

---

## Stage 2: Process DAG

### Latent states
Each farm j on day t can be in one of:
- **S**: Susceptible (has birds, not yet infected)
- **E**: Exposed/latent (infected but not yet infectious)
- **I**: Infectious (shedding virus, can transmit)
- **D**: Detected (suspicion raised, awaiting confirmation)
- **C**: Confirmed (lab-confirmed)
- **R**: Removed (culled)

### Transitions
```
S → E: infection event (spillover + spatial + movement hazard)
E → I: after latent period
I → D: after detection delay
D → C: after D→C delay
C → R: after C→R delay (reactive culling)
S → R: preventive culling (if within 1km of confirmed case, from 1 Jan)
```

### Key process components

**a) Wild-bird spillover**
External force of infection from wild birds. Needs:
- Spatial heterogeneity (HRZ vs non-HRZ)
- Temporal profile (onset, peak, decline)
- Species-specific susceptibility

**b) Farm-to-farm spatial transmission**
Local spread via shared resources, aerosol, fomites, etc. Needs:
- Distance-dependent kernel (form TBD)
- Farm-level infectiousness profile (function of time since infection)

**c) Movement-mediated transmission**
Direct transport of infected birds. Needs:
- Movement network (from data)
- Probability of transmission per movement
- Testing/zone blocking

**d) Within-farm dynamics**
Infectiousness builds over time as within-farm prevalence grows.
- Informed by mortality ledgers: exponential growth, r ≈ 1.0/day
- Open question: smooth ramp-up vs hard latent period cutoff

**e) Interventions**
- Reactive culling (C→R delay from data)
- Regulated zones (10km, 28 days, reduce hazard)
- Preventive culling (1km, from 1 Jan)
- Pre-shipment testing (HRZ, sensitivity 0.9)

### Open modelling decisions (for Stage 5)
1. Spillover temporal profile: what shape? How to parameterise onset?
2. Spatial kernel: exponential, Gaussian, power-law?
3. Infectiousness profile: hard latent cutoff or smooth ramp-up?
4. Which delays to fix vs estimate?
5. Movement transmission: fixed or estimated?

---

## Stage 3: Data source selection

For the process model, we use:
- **cases.csv**: Primary data for inference (infection times, confirmation dates, spatial locations)
- **population.csv**: Susceptible population structure
- **activity.csv**: Time-varying susceptibility (which farms have birds)
- **movement.csv**: Movement transmission pathway
- **prev_culls.csv**: Removal of farms from susceptible pool
- **hrz_32626.geojson**: HRZ classification for spillover

The mortality ledgers (3 farms) inform the fixed within-farm growth rate but are too few for formal inference.

---

## Stage 4: Observation DAG

### How observed data relates to latent process

The key observation model: we observe *confirmations* (date_confirmed), not *infections*. The mapping:

```
Infection (latent) → [latent period] → Infectious → [detection delay] → Suspicion
→ [D→C delay] → Confirmation (observed) → [C→R delay] → Removal (observed)
```

For passive cases (96/103): we observe date_suspicious AND date_confirmed.
For preshipment cases (6): we observe date_confirmed only.
For contact tracing cases (1): we observe date_confirmed only.

**Observation biases:**
- Only confirmed cases are observed — subclinical infections missed
- Detection probability may vary by:
  - Species (chickens show clinical signs earlier?)
  - Farm size
  - Proximity to other cases (heightened surveillance)
  - Time (surveillance intensity increases during outbreak)
- Preventive culls remove farms that might have been infected but never confirmed

**For our model:**
The FOI likelihood conditions on the observed infection times (back-calculated from confirmation dates). This means:
- We assume all confirmed cases are true cases (no false positives)
- We assume no other farms were infected but undetected (strong assumption)
- The back-calculation depends on the fixed delay parameters

---

## Stage 5: Refine the model DAGs

The goal: reconcile the process DAG (Stage 2) with what the data can support (Stages 3-4). Resolve the open modelling decisions with domain knowledge and simplicity-first reasoning.

### Decision 1: Spillover temporal profile

**Domain knowledge:**
- Wild birds "migrate in early winter with stragglers seen through February"
- First confirmed case: 22 Dec 2025 (back-calculated infection ~12-13 Dec)
- Index case detected ~30 days before narrative was written
- Virus believed introduced by wild birds near the Southern region (HRZ)

**Options:**
a) Fixed ramp-up rate (iteration 1: 0.02/day) — too gentle, ψ already 66% of peak on day 1
b) Logistic onset with estimated midpoint and steepness — flexible, but 2 extra params
c) Piecewise: zero before estimated onset day t₀, then quick ramp to peak, then decay
d) Step function: zero before t₀, constant after — simplest

**Decision:** Start with option (c) — a piecewise profile with an estimated onset day t₀ and a fast ramp-up. This captures the key feature (no spillover before migration arrives) while remaining simple. Specifically:
- Before t₀: ψ(t) = 0
- t₀ to t_peak: linear or exponential ramp-up (short window, ~5-10 days)
- After t_peak: exponential decay with rate δ

The onset day t₀ is the critical missing piece from iteration 1. It should be estimated with an informative prior (say Normal(15, 5) on the day scale, centred around mid-December).

Alternative: could use a log-normal pulse shape with estimated peak time and width — single parametric form, no piecewise discontinuity.

### Decision 2: Spatial kernel form

**Domain knowledge:**
- HPAI local spread via shared equipment, personnel, aerosol, fomites
- Generally short-range (<5 km is most important), but occasional longer-range jumps
- Dutch 2003 and Italian 1999-2000 HPAI outbreaks used exponential or power-law kernels

**Decision:** Start with exponential kernel h₀·exp(-d/α). Simplest, well-established. If posterior predictive checks show spatial patterns not captured, consider power-law (heavier tail) in a later iteration. The workflow says: start simple, iterate.

### Decision 3: Infectiousness profile

**Domain knowledge:**
- Mortality ledgers show exponential growth in within-farm deaths (r ≈ 1.0/day)
- Farm-level infectiousness should reflect viral output, which correlates with prevalence
- The latent period represents time from virus introduction to sufficient within-farm prevalence to transmit — not a sharp biological boundary

**Options:**
a) Hard latent cutoff + exponential growth (iteration 1) — simple but discontinuous
b) Gamma CDF ramp-up + exponential growth — smooth, biologically motivated
c) Logistic ramp-up — smooth but different functional form
d) Full Gamma-distributed incubation with infectiousness integral — most rigorous

**Decision:** Keep the hard latent cutoff for now (option a). The cutoff is at the farm level (not individual bird) and represents when enough birds are infected to generate detectable viral shedding. The exact shape matters less than getting the delay parameters right. Revisit if needed.

### Decision 4: Delay parameters

**Domain knowledge:**
- Farm-level latent period: virus introduction → enough infected birds to transmit. Literature: 2-5 days for HPAI at farm level.
- Detection delay: onset of detectable mortality → farmer raises suspicion. Narrative: "sudden increase in mortality" triggers action. Could be 3-7 days.
- D→C delay: sampling → lab confirmation. Narrative: "~2 days". Data confirms median 2 days.
- C→R delay: confirmation → culling. Narrative: "within subsequent few days"; capacity constraints from 6 Jan. Data-driven per case.

**Decision:** Fix μ_E and μ_ID because they enter through integer rounding of infection day assignments, giving zero gradient for HMC/NUTS. This is a well-known problem in spatial FOI models. Use literature-informed values:
- μ_E = 3.5 days (midpoint of 2-5 day range)
- μ_ID = 5.0 days (midpoint of 3-7 day range)
- D→C = 2 days (from data)
- C→R = per case from data (or Gamma(2.0, 1.5) for predictions)

Conduct sensitivity analysis over μ_E ∈ {2, 3, 4, 5} and μ_ID ∈ {3, 4, 5, 6, 7} in a later iteration.

### Decision 5: Movement transmission

**Domain knowledge:**
- Only 6 preshipment detections out of 103 cases
- 7,187 movements recorded — most between non-infected farms
- Pre-shipment testing (sensitivity 0.9) in HRZ should block most infectious movements
- Movements from regulated zones blocked entirely

**Decision:** Fix p_mov = 0.01 initially. The signal is weak (6 cases) and movements are heavily confounded with spatial proximity. Include the movement pathway for mechanistic completeness but don't try to estimate p_mov — it would be poorly identified. Revisit if posterior predictive checks suggest missing transmission events that movements could explain.

### Decision 6: Zone biosecurity effectiveness

**Domain knowledge:**
- Surveillance zones (10 km, 28 days) from narrative
- Heightened biosecurity, movement restrictions in zones
- No direct data on effectiveness

**Decision:** Fix ε = 0.5 (zones reduce hazard by 50%). This is a common assumption in the literature. Could do sensitivity analysis later.

### Refined process model summary

**Estimated parameters (target: 5-7):**
1. h₀ — spatial transmission rate
2. α — spatial range
3. β_duck — duck susceptibility relative to chicken
4. φ_hrz — HRZ spillover rate at peak
5. φ_non — non-HRZ spillover rate at peak
6. δ — post-peak spillover decay rate
7. t₀ — spillover onset day (NEW — addresses the key iteration 1 problem)

**Fixed parameters:** μ_E, μ_ID, D→C delay, r, p_mov, ε, zone radius/duration, preventive cull radius/start date, σ_test

### What changes from iteration 1?
- **Added:** spillover onset parameter t₀ (previously hardcoded ramp from day 1)
- **Changed:** spillover profile from gentle exp ramp-up to zero-before-onset + ramp
- Everything else stays the same for this iteration — start simple, iterate

---

## Stage 6: Modularisation

Following Abbott et al., decompose the model into independently verifiable submodels. Validate each before combining.

### Module structure

**Module 1: Spillover-only model**
- Fix h₀ = 0, p_mov = 0 (no farm-to-farm transmission)
- Estimate: φ_hrz, φ_non, δ, t₀, β_duck
- Data: cases.csv, population.csv, activity.csv, hrz_32626.geojson, prev_culls.csv
- Purpose: Can spillover alone explain the observed case pattern? What are the spillover parameters in isolation?
- Expected: will over-predict late cases (can't explain clustering without spatial kernel) or under-predict (needs spatial amplification). Either way, establishes baseline spillover rates.

**Module 2: Spillover + spatial transmission**
- Add h₀, α to Module 1
- Estimate: h₀, α, φ_hrz, φ_non, δ, t₀, β_duck (7 params)
- Data: same as Module 1
- Purpose: Full transmission model without movement pathway
- Expected: should capture most dynamics; movements contribute <6% of cases

**Module 3: Full model (add movements)**
- Add p_mov (fixed or estimated) to Module 2
- Purpose: Complete model for prediction and counterfactuals
- May not need separate fitting if p_mov is fixed

### Validation at each module
For each module, apply Stage 9 checks:
1. Prior predictive check: do simulated epidemics look plausible?
2. Synthetic data recovery: can the model recover known parameters?
3. Fit to real data
4. Posterior predictive check: does the fitted model reproduce observations?

---

## Stage 7: Inference and computation choices

### Model complexity assessment
- 5-7 continuous parameters (low-dimensional)
- No latent states to integrate out (infection times back-calculated from data)
- Likelihood is analytical (FOI survival — no simulation)
- Likelihood is differentiable w.r.t. all estimated parameters (verified in iteration 1 for the 6-param model)

### Inference method
**NUTS (No-U-Turn Sampler)** via Turing.jl — the natural choice because:
- Low parameter dimension (5-7)
- Analytical, differentiable likelihood
- Iteration 1 demonstrated excellent performance (0.19s/sample, ESS 272-595 for 6 params)
- ForwardDiff automatic differentiation (forward-mode appropriate for ≤10 parameters)

### Implementation considerations
- Chunk size for ForwardDiff: match parameter dimension (6 or 7)
- Warmup: 500 iterations, sampling: 500-1000 iterations
- Multiple chains for convergence diagnostics (R-hat)
- Grid search for MAP initialisation to reduce warmup time

### Diagnostics
- R-hat < 1.01 for all parameters
- ESS > 100 (ideally > 400)
- No divergent transitions
- Tree depth not maxing out (target max 10)
- Visual trace plot inspection

---

## Stage 8: Implementation

### Language and framework
- Julia + Turing.jl (probabilistic programming language)
- Custom likelihood via `@addlogprob!`
- Generative model capability: same code structure for prior predictive, fitting, and prediction

### Code structure
- `prepare_model_data()`: load and preprocess all data, compute neighbours
- `foi_loglik()`: evaluate FOI survival log-likelihood for given parameters
- `@model hpai_model()`: Turing model specifying priors + likelihood
- `simulate_forward()`: generative simulation for prediction/PPC
- `run_predictions()`: batch forward simulations across posterior draws
- Modular: spillover/spatial/movement components in separate functions for easy on/off

### Computational optimisation
- Bulk spillover survival aggregation (~9000 farms → 8 category counts)
- Sparse spatial hazard tracking (only farms near infectious cases)
- Flat neighbour arrays for cache-friendly access
- Target: <2ms per likelihood evaluation

---

## Stage 9: Model specification and validation

### 9.1 Model family
Discrete-time (daily) stochastic farm-level model with:
- Discrete states (S, E, I, D, C, R) per farm
- Continuous transmission parameters
- FOI survival likelihood (point-process)

### 9.2 Prior specification

**Prior choices (to be refined by prior predictive checks):**

| Parameter | Prior | Rationale |
|-----------|-------|-----------|
| h₀ | LogNormal(log(1e-4), 1.5) | Weak spatial transmission; wide range |
| α | LogNormal(log(3500), 0.5) | ~3.5 km range, informed by HPAI kernel literature |
| β_duck | Beta(2, 8) | Mean 0.2; ducks ~5× less susceptible than chickens (observed attack rates) |
| φ_hrz | LogNormal(log(1e-3), 1.0) | Daily per-farm spillover in HRZ |
| φ_non | LogNormal(log(1e-4), 1.0) | ~10× lower outside HRZ |
| δ | Exponential(mean=0.02) | Slow post-peak decay; wide right tail |
| t₀ | Normal(15, 5) | Mid-December onset; ±10 days at 95% |

**Prior predictive check protocol:**
1. Draw N=500 parameter sets from priors
2. For each, simulate an epidemic from day 1 (spillover only at first, then full model)
3. Check: total cases, timing of first case, spatial distribution, species ratio
4. Criteria for acceptable priors:
   - 95% of simulations produce 10-500 total cases (plausible epidemic range)
   - First case between day 10 and day 30 (early-to-mid December)
   - Duck:chicken case ratio roughly 1:3 to 1:5

### 9.3 Computational validation

**Step 1: Synthetic data recovery**
- Choose "true" parameter values (e.g., posterior medians from iteration 1, or domain-informed values)
- Simulate one synthetic epidemic from these parameters
- Fit the model to this synthetic data
- Check: can we recover the true parameters within 95% CrI?

**Step 2: Simulation-based calibration (SBC)**
- Draw parameters from prior, simulate data, fit, check coverage
- 95% CrI should contain the true value ~95% of the time
- Computationally expensive — do for Module 1 (5 params) first

**Step 3: Fit to real data**
- After passing Steps 1-2, fit to observed cases
- Check diagnostics (R-hat, ESS, divergences, tree depth)

### 9.4 Posterior evaluation

**Posterior predictive checks:**
- Simulate 200 epidemics from posterior draws (same as prediction, but over observation period)
- Compare to observed data:
  1. Total case count: posterior median vs 103
  2. Epidemic curve shape: timing, peak, decline
  3. Spatial distribution: where do simulated cases cluster?
  4. Species ratio: chicken vs duck case proportions
  5. First case timing: when do simulated epidemics start?
- Acceptable if observed data falls within 95% CrI for key summaries

**Sensitivity analysis:**
- Vary fixed parameters (μ_E, μ_ID, p_mov, ε) and re-fit
- Check if conclusions (Q3-Q5) change materially

---

## Stage 10: Data integration

Not directly applicable here — we have a single joint model rather than multiple independent modules requiring ensembling. However, the sequential module-building (Module 1 → 2 → 3) serves the same purpose: validate simpler models before adding complexity.

---

## Implementation plan (ordered tasks)

1. ☐ Implement spillover onset parameter t₀ in the model
2. ☐ Module 1: spillover-only model — prior predictive check
3. ☐ Module 1: synthetic data recovery test
4. ☐ Module 1: fit to real data + posterior predictive check
5. ☐ Module 2: add spatial transmission — prior predictive check
6. ☐ Module 2: synthetic data recovery test
7. ☐ Module 2: fit to real data + posterior predictive check
8. ☐ Module 3: add movements (if needed) + validate
9. ☐ Generate final Q2-Q5 outputs from validated model
10. ☐ Sensitivity analysis on fixed parameters
