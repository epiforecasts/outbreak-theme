# HPAI Spatial Transmission Model — Jolly Island

## Overview

A Bayesian farm-level spatial transmission model for highly pathogenic avian
influenza (HPAI). Each farm is a single epidemiological unit — either
susceptible, infected (with a time-varying infectiousness), or removed
(culled). Two transmission pathways are modelled: wild-bird spillover and
local spatial spread (including movement-mediated transmission as a
sensitivity analysis).

The model is implemented in Julia using Turing.jl for Bayesian inference
and a custom forward simulation engine for prediction.

- **Observation period**: 1 Dec 2025 – 13 Jan 2026 (44 days)
- **Prediction horizon**: 14 Jan – 10 Feb 2026 (28 days)

---

## 1. Farm-level infectiousness

Once infected, a farm's infectiousness rises over time as within-farm
prevalence builds, then is truncated by culling. The infectiousness profile
w(τ), where τ is time since infection, follows a saturating form:

    w(τ) = 0                              if τ < τ_min
    w(τ) = 1 − exp(−r · (τ − τ_min))     if τ ≥ τ_min

where r = 1.0/day is the within-farm growth rate and τ_min = 1 day is the
minimum delay before a farm contributes to onward transmission. This is
informed by mortality ledgers showing consistent exponential mortality
growth across the three example farms.

Note: τ_min and μ_E (§2) serve different purposes. τ_min is the hard
latent period in the transmission kernel — farms cannot transmit before
τ_min days post-infection. μ_E = 3.5 days is the mean incubation period
used to back-calculate infection times from observed suspicion dates.
Both are fixed; μ_E is not used in the hazard computation.

Because w(τ) is increasing, an extra day before culling adds
disproportionately to transmission.

---

## 2. Timeline of a single case

For each infected farm i, the timeline is:

    t_inf ──── t_inf + μ_E ──── t_susp ──── t_conf ──── t_cull
      │            │               │            │           │
    infection   infectious     suspicion   confirmation  removal
                  starts        raised       (lab)       (culled)

Key intervals:
- **Latent period** μ_E = 3.5 days (fixed). Infection → infectious.
- **Detection delay** μ_ID = 5.0 days (fixed). Infectious → suspicion raised.
  Together μ_E + μ_ID = 8.5 days from infection to suspicion.
- **D→C delay** = 2 days (fixed). Suspicion → confirmation.
- **C→R delay**: from data per case (cull_start − date_confirmed); default
  3 days if missing.
- **Compound delay** = ⌈μ_E + μ_ID + D→C⌉ = 11 days (infection → confirmation).

Both μ_E and μ_ID are fixed because they enter the likelihood through
integer rounding of infection day assignments, which gives zero gradient
for NUTS. Fixing them at literature-informed values is standard practice
for spatial FOI models where delay parameters are confounded with
transmission parameters.

For observed cases, we know t_susp (for passive cases), t_conf, and t_cull
from the data. Given the fixed delays, we back-calculate:
- t_inf = t_susp − μ_E − μ_ID  (passive cases, using observed suspicion date)
- t_inf = t_conf − μ_E − μ_ID − 2  (other detection methods)
- t_infectious_start = t_inf + μ_E

The infectious window runs from t_infectious_start to t_cull.

---

## 3. Force of infection

The daily hazard for susceptible farm j on day t:

    λ_j(t) = β_sj · [ spillover_j(t) + spatial_j(t) + movement_j(t) ] · zone_j(t)

where zone_j(t) = (1 − ε) if farm j is inside a regulated zone on day t,
and 1.0 otherwise (see §4).

### 3.1 Species-specific susceptibility

    β_sj = 1.0 if chicken, β_duck if duck

β_duck is estimated (prior mean ~0.2, reflecting the ~5× lower observed
attack rate in ducks vs chickens).

### 3.2 Wild-bird spillover

    spillover_j(t) = φ_j · ψ(t)

where:
- φ_j = φ_hrz if farm j is in the high-risk zone, φ_non otherwise
- ψ(t) is a Bateman function capturing the rise and decay of wild-bird
  pressure:

      ψ(t) = 0                                           if t < t₀
      ψ(t) = [(1 − exp(−σ·τ)) · exp(−δ·τ)] / ψ_peak    if t ≥ t₀

  where τ = t − t₀, and ψ_peak normalises the profile to peak at 1.0:

      ψ_peak = [σ/(σ+δ)] · [δ/(σ+δ)]^(δ/σ)

The Bateman function models a process that rises with rate σ (onset of
wild-bird migration) and decays with rate δ (waning environmental
pressure). The onset time t₀ and both rates σ and δ are estimated from
data.

### 3.3 Local spatial transmission

    spatial_j(t) = Σ_{i ∈ infectious(t)}  β · exp(−d_ij / α) · w(t − t_inf_i)

Sum over all farms i that are infectious at time t. The exponential kernel
β · exp(−d/α) captures distance-dependent local spread (shared equipment,
personnel, aerosol, fomites). The infectiousness profile w(τ) weights
each infectious farm's contribution by its current viral output.

Only farms within 50 km are considered (computational cutoff).

### 3.4 Movement-mediated transmission

For each movement on day t where the source farm is infectious and the
movement is not blocked by zone restrictions:

    movement_j(t) = p_mov · w(t − t_inf_source)

Pre-shipment testing in the HRZ intercepts infectious source farms with
probability σ_test = 0.9, reducing the effective movement hazard:
- p_eff = p_mov × (1 − σ_test) for HRZ source farms
- p_eff = p_mov for non-HRZ source farms

Movements from farms inside a regulated zone are blocked entirely.

During the prediction period, movements are generated from a fitted
gravity model (see §7.3).

### 3.5 Infection probability

    P(farm j infected on day t) = 1 − exp(−λ_j(t))

---

## 4. Interventions

All intervention parameters are fixed from the narrative:

- **Reactive culling**: confirmed farms culled after observed C→R delay
  (from data). For new cases in the prediction period, a fixed buffer
  of 3 days after confirmation is used.
- **Regulated zones**: 10 km surveillance zones around confirmed farms for
  28 days after confirmation. Farms inside zones have reduced hazard
  (multiplied by 1 − ε = 0.5) and cannot ship movements.
- **Preventive culling** (from 1 Jan 2026): active susceptible farms within
  1 km of confirmed cases are removed.
- **Pre-shipment testing** (HRZ only): infectious source farms detected
  with sensitivity σ_test = 0.9, blocking the movement.

For the inference likelihood, interventions affect the hazard computation
(zone effects reduce λ_j) and the infectious window (culling truncates it).

---

## 5. Parameters

### Spillover model — estimated (6 parameters)

The primary model used for inference and Q2 prediction.

| Parameter | Symbol | Description | Prior |
|-----------|--------|-------------|-------|
| Spillover onset | t₀ | Day wild-bird pressure begins | truncated Normal(15, 5) on [1, 44] |
| Spillover rate (HRZ) | φ_hrz | Daily spillover hazard in HRZ | LogNormal(log(1e-3), 1.0) |
| Spillover rate (non-HRZ) | φ_non | Daily spillover hazard outside HRZ | LogNormal(log(1e-4), 1.0) |
| Spillover decay | δ | Post-peak decay rate (/day) | Exponential(mean=0.02) |
| Duck susceptibility | β_duck | Relative to chicken (=1) | Beta(2, 8) |
| Spillover rise rate | σ | Bateman onset rate (/day) | LogNormal(log(0.3), 1.0) |

### Full model — additional parameters (3)

The full model adds spatial and movement transmission. These parameters
are not independently identifiable from spillover in this dataset (see §7.2),
so they are used only for sensitivity analysis with values fixed at the
full model's posterior medians.

| Parameter | Symbol | Description | Prior / Fixed value |
|-----------|--------|-------------|---------------------|
| Spatial transmission rate | β | Kernel amplitude | LogNormal(log(1e-4), 1.5) / fixed: 0.006 |
| Spatial range | α | Kernel decay distance (metres) | LogNormal(log(3500), 0.5) / fixed: 4572 |
| Movement transmission | p_mov | Per-movement infection probability | Beta(2, 20) / fixed: 0.076 |

### Fixed

| Parameter | Value | Source |
|-----------|-------|--------|
| Latent period μ_E | 3.5 days | HPAI literature; not identifiable (zero gradient through integer rounding) |
| Detection delay μ_ID | 5.0 days | HPAI literature; not identifiable (zero gradient through integer rounding) |
| β_chicken | 1.0 | Reference category |
| Infectiousness growth rate r | 1.0/day | Mortality ledgers (3 farms) |
| Hard latent period τ_min | 1 day | Minimum time before infectiousness |
| D→C delay | 2 days | Data (median) |
| C→R delay (observed cases) | Per case | Data (cull_start − date_confirmed); default 3 days if missing |
| C→R delay (predicted cases) | 3 days | Fixed buffer (REMOVAL_BUFFER) |
| Zone radius | 10 km | Narrative (surveillance zone) |
| Zone duration | 28 days | Narrative |
| Zone biosecurity ε | 0.5 | Assumption |
| Preventive cull radius | 1 km | Narrative |
| Preventive cull start | 1 Jan 2026 | Narrative |
| Pre-shipment sensitivity σ_test | 0.9 | Assumption |
| Neighbour cutoff | 50 km | Computational |

---

## 6. Likelihood

A force-of-infection (survival) likelihood, evaluated over the observation
period (1 Dec 2025 – 13 Jan 2026, T = 44 days).

For each day t = 1, …, T:

1. Back-calculate infection times for all 103 observed cases from their
   suspicion/confirmation dates using the fixed delay parameters.

2. Determine which cases are infectious on day t (between t_infectious
   and t_cull).

3. For each active susceptible farm j, compute hazard λ_j(t) from all
   active transmission pathways plus zone effects.

4. **Case farms infected on day t** contribute:
       log(1 − exp(−λ_j(t)))
   This is the log-probability of the observed infection event.

5. **Susceptible farms that survived day t** contribute:
       −λ_j(t)
   This is the log-survival probability.

The total log-likelihood:

    ℓ(θ) = Σ_t [ Σ_{j ∈ infected(t)} log(1 − e^{−λ_j(t)})
                + Σ_{j ∈ survived(t)} (−λ_j(t)) ]

This is a standard spatial point-process likelihood — fully analytical,
no simulation required.

### Computational optimisation

The likelihood evaluation exploits the sparsity of the problem:

1. **Bulk spillover survival** (~9,000 non-case, non-prev-culled farms):
   aggregated into 8 category counts per day (species × HRZ × zone status),
   avoiding per-farm loops.

2. **Per-case and per-prev-cull corrections**: only ~115 farms need
   individual spillover accounting (for their pre-infection/pre-cull days).

3. **Sparse spatial hazard**: only farms within 50 km of an infectious
   case receive non-zero spatial hazard on any given day; tracked via
   a "touched" list to avoid full-array operations.

4. **Flat neighbour arrays**: case-farm neighbour indices and distances
   stored contiguously for cache-friendly access.

This gives ~1.4 ms per likelihood evaluation, enabling efficient NUTS
sampling with ForwardDiff automatic differentiation (chunk size = 6).

---

## 7. Inference and prediction

### 7.1 Inference

The spillover model (6 parameters) is used for inference. The Turing
`@model` specifies priors and evaluates the FOI log-likelihood via
`@addlogprob!`. Sampled with NUTS (No-U-Turn Sampler) using ForwardDiff
for gradient computation.

Initialisation: a coarse-then-fine grid search over parameter space
identifies an approximate MAP estimate, which is passed to NUTS as the
starting point to reduce warmup time.

### 7.2 Model selection

Both the spillover-only and full (9-parameter) models were fitted. The
full model includes spatial and movement transmission but its additional
parameters (β, α, p_mov) are not independently identifiable from the
spillover parameters in this dataset: adding spatial transmission pulls
down the spillover rates, and both pathways explain the same observed
infections.

Posterior predictive checks confirm this: the spillover model produces
forward simulations consistent with the observed case count (median 82,
95% CrI [62–106] vs observed 103), while the full model under-predicts
(median 44) because the competing transmission pathways each receive
less weight.

### 7.3 Gravity model for movements

For counterfactual scenarios (Q4, Q5) that include farm-to-farm
transmission, movements during the prediction period are generated from
a gravity model fitted to observed movement data.

Observed movements are exclusively broiler_1 (source) → broiler_2
(destination), representing the day-old chick supply chain. The gravity
model is:

    P(i → j) ∝ cap_i · cap_j · exp(−d_ij / κ)

where cap is farm capacity and d_ij is Euclidean distance. The kernel
scale κ is estimated by MLE with a 2D pair density correction:
κ = mean(d) / 2. The daily movement count is Poisson-distributed at the
observed rate.

### 7.4 Prediction (Q2)

For each posterior draw, run the spillover model generatively: step
through time day by day from the end of the observation period, compute
hazards, draw new infections from Bernoulli(1 − exp(−λ)), and propagate
newly infected farms through the confirmation and culling pipeline.

For each new infection on day t:
- Confirmation at t + 11 (compound delay)
- Removal at t + 11 + 3 (confirmation + removal buffer)
- Surveillance zone established at confirmation (10 km, 28 days)
- Preventive culling triggered at confirmation (1 km)

500 forward trajectories are drawn per scenario.

### 7.5 Counterfactual scenarios (Q4, Q5)

Each counterfactual is run under two transmission assumptions:

**Scenario A — spillover only** (same as Q2 baseline): captures the effect
of culling policy on removing susceptible farms from the spillover-exposed
population.

**Scenario B — spillover + transmission**: adds spatial and movement
transmission using fixed parameters (β = 0.006, α = 4572 m, p_mov = 0.076)
from the full model's posterior medians, with movements generated by the
gravity model. This provides a sensitivity analysis for how culling policy
interacts with farm-to-farm spread.

The two scenarios bracket the range of possible outcomes:
- **Q4** (chicken-only preventive culling): restrict preventive culling
  to chicken farms only; duck farms within 1 km of cases are not culled.
- **Q5** (faster reactive culling): remove preventive culling entirely,
  reduce removal buffer by 1 day (from 3 to 2).

---

## 8. Data inputs

| File | Use |
|------|-----|
| population.csv | Farm locations (UTM), species, production type, capacity |
| activity.csv | Which farms have birds on each day |
| movement.csv | Movement schedule between farms |
| cases.csv | Observed cases: suspicion dates, confirmation dates, detection method, cull start |
| prev_culls.csv | Preventive culling records (farm, cull start date) |
| hrz_32626.geojson | High-risk zone polygon (UTM EPSG:32626) for spillover and testing classification |

### Data summary (Phase 1)

- 9,160 farms (4,435 chicken, 4,725 duck); 1,962 in HRZ
- 103 confirmed cases over 24 days (22 Dec 2025 – 13 Jan 2026)
- 12 preventive culls with recorded cull dates
- 81 chicken cases (1.83% attack rate), 22 duck cases (0.47% attack rate)
