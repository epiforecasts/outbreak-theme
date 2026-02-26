# Observation DAG

This document describes the observation processes that map latent states (from the process DAG) to observed data.

The observation DAG separates **what actually happens** (process) from **what we see** (observation). This distinction is critical for inference: we fit the model to observations, but the likelihood depends on how observations relate to latent states.

---

## Overview

| Latent quantity | Observed quantity | Observation process |
|-----------------|-------------------|---------------------|
| Infection time $T_j^I$ | Confirmation time $T_j^C$ | Detection delay |
| Infection time $T_j^I$ | Suspicion time $T_j^S$ | Farmer reporting |
| Infectious state $I_j(t)$ | Case record | Surveillance + confirmation |
| Removal time $T_j^R$ | Cull dates | Administrative recording |
| Transmission source | — | Not observed (spillover vs which farm) |
| Flock size at infection | Registered capacity | Proxy (mortality creates discrepancy) |
| Non-infection | Absence of case record | Survival likelihood |

---

## Detection Process

### From infection to confirmation

A farm progresses through stages:

```
Infection (T^I) → Detectable (T^D) → Suspected (T^S) → Confirmed (T^C) → Culled (T^R)
```

**Delays**:
- $\delta_{\text{amplification}} = T^D - T^I$: Incubation/amplification — time until detectable signs. Includes within-farm epidemic build-up (depends on within-farm growth rate $r$; see `step02_process_dag.md`)
- $\delta_{\text{recognition}} = T^S - T^D$: Recognition — time until farmer reports suspicion
- $\delta_{\text{confirmation}} = T^C - T^S$: Confirmation — time from suspicion to lab confirmation

**Observed**: $T^S$ (99/103 cases), $T^C$ (all cases)
**Latent**: $T^I$, $T^D$

### Detection methods

From `cases.csv`, detection occurs via:

| Method | Count | Description |
|--------|-------|-------------|
| passive | 96 | Farmer reports clinical signs |
| pre_shipment | 6 | Testing before broiler_1 → broiler_2 movement (HRZ only) |
| contact_tracing | 1 | Backward/forward tracing from confirmed case |

**Implications for observation model**:
- Passive detection: delay depends on clinical signs (species-dependent)
- Pre-shipment: may detect earlier (before clinical signs); only applies to HRZ farms with scheduled movements
- Contact tracing: detection triggered by epidemiological link, not clinical signs

### Species differences in detection

Chickens and ducks may differ in both detection probability and delay:

**Detection probability** $p_{\text{detect}}$:
- **Chickens**: High mortality → nearly all infections detected
- **Ducks**: Lower mortality, subclinical shedding → some infections may never be detected

**Detection delay** $\delta$:
- **Chickens**: Rapid onset of clinical signs → short delay
- **Ducks**: Subtler signs, longer time to recognition → longer delay

These are distinct effects:
- Lower $p_{\text{detect}}$ → under-ascertainment (missing cases in data)
- Longer $\delta$ → bias in inferred infection times (but case still observed)

Both contribute to why $\beta_{\text{duck}}$ is not identifiable from case data alone.

---

## Observation Model Components

### 1. Infection-to-confirmation delay

**Model**: $T^C = T^I + \delta$

where $\delta$ is the total delay from infection to confirmation.

**Options for $\delta$**:
1. **Fixed**: Set $\delta$ = median($T^C - T^S$) + assumed incubation
2. **Estimated**: Infer distribution of $\delta$ from data
3. **Species-specific**: $\delta_{\text{chicken}} \neq \delta_{\text{duck}}$

**Data available**:
- $T^C - T^S$ observed for 99/103 cases (median ~3 days based on earlier analysis)
- $T^I$ not observed — must infer

**Recommendation**: Model $\delta$ as sum of:
- $\delta_{\text{amplification}}$: time from farm infection to detectable mortality/morbidity. Includes within-farm epidemic build-up to detection threshold. It depends on flock size, mortality threshold, and within-farm growth rate $r$ (fixed at 1.0/day from mortality ledgers; see `step02_process_dag.md`).
- $\delta_{\text{recognition}} + \delta_{\text{confirmation}}$: estimable from $T^C - T^S$ distribution (observed for 99/103 cases)

### 2. Case ascertainment

Not all infections may be detected:
- Farms culled preventively before detection
- Farms with subclinical infection (especially ducks)
- Farms outside surveillance zones with mild presentation

**Model options**:
1. **Assume complete ascertainment**: All infected farms eventually confirmed (simplest)
2. **Species-specific ascertainment**: $p_{\text{detect,chicken}} = 1$, $p_{\text{detect,duck}} < 1$
3. **Scenario analysis**: Run model with different assumed duck ascertainment rates

**Recommendation**: Start with complete ascertainment (option 1). For Q3 interpretation, use scenario analysis (option 3) to explore how different duck ascertainment rates affect β_duck estimates. This fits the bracketing approach:
- Scenario A: Complete ascertainment ($p_{\text{detect,duck}} = 1$)
- Scenario B-D: Partial duck ascertainment ($p_{\text{detect,duck}} \in \{0.7, 0.5, 0.3\}$)

The range should be explored as sensitivity analysis rather than point estimates, since we lack external data to calibrate duck ascertainment.

Option 2 (estimating $p_{\text{detect,duck}}$) is not identifiable without external data on true infection rates.

**Identifiability note**: Under partial ascertainment, β_duck and $p_{\text{detect,duck}}$ trade off — higher ascertainment implies lower susceptibility and vice versa. Only the product $\beta_{\text{duck}} \times p_{\text{detect,duck}}$ is identified from case data.

### 3. Removal observation

**Reactive culls** (from `cases.csv`):
- `cull_start`, `cull_end` observed for confirmed cases
- Delay $\delta_{\text{reactive}} = T^R - T^C \geq 0$ from confirmation to cull start is observable (median ~2 days)

**Preventive culls** (from `prev_culls.csv`):
- 77% (40/52) missing cull dates
- Must impute or model

**Imputation approach** (see `step02_process_dag.md` §Preventive Culling and `step03_data_source_selection.md` §5):
For farms with missing preventive cull dates:
1. Identify trigger case (confirmed farm within 1 km that triggered preventive cull)
2. Estimate $\delta_{\text{prev}} \geq 0$ from the 12 complete records (trigger-to-cull delay)
3. Impute: $T^R_{\text{prev}} = T_{\text{trigger}} + \delta_{\text{prev}}$
4. **Fallback**: if no confirmed farm within 1 km, use the nearest confirmed case regardless of distance; if multiple equidistant, use the earliest confirmation date

### 4. Transmission source (latent)

Who infected whom is not observed. For each case, the source could be:
- **Spillover**: Wild bird introduction (if farm in HRZ)
- **Local transmission**: Infection from another farm

We observe the spatial and temporal pattern of cases, but not the transmission tree. The likelihood marginalises over possible sources:

$$P(T_j^I | \text{process}) = P(\text{spillover}) + \sum_{i \neq j} P(\text{infected by } i)$$

This is handled implicitly by the force of infection formulation — we don't need to identify specific transmission pairs.

**Note**: Genomic data (if available) could inform transmission pairs, but we don't have this.

### 5. Flock size observation

Flock size ($N_j$) from `population.csv` is the registered capacity, not necessarily the size at infection or detection. Mortality between infection and detection (especially for chickens with high HPAI mortality) means:
- True size at infection may differ from registered capacity
- Size at detection is reduced by deaths

**Current treatment**: Use registered capacity as proxy for flock size. This is standard practice but introduces error, particularly for high-mortality farms.

### 6. Non-case observation

Farms that were never confirmed contribute to the likelihood via survival:

$$P(\text{farm } j \text{ not infected by } T) = \exp\left(-\int_0^T \lambda_j(t) dt\right)$$

This requires knowing which farms were at risk (active, not yet culled) at each time point.

**Note**: We do not have data on farms that were tested but negative. If active surveillance tested farms in zones around cases, this would provide additional information constraining undetected infections. Currently we assume no such data are available.

---

## Likelihood Structure

The likelihood connects observations to the latent process:

$$L(\theta | \text{data}) = \prod_j L_j(\theta | T_j^C, T_j^R, \ldots)$$

For each observed case $j$:

$$L_j \propto P(T_j^C | T_j^I, \delta) \times P(T_j^I | \lambda_j(t), \text{process})$$

where:
- First term: observation model (delay distribution)
- Second term: process model (force of infection)

For non-cases (farms that remained susceptible):

$$L_j \propto P(\text{no infection by end}) = \exp\left(-\int_0^T \lambda_j(t) dt\right)$$

---

## DAG Representation

```
LATENT (Process DAG)              OBSERVED (Data)
─────────────────────             ──────────────────

[Spillover φ·ψ(t)]
     |
     v
[Force of infection λ_j(t)]
     |
     v
[Infection time T^I_j]  ──────>  {Suspicion time T^S_j}
     |                                    |
     | (delay δ)                          |
     v                                    v
[Detectable state]      ──────>  {Confirmation time T^C_j}
     |                                    |
     v                                    v
[Removal time T^R_j]    ──────>  {Cull dates}


OBSERVATION PROCESSES:
─────────────────────
(1) Detection delay: T^I → T^C via δ = δ_amplification + δ_recognition + δ_confirmation
(2) Case ascertainment: infected farms → confirmed cases
(3) Removal recording: true cull time → recorded dates
```

**Legend**:
- `[Square brackets]`: latent states/parameters
- `{Curly braces}`: observed data
- `→`: deterministic or stochastic mapping

---

## Key Assumptions

1. **All infections eventually confirmed** (no under-ascertainment) — relaxed via scenario analysis for ducks
2. **Detection delay is additive**: $T^C = T^I + \delta$
3. **Delay distribution is stationary** (doesn't change over outbreak) — may not hold; detection likely faster in later phases as surveillance intensifies
4. **Cull dates accurate for reactive culls** (recorded correctly)
5. **Preventive cull dates imputable** from trigger case timing
6. **Farm locations observed without error** — geocoding uncertainty not modelled
7. **Species classification accurate** — no mixed-species farms or misclassification

---

## Data Quality Issues

| Issue | Impact | Mitigation |
|-------|--------|------------|
| 4/103 missing $T^S$ | Minor | Use $T^C$ only for these cases |
| 77% missing preventive cull dates | Moderate | Impute from trigger + delay |
| Detection method variation | Moderate | Consider method-specific delays |
| Species-specific detection | High | Scenario analysis for ascertainment |
| Flock size vs actual size | Moderate | Use registered capacity as proxy |
| No negative test data | Minor | Assume no active surveillance data |
| Stationarity of delays | Unknown | Sensitivity analysis; may be faster late in outbreak |

---

## Link to Inference

The observation DAG determines:

1. **Likelihood function**: How to compute $P(\text{data} | \theta)$
2. **Latent variable augmentation**: Which latent states to sample (e.g., $T^I$)
3. **Missing data handling**: How to treat missing cull dates

**Inference approach options**:
- **Data augmentation MCMC**: Sample latent $T^I$ for each case
- **Marginalisation**: Integrate out $T^I$ analytically if delay distribution is tractable
- **Approximate**: Use $T^S$ or $T^C$ as proxy for $T^I$ with fixed offset

---

## Parameters in Observation Model

| Parameter | Type | Notes |
|-----------|------|-------|
| $\delta_{\text{amplification}}$ | Fix or estimate | Time from infection to detectable signs (1–3 days); depends on within-farm growth rate $r$ (fixed at 1.0/day in process model) |
| $\delta_{\text{recognition}}$ | Estimate | Time from detectable signs to farmer reporting suspicion |
| $\delta_{\text{confirmation}}$ | Estimate | Time from suspicion to lab confirmation (~2 days from $T^C - T^S$ data) |
| $\sigma_\delta$ | Estimate | Variance in total detection delay $\delta$ |
| $p_{\text{detect}}$ | Fix = 1 | Ascertainment probability (assume complete); relaxed via scenario analysis for ducks |
| $\delta_{\text{reactive}}$ | Impute from data | Confirmation-to-removal delay ($\geq 0$; median ~2 days); see `step02_process_dag.md` |
| $\delta_{\text{prev}}$ | Impute from data | Trigger-to-removal delay for preventive culling ($\geq 0$; from 12 complete records); see `step02_process_dag.md` |

**Note**: Process model parameters ($\phi$, $\eta$, $t_{\text{change}}$, $\rho$, $\beta$, $\alpha$, $\beta_{\text{duck}}$, and fixed parameters $\tau_{\min}$, $r$, $p_{\text{mov}}$, $\sigma_{\text{test}}$) are defined in `step02_process_dag.md`.

---

## Link to Research Questions

| Question | Observation model requirements |
|----------|-------------------------------|
| Q1 | Minimal — descriptive uses observed dates directly |
| Q2 | Full observation model for forecasting |
| Q3 | Species-specific delays inform β_duck interpretation |
| Q4 | Observation model unchanged; process model modified |
| Q5 | Observation model unchanged; process model modified |

---
