# Research Questions and Target Estimands

This document defines the research questions from the HPAI Modelling Challenge Phase 1 and translates them into target estimands.

---

## Overview

The challenge poses 5 questions (page 3 of narrative):

| Question | Type | Priority |
|----------|------|----------|
| Q1 | Descriptive | High |
| Q2 | Prediction | High |
| Q3 | Mechanistic inference | High |
| Q4 | Counterfactual | Bonus |
| Q5 | Counterfactual | Bonus |

**Deadline**: 28 February 2026

---

## Q1: Epidemic Description (Descriptive)

### Research Question (verbatim from narrative)

> "A general description of the ongoing epidemic, including a table describing the distribution of the outbreaks per species and production types, a timeline of the incidence and a visual representation of the spatial distribution of the outbreaks."

### Target Estimands

This is purely descriptive — quantities to compute from data:

1. **Distribution table**: Case counts by species (chicken/duck) × production type (layer, broiler_1, broiler_2, organic, conventional)

2. **Timeline**: Epidemic curve showing incidence over time (22 Dec 2025 – 13 Jan 2026)

3. **Spatial distribution**: Map of outbreak locations

### Data Requirements

- `cases.csv`: farm_id, date_confirmed
- `population.csv`: species, production_type, coordinates

### Notes

- Describes confirmed cases only
- No model required

---

## Q2: Epidemic Prediction (Prediction)

### Research Question (verbatim)

> "A prediction of the likely temporal and spatial evolution of the epidemic over the following four weeks, assuming the management strategy remains as it is today; if you can have it by farm production type, that would be even better."

### Target Estimands

1. **Temporal prediction**: Expected number of new confirmed cases per day for 14 Jan – 10 Feb 2026 (28 days)
   - With uncertainty (credible intervals or distribution)
   - Stratified by production type if possible

2. **Spatial prediction**: Per-farm confirmation probability $P(\text{confirmed by } t_{\text{end}} \mid \text{data})$ for each susceptible farm, using farm coordinates from `population.csv`. For district-level reporting, aggregate by summing expected cases within each district polygon.

### Formal Definition

Let $C(t)$ = number of farms confirmed on day $t$.

**Primary estimand**: Posterior predictive distribution
$$P(C(t) \mid \text{data}), \quad t \in \{14\text{ Jan}, \ldots, 10\text{ Feb } 2026\}$$

**Secondary**: Cumulative cases over horizon, per-farm infection probabilities

### Required Outputs (from narrative p.3-4)

- Raw simulated trajectories as CSV, or
- Synthetic metrics with distributions/credible intervals
- Spatial risk as CSV or shapefile

### Assumptions

- Management strategy unchanged from current (reactive culling, preventive culling within 1km, pre-movement testing in HRZ)
- Predicting confirmations, not infections (infections are latent)

---

## Q3: Species Contribution (Mechanistic Inference)

### Research Question (verbatim)

> "A characterization of the relative contribution of the galliform farms to the virus spread as compared to the palmiped farms"

### Context from Narrative

The narrative states (p.2):
> "early raw epidemic data suggest that galliform (chicken) farms may have been more frequently infected as compared to palmiped (duck) farms, suggesting galliform farms could play a more important role in the spread"

### Target Estimands

The question asks about "contribution to spread" — this could mean:

1. **Susceptibility**: Are chickens more likely to become infected given exposure?
2. **Infectiousness**: Do infected chicken farms transmit more than infected duck farms?
3. **Overall contribution**: What fraction of transmission events involve chicken farms?

### Possible Formal Definitions

**Option A — Relative susceptibility**:
$$\beta_{\text{ratio}} = \frac{\beta_{\text{duck}}}{\beta_{\text{chicken}}}$$
where $\beta_s$ modifies force of infection for species $s$

**Option B — Attributable transmission**:
Fraction of infections caused by chicken farms vs duck farms (requires reconstruction of transmission tree)

**Option C — Descriptive**:
Attack rate ratio = (cases in chickens / chicken farms) ÷ (cases in ducks / duck farms)

**Causal vs associational**: Options A and B are associational estimates from observational data — they describe patterns but do not identify causal effects without additional assumptions (no unmeasured confounding, no selection bias). Option C is purely descriptive. We treat Q3 as associational throughout.

### Interpretation Challenges

From data characterisation (Gap 2): We cannot distinguish "chickens more susceptible" from "chickens more detectable". The observed 4× attack rate difference (chicken 1.83% vs duck 0.47%) is consistent with either explanation.

Any estimate of $\beta_{\text{duck}}$ conflates:
- True biological susceptibility
- Detection probability (passive surveillance triggers on mortality; ducks may have lower mortality)

### Proposed Approach

Given non-identifiability, we propose three complementary analyses:

1. **Estimate from data**: Fit model with species-specific susceptibility parameter. Report estimate with caveat that it reflects "apparent susceptibility" (true susceptibility × detection probability).

2. **Scenario A — All susceptibility**: Assume detection probability equal for both species. The 4× attack rate difference reflects true biological susceptibility (chickens 4× more susceptible).

3. **Scenario B — All detectability**: Assume true susceptibility equal. The 4× difference reflects detection bias (ducks 4× less likely to be detected when infected).

This bracketing approach:
- Shows sensitivity of conclusions to interpretation
- If Q4/Q5 counterfactuals differ substantially between scenarios → interpretation matters for policy
- If similar → ambiguity may not matter for practical decisions

**This limitation and the scenario analysis should be explicitly reported.**

---

## Q4: Chicken-Only Preventive Culling (Counterfactual)

### Research Question (verbatim)

> "Given we have reached the culling capacity and that chicken farms likely are key contributor to the epidemic, how epidemiologically-relevant would it be to focus the preventive culling actions on the chicken farms only?"

### Context from Narrative

- Preventive culling started 1 Jan 2026
- Targets all farms within 1km of confirmed cases (irrespective of species)
- Culling capacity reached since 6 Jan — delays occurring

### Target Estimand

Compare predicted outcomes under two scenarios:

**Baseline**: Current policy — preventive cull all farms within 1km

**Alternative**: Preventive cull only chicken farms within 1km

**Estimand**: Difference in cumulative cases over 14 Jan – 10 Feb 2026 (same horizon as Q2)
$$\Delta C = C_{\text{chicken-only}} - C_{\text{baseline}}$$

### Interpretation

- If $\Delta C < 0$: chicken-only policy performs *better* (fewer cases)
- If $\Delta C > 0$: chicken-only policy performs *worse*
- If $\Delta C \approx 0$: no meaningful difference

### Assumptions Required

- Same total culling capacity under both scenarios
- Model correctly captures species-specific transmission
- Forward simulation under modified intervention rules

---

## Q5: Faster Reactive vs Preventive Culling (Counterfactual)

### Research Question (verbatim)

> "Given that outbreak farms are likely more infectious than preventively-culled farms, how epidemiologically-relevant would it be to ignore the preventive culling actions and spend all efforts possible to reduce by 1 day the start of the reactive culling actions in outbreak farms?"

### Context from Narrative

- Reactive culling has priority over preventive culling
- Culling capacity reached — delays occurring
- Question assumes outbreak farms (confirmed infected, actively shedding) are more infectious than farms preventively culled (may not be infected)

### Target Estimand

Compare predicted outcomes under two scenarios:

**Baseline**: Current policy — reactive culling + preventive culling within 1km

**Alternative**: No preventive culling; confirmation-to-culling delay reduced by 1 day

**Estimand**: Difference in cumulative cases over 14 Jan – 10 Feb 2026 (same horizon as Q2)
$$\Delta C = C_{\text{faster-reactive}} - C_{\text{baseline}}$$

### Biological Rationale

The question implicitly assumes:
- Outbreak farms are more infectious than not-yet-infected farms
- Reducing infectious period of confirmed farms may have larger impact than removing potentially uninfected neighbours
- Resources spent on preventive culling could instead speed up reactive culling

### Assumptions Required

- 1-day reduction in reactive culling delay is achievable
- Model tracks culling delays
- Preventive culling resources can be reallocated to reactive culling

---

## Summary

| Question | Estimand | Type | Model Required |
|----------|----------|------|----------------|
| Q1 | Case counts, epi curve, map | Descriptive | No |
| Q2 | $P(C(t) \mid \text{data})$ for 4 weeks | Predictive | Yes (generative) |
| Q3 | Relative contribution of species | Inferential | Yes (but interpretation limited) |
| Q4 | $\Delta C$ under chicken-only culling | Counterfactual | Yes (simulation) |
| Q5 | $\Delta C$ under faster reactive | Counterfactual | Yes (simulation) |

---

## Model Requirements (for Q2-Q5)

To address the non-descriptive questions, a model must:

1. Capture spatial transmission (farm-to-farm spread)
2. Include wild bird spillover pathway
3. Allow species-specific parameters (for Q3)
4. Track intervention effects: zones, preventive culling, reactive culling
5. Be generative (forward simulation for Q2, Q4, Q5)
6. Quantify uncertainty

---

## Data Gaps Affecting These Questions

From `step00_data_source_characterisation.md`:

| Gap | Affects | Implication |
|-----|---------|-------------|
| Infection times latent | Q2, Q3, Q4, Q5 | Must infer infection times from confirmation dates |
| Species detection difference uncertain | Q3 | β estimates confounded with detectability |
| Spillover not observed | Q2, Q3 | Must infer from HRZ/spatial pattern |
| 77% preventive cull dates missing | Q4, Q5 | Must impute or assume removal timing |
| Movement data partial | Q2, Q3 | Movement transmission may be underestimated |

---

