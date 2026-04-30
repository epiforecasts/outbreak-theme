# Predictions (Workflow Step 10)

This step documents the predictions and scenario analyses for the phase 2 submission.

---

## Phase 2 iteration

### Q2: Baseline prediction

4-week forward prediction from 14 Feb 2026, assuming current management continues.

**Result**: 382 [261–525] cumulative new cases.

### Q3: Drop preventive culling + faster reactive culling

No preventive culling of chicken farms; reactive culling delay reduced by 1 day.

**Result**: 451 [274–628] cumulative new cases, approximately +18% relative to baseline.

Dropping preventive culling increases cases more than the 1-day reduction in reactive delay can offset.

### Q4: Stop confinement

#### Q4a: Stop confinement for organic duck farms only

**Result**: 415 [286–590], approximately +9% relative to baseline.

Organic duck farms make up a small fraction of confined farms, so the effect is modest.

#### Q4b: Stop confinement for all confined farms

**Result**: 773 [550–1035], approximately +102% relative to baseline.

Removing confinement from all broiler\_2 and organic duck farms roughly doubles expected cases. This is the largest effect of any scenario tested.

### Q5: Prohibit restocking in surveillance zones

**Result**: 353 [237–490], approximately -8% relative to baseline.

A small reduction. Few restocked farms become infected during the prediction window — restocking is slow (14-day delay) and surveillance zones cover a limited area.

---

## Phase 3 iteration

### Q1: Full epidemic description

Descriptive summary of the complete epidemic: 560 confirmed cases, last case on 7 Apr. No forward simulation required — all data are observed. Summary statistics (epidemic curve, spatial distribution, farm-type breakdown) are computed directly from the phase 3 dataset.

### Q2: Restocking safety sweep

For each candidate restocking start date, 500 forward simulations are run from the posterior. The earliest date for which the probability of at least one rebound case (P(rebound)) falls below 5% is identified as the recommended safe restocking date.

### Q3: Counterfactual — preventive culling

Pairs of `simulate_epidemic` runs are drawn from the same posterior samples, one with preventive culling enabled (as observed) and one with `:no_prev_cull`. The difference in total cases between the two arms gives the distribution of outbreaks averted by the preventive culling policy.

### Q4: Capacity sweep

For each candidate daily cull capacity level, 500 `simulate_epidemic` runs are carried out to day 75. The minimum capacity for which 95% of runs achieve eradication (zero active farms) by day 75 is identified.

### Gravity model for predicted movements

Future movements (beyond the data cut-off) are needed for forward simulation. We fit a gravity model to observed broiler\_1 $\to$ broiler\_2 movements:

- Distance decay: $\kappa = 36.8$ km
- Rate: 48.5 movements/day

Predicted movements are sampled from this model during forward simulation.

### Summary

| Scenario | Median [95% CrI] | Change vs baseline |
|----------|-------------------|--------------------|
| Q2 baseline | 382 [261–525] | — |
| Q3 no prev cull + faster reactive | 451 [274–628] | +18% |
| Q4a stop organic duck confinement | 415 [286–590] | +9% |
| Q4b stop all confinement | 773 [550–1035] | +102% |
| Q5 prohibit restocking | 353 [237–490] | -8% |

---
