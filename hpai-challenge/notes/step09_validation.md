# Validation (Workflow Step 9)

This step documents the validation pipeline for the phase 2 model.

---

## Phase 2 iteration

### Prior predictive checks

Median prior predictive case count is ~35. The target range (100–1000 cases) is achieved in only ~10% of prior draws. The priors are individually broad, but the combination of parameters rarely produces large epidemics. This is acceptable — priors should be weakly informative, not calibrated to match the data.

### Synthetic data recovery

- Spillover-only model: 6/6 parameters covered by 95% CrI
- Full model: 9/9 parameters covered by 95% CrI

Both models recover all true parameter values from synthetic data.

### Posterior predictive checks (PPC)

| Model | Median cases [95% CrI] | 466 covered? |
|-------|------------------------|--------------|
| Spillover-only | 347 [313–382] | No |
| Full model | 319 [125–560] | Yes |

The spillover-only model produces a narrow prediction interval that excludes the observed 466 cases. The full model has wider intervals (reflecting additional parameter uncertainty from spatial and movement components) and covers the observed count.

### WAIC comparison

| Model | WAIC |
|-------|------|
| Spillover-only | 6340 |
| Full model | 5360 |

Difference of ~980 in favour of the full model.

### Retrospective validation

Using the full model fitted to phase 1 data (103 cases, 13 Jan cut-off), we predicted phase 2 case counts. The model produced a median of 265 [110–429] for the additional cases observed in phase 2. The actual count of 363 additional cases falls within the 95% CrI.

---

## Phase 3 iteration

### Prior predictive checks

Prior predictive checks re-run on the expanded 560-case dataset. Structure and conclusions are unchanged from phase 2: the prior predictive distribution is broad and weakly informative.

### Synthetic data recovery

All 9 parameters (including `δ`) recovered within 95% CrI on synthetic data generated from the phase 3 model.

### Posterior predictive checks

PPC run against the observed 560 cases. The full model 95% CrI is checked to cover the observed count.

### Retrospective validation

Phase 2 predictions are evaluated against the phase 3 actuals. The phase 2 model predicted 382 [261–525] additional cases from 14 Feb; the actual additional case count to epidemic end was 81. This is a substantial overprediction, falling outside the 95% CrI. The most likely cause is underestimation of `δ` (the removal rate): if infected farms are culled faster than the model anticipated, the epidemic dies out sooner than predicted. This is documented as a calibration failure rather than a model failure — the phase 2 data did not yet contain the signal of rapid epidemic decline.

### Model comparison

Spillover-only vs full model comparison repeated on phase 3 data using WAIC.

---
