# Implementation (Workflow Step 8)

This step documents the software implementation of the model specified in steps 05-07.

---

## Phase 2 iteration

### Software stack

- Julia 1.12
- Turing.jl 0.43 for probabilistic programming
- ForwardDiff for automatic differentiation (forward-mode, chunk size = number of parameters)

### MCMC configuration

`MCMCSerial()` is used instead of `MCMCThreads()`. The threaded sampler causes `NumericalError` due to non-thread-safe state in the likelihood evaluation. Serial chains are run sequentially on the same core, or distributed across HPC jobs.

### HPC setup

- 4 cores, 32 GB RAM per job
- 24–48 hour wall time
- 4 chains run as separate jobs, combined post-hoc

### Forward simulation

Forward simulation starts from $T_{\text{LIK}}$ (day 64), not $T_{\text{SIM\_END}}$ (day 75). At day 64, approximately 73 farms are infectious. At day 75, only ~9 remain — most have been culled in the intervening 11 days. Starting from day 64 gives the forward simulation a realistic number of active sources and lets it generate unobserved recent infections through the censored gap.

---

## Phase 3 iteration

### config.jl

- `SIM_END` extended to 131 (covering the full epidemic to 7 Apr).
- Data paths updated to `phase-3/`.
- `CONFINEMENT_END_DAY_ORGANIC_DUCK` added to mark the date organic duck confinement was lifted.

### data.jl

- `ModelData` gains an `is_organic_duck` field to distinguish organic duck farms from other confined farm types.
- `bulk_counts` is now time-varying, updated to reflect restocking activity during the epidemic.

### likelihood.jl

- `is_confined_at(farm, day)` helper added, handling both the standard confinement period and the organic duck confinement end day.
- Applied across 8 confinement check sites in the likelihood.

### simulation.jl

- `confinement_end_day` for organic duck farms is threaded through the simulation state.
- `restocking_start_day` and `daily_cull_capacity` added as keyword arguments to enable scenario sweeps without recompiling.
- New `:no_prev_cull` scenario disables preventive culling for the counterfactual analysis.

### generative.jl

- `disable_prev_cull` and `daily_cull_capacity` forwarded to `simulate_forward` for counterfactual and capacity-sweep simulations.

### predict.jl

Restructured around four questions:

- **Q1**: Descriptive summary of the full epidemic — no forward simulation required.
- **Q2**: Restocking safety sweep — forward simulations over a grid of candidate start dates.
- **Q3**: Counterfactual — paired runs with and without preventive culling to estimate outbreaks averted.
- **Q4**: Capacity sweep — forward simulations over a grid of daily cull capacity values.

---
