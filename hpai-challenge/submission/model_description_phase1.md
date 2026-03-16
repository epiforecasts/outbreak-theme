# HPAI modelling challenge model description — first phase

- **Team / institution**: LSHTM
- **Authors**: Sebastian Funk
- **Date**: 2026-01-14
- **Model name**: Jolly Island HPAI Spatial Transmission Model
- **Version**: 1.0
- **Contact**: sebastian.funk@lshtm.ac.uk
- **Repository / code link**: https://github.com/epiforecasts/outbreak-theme

---

## Model description

### Model overview

- **Type**: Mechanistic (Bayesian farm-level spatial force-of-infection model)
- **Stochastic vs deterministic**: Stochastic (Bernoulli infection events per farm per day; Bayesian posterior with MCMC)
- **Timestep**: Discrete, 1 day
- **Entities and levels considered**: Farm (single epidemiological unit per farm)
- **Spatial scale**: Metapopulation of 9,160 individual farms with explicit coordinates (UTM EPSG:32626)
- **Interaction mechanisms**: Wild-bird spillover (environmental), spatial proximity (exponential kernel), livestock movements (gravity model for prediction period)
- **Implementation**: Julia 1.11, Turing.jl (Bayesian inference), custom forward simulation engine

#### Narrative description

Each farm is a single epidemiological unit in one of three states: susceptible, infected, or removed (culled). Two main transmission pathways are modelled:

1. **Wild-bird spillover**: A time-varying environmental pressure, modelled as a Bateman function (rise then decay), with separate hazard rates for farms inside and outside the high-risk zone (HRZ). This captures the introduction of HPAI from wild bird populations.

2. **Local spatial transmission**: An exponential distance kernel β·exp(−d/α) weighted by farm infectiousness since infection. Only farms within 50 km contribute (computational cutoff).

Species-specific susceptibility is captured via a duck susceptibility parameter β_duck relative to chickens (β_chicken = 1.0). Disease management measures (reactive culling, surveillance zones, preventive culling, pre-shipment testing) are modelled mechanistically based on the narrative.

Inference uses a survival likelihood (force-of-infection model) fitted to the 103 observed cases. The spillover-only model (6 parameters), which fits the observed data adequately, is used for prediction; a full model adding spatial and movement transmission is used for sensitivity analyses in counterfactual scenarios. Predictions from the spillover model should be interpreted as a lower bound on future case counts if local transmission is occurring.

### State variables

| Name | Symbol | Description |
|------|--------|-------------|
| Susceptible | S | Active farms with birds that have not been infected or preventively culled |
| Infected | I | Farms that have been infected; infectiousness increases over time following w(τ) = 1 − exp(−r·(τ − τ_min)) |
| Removed | R | Farms that have been culled (reactive or preventive) and no longer contribute to transmission |

---

## Processes and parameter values

### Demographic processes

**Entry/exit processes**: No births or natural mortality are modelled. Farms exit the susceptible pool through infection (→ Infected) or preventive culling (→ Removed). Infected farms exit through reactive culling.

**Production cycles / occupancy / downtime**: Farm activity status (whether birds are present) is read from the activity dataset and varies by day. Only active farms can be infected.

**Movements or exchanges**: Livestock movements are read from the movement dataset during the observation period. For the prediction period, movements are generated from a fitted gravity model: P(i→j) ∝ cap_i · cap_j · exp(−d_ij/κ), with daily movement count drawn from a Poisson distribution at the observed rate. Movements from farms inside regulated zones are blocked.

### Epidemiological processes

**Health states**: Susceptible → Infected → Removed (SIR at the farm level). Once infected, a farm's infectiousness follows a saturating profile w(τ) = 1 − exp(−r·(τ − τ_min)) for τ ≥ τ_min, and w(τ) = 0 otherwise.

**Transmission**: The daily hazard for susceptible farm j on day t is:

    λ_j(t) = β_sj · [spillover_j(t) + spatial_j(t) + movement_j(t)] · zone_j(t)

where:
- β_sj = 1.0 (chicken) or β_duck (duck)
- spillover_j(t) = φ_zone · ψ(t), with ψ(t) a Bateman function
- spatial_j(t) = Σ_i β·exp(−d_ij/α)·w(t − t_inf_i) over infectious farms i
- movement_j(t) = p_mov · w(t − t_inf_source) for each movement from an infectious source
- zone_j(t) = (1 − ε) if inside a regulated zone, 1.0 otherwise

Infection probability: P(infected) = 1 − exp(−λ_j(t))

**Progression and recovery**: Fixed timeline per case:
- Latent period (μ_E = 3.5 days): infection → infectious
- Detection delay (μ_ID = 5.0 days): infectious → suspicion raised
- D→C delay (2 days): suspicion → confirmation
- C→R delay: from data per case (median 3 days)
- Compound delay: 11 days from infection to confirmation

**Detection and surveillance**: Not explicitly modelled as a stochastic process. Detection times are determined by the fixed delay pipeline. Surveillance zones enhance biosecurity (reducing hazard by factor 1 − ε = 0.5) and block movements.

**Disease-induced mortality / production impact**: Not modelled explicitly. Within-farm dynamics are captured implicitly through the infectiousness profile w(τ).

### Disease management (control measures)

| Management measure | Included? | Implementation description |
|---|---|---|
| National standstill | N | Not modelled |
| Suspicion management | N | Detection delays are fixed parameters |
| Reactive culling | Y | Confirmed farms culled after observed C→R delay (from data). Predicted cases: 3-day buffer after confirmation. Culling removes farm from infectious pool. |
| Contact tracing | N | Not modelled |
| Zoning: Movement bans in regulated zones | Y | 10 km surveillance zones around confirmed farms for 28 days. Movements from/to zone farms are blocked. |
| Zoning: Enhanced biosecurity in regulated zones | Y | Farms inside zones have hazard multiplied by (1 − ε) = 0.5 |
| Preventive culling | Y | From 1 Jan 2026: active susceptible farms within 1 km of confirmed cases are removed |
| Pre-movement testing in the HRZ | Y | Infectious source farms in HRZ detected with sensitivity σ_test = 0.9, blocking the movement |

---

## Parameter table

### Estimated parameters (spillover model)

| Symbol | Description | Unit | Default value | Range / distribution | Source |
|--------|-------------|------|---------------|---------------------|--------|
| t₀ | Spillover onset day | day | — | truncated Normal(15, 5) on [1, 44] | Estimated |
| φ_hrz | Spillover rate in HRZ | /farm/day | 0.0017 | LogNormal(log(10⁻³), 1.0); posterior median 0.0017 [0.0012–0.0024] | Estimated |
| φ_non | Spillover rate outside HRZ | /farm/day | 0.0002 | LogNormal(log(10⁻⁴), 1.0); posterior median 0.0002 [0.0001–0.0004] | Estimated |
| δ | Spillover post-peak decay rate | /day | 0.011 | Exponential(mean=0.02); posterior median 0.011 [0.0001–0.041] | Estimated |
| β_duck | Duck susceptibility (relative to chicken) | — | 0.24 | Beta(2, 8); posterior median 0.24 [0.13–0.37] | Estimated |
| σ | Spillover rise rate | /day | — | LogNormal(log(0.3), 1.0) | Estimated |

### Fixed parameters (spatial/movement transmission, used in sensitivity analyses)

| Symbol | Description | Unit | Value | Source |
|--------|-------------|------|-------|--------|
| β | Spatial transmission rate | — | 0.006 | Full model posterior median |
| α | Spatial kernel range | m | 4,572 | Full model posterior median |
| p_mov | Per-movement infection probability | — | 0.076 | Full model posterior median |

### Fixed parameters (epidemiology and interventions)

| Symbol | Description | Unit | Value | Source |
|--------|-------------|------|-------|--------|
| μ_E | Latent period | days | 3.5 | HPAI literature |
| μ_ID | Detection delay | days | 5.0 | HPAI literature |
| β_chicken | Chicken susceptibility | — | 1.0 | Reference category |
| r | Infectiousness growth rate | /day | 1.0 | Mortality ledgers (3 farms) |
| τ_min | Hard latent period | day | 1 | Minimum before infectiousness |
| D→C delay | Suspicion to confirmation | days | 2 | Data (median) |
| C→R delay (observed) | Confirmation to culling | days | Per case | Data |
| C→R delay (predicted) | Confirmation to culling | days | 3 | Fixed buffer |
| Zone radius | Surveillance zone radius | km | 10 | Narrative |
| Zone duration | Surveillance zone duration | days | 28 | Narrative |
| ε | Zone biosecurity effect | — | 0.5 | Assumption |
| Preventive cull radius | — | km | 1 | Narrative |
| Preventive cull start | — | date | 1 Jan 2026 | Narrative |
| σ_test | Pre-shipment test sensitivity | — | 0.9 | Assumption |
| Neighbour cutoff | Spatial kernel cutoff | km | 50 | Computational |

---

## Parameter estimation

The spillover model (6 parameters) is fitted using Bayesian inference with MCMC. The likelihood is a force-of-infection survival model evaluated over the observation period (1 Dec 2025 – 13 Jan 2026, T = 44 days):

    ℓ(θ) = Σ_t [ Σ_{j ∈ infected(t)} log(1 − exp(−λ_j(t))) + Σ_{j ∈ survived(t)} (−λ_j(t)) ]

This is a standard spatial point-process likelihood — fully analytical, no simulation required.

Computational optimisation exploits sparsity: bulk spillover survival for ~9,000 non-case farms is aggregated into 8 category counts per day (species × HRZ × zone status); only ~115 farms need individual accounting. This gives ~1.4 ms per likelihood evaluation.

Sampling uses NUTS (No-U-Turn Sampler) with ForwardDiff automatic differentiation (chunk size = 6). A coarse-then-fine grid search over parameter space identifies an approximate MAP estimate as the NUTS starting point.

Both spillover-only and full (9-parameter) models were fitted. The full model's spatial and movement parameters are not independently identifiable from spillover in this dataset. Posterior predictive checks confirm the spillover model produces forward simulations consistent with observed cases (median 82, 95% CrI [62–106] vs observed 103).

---

## Initial conditions

- **Simulation start date**: 1 Dec 2025 (day 1)
- **Initial infection seeding**: All 103 observed cases are seeded with back-calculated infection times derived from their suspicion/confirmation dates and the fixed delay parameters. No additional index cases are assumed.

---

## Simulations and outputs

- **Number of runs**: 500 forward trajectories per scenario, each drawn from a different posterior sample
- **Random seed handling**: Each trajectory uses a unique random seed
- **Outputs**: Daily new confirmed cases (total, by species), cumulative cases, spatial distribution by farm and by district
- **Output format**: CSV files

---

## First period

### Q1. General description of the ongoing epidemic

The HPAI epidemic on Jolly Island comprises 103 confirmed cases over a 23-day period (22 Dec 2025 – 13 Jan 2026). The epidemic is strongly concentrated in chicken farms and in the southern high-risk zone.

**Distribution by species and production type:**

| Species | Production type | Cases | Active farms | Attack rate (%) |
|---------|----------------|-------|--------------|----------------|
| Chicken | broiler_1 | 18 | 517 | 3.5 |
| Chicken | broiler_2 | 56 | 1,411 | 4.0 |
| Chicken | layer | 7 | 415 | 1.7 |
| Duck | conventional | 11 | 2,684 | 0.4 |
| Duck | organic | 11 | 822 | 1.3 |
| **Total** | | **103** | **5,849** | **1.8** |

Chickens account for 81/103 cases (79%), with a chicken attack rate of 3.5% (81/2,343 active farms) vs 0.6% for ducks (22/3,506 active farms). Broiler farms (types 1 and 2) are disproportionately affected.

**Geographic distribution**: Cases are concentrated in the eastern high-risk zone. Berks county has the most cases (44, attack rate 3.8%), followed by Susquehanna (22, 3.9%) and Indiana (20, 3.9%). Additional cases appear in Allegheny (5), Luzerne (4), Lancaster (3), Lycoming (3), and Cumberland (2).

**Temporal dynamics**: The epidemic curve shows an initial cluster of cases from 22 Dec, building to a peak of approximately 8–10 new confirmed cases per day in the first week of January 2026, with a slight decline towards the end of the observation period. This temporal pattern is consistent with the Bateman spillover function peaking and beginning to decay.

See `epicurve_confirmed.csv` for the full daily epicurve, `cases_by_county.csv` for geographic distribution, and `q1_spatial_map.pdf` for a visual representation of the spatial distribution.

### Q2. Prediction of temporal and spatial evolution

Using the spillover-only model with 500 forward simulations from the posterior, we predict:

- **Cumulative new confirmed cases over 4 weeks (14 Jan – 10 Feb 2026)**: median 45, 95% CrI [22–79]
- **Daily new confirmed cases**: median approximately 2 per day, declining slightly over the period
- **Species breakdown**: approximately 78% chicken, 22% duck (consistent with observed proportions)

The epidemic is predicted to continue at a declining rate as the wild-bird spillover pressure wanes (captured by the decay parameter δ). The spatial distribution of new cases is expected to remain concentrated in the high-risk zone, particularly in Berks and surrounding counties.

See `q2_temporal_summary.csv` for daily summary statistics, `q2_trajectories.csv` for individual trajectories, `q2_spatial_farm.csv` and `q2_spatial_district.csv` for spatial predictions.

### Q3. Relative contribution of chicken vs duck farms

The posterior distribution of β_duck (duck susceptibility relative to chicken) has:
- **Median**: 0.24 (95% CrI: 0.13–0.37)
- **Implied chicken relative risk**: median 4.2× (95% CrI: 2.7–7.5×)

This means chickens are approximately 4 times more susceptible to HPAI infection than ducks. The susceptibility-weighted farm share is computed as chicken_share = n_chicken / (n_chicken + β_duck × n_duck), where n is the number of active farms:
- **Chicken share**: ~80% of total susceptibility-weighted exposure
- **Duck share**: ~20% of total susceptibility-weighted exposure

Despite ducks comprising ~60% of active farms, their lower susceptibility means they contribute only about one-fifth of the effective susceptible population. Chickens bear a disproportionate burden of infection.

See `q3_species_contribution.csv` for the full posterior distribution of species contributions.

### Q4. Chicken-only preventive culling

We evaluated the impact of restricting preventive culling to chicken farms only (i.e., not culling duck farms within 1 km of confirmed cases) under two transmission assumptions:

**Spillover-only scenario**: Cumulative new cases over 4 weeks: median 45 [21–79], essentially identical to the baseline (45 [22–79]). Under spillover-dominated dynamics, removing susceptible duck farms from preventive culling has negligible impact because (a) duck farms have low susceptibility (β_duck ≈ 0.24) and (b) spillover acts independently on each farm regardless of neighbours' status.

**Spillover + transmission scenario**: Cumulative new cases: median 103 [54–178], compared to a baseline with full preventive culling. The addition of spatial transmission means that un-culled duck farms near outbreaks could serve as stepping stones for onward spread, though the effect is modest given ducks' lower susceptibility.

**Conclusion**: Under the spillover-dominated model, which fits the observed data adequately, restricting preventive culling to chickens only would have minimal epidemiological impact. If farm-to-farm transmission plays a larger role than estimated, there could be a moderate increase in cases, but this depends on the spatial clustering of duck farms near outbreak sites.

See `q4_spillover_temporal_summary.csv` and `q4_transmission_temporal_summary.csv` for full results.

### Q5. Faster reactive culling vs preventive culling

We evaluated replacing preventive culling with faster reactive culling (reducing the confirmation-to-removal delay by 1 day, from 3 to 2 days) under two transmission assumptions:

**Spillover-only scenario**: Cumulative new cases: median 44 [21–81], essentially identical to the baseline (45 [22–79]). Under spillover-dominated dynamics, neither preventive culling nor reactive culling timing substantially affects the epidemic because each new case arises independently from wild-bird pressure, not from farm-to-farm transmission.

**Spillover + transmission scenario**: Cumulative new cases: median 105 [55–170]. This is comparable to the chicken-only preventive culling scenario (103 [54–178]), suggesting that the benefit of removing the most infectious day (when w(τ) is highest, just before culling) is approximately offset by the loss of preventive culling's protection of nearby susceptible farms.

**Conclusion**: Under the spillover-dominated model, switching to faster reactive culling without preventive culling makes negligible difference. Under a transmission scenario, the trade-off between faster removal of infected farms and loss of preventive protection is approximately neutral. The current combined strategy (reactive + preventive culling) is therefore a reasonable approach.

See `q5_spillover_temporal_summary.csv` and `q5_transmission_temporal_summary.csv` for full results.

---

## Strengths and limitations

### Strengths

- **Mechanistic framework**: Explicit modelling of transmission pathways allows scenario comparison and counterfactual analysis
- **Bayesian inference**: Full posterior uncertainty propagated through predictions
- **Efficient likelihood**: Analytical FOI likelihood (~1.4 ms/evaluation) enables thorough MCMC exploration
- **Species-specific susceptibility**: Estimated from data rather than assumed
- **Comprehensive intervention modelling**: All major control measures from the narrative are represented

### Limitations

- **Identifiability**: Spatial and movement transmission parameters are not independently identifiable from spillover in this dataset. The full model is used only for sensitivity analysis.
- **Farm-level granularity**: Within-farm dynamics (bird-level transmission, partial depopulation) are not modelled; each farm is a single unit.
- **Fixed delay parameters**: Latent period and detection delay are fixed at literature values because they cannot be estimated from the data (zero gradient through integer rounding).
- **Spillover model**: Wild-bird pressure is modelled as a simple parametric function (Bateman) rather than driven by explicit wild-bird surveillance data.
- **No between-species transmission differences**: The model estimates susceptibility differences but assumes the same infectiousness profile for both species once infected.

---

## Effort estimate

Approximately 8–10 person-days over the challenge period, including model development, data processing, inference, prediction, and scenario analysis.

---

## References

1. Boender GJ, et al. (2007). Risk maps for the spread of highly pathogenic avian influenza in poultry. PLoS Computational Biology, 3(4), e71.
2. Tildesley MJ, et al. (2006). Optimal reactive vaccination strategies for a foot-and-mouth outbreak in the UK. Nature, 440(7080), 83–86.
3. Keeling MJ, et al. (2001). Dynamics of the 2001 UK foot and mouth epidemic. Science, 294(5543), 813–817.
