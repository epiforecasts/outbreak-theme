# HPAI modelling challenge model description — third phase

**Team:** LSHTM
**Authors:** Sebastian Funk
**Date:** 2026-04-20

**Model name:** Spatial FOI model
**Version:** 3.0 (phase 3 iteration)
**Contact:** sebastian.funk@lshtm.ac.uk

## Model description

### Changes since phase 2

The core model structure is unchanged: a continuous-time force-of-infection (FOI) model fitted to farm-level outbreak data via Bayesian inference. Two modifications were made for phase 3:

1. **Organic duck confinement lifted.** From 14 February 2026, confinement of organic duck farms was lifted. We represent this as a time-limited confinement window (14 January – 13 February) for organic duck farms, while stage 2 broiler confinement continues indefinitely. The confinement factor γ = 0.5 is applied piecewise by farm type and time.

2. **Extended observation window.** The simulation window was extended from 75 days (13 February) to 131 days (10 April), covering the entire epidemic. The right-censoring correction (truncating the survival likelihood at T − 11 days) remains in place for structural consistency but is vacuous, since no infections in the final 11 days remain unconfirmed.

No changes were made to the model structure, transmission mechanisms, or inference approach.

### Model overview

- **Type:** Mechanistic, continuous-time FOI
- **Stochastic vs deterministic:** Deterministic likelihood; stochastic forward simulation
- **Timestep:** Discrete (1 day)
- **Entities:** Individual farms
- **Spatial scale:** Island-wide, farm-level coordinates (EPSG:32626)
- **Interaction mechanisms:** Environmental spillover, spatial proximity (exponential kernel), live animal movements
- **Implementation:** Julia 1.12, Turing.jl (NUTS sampler), ForwardDiff for automatic differentiation

### Model structure

Each susceptible farm *j* on day *t* experiences a total hazard rate:

$$\lambda_j(t) = s_j \cdot [\varphi_j \cdot \psi(t) + \beta \sum_{i \in \mathcal{I}(t)} w(\tau_{ij}) \cdot e^{-d_{ij}/\alpha} + \sum_{m \in \mathcal{M}_j(t)} p_{\mathrm{eff}} \cdot w(\tau_m)] \cdot z_j(t) \cdot \gamma_j(t)$$

where:
- $s_j$ = species susceptibility ($\beta_{\mathrm{duck}}$ for duck, 1 for chicken)
- $\varphi_j$ = spillover rate ($\varphi_{\mathrm{HRZ}}$ in the high-risk zone, $\varphi_{\mathrm{non}}$ elsewhere)
- $\psi(t)$ = Bateman-function spillover profile, rising from $t_0$ with rate $\sigma$ and decaying with rate $\delta$
- $\beta$ = spatial transmission rate, $\alpha$ = kernel range, $d_{ij}$ = distance between farms
- $w(\tau)$ = within-farm infectiousness, saturating from the latent period onward
- $p_{\mathrm{eff}}$ = effective movement transmission probability (reduced by pre-shipment testing in the HRZ, blocked in surveillance zones)
- $z_j(t)$ = surveillance zone modifier ($1 - \varepsilon$ if in zone, 1 otherwise)
- $\gamma_j(t)$ = confinement modifier (0.5 for broiler\_2 farms from 14 Jan onward; 0.5 for organic duck farms from 14 Jan to 13 Feb only; 1 otherwise)

### State variables

| Name | Symbol | Description |
|------|--------|-------------|
| Susceptible | S | Farm with no current infection, can become infected |
| Latent/Infectious | I | Infected farm, becomes infectious after latent period |
| Confirmed | C | Farm with confirmed outbreak, undergoing reactive culling |
| Removed | R | Farm depopulated (reactive or preventive cull) |

### Processes and parameter values

**Demographic processes**

- *Entry/exit:* Farm activity periods from the activity data set determine when farms have birds. No births or natural mortality are modelled explicitly.
- *Production cycles:* Captured through the activity matrix (date ranges of flock presence).
- *Movements:* Observed broiler\_1 to broiler\_2 movements during the fitting period; a gravity model (exponential distance kernel fitted to observed movements) generates synthetic movements for the prediction period.

**Epidemiological processes**

- *Transmission:* Three routes as described above — environmental spillover, spatial proximity, and live animal movements.
- *Progression:* Fixed compound delay of 11 days from infection to confirmation (3.5 days amplification + 5 days detection delay + 2 days suspicion-to-confirmation).
- *Detection:* Passive surveillance (84%), post-culling detection (13%), pre-shipment testing (3%), contact tracing (1%).

**Disease management**

| Measure | Included? | Implementation |
|---------|-----------|----------------|
| National standstill | N | — |
| Suspicion management | N | — |
| Reactive culling | Y | Farm removed at cull\_start (from data) or confirmation + 3 days |
| Contact tracing | N | Not modelled explicitly |
| Zoning: Movement bans | Y | No movements from farms in surveillance zones |
| Zoning: Enhanced biosecurity | Y | Susceptibility reduced by factor ε = 0.5 in surveillance zones |
| Preventive culling | Y | Phase 1: 1 km all species; Phase 2: 3 km chicken only |
| Pre-movement testing in HRZ | Y | Movement probability reduced by test sensitivity (90%) |
| Confinement | Y | 50% susceptibility reduction; broiler\_2 from 14 Jan onward; organic duck 14 Jan – 13 Feb only |
| Restocking ban | N | Not in baseline; explored in Q2 |

### Parameter table

| Symbol | Description | Unit | Value / Distribution | Source |
|--------|-------------|------|---------------------|--------|
| $t_0$ | Spillover onset day | day | Normal(15, 5), truncated [1, 75] | Estimated |
| $\varphi_{\mathrm{HRZ}}$ | Spillover rate, HRZ | day⁻¹ | LogNormal(log(10⁻³), 1) | Estimated |
| $\varphi_{\mathrm{non}}$ | Spillover rate, non-HRZ | day⁻¹ | LogNormal(log(10⁻⁴), 1) | Estimated |
| $\delta$ | Spillover decay rate | day⁻¹ | Exponential(1/50) | Estimated |
| $\sigma$ | Spillover rise rate | day⁻¹ | LogNormal(log(0.3), 1) | Estimated |
| $\beta_{\mathrm{duck}}$ | Duck susceptibility (relative) | — | Beta(2, 8) | Estimated |
| $\beta$ | Spatial transmission rate | day⁻¹ | LogNormal(log(10⁻⁴), 1.5) | Estimated |
| $\alpha$ | Spatial kernel range | m | LogNormal(log(3500), 0.5) | Estimated |
| $p_{\mathrm{mov}}$ | Movement transmission probability | — | Beta(2, 20) | Estimated |
| $\mu_E$ | Amplification period | day | 3.5 | Fixed (literature) |
| $\mu_{ID}$ | Detection delay | day | 5.0 | Fixed (literature) |
| $D_{TC}$ | Suspicion to confirmation | day | 2.0 | Fixed (data) |
| $\sigma_{\mathrm{test}}$ | Pre-shipment test sensitivity | — | 0.9 | Fixed (assumption) |
| $\varepsilon$ | Zone biosecurity reduction | — | 0.5 | Fixed (assumption) |
| $\gamma$ | Confinement factor | — | 0.5 | Fixed (assumption) |

### Parameter estimation

We fit the model using Hamiltonian Monte Carlo (NUTS) via Turing.jl, with 4 chains of 1000 post-warmup samples each (500 warmup). Initial values come from a random prior search (MAP approximation). Gradients are computed via ForwardDiff with chunk size 9.

Convergence was assessed by R-hat (all < 1.01) and effective sample size (all > 400).

### Posterior parameter estimates

| Parameter | Median | 95% CrI |
|-----------|--------|---------|
| $t_0$ | 9.5 | [5.9, 10.8] |
| $\varphi_{\mathrm{HRZ}}$ | 1.24 × 10⁻³ | [6.3 × 10⁻⁴, 2.1 × 10⁻³] |
| $\varphi_{\mathrm{non}}$ | 1.08 × 10⁻⁴ | [4.4 × 10⁻⁵, 2.1 × 10⁻⁴] |
| $\delta$ | 0.018 | [0.005, 0.036] |
| $\sigma$ | 0.135 | [0.027, 1.14] |
| $\beta_{\mathrm{duck}}$ | 0.173 | [0.142, 0.210] |
| $\beta$ | 0.0106 | [0.008, 0.014] |
| $\alpha$ | 5167 m | [4621, 5787] |
| $p_{\mathrm{mov}}$ | 0.079 | [0.011, 0.235] |

### Initial conditions

- **Simulation start date:** 1 December 2025 (day 1)
- **Initial infection seeding:** No explicit seed — spillover begins at estimated $t_0$ (posterior median ~9.5, corresponding to ~10 December)

### Simulations and outputs

- **Number of runs:** 500 posterior draws, each producing one forward trajectory
- **Random seed:** Fixed (MersenneTwister(2025)) for reproducibility
- **Outputs:** Daily new confirmed cases; per-farm infection probability; restocking safety profile; counterfactual outbreak totals; capacity sweep
- **Output format:** CSV files

---

## Third period results

### Q1: General description of the resolved epidemic

A total of 560 outbreaks were confirmed between 22 December 2025 and 7 April 2026 (107 days). The last case was confirmed on 7 April; no further outbreaks have been reported since.

**Distribution by species and production type:**

| Species | Production | Cases | % |
|---------|-----------|-------|---|
| Chicken | broiler\_2 | 255 | 45.5 |
| Chicken | broiler\_1 | 101 | 18.0 |
| Duck | conventional | 78 | 13.9 |
| Chicken | layer | 64 | 11.4 |
| Duck | organic | 62 | 11.1 |
| **Total** | | **560** | |

Chicken farms accounted for 75% of all cases (420/560), with stage 2 broilers alone making up 46%. Duck farms (conventional and organic combined) accounted for 25%.

**Temporal pattern:** Incidence peaked in the week of 18–24 January at approximately 18 cases per day, then declined steadily to 1–3 cases per day by early March. The final case was confirmed on 7 April. Several factors likely contributed to the decline: depletion of susceptible farms near existing clusters, preventive culling, confinement of broiler\_2 farms, and the end of the wildlife migration period.

**Spatial distribution:** Berks county accounted for 182 cases (33%), followed by Fulton (57), Cumberland (55), Susquehanna (49), Allegheny (39), and Indiana (39). These six counties together contained 75% of all cases. Cases were initially concentrated in the high-risk zone but spread island-wide during January.

**Detection methods:** Passive surveillance detected 84% of cases (468/560). Post-culling detection identified 74 cases (13%), pre-shipment testing 14 cases (3%), and contact tracing 4 cases (1%).

**Preventive culling:** 691 farms were preventively culled over the course of the epidemic, all of which had been completed by the data cut-off. 74 of these turned out to be infected (detected as post-culling cases), which means preventive culling was removing farms that were already incubating infection.

*(Spatial map provided as separate PDF)*

### Q2: When can restocking begin safely?

We simulated the epidemic forward from the end of the observation window (10 April, day 131) under a range of restocking start dates, from 10 April to 7 August 2026. For each candidate date, we ran 500 forward simulations (one per posterior draw) for 120 days and recorded whether any new infections occurred after restocking.

Under the fitted model, P(rebound) remains above 95% for all candidate dates within the sweep window. The predicted number of new cases declines with later restocking (median 982 at 10 April, 355 at 7 August) but does not approach zero.

This result is driven by the fitted spillover decay rate (δ = 0.018 day⁻¹, half-life ~39 days). The Bateman function spillover profile has not fully decayed by the end of the observation period: at day 131, residual spillover intensity is approximately 11% of its peak value. As long as any spillover persists, restocked farms face infection risk, and the estimated spatial transmission rate (β = 0.011, kernel range α ≈ 5.2 km) is sufficient to amplify isolated infections into local outbreaks.

**Interpretation and caveats.** The model lacks a mechanism for the abrupt end of spillover that the narrative describes ("the migration period has progressively come to an end"). The Bateman function produces a gradual exponential decay rather than a hard cutoff. In practice, if no new cases have been observed for more than two weeks after the last confirmed case, and if wildlife migration has ended, the actual rebound risk is likely much lower than the model predicts.

Surveillance is likely more informative than modelling here: if active monitoring of farms near former outbreak clusters detects no new infections through May, restocking could begin in June with low risk, starting with farms furthest from former outbreak foci. If the spillover profile were truncated to zero after day 131, farm-to-farm transmission alone would not sustain an epidemic, because all known infection chains have been broken.

### Q3: Outbreaks averted by preventive culling

We estimated the number of outbreaks averted by comparing paired counterfactual simulations from each posterior draw: one with preventive culling as implemented (phased: 1 km all-species through 13 January, then 3 km chicken-only), and one without any preventive culling. Both used the generative epidemic simulator starting from day 1 with fitted parameters.

Preventive culling averted a median of **130 outbreaks** (95% CrI: −195 to 503).

The wide credible interval includes negative values. In some posterior draws, the scenario without preventive culling produces fewer total cases because the susceptible farm population depletes faster through infection — more farms get infected and are removed through reactive culling instead, which can lead to faster local burnout. The median effect is positive: in most posterior draws, preventive culling reduced total case counts.

The with-culling simulations produced a median of 393 cases [187–644], consistent with the posterior predictive check (observed 560, PPC median 394 [204–636]). The without-culling simulations produced a median of 514 cases [233–926].

### Q4: Daily culling capacity for 95% eradication by end of phase 2

We swept daily culling capacity from 5 to 100 farms per day (plus unlimited) and asked: for each capacity level, what fraction of simulated epidemics are eradicated (no active infections remaining) by 13 February 2026 (day 75)?

At no capacity level did the model achieve eradication. P(eradicated by day 75) = 0 for all capacity levels tested, including unlimited.

This result follows from the model structure: the environmental spillover process continues to introduce new infections from the wildlife reservoir regardless of culling capacity. Even if every confirmed case were culled instantly, new farms would continue to be infected through spillover. The model estimates that spillover at day 75 is still at approximately 30% of its peak intensity, generating roughly 1–3 new spillover cases per day in the high-risk zone.

Culling alone, however fast, could not have eradicated the epidemic by February. Eradication would have required eliminating the wildlife source, either through the natural end of the migration period or through direct intervention in wild bird populations. Neither was available as a policy lever.

## Phase 2 prediction retrospective

Our phase 2 baseline prediction (382 new cases, 95% CrI: 261–525, over 28 days from 14 February) overpredicted by a factor of nearly five: the actual outcome was 81 cases, below the lower bound of our credible interval.

The primary cause was underestimation of the spillover decay rate δ. The phase 2 fit (on 75 days of data) estimated the epidemic trajectory at a point where incidence was declining but had not yet reached the tail. With only 75 days, the model could not distinguish between a slowly decaying epidemic and one that was about to collapse. The phase 3 fit (131 days) estimates δ = 0.018 [0.005–0.036], somewhat larger than the phase 2 estimate, reflecting the long tail of declining incidence.

The overprediction also reflects a structural limitation: the right-censoring correction treated recent infections as potentially unobserved, biasing the forward simulation toward higher initial infection counts.

## Strengths and limitations

The farm-level spatial resolution and three transmission routes allow scenario comparisons that would not be possible with a simpler model. The Bayesian framework propagates parameter uncertainty through to predictions.

The main limitation for phase 3 is the Bateman function spillover profile, which decays exponentially rather than abruptly. The fitted decay rate implies persistent spillover well beyond the observed epidemic, leading to overly pessimistic restocking assessments. A model with an explicit migration end date (informed by ornithological data) would give more practical answers to Q2.

The confinement factor (50% reduction) remains assumed rather than estimated. The Q3 and Q4 counterfactuals depend on the generative model reproducing observed outbreak dynamics, which it does only approximately (median 393 simulated vs 560 observed cases, though within the 95% CrI).

## Effort estimate

Approximately 6 person-days for phase 3, including model updates (confinement end date, counterfactual simulation capabilities, capacity-constrained culling), HPC runs, and report preparation.

## References

1. Abbott S, Sherratt K, Funk S. Infectious disease modelling with multiple data sources. 2025.
2. Carpenter B et al. Stan: A probabilistic programming language. J Stat Softw. 2017;76(1).
3. Ge H, Xu K, Ghahramani Z. Turing: a language for flexible probabilistic inference. AISTATS. 2018.
