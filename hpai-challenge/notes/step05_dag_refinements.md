# Refining the Model DAGs (Workflow Step 5)

Following Abbott et al., this step revisits the process DAG (step 02) and observation DAG (step 04) to ensure alignment between what we want to model and what the data can support. We adjust complexity — both up and down — based on identifiability analysis, data constraints, and domain knowledge.

---

## Process DAG Refinements

### 1. Spillover temporal profile

**Step 02 specification**: piecewise constant $\psi(t)$ with estimated changepoint $t_{\text{change}}$ and post-change decline factor $\rho$. Four spillover parameters: $\phi$, $\eta$, $t_{\text{change}}$, $\rho$.

**Problem**: the piecewise constant profile assumes spillover is active from the start of the simulation period ($\psi = 1$ for $t \leq t_{\text{change}}$). This means any farm can be infected from day 1 — there is no mechanism for zero spillover before wild bird arrival. Since the simulation window starts 1 Dec 2025 but the first case is not confirmed until 22 Dec 2025 (back-calculated infection ~12–13 Dec 2025), the model will attribute ~2 weeks of accumulated spillover hazard to all farms before any cases occur. This inflates the baseline spillover exposure and risks confounding spillover with local transmission.

**Refinement**: replace the piecewise constant profile with onset + exponential decay:

$$\psi(t) = \begin{cases} 0 & \text{if } t < t_0 \\ \exp(-\delta \cdot (t - t_0)) & \text{if } t \geq t_0 \end{cases}$$

where $t_0$ is the onset day (estimated) and $\delta$ is the post-onset decay rate (estimated). Spillover pressure is at its maximum on day $t_0$ and declines thereafter.

**Rationale**: the narrative states birds "migrate in early winter with stragglers seen through February" — consistent with a relatively sharp arrival followed by declining pressure. A hard onset at $t_0$ with exponential decay is simpler than the piecewise constant specification (2 spillover-shape parameters instead of 2) and better captures the key feature: zero pressure before migration arrival.

**Spillover parameterisation**: step 02 used $\phi$ (HRZ rate) and $\eta$ (non-HRZ reduction factor, $\phi_{\text{non}} = \phi \times \eta$). An equivalent parameterisation estimates $\phi_{\text{hrz}}$ and $\phi_{\text{non}}$ directly, avoiding the implicit constraint $\phi_{\text{non}} < \phi_{\text{hrz}}$ imposed by $\eta \in (0, 1)$. The direct parameterisation is preferred — the data should determine the relative rates.

**Net effect on estimated parameters**: step 02 had $\phi$, $\eta$, $t_{\text{change}}$, $\rho$ (4 spillover-shape parameters). The refined model has $\phi_{\text{hrz}}$, $\phi_{\text{non}}$, $t_0$, $\delta$ (4 parameters). Same count, different structure.

### 2. Spatial kernel form

**Step 02 specification**: Cauchy kernel $K(d) = 1/(1 + (d/\alpha)^2)$, motivated by HPAI literature favouring fat-tailed kernels.

**Refinement**: start with the exponential kernel $K(d) = \exp(-d/\alpha)$. The exponential is simpler, computationally cheaper, and well-established in spatial epidemic models. The Cauchy kernel's advantage — heavier tails capturing occasional long-range jumps — matters most when the spatial kernel carries substantial transmission signal. With 103 cases over 24 days on an island with ~9000 farms, the data may not resolve kernel tail behaviour. Movement transmission (included separately) already provides a mechanism for long-range jumps.

**Action**: retain Cauchy and power-law kernels as alternatives for model comparison via posterior predictive checks. If the exponential kernel produces spatial patterns inconsistent with observed case clustering, switch to Cauchy.

### 3. Identifiability analysis

Step 02 identified key identifiability concerns. We expand that analysis here with proposed resolutions.

| Parameter pair | Identifiability concern | Resolution |
|---|---|---|
| $\phi$ vs $\beta$ (spillover vs local) | Early outbreak: few infected farms → low local pressure → early infections inform spillover + $\beta$ jointly | HRZ/non-HRZ contrast (~21% of farms in HRZ) separates the two; $t_0$ onset removes pre-arrival confounding |
| $\beta$ vs $\alpha$ (rate vs scale) | Higher $\beta$ with smaller $\alpha$ mimics lower $\beta$ with larger $\alpha$ | Reparameterise: estimate $\beta_0 = \beta \cdot K(d_0)$ and $\alpha$ separately (see §4) |
| $\beta_{\text{duck}}$ | Conflates susceptibility, infectiousness, detectability | Scenario analysis fallback: if posterior 95% CrI covers >80% of the prior range, fix $\beta_{\text{duck}} \in \{0.5, 1.0, 1.5\}$ per run |
| $t_0$ vs $\delta$ (onset vs decay) | Later $t_0$ with slower $\delta$ can produce similar cumulative hazard to earlier $t_0$ with faster $\delta$ | Informative priors on both; $t_0$ is anchored by first-case timing |
| $p_{\text{mov}}$ | Confounded with spatial kernel for nearby farm pairs | Fix at 0.01; include pathway for mechanistic completeness |

### 4. Reparameterisation ($\beta_0 / \alpha$)

Step 02 proposed estimating $\beta_0 = \beta \cdot K(d_0)$ at reference distance $d_0$ to resolve the $\beta$–$\alpha$ trade-off. This separates "how much transmission" ($\beta_0$, identified from overall transmission intensity) from "how far" ($\alpha$, identified from spatial decay pattern).

**Decision**: implement the $\beta_0$ reparameterisation from the start. Set $d_0$ to the median inter-farm distance among observed case-farm neighbour pairs. For the exponential kernel: $\beta = \beta_0 / \exp(-d_0 / \alpha)$.

### 5. Complexity adjustments

**Complexity retained from step 02**:
- Three transmission pathways (spillover, local, movement)
- Species-specific susceptibility ($\beta_{\text{duck}}$)
- Zone effects on susceptibility
- Preventive and reactive culling
- Within-farm infectiousness ramp-up

**Complexity reduced**:
- Spillover profile: onset + decay instead of piecewise constant (same parameter count, simpler structure)
- Spatial kernel: exponential instead of Cauchy (simpler, defer comparison)

**Complexity deferred** (for later iterations if base model is insufficient):
- Continuous spillover risk (distance to wetlands instead of binary HRZ)
- Pathway-specific species effects
- Culling capacity model (time-varying $\delta_{\text{reactive}}$)
- Farm-level random effects
- Estimating $p_{\text{mov}}$

---

## Observation DAG Refinements

### 1. Compound delay structure

Step 04 committed to a single compound delay $\delta = T^C - T^I$ rather than decomposing into amplification ($T^D - T^I$), recognition ($T^S - T^D$), and confirmation ($T^C - T^S$) components. Individual components are not separately identifiable from the available data.

**Refinement**: fix the compound delay at a literature-informed value. The total delay from infection to confirmation combines:
- Latent/amplification period: 2–5 days (literature)
- Detection delay: 3–7 days (clinical recognition)
- Suspicion to confirmation: ~2 days (from data, median $T^C - T^S$)

These enter the likelihood through integer rounding of back-calculated infection days, giving zero gradient for HMC/NUTS — they cannot be estimated by gradient-based samplers. Sensitivity analysis over delay values is required (see below).

### 2. Inference approach

Step 04 committed to data augmentation MCMC — sampling latent infection times $T^I$ jointly with model parameters. This properly propagates uncertainty in infection times through to parameter estimates and avoids the zero-gradient problem of fixed delays.

**Refinement**: implement data augmentation as the target approach. As a pragmatic first step, deterministic back-calculation with fixed delays is acceptable for initial model development, provided sensitivity analysis is conducted. The transition to data augmentation should follow once the base model structure is validated.

### 3. Stationarity of detection delays

Step 04 flagged that detection delays may shorten as the outbreak progresses (heightened surveillance, faster farmer recognition after initial awareness). This violates the stationarity assumption.

**Refinement**: conduct sensitivity analysis with two delay regimes:
- **Stationary**: constant delays throughout
- **Non-stationary**: shorter delays after 1 Jan 2026 (when surveillance zones and heightened awareness are in place)

If key conclusions ($\beta_{\text{duck}}$, spillover/local partitioning) are robust to this change, the stationary assumption is acceptable.

### 4. Delay sensitivity analysis

Fixed delay values directly determine back-calculated infection times, propagating to all downstream parameter estimates. This is the most consequential fixed-parameter choice in the model.

**Refinement**: conduct a grid sensitivity analysis over plausible delay values and report how posteriors for transmission and spillover parameters vary.

---

## Candidate DAGs

Following the workflow's recommendation to proceed with multiple candidate DAGs when mechanisms are uncertain, we identify two candidates:

**DAG A (base)**: onset + exponential decay spillover, exponential spatial kernel, $\beta_0/\alpha$ reparameterisation, 7 estimated parameters.

**DAG B (fat-tailed)**: as DAG A but with Cauchy spatial kernel. Same parameter count; different tail behaviour for local transmission.

Both candidates proceed through the remaining workflow steps (modularisation, inference, validation). Model comparison via posterior predictive checks determines which better reproduces observed spatial patterns.

---

## Consolidated Parameter Table

### Estimated parameters (7)

| Parameter | Symbol | Prior | Role |
|---|---|---|---|
| Spillover onset day | $t_0$ | $\text{Normal}(15, 5)$ | Day spillover begins ($t = 1$ is 1 Dec) |
| HRZ spillover rate | $\phi_{\text{hrz}}$ | $\text{LogNormal}(\log(10^{-3}), 1.0)$ | Daily per-farm spillover in HRZ at onset |
| Non-HRZ spillover rate | $\phi_{\text{non}}$ | $\text{LogNormal}(\log(10^{-4}), 1.0)$ | Daily per-farm spillover outside HRZ |
| Spillover decay rate | $\delta$ | $\text{Exponential}(\text{rate} = 50\ \text{day}^{-1})$ (mean $= 0.02\ \text{day}^{-1}$) | Post-onset decline in spillover |
| Transmission rate at reference distance | $\beta_0$ | $\text{LogNormal}(\log(10^{-4}), 1.5)$ | Farm-to-farm transmission intensity at $d_0$ |
| Spatial kernel scale | $\alpha$ | $\text{LogNormal}(\log(3500), 0.5)$ | Characteristic distance (metres) |
| Duck susceptibility | $\beta_{\text{duck}}$ | $\text{Beta}(2, 8)$ | Relative susceptibility (chicken = 1) |

### Fixed parameters

| Parameter | Symbol | Value | Source |
|---|---|---|---|
| Within-farm growth rate | $r$ | 1.0/day | Mortality ledgers (3 farms) |
| Hard latent period | $\tau_{\min}$ | 1 day | Literature; effective 2-day discrete-time latent period |
| Movement transmission probability | $p_{\text{mov}}$ | 0.01 | Expert assumption; not identifiable |
| Pre-shipment testing sensitivity | $\sigma_{\text{test}}$ | 0.9 | Assumed |
| Zone biosecurity reduction | $\varepsilon$ | 0.5 | Assumed |
| Surveillance zone radius | — | 10 km | Regulatory |
| Zone duration | — | 28 days | Regulatory |
| Preventive cull radius | — | 1 km | Regulatory (from 1 Jan) |
| Reference distance | $d_0$ | Computed at data preparation time | Median inter-farm distance among case-farm neighbour pairs (see §4) |

### Fixed delays (subject to sensitivity analysis)

| Parameter | Symbol | Default | Range for sensitivity |
|---|---|---|---|
| Latent period | $\mu_E$ | 3.5 days | {2, 3, 4, 5} |
| Detection delay | $\mu_{ID}$ | 5.0 days | {3, 4, 5, 6, 7} |
| Suspicion to confirmation | $D \to C$ | 2 days | Fixed from data |

### Imputed from data

| Parameter | Symbol | Source | Notes |
|---|---|---|---|
| Reactive cull delay | $\delta_{\text{reactive}}$ | $T^R - T^C$ from cases.csv | Per-case from data; median ~2 days |
| Preventive cull delay | $\delta_{\text{prev}}$ | 12 complete records | Imputed for 40/52 missing records |

---

## Refined Process DAG (Pseudocode)

Changes from step 02: (1) spillover uses onset $t_0$ + decay $\delta$ instead of piecewise constant, (2) $\phi_{\text{hrz}}$/$\phi_{\text{non}}$ estimated directly instead of $\phi$/$\eta$, (3) $\beta_0$ reparameterisation, (4) exponential kernel as default.

```text
Parameters (estimated):
  t₀       = spillover onset day
  φ_hrz    = HRZ spillover rate (at onset)
  φ_non    = non-HRZ spillover rate (at onset)
  δ        = spillover decay rate (post-onset)
  β₀       = transmission rate at reference distance d₀
  α        = spatial kernel scale (metres)
  β_duck   = relative duck susceptibility

Parameters (fixed):
  τ_min    = hard latent period (= 1 day)
  r        = within-farm growth rate (= 1.0/day)
  p_mov    = per-movement transmission probability (= 0.01)
  σ_test   = pre-shipment testing sensitivity (= 0.9)
  ε        = zone biosecurity reduction (= 0.5)
  d₀       = reference distance for β₀

Process:
  For each susceptible farm j at time t:

    # Spillover (onset + exponential decay)
    ψ(t) = 0 if t < t₀, else exp(-δ × (t - t₀))
    hazard_spillover = φ_j × ψ(t)
    # where φ_j = φ_hrz if j ∈ HRZ, else φ_non

    # Within-farm infectiousness (hard latent period, then saturating ramp-up)
    w(τ) = 0 if τ < τ_min, else 1 - exp(-r × (τ - τ_min))

    # Local transmission (exponential kernel; β derived from β₀)
    β = β₀ / exp(-d₀ / α)
    hazard_local = β × Σ_{i: I_i(t)=1} w(t - T_i^I) × exp(-d_ij / α)

    # Movement transmission
    hazard_movement = Σ_{i: I_i(t)=1} M_i→j(t) × p_eff(i,t) × w(t - T_i^I)

    # Total hazard (species modifier applies to all pathways)
    λ_j(t) = β_species[j] × (hazard_spillover + hazard_local + hazard_movement)

    # Zone modifier
    if j in surveillance zone: λ_j(t) *= (1 - ε)

    # Infection event
    P(infection on day t) = 1 - exp(-λ_j(t))
```

---
