---
bibliography: references.bib
---

# Process DAG

This document develops the process DAG — a representation of the latent transmission process that generates HPAI infections on Jolly Island farms.

The process DAG describes **how transmission occurs**, independent of how we observe it. Observation processes (surveillance, confirmation, reporting) are handled separately.

---

## Key Process Components

Based on the epidemic narrative and research questions, the transmission process involves:

1. **Wild bird spillover** — external introduction from infected migratory birds
2. **Farm-to-farm transmission** — local spread between poultry farms
3. **Within-farm infection dynamics** — progression from introduction to detectable outbreak
4. **Farm removal** — culling removes farms from susceptible/infectious pools
5. **Species heterogeneity** — potential differences in chicken vs duck susceptibility/infectiousness

---

## State Variables

### Population filter

Farms without birds present (downtime between production cycles) are excluded from the modelled population as a preprocessing step. Only farms active during the study period enter the model. Restocking is not modelled.

### Farm-level states

Each farm $j$ at time $t$ is in exactly one state $\text{State}_j(t) \in \{S, I, R\}$:

| State | Description |
|-------|-------------|
| $S$ (Susceptible) | Farm has birds present, not infected |
| $I$ (Infected) | Farm infected, potentially infectious to others |
| $R$ (Removed) | Farm culled or depopulated |

For convenience, we also define indicator variables: $S_j(t) = \mathbf{1}[\text{State}_j(t) = S]$, etc.

**Note**: "Infected" here means the farm has HPAI infection, not that it has been detected. Detection is an observation process.

### Aggregate states (derived)

- Total susceptible farms: $S(t) = \sum_j S_j(t)$
- Total infected farms: $I(t) = \sum_j I_j(t)$
- Total removed farms: $R(t) = \sum_j R_j(t)$

---

## Transmission Processes

**Time convention**: We use discrete time with $\Delta t = 1$ day. All hazards $\lambda_j(t)$ represent daily rates, and infection probability on day $t$ is $P(\text{infection}) = 1 - \exp(-\lambda_j(t))$.

### 1. Wild Bird Spillover

Farms in the high-risk zone (HRZ) can be infected by contact with infected wild birds.

**Temporal profile**: Wild bird migration is seasonal — the narrative states birds "migrate in early winter with stragglers seen through February". Spillover pressure is not constant but follows a temporal profile:

$$\text{hazard}_{\text{spillover},j}(t) = \begin{cases} \phi_{\text{hrz}} \cdot \psi(t) & \text{if } j \in \text{HRZ} \\ \phi_{\text{non}} \cdot \psi(t) & \text{otherwise} \end{cases}$$

where $\phi_{\text{hrz}}$ and $\phi_{\text{non}}$ are the spillover rate scalars for HRZ and non-HRZ farms respectively, and $\psi(t)$ is a flexible temporal profile estimated from the data.

**Temporal profile — piecewise constant**: We model $\psi(t)$ as a two-piece step function with an estimated changepoint:

$$\psi(t) = \begin{cases} 1 & \text{if } t \leq t_{\text{change}} \\ \rho & \text{if } t > t_{\text{change}} \end{cases}$$

where $t_{\text{change}}$ is the changepoint (estimated) and $\rho \in (0, 1]$ is the relative spillover intensity in the second period. The first period is the reference ($\psi = 1$), so $\phi_{\text{hrz}}$ and $\phi_{\text{non}}$ represent the spillover rates during the early high-pressure phase.

**Rationale**: No external wild bird surveillance data are available for Jolly Island. With 103 cases over 44 days, the data cannot resolve a detailed temporal profile. A two-piece model captures the key feature — spillover pressure was higher early in the outbreak (migration arrival) and declined later — without over-parameterising. The narrative ("wild birds migrate in early winter with stragglers seen through February") supports a declining profile but not a specific shape.

**Identifiability**: $t_{\text{change}}$ is identified from the changing proportion of HRZ vs non-HRZ cases over time. $\rho$ is identified from the HRZ excess in the second period relative to the first. The HRZ/non-HRZ contrast (~21% of farms in HRZ) separates spillover from local transmission.

**Spatial alternatives** (for later consideration):
- Continuous distance-to-wetland using `clc_32626.geojson`: $\phi \cdot \exp(-d_j / \rho)$

### 2. Farm-to-Farm Transmission

Infected farms can transmit to nearby susceptible farms.

**Spatial kernel**: For each susceptible farm $j$ at time $t$, the force of infection from farm-to-farm spread is:
$$\lambda_j^{\text{local}}(t) = \beta \cdot \sum_{i: I_i(t)=1} w_i(\tau_i) \cdot K(d_{ij})$$

where:
- $\beta$ is the baseline transmission rate
- $w_i(\tau_i)$ is the infectiousness of farm $i$ at time $\tau_i$ since infection
- $K(d_{ij})$ is a spatial kernel depending on distance $d_{ij}$ between farms
- Sum is over all currently infected farms

**Kernel options**:

- Exponential: $K(d) = \exp(-d/\alpha)$
- Power law: $K(d) = (1 + d/\alpha)^{-\gamma}$
- Cauchy: $K(d) = 1/(1 + (d/\alpha)^2)$

**Decision**: Use **Cauchy** kernel $K(d) = 1/(1 + (d/\alpha)^2)$. HPAI transmission literature consistently favours fat-tailed kernels: @Boender2007 fitted a Cauchy-type kernel ($\alpha \approx 2.1$ km) to the Netherlands H7N7 outbreak; @Boender2023 found common features in spatial livestock disease transmission parameters; @Seymour2021 used Bayesian nonparametric kernel estimation for the Netherlands H7N7 outbreak and found sub-exponential tails. Fat tails better capture occasional long-range jumps via unrecorded pathways (wind, shared services) absorbed into the spatial kernel. The reparameterisation $\beta_0 = \beta \cdot K(d_0)$ applies equally. Exponential and power-law kernels are retained as alternatives for model comparison.

**Infectiousness profile**: Farm infectiousness increases over time as within-farm prevalence grows, with a hard latent period before any between-farm transmission is possible:
$$w(\tau) = \begin{cases} 0 & \text{if } \tau < \tau_{\min} \\ 1 - \exp(-r \cdot (\tau - \tau_{\min})) & \text{if } \tau \geq \tau_{\min} \end{cases}$$
where $\tau_{\min} = 1$ day is the minimum latent period and $r \approx 1.0$/day is the within-farm growth rate (from mortality ledgers). The hard latent period reflects that within-flock spread must occur before a farm generates sufficient environmental contamination for between-farm transmission (latent periods consistently <2 days across HPAI subtypes; @Guinat2023). After $\tau_{\min}$, infectiousness ramps up and saturates to 1.

### 3. Movement-Based Transmission

Broiler_1 → broiler_2 movements can transmit infection via transport of infected birds.

**Process**: Let $M_{i \to j}(t) \in \{0,1,2,...\}$ be the number of movements from farm $i$ to farm $j$ on day $t$. The movement hazard for farm $j$ is:
$$\lambda_{\text{movement},j}(t) = \sum_{i: I_i(t)=1} M_{i \to j}(t) \cdot p_{\text{eff}}(i,t) \cdot w(t - T_i^I)$$

where $p_{\text{eff}}(i,t)$ incorporates both pre-shipment testing and zone-based blocking:

- If farm $i$ is in a regulated zone at time $t$: $p_{\text{eff}}(i,t) = 0$ (movements blocked)
- Else if farm $i$ is in HRZ: $p_{\text{eff}}(i,t) = p_{\text{mov}} \cdot (1 - \sigma_{\text{test}})$ (testing intercepts with probability $\sigma_{\text{test}} = 0.9$)
- Else: $p_{\text{eff}}(i,t) = p_{\text{mov}}$

**Data**: 7,187 recorded movements (broiler_1 → broiler_2 only). Other movement types (to slaughter, equipment, personnel) are not recorded.

**Evidence for inclusion**:

- 6/103 cases detected via pre-shipment testing — these are farms *intercepted before transmitting*, indicating movement is a real transmission pathway
- Pre-shipment testing sensitivity <100%, so some infected movements may get through
- @Yoo2021 found ~30% of Korean H5N6 transmission via vehicle movements
- Movement provides a mechanistically distinct pathway from spatial proximity (can explain long-range jumps)

**Modifiers**:

- Pre-shipment testing in HRZ intercepts infectious farms with probability $\sigma_{\text{test}}$ (assume 0.9)
- Movements from farms in regulated zones are blocked
- Effective probability: $p_{\text{eff}} = p_{\text{mov}} \cdot (1 - \sigma_{\text{test}})$ for HRZ sources
- Otherwise (non-HRZ, non-regulated source): $p_{\text{eff}} = p_{\text{mov}}$ (no modification)

**Decision**: Include movement transmission with **fixed $p_{\text{mov}} = 0.01$** (expert assumption — calibrated so that at typical movement volumes the movement pathway contributes a share broadly consistent with @Yoo2021's ~30% attribution for Korean H5N6; the cited study reports population-level attribution, not a per-movement probability). Sensitivity analysis should vary $p_{\text{mov}}$ over [0.001, 0.05]. The parameter may not be identifiable (confounded with spatial kernel for nearby farm pairs), but including it avoids misattributing movement transmission to the spatial kernel.

**Limitation**: Only broiler_1 → broiler_2 movements are recorded. Other transmission via shared equipment, personnel, or feed trucks is absorbed into the spatial kernel.

### 4. Species-Specific Susceptibility

Chickens and ducks may have different susceptibility to infection.

**Process**: Modify total force of infection by species:
$$\lambda_j(t) = \beta_s \cdot (\lambda_j^{\text{spillover}}(t) + \lambda_j^{\text{local}}(t) + \lambda_j^{\text{movement}}(t))$$

where $\beta_s$ (≡ `β_species[j]` in pseudocode) is species-specific:

- $\beta_{\text{chicken}} = 1$ (reference)
- $\beta_{\text{duck}} \in (0, 2]$ (to be estimated)

**Rationale for bidirectional bound**: Ducks may be either more or less susceptible to infection than chickens depending on the pathway. Ducks are generally more resistant to HPAI *mortality* (@Smith2015 — chickens lack the RIG-I innate immune receptor), but may be *more* susceptible to *infection*, particularly via wild-bird spillover since they share aquatic habitats with migratory waterfowl. The (0, 2] bound allows the data to inform direction. See @PantinJackwood2013 for strain-dependent species differences.

**Interpretation caveat**: $\beta_{\text{duck}}$ conflates true susceptibility with detectability (see `step01_research_questions.md`, Q3).

### 5. Infection Dynamics

Once infected, a farm progresses through infection until removed.

**Infection time**: Farm $j$ becomes infected at time $T_j^I$.

**Removal time**: Farm $j$ is removed (culled) at time $T_j^R$.

**Infectious period**: Farm is infectious from $T_j^I$ to $T_j^R$.

The removal process is partially observed (cull dates in `cases.csv`) but is affected by interventions.

---

## Intervention Processes

### Reactive Culling

Confirmed farms are culled.

**Modelling choice**: We treat confirmation times $T_j^C$ as observed inputs (from the observation model). Removal is then a deterministic process:
$$T_j^R = T_j^C + \delta_{\text{reactive}}$$

where $\delta_{\text{reactive}}$ is a fixed or estimated delay (median ~2 days from data, but longer during capacity constraints from 6 Jan).

**Limitation**: A single $\delta_{\text{reactive}}$ averages over pre- and post-capacity-constraint periods, potentially underestimating infectious periods during high incidence (post-6 Jan) and inflating $\beta$. A time-varying $\delta_{\text{reactive}}(t)$ or explicit capacity model (analogous to the preventive culling treatment) could address this; see Complexity Option 4.

### Preventive Culling

From 1 Jan, farms within 1km of confirmed cases are preventively culled.

**Process**: If farm $j$ is within 1km of a confirmed farm and $j$ is susceptible:
$$T_j^R = \max(T_{\text{trigger}}, T_{\text{capacity}})$$

where:

- $T_{\text{trigger}}$ is when the trigger case was confirmed
- $T_{\text{capacity}} = T_{\text{trigger}} + \delta_{\text{prev}}$, i.e. the capacity delay is modelled via $\delta_{\text{prev}}$ (estimated from the 12 complete records). Since $T_j^R = \max(T_{\text{trigger}}, T_{\text{trigger}} + \delta_{\text{prev}}) = T_{\text{trigger}} + \delta_{\text{prev}}$ when $\delta_{\text{prev}} \geq 0$, the process simplifies to $T_j^R = T_{\text{trigger}} + \delta_{\text{prev}}$

**Note**: 77% of preventive cull dates are missing — imputed using the same $\delta_{\text{prev}}$ parameter.

**Limitation**: $\delta_{\text{prev}}$ estimated from only 12 records carries substantial uncertainty that may bias estimates of β and culling effectiveness. Consider: (a) informative prior borrowing from reactive culling distribution, or (b) sensitivity analysis over plausible $\delta_{\text{prev}}$ range.

### Zones

Regulated zones (3km protection, 10km surveillance) affect:

- Movement restrictions
- Enhanced detection (observation process, not transmission)

For the process DAG, zones primarily affect which farms can be preventively culled.

---

## Simplified Process DAG (Initial)

For initial model development, we propose a simplified process:

```text
Parameters (estimated):
  φ_hrz   = spillover rate scalar (HRZ farms)
              Prior: LogNormal(log(0.01), 2)  (median 0.01/day; wide — no direct literature estimates for per-farm spillover rates; h₀ ≈ 10⁻³–10⁻² from farm-to-farm literature as rough anchor)
  φ_non   = spillover rate scalar (non-HRZ farms)
              Prior: LogNormal(log(0.01), 2)  (same as φ_hrz; data separates HRZ from non-HRZ)
  t_change = spillover changepoint (day)
              Prior: Normal(1 Jan 2026, SD = 10 days)  (narrative: "migrate in early winter with stragglers through February")
  ρ       = relative spillover intensity after changepoint
              Prior: Beta(2, 2)  (centred at 0.5; narrative supports decline but not specific magnitude)
  β       = farm-to-farm transmission rate
              Prior: Exponential(mean = 0.1)  (weakly informative; see α–β identifiability note)
  α       = spatial kernel scale (km)
              Prior: LogNormal(log(2), 1)  (median ≈ 2 km; broad enough for data to inform)
  β_duck  = relative susceptibility of ducks
              Prior: LogNormal(0, 0.5) truncated to (0,2] (centred at 1, allows either direction)
              Fallback: if posterior is prior-dominated, treat as scenario parameter (fixed per run: 0.5, 1.0, 1.5)

Parameters (fixed):
  τ_min   = hard latent period (= 1 day; within-flock spread before between-farm transmission)
  r       = within-farm growth rate (= 1.0/day from mortality ledgers)
  p_mov   = per-movement transmission probability (= 0.01)
  σ_test  = pre-shipment testing sensitivity (= 0.9)

Parameters (imputed / estimated from data summaries):
  δ_reactive = confirmation-to-removal delay for reactive culling (median ~2 days from cases.csv)
  δ_prev     = trigger-to-removal delay for preventive culling (estimated from 12 complete records; see limitation note)

States:
  S_j(t)  = farm j susceptible at t
  I_j(t)  = farm j infected at t
  R_j(t)  = farm j removed at t
  T_j^I   = infection time of farm j

Process:
  For each susceptible farm j at time t:

    # Spillover (piecewise constant temporal profile)
    ψ(t) = 1 if t ≤ t_change, else ρ
    hazard_spillover = φ_hrz × ψ(t)  if j ∈ HRZ, else φ_non × ψ(t)

    # Within-farm infectiousness (shared by local and movement pathways)
    w(τ) = 0 if τ < τ_min, else 1 - exp(-r × (τ - τ_min))   # hard latent period, then saturates to 1

    # Local transmission (spatial kernel, sum over infected farms)
    hazard_local = β × Σ_{i: I_i(t)=1} w(t - T_i^I) × K(d_ij)

    # Movement transmission (broiler_1 → broiler_2, sum over infected sources)
    hazard_movement = Σ_{i: I_i(t)=1} M_i→j(t) × p_eff(i,t) × w(t - T_i^I)
    # where M_i→j(t) = number of movements from i to j on day t
    # p_eff(i,t) = 0 if i in regulated zone, else p_mov×(1-σ_test) if i in HRZ, else p_mov

    # Total hazard (species modifier applies to all pathways)
    λ_j(t) = β_species[j] × (hazard_spillover + hazard_local + hazard_movement)

    # Infection event
    P(infection at t) = 1 - exp(-λ_j(t))

  For each confirmed infected farm j (reactive culling):
    T_j^R = T_j^C + δ_reactive

  For each susceptible farm j within 1 km of a confirmed farm (preventive culling, from 1 Jan):
    # S → R directly, bypassing I
    T_j^R = T_trigger + δ_prev
    # where T_trigger = confirmation time of trigger case
    # Precedence: if farm j transitions S→I before scheduled preventive T_j^R,
    #   cancel preventive schedule; treat as reactive case (T_j^R = T_j^C + δ_reactive)
```

---

## DAG Representation

```text
    PARAMETERS                                          FIXED                COVARIATES                       DATA
    [φ_hrz, φ_non, t_change, ρ, β, α, β_duck]    {p_mov, r, τ_min, σ_test}      {HRZ}  {species}  {location}   {M_{i→j}(t): movements}
              |              |        |                  |                    |        |          |                  |
              v              v        v                  v                    v        v          v                  v
    (Spillover) + (Local transmission) + (Movement transmission) ←---------
          \              |                    /
           \             |                   /
            \            v                  /
             \--→ [Force of infection λ_j(t)] ←--/
                             |
                             v
                   [Farm infection state]
                        S → I → R
                             ^
                             |
                  [Interventions: culling]
```

**Nodes**:

- `[Square brackets]`: estimated parameters or latent states
- `{Curly braces}`: fixed covariates or fixed parameters (observed or set from literature)
- `(Parentheses)`: processes/transformations

**Key distinction**: HRZ and species are covariates that modify how parameters act, not parameters themselves. We estimate φ_hrz and φ_non (spillover rate scalars), t_change and ρ (spillover temporal profile), β (transmission rate), α (kernel scale), and β_duck (relative susceptibility). Movement transmission probability p_mov is fixed at 0.01 based on literature (mechanistically important but not identifiable from these data). Other fixed inputs: τ_min (hard latent period), r (within-farm growth rate), and σ_test (pre-shipment testing sensitivity).

---

## Unobserved Quantities

The process DAG involves latent quantities that are not directly observed:

| Quantity | Symbol | Status |
|----------|--------|--------|
| Infection times | $T_j^I$ | Latent — must infer |
| Spillover events | — | Latent — inferred from HRZ/non-HRZ pattern |
| Transmission pairs | — | Latent — not identifiable without genomics |
| Within-farm prevalence | — | Not modelled (farm-level abstraction) |

---

## Parameter Identifiability

| Parameter(s) | Identifiable? | Notes |
|--------------|---------------|-------|
| φ_hrz (HRZ spillover) | Conditional | Identified from HRZ/non-HRZ case contrast; early outbreak confounded with β |
| φ_non (non-HRZ spillover) | Weak | Identified from non-HRZ cases not explained by local transmission; expected to be small |
| t_change (changepoint) | Conditional | Identified from changing HRZ/non-HRZ case ratio over time |
| ρ (post-change intensity) | Conditional | Identified from HRZ excess in second period relative to first |
| β (transmission) | Conditional | Identified from non-HRZ cases + spatial-temporal clustering |
| α (kernel scale) | Weak | Reparameterise: estimate β₀ = β·K(d₀) and α separately |
| β_duck (species) | **Weak** | Conflates susceptibility, infectiousness, and detectability; LogNormal(0, 0.5) prior on (0,2]. If posterior is prior-dominated, fall back to scenario analysis (fixed β_duck ∈ {0.5, 1.0, 1.5}) |
| r (growth rate) | Fixed | Set from mortality ledger data (= 1.0/day) |
| p_mov (movement) | Fixed | Set from literature (= 0.01); not identifiable, confounded with spatial kernel |

### φ_hrz / φ_non vs β

Structurally identifiable because:

- φ_hrz/φ_non: hazard depending on HRZ membership, modulated by the piecewise constant temporal profile ψ(t)
- β: proximity-dependent hazard varying with distance to specific infected farms

**Key diagnostic**: Infections outside the HRZ can arise from local transmission (β) or non-HRZ spillover (φ_non). Since φ_non is expected to be small, non-HRZ infections primarily identify β. The HRZ excess then identifies φ_hrz.

**Practical concerns**:

- Early outbreak: few infected farms means low local pressure everywhere, so early infections inform spillover + β jointly rather than separately
- Culling delays: capacity constraints extend infectious periods, inflating apparent β
- Undetected farms: if surveillance misses cases, apparent spillover may be local transmission from undetected sources

Simulation-based identifiability checks recommended before drawing conclusions.

### α vs β

Kernel scale and transmission rate trade off: higher β with smaller α can produce similar patterns to lower β with larger α.

**Resolution**: Reparameterise to a compound parameter that is well-identified:

- Estimate $\beta_0 = \beta \cdot K(d_0)$ — transmission rate at reference distance $d_0$ (e.g., median inter-farm distance)
- Estimate $\alpha$ — kernel decay scale
- Derive $\beta = \beta_0 / K(d_0)$

This separates "how much transmission" ($\beta_0$) from "how far" ($\alpha$). The compound $\beta_0$ is identified from overall transmission intensity; $\alpha$ is identified from the spatial decay pattern.

### β_duck

The apparent species difference conflates:

- True biological susceptibility
- Viral shedding intensity/duration (ducks may shed longer)
- Detection probability (lower mortality → delayed detection)

Ducks may be *more* infectious despite lower apparent susceptibility. The research questions address this via scenario analysis.

---

## Simplifying Assumptions

1. **Homogeneous mixing within species/location** — no farm-level heterogeneity beyond species and location
2. **Fixed movement transmission probability** — $p_{\text{mov}} = 0.01$ not estimated (confounded with spatial kernel)
3. **Piecewise constant spillover profile** — two periods separated by an estimated changepoint; captures declining spillover without over-parameterising
4. **Binary HRZ** — spillover risk is HRZ vs non-HRZ (continuous distance-to-water preferred if feasible)
5. **Hard latent period via infectiousness profile** — $w(\tau) = 0$ for $\tau < \tau_{\min}$ (1 day), then ramp-up; no separate S→E→I compartment
6. **Fixed within-farm growth rate** — $r = 1.0$/day from mortality ledgers
7. **Single species modifier across pathways** — $\beta_{\text{duck}} \in (0, 2]$ applied uniformly; pathway-specific modifiers deferred

---

## Complexity Options (for refinement)

If the simplified model is insufficient:

**High priority** (recommended before fitting):
1. **Longer latent period** — increase $\tau_{\min}$ from 1 to 2 days, or estimate $\tau_{\min}$ from data if sufficient information exists.
2. **Pathway-specific species effect** — separate $\beta_{\text{duck}}^{\text{spillover}}$, $\beta_{\text{duck}}^{\text{local}}$, $\beta_{\text{duck}}^{\text{movement}}$ if a single modifier is insufficient

**Medium priority**:
3. **Continuous spillover risk** — distance to wetlands instead of binary HRZ (using `clc_32626.geojson`)
4. **Culling capacity model** — explicitly model queue and delays as function of outbreak intensity
5. **Estimate $p_{\text{mov}}$** — if movement transmission appears important, consider estimating rather than fixing

**Lower priority**:
6. **Shared services** — feed trucks, veterinary personnel, equipment (if data available)
7. **Farm-level heterogeneity** — random effects on susceptibility/infectiousness
8. **Wind-borne transmission** — directional kernel if spatial patterns suggest it
9. **More flexible spillover profile** — spline or Gaussian process if the piecewise constant is too rigid and more data become available

---

## Link to Research Questions

| Question | Process DAG components required |
|----------|-------------------------------|
| Q1 | None (descriptive) |
| Q2 | Full process: spillover + local transmission + removal |
| Q3 | Species-specific susceptibility ($\beta_{\text{duck}}$): posterior estimate if identifiable, otherwise scenario comparison |
| Q4 | Preventive culling process (modify intervention rules) |
| Q5 | Reactive culling delays (modify $\delta_{\text{reactive}}$) |

---

