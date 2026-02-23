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

### Farm-level states

Each farm $j$ at time $t$ is in exactly one state $\text{State}_j(t) \in \{S, I, R, A\}$:

| State | Description |
|-------|-------------|
| $S$ (Susceptible) | Farm has birds present, not infected |
| $I$ (Infected) | Farm infected, potentially infectious to others |
| $R$ (Removed) | Farm culled or depopulated |
| $A$ (Inactive) | No birds present (downtime) |

For convenience, we also define indicator variables: $S_j(t) = \mathbf{1}[\text{State}_j(t) = S]$, etc.

**Transition note**: $A$ is treated as a static covariate — farms in downtime at the start of the modelled period remain inactive throughout and are excluded from $S \to I \to R$ dynamics. Restocking (A → S) is not modelled.

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

$$\text{hazard}_{\text{spillover},j}(t) = \begin{cases} \epsilon \cdot \psi(t) & \text{if } j \in \text{HRZ} \\ 0 & \text{otherwise} \end{cases}$$

where $\epsilon$ is the peak daily spillover hazard and $\psi(t)$ captures the temporal profile.

**Temporal profile options**:

| Profile | Parameters | Form |
|---------|------------|------|
| Constant | ε | $\psi(t) = 1$ — ignores seasonality |
| Step onset | ε, t₀ | $\psi(t) = 0$ for $t < t_0$, $= 1$ for $t \geq t_0$ |
| Onset + decay | ε, t₀, δ | $\psi(t) = 0$ for $t < t_0$, $= \exp(-\delta(t - t_0))$ for $t \geq t_0$ |
| Pulse | ε, t₀, σ | $\psi(t) \propto \exp(-(t - t_0)^2 / 2\sigma^2)$ — symmetric around peak |

**Decision**: Use **onset + decay** profile. This captures: (1) zero spillover before migration arrives (t₀), (2) peak at onset, (3) gradual decline as birds move on (δ). The first confirmed case (22 Dec) with back-calculated infection ~12-14 Dec suggests t₀ ≈ early December.

**Prior for t₀**: Normal centred on 10–15 December 2025 (i.e., days 10–15 counting from 1 December), with SD ~5 days to allow data to inform onset timing.

**Identifiability note**: With binary HRZ, we estimate $\epsilon$ (peak spillover rate for HRZ farms). Non-HRZ spillover is assumed negligible and fixed at 0 for the initial model (narrative mentions "stragglers" but we treat this as second-order). The onset t₀ is identified from the timing of first infections; the decay δ is identified from the declining proportion of spillover-attributable cases over time.

**Spatial alternatives** (for later consideration):
- Binary HRZ + non-HRZ spillover: $\epsilon_{\text{HRZ}}$ and $\epsilon_{\text{non}}$ (narrative mentions "stragglers" outside HRZ)
- Continuous distance-to-wetland using `clc_32626.geojson`: $\epsilon \cdot \exp(-d_j / \rho)$

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

**Decision**: Use **exponential** kernel $K(d) = \exp(-d/\alpha)$. Parsimonious (single scale parameter $\alpha$), and the reparameterisation $\beta_0 = \beta \cdot K(d_0)$ is straightforward. Power-law and Cauchy kernels are retained as complexity options.

**Infectiousness profile**: Farm infectiousness increases over time as within-farm prevalence grows:
$$w(\tau) = 1 - \exp(-r \cdot \tau)$$
where $r \approx 1.0$/day is the within-farm growth rate (from mortality ledgers). This starts at 0 when $\tau=0$ and saturates to 1. Note: $w(0) = 0$ provides a "soft" de facto latency — a farm infected on day $t$ contributes negligible hazard on that day, ramping up continuously thereafter. This is distinct from a "hard" latent period (see Complexity Option 1).

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
- $\beta_{\text{duck}} \in (0, 1]$ (to be estimated)

**Provisional assumption for spillover**: Applying a single $\beta_{\text{duck}} \in (0, 1]$ uniformly to all pathways including $\lambda_j^{\text{spillover}}$ is provisional — ducks (waterfowl) might be *more* susceptible to wild-bird spillover than chickens, violating the (0, 1] constraint for that pathway. This simplification is adopted for the initial model; see Complexity Option 2 for pathway-specific modifiers or relaxed bounds.

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
- $T_{\text{capacity}}$ is the earliest time culling resources are available, modelled as a deterministic queue delay based on backlog (fixed offset from 6 Jan when capacity was reached)

**Note**: 77% of preventive cull dates are missing — impute as $T_{\text{trigger}} + \delta_{\text{prev}}$ where $\delta_{\text{prev}}$ is estimated from the 12 complete records.

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
  ε       = peak spillover hazard (HRZ farms)
  t₀      = spillover onset time (early December)
  δ       = spillover decay rate after onset
  β       = farm-to-farm transmission rate
  α       = spatial kernel scale
  β_duck  = relative susceptibility of ducks
              Prior: Beta(2,2) on (0,1] (weakly informative, symmetric)
              Fallback: if posterior is prior-dominated, treat as scenario parameter (fixed per run: 1.0, 0.5, 0.25)

Parameters (fixed):
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

    # Spillover (HRZ only, with temporal profile)
    ψ(t) = 0 if t < t₀, else exp(-δ(t - t₀))
    hazard_spillover = ε × ψ(t)  if j ∈ HRZ, else 0

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

  For each infected farm j:
    # Removal (culling) - modelled separately
    T_j^R = T_j^C + delay  (depends on observation)
```

---

## DAG Representation

```text
    PARAMETERS                      FIXED           COVARIATES                       DATA
    [ε, t₀, δ]  [β, α]  [β_duck]    {p_mov, r}      {HRZ}  {species}  {location}   {M_{i→j}(t): movements}
         |         |        |           |             |        |          |                  |
         v         v        v           v             v        v          v                  v
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

**Key distinction**: HRZ and species are covariates that modify how parameters act, not parameters themselves. We estimate ε (peak spillover rate), t₀ (spillover onset), δ (spillover decay), β (transmission rate), α (kernel scale), and β_duck (relative susceptibility). Movement transmission probability p_mov is fixed at 0.01 based on literature (mechanistically important but not identifiable from these data).

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
| ε (spillover) | Conditional | Requires non-HRZ infections; early outbreak confounded with β |
| t₀ (onset) | Conditional | Identified from timing of first infections; prior centred on early December |
| δ (decay) | Weak | Identified from declining proportion of spillover cases over time; jointly correlated with t₀ (see note below) |
| β (transmission) | Conditional | Identified from non-HRZ cases + spatial-temporal clustering |
| α (kernel scale) | Weak | Reparameterize: estimate β₀ = β·K(d₀) and α separately |
| β_duck (species) | **Weak** | Conflates susceptibility, infectiousness, and detectability; weakly informative prior Beta(2,2) on (0,1]. If posterior is prior-dominated, fall back to scenario analysis (fixed β_duck ∈ {1.0, 0.5, 0.25}) |
| r (growth rate) | Fixed | Set from mortality ledger data (= 1.0/day) |
| p_mov (movement) | Fixed | Set from literature (= 0.01); not identifiable, confounded with spatial kernel |

### t₀ vs δ (joint identifiability)

Under the onset+decay profile $\psi(t) = \exp(-\delta(t - t_0))$, the parameters t₀ and δ are inherently correlated: a later t₀ with slower δ can produce a similar spillover time-series to an earlier t₀ with faster δ, creating a ridge-shaped likelihood surface. With limited early-outbreak observations, this may cause poor MCMC mixing and inflated marginal uncertainties.

**Mitigations**:
- Reparameterise to $(t_0, \tau_{\text{half}} = \ln 2 / \delta)$, which is more interpretable and partially orthogonalises the parameterisation
- Profile the joint likelihood for $(t_0, \delta)$ before full inference to diagnose the ridge empirically
- Informative priors on t₀ (centred on early December) help anchor one end of the ridge

### ε vs β

Structurally identifiable because:
- ε: HRZ-only hazard with fixed spatial footprint (ε is constant amplitude; actual hazard = ε·ψ(t) where ψ(t) encodes temporal onset/decay via t₀ and δ)
- β: proximity-dependent hazard varying with distance to infected farms

**Key diagnostic**: Infections outside the HRZ can arise from local transmission (β) or movement from infected farms (including HRZ sources). Since movement transmission uses fixed $p_{\text{mov}}$, non-HRZ infections primarily identify β, and ε is then identified from the excess HRZ risk.

**Conditions for identification**:
1. Some infections outside HRZ
2. Some HRZ infections when local infectious pressure is low
3. Sufficient temporal variation

**Practical concerns**:
- Early outbreak: few infected farms means low local pressure everywhere, so early infections inform ε + β jointly rather than separately
- Culling delays: capacity constraints extend infectious periods, inflating apparent β
- Undetected farms: if surveillance misses cases, apparent spillover may be local transmission from undetected sources

Simulation-based identifiability checks recommended before drawing conclusions.

### α vs β

Kernel scale and transmission rate trade off: higher β with smaller α can produce similar patterns to lower β with larger α.

**Resolution**: Reparameterize to a compound parameter that is well-identified:
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
3. **Onset + decay spillover profile** — $\psi(t) = 0$ for $t < t_0$, then $\exp(-\delta(t - t_0))$; captures migration arrival and departure but not complex seasonal patterns
4. **Binary HRZ** — spillover risk is HRZ vs non-HRZ (continuous distance-to-water preferred if feasible)
5. **No explicit latent state** — no separate S→E→I transition; the infectiousness profile $w(\tau)$ provides soft latency ($w(0)=0$, ramping up continuously)
6. **Fixed within-farm growth rate** — $r = 1.0$/day from mortality ledgers
7. **Species effect constrained to (0,1]** — assumes ducks less susceptible; literature is mixed

---

## Complexity Options (for refinement)

If the simplified model is insufficient:

**High priority** (recommended before fitting):
1. **Hard latent period** — add explicit delay before farm becomes infectious: $w(\tau) = 0$ for $\tau < \tau_{\min}$ (e.g., $\tau_{\min} = 1$–2 days). This supplements the soft latency from $w(0)=0$ with a hard floor.
2. **Bidirectional species effect** — allow $\beta_{\text{duck}} \in (0, 2]$ or use log-scale

**Medium priority**:
3. **Continuous spillover risk** — distance to wetlands instead of binary HRZ (using `clc_32626.geojson`)
4. **Culling capacity model** — explicitly model queue and delays as function of outbreak intensity
5. **Estimate $p_{\text{mov}}$** — if movement transmission appears important, consider estimating rather than fixing

**Lower priority**:
6. **Shared services** — feed trucks, veterinary personnel, equipment (if data available)
7. **Farm-level heterogeneity** — random effects on susceptibility/infectiousness
8. **Wind-borne transmission** — directional kernel if spatial patterns suggest it
9. **More flexible spillover profile** — e.g., Gaussian pulse or spline if onset+decay is too restrictive

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

