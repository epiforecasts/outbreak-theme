# Data Source Selection

This document selects which data sources to use for each model component, based on the data characterisation and process DAG.

---

## Selection Criteria

For each data source, consider:
1. **Information content** — what parameters/states does it inform?
2. **Quality** — completeness, accuracy, bias
3. **Complexity cost** — what observation model is required?
4. **Redundancy** — do multiple sources inform the same quantity?

---

## Data Sources and Their Use

### Primary Data Sources (Required)

| Data Source | Use | Informs |
|-------------|-----|---------|
| `cases.csv` | Fit transmission model | Infection times (latent), spatial pattern, species distribution |
| `population.csv` | Define susceptible population | Farm locations, species, production type |
| `activity.csv` | Time-varying susceptibility | Which farms are at risk on each day |
| `hrz_32626.geojson` | Spillover risk | Binary HRZ indicator for each farm |

**Rationale**: These four sources are essential. Cases are the primary outcome; population defines the at-risk set; activity determines time-varying exposure; HRZ informs spillover.

### Secondary Data Sources (Conditional Use)

| Data Source | Use | Decision |
|-------------|-----|----------|
| `prev_culls.csv` | Track removed farms | **Include** — needed for susceptible population, despite 77% missing dates |
| `movement.csv` | Movement transmission | **Include (fixed parameter)** — mechanistically important; 6/103 cases detected via pre-shipment testing indicate movement is a real pathway |
| `clc_32626.geojson` | Continuous spillover risk | **Exclude initially** — binary HRZ is simpler; revisit if HRZ insufficient |

### Fixed Inputs (Not Fitted)

| Data Source | Use | Decision |
|-------------|-----|----------|
| Mortality ledgers | Within-farm growth rate | **Fix r ≈ 1.0/day** — n=3 insufficient to estimate, but consistent pattern |

---

## Detailed Selection Decisions

### 1. Cases (`cases.csv`) — PRIMARY

**Fields used**:
- `farm_id` — link to population
- `date_confirmed` — primary time variable (observation time)
- `date_suspicious` — informative for detection delay (99/103 have this)
- `cull_start`, `cull_end` — removal timing

**Fields not directly used**:
- `detection_method` — could stratify detection model, but adds complexity

**Observation model needed**: Map latent infection times → confirmation times

### 2. Population (`population.csv`) — PRIMARY

**Fields used**:
- `farm_id` — link to cases
- `x`, `y` — spatial coordinates for kernel
- `species` — chicken/duck for β_duck
- `production_type` — for stratified predictions (Q2)
- `capacity` — could weight infectiousness (not in initial model)

**Derived quantities**:
- Distance matrix between all farms
- HRZ membership (spatial join with `hrz_32626.geojson`)

### 3. Activity (`activity.csv`) — PRIMARY

**Fields used**:
- `farm_id` — link to population
- `date_start`, `date_end` — activity period
- `volume` — batch size (not used initially)

**Use**: Farm $j$ is only susceptible on day $t$ if active: $A_j(t) = 1$

**Missing data**: `date_end` missing for currently active farms — treat as ongoing

### 4. High-Risk Zone (`hrz_32626.geojson`) — PRIMARY

**Use**: Spatial join to assign HRZ indicator to each farm

**Alternative considered**: Use `clc_32626.geojson` for continuous distance-to-wetland. Decision: start with binary HRZ; less complex observation model.

### 5. Preventive Culls (`prev_culls.csv`) — SECONDARY

**Problem**: 77% of records (40/52) missing cull dates

**Options**:
1. Exclude — assume preventive culls don't affect transmission (unrealistic)
2. Impute dates based on trigger case timing + delay
3. Model culling queue explicitly

**Decision**: Option 2 — impute missing dates as $T_{\text{trigger}} + \delta_{\text{prev}}$, where $T_{\text{trigger}}$ is the confirmation date of the trigger case (confirmed farm within 1 km) and $\delta_{\text{prev}}$ is the median trigger-to-cull delay estimated from the 12 complete records.

### 6. Movements (`movement.csv`) — INCLUDED (fixed parameter)

**Rationale for inclusion**:
- 6/103 cases detected via pre-shipment testing — these are farms *intercepted before transmitting*, indicating movement is a real transmission pathway
- Movement provides a mechanistically distinct pathway from spatial proximity (can explain long-range jumps)
- @Yoo2021 found ~30% of Korean H5N6 transmission via vehicle movements
- Including with fixed $p_{\text{mov}}$ avoids misattributing movement transmission to the spatial kernel

**Fields used**:
- `source_farm_id`, `destination_farm_id` — define transmission pathway
- `date` — timing of movement event
- HRZ status of source farm — determines whether pre-shipment testing applies

**Fixed parameter**: $p_{\text{mov}} = 0.01$ per movement. Not estimated because:
- Confounded with spatial kernel for nearby farm pairs
- Pre-shipment testing sensitivity unknown
- Limited number of observed movement-related infections

**Modifiers**:
- Pre-shipment testing in HRZ intercepts infectious farms with probability $\sigma_{\text{test}} = 0.9$
- Movements from farms in regulated zones are blocked

### 7. Mortality Ledgers — FIXED PARAMETER

**Use**: Fix within-farm growth rate $r = 1.0$/day

**Rationale**:
- n=3 insufficient to estimate with uncertainty
- Consistent pattern across 3 farms (2 chicken, 1 duck)
- Affects infectiousness profile $w(\tau) = 1 - \exp(-r \cdot (\tau - \tau_{\min}))$ for $\tau \geq \tau_{\min}$

---

## Summary: Data Flow

```text
┌─────────────────────────────────────────────────────────────────┐
│                        PRIMARY DATA                              │
├─────────────────────────────────────────────────────────────────┤
│  population.csv  ──┬──> Farm locations, species                 │
│                    │                                             │
│  activity.csv    ──┼──> Time-varying susceptibility             │
│                    │                                             │
│  hrz_32626.geojson─┼──> Spillover risk (HRZ indicator)          │
│                    │                                             │
│  cases.csv       ──┴──> Observed confirmations (fit target)     │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      SECONDARY DATA                              │
├─────────────────────────────────────────────────────────────────┤
│  prev_culls.csv  ────> Removed farms (impute missing dates)     │
│  movement.csv    ────> Movement transmission (fixed p_mov)      │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      FIXED INPUTS                                │
├─────────────────────────────────────────────────────────────────┤
│  mortality_ledgers ──> r = 1.0/day (within-farm growth)         │
│  movement.csv      ──> p_mov = 0.01 (movement transmission)     │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                      EXCLUDED (initially)                        │
├─────────────────────────────────────────────────────────────────┤
│  clc_32626.geojson ──> Continuous spillover (revisit later)     │
└─────────────────────────────────────────────────────────────────┘
```

---

## Pre-processing Required

Before model fitting:

1. **Spatial join**: Assign HRZ indicator to each farm using `hrz_32626.geojson`
2. **Distance matrix**: Compute pairwise distances between all 9,160 farms
3. **Activity lookup**: Build function to query $A_j(t)$ for any farm/day
4. **Impute prev_cull dates**: For 40 farms with missing dates:
   - Identify trigger case (confirmed farm within 1 km)
   - Assign cull date = $T_{\text{trigger}} + \delta_{\text{prev}}$ (where $\delta_{\text{prev}}$ is the median delay from 12 complete records)
   - **Fallback** (no confirmed farm within 1 km): use the nearest confirmed case regardless of distance; if multiple equidistant, use the earliest confirmation date
5. **Merge removals**: Combine reactive culls (from `cases.csv`) and preventive culls into single removal timeline
6. **Movement network**: Build lookup for movements by date/source/destination
   - Filter to broiler_1 → broiler_2 movements only (these are the only live-bird transfers between production farms in the data; movements to slaughter do not create new farm infections)
   - Assign HRZ status to source farms for pre-shipment testing modifier
   - Identify movements blocked by regulated zones

---

## Parameters to Estimate vs Fix

| Parameter | Estimate/Fix | Rationale |
|-----------|--------------|-----------|
| φ (HRZ spillover rate) | Estimate | Key process parameter; LogNormal(log(0.01), 2) prior |
| η (non-HRZ reduction) | Estimate | φ_non = φ × η; Beta(1, 5) prior |
| t_change (spillover changepoint) | Estimate | Identified from changing HRZ/non-HRZ case ratio; Normal(31, SD=10) prior (t=0 ≡ 1 Dec 2025) |
| ρ (post-change intensity) | Estimate | Identified from HRZ excess in second period; Beta(2, 5) prior |
| β₀ (transmission at d₀) | Estimate | Compound parameter β·K(d₀); well-identified |
| α (kernel scale) | Estimate | Decay scale; identified from spatial pattern |
| β_duck (species susceptibility) | Estimate | Q3 target (interpret with caution); LogNormal(0, 0.5) on (0,2] |
| r (within-farm growth) | Fix = 1.0/day | n=3 data, consistent pattern |
| τ_min (hard latent period) | Fix = 1 day | Within-flock spread before between-farm transmission |
| p_mov (movement transmission) | Fix = 0.01 | Confounded with spatial kernel; @Yoo2021 suggests ~30% via movement |
| σ_test (pre-shipment sensitivity) | Fix = 0.9 | Pre-shipment testing intercept probability |
| Detection delay distribution | Estimate | Infection → confirmation delay; informed by date_suspicious → date_confirmed interval. Deferred to observation model (step 04). |

**Note**: β is derived as β = β₀/K(d₀) where d₀ is a reference distance (e.g., median inter-farm distance). See `step02_process_dag.md` for full parameter specifications and priors.

---

## Link to Questions

| Question | Data sources required |
|----------|----------------------|
| Q1 | cases.csv, population.csv (descriptive only) |
| Q2 | All primary + prev_culls + movement.csv (full model) |
| Q3 | All primary + prev_culls + movement.csv (estimate β_duck) |
| Q4 | All primary + prev_culls + movement.csv (simulate modified culling) |
| Q5 | All primary + prev_culls + movement.csv (simulate modified delays) |

*Note*: Q2–Q5 additionally require all fixed inputs (mortality ledgers for r, movement.csv for p_mov, σ_test); these are omitted from rows for brevity.

---
