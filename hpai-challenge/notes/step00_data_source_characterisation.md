---
title: "Data Source Characterisation — HPAI Modelling Challenge"
bibliography: references.bib
---

This document characterises each data source available for the Jolly Island HPAI outbreak using a six-category checklist: metadata, scope, resolution, data quality, data utility, and practical considerations.

---

## Data Source 1: Confirmed Case Notifications (`cases.csv`)

### 1.1 Basic Metadata

| Attribute | Value |
|-----------|-------|
| **Classification** | Farm-level outbreak time series |
| **Study design** | Routine passive surveillance + targeted testing (pre-shipment, contact tracing) |
| **Description** | All laboratory-confirmed HPAI H5N0 outbreaks on Jolly Island poultry farms, with dates of suspicion, confirmation, and culling |
| **Primary purpose** | Disease control and notification (mandatory reporting) |

### 1.2 Scope

| Attribute | Value |
|-----------|-------|
| **Source population** | All poultry farms on Jolly Island (n=9,160) |
| **Target population** | Farms with confirmed HPAI infection |
| **Stratifications** | Species (chicken/duck), production type, detection method, HRZ status (derivable) |
| **Collection type** | Routine (passive surveillance) + triggered (contact tracing, pre-shipment testing) |
| **Potential triggers** | Farmer reports sudden mortality increase; pre-shipment testing in HRZ; contact tracing from stage 2 broiler cases |

### 1.3 Resolution

| Attribute | Value |
|-----------|-------|
| **Data aggregation** | Individual farm-level |
| **Temporal data** | Yes |
| **Collection frequency** | Continuous (event-based) |
| **Reporting frequency** | Daily |
| **Time period covered** | 22 Dec 2025 – 13 Jan 2026 (24 days, 103 cases: 81 chicken, 22 duck) |
| **Spatial data** | Yes (via linkage to population.csv) |
| **Spatial resolution** | Farm-level (exact coordinates in metres) |
| **Geographical coverage** | Complete (all confirmed cases on the island) |

### 1.4 Data Quality

#### Quality of case definitions
**Assessment**: High
**Reasoning**: Case definition is standardised and objective — RT-PCR confirmation at the national reference laboratory. The narrative states "results are usually provided in around 2 days" and the surveillance system is described as "efficient" with "fair and timely compensation", reducing incentives for non-reporting.

#### Test sensitivity and specificity
**Assessment**: High specificity, variable sensitivity
**Reasoning**: RT-PCR has very high specificity (near 100% — false positives extremely rare). Sensitivity depends on sampling: 20 birds sampled per suspicion, which should detect infection once within-farm prevalence is substantial. However, **sensitivity for subclinical infections is unknown** — the test is only applied when suspicion is raised.

#### Reporting delays
**Assessment**: ~2 days (suspicion → confirmation), but infection → suspicion delay is unobserved
**Reasoning**: Narrative states "results are usually provided in around 2 days". The critical unmeasured delay is infection → detectable mortality → farmer suspicion, estimated at 5.9–7.4 days (2016) to 9.8–14.8 days (2014) for chickens in Dutch H5N8 outbreaks [@Hobbelen2020], and potentially longer or never for ducks due to subclinical infection [@Elbers2021].

#### Missingness
**Assessment**: Low for core fields, moderate for ancillary
**Reasoning from data**:
- `date_suspicious`: Missing for 4/103 cases (non-passive detections: 6 preshipment + 1 contact_tracing, but only 4 missing — some non-passive cases also have suspicion dates)
- `cull_start`/`cull_end`: Some missing — culls in progress at data cutoff

#### Sources of bias

| Bias type | Assessment | Reasoning |
|-----------|------------|-----------|
| **Selection bias (species)** | **Possible — informed by literature** | Passive surveillance triggers on "sudden increase in mortality". Ducks may show lower mortality with HPAI than chickens [@Elbers2021; @MSDVetManual], which could mean duck infections are less likely to trigger suspicion. However, this is an assumption from external literature, not validated in this dataset. The 4× attack rate difference (chicken 1.83% vs duck 0.47%) is *consistent with* detection bias but also consistent with true susceptibility differences. |
| **Selection bias (detection method)** | Moderate | 96/103 cases detected by passive surveillance, 6 by pre-shipment testing, 1 by contact tracing. Pre-shipment testing only applies to HRZ broiler_1 farms, creating differential surveillance intensity. |
| **Left truncation** | Present | Infections before 22 Dec (first detection) are unobserved. The true epidemic start is unknown. |
| **Right censoring** | Present | Some culls incomplete at 13 Jan cutoff. |
| **Time-varying detection** | Likely | Narrative mentions "enhanced biosecurity and detectability in regulated zones" — farms in zones may report suspicions faster. Early outbreak detection delays may be longer (lower awareness) than later. |

#### Direction of overall bias
**Assessment**: Likely under-ascertainment; species difference uncertain
**Reasoning**: The data likely misses some infections that never reach the mortality threshold triggering suspicion. Whether this affects ducks more than chickens is suggested by domain knowledge but not validated here. Observed attack rates (chicken: 81/4,435 = 1.83%; duck: 22/4,725 = 0.47%) show 4× higher detection in chickens — this could reflect true susceptibility differences, differential detection, or both. The data alone cannot distinguish these explanations.

### 1.5 Data Utility

| Quantity | Informed? | Reasoning |
|----------|-----------|-----------|
| Infection incidence | Indirect | Must back-calculate from confirmation dates using assumed delays |
| Transmission rate | Indirect | Inferred from spatial-temporal patterns of cases |
| Species susceptibility | Direct but biased | Observed case counts by species, but confounded by detection differences |
| Spatial kernel | Indirect | Inter-farm distances between cases inform local spread, but confounded with spillover |
| Detection delays (D→C) | Direct | `date_suspicious` to `date_confirmed` observed for passive cases |
| Infection-to-detection delay | No | `date_suspicious` is NOT infection date; the delay from infection to detectable mortality is unobserved |

### 1.6 Practical Considerations

| Attribute | Assessment | Reasoning |
|-----------|------------|-----------|
| **Scalability** | Good | Event-based collection scales with outbreak |
| **Accessibility** | Full | Complete dataset provided |
| **Linkage** | Excellent | `farm_id` links to all other datasets |
| **Format** | Good | Clean CSV |

### 1.7 Key Gaps Identified

1. **Infection times are latent** — only confirmation (and sometimes suspicion) dates observed
2. **Species ascertainment difference unknown** — cannot distinguish "ducks less susceptible" from "ducks less detectable"; domain knowledge suggests detection bias is plausible but this is not validated
3. **No measure of surveillance intensity** — cannot model time-varying detection probability
4. **Early epidemic truncated** — infections before 22 Dec not captured

---

## Data Source 2: Farm Population Register (`population.csv`)

### 2.1 Basic Metadata

| Attribute | Value |
|-----------|-------|
| **Classification** | Static population census |
| **Study design** | Complete census (country-wide) |
| **Description** | Exhaustive register of all poultry farms on Jolly Island with location, species, production type, capacity |
| **Primary purpose** | Regulatory/administrative (farm registration) |

### 2.2 Scope

| Attribute | Value |
|-----------|-------|
| **Source population** | All registered poultry farms |
| **Target population** | All farms (complete) |
| **Stratifications** | Species (chicken/duck), production type (layer, broiler_1, broiler_2, organic, conventional), county, district |
| **Collection type** | One-off census (established when HPAI 2.3.4.4b emerged in Europe) |

### 2.3 Resolution

| Attribute | Value |
|-----------|-------|
| **Data aggregation** | Individual farm-level |
| **Temporal data** | No (snapshot) |
| **Spatial data** | Yes |
| **Spatial resolution** | Farm-level (UTM coordinates in metres, EPSG:32626) |
| **Geographical coverage** | Complete (9,160 farms across 14 counties, 556 districts) |
| **Species breakdown** | 4,435 chicken (48%), 4,725 duck (52%) |
| **Production types** | 1,248 broiler_1, 2,746 broiler_2, 441 layer (chicken); 3,691 conventional, 1,034 organic (duck) |

### 2.4 Data Quality

#### Completeness
**Assessment**: High (stated as exhaustive)
**Reasoning**: Documentation states "the provided list is considered as exhaustive and exact as possible" — census conducted in collaboration with poultry organisations when HPAI 2.3.4.4b emerged. Regulatory requirement for farm registration.

#### Accuracy of attributes
**Assessment**: High for categorical variables, uncertain for capacity
**Reasoning**: Species and production type are fundamental to farm operations and registration — unlikely to be misclassified. Capacity is "maximum number of birds that can be housed" — may not reflect actual flock size, which varies with activity periods.

#### Spatial accuracy
**Assessment**: High (to be verified)
**Reasoning**: UTM coordinates in metres suggest precise geocoding. Need to verify no duplicate coordinates or obvious errors.

#### Temporal validity
**Assessment**: Concern — static snapshot
**Reasoning**: Census is static. New farms established or farms closed since census would not be captured. For a short outbreak period (44 days), this is likely acceptable.

### 2.5 Data Utility

| Quantity | Informed? | Reasoning |
|----------|-----------|-----------|
| Susceptible population | Direct | Denominators for attack rates |
| Spatial structure | Direct | Farm locations enable distance calculations |
| Species/production stratification | Direct | Core attributes |
| HRZ classification | Indirect | Requires spatial join with HRZ polygon |
| Actual flock size | No | Capacity ≠ current flock size; need activity.csv for presence |

### 2.6 Practical Considerations

| Attribute | Assessment | Reasoning |
|-----------|------------|-----------|
| **Scalability** | N/A | Static census |
| **Accessibility** | Full | Complete dataset provided |
| **Linkage** | Excellent | `farm_id` links to all other datasets |
| **Format** | Good | Clean CSV with UTM coordinates |

### 2.7 Key Gaps Identified

1. **No temporal dimension** — cannot track farm establishment/closure
2. **Capacity vs actual population** — capacity is maximum, not current flock size
3. **No information on biosecurity practices** — all farms treated equally within production type

---

## Data Source 3: Farm Activity Periods (`activity.csv`)

### 3.1 Basic Metadata

| Attribute | Value |
|-----------|-------|
| **Classification** | Time-varying farm status |
| **Study design** | Routine administrative records |
| **Description** | Periods when birds are present on each farm (date_start to date_end); absence indicates downtime |
| **Primary purpose** | Farm management / regulatory compliance |

### 3.2 Scope

| Attribute | Value |
|-----------|-------|
| **Source population** | All poultry farms with recorded batch activity |
| **Target population** | Active farms (birds present) |
| **Temporal coverage** | Full outbreak period (batch start/end dates) |

### 3.3 Resolution

| Attribute | Value |
|-----------|-------|
| **Data aggregation** | Farm-level batches |
| **Temporal resolution** | Daily (date_start, date_end) |
| **Spatial resolution** | N/A (links to population.csv via farm_id) |

### 3.4 Data Quality

#### Completeness
**Assessment**: High
**Reasoning**: Movement registration is mandatory and automated ("Automated systems and a long-established trust between farmers and the authorities enabled good compliance"). Activity periods are derived from batch arrivals/departures which are administratively tracked.

#### Accuracy
**Assessment**: High for dates, volume may differ from actual
**Reasoning**: Dates tied to regulatory movement records. Volume is "at most the capacity of the farm" — actual mortality during the period would reduce this, but mortality is not tracked in activity data.

#### Missingness
**Assessment**: Low, by design
**Reasoning**: `date_end` missing for currently active farms is intentional — indicates ongoing activity at data cutoff.

### 3.5 Data Utility

| Quantity | Informed? | Reasoning |
|----------|-----------|-----------|
| Time-varying susceptibility | Direct | Farms only at risk when active (birds present) |
| Susceptible farm-days | Direct | Can compute exposure time per farm |
| Actual flock size | Partial | Volume gives batch size, but not adjusted for mortality |

### 3.6 Practical Considerations

| Attribute | Assessment | Reasoning |
|-----------|------------|-----------|
| **Scalability** | Good | Routine administrative records |
| **Accessibility** | Full | Complete dataset provided |
| **Linkage** | Excellent | `farm_id` links to population and cases |
| **Format** | Good | Clean CSV |

### 3.7 Key Gaps Identified

1. **No within-period mortality adjustment** — volume is starting batch size, not current population
2. **Downtime interpretation** — assumes empty farms are truly at zero risk (reasonable for cleaned facilities)

---

## Data Source 4: Animal Movement Records (`movement.csv`)

### 4.1 Basic Metadata

| Attribute | Value |
|-----------|-------|
| **Classification** | Movement network data |
| **Study design** | Mandatory movement notification |
| **Description** | All movements between stage 1 and stage 2 chicken broiler farms |
| **Primary purpose** | Regulatory (traceability) |

### 4.2 Scope — Critical Limitation

**Assessment**: Partial coverage
**Reasoning**: Movement data **only captures broiler_1 → broiler_2 transfers**. This excludes:
- Movements to slaughterhouses (all production types)
- Layer farm movements
- Duck farm movements (organic, conventional)
- Equipment/personnel movements (not recorded)

This means movement-mediated transmission can only be assessed for one specific pathway.

### 4.3 Resolution

| Attribute | Value |
|-----------|-------|
| **Data aggregation** | Individual movements |
| **Temporal resolution** | Daily (movement date) |
| **Spatial resolution** | Farm-to-farm (via farm_id linkage) |
| **Network size** | 7,187 movements between broiler farms |

### 4.4 Data Quality

#### Completeness (within scope)
**Assessment**: High
**Reasoning**: "Automated systems and a long-established trust between farmers and the authorities enabled a good compliance, making the list of movements exhaustive."

#### Selection bias
**Assessment**: High (by design)
**Reasoning**: Only broiler movements recorded. If movement transmission occurs via other pathways (shared equipment, personnel, vehicles visiting multiple farms), this is entirely unobserved.

### 4.5 Data Utility

| Quantity | Informed? | Reasoning |
|----------|-----------|-----------|
| Movement transmission (broiler pathway) | Direct | Can assess if infected broiler_1 farms sent birds to farms that became cases |
| Movement transmission (other pathways) | No | No data on equipment, personnel, or non-broiler movements |
| Network structure | Partial | Only broiler_1 → broiler_2 edges |

#### Empirical signal assessment
**Assessment**: Weak signal
**Reasoning**: Only 6/103 cases detected via pre-shipment testing. This could mean: (a) movement transmission is rare, (b) pre-shipment testing is effective at blocking it, or (c) most movement transmission occurs via unrecorded pathways. Cannot distinguish these explanations.

### 4.6 Practical Considerations

| Attribute | Assessment | Reasoning |
|-----------|------------|-----------|
| **Scalability** | Good | Routine administrative records |
| **Accessibility** | Full | Complete dataset provided |
| **Linkage** | Excellent | `source_farm_id` and `destination_farm_id` link to population |
| **Format** | Good | Clean CSV |

### 4.7 Key Gaps Identified

1. **Partial coverage** — only broiler movements, not other farm types or non-animal contacts
2. **Cannot separate transmission routes** — movement vs local spread vs spillover confounded
3. **Pre-shipment testing confounds interpretation** — blocked movements not observed

---

## Data Source 5: Preventive Culling Records (`prev_culls.csv`)

### 5.1 Basic Metadata

| Attribute | Value |
|-----------|-------|
| **Classification** | Intervention records |
| **Study design** | Operational records from outbreak response |
| **Description** | Farms culled preventively (within 1km of confirmed case, from 1 Jan) |
| **Primary purpose** | Disease control operations |

### 5.2 Scope

| Attribute | Value |
|-----------|-------|
| **Source population** | Farms within 1km of confirmed cases |
| **Target population** | All susceptible farms in proximity zone (after 1 Jan policy) |
| **Temporal coverage** | 1 Jan onwards (policy implementation date) |
| **Total records** | 52 farms |

### 5.3 Resolution

| Attribute | Value |
|-----------|-------|
| **Data aggregation** | Farm-level |
| **Temporal resolution** | Daily (when recorded) |
| **Key limitation** | 77% (40/52) of cull dates missing |

### 5.4 Data Quality

#### Completeness
**Assessment**: Low for temporal data
**Reasoning**: Documentation and empirical check needed — narrative mentions "52 entries, 12 with recorded cull dates". This means **77% of preventive culls lack timing information**.

#### Why dates are missing
**Assessment**: Unknown — possibly operational
**Reasoning**: Could be: (a) culls planned but not yet executed at data cutoff, (b) recording gaps during operational pressure, (c) data extraction issue. The narrative mentions "culling capacity has been reached since 6 January" — backlogs may explain missing dates.

### 5.5 Data Utility

| Quantity | Informed? | Reasoning |
|----------|-----------|-----------|
| Which farms were preventively culled | Direct | Farm IDs recorded |
| When they were removed from susceptible pool | Partial | Only 12/52 have dates |
| Infection status of culled farms | Indirect | RT-PCR applied to all; positives appear in cases.csv |

### 5.6 Practical Considerations

| Attribute | Assessment | Reasoning |
|-----------|------------|-----------|
| **Scalability** | N/A | Outbreak-specific |
| **Accessibility** | Full | Complete dataset provided |
| **Linkage** | Good | `farm_id` links to population |
| **Format** | Good | Clean CSV |

### 5.7 Key Gaps Identified

1. **Missing cull dates (77%)** — severely limits ability to model time-varying susceptible population
2. **No information on reason for delay** — cannot model culling capacity constraints

---

## Data Source 6: Mortality Ledgers (3 example farms)

### 6.1 Basic Metadata

| Attribute | Value |
|-----------|-------|
| **Classification** | Within-farm surveillance data |
| **Study design** | Convenience sample of routine farm records |
| **Description** | Daily mortality counts from 3 infected farms |
| **Primary purpose** | Farm management / regulatory compliance |

### 6.2 Scope

| Attribute | Value |
|-----------|-------|
| **Source population** | 3 of 103 confirmed outbreak farms |
| **Target population** | Representative infected farms (uncertain) |
| **Temporal coverage** | Days leading up to detection for each farm |
| **Selection criteria** | Unknown — may not be representative |

### 6.3 Resolution

| Attribute | Value |
|-----------|-------|
| **Data aggregation** | Farm-level daily counts |
| **Temporal resolution** | Daily |
| **Sample size** | n=3 farms (2 chicken broiler_2, 1 duck conventional) |

### 6.4 Data Quality

#### Internal quality
**Assessment**: High for recorded data
**Reasoning**: Official farm records, mandated to be kept "in pristine condition". Daily mortality counts with farm ID, capacity, date.

#### External validity (representativeness)
**Assessment**: **LOW — major concern**
**Reasoning**: Only 3 of 103 outbreak farms. Selection criteria unknown — were these chosen because they had good records? Because they showed "typical" dynamics? Because they were accessible? The sample may not be representative.

#### Sample size for inference
**Assessment**: Adequate for growth rate estimation; limited for heterogeneity
**Reasoning**: n=3 farms but each provides a detailed daily time series (multiple data points per farm). The exponential growth pattern is consistent across all 3 farms, suggesting the within-farm growth rate r ≈ 1.0/day is robust. This is sufficient to estimate a single growth rate parameter with reasonable confidence, but insufficient to characterise between-farm variability or production-type-specific differences.

### 6.5 Observed Pattern

The mortality ledgers show consistent exponential growth in daily mortality across all 3 farms, with estimated growth rate r ≈ 1.0/day (doubling time ~0.7 days). This consistency across farms suggests the growth rate is driven by HPAI biology rather than farm-specific factors.

**Validation against literature**: Published within-flock R₀ estimates vary widely — @Savill2006 report R₀ ranging from 2.5 (H5N1 in vaccinated flocks) to 208 (H7N7 in unvaccinated flocks). Converting R₀ to exponential growth rate r requires assumptions about generation time: r = (R₀ - 1) / T_g. With T_g ≈ 1–2 days for HPAI, R₀ = 2.5 gives r ≈ 0.75–1.5/day, which brackets our estimate of r ≈ 1.0/day. The consistency suggests our estimate is biologically plausible, though the small sample (n=3) limits confidence.

### 6.6 Species and Production Type Coverage

**Verified from data**: The 3 mortality ledger farms are:
- Farm 2395: broiler_2, chicken
- Farm 3013: broiler_2, chicken
- Farm 8120: conventional, duck

This provides representation of both species but limited production type coverage (2 broiler_2, 1 conventional; no layer, broiler_1, or organic).

### 6.7 Practical Considerations

| Attribute | Assessment | Reasoning |
|-----------|------------|-----------|
| **Scalability** | Limited | Manual records, small sample |
| **Accessibility** | Partial | Only 3 farms provided |
| **Linkage** | Good | `farm_id` links to population |
| **Format** | Good | Clean daily records |

### 6.8 Key Gaps Identified

1. **Cannot characterise heterogeneity** — n=3 sufficient for point estimate but not for variance across farms
2. **Limited production type coverage** — 2 broiler_2 chicken, 1 conventional duck; missing layer, broiler_1, organic
3. **Selection criteria unknown** — were these farms chosen because they had good records, or are they representative?

---

## Data Source 7: High-Risk Zone Boundary (`hrz_32626.geojson`)

### 7.1 Basic Metadata

| Attribute | Value |
|-----------|-------|
| **Classification** | Regulatory spatial boundary |
| **Study design** | Administrative definition |
| **Description** | Polygon defining areas with elevated wild bird contact risk |
| **Primary purpose** | Pre-shipment testing requirement (regulatory) |

### 7.2 Scope

| Attribute | Value |
|-----------|-------|
| **Source population** | N/A (boundary, not population) |
| **Spatial coverage** | Areas of wild bird stopover / high gull concentration |
| **Farms in HRZ** | 1,962 of 9,160 (21%) |

### 7.3 Resolution

| Attribute | Value |
|-----------|-------|
| **Data aggregation** | Single polygon boundary |
| **Spatial resolution** | Precise boundary (EPSG:32626) |
| **Temporal resolution** | Static (no temporal dimension) |

### 7.4 Data Quality

**Assessment**: High (definitional)
**Reasoning**: This is a regulatory boundary, not a measured quantity. It defines where pre-shipment testing is required and where spillover risk is administratively considered elevated. The boundary itself is not subject to measurement error.

### 7.5 Validity as Spillover Proxy

**Assessment**: Uncertain
**Reasoning**: The HRZ is defined as areas "where migratory birds stop over during their migration... or where high concentrations of Laridae (gulls) are found." This is based on ecological knowledge, but:
- No wild bird surveillance data to validate
- Spillover risk may vary within HRZ (distance to water bodies, etc.)
- Non-HRZ areas still have some spillover risk ("stragglers seen through February")

### 7.6 Data Utility

| Quantity | Informed? | Reasoning |
|----------|-----------|-----------|
| Binary spillover risk | Direct | HRZ vs non-HRZ classification |
| Continuous spillover risk | No | Need distance-to-wetland from CLC data |
| Temporal spillover pattern | No | Static boundary |

### 7.7 Practical Considerations

| Attribute | Assessment | Reasoning |
|-----------|------------|-----------|
| **Scalability** | N/A | Static boundary |
| **Accessibility** | Full | Complete dataset provided |
| **Linkage** | Via spatial join | Point-in-polygon with farm locations |
| **Format** | Good | GeoJSON (EPSG:32626) |

### 7.8 Key Gaps Identified

1. **Binary classification may be too coarse** — spillover risk likely varies continuously
2. **No temporal dimension** — migration timing not encoded in spatial boundary
3. **No wild bird data to validate** — HRZ is an assumption, not observation

---

## Data Source 8: Supporting Spatial Data

### 8.1 Land Cover (`clc_32626.geojson`)

| Attribute | Value |
|-----------|-------|
| **Classification** | CORINE Land Cover (CLC) — European standard |
| **Description** | Wetlands and water bodies across Jolly Island |
| **Contents** | 647 polygon features |
| **Categories** | Wetlands (4xx): marshes, bogs, salt marshes, intertidal; Water (5xx): lakes, lagoons, estuaries |
| **Data quality** | High — standardised classification, precise boundaries (EPSG:32626) |

**Potential utility**: Could inform continuous spillover risk model based on distance to wetlands/water, rather than binary HRZ classification.

**Current status**: Excluded from initial model in favour of binary HRZ classification (simpler to implement and interpret). Revisit if binary HRZ proves insufficient.

### 8.2 Administrative Boundaries

| File | Description | Use |
|------|-------------|-----|
| `counties_32626.geojson` | 14 NUTS-3 level boundaries | Aggregation for reporting |
| `districts_32626.geojson` | 556 LAU level boundaries | Aggregation for reporting |

**Utility**: County/district-level summaries for descriptive epidemiology (Q1), not transmission modelling.

### 8.3 Practical Considerations (all supporting spatial data)

| Attribute | Assessment | Reasoning |
|-----------|------------|-----------|
| **Scalability** | N/A | Static boundaries |
| **Accessibility** | Full | Complete datasets provided |
| **Linkage** | Via spatial operations | Point-in-polygon, distance calculations |
| **Format** | Good | GeoJSON (EPSG:32626) |

---

## Summary: Quality Assessment Reasoning

| Data Source | Quality | Key Reasoning |
|-------------|---------|---------------|
| cases.csv | High for confirmed cases; species bias uncertain | RT-PCR confirmation is definitive; passive surveillance may miss subclinical infections (domain knowledge suggests ducks more affected, but not validated) |
| population.csv | High | Complete census, regulatory requirement, unlikely misclassification |
| activity.csv | High | Tied to mandatory movement records |
| movement.csv | High within scope; partial coverage | Only broiler movements captured; other transmission routes unobserved |
| prev_culls.csv | **Low for timing** | 77% of cull dates missing |
| mortality ledgers | High internally; limited external validity | n=3 with detailed time series — adequate for growth rate, limited for heterogeneity |
| hrz_32626.geojson | High (definitional) | Regulatory boundary, not measured |
| clc_32626.geojson | High | CORINE land cover — 647 wetland/water features; potential for continuous spillover risk |

---

## Phase 2 iteration

Phase 2 data were released with a later cut-off (13 Feb 2026 vs 13 Jan 2026). No structural changes to CSV schemas. The key differences are in volume and coverage.

### Updated data volumes

| Data source | Phase 1 | Phase 2 | Notes |
|-------------|---------|---------|-------|
| `cases.csv` | 103 cases | 466 cases | New detection method `post-culling` (61 cases) — farms found positive during preventive culling |
| `prev_culls.csv` | 52 culls | 539 culls | Now heavily chicken-focused (506/539) |
| `activity.csv` | 23,969 records | 27,506 records | Extended to 13 Feb |
| `movement.csv` | 7,187 records | 7,959 records | Extended to 13 Feb |

### Post-culling detection

The new `post-culling` detection method (61/466 cases) represents farms found positive during preventive culling operations. For these farms, `cull_start` precedes `date_confirmed`. The back-calculation still works: $T_j^I = T_j^C - 11$ gives an infection date before the cull start.

### Preventive culling coverage

The jump from 52 to 539 preventive culls reflects both the expanded observation window and the phase 2 policy shift to chicken-only culling within 3 km (from all-species within 1 km). Missing cull date rates improved — most phase 2 records have timing data.

### Cases with incomplete culls

14 cases have `cull_status` of "planned" or "in progress" at the data cut-off. These have `cull_start` but no `cull_end`. We only use `cull_start` for removal timing, so this does not affect the model.

---

## Key Data Gaps and Modelling Implications

### Gap 1: Infection times are latent
**What's missing**: True infection dates for cases
**Why it matters**: Any transmission analysis requires infection times to understand who-infected-whom
**Implication**: Must back-calculate using assumed delay distributions (infection → detection → confirmation)
**Domain knowledge input**: Literature suggests 5-15 days from introduction to detection for chickens

### Gap 2: Species difference in ascertainment (uncertain)
**What's missing**: True infection status of duck farms; validation of detection differences
**Why it matters**: Cannot distinguish "ducks less susceptible" from "ducks less detectable"
**Implication**: Any species effect estimate will conflate susceptibility and detectability
**Domain knowledge input**: Elbers et al. (2021) suggests lower mortality in duck HPAI infections, but this is external literature, not validated in this outbreak
**Epistemic status**: Hypothesis informed by domain knowledge, not established fact

### Gap 3: Spillover not directly observed
**What's missing**: Wild bird infection data, spillover event times
**Why it matters**: Cannot separate spillover from local transmission
**Implication**: Must infer spillover from spatial pattern (HRZ vs non-HRZ cases)
**Domain knowledge input**: Migration timing ("early winter with stragglers through February")

### Gap 4: Within-farm heterogeneity unknown
**What's missing**: Mortality data from broader sample of farms across production types
**Why it matters**: Farm-level infectiousness profile may vary by production type
**Implication**: Can estimate mean growth rate r ≈ 1.0/day from n=3 sample (consistent across 2 chicken broiler_2, 1 duck conventional), but cannot characterise variance or production-type differences

### Gap 5: Incomplete intervention timing
**What's missing**: Cull dates for 77% of preventive culls
**Why it matters**: Affects susceptible population calculation
**Implication**: Must impute or make assumptions about removal timing

### Gap 6: Movement data partial
**What's missing**: Non-broiler movements, equipment/personnel contacts
**Why it matters**: Movement transmission may occur via unrecorded pathways
**Implication**: Movement transmission effect may be underestimated or confounded with proximity

---

## Phase 3 iteration

Phase 3 data cover the full resolved epidemic (last case 7 April 2026). No structural changes to CSV schemas.

### Updated data volumes

| Data source | Phase 2 | Phase 3 | Notes |
|-------------|---------|---------|-------|
| `cases.csv` | 466 cases | 560 cases (+94) | Last case 7 Apr. No new detection methods |
| `prev_culls.csv` | 539 (21 incomplete) | 691 (all completed) | +152 culls, all finished |
| `activity.csv` | 27,506 records | 33,420 records | Extended to mid-April |
| `movement.csv` | 7,959 records | 9,341 records | Extended to 11 Apr |

### Preventive cull completion

All 691 preventive culls now have recorded completion dates. The right-censoring gap identified in phase 1 (77% missing cull dates) and residual incompleteness from phase 2 are fully resolved. The imputation approach used in earlier phases is no longer needed for retrospective analysis.

### No schema changes

All six data sources retain their phase 1/2 structure. The quality issues identified for earlier phases (species detection bias, movement data partial coverage, n=3 mortality ledgers) remain unchanged — additional data were not collected to address them.

---

