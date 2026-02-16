---
title: "Domain Knowledge: HPAI Biology and Farm-Level Transmission"
bibliography: references.bib
---

This document compiles domain knowledge from HPAI literature and the Jolly Island
narrative to inform the process DAG. It is independent of any model fitting results.

---

## 1. Virus biology (individual bird level)

### Latent and infectious periods
- **Latent period** (infection → virus shedding): shedding in faeces and respiratory
  secretions as early as 1–2 days post-inoculation in experimentally infected birds
  [@Spickler2008]
- **Infectious period**: mean 4.2 days (95% CI 3.3–5.5) estimated for H5N2 in
  turkey flocks, with a latent period of 1.4 days (95% CI 0.9–2.6)
  [@Ssematimba2019]
- **Incubation period** (infection → clinical signs): "highly variable, ranging from
  a few days in individual birds to 2 weeks in a flock" [@MSDVetManual]. The WOAH
  Terrestrial Code specifies 14 days as the flock-level incubation period for control
  programmes [@WOAHCode].

### Species differences
- **Chickens/turkeys (galliformes)**: highly susceptible; "HPAI viruses cause severe,
  systemic disease with high mortality rates in chickens, turkeys, and other gallinaceous
  poultry; mortality rates can be as high as 100% in a few days after exposure in
  unvaccinated birds" [@MSDVetManual].
- **Ducks (palmipedes)**: "wild waterfowl and shorebirds are often subclinically
  affected carriers" [@MSDVetManual]. Detection is much harder in domestic duck
  flocks — "mild clinical signs and low mortality levels induced by HPAIV infection
  in farmed ducks pose an increased risk of late diagnosis or misdiagnosis, leading
  to delayed notification" [@Elbers2021]. A daily mortality threshold >0.3% (after
  the first week of production) can detect infection with sensitivity and specificity
  >80%, and drops in feed or water intake may signal infection 1–2 days before
  mortality does [@Elbers2021].
- **Implication**: ducks have lower *apparent* susceptibility (fewer detected cases)
  but this conflates true susceptibility with detectability. A duck farm may be infected
  but never detected because mortality never triggers suspicion.

### Virus shedding and environmental persistence
- Faecal shedding is a major route; respiratory less important for farm-to-farm
- High virus titres in litter, water, on surfaces
- Cold temperatures and wet conditions favour environmental survival of AIV

---

## 2. Farm-level dynamics

### Infection to detection timeline
The within-farm epidemic unfolds as a chain of events:

```text
Virus introduction → initial bird infections → exponential within-farm spread
→ detectable mortality increase → farmer raises suspicion → sampling → lab confirmation
→ culling
```

**Key time intervals (from literature):**

| Interval | Estimate | Source |
|----------|----------|-------|
| Introduction → first shedding | 1–2 days | [@Spickler2008] |
| Introduction → detectable mortality in chickens | 5.9–14.8 days | Dutch H5N8: 5.9–7.4d (2016), 9.8–14.8d (2014) [@Hobbelen2020] |
| Introduction → detection in ducks | Longer, highly variable | May never reach suspicion threshold [@Elbers2021] |
| Detection → lab confirmation | ~2 days | Narrative |
| Confirmation → culling | Days to >1 week | Capacity-dependent; narrative |

Note: the farm-level incubation period (introduction → detection) was estimated
in only 4 of 24 between-farm transmission modelling studies identified in a
systematic review [@Guinat2023], reflecting the difficulty of observing this
quantity. Estimates are sensitive to assumptions about the bird-level latent
period [@Ssematimba2019].

**Netherlands 2003 (H7N7):** farm infectious period estimated at 7.47 days
(95% CI 7.2–7.8), Gamma-distributed with shape 11.1 [@Boender2007].

**Minnesota 2015 (H5N2):** assumed infection 8 days before detection, 3-day latent
period, then infectious until disposal. Average infectious period 17.2 days — much
longer due to slower US depopulation timelines [@Bonney2018].

### Within-farm infectiousness growth
- Mortality ledgers from the Jolly Island narrative show consistent exponential growth
  in daily mortality across 3 example farms
- Doubling time ~0.6–0.7 days (growth rate r ≈ 1.0/day)
- This growth reflects the within-farm epidemic: more infected birds = more virus
  output = higher between-farm transmission risk
- Farm-level infectiousness should increase over time until culling truncates it

---

## 3. Between-farm transmission

### Spatial kernel (farm-to-farm)

The standard kernel form in HPAI literature is a **logistic (power-law-like)** function
[@Boender2007]:

```text
h(r) = h₀ / (1 + (r/r₀)^α)
```

**Estimated parameters from previous outbreaks:**

| Outbreak | h₀ (day⁻¹) | r₀ (km) | α | Source |
|----------|-----------|---------|---|--------|
| Netherlands 2003 (H7N7) | 0.0020 | 5.6 | 3.38 | [@Boender2007] |
| Minnesota 2015 (H5N2) | 0.0061 | 7.02 | 2.46 | [@Bonney2018] |

Key features:

- h₀ ≈ 10⁻³ to 10⁻² per day (maximum hazard at zero distance)
- r₀ ≈ 5–7 km (half-maximum distance)
- α ≈ 2.5–3.5 (steepness of decay; α=2 is Cauchy-like, α→∞ is step function)
- Majority of transmission within 25 km
- Minnesota had longer-range transmission than Netherlands

**Important:** Four kernel formulations were compared for Minnesota; AIC differences
were <2, so kernel form choice may not be critical [@Bonney2018]. However, the logistic
form is the standard in the HPAI literature.

### Transmission mechanisms
- Shared equipment, vehicles, personnel between nearby farms
- Aerosol/dust (short range, <1 km)
- Fomites on clothing, vehicles
- Wind-borne transmission from contaminated farm dust, estimated to explain ~24% of
  transmission within 25 km in the Dutch 2003 epidemic [@Ssematimba2012]
- The kernel captures ALL these mechanisms aggregated over distance

### Movement-mediated transmission
- Direct transport of birds between farms
- Can cause long-range jumps (not captured by distance kernel)
- Vehicle movements estimated to contribute ~30% of H5N6 inter-farm infections
  in the 2016–2017 Korean outbreak [@Yoo2021]
- Movements are blocked in regulated zones
- Pre-shipment testing can intercept (sensitivity depends on test and timing)

---

## 4. Wild bird spillover

### Seasonal pattern
- Migratory waterfowl (natural AIV reservoirs) move south in autumn, north in spring
- Global HPAI seasonality: typically lowest in September, rising from October,
  peaking around January–February in the Northern Hemisphere
- Cold and wet conditions favour environmental virus survival
- Spillover risk varies with stopover duration and local bird density

### Jolly Island narrative specifics
- "Wild birds migrate in early winter with stragglers seen through February"
- Virus believed introduced by wild birds, originating near the Southern region (HRZ)
- HRZ defined around areas of high wild bird concentration
- First case detected 22 Dec; back-calculated infection ~12–14 Dec
- This suggests migration/spillover onset in early-to-mid December

### Spatial heterogeneity
- HRZ (1,962 of 9,160 farms) has higher wild bird density/contact
- Non-HRZ farms still at some risk (stragglers, wider bird movements)
- Ratio: observed attack rate in HRZ >> non-HRZ (to be quantified from data)

### Temporal profile
- Not a constant force — builds up as birds arrive, peaks, then declines
- The precise shape is uncertain:
  - Could be gradual build-up over weeks
  - Could be relatively sudden arrival/departure
  - "Stragglers through February" suggests slow tail
- Key biological question: when did spillover pressure *start*? This constrains
  when the first infections could have occurred.

---

## 5. Detection and surveillance

### Passive surveillance (96/103 cases)
- Farmer notices "sudden increase in mortality" (narrative)
- Triggers: mortality rate exceeding baseline expectations
- More reliable in chickens (dramatic mortality) than ducks (subclinical)
  [@Elbers2021]
- Reporting delay: farmer → vet → sampling → lab

### Pre-shipment testing (6/103 cases)
- Required in HRZ before bird movements (narrative)
- Batch test of 20 birds, 3–7 days before shipment (narrative)
- Test sensitivity: not explicitly stated but batch testing of 20 birds from a
  flock with rising prevalence should have reasonable sensitivity once within-farm
  epidemic is established

### Contact tracing (1/103 cases)
- Backward tracing from confirmed case to source farms (narrative)
- Forward tracing from source to destination farms
- Only yielded 1 case — suggesting limited movement-mediated transmission

### Detection biases
- Ducks: under-detected (subclinical carriage, lower mortality) [@Elbers2021]
- Small farms: may have less monitoring
- Early epidemic: lower awareness, longer detection delays
- Later epidemic: heightened awareness, possibly faster detection
- Preventive culls remove farms that might have been infected but never confirmed

---

## 6. Interventions (from narrative)

### Reactive culling
- All birds on confirmed farm killed and disposed of
- Timeline: confirmation → culling start (data-driven, variable)
- Capacity constraints from 6 January (backlog)

### Regulated zones
- 10 km surveillance zones around confirmed farms
- Duration: 28 days after confirmation
- Effects: heightened biosecurity, movement restrictions
- Effectiveness: not directly measured; assumed to reduce hazard

### Preventive culling (from 1 January)
- All farms within 1 km of confirmed case, irrespective of species
- Removes susceptible farms from the at-risk pool
- Also removes potentially infected-but-undetected farms

### National standstill
- 3-day national standstill after index case (early outbreak)
- Blocks all movements temporarily

### Pre-shipment testing (HRZ)
- Testing required before movements from HRZ farms
- Infectious farms intercepted with some probability

---

## 7. Key uncertainties and modelling implications

1. **Duck detection**: are we missing duck infections? The "complete case ascertainment"
   assumption may be violated for ducks specifically [@Elbers2021].

2. **Spillover onset**: when did wild bird pressure begin? This is not directly
   observed and critically determines the epidemic's time zero.

3. **Kernel form**: exponential vs logistic/power-law — literature favours
   logistic [@Boender2007], but kernel form may not matter much
   (AIC evidence from [@Bonney2018]).

4. **Farm-level delays**: the chain infection → shedding → within-farm spread →
   mortality → detection involves multiple poorly-observed steps. Aggregating into
   "infection → detection delay" is pragmatic but hides biological complexity;
   only a minority of modelling studies estimate this directly [@Guinat2023].

5. **Spatial vs spillover attribution**: without genomic data, we cannot distinguish
   "this farm was infected by its neighbour" from "both farms received independent
   spillover". The model must partition observed cases between pathways using spatial
   and temporal patterns alone.

## References
