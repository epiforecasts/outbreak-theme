# Output file descriptions

All output CSV files are in the `output/` directory.

## Descriptive outputs (Q1)

| File | Description |
|------|-------------|
| `spatial_cases.csv` | All 103 confirmed cases with farm coordinates (UTM), county, district, species, production type, detection method, and dates |
| `epicurve_confirmed.csv` | Daily epicurve of confirmed cases (total, chicken, duck, cumulative) |
| `epicurve_by_production.csv` | Daily epicurve broken down by production type |
| `cases_by_county.csv` | Case counts and attack rates by county |
| `table_cases_by_type.csv` | Case counts and attack rates by species and production type |
| `q1_spatial_map.pdf` | Spatial map of confirmed cases by species, with county boundaries and high-risk zone (in `submission/` directory) |

## Temporal predictions (Q2)

| File | Description |
|------|-------------|
| `q2_temporal_summary.csv` | Daily summary statistics (mean, median, 2.5th/97.5th percentiles) for new and cumulative confirmed cases over the 4-week prediction period (14 Jan – 10 Feb 2026), plus daily means by species |
| `q2_trajectories.csv` | Individual forward simulation trajectories (500 runs) with daily new confirmed cases |

## Spatial predictions (Q2)

| File | Description |
|------|-------------|
| `q2_spatial_farm.csv` | Predicted infection probability per farm over the 4-week prediction period |
| `q2_spatial_district.csv` | Predicted new cases aggregated by district (summary statistics across simulations) |

## Species contribution (Q3)

| File | Description |
|------|-------------|
| `q3_species_contribution.csv` | Posterior samples of β_duck, chicken/duck susceptibility shares, and implied chicken relative risk |
| `posterior_samples.csv` | MCMC posterior samples from an earlier model iteration; columns `h0` and `α` are from a previous parameterisation. The β_duck, φ_hrz, φ_non, δ columns remain valid. For the β_duck posterior used in Q3 analyses, see `q3_species_contribution.csv` |

## Counterfactual: chicken-only preventive culling (Q4)

| File | Description |
|------|-------------|
| `q4_spillover_temporal_summary.csv` | Daily summary statistics under spillover-only dynamics with chicken-only preventive culling |
| `q4_spillover_trajectories.csv` | Individual trajectories (spillover-only) |
| `q4_spillover_spatial_farm.csv` | Per-farm infection probabilities (spillover-only) |
| `q4_spillover_spatial_district.csv` | District-level predictions (spillover-only) |
| `q4_transmission_temporal_summary.csv` | Daily summary statistics with spatial + movement transmission, chicken-only preventive culling |
| `q4_transmission_trajectories.csv` | Individual trajectories (spillover + transmission) |
| `q4_transmission_spatial_farm.csv` | Per-farm infection probabilities (spillover + transmission) |
| `q4_transmission_spatial_district.csv` | District-level predictions (spillover + transmission) |

## Counterfactual: faster reactive culling (Q5)

| File | Description |
|------|-------------|
| `q5_spillover_temporal_summary.csv` | Daily summary statistics under spillover-only dynamics with faster reactive culling (no preventive culling) |
| `q5_spillover_trajectories.csv` | Individual trajectories (spillover-only) |
| `q5_spillover_spatial_farm.csv` | Per-farm infection probabilities (spillover-only) |
| `q5_spillover_spatial_district.csv` | District-level predictions (spillover-only) |
| `q5_transmission_temporal_summary.csv` | Daily summary statistics with spatial + movement transmission, faster reactive culling |
| `q5_transmission_trajectories.csv` | Individual trajectories (spillover + transmission) |
| `q5_transmission_spatial_farm.csv` | Per-farm infection probabilities (spillover + transmission) |
| `q5_transmission_spatial_district.csv` | District-level predictions (spillover + transmission) |

## Validation

| File | Description |
|------|-------------|
| `validation/` | Directory containing posterior predictive check outputs and model validation results |
