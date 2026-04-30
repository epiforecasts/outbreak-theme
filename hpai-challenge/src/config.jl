# ── 0. Configuration ──────────────────────────────────────────────────────────
# Constants, date helpers, and fixed parameters.

using Dates

# Paths (relative to hpai-challenge/)
const DATA_DIR = joinpath(@__DIR__, "..", "data")
const POP_CSV = joinpath(DATA_DIR, "canonical", "population.csv")
const HRZ_GEOJSON = joinpath(DATA_DIR, "canonical", "hrz_32626.geojson")
const CASES_CSV = joinpath(DATA_DIR, "phase-3", "cases.csv")
const ACTIVITY_CSV = joinpath(DATA_DIR, "phase-3", "activity.csv")
const MOVEMENT_CSV = joinpath(DATA_DIR, "phase-3", "movement.csv")
const PREV_CULLS_CSV = joinpath(DATA_DIR, "phase-3", "prev_culls.csv")

# Reference date: day 1 = 1 Dec 2025
const REF_DATE = Date(2025, 11, 30)  # day 0, so day 1 = 1 Dec

date_to_day(d::Date) = Dates.value(d - REF_DATE)
day_to_date(t::Int) = REF_DATE + Day(t)

# Simulation window
const SIM_START = 1    # 1 Dec 2025
const SIM_END = 131    # 10 Apr 2026

# Fixed delay parameters
const MU_E = 3.5       # latent/amplification period (days)
const MU_ID = 5.0      # detection delay (days)
const D_TO_C = 2.0     # suspicion to confirmation (days)
const COMPOUND_DELAY = ceil(Int, MU_E + MU_ID + D_TO_C)  # = 11

# Right-censoring: last day where infections are fully observable
const T_LIK = SIM_END - COMPOUND_DELAY
const REMOVAL_BUFFER = 3  # additional days from confirmation to farm removal

# Within-farm dynamics
const R_GROWTH = 1.0   # within-farm growth rate (per day)
const TAU_MIN = 1      # hard latent period (days)

# Movement and testing
const SIGMA_TEST = 0.9  # pre-shipment testing sensitivity
const EPSILON = 0.5     # zone biosecurity reduction

# Zones
const SURV_ZONE_RADIUS = 10_000.0  # surveillance zone radius (metres)
const SURV_ZONE_DURATION = 28       # surveillance zone duration (days)
const PREV_CULL_START_DAY = date_to_day(Date(2026, 1, 1))  # preventive culls from 1 Jan 2026

# Phased preventive culling
const PREV_CULL_RADIUS_P1 = 1_000.0   # phase 1: 1 km, all species (1–13 Jan)
const PREV_CULL_RADIUS_P2 = 3_000.0   # phase 2: 3 km, chicken only (14 Jan+)
const PREV_CULL_POLICY_CHANGE_DAY = date_to_day(Date(2026, 1, 14))

# Confinement (broiler_2 and organic duck farms)
const CONFINEMENT_START_DAY = date_to_day(Date(2026, 1, 14))
const CONFINEMENT_END_DAY_ORGANIC_DUCK = date_to_day(Date(2026, 2, 14))  # organic duck confinement lifted
const CONFINEMENT_FACTOR = 0.5  # susceptibility reduction for confined farms

# Restocking
const RESTOCKING_DELAY = 14  # days after culling before farm can restock
