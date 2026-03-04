# ── 3–4. Data structure and loading ───────────────────────────────────────────

using CSV, DataFrames, Dates

# ── 3. Data structure ────────────────────────────────────────────────────────

"""
Precomputed data for MCMC evaluation. All arrays use internal 1-based farm
indices (not raw farm_id).
"""
struct ModelData
    # Farm attributes (length N)
    N::Int
    is_duck::BitVector        # true if species == "duck"
    is_hrz::BitVector         # true if farm in HRZ
    x::Vector{Float64}
    y::Vector{Float64}

    # Activity: active[i, t] = true if farm i has birds on day t
    active::BitMatrix         # N × T

    # Spatial neighbours (full population)
    nbr_idx::Vector{Vector{Int}}
    nbr_dist::Vector{Vector{Float64}}

    # Movements: mov_by_day[t] = Vector of (src_idx, dst_idx)
    mov_by_day::Vector{Vector{Tuple{Int,Int}}}

    # Cases (sorted by confirmation day)
    n_cases::Int
    case_idx::Vector{Int}        # farm index for each case
    case_confirm_day::Vector{Int} # T^C
    case_infect_day::Vector{Int}  # T^I = T^C - ceil(d)
    case_removal_day::Vector{Int} # day farm is removed (cull start)

    # Preventive culls (non-case farms culled preventively)
    n_prev_culls::Int
    prev_cull_idx::Vector{Int}
    prev_cull_day::Vector{Int}

    # Zone membership: in_surv_zone[i, t] = true if farm i is in a
    # surveillance zone on day t (precomputed from case locations)
    in_surv_zone::BitMatrix   # N × T

    # Bulk spillover bins for non-case, non-prev-culled farms
    # 8 bins: (duck/chicken) × (hrz/non) × (in_zone/out_zone) per day
    # bulk_counts[bin, t] = number of active farms in that bin on day t
    bulk_counts::Matrix{Int}  # 8 × T

    # Flat neighbour arrays for case farms (cache-efficient)
    # For case c, neighbours are at flat_nbr_idx[flat_nbr_offset[c]:flat_nbr_offset[c+1]-1]
    case_farm_set::BitVector     # N-length: true if farm is a case
    prev_cull_set::BitVector     # N-length: true if farm is preventively culled
    flat_nbr_farm::Vector{Int}   # neighbour farm index
    flat_nbr_dist_val::Vector{Float64}  # distance to neighbour
    flat_nbr_offset::Vector{Int} # offset into flat arrays per case farm

    # Which non-case farms receive spatial hazard from at least one case
    # (for the full model: need per-farm survival tracking)
    susceptible_with_nbrs::Vector{Int}  # farm indices
    susc_flat_case_idx::Vector{Int}     # which case is the neighbour
    susc_flat_dist::Vector{Float64}     # distance
    susc_flat_offset::Vector{Int}       # offsets into flat arrays
end

# ── 4. Data loading ──────────────────────────────────────────────────────────

"""
    bin_index(is_duck, is_hrz, in_zone) -> Int

Map (species, HRZ, zone) to a bin index 1–8.
"""
function bin_index(is_duck::Bool, is_hrz::Bool, in_zone::Bool)::Int
    return 1 + Int(is_duck) + 2 * Int(is_hrz) + 4 * Int(in_zone)
end

"""
    prepare_model_data(; max_dist=50_000.0, verbose=true) -> ModelData

Load all data files, precompute spatial structures, and return a ModelData
struct ready for likelihood evaluation.
"""
function prepare_model_data(; max_dist::Float64 = 50_000.0, verbose::Bool = true)
    T = SIM_END

    # ── Load population ──
    pop = CSV.read(POP_CSV, DataFrame)
    N = nrow(pop)
    id_to_idx = Dict{Int,Int}(id => i for (i, id) in enumerate(pop.farm_id))
    x = Float64.(pop.x)
    y = Float64.(pop.y)
    is_duck = pop.species .== "duck"

    verbose && println("Loaded $N farms ($(sum(is_duck)) duck, $(sum(.!is_duck)) chicken)")

    # ── HRZ membership ──
    is_hrz = identify_hrz_farms(x, y)
    verbose && println("HRZ farms: $(sum(is_hrz))")

    # ── Activity matrix ──
    act_df = CSV.read(ACTIVITY_CSV, DataFrame)
    active = falses(N, T)
    for row in eachrow(act_df)
        idx = get(id_to_idx, row.farm_id, nothing)
        idx === nothing && continue
        d_start = date_to_day(row.date_start)
        d_end = ismissing(row.date_end) ? T : date_to_day(row.date_end)
        for t in max(d_start, SIM_START):min(d_end, T)
            active[idx, t] = true
        end
    end
    verbose && println("Activity loaded")

    # ── Neighbours ──
    nbr_idx, nbr_dist = compute_neighbours(x, y; max_dist)
    verbose && println("Neighbours computed (max_dist=$(max_dist)m)")

    # ── Movements ──
    mov_df = CSV.read(MOVEMENT_CSV, DataFrame)
    mov_by_day = [Tuple{Int,Int}[] for _ in 1:T]
    n_mov_loaded = 0
    for row in eachrow(mov_df)
        t = date_to_day(row.date)
        (t < SIM_START || t > T) && continue
        src = get(id_to_idx, row.source_farm, nothing)
        dst = get(id_to_idx, row.dest_farm, nothing)
        (src === nothing || dst === nothing) && continue
        push!(mov_by_day[t], (src, dst))
        n_mov_loaded += 1
    end
    verbose && println("Movements loaded: $n_mov_loaded in simulation window")

    # ── Cases ──
    cases_df = CSV.read(CASES_CSV, DataFrame)
    # Sort by confirmation date
    sort!(cases_df, :date_confirmed)
    n_cases = nrow(cases_df)

    case_idx = Int[]
    case_confirm_day = Int[]
    case_infect_day = Int[]
    case_removal_day = Int[]

    for row in eachrow(cases_df)
        idx = id_to_idx[row.farm_id]
        tc = date_to_day(row.date_confirmed)
        ti = tc - COMPOUND_DELAY
        if ti < SIM_START
            error("Case farm_id=$(row.farm_id): back-calculated infection day $ti " *
                  "is before simulation start (day $SIM_START). " *
                  "Increase simulation window or reduce compound delay.")
        end
        # Removal = cull_start if available, else confirmation day + 3 (fallback)
        tr = if !ismissing(row.cull_start)
            date_to_day(row.cull_start)
        else
            tc + 3
        end
        push!(case_idx, idx)
        push!(case_confirm_day, tc)
        push!(case_infect_day, ti)
        push!(case_removal_day, tr)
    end

    verbose && println("Cases: $n_cases (infection days $(minimum(case_infect_day))–$(maximum(case_infect_day)))")

    # ── Preventive culls ──
    prev_df = CSV.read(PREV_CULLS_CSV, DataFrame)
    case_farm_ids = Set(cases_df.farm_id)
    prev_cull_idx = Int[]
    prev_cull_day = Int[]

    # Median delay for imputation of missing cull dates
    completed = filter(row -> !ismissing(row.cull_start), prev_df)
    # For prev culls without dates, impute using median of known dates
    # relative to a reference (use the latest case confirmation before each cull)
    median_prev_day = if nrow(completed) > 0
        days = sort!([date_to_day(row.cull_start) for row in eachrow(completed)])
        days[(length(days) + 1) ÷ 2]
    else
        PREV_CULL_START_DAY + 3
    end

    for row in eachrow(prev_df)
        # Skip farms that are also confirmed cases
        row.farm_id in case_farm_ids && continue
        idx = get(id_to_idx, row.farm_id, nothing)
        idx === nothing && continue
        cull_day = if !ismissing(row.cull_start)
            date_to_day(row.cull_start)
        else
            median_prev_day
        end
        (cull_day < SIM_START || cull_day > T) && continue
        push!(prev_cull_idx, idx)
        push!(prev_cull_day, cull_day)
    end
    n_prev_culls = length(prev_cull_idx)
    verbose && println("Preventive culls (non-case): $n_prev_culls")

    # ── Surveillance zones ──
    in_surv_zone = falses(N, T)
    for c in 1:n_cases
        ci = case_idx[c]
        tc = case_confirm_day[c]
        zone_end = min(tc + SURV_ZONE_DURATION - 1, T)
        for i in 1:N
            dx = x[i] - x[ci]
            dy = y[i] - y[ci]
            if dx * dx + dy * dy <= SURV_ZONE_RADIUS^2
                for t in tc:zone_end
                    in_surv_zone[i, t] = true
                end
            end
        end
    end
    verbose && println("Surveillance zones computed")

    # ── Identify case and prev-cull farm sets ──
    case_farm_set = falses(N)
    case_farm_set[case_idx] .= true
    prev_cull_set = falses(N)
    prev_cull_set[prev_cull_idx] .= true

    # ── Bulk spillover bins ──
    # For farms that are neither cases nor preventively culled
    bulk_counts = zeros(Int, 8, T)
    for i in 1:N
        (case_farm_set[i] || prev_cull_set[i]) && continue
        for t in 1:T
            active[i, t] || continue
            bin = bin_index(is_duck[i], is_hrz[i], in_surv_zone[i, t])
            bulk_counts[bin, t] += 1
        end
    end
    verbose && println("Bulk spillover bins computed")

    # ── Flat neighbour arrays for case farms ──
    # For each susceptible farm that is a neighbour of at least one case farm,
    # store which cases are neighbours and at what distance.
    flat_nbr_farm = Int[]
    flat_nbr_dist_val = Float64[]
    flat_nbr_offset = Int[1]

    for c in 1:n_cases
        ci = case_idx[c]
        for (k, ni) in enumerate(nbr_idx[ci])
            push!(flat_nbr_farm, ni)
            push!(flat_nbr_dist_val, nbr_dist[ci][k])
        end
        push!(flat_nbr_offset, length(flat_nbr_farm) + 1)
    end
    verbose && println("Case-farm neighbour arrays: $(length(flat_nbr_farm)) entries")

    # ── Flat neighbour arrays for susceptible farms ──
    # For each non-case farm, which case farms are its neighbours
    susc_nbr_map = Dict{Int, Vector{Tuple{Int,Float64}}}()  # farm_idx => [(case_num, dist)]
    for c in 1:n_cases
        ci = case_idx[c]
        for (k, ni) in enumerate(nbr_idx[ci])
            case_farm_set[ni] && continue  # skip case-to-case
            if !haskey(susc_nbr_map, ni)
                susc_nbr_map[ni] = Tuple{Int,Float64}[]
            end
            push!(susc_nbr_map[ni], (c, nbr_dist[ci][k]))
        end
    end

    susceptible_with_nbrs = sort(collect(keys(susc_nbr_map)))
    susc_flat_case_idx = Int[]
    susc_flat_dist = Float64[]
    susc_flat_offset = Int[1]
    for si in susceptible_with_nbrs
        for (c, d) in susc_nbr_map[si]
            push!(susc_flat_case_idx, c)
            push!(susc_flat_dist, d)
        end
        push!(susc_flat_offset, length(susc_flat_case_idx) + 1)
    end
    verbose && println("Susceptible-farm neighbour arrays: $(length(susceptible_with_nbrs)) farms")

    data = ModelData(
        N, is_duck, is_hrz, x, y,
        active,
        nbr_idx, nbr_dist,
        mov_by_day,
        n_cases, case_idx, case_confirm_day, case_infect_day, case_removal_day,
        n_prev_culls, prev_cull_idx, prev_cull_day,
        in_surv_zone,
        bulk_counts,
        case_farm_set, prev_cull_set,
        flat_nbr_farm, flat_nbr_dist_val, flat_nbr_offset,
        susceptible_with_nbrs, susc_flat_case_idx, susc_flat_dist, susc_flat_offset,
    )

    verbose && println("\nModelData ready: $N farms, $n_cases cases, T=$T days")
    return data
end
