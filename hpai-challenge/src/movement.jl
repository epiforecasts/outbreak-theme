# ── Gravity model for poultry movements ──────────────────────────────────────
# Movements are exclusively broiler_1 (source) → broiler_2 (destination).
# Model: P(i→j) ∝ cap_i · cap_j · exp(-d_ij/κ)
# Fitted from observed movement data.

using CSV, DataFrames, Statistics, Distributions, Random

struct GravityModel
    daily_rate::Float64           # Poisson rate (movements/day)
    κ::Float64                    # distance decay scale (metres)
    src_idx::Vector{Int}          # farm indices of broiler_1 sources
    dst_idx::Vector{Int}          # farm indices of broiler_2 destinations
    src_cap::Vector{Float64}      # source capacities
    dst_cap::Vector{Float64}      # destination capacities
    src_cumw::Vector{Float64}     # cumulative source weights for sampling
    x::Vector{Float64}            # all farm x coordinates
    y::Vector{Float64}            # all farm y coordinates
end

"""
    fit_gravity_model(pop, mov_df; verbose) -> GravityModel

Fit a gravity model from observed movement data.
Kernel: exponential exp(-d/κ). MLE with 2D pair density correction: κ = mean(d)/2.
"""
function fit_gravity_model(pop::DataFrame, mov_df::DataFrame; verbose::Bool = true)
    id_to_idx = Dict(id => i for (i, id) in enumerate(pop.farm_id))
    N = nrow(pop)

    # Identify source (broiler_1) and destination (broiler_2) pools
    src_idx = findall(pop.production .== "broiler_1")
    dst_idx = findall(pop.production .== "broiler_2")

    # Filter to broiler_1 → broiler_2 movements only
    src_set = Set(pop.farm_id[src_idx])
    dst_set = Set(pop.farm_id[dst_idx])

    # Compute movement distances
    distances = Float64[]
    n_filtered = 0
    filtered_dates = Set{eltype(mov_df.date)}()
    for row in eachrow(mov_df)
        src = get(id_to_idx, row.source_farm, nothing)
        dst = get(id_to_idx, row.dest_farm, nothing)
        (src === nothing || dst === nothing) && continue
        (row.source_farm in src_set && row.dest_farm in dst_set) || continue
        n_filtered += 1
        push!(filtered_dates, row.date)
        dx = pop.x[src] - pop.x[dst]
        dy = pop.y[src] - pop.y[dst]
        push!(distances, sqrt(dx^2 + dy^2))
    end

    isempty(distances) && error("No valid broiler_1→broiler_2 movements found; cannot fit gravity model.")

    # MLE for exponential kernel with 2D pair density: κ = mean(d) / 2
    κ = mean(distances) / 2.0
    (!isfinite(κ) || κ <= 0) && error("Fitted κ is invalid ($κ); check movement distances.")

    # Daily rate
    n_dates = length(filtered_dates)
    daily_rate = n_filtered / Float64(n_dates)

    # Source weights ∝ capacity
    src_cap = Float64.(pop.capacity[src_idx])
    src_cumw = cumsum(src_cap)

    # Destination capacities
    dst_cap = Float64.(pop.capacity[dst_idx])

    if verbose
        println("Gravity model: κ=$(round(κ/1000, digits=1))km, " *
                "rate=$(round(daily_rate, digits=1))/day")
        println("  Sources: $(length(src_idx)) broiler_1 farms")
        println("  Destinations: $(length(dst_idx)) broiler_2 farms")
        println("  Movement distances: median=$(round(median(distances)/1000, digits=1))km, " *
                "mean=$(round(mean(distances)/1000, digits=1))km")
    end

    return GravityModel(daily_rate, κ, src_idx, dst_idx,
                         src_cap, dst_cap, src_cumw,
                         Float64.(pop.x), Float64.(pop.y))
end

"""
    generate_movements(grav, n_days; rng) -> Vector{Vector{Tuple{Int,Int}}}

Generate synthetic movements for `n_days` using the fitted gravity model.
Returns movements indexed by farm index (not farm_id).
"""
function generate_movements(
    grav::GravityModel, n_days::Int;
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    n_src = length(grav.src_idx)
    n_dst = length(grav.dst_idx)
    (n_src == 0 || n_dst == 0) && return [Tuple{Int,Int}[] for _ in 1:n_days]
    src_total = grav.src_cumw[end]
    inv_κ = 1.0 / grav.κ

    mov_by_day = [Tuple{Int,Int}[] for _ in 1:n_days]

    for d in 1:n_days
        n_mov = rand(rng, Poisson(grav.daily_rate))

        for _ in 1:n_mov
            # Sample source ∝ capacity
            u = rand(rng) * src_total
            si = searchsortedfirst(grav.src_cumw, u)
            si > n_src && continue
            src_farm = grav.src_idx[si]

            # Compute destination weights: cap_j * exp(-d_ij / κ)
            sx = grav.x[src_farm]
            sy = grav.y[src_farm]
            weights = Vector{Float64}(undef, n_dst)
            @inbounds for k in 1:n_dst
                df = grav.dst_idx[k]
                dx = sx - grav.x[df]
                dy = sy - grav.y[df]
                dist = sqrt(dx * dx + dy * dy)
                weights[k] = grav.dst_cap[k] * exp(-dist * inv_κ)
            end

            # Sample destination from weights
            total_w = sum(weights)
            total_w <= 0.0 && continue
            u2 = rand(rng) * total_w
            cum = 0.0
            dst_farm = grav.dst_idx[1]
            @inbounds for k in 1:n_dst
                cum += weights[k]
                if u2 <= cum
                    dst_farm = grav.dst_idx[k]
                    break
                end
            end

            push!(mov_by_day[d], (src_farm, dst_farm))
        end
    end

    return mov_by_day
end
