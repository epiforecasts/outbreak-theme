# ── Generative simulation engine ──────────────────────────────────────────────
# Functions for prior/posterior predictive checks and synthetic data recovery.
# simulate_epidemic mirrors the likelihood (same FOI formulae) but generates
# epidemics from day 1 with no observed data.

using Random, Distributions

# ── Prior sampling ────────────────────────────────────────────────────────────

"""
    sample_prior_params(rng; model=:full) -> NamedTuple

Draw one parameter set from the prior distributions (matching models.jl).
"""
function sample_prior_params(rng::AbstractRNG; model::Symbol = :full)
    t₀ = rand(rng, truncated(Normal(15, 5), 1, 44))
    φ_hrz = rand(rng, LogNormal(log(1e-3), 1.0))
    φ_non = rand(rng, LogNormal(log(1e-4), 1.0))
    δ = rand(rng, Exponential(1 / 50))
    β_duck = rand(rng, Beta(2, 8))
    σ = rand(rng, LogNormal(log(0.3), 1.0))

    if model == :spillover
        return (; t₀, φ_hrz, φ_non, δ, β_duck, σ)
    else
        β = rand(rng, LogNormal(log(1e-4), 1.5))
        α = rand(rng, LogNormal(log(3500), 0.5))
        p_mov = rand(rng, Beta(2, 20))
        return (; t₀, φ_hrz, φ_non, δ, β_duck, σ, β, α, p_mov)
    end
end

# ── Posterior extraction ──────────────────────────────────────────────────────

"""
    extract_posterior_params(chains; n_draws=200, rng) -> Vector{NamedTuple}

Sample parameter sets from posterior chains.
"""
function extract_posterior_params(chains; n_draws::Int = 200,
                                  rng::AbstractRNG = Random.GLOBAL_RNG)
    param_names = names(chains, :parameters)
    n_samples = size(chains, 1)
    n_chains = size(chains, 3)
    NT = NamedTuple{Tuple(param_names)}
    draws = Vector{NT}(undef, n_draws)

    for i in 1:n_draws
        si = rand(rng, 1:n_samples)
        ci = rand(rng, 1:n_chains)
        draws[i] = NT(Tuple(chains[si, pn, ci] for pn in param_names))
    end
    return draws
end

# ── Generative epidemic simulation ───────────────────────────────────────────

"""
    simulate_epidemic(data, params; rng, T_end=SIM_END) -> NamedTuple

Simulate an entire epidemic from day 1, purely generatively.
Uses the same FOI formulae as the likelihood functions:
- Spillover: sp * φ * ψ(t), zone-modulated
- Spatial: β * Σ w(τ) * exp(-d/α) over infectious neighbours
- Movement: p_eff * w(τ) with zone/testing rules
Dynamic surveillance zones and preventive culling are computed on-the-fly.
"""
function simulate_epidemic(
    data::ModelData, params::NamedTuple;
    rng::AbstractRNG = Random.GLOBAL_RNG,
    T_end::Int = SIM_END,
)
    N = data.N

    # Extract parameters
    t₀ = params.t₀
    φ_hrz = params.φ_hrz
    φ_non = params.φ_non
    δ_val = params.δ
    β_duck_val = params.β_duck
    σ_val = params.σ

    has_transmission = haskey(params, :β)
    β_val = has_transmission ? params.β : 0.0
    α_val = has_transmission ? params.α : 1.0
    p_mov_val = has_transmission ? params.p_mov : 0.0

    inv_α = 1.0 / α_val

    # Precompute spillover profile ψ(t) — Bateman function (rise then decay)
    ψ_peak_val = (σ_val / (σ_val + δ_val)) *
                 exp((δ_val / σ_val) * log(δ_val / (σ_val + δ_val)))
    ψ_profile = Vector{Float64}(undef, T_end)
    for t in 1:T_end
        if t < t₀
            ψ_profile[t] = 0.0
        else
            τ = t - t₀
            ψ_profile[t] = (1.0 - exp(-σ_val * τ)) * exp(-δ_val * τ) / ψ_peak_val
        end
    end

    # Precompute infectiousness profile w(τ)
    w = Vector{Float64}(undef, T_end)
    for τ in 1:T_end
        w[τ] = τ < TAU_MIN ? 0.0 : 1.0 - exp(-R_GROWTH * (τ - TAU_MIN))
    end

    # Track infection status
    infected = falses(N)
    infect_day = zeros(Int, N)
    confirm_day = zeros(Int, N)
    removal_day = fill(typemax(Int), N)
    culled = falses(N)  # preventive culls

    # Dynamic surveillance zones: zone_end_day[i] = last day farm i is in zone
    zone_end_day = zeros(Int, N)

    # Outputs
    daily_infections = zeros(Int, T_end)
    case_farms = Int[]
    case_days = Int[]

    for t in 1:T_end
        ψt = ψ_profile[t]

        # Check for newly confirmed cases (confirmation = infection + COMPOUND_DELAY)
        # and trigger surveillance zones + preventive culling
        for ci in eachindex(case_farms)
            tc = case_days[ci] + COMPOUND_DELAY
            tc != t && continue
            cf = case_farms[ci]
            confirm_day[cf] = tc

            # Establish surveillance zone around confirmed case
            for i in 1:N
                dx = data.x[i] - data.x[cf]
                dy = data.y[i] - data.y[cf]
                if dx * dx + dy * dy <= SURV_ZONE_RADIUS^2
                    zone_end_day[i] = max(zone_end_day[i], tc + SURV_ZONE_DURATION - 1)
                end
            end

            # Preventive culling (from PREV_CULL_START_DAY onwards)
            if tc >= PREV_CULL_START_DAY
                for i in 1:N
                    (infected[i] || culled[i]) && continue
                    dx = data.x[i] - data.x[cf]
                    dy = data.y[i] - data.y[cf]
                    if dx * dx + dy * dy <= PREV_CULL_RADIUS^2
                        culled[i] = true
                    end
                end
            end
        end

        # Simulate new infections
        T_active = size(data.active, 2)
        for j in 1:N
            (infected[j] || culled[j]) && continue
            # Activity check: use activity matrix if within range, assume active beyond
            if t <= T_active
                data.active[j, t] || continue
            end

            sp = data.is_duck[j] ? β_duck_val : 1.0
            φ = data.is_hrz[j] ? φ_hrz : φ_non
            in_zone = zone_end_day[j] >= t

            # Spillover hazard
            λ = φ * ψt

            # Spatial transmission from all infectious neighbours
            if has_transmission
                for (k, ni) in enumerate(data.nbr_idx[j])
                    infected[ni] || continue
                    τ = t - infect_day[ni]
                    (τ < TAU_MIN || t >= removal_day[ni]) && continue
                    τ > T_end && continue
                    λ += β_val * w[τ] * exp(-data.nbr_dist[j][k] * inv_α)
                end
            end

            # Movement hazard
            if has_transmission && t <= length(data.mov_by_day)
                for (src, dst) in data.mov_by_day[t]
                    dst != j && continue
                    infected[src] || continue
                    τ = t - infect_day[src]
                    (τ < TAU_MIN || t >= removal_day[src]) && continue
                    τ > T_end && continue

                    # Effective movement probability with zone/testing rules
                    p_eff = if zone_end_day[src] >= t
                        0.0  # no movements from farms in surveillance zone
                    elseif data.is_hrz[src]
                        p_mov_val * (1.0 - SIGMA_TEST)
                    else
                        p_mov_val
                    end
                    p_eff == 0.0 && continue
                    λ += p_eff * w[τ]
                end
            end

            λ *= sp
            if in_zone
                λ *= (1.0 - EPSILON)
            end

            # Stochastic infection
            if rand(rng) < 1.0 - exp(-λ)
                infected[j] = true
                infect_day[j] = t
                removal_day[j] = t + COMPOUND_DELAY + REMOVAL_BUFFER
                daily_infections[t] += 1
                push!(case_farms, j)
                push!(case_days, t)
            end
        end
    end

    total = sum(daily_infections)
    cum_infections = cumsum(daily_infections)
    first_day = total > 0 ? findfirst(>(0), daily_infections) : 0

    # Collect species info for cases
    case_species = [data.is_duck[f] ? "duck" : "chicken" for f in case_farms]
    case_x = [data.x[f] for f in case_farms]
    case_y = [data.y[f] for f in case_farms]

    return (;
        daily_infections, cum_infections, total_cases=total,
        first_case_day=first_day, case_farms, case_days,
        case_species, case_x, case_y,
        infected, infect_day, confirm_day, removal_day, culled,
    )
end

# ── Synthetic data preparation ────────────────────────────────────────────────

"""
    prepare_synthetic_data(data, sim_result) -> ModelData

Construct a modified ModelData from a simulated epidemic result.
Keeps farm attributes, neighbours, activity, and movements from the original
data, but replaces case-related fields with the simulated epidemic.
"""
function prepare_synthetic_data(data::ModelData, sim::NamedTuple)
    T = SIM_END
    N = data.N

    all_case_idx = sim.case_farms
    all_infect_day = sim.case_days
    all_confirm_day = all_infect_day .+ COMPOUND_DELAY

    # Filter to cases whose infection day falls within the simulation window
    # (matching prepare_model_data which requires case_infect_day >= SIM_START)
    valid = findall(d -> SIM_START <= d <= T, all_infect_day)
    case_idx = all_case_idx[valid]
    case_infect_day = all_infect_day[valid]
    case_confirm_day = all_confirm_day[valid]

    # Sort by confirmation day (matching prepare_model_data behaviour)
    order = sortperm(case_confirm_day)
    case_idx = case_idx[order]
    case_infect_day = case_infect_day[order]
    case_confirm_day = case_confirm_day[order]
    n_cases = length(case_idx)

    # Removal day: use confirm + REMOVAL_BUFFER, clamped to be usable by likelihood
    case_removal_day = case_confirm_day .+ REMOVAL_BUFFER

    # Recompute surveillance zones from simulated cases
    in_surv_zone = falses(N, T)
    for c in 1:n_cases
        ci = case_idx[c]
        tc = case_confirm_day[c]
        tc > T && continue
        zone_end = min(tc + SURV_ZONE_DURATION - 1, T)
        for i in 1:N
            dx = data.x[i] - data.x[ci]
            dy = data.y[i] - data.y[ci]
            if dx * dx + dy * dy <= SURV_ZONE_RADIUS^2
                for t in tc:zone_end
                    t > T && break
                    in_surv_zone[i, t] = true
                end
            end
        end
    end

    # Recompute preventive culls from simulated cases
    prev_cull_idx = Int[]
    prev_cull_day = Int[]
    case_farm_set = falses(N)
    case_farm_set[case_idx] .= true
    prev_cull_set = falses(N)

    for c in 1:n_cases
        ci = case_idx[c]
        tc = case_confirm_day[c]
        (tc < PREV_CULL_START_DAY || tc > T) && continue
        for i in 1:N
            (case_farm_set[i] || prev_cull_set[i]) && continue
            dx = data.x[i] - data.x[ci]
            dy = data.y[i] - data.y[ci]
            if dx * dx + dy * dy <= PREV_CULL_RADIUS^2
                push!(prev_cull_idx, i)
                push!(prev_cull_day, tc)
                prev_cull_set[i] = true
            end
        end
    end
    n_prev_culls = length(prev_cull_idx)

    # Recompute bulk counts
    bulk_counts = zeros(Int, 8, T)
    for i in 1:N
        (case_farm_set[i] || prev_cull_set[i]) && continue
        for t in 1:T
            data.active[i, t] || continue
            bin = bin_index(data.is_duck[i], data.is_hrz[i], in_surv_zone[i, t])
            bulk_counts[bin, t] += 1
        end
    end

    # Recompute flat neighbour arrays for new case farms
    flat_nbr_farm = Int[]
    flat_nbr_dist_val = Float64[]
    flat_nbr_offset = Int[1]

    for c in 1:n_cases
        ci = case_idx[c]
        for (k, ni) in enumerate(data.nbr_idx[ci])
            push!(flat_nbr_farm, ni)
            push!(flat_nbr_dist_val, data.nbr_dist[ci][k])
        end
        push!(flat_nbr_offset, length(flat_nbr_farm) + 1)
    end

    # Susceptible-farm flat arrays
    susc_nbr_map = Dict{Int, Vector{Tuple{Int,Float64}}}()
    for c in 1:n_cases
        ci = case_idx[c]
        for (k, ni) in enumerate(data.nbr_idx[ci])
            case_farm_set[ni] && continue
            if !haskey(susc_nbr_map, ni)
                susc_nbr_map[ni] = Tuple{Int,Float64}[]
            end
            push!(susc_nbr_map[ni], (c, data.nbr_dist[ci][k]))
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

    return ModelData(
        N, data.is_duck, data.is_hrz, data.x, data.y,
        data.active,
        data.nbr_idx, data.nbr_dist,
        data.mov_by_day,
        n_cases, case_idx, case_confirm_day, case_infect_day, case_removal_day,
        n_prev_culls, prev_cull_idx, prev_cull_day,
        in_surv_zone,
        bulk_counts,
        case_farm_set, prev_cull_set,
        flat_nbr_farm, flat_nbr_dist_val, flat_nbr_offset,
        susceptible_with_nbrs, susc_flat_case_idx, susc_flat_dist, susc_flat_offset,
    )
end
