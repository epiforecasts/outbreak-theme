# ── 11. Forward simulation ────────────────────────────────────────────────────

using Random, DataFrames, Dates

"""
    simulate_forward(data, params; scenario, n_days_pred, rng) -> NamedTuple

Simulate new infections forward from SIM_END using posterior parameter draws.
Initialises from observed cases and carries forward surveillance zones.

Scenarios:
- `:baseline` — current policy (Q2)
- `:chicken_only_cull` — preventive cull only chicken farms (Q4)
- `:no_prev_cull_faster_reactive` — no preventive culling, REMOVAL_BUFFER-1 (Q5)

Returns NamedTuple with `new_infections_by_day`, `new_confirmed_by_day`,
`case_farms`, `case_days`, and per-farm `infected` status.
"""
function simulate_forward(
    data::ModelData, params::NamedTuple;
    scenario::Symbol = :baseline,
    n_days_pred::Int = 28,
    pred_movements::Union{Nothing, Vector{Vector{Tuple{Int,Int}}}} = nothing,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    scenario in (:baseline, :chicken_only_cull, :no_prev_cull_faster_reactive) ||
        error("Unknown scenario: $scenario")
    if pred_movements !== nothing && length(pred_movements) < n_days_pred
        error("pred_movements must have at least $n_days_pred days, got $(length(pred_movements))")
    end
    N = data.N
    T_end = SIM_END + n_days_pred

    # Scenario-specific parameters
    removal_buffer = scenario == :no_prev_cull_faster_reactive ? REMOVAL_BUFFER - 1 : REMOVAL_BUFFER
    do_prev_cull = scenario != :no_prev_cull_faster_reactive
    chicken_only_cull = scenario == :chicken_only_cull

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

    (σ_val <= 0 || δ_val <= 0) && error("σ and δ must be positive (got σ=$σ_val, δ=$δ_val)")
    has_transmission && α_val <= 0 && error("α must be positive (got $α_val)")
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

    # Track infection status — initialise from observed cases
    infected = falses(N)
    infect_day = zeros(Int, N)
    removal_day_arr = fill(typemax(Int), N)

    for c in 1:data.n_cases
        ci = data.case_idx[c]
        infected[ci] = true
        infect_day[ci] = data.case_infect_day[c]
        removal_day_arr[ci] = data.case_removal_day[c]
    end

    # Track preventive culls from observed data
    culled = falses(N)
    for p in 1:data.n_prev_culls
        culled[data.prev_cull_idx[p]] = true
    end

    # Carry forward surveillance zones from observed case confirmations
    zone_end_day = zeros(Int, N)
    for c in 1:data.n_cases
        ci = data.case_idx[c]
        tc = data.case_confirm_day[c]
        for i in 1:N
            dx = data.x[i] - data.x[ci]
            dy = data.y[i] - data.y[ci]
            if dx * dx + dy * dy <= SURV_ZONE_RADIUS^2
                zone_end_day[i] = max(zone_end_day[i], tc + SURV_ZONE_DURATION - 1)
            end
        end
    end

    # Outputs (prediction period only)
    new_infections_by_day = zeros(Int, n_days_pred)
    new_confirmed_by_day = zeros(Int, n_days_pred)
    case_farms = Int[]
    case_days = Int[]

    # Track new infections for deferred confirmation
    pending_confirm = Tuple{Int,Int}[]  # (farm_idx, confirm_day)

    for t in (SIM_END + 1):T_end
        pred_idx = t - SIM_END

        # Process pending confirmations that trigger on this day
        for (cf, tc) in pending_confirm
            tc != t && continue

            # Record confirmed case
            if tc <= T_end
                ci_pred = tc - SIM_END
                if 1 <= ci_pred <= n_days_pred
                    new_confirmed_by_day[ci_pred] += 1
                end
            end

            # Establish surveillance zone
            for i in 1:N
                dx = data.x[i] - data.x[cf]
                dy = data.y[i] - data.y[cf]
                if dx * dx + dy * dy <= SURV_ZONE_RADIUS^2
                    zone_end_day[i] = max(zone_end_day[i], tc + SURV_ZONE_DURATION - 1)
                end
            end

            # Preventive culling
            if do_prev_cull && tc >= PREV_CULL_START_DAY
                for i in 1:N
                    (infected[i] || culled[i]) && continue
                    if chicken_only_cull && data.is_duck[i]
                        continue
                    end
                    dx = data.x[i] - data.x[cf]
                    dy = data.y[i] - data.y[cf]
                    if dx * dx + dy * dy <= PREV_CULL_RADIUS^2
                        culled[i] = true
                    end
                end
            end
        end

        # Simulate new infections
        ψt = ψ_profile[t]

        for j in 1:N
            (infected[j] || culled[j]) && continue
            # All farms assumed active in prediction period

            sp = data.is_duck[j] ? β_duck_val : 1.0
            φ = data.is_hrz[j] ? φ_hrz : φ_non
            in_zone = zone_end_day[j] >= t

            # Spillover hazard
            λ = φ * ψt

            # Spatial transmission from infectious neighbours
            if has_transmission
                for (k, ni) in enumerate(data.nbr_idx[j])
                    infected[ni] || continue
                    τ = t - infect_day[ni]
                    (τ < TAU_MIN || t >= removal_day_arr[ni]) && continue
                    τ > T_end && continue
                    λ += β_val * w[τ] * exp(-data.nbr_dist[j][k] * inv_α)
                end
            end

            # Movement hazard (from gravity model)
            if has_transmission && pred_movements !== nothing
                for (src, dst) in pred_movements[pred_idx]
                    dst != j && continue
                    (infected[src] && !culled[src]) || continue
                    τ = t - infect_day[src]
                    (τ < TAU_MIN || t >= removal_day_arr[src]) && continue
                    τ > T_end && continue

                    p_eff = if zone_end_day[src] >= t
                        0.0
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
                removal_day_arr[j] = t + COMPOUND_DELAY + removal_buffer
                new_infections_by_day[pred_idx] += 1
                push!(case_farms, j)
                push!(case_days, t)

                # Schedule confirmation
                tc = t + COMPOUND_DELAY
                push!(pending_confirm, (j, tc))
            end
        end
    end

    return (;
        new_infections_by_day,
        new_confirmed_by_day,
        case_farms,
        case_days,
        infected,
    )
end
