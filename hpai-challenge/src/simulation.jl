# ── 11. Forward simulation ────────────────────────────────────────────────────

using Random, DataFrames, Dates

"""
    simulate_forward(data, params; scenario, n_days_pred, rng) -> NamedTuple

Simulate new infections forward from SIM_END using posterior parameter draws.
Initialises from observed cases and carries forward surveillance zones.

Scenarios:
- `:baseline` — current policy (Q2)
- `:no_prev_cull_faster_reactive` — no preventive culling, REMOVAL_BUFFER-1 (Q3 phase 2)
- `:stop_organic_duck_confinement` — remove confinement for organic ducks (Q4a phase 2)
- `:stop_all_confinement` — remove all confinement (Q4b phase 2)
- `:prohibit_restocking` — farms in surveillance zones cannot restock (Q5 phase 2)
- `:no_prev_cull` — full no-preventive-culling counterfactual (Q3 phase 3)

Returns NamedTuple with `new_infections_by_day`, `new_confirmed_by_day`,
`case_farms`, `case_days`, and per-farm `infected` status.
"""
function simulate_forward(
    data::ModelData, params::NamedTuple;
    scenario::Symbol = :baseline,
    n_days_pred::Int = 28,
    pred_start::Int = T_LIK,
    pred_movements::Union{Nothing, Vector{Vector{Tuple{Int,Int}}}} = nothing,
    restocking_start_day::Int = 0,
    daily_cull_capacity::Int = typemax(Int),
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    valid_scenarios = (:baseline, :no_prev_cull_faster_reactive,
                       :stop_organic_duck_confinement, :stop_all_confinement,
                       :prohibit_restocking, :no_prev_cull)
    scenario in valid_scenarios ||
        error("Unknown scenario: $scenario")
    if pred_movements !== nothing && length(pred_movements) < n_days_pred
        error("pred_movements must have at least $n_days_pred days, got $(length(pred_movements))")
    end
    N = data.N
    T_end = SIM_END + n_days_pred

    # Scenario-specific parameters
    removal_buffer = scenario == :no_prev_cull_faster_reactive ? REMOVAL_BUFFER - 1 : REMOVAL_BUFFER
    do_prev_cull = !(scenario in (:no_prev_cull_faster_reactive, :no_prev_cull))
    block_restocking = scenario == :prohibit_restocking

    # Confinement: which farms are confined under this scenario
    is_confined_eff = if scenario == :stop_all_confinement
        falses(N)
    elseif scenario == :stop_organic_duck_confinement
        # Keep broiler_2 confined, remove organic duck confinement
        BitVector(data.production[i] == "broiler_2" for i in 1:N)
    else
        data.is_confined
    end

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

    # Track infection status — initialise from cases infected before pred_start
    infected = falses(N)
    infect_day = zeros(Int, N)
    removal_day_arr = fill(typemax(Int), N)

    for c in 1:data.n_cases
        data.case_infect_day[c] <= pred_start || continue
        ci = data.case_idx[c]
        infected[ci] = true
        infect_day[ci] = data.case_infect_day[c]
        removal_day_arr[ci] = data.case_removal_day[c]
    end

    # Track preventive culls before pred_start
    culled = falses(N)
    cull_day_arr = fill(typemax(Int), N)
    for p in 1:data.n_prev_culls
        data.prev_cull_day[p] <= pred_start || continue
        pi = data.prev_cull_idx[p]
        culled[pi] = true
        cull_day_arr[pi] = data.prev_cull_day[p]
    end

    # Carry forward surveillance zones from confirmations before pred_start
    zone_end_day = zeros(Int, N)
    for c in 1:data.n_cases
        ci = data.case_idx[c]
        tc = data.case_confirm_day[c]
        tc <= pred_start || continue
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

    # Capacity-limited culling queue: (farm_idx, queued_day)
    cull_queue = Tuple{Int,Int}[]

    for t in (pred_start + 1):T_end
        # pred_idx counts from SIM_END (only positive indices are in the output window)
        pred_idx = t - SIM_END

        # Process capacity-limited culling queue
        if daily_cull_capacity < typemax(Int) && !isempty(cull_queue)
            n_to_cull = min(daily_cull_capacity, length(cull_queue))
            for _ in 1:n_to_cull
                (qi, _) = popfirst!(cull_queue)
                if !culled[qi] && !infected[qi]
                    culled[qi] = true
                    cull_day_arr[qi] = t
                end
            end
        end

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

            # Phased preventive culling
            if do_prev_cull && tc >= PREV_CULL_START_DAY
                if tc < PREV_CULL_POLICY_CHANGE_DAY
                    # Phase 1: 1 km, all species
                    radius_sq = PREV_CULL_RADIUS_P1^2
                    for i in 1:N
                        (infected[i] || culled[i]) && continue
                        dx = data.x[i] - data.x[cf]
                        dy = data.y[i] - data.y[cf]
                        if dx * dx + dy * dy <= radius_sq
                            if daily_cull_capacity < typemax(Int)
                                push!(cull_queue, (i, t))
                            else
                                culled[i] = true
                                cull_day_arr[i] = t
                            end
                        end
                    end
                else
                    # Phase 2: 3 km, chicken only
                    radius_sq = PREV_CULL_RADIUS_P2^2
                    for i in 1:N
                        (infected[i] || culled[i]) && continue
                        data.is_duck[i] && continue
                        dx = data.x[i] - data.x[cf]
                        dy = data.y[i] - data.y[cf]
                        if dx * dx + dy * dy <= radius_sq
                            if daily_cull_capacity < typemax(Int)
                                push!(cull_queue, (i, t))
                            else
                                culled[i] = true
                                cull_day_arr[i] = t
                            end
                        end
                    end
                end
            end
        end

        # Restocking: culled farms can restock after RESTOCKING_DELAY
        # (unless in surveillance zone with prohibit_restocking scenario,
        #  or before restocking_start_day)
        for i in 1:N
            culled[i] || continue
            (t - cull_day_arr[i]) >= RESTOCKING_DELAY || continue
            if restocking_start_day > 0 && t < restocking_start_day
                continue
            end
            if block_restocking && zone_end_day[i] >= t
                continue
            end
            culled[i] = false
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

            # Movement hazard: use observed data for t ≤ SIM_END, gravity model after
            if has_transmission && t <= SIM_END && t <= length(data.mov_by_day)
                for (src, dst) in data.mov_by_day[t]
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
            elseif has_transmission && pred_movements !== nothing && t > SIM_END
                mov_idx = t - SIM_END
                for (src, dst) in pred_movements[mov_idx]
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
            # Confinement reduces susceptibility (organic duck lifted on day 76)
            if is_confined_eff[j] && t >= CONFINEMENT_START_DAY &&
               !(data.is_organic_duck[j] && t >= CONFINEMENT_END_DAY_ORGANIC_DUCK)
                λ *= CONFINEMENT_FACTOR
            end

            # Stochastic infection
            if rand(rng) < 1.0 - exp(-λ)
                infected[j] = true
                infect_day[j] = t
                removal_day_arr[j] = t + COMPOUND_DELAY + removal_buffer
                if 1 <= pred_idx <= n_days_pred
                    new_infections_by_day[pred_idx] += 1
                end
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
