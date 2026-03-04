# ── 5–7. Likelihood functions ─────────────────────────────────────────────────
# All functions are AD-compatible: no Float64 type annotations on parameters
# or return types, and all intermediate arrays use the promoted type from
# the input parameters (to support ForwardDiff Dual numbers).

# ── 5. Numerics ──────────────────────────────────────────────────────────────

"""
    safe_log1mexp(x)

Compute log(1 - exp(-x)) in a numerically stable way. Assumes x > 0.
"""
function safe_log1mexp(x)
    if x < 0.6931471805599453  # log(2)
        return log(-expm1(-x))
    else
        return log1p(-exp(-x))
    end
end

# ── 6. Spillover-only likelihood (Module A + C) ─────────────────────────────

"""
    foi_loglik_spillover(data, in_zone, t₀, φ_hrz, φ_non, δ, β_duck)

Log-likelihood for the spillover-only model (no farm-to-farm transmission).
Uses bulk aggregation for non-case farm survival.
"""
function foi_loglik_spillover(
    data::ModelData, in_zone::BitMatrix,
    t₀, φ_hrz, φ_non, δ, β_duck,
)
    T = SIM_END
    RT = typeof(t₀ + φ_hrz + φ_non + δ + β_duck)
    ll = zero(RT)

    # Precompute spillover profile ψ(t) for each day
    ψ = Vector{RT}(undef, T)
    @inbounds for t in 1:T
        if t < t₀
            ψ[t] = zero(RT)
        else
            ψ[t] = exp(-δ * (t - t₀))
        end
    end

    zone_mult = RT(1.0 - EPSILON)

    # ── Bulk non-case survival ──
    @inbounds for t in 1:T
        ψt = ψ[t]
        ψt == zero(RT) && continue
        for bin in 1:8
            cnt = data.bulk_counts[bin, t]
            cnt == 0 && continue

            is_duck_b = ((bin - 1) & 1) != 0
            is_hrz_b = ((bin - 1) & 2) != 0
            in_zone_b = ((bin - 1) & 4) != 0

            sp = is_duck_b ? β_duck : one(RT)
            φ = is_hrz_b ? φ_hrz : φ_non
            λ = sp * φ * ψt
            if in_zone_b
                λ *= zone_mult
            end
            ll -= λ * cnt
        end
    end

    # ── Case spillover survival (days 1 to T_j^I - 1) + infection contribution ──
    @inbounds for c in 1:data.n_cases
        ci = data.case_idx[c]
        ti = data.case_infect_day[c]
        sp = data.is_duck[ci] ? β_duck : one(RT)
        φ = data.is_hrz[ci] ? φ_hrz : φ_non

        # Survival before infection day
        for t in 1:(ti - 1)
            ψt = ψ[t]
            ψt == zero(RT) && continue
            data.active[ci, t] || continue
            λ = sp * φ * ψt
            if in_zone[ci, t]
                λ *= zone_mult
            end
            ll -= λ
        end

        # Infection on day T_j^I
        data.active[ci, ti] || continue
        λ_inf = sp * φ * ψ[ti]
        if in_zone[ci, ti]
            λ_inf *= zone_mult
        end
        if λ_inf <= zero(RT)
            return RT(-Inf)
        end
        ll += safe_log1mexp(λ_inf)
    end

    # ── Preventive-cull farm survival (days 1 to cull day - 1) ──
    @inbounds for p in 1:data.n_prev_culls
        pi = data.prev_cull_idx[p]
        cull_t = data.prev_cull_day[p]
        sp = data.is_duck[pi] ? β_duck : one(RT)
        φ = data.is_hrz[pi] ? φ_hrz : φ_non

        for t in 1:(cull_t - 1)
            ψt = ψ[t]
            ψt == zero(RT) && continue
            data.active[pi, t] || continue
            λ = sp * φ * ψt
            if in_zone[pi, t]
                λ *= zone_mult
            end
            ll -= λ
        end
    end

    return ll
end

# ── 7. Full likelihood (Module A + B + C) ────────────────────────────────────

"""
    foi_loglik(data, in_zone, t₀, φ_hrz, φ_non, δ, β_duck, β, α, p_mov)

Full log-likelihood including spillover, spatial transmission, and movement.

Strategy:
- Precompute ψ, w, and kernel values with the promoted Dual type
- Accumulate spatial hazard per susceptible farm per day
- Movement hazard via sparse dict
- Bulk spillover survival for ~9000 non-case farms (8 bins × T)
- Per-case corrections and infection-day contributions
"""
function foi_loglik(
    data::ModelData, in_zone::BitMatrix,
    t₀, φ_hrz, φ_non, δ, β_duck, β_val, α, p_mov,
)
    T = SIM_END
    RT = typeof(t₀ + φ_hrz + φ_non + δ + β_duck + β_val + α + p_mov)
    ll = zero(RT)
    zone_mult = RT(1.0 - EPSILON)

    # Precompute spillover profile
    ψ = Vector{RT}(undef, T)
    @inbounds for t in 1:T
        ψ[t] = t < t₀ ? zero(RT) : exp(-δ * (t - t₀))
    end

    # Precompute infectiousness w(τ) for each possible delay
    w = Vector{RT}(undef, T)
    @inbounds for τ in 1:T
        w[τ] = τ < TAU_MIN ? zero(RT) : one(RT) - exp(-RT(R_GROWTH) * (τ - TAU_MIN))
    end

    # Precompute kernel values for susceptible-farm ↔ case-farm pairs
    inv_α = one(RT) / α
    kernel_vals = Vector{RT}(undef, length(data.susc_flat_dist))
    @inbounds for k in eachindex(data.susc_flat_dist)
        kernel_vals[k] = exp(-data.susc_flat_dist[k] * inv_α)
    end

    # Kernel values for case-farm flat neighbours (case survival/infection)
    case_kernel = Vector{RT}(undef, length(data.flat_nbr_dist_val))
    @inbounds for k in eachindex(data.flat_nbr_dist_val)
        case_kernel[k] = exp(-data.flat_nbr_dist_val[k] * inv_α)
    end

    # ── Spatial hazard for susceptible non-case farms ──
    n_susc = length(data.susceptible_with_nbrs)
    susc_total_hazard = zeros(RT, n_susc, T)

    @inbounds for s in 1:n_susc
        si = data.susceptible_with_nbrs[s]
        off_start = data.susc_flat_offset[s]
        off_end = data.susc_flat_offset[s + 1] - 1

        for t in 1:T
            data.active[si, t] || continue
            h_spatial = zero(RT)
            for k in off_start:off_end
                c = data.susc_flat_case_idx[k]
                ti = data.case_infect_day[c]
                tr = data.case_removal_day[c]
                τ = t - ti
                (τ < TAU_MIN || t >= tr) && continue
                τ > T && continue
                wτ = w[τ]
                wτ == zero(RT) && continue
                h_spatial += wτ * kernel_vals[k]
            end
            susc_total_hazard[s, t] = β_val * h_spatial
        end
    end

    # ── Movement hazard ──
    # Build a reverse lookup: farm_idx → case_number for movement source matching
    farm_to_case = Dict{Int,Int}()
    for c in 1:data.n_cases
        farm_to_case[data.case_idx[c]] = c
    end

    mov_hazard = Dict{Tuple{Int,Int}, RT}()  # (farm_idx, day) → hazard

    @inbounds for t in 1:T
        for (src, dst) in data.mov_by_day[t]
            src_case = get(farm_to_case, src, 0)
            src_case == 0 && continue
            data.active[dst, t] || continue

            ti = data.case_infect_day[src_case]
            tr = data.case_removal_day[src_case]
            τ = t - ti
            (τ < TAU_MIN || t >= tr) && continue
            τ > T && continue
            wτ = w[τ]
            wτ == zero(RT) && continue

            # Effective movement probability
            p_eff = if in_zone[src, t]
                zero(RT)
            elseif data.is_hrz[src]
                p_mov * RT(1.0 - SIGMA_TEST)
            else
                p_mov
            end
            p_eff == zero(RT) && continue

            key = (dst, t)
            mov_hazard[key] = get(mov_hazard, key, zero(RT)) + p_eff * wτ
        end
    end

    # ── Susceptible farm index lookup ──
    susc_idx_map = Dict{Int,Int}()
    for s in 1:n_susc
        susc_idx_map[data.susceptible_with_nbrs[s]] = s
    end

    # ── Bulk non-case survival (spillover) ──
    @inbounds for t in 1:T
        ψt = ψ[t]
        ψt == zero(RT) && continue
        for bin in 1:8
            cnt = data.bulk_counts[bin, t]
            cnt == 0 && continue
            is_duck_b = ((bin - 1) & 1) != 0
            is_hrz_b = ((bin - 1) & 2) != 0
            in_zone_b = ((bin - 1) & 4) != 0
            sp = is_duck_b ? β_duck : one(RT)
            φ = is_hrz_b ? φ_hrz : φ_non
            λ_spill = sp * φ * ψt
            if in_zone_b
                λ_spill *= zone_mult
            end
            ll -= λ_spill * cnt
        end
    end

    # ── Non-case, non-prev-cull spatial + movement survival ──
    # (prev-cull farms handled separately below with truncated time range)
    @inbounds for s in 1:n_susc
        si = data.susceptible_with_nbrs[s]
        data.prev_cull_set[si] && continue
        sp = data.is_duck[si] ? β_duck : one(RT)
        for t in 1:T
            h = susc_total_hazard[s, t]
            h == zero(RT) && continue
            λ_trans = sp * h
            if in_zone[si, t]
                λ_trans *= zone_mult
            end
            ll -= λ_trans
        end
    end

    for ((dst, t), h_mov) in mov_hazard
        data.case_farm_set[dst] && continue
        data.prev_cull_set[dst] && continue
        sp = data.is_duck[dst] ? β_duck : one(RT)
        λ_mov = sp * h_mov
        if in_zone[dst, t]
            λ_mov *= zone_mult
        end
        ll -= λ_mov
    end

    # ── Preventive-cull farm survival (spillover + transmission) ──
    @inbounds for p in 1:data.n_prev_culls
        pi = data.prev_cull_idx[p]
        cull_t = data.prev_cull_day[p]
        sp = data.is_duck[pi] ? β_duck : one(RT)
        φ = data.is_hrz[pi] ? φ_hrz : φ_non

        for t in 1:(cull_t - 1)
            data.active[pi, t] || continue
            λ = φ * ψ[t]

            s_pos = get(susc_idx_map, pi, 0)
            if s_pos > 0
                λ += susc_total_hazard[s_pos, t]
            end

            h_mov = get(mov_hazard, (pi, t), zero(RT))
            λ += h_mov

            λ *= sp
            if in_zone[pi, t]
                λ *= zone_mult
            end
            ll -= λ
        end
    end

    # ── Case contributions ──
    @inbounds for c in 1:data.n_cases
        ci = data.case_idx[c]
        ti = data.case_infect_day[c]
        sp = data.is_duck[ci] ? β_duck : one(RT)
        φ = data.is_hrz[ci] ? φ_hrz : φ_non

        # Survival before infection day
        for t in 1:(ti - 1)
            data.active[ci, t] || continue
            λ = φ * ψ[t]

            # Spatial hazard from other infectious cases
            off_start = data.flat_nbr_offset[c]
            off_end = data.flat_nbr_offset[c + 1] - 1
            h_spatial = zero(RT)
            for k in off_start:off_end
                ni = data.flat_nbr_farm[k]
                data.case_farm_set[ni] || continue
                c2 = get(farm_to_case, ni, 0)
                c2 == 0 && continue
                ti2 = data.case_infect_day[c2]
                tr2 = data.case_removal_day[c2]
                τ2 = t - ti2
                if τ2 >= TAU_MIN && t < tr2 && τ2 <= T
                    h_spatial += w[τ2] * case_kernel[k]
                end
            end
            λ += β_val * h_spatial

            h_mov = get(mov_hazard, (ci, t), zero(RT))
            λ += h_mov

            λ *= sp
            if in_zone[ci, t]
                λ *= zone_mult
            end
            ll -= λ
        end

        # Infection on day T_j^I
        data.active[ci, ti] || continue
        λ_inf = φ * ψ[ti]

        # Spatial hazard on infection day
        off_start = data.flat_nbr_offset[c]
        off_end = data.flat_nbr_offset[c + 1] - 1
        h_spatial = zero(RT)
        for k in off_start:off_end
            ni = data.flat_nbr_farm[k]
            data.case_farm_set[ni] || continue
            c2 = get(farm_to_case, ni, 0)
            c2 == 0 && continue
            ti2 = data.case_infect_day[c2]
            tr2 = data.case_removal_day[c2]
            τ2 = ti - ti2
            if τ2 >= TAU_MIN && ti < tr2 && τ2 <= T
                h_spatial += w[τ2] * case_kernel[k]
            end
        end
        λ_inf += β_val * h_spatial

        h_mov = get(mov_hazard, (ci, ti), zero(RT))
        λ_inf += h_mov

        λ_inf *= sp
        if in_zone[ci, ti]
            λ_inf *= zone_mult
        end

        if λ_inf <= zero(RT)
            return RT(-Inf)
        end
        ll += safe_log1mexp(λ_inf)
    end

    return ll
end
