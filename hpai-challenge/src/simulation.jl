# ── 11. Forward simulation ────────────────────────────────────────────────────

using Random, DataFrames, Dates

"""
    simulate_forward(data, chains; n_days_pred=14, n_sims=100, rng=Random.GLOBAL_RNG)
        -> DataFrame

Sample from the posterior and simulate new infections forward in time.
Returns a DataFrame with columns: sim, day, date, new_infections, cum_infections.

Used for Q2 (forecasting) and posterior predictive checks.
"""
function simulate_forward(
    data::ModelData, chains;
    n_days_pred::Int = 14,
    n_sims::Int = 100,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    T_end = SIM_END + n_days_pred
    param_names = names(chains, :parameters)
    n_draws = size(chains, 1) * size(chains, 3)

    results = DataFrame(sim=Int[], day=Int[], date=Date[],
                        new_infections=Int[], cum_infections=Int[])

    for sim in 1:n_sims
        # Sample a random posterior draw
        sample_idx = rand(rng, 1:size(chains, 1))
        chain_idx = rand(rng, 1:size(chains, 3))

        # Extract parameters
        params = Dict(pn => chains[sample_idx, pn, chain_idx] for pn in param_names)

        t₀ = params[:t₀]
        φ_hrz = params[:φ_hrz]
        φ_non = params[:φ_non]
        δ_val = params[:δ]
        β_duck = params[:β_duck]

        has_transmission = haskey(params, :β)
        β_val = has_transmission ? params[:β] : 0.0
        α_val = has_transmission ? params[:α] : 1.0
        p_mov = has_transmission ? params[:p_mov] : 0.0

        inv_α = 1.0 / α_val

        # Infectiousness function
        w = [τ < TAU_MIN ? 0.0 : 1.0 - exp(-R_GROWTH * (τ - TAU_MIN)) for τ in 1:T_end]

        # Track infection status
        infected = falses(data.N)
        infect_day = zeros(Int, data.N)
        removal_day = fill(typemax(Int), data.N)

        # Initialise with observed cases
        for c in 1:data.n_cases
            ci = data.case_idx[c]
            infected[ci] = true
            infect_day[ci] = data.case_infect_day[c]
            removal_day[ci] = data.case_removal_day[c]
        end

        cum = data.n_cases
        for t in (SIM_END + 1):T_end
            new_today = 0
            ψt = t < t₀ ? 0.0 : exp(-δ_val * (t - t₀))

            for j in 1:data.N
                infected[j] && continue
                # Assume all farms active in prediction period
                sp = data.is_duck[j] ? β_duck : 1.0
                φ = data.is_hrz[j] ? φ_hrz : φ_non

                # Spillover
                λ = φ * ψt

                # Spatial transmission
                if has_transmission
                    for (k, ni) in enumerate(data.nbr_idx[j])
                        infected[ni] || continue
                        τ = t - infect_day[ni]
                        (τ < 1 || t >= removal_day[ni]) && continue
                        τ > T_end && continue
                        λ += β_val * w[τ] * exp(-data.nbr_dist[j][k] * inv_α)
                    end
                end

                # Movement hazard (use observed movements if within data range,
                # otherwise no movements in prediction period)
                if has_transmission && t <= SIM_END
                    for (src, dst) in data.mov_by_day[t]
                        dst != j && continue
                        infected[src] || continue
                        τ = t - infect_day[src]
                        (τ < TAU_MIN || t >= removal_day[src]) && continue
                        τ > T_end && continue
                        λ += p_mov * w[τ]
                    end
                end

                λ *= sp
                # No zone effect in prediction period (conservative)

                # Stochastic infection
                if rand(rng) < 1.0 - exp(-λ)
                    infected[j] = true
                    infect_day[j] = t
                    removal_day[j] = t + COMPOUND_DELAY + 3  # rough removal
                    new_today += 1
                end
            end

            cum += new_today
            push!(results, (sim=sim, day=t, date=day_to_date(t),
                            new_infections=new_today, cum_infections=cum))
        end
    end

    return results
end
