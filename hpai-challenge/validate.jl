# ── Model validation ──────────────────────────────────────────────────────────
# Prior predictive checks, posterior predictive checks, MCMC diagnostics,
# and synthetic data recovery for the HPAI spatial transmission model.
#
# Usage:
#   julia --threads=4 validate.jl                # full validation
#   julia --threads=4 validate.jl --prior-only   # prior predictive only
#   julia --threads=4 validate.jl --short        # reduced samples for testing

# ── Source files ──
include("src/config.jl")
include("src/spatial.jl")
include("src/data.jl")
include("src/likelihood.jl")
include("src/models.jl")
include("src/inference.jl")
include("src/simulation.jl")
include("src/generative.jl")

using CSV, DataFrames, Statistics, Random, MCMCChains

const VALIDATION_DIR = joinpath(@__DIR__, "output", "validation")

# ── §1 Prior predictive checks ───────────────────────────────────────────────

"""
    run_prior_predictive(data; n_sims, model_type, rng) -> DataFrame

Run prior predictive simulations and return per-simulation summaries.
"""
function run_prior_predictive(
    data::ModelData;
    n_sims::Int = 500,
    model_type::Symbol = :full,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    results = DataFrame(
        sim=Int[], total_cases=Int[], first_case_day=Int[],
        n_duck=Int[], n_chicken=Int[], duck_ratio=Float64[],
        t₀=Float64[], φ_hrz=Float64[], φ_non=Float64[],
        δ=Float64[], β_duck=Float64[],
    )
    if model_type == :full
        results.β = Float64[]
        results.α = Float64[]
        results.p_mov = Float64[]
    end

    for i in 1:n_sims
        params = sample_prior_params(rng; model=model_type)
        sim = simulate_epidemic(data, params; rng)

        n_duck = count(==("duck"), sim.case_species)
        n_chicken = count(==("chicken"), sim.case_species)
        total = sim.total_cases
        duck_ratio = total > 0 ? n_duck / total : 0.0

        row = Dict{Symbol,Any}(
            :sim => i, :total_cases => total,
            :first_case_day => sim.first_case_day,
            :n_duck => n_duck, :n_chicken => n_chicken,
            :duck_ratio => duck_ratio,
            :t₀ => params.t₀, :φ_hrz => params.φ_hrz,
            :φ_non => params.φ_non, :δ => params.δ,
            :β_duck => params.β_duck,
        )
        if model_type == :full
            row[:β] = params.β
            row[:α] = params.α
            row[:p_mov] = params.p_mov
        end
        push!(results, row)

        if i % 100 == 0
            println("  Prior predictive ($model_type): $i/$n_sims")
        end
    end

    return results
end

"""
    check_prior_predictive(df, model_type) -> Bool

Print acceptance checks against workflow criteria and return pass/fail.
"""
function check_prior_predictive(df::DataFrame, model_type::Symbol)
    n = nrow(df)
    println("\n── Prior predictive checks ($model_type, n=$n) ──")

    # Criterion 1: 95% of sims produce 10–500 cases
    in_range = count(r -> 10 <= r.total_cases <= 500, eachrow(df))
    pct_range = round(100 * in_range / n, digits=1)
    pass_range = pct_range >= 95
    println("  Cases 10–500: $pct_range% (target ≥95%) $(pass_range ? "✓" : "✗")")

    # Criterion 2: first case between day 10–30
    nonzero = filter(r -> r.first_case_day > 0, df)
    if nrow(nonzero) > 0
        in_window = count(r -> 10 <= r.first_case_day <= 30, eachrow(nonzero))
        pct_window = round(100 * in_window / nrow(nonzero), digits=1)
    else
        pct_window = 0.0
    end
    pass_window = pct_window >= 80
    println("  First case day 10–30: $pct_window% (target ≥80%) $(pass_window ? "✓" : "✗")")

    # Criterion 3: duck:chicken ratio ~1:3 to 1:5 (duck_ratio ~0.17–0.25)
    has_cases = filter(r -> r.total_cases > 10, df)
    if nrow(has_cases) > 0
        med_ratio = median(has_cases.duck_ratio)
        pass_ratio = 0.10 <= med_ratio <= 0.40
        println("  Median duck ratio: $(round(med_ratio, digits=3)) (target 0.10–0.40) $(pass_ratio ? "✓" : "✗")")
    else
        pass_ratio = false
        println("  Median duck ratio: N/A (no sims with >10 cases) ✗")
    end

    # Summary statistics
    println("\n  Total cases: median=$(median(df.total_cases)), " *
            "IQR=$(quantile(df.total_cases, 0.25))–$(quantile(df.total_cases, 0.75))")
    nonzero_days = filter(>(0), df.first_case_day)
    if !isempty(nonzero_days)
        println("  First case day: median=$(median(nonzero_days)), " *
                "IQR=$(quantile(nonzero_days, 0.25))–$(quantile(nonzero_days, 0.75))")
    end

    return pass_range && pass_window && pass_ratio
end

# ── §3 Posterior predictive checks ────────────────────────────────────────────

"""
    run_posterior_predictive(data, chains; n_sims, rng) -> (summary_df, epicurve_df)

Run posterior predictive simulations and return summaries + epicurve quantiles.
"""
function run_posterior_predictive(
    data::ModelData, chains;
    n_sims::Int = 200,
    rng::AbstractRNG = Random.GLOBAL_RNG,
    convert_params::Union{Function, Nothing} = nothing,
)
    draws = extract_posterior_params(chains; n_draws=n_sims, rng)
    if convert_params !== nothing
        draws = [convert_params(d) for d in draws]
    end

    # Per-sim summaries
    summary_df = DataFrame(
        sim=Int[], total_cases=Int[], first_case_day=Int[],
        n_duck=Int[], n_chicken=Int[], duck_ratio=Float64[],
    )

    # Daily infection counts for epicurve
    all_daily = Matrix{Int}(undef, n_sims, SIM_END)

    for i in 1:n_sims
        sim = simulate_epidemic(data, draws[i]; rng)

        n_duck = count(==("duck"), sim.case_species)
        n_chicken = count(==("chicken"), sim.case_species)
        total = sim.total_cases
        duck_ratio = total > 0 ? n_duck / total : 0.0

        push!(summary_df, (
            sim=i, total_cases=total,
            first_case_day=sim.first_case_day,
            n_duck=n_duck, n_chicken=n_chicken,
            duck_ratio=duck_ratio,
        ))

        all_daily[i, :] = sim.daily_infections

        if i % 50 == 0
            println("  Posterior predictive: $i/$n_sims")
        end
    end

    # Compute epicurve quantiles
    epicurve_df = DataFrame(
        day=Int[], date=Date[],
        median=Float64[], lower=Float64[], upper=Float64[],
    )
    for t in 1:SIM_END
        vals = Float64.(all_daily[:, t])
        push!(epicurve_df, (
            day=t, date=day_to_date(t),
            median=median(vals),
            lower=quantile(vals, 0.025),
            upper=quantile(vals, 0.975),
        ))
    end

    # Print comparison with observed data
    println("\n── Posterior predictive summary ──")
    obs_cases = data.n_cases
    med_cases = median(summary_df.total_cases)
    q025 = quantile(summary_df.total_cases, 0.025)
    q975 = quantile(summary_df.total_cases, 0.975)
    in_cri = q025 <= obs_cases <= q975
    println("  Total cases: observed=$obs_cases, simulated median=$med_cases " *
            "(95% CrI: $(q025)–$(q975)) $(in_cri ? "✓" : "✗")")

    obs_first = minimum(data.case_infect_day)
    med_first = median(filter(>(0), summary_df.first_case_day))
    println("  First case day: observed=$obs_first, simulated median=$med_first")

    has_cases = filter(r -> r.total_cases > 0, summary_df)
    if nrow(has_cases) > 0
        med_ratio = median(has_cases.duck_ratio)
        obs_duck = count(data.is_duck[data.case_idx])
        obs_ratio = obs_duck / obs_cases
        println("  Duck ratio: observed=$(round(obs_ratio, digits=3)), " *
                "simulated median=$(round(med_ratio, digits=3))")
    end

    return summary_df, epicurve_df
end

# ── §4 MCMC diagnostics ──────────────────────────────────────────────────────

"""
    save_mcmc_diagnostics(chains, model_name)

Save parameter summaries and raw trace to CSV.
"""
function save_mcmc_diagnostics(chains, model_name::String)
    param_names = names(chains, :parameters)

    # Parameter summaries
    diag_df = DataFrame(
        parameter=String[], mean=Float64[], std=Float64[],
        q025=Float64[], q25=Float64[], q50=Float64[],
        q75=Float64[], q975=Float64[],
        ess=Float64[], rhat=Float64[],
    )

    for pn in param_names
        vals = vec(chains[:, pn, :].data)
        # ESS and R-hat from MCMCChains
        summ = summarystats(chains[:, [pn], :])
        ess_val = summ[1, :ess_tail]
        rhat_val = summ[1, :rhat]

        push!(diag_df, (
            parameter=string(pn),
            mean=mean(vals), std=std(vals),
            q025=quantile(vals, 0.025), q25=quantile(vals, 0.25),
            q50=quantile(vals, 0.5), q75=quantile(vals, 0.75),
            q975=quantile(vals, 0.975),
            ess=ess_val, rhat=rhat_val,
        ))
    end

    diag_path = joinpath(VALIDATION_DIR, "mcmc_diagnostics_$(model_name).csv")
    CSV.write(diag_path, diag_df)
    println("  Saved: $diag_path")

    # Raw trace
    trace_df = DataFrame()
    for pn in param_names
        for ch in 1:size(chains, 3)
            col_name = Symbol("$(pn)_chain$(ch)")
            trace_df[!, col_name] = vec(chains[:, pn, ch].data)
        end
    end
    trace_path = joinpath(VALIDATION_DIR, "mcmc_trace_$(model_name).csv")
    CSV.write(trace_path, trace_df)
    println("  Saved: $trace_path")

    # Print summary
    println("\n── MCMC diagnostics ($model_name) ──")
    for row in eachrow(diag_df)
        rhat_ok = row.rhat <= 1.05
        ess_ok = row.ess >= 100
        println("  $(row.parameter): R̂=$(round(row.rhat, digits=3)) " *
                "$(rhat_ok ? "✓" : "✗"), " *
                "ESS=$(round(row.ess, digits=0)) " *
                "$(ess_ok ? "✓" : "✗")")
    end
end

# ── §5 Synthetic data recovery ────────────────────────────────────────────────

"""
    run_synthetic_recovery(data, chains, model_type; rng, short) -> DataFrame

Generate synthetic data from posterior medians, fit the model, and check
whether 95% CrIs cover the true values.
"""
function run_synthetic_recovery(
    data::ModelData, chains, model_type::Symbol;
    rng::AbstractRNG = Random.GLOBAL_RNG,
    short::Bool = false,
    n_chains::Int = short ? 1 : 4,
    n_samples::Int = short ? 100 : 1000,
    n_warmup::Int = short ? 50 : 500,
)
    param_names = names(chains, :parameters)

    # True parameters = posterior medians from initial fit
    true_params = Dict{Symbol,Float64}()
    for pn in param_names
        true_params[pn] = median(vec(chains[:, pn, :].data))
    end
    true_nt = NamedTuple{Tuple(keys(true_params))}(values(true_params))

    println("\n── Synthetic recovery ($model_type) ──")
    println("  True parameters:")
    for (k, v) in pairs(true_nt)
        println("    $k = $(round(v, sigdigits=4))")
    end

    # Convert reparameterised params for simulation
    sim_params = model_type == :full_reparam ? reparam_to_sim_params(true_nt) : true_nt

    # Simulate synthetic epidemic
    println("  Simulating synthetic epidemic...")
    sim = simulate_epidemic(data, sim_params; rng)
    println("  Synthetic epidemic: $(sim.total_cases) cases, first day $(sim.first_case_day)")

    if sim.total_cases < 5
        println("  WARNING: too few synthetic cases ($(sim.total_cases)). Skipping recovery.")
        return DataFrame()
    end

    # Prepare synthetic ModelData (filters cases to simulation window internally)
    println("  Preparing synthetic data...")
    synth_data = prepare_synthetic_data(data, sim)
    synth_zone = synth_data.in_surv_zone
    println("  Synthetic data: $(synth_data.n_cases) cases")

    if synth_data.n_cases < 5
        println("  WARNING: too few cases within simulation window. Skipping recovery.")
        return DataFrame()
    end

    # Fit model to synthetic data
    println("  Fitting model to synthetic data...")
    if model_type == :spillover
        map_params = find_map_spillover(synth_data, synth_zone)
        m = spillover_model(synth_data, synth_zone)
    elseif model_type == :full_reparam
        map_params = find_map_full_reparam(synth_data, synth_zone)
        m = full_reparam_model(synth_data, synth_zone)
    else
        map_params = find_map_full(synth_data, synth_zone)
        m = full_model(synth_data, synth_zone)
    end

    synth_chains = run_inference(m, map_params;
                                 n_chains, n_samples, n_warmup)

    # Check coverage
    recovery_df = DataFrame(
        parameter=String[], true_value=Float64[],
        posterior_median=Float64[], q025=Float64[], q975=Float64[],
        covered=Bool[],
    )

    println("\n  Recovery results:")
    all_covered = true
    for pn in param_names
        vals = vec(synth_chains[:, pn, :].data)
        q025 = quantile(vals, 0.025)
        q975 = quantile(vals, 0.975)
        tv = true_params[pn]
        covered = q025 <= tv <= q975
        all_covered = all_covered && covered

        push!(recovery_df, (
            parameter=string(pn),
            true_value=tv,
            posterior_median=median(vals),
            q025=q025, q975=q975,
            covered=covered,
        ))

        println("    $(pn): true=$(round(tv, sigdigits=4)), " *
                "95% CrI=[$(round(q025, sigdigits=4)), $(round(q975, sigdigits=4))] " *
                "$(covered ? "✓" : "✗")")
    end

    n_covered = count(recovery_df.covered)
    n_total = nrow(recovery_df)
    println("  Coverage: $n_covered/$n_total parameters covered by 95% CrI " *
            "$(all_covered ? "✓" : "✗")")

    return recovery_df
end

# ── Main ──────────────────────────────────────────────────────────────────────

function main(;
    prior_only::Bool = false,
    short::Bool = false,
    n_prior_sims::Int = short ? 50 : 500,
    n_posterior_sims::Int = short ? 50 : 200,
    n_chains::Int = short ? 1 : 4,
    n_samples::Int = short ? 100 : 1000,
    n_warmup::Int = short ? 50 : 500,
)
    mkpath(VALIDATION_DIR)
    rng = Random.MersenneTwister(2024)

    println("═══ HPAI model validation ═══")
    println("Mode: $(prior_only ? "prior-only" : "full")" *
            "$(short ? " (short)" : "")")
    println()

    # ── Load data ──
    data = prepare_model_data()
    in_zone = data.in_surv_zone

    # ── §1 Prior predictive checks ──
    for model_type in [:spillover, :full, :full_reparam]
        println("\n═══ Prior predictive: $model_type ═══")
        df = run_prior_predictive(data; n_sims=n_prior_sims,
                                  model_type, rng)
        check_prior_predictive(df, model_type)

        path = joinpath(VALIDATION_DIR, "prior_predictive_$(model_type).csv")
        CSV.write(path, df)
        println("  Saved: $path")
    end

    prior_only && return println("\n── Prior-only mode: done ──")

    # ── §2 Fit models ──
    println("\n═══ Stage 1: Spillover-only model ═══")
    map_spill = find_map_spillover(data, in_zone)
    m_spill = spillover_model(data, in_zone)
    chains_spill = run_inference(m_spill, map_spill;
                                 n_chains, n_samples, n_warmup)

    println("\n═══ Stage 2: Full model ═══")
    map_full = find_map_full(data, in_zone)
    m_full = full_model(data, in_zone)
    chains_full = run_inference(m_full, map_full;
                                n_chains, n_samples, n_warmup)

    println("\n═══ Stage 3: Reparameterised full model ═══")
    map_reparam = find_map_full_reparam(data, in_zone)
    m_reparam = full_reparam_model(data, in_zone)
    chains_reparam = run_inference(m_reparam, map_reparam;
                                    n_chains, n_samples, n_warmup)

    # ── §3 Posterior predictive checks ──
    for (model_name, chains, conv) in [
        ("spillover", chains_spill, nothing),
        ("full", chains_full, nothing),
        ("full_reparam", chains_reparam, reparam_to_sim_params),
    ]
        println("\n═══ Posterior predictive: $model_name ═══")
        summary_df, epicurve_df = run_posterior_predictive(
            data, chains; n_sims=n_posterior_sims, rng,
            convert_params=conv)

        CSV.write(joinpath(VALIDATION_DIR,
                           "posterior_predictive_$(model_name)_summary.csv"), summary_df)
        CSV.write(joinpath(VALIDATION_DIR,
                           "posterior_predictive_$(model_name)_epicurve.csv"), epicurve_df)
        println("  Saved posterior predictive CSVs")
    end

    # ── §4 MCMC diagnostics ──
    println("\n═══ MCMC diagnostics ═══")
    save_mcmc_diagnostics(chains_spill, "spillover")
    save_mcmc_diagnostics(chains_full, "full")
    save_mcmc_diagnostics(chains_reparam, "full_reparam")

    # ── §5 Synthetic data recovery ──
    for (model_name, model_type, chains) in [
        ("spillover", :spillover, chains_spill),
        ("full", :full, chains_full),
        ("full_reparam", :full_reparam, chains_reparam),
    ]
        println("\n═══ Synthetic recovery: $model_name ═══")
        recovery_df = run_synthetic_recovery(
            data, chains, model_type; rng, short, n_chains, n_samples, n_warmup)

        if nrow(recovery_df) > 0
            path = joinpath(VALIDATION_DIR, "synthetic_recovery_$(model_name).csv")
            CSV.write(path, recovery_df)
            println("  Saved: $path")
        end
    end

    println("\n═══ Validation complete ═══")
    println("Results in: $VALIDATION_DIR")
end

# Run from command line
if abspath(PROGRAM_FILE) == @__FILE__
    short = "--short" in ARGS
    prior_only = "--prior-only" in ARGS
    main(; prior_only, short)
end
