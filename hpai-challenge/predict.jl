# ── Predictions (Phase 3) ─────────────────────────────────────────────────────
# Generate predictions for Q1–Q4 of the HPAI challenge, phase 3.
#
# Q1: Full epidemic description (descriptive — no modelling)
# Q2: When can restocking begin safely?
# Q3: How many outbreaks did preventive culling avert? (bonus)
# Q4: Daily culling capacity for 95% eradication by end of phase 2? (bonus)
#
# Usage:
#   julia --threads=4 predict.jl                # full run (500 draws)
#   julia --threads=4 predict.jl --short        # reduced samples for testing

# ── Source files ──
include("src/config.jl")
include("src/spatial.jl")
include("src/data.jl")
include("src/likelihood.jl")
include("src/models.jl")
include("src/inference.jl")
include("src/simulation.jl")
include("src/generative.jl")
include("src/movement.jl")

using CSV, DataFrames, Statistics, Random, MCMCChains, Dates

const OUTPUT_DIR = joinpath(@__DIR__, "output")

# ── Helpers ──────────────────────────────────────────────────────────────────

function load_population_metadata()
    return CSV.read(POP_CSV, DataFrame)
end

function load_movement_data()
    return CSV.read(MOVEMENT_CSV, DataFrame)
end

"""
    print_scenario_summary(name, totals)

Print median and 95% CrI for a vector of total case counts.
"""
function print_scenario_summary(name::String, totals::Vector)
    println("  $name: median=$(round(median(totals), digits=1)), " *
            "95% CrI [$(round(quantile(totals, 0.025), digits=1))–" *
            "$(round(quantile(totals, 0.975), digits=1))]")
end

# ── Q1: Descriptive summary ─────────────────────────────────────────────────

function run_q1(data::ModelData, pop::DataFrame)
    println("\n═══ Q1: Epidemic description ═══")

    cases_df = CSV.read(CASES_CSV, DataFrame)

    # Epicurve
    epicurve = combine(groupby(cases_df, :date_confirmed), nrow => :n_cases)
    sort!(epicurve, :date_confirmed)
    path = joinpath(OUTPUT_DIR, "q1_epicurve.csv")
    CSV.write(path, epicurve)
    println("  Saved: $path ($(nrow(epicurve)) days)")

    # Species / production breakdown
    cases_full = leftjoin(cases_df, pop; on=:farm_id)
    breakdown = combine(
        groupby(cases_full, [:species, :production]),
        nrow => :n_cases,
    )
    sort!(breakdown, :n_cases, rev=true)
    path = joinpath(OUTPUT_DIR, "q1_species_production.csv")
    CSV.write(path, breakdown)
    println("  Saved: $path")

    # Detection method breakdown
    det = combine(groupby(cases_df, :detection_method), nrow => :n_cases)
    sort!(det, :n_cases, rev=true)
    path = joinpath(OUTPUT_DIR, "q1_detection_methods.csv")
    CSV.write(path, det)
    println("  Saved: $path")

    # Spatial summary by county
    spatial = combine(
        groupby(cases_full, :county),
        nrow => :n_cases,
    )
    sort!(spatial, :n_cases, rev=true)
    path = joinpath(OUTPUT_DIR, "q1_cases_by_county.csv")
    CSV.write(path, spatial)
    println("  Saved: $path")

    total = nrow(cases_df)
    dates = sort(cases_df.date_confirmed)
    println("  Total cases: $total")
    println("  Date range: $(dates[1]) to $(dates[end])")
    println("  Species: $(count(cases_full.species .== "chicken")) chicken, " *
            "$(count(cases_full.species .== "duck")) duck")
end

# ── Q2: Restocking safety ───────────────────────────────────────────────────

"""
    run_q2_restocking(data, draws, pop; n_days_pred, rng)

Sweep restocking start dates and find when restocking is safe.
"""
function run_q2_restocking(
    data::ModelData, draws::Vector, pop::DataFrame;
    n_days_pred::Int = 120,
    pred_movements = nothing,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    println("\n═══ Q2: Restocking safety sweep ═══")

    n_sims = length(draws)

    # Candidate restocking days: weekly from SIM_END to SIM_END + 120
    restock_days = collect(SIM_END:7:(SIM_END + n_days_pred))

    results = DataFrame(
        restocking_day=Int[], restocking_date=Date[],
        p_rebound=Float64[],
        mean_new_cases=Float64[], median_new_cases=Float64[],
        q025_new_cases=Float64[], q975_new_cases=Float64[],
    )

    for r_day in restock_days
        total_cases = zeros(Int, n_sims)

        for i in 1:n_sims
            # Start from SIM_END (not T_LIK): the epidemic is over,
            # so no right-censoring fill-in needed
            sim = simulate_forward(data, draws[i];
                                    scenario=:baseline, n_days_pred,
                                    pred_start=SIM_END,
                                    pred_movements, rng,
                                    restocking_start_day=r_day)
            total_cases[i] = sum(sim.new_confirmed_by_day)
        end

        p_rebound = count(total_cases .> 0) / n_sims
        push!(results, (
            restocking_day=r_day,
            restocking_date=day_to_date(r_day),
            p_rebound=round(p_rebound, digits=4),
            mean_new_cases=round(mean(total_cases), digits=2),
            median_new_cases=round(median(total_cases), digits=1),
            q025_new_cases=round(quantile(Float64.(total_cases), 0.025), digits=1),
            q975_new_cases=round(quantile(Float64.(total_cases), 0.975), digits=1),
        ))

        println("  Restocking $(day_to_date(r_day)): P(rebound)=$(round(p_rebound, digits=3)), " *
                "median cases=$(round(median(total_cases), digits=1))")
    end

    path = joinpath(OUTPUT_DIR, "q2_restocking_safety.csv")
    CSV.write(path, results)
    println("  Saved: $path")

    # Find safe date
    safe_rows = filter(row -> row.p_rebound < 0.05, results)
    if nrow(safe_rows) > 0
        safe_date = safe_rows[1, :restocking_date]
        println("  Safe restocking date (P(rebound) < 5%): $safe_date")
    else
        println("  No date found with P(rebound) < 5% within sweep window")
    end

    return results
end

# ── Q3: Outbreaks averted by preventive culling ─────────────────────────────

"""
    run_q3_averted(data, draws; rng)

Counterfactual: simulate full epidemic with and without preventive culling.
"""
function run_q3_averted(
    data::ModelData, draws::Vector;
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    println("\n═══ Q3: Outbreaks averted by preventive culling ═══")

    n_sims = length(draws)
    T_end = SIM_END

    total_with = zeros(Int, n_sims)
    total_without = zeros(Int, n_sims)

    for i in 1:n_sims
        # Baseline (with preventive culling)
        sim_with = simulate_epidemic(data, draws[i]; rng, T_end)
        total_with[i] = sim_with.total_cases

        # Counterfactual (no preventive culling)
        sim_without = simulate_epidemic(data, draws[i]; rng, T_end,
                                         disable_prev_cull=true)
        total_without[i] = sim_without.total_cases

        if i % 50 == 0
            println("  Q3: $i/$n_sims draws")
        end
    end

    averted = total_without .- total_with

    results = DataFrame(
        trajectory=1:n_sims,
        total_with_cull=total_with,
        total_without_cull=total_without,
        averted=averted,
    )
    path = joinpath(OUTPUT_DIR, "q3_averted_trajectories.csv")
    CSV.write(path, results)
    println("  Saved: $path")

    summary = DataFrame(
        metric=["with_prev_cull", "without_prev_cull", "averted"],
        median=[median(total_with), median(total_without), median(averted)],
        q025=[quantile(Float64.(total_with), 0.025),
              quantile(Float64.(total_without), 0.025),
              quantile(Float64.(averted), 0.025)],
        q975=[quantile(Float64.(total_with), 0.975),
              quantile(Float64.(total_without), 0.975),
              quantile(Float64.(averted), 0.975)],
    )
    path = joinpath(OUTPUT_DIR, "q3_averted_summary.csv")
    CSV.write(path, summary)
    println("  Saved: $path")

    println("  With culling: median=$(round(median(total_with), digits=1))")
    println("  Without culling: median=$(round(median(total_without), digits=1))")
    println("  Averted: median=$(round(median(averted), digits=1)) " *
            "[$(round(quantile(Float64.(averted), 0.025), digits=1))–" *
            "$(round(quantile(Float64.(averted), 0.975), digits=1))]")

    return averted
end

# ── Q4: Daily culling capacity for 95% eradication ──────────────────────────

"""
    run_q4_capacity(data, draws; rng)

Sweep daily culling capacity and find threshold for 95% eradication by day 75.
"""
function run_q4_capacity(
    data::ModelData, draws::Vector;
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    println("\n═══ Q4: Culling capacity sweep ═══")

    n_sims = length(draws)
    T_target = date_to_day(Date(2026, 2, 13))  # end of phase 2

    capacities = [5, 10, 15, 20, 25, 30, 40, 50, 75, 100, typemax(Int)]

    results = DataFrame(
        daily_capacity=Int[],
        capacity_label=String[],
        p_eradicated=Float64[],
        median_cases=Float64[],
        q025_cases=Float64[], q975_cases=Float64[],
    )

    for cap in capacities
        cap_label = cap == typemax(Int) ? "unlimited" : string(cap)
        n_eradicated = 0
        total_cases = zeros(Int, n_sims)

        for i in 1:n_sims
            sim = simulate_epidemic(data, draws[i]; rng,
                                     T_end=T_target,
                                     daily_cull_capacity=cap)
            total_cases[i] = sim.total_cases

            # Eradicated = no active infections at T_target
            any_active = false
            for j in 1:data.N
                if sim.infected[j] && sim.removal_day[j] > T_target
                    any_active = true
                    break
                end
            end
            if !any_active
                n_eradicated += 1
            end
        end

        p_erad = n_eradicated / n_sims
        push!(results, (
            daily_capacity=cap == typemax(Int) ? 9999 : cap,
            capacity_label=cap_label,
            p_eradicated=round(p_erad, digits=4),
            median_cases=round(median(total_cases), digits=1),
            q025_cases=round(quantile(Float64.(total_cases), 0.025), digits=1),
            q975_cases=round(quantile(Float64.(total_cases), 0.975), digits=1),
        ))

        println("  Capacity $cap_label: P(eradicated)=$(round(p_erad, digits=3)), " *
                "median cases=$(round(median(total_cases), digits=1))")
    end

    path = joinpath(OUTPUT_DIR, "q4_capacity_sweep.csv")
    CSV.write(path, results)
    println("  Saved: $path")

    # Find threshold
    safe_rows = filter(row -> row.p_eradicated >= 0.95, results)
    if nrow(safe_rows) > 0
        threshold = safe_rows[1, :capacity_label]
        println("  Minimum capacity for 95% eradication: $threshold farms/day")
    else
        println("  95% eradication not achieved at any tested capacity")
    end

    return results
end

# ── Main ──────────────────────────────────────────────────────────────────────

function main(;
    short::Bool = false,
    n_draws::Int = short ? 100 : 500,
    n_chains::Int = short ? 1 : 4,
    n_samples::Int = short ? 100 : 1000,
    n_warmup::Int = short ? 50 : 500,
)
    mkpath(OUTPUT_DIR)
    rng = Random.MersenneTwister(2025)

    println("═══ HPAI predictions (Phase 3) ═══")
    println("Mode: $(short ? "short" : "full")")
    println("Draws: $n_draws")
    println()

    # ── Load data ──
    data = prepare_model_data()
    pop = load_population_metadata()
    mov_df = load_movement_data()
    in_zone = data.in_surv_zone

    # ── Q1: Descriptive summary (no model fitting needed) ──
    run_q1(data, pop)

    # ── Fit full model ──
    println("\n═══ Fitting full model ═══")
    map_full = find_map_full(data, in_zone)
    m_full = full_model(data, in_zone)
    chains_full = run_inference(m_full, map_full;
                                n_chains, n_samples, n_warmup)
    draws_full = extract_posterior_params(chains_full; n_draws, rng)
    println("Extracted $n_draws posterior draws (full model)")

    # ── Fit gravity model for movements (needed for Q2) ──
    println("\n═══ Gravity model ═══")
    grav = fit_gravity_model(pop, mov_df)

    # ── Q2: Restocking safety ──
    n_days_pred_q2 = 120
    movs_q2 = generate_movements(grav, n_days_pred_q2; rng)
    run_q2_restocking(data, draws_full, pop;
                      n_days_pred=n_days_pred_q2,
                      pred_movements=movs_q2, rng)

    # ── Q3: Outbreaks averted (bonus) ──
    run_q3_averted(data, draws_full; rng)

    # ── Q4: Culling capacity (bonus) ──
    run_q4_capacity(data, draws_full; rng)

    println("\n═══ All predictions complete ═══")
    println("Results in: $OUTPUT_DIR")
end

# Run from command line
if abspath(PROGRAM_FILE) == @__FILE__
    short = "--short" in ARGS
    main(; short)
end
