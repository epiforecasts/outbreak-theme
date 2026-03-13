# ── Predictions ──────────────────────────────────────────────────────────────
# Generate predictions for Q2–Q5 of the HPAI challenge.
#
# Q2: 4-week baseline prediction (spillover only)
# Q3: Species contribution (posterior β_duck)
# Q4: Chicken-only preventive culling (two scenarios: spillover / with transmission)
# Q5: No preventive culling + faster reactive (two scenarios)
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

# Fixed transmission parameters from full model posterior medians.
# Not independently identifiable from spillover, but used as a sensitivity
# scenario to bracket the effect of culling policy on farm-to-farm spread.
const FIXED_TRANSMISSION = (β = 0.006, α = 4572.0, p_mov = 0.076)

# ── Helpers ──────────────────────────────────────────────────────────────────

"""
Load population metadata for output enrichment.
"""
function load_population_metadata()
    return CSV.read(POP_CSV, DataFrame)
end

"""
Load raw movement data for gravity model fitting.
"""
function load_movement_data()
    return CSV.read(MOVEMENT_CSV, DataFrame)
end

"""
    add_transmission_params(spillover_draw) -> NamedTuple

Merge spillover posterior draw with fixed transmission parameters.
"""
function add_transmission_params(draw::NamedTuple)
    return merge(draw, FIXED_TRANSMISSION)
end

"""
    run_prediction_sims(data, draws, pop; scenario, n_days_pred,
                        pred_movements, rng)

Run forward simulations for each posterior draw and collect results.
"""
function run_prediction_sims(
    data::ModelData, draws::Vector, pop::DataFrame;
    scenario::Symbol = :baseline,
    n_days_pred::Int = 28,
    pred_movements::Union{Nothing, Vector{Vector{Tuple{Int,Int}}}} = nothing,
    rng::AbstractRNG = Random.GLOBAL_RNG,
)
    n_sims = length(draws)
    N = data.N

    pred_dates = [day_to_date(SIM_END + d) for d in 1:n_days_pred]

    all_confirmed = Matrix{Int}(undef, n_sims, n_days_pred)
    all_chicken = Matrix{Int}(undef, n_sims, n_days_pred)
    all_duck = Matrix{Int}(undef, n_sims, n_days_pred)
    all_broiler_1 = Matrix{Int}(undef, n_sims, n_days_pred)
    all_broiler_2 = Matrix{Int}(undef, n_sims, n_days_pred)
    all_conventional = Matrix{Int}(undef, n_sims, n_days_pred)
    all_layer = Matrix{Int}(undef, n_sims, n_days_pred)
    all_organic = Matrix{Int}(undef, n_sims, n_days_pred)
    farm_infected_count = zeros(Int, N)

    for i in 1:n_sims
        sim = simulate_forward(data, draws[i];
                                scenario, n_days_pred, pred_movements, rng)

        for f in sim.case_farms
            farm_infected_count[f] += 1
        end

        daily_chicken = zeros(Int, n_days_pred)
        daily_duck = zeros(Int, n_days_pred)
        daily_broiler_1 = zeros(Int, n_days_pred)
        daily_broiler_2 = zeros(Int, n_days_pred)
        daily_conventional = zeros(Int, n_days_pred)
        daily_layer = zeros(Int, n_days_pred)
        daily_organic = zeros(Int, n_days_pred)

        for (fi, fd) in zip(sim.case_farms, sim.case_days)
            confirm_d = fd + COMPOUND_DELAY
            ci = confirm_d - SIM_END
            (1 <= ci <= n_days_pred) || continue
            sp = pop.species[fi]
            pr = pop.production[fi]
            if sp == "chicken"
                daily_chicken[ci] += 1
            else
                daily_duck[ci] += 1
            end
            if pr == "broiler_1"
                daily_broiler_1[ci] += 1
            elseif pr == "broiler_2"
                daily_broiler_2[ci] += 1
            elseif pr == "conventional"
                daily_conventional[ci] += 1
            elseif pr == "layer"
                daily_layer[ci] += 1
            elseif pr == "organic"
                daily_organic[ci] += 1
            end
        end

        all_confirmed[i, :] = sim.new_confirmed_by_day
        all_chicken[i, :] = daily_chicken
        all_duck[i, :] = daily_duck
        all_broiler_1[i, :] = daily_broiler_1
        all_broiler_2[i, :] = daily_broiler_2
        all_conventional[i, :] = daily_conventional
        all_layer[i, :] = daily_layer
        all_organic[i, :] = daily_organic

        if i % 50 == 0
            println("  Scenario $scenario: $i/$n_sims simulations")
        end
    end

    return (;
        all_confirmed, all_chicken, all_duck,
        all_broiler_1, all_broiler_2, all_conventional,
        all_layer, all_organic,
        farm_infected_count, pred_dates,
    )
end

# ── Output writers ───────────────────────────────────────────────────────────

function write_trajectories(sims::NamedTuple, prefix::String)
    n_sims = size(sims.all_confirmed, 1)
    n_days = size(sims.all_confirmed, 2)

    rows = DataFrame(
        trajectory=Int[], date=Date[],
        new_confirmed=Int[], cum_new_confirmed=Int[],
        new_chicken=Int[], new_duck=Int[],
        new_broiler_1=Int[], new_broiler_2=Int[],
        new_conventional=Int[], new_layer=Int[], new_organic=Int[],
    )

    for i in 1:n_sims
        cum = 0
        for d in 1:n_days
            nc = sims.all_confirmed[i, d]
            cum += nc
            push!(rows, (
                trajectory=i, date=sims.pred_dates[d],
                new_confirmed=nc, cum_new_confirmed=cum,
                new_chicken=sims.all_chicken[i, d],
                new_duck=sims.all_duck[i, d],
                new_broiler_1=sims.all_broiler_1[i, d],
                new_broiler_2=sims.all_broiler_2[i, d],
                new_conventional=sims.all_conventional[i, d],
                new_layer=sims.all_layer[i, d],
                new_organic=sims.all_organic[i, d],
            ))
        end
    end

    path = joinpath(OUTPUT_DIR, "$(prefix)_trajectories.csv")
    CSV.write(path, rows)
    println("  Saved: $path")
end

function write_temporal_summary(sims::NamedTuple, prefix::String)
    n_days = size(sims.all_confirmed, 2)
    cum_matrix = cumsum(sims.all_confirmed, dims=2)

    rows = DataFrame(
        date=Date[],
        new_mean=Float64[], new_median=Float64[],
        new_q025=Float64[], new_q975=Float64[],
        cum_mean=Float64[], cum_median=Float64[],
        cum_q025=Float64[], cum_q975=Float64[],
        new_chicken_mean=Float64[], new_duck_mean=Float64[],
    )

    for d in 1:n_days
        new_vals = Float64.(sims.all_confirmed[:, d])
        cum_vals = Float64.(cum_matrix[:, d])
        push!(rows, (
            date=sims.pred_dates[d],
            new_mean=round(mean(new_vals), digits=2),
            new_median=round(median(new_vals), digits=2),
            new_q025=round(quantile(new_vals, 0.025), digits=2),
            new_q975=round(quantile(new_vals, 0.975), digits=2),
            cum_mean=round(mean(cum_vals), digits=2),
            cum_median=round(median(cum_vals), digits=2),
            cum_q025=round(quantile(cum_vals, 0.025), digits=2),
            cum_q975=round(quantile(cum_vals, 0.975), digits=2),
            new_chicken_mean=round(mean(Float64.(sims.all_chicken[:, d])), digits=2),
            new_duck_mean=round(mean(Float64.(sims.all_duck[:, d])), digits=2),
        ))
    end

    path = joinpath(OUTPUT_DIR, "$(prefix)_temporal_summary.csv")
    CSV.write(path, rows)
    println("  Saved: $path")
end

function write_spatial_farm(
    data::ModelData, sims::NamedTuple, pop::DataFrame, prefix::String,
)
    n_sims = size(sims.all_confirmed, 1)
    N = data.N

    rows = DataFrame(
        farm_id=Int[], x=Float64[], y=Float64[],
        county=String[], district=String[],
        species=String[], production=String[],
        observed_case=Bool[], p_new_case=Float64[],
    )

    for i in 1:N
        p_new = sims.farm_infected_count[i] / n_sims
        push!(rows, (
            farm_id=pop.farm_id[i],
            x=pop.x[i], y=pop.y[i],
            county=pop.county[i], district=pop.district[i],
            species=pop.species[i], production=pop.production[i],
            observed_case=data.case_farm_set[i],
            p_new_case=round(p_new, digits=4),
        ))
    end

    path = joinpath(OUTPUT_DIR, "$(prefix)_spatial_farm.csv")
    CSV.write(path, rows)
    println("  Saved: $path")
end

function write_spatial_district(
    data::ModelData, sims::NamedTuple, pop::DataFrame, prefix::String,
)
    n_sims = size(sims.all_confirmed, 1)
    N = data.N
    farm_risk = sims.farm_infected_count ./ n_sims

    district_data = Dict{Tuple{String,String}, @NamedTuple{
        n_farms::Int, n_existing_cases::Int,
        total_risk::Float64, max_risk::Float64, n_at_risk::Int,
    }}()

    for i in 1:N
        key = (pop.county[i], pop.district[i])
        d = get(district_data, key, (n_farms=0, n_existing_cases=0,
                                      total_risk=0.0, max_risk=0.0, n_at_risk=0))
        district_data[key] = (
            n_farms=d.n_farms + 1,
            n_existing_cases=d.n_existing_cases + Int(data.case_farm_set[i]),
            total_risk=d.total_risk + farm_risk[i],
            max_risk=max(d.max_risk, farm_risk[i]),
            n_at_risk=d.n_at_risk + (farm_risk[i] > 0.0 ? 1 : 0),
        )
    end

    rows = DataFrame(
        county=String[], district=String[],
        n_farms=Int[], n_existing_cases=Int[],
        expected_new_cases=Float64[], max_farm_risk=Float64[],
        n_farms_at_risk=Int[],
    )

    for ((county, district), d) in sort(collect(district_data),
                                         by=x -> -x[2].total_risk)
        push!(rows, (
            county=county, district=district,
            n_farms=d.n_farms, n_existing_cases=d.n_existing_cases,
            expected_new_cases=round(d.total_risk, digits=2),
            max_farm_risk=round(d.max_risk, digits=3),
            n_farms_at_risk=d.n_at_risk,
        ))
    end

    path = joinpath(OUTPUT_DIR, "$(prefix)_spatial_district.csv")
    CSV.write(path, rows)
    println("  Saved: $path")
end

function write_all_outputs(
    data::ModelData, sims::NamedTuple, pop::DataFrame, prefix::String,
)
    write_trajectories(sims, prefix)
    write_temporal_summary(sims, prefix)
    write_spatial_farm(data, sims, pop, prefix)
    write_spatial_district(data, sims, pop, prefix)
end

function write_species_contribution(draws::Vector, pop::DataFrame, data::ModelData)
    n_duck = count(data.is_duck)
    n_chicken = data.N - n_duck
    f_duck = n_duck / data.N
    f_chicken = n_chicken / data.N

    rows = DataFrame(
        β_duck=Float64[],
        chicken_susceptibility_share=Float64[],
        duck_susceptibility_share=Float64[],
        implied_chicken_rr=Float64[],
    )

    for params in draws
        bd = params.β_duck
        total_susc = f_chicken * 1.0 + f_duck * bd
        chicken_share = (f_chicken * 1.0) / total_susc
        duck_share = (f_duck * bd) / total_susc
        implied_rr = 1.0 / bd

        push!(rows, (
            β_duck=round(bd, digits=4),
            chicken_susceptibility_share=round(chicken_share, digits=4),
            duck_susceptibility_share=round(duck_share, digits=4),
            implied_chicken_rr=round(implied_rr, digits=2),
        ))
    end

    path = joinpath(OUTPUT_DIR, "q3_species_contribution.csv")
    CSV.write(path, rows)
    println("  Saved: $path")
end

"""
    print_scenario_summary(name, sims)

Print median and 95% CrI for total confirmed cases.
"""
function print_scenario_summary(name::String, sims::NamedTuple)
    total = vec(sum(sims.all_confirmed, dims=2))
    println("  $name: median=$(round(median(total), digits=1)), " *
            "95% CrI [$(round(quantile(total, 0.025), digits=1))–" *
            "$(round(quantile(total, 0.975), digits=1))]")
end

# ── Main ──────────────────────────────────────────────────────────────────────

function main(;
    short::Bool = false,
    n_draws::Int = short ? 100 : 500,
    n_chains::Int = short ? 1 : 4,
    n_samples::Int = short ? 100 : 1000,
    n_warmup::Int = short ? 50 : 500,
    n_days_pred::Int = 28,
)
    mkpath(OUTPUT_DIR)
    rng = Random.MersenneTwister(2025)

    println("═══ HPAI predictions ═══")
    println("Mode: $(short ? "short" : "full")")
    println("Draws: $n_draws, prediction days: $n_days_pred")
    println()

    # ── Load data ──
    data = prepare_model_data()
    pop = load_population_metadata()
    mov_df = load_movement_data()
    in_zone = data.in_surv_zone

    # ── Fit gravity model for movements ──
    println("\n═══ Gravity model ═══")
    grav = fit_gravity_model(pop, mov_df)

    # ── Fit spillover model ──
    println("\n═══ Fitting spillover model ═══")
    map_spill = find_map_spillover(data, in_zone)
    m_spill = spillover_model(data, in_zone)
    chains_spill = run_inference(m_spill, map_spill;
                                 n_chains, n_samples, n_warmup)

    # ── Extract posterior draws ──
    draws_spillover = extract_posterior_params(chains_spill; n_draws, rng)
    println("Extracted $n_draws posterior draws")

    # Augment draws with fixed transmission parameters
    draws_transmission = [add_transmission_params(d) for d in draws_spillover]

    # ── Q2: Baseline prediction (spillover only) ──
    println("\n═══ Q2: 4-week baseline prediction ═══")
    sims_q2 = run_prediction_sims(data, draws_spillover, pop;
                                   scenario=:baseline, n_days_pred, rng)
    write_all_outputs(data, sims_q2, pop, "q2")

    # ── Q3: Species contribution ──
    println("\n═══ Q3: Species contribution ═══")
    write_species_contribution(draws_spillover, pop, data)

    # ── Q4: Chicken-only preventive culling ──
    # Scenario A: spillover only
    println("\n═══ Q4a: Chicken-only cull (spillover only) ═══")
    sims_q4a = run_prediction_sims(data, draws_spillover, pop;
                                    scenario=:chicken_only_cull, n_days_pred, rng)
    write_all_outputs(data, sims_q4a, pop, "q4_spillover")

    # Scenario B: with transmission + gravity model movements
    println("\n═══ Q4b: Chicken-only cull (with transmission) ═══")
    movs_q4 = generate_movements(grav, n_days_pred; rng)
    sims_q4b = run_prediction_sims(data, draws_transmission, pop;
                                    scenario=:chicken_only_cull, n_days_pred,
                                    pred_movements=movs_q4, rng)
    write_all_outputs(data, sims_q4b, pop, "q4_transmission")

    # ── Q5: No preventive culling + faster reactive ──
    # Scenario A: spillover only
    println("\n═══ Q5a: No prev cull + faster reactive (spillover only) ═══")
    sims_q5a = run_prediction_sims(data, draws_spillover, pop;
                                    scenario=:no_prev_cull_faster_reactive,
                                    n_days_pred, rng)
    write_all_outputs(data, sims_q5a, pop, "q5_spillover")

    # Scenario B: with transmission + gravity model movements
    println("\n═══ Q5b: No prev cull + faster reactive (with transmission) ═══")
    movs_q5 = generate_movements(grav, n_days_pred; rng)
    sims_q5b = run_prediction_sims(data, draws_transmission, pop;
                                    scenario=:no_prev_cull_faster_reactive,
                                    n_days_pred, pred_movements=movs_q5, rng)
    write_all_outputs(data, sims_q5b, pop, "q5_transmission")

    # ── Summary ──
    println("\n═══ Prediction summary ═══")
    print_scenario_summary("Q2 baseline", sims_q2)
    println()
    print_scenario_summary("Q4a chicken-only cull (spillover)", sims_q4a)
    print_scenario_summary("Q4b chicken-only cull (transmission)", sims_q4b)
    println()
    print_scenario_summary("Q5a no prev cull (spillover)", sims_q5a)
    print_scenario_summary("Q5b no prev cull (transmission)", sims_q5b)

    println("\n═══ All predictions complete ═══")
    println("Results in: $OUTPUT_DIR")
end

# Run from command line
if abspath(PROGRAM_FILE) == @__FILE__
    short = "--short" in ARGS
    main(; short)
end
