# ── HPAI spatial transmission model ───────────────────────────────────────────
# Main entry point. Includes all source files and provides a `main()` function
# for running inference.
#
# Usage:
#   julia --threads=4 model.jl [--spillover-only] [--short]

# ── Source files ──
include("src/config.jl")
include("src/spatial.jl")
include("src/data.jl")
include("src/likelihood.jl")
include("src/models.jl")
include("src/inference.jl")
include("src/simulation.jl")

# ── 12. Main entry point ─────────────────────────────────────────────────────

function main(;
    spillover_only::Bool = false,
    short::Bool = false,
    n_chains::Int = short ? 1 : 4,
    n_samples::Int = short ? 10 : 1000,
    n_warmup::Int = short ? 10 : 500,
)
    println("═══ HPAI spatial transmission model ═══")
    println("Compound delay: $COMPOUND_DELAY days (ceil of $(MU_E + MU_ID + D_TO_C))")
    println("Simulation window: $(day_to_date(SIM_START)) to $(day_to_date(SIM_END))")
    println()

    # ── Load data ──
    data = prepare_model_data()
    in_zone = data.in_surv_zone

    # ── Stage 1: spillover-only (A + C) ──
    println("\n═══ Stage 1: Spillover-only model (5 parameters) ═══")
    map_spill = find_map_spillover(data, in_zone)

    m_spill = spillover_model(data, in_zone)
    chains_spill = run_inference(m_spill, map_spill;
                                n_chains, n_samples, n_warmup)

    if spillover_only
        println("\n── Spillover-only mode: done ──")
        return (; data, chains_spillover=chains_spill)
    end

    # ── Stage 2: full model (A + B + C) ──
    println("\n═══ Stage 2: Full model (8 parameters) ═══")
    map_full = find_map_full(data, in_zone)

    m_full = full_model(data, in_zone)
    chains_full = run_inference(m_full, map_full;
                                n_chains, n_samples, n_warmup)

    return (; data, chains_spillover=chains_spill, chains_full)
end

# Run from command line
if abspath(PROGRAM_FILE) == @__FILE__
    short = "--short" in ARGS
    spill_only = "--spillover-only" in ARGS
    main(; spillover_only=spill_only, short)
end
