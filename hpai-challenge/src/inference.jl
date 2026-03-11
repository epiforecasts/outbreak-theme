# ── 9–10. MAP initialisation and MCMC runner ──────────────────────────────────

using Turing, Distributions, MCMCChains
using Random

# ── 9. MAP initialisation ────────────────────────────────────────────────────

"""
    find_map_spillover(data, in_zone; n_samples=2000) -> NamedTuple

Random search from prior for the 6 spillover parameters to find an approximate MAP.
Returns a NamedTuple suitable for chain initialisation.
"""
function find_map_spillover(data::ModelData, in_zone::BitMatrix; n_samples::Int = 2000)
    best_ll = -Inf
    best_params = nothing

    rng = Random.MersenneTwister(42)

    d_t₀ = truncated(Normal(15, 5), 1, 44)
    d_φ_hrz = LogNormal(log(1e-3), 1.0)
    d_φ_non = LogNormal(log(1e-4), 1.0)
    d_δ = Exponential(1/50)
    d_β_duck = Beta(2, 8)
    d_σ = LogNormal(log(0.3), 1.0)

    for _ in 1:n_samples
        t₀ = rand(rng, d_t₀)
        φ_hrz = rand(rng, d_φ_hrz)
        φ_non = rand(rng, d_φ_non)
        δ = rand(rng, d_δ)
        β_duck = rand(rng, d_β_duck)
        σ = rand(rng, d_σ)

        ll = foi_loglik_spillover(data, in_zone, t₀, φ_hrz, φ_non, δ, β_duck, σ)
        isinf(ll) && continue
        lp = logpdf(d_t₀, t₀) + logpdf(d_φ_hrz, φ_hrz) +
             logpdf(d_φ_non, φ_non) + logpdf(d_δ, δ) +
             logpdf(d_β_duck, β_duck) + logpdf(d_σ, σ)
        total = ll + lp
        if total > best_ll
            best_ll = total
            best_params = (t₀=t₀, φ_hrz=φ_hrz, φ_non=φ_non, δ=δ, β_duck=β_duck, σ=σ)
        end
    end

    if best_params === nothing
        error("MAP spillover search found no finite posterior evaluations. " *
              "Check data and priors.")
    end
    println("MAP spillover: log-posterior = $(round(best_ll, digits=2))")
    println("  t₀=$(round(best_params.t₀, digits=1)), " *
            "φ_hrz=$(round(best_params.φ_hrz, sigdigits=3)), " *
            "φ_non=$(round(best_params.φ_non, sigdigits=3)), " *
            "δ=$(round(best_params.δ, sigdigits=3)), " *
            "β_duck=$(round(best_params.β_duck, digits=3)), " *
            "σ=$(round(best_params.σ, sigdigits=3))")
    return best_params
end

"""
    find_map_full(data, in_zone; n_samples=1000) -> NamedTuple

Random search from prior for the full 8-parameter model.
Returns a NamedTuple suitable for chain initialisation.
"""
function find_map_full(data::ModelData, in_zone::BitMatrix; n_samples::Int = 1000)
    best_ll = -Inf
    best_params = nothing

    rng = Random.MersenneTwister(42)

    d_t₀ = truncated(Normal(15, 5), 1, 44)
    d_φ_hrz = LogNormal(log(1e-3), 1.0)
    d_φ_non = LogNormal(log(1e-4), 1.0)
    d_δ = Exponential(1/50)
    d_β_duck = Beta(2, 8)
    d_σ = LogNormal(log(0.3), 1.0)
    d_β = LogNormal(log(1e-4), 1.5)
    d_α = LogNormal(log(3500), 0.5)
    d_p_mov = Beta(2, 20)

    for _ in 1:n_samples
        t₀ = rand(rng, d_t₀)
        φ_hrz = rand(rng, d_φ_hrz)
        φ_non = rand(rng, d_φ_non)
        δ = rand(rng, d_δ)
        β_duck = rand(rng, d_β_duck)
        σ = rand(rng, d_σ)
        β_val = rand(rng, d_β)
        α = rand(rng, d_α)
        p_mov = rand(rng, d_p_mov)

        ll = foi_loglik(data, in_zone, t₀, φ_hrz, φ_non, δ, β_duck, σ, β_val, α, p_mov)
        isinf(ll) && continue

        lp = logpdf(d_t₀, t₀) + logpdf(d_φ_hrz, φ_hrz) +
             logpdf(d_φ_non, φ_non) + logpdf(d_δ, δ) +
             logpdf(d_β_duck, β_duck) + logpdf(d_σ, σ) +
             logpdf(d_β, β_val) + logpdf(d_α, α) + logpdf(d_p_mov, p_mov)
        total = ll + lp
        if total > best_ll
            best_ll = total
            best_params = (t₀=t₀, φ_hrz=φ_hrz, φ_non=φ_non, δ=δ,
                           β_duck=β_duck, σ=σ, β=β_val, α=α, p_mov=p_mov)
        end
    end

    if best_params === nothing
        error("MAP full search found no finite posterior evaluations. " *
              "Check data and priors.")
    end
    println("MAP full: log-posterior = $(round(best_ll, digits=2))")
    println("  t₀=$(round(best_params.t₀, digits=1)), " *
            "φ_hrz=$(round(best_params.φ_hrz, sigdigits=3)), " *
            "φ_non=$(round(best_params.φ_non, sigdigits=3)), " *
            "δ=$(round(best_params.δ, sigdigits=3)), " *
            "β_duck=$(round(best_params.β_duck, digits=3)), " *
            "σ=$(round(best_params.σ, sigdigits=3)), " *
            "β=$(round(best_params.β, sigdigits=3)), " *
            "α=$(round(best_params.α, digits=0)), " *
            "p_mov=$(round(best_params.p_mov, digits=4))")
    return best_params
end

# ── 10. MCMC runner ──────────────────────────────────────────────────────────

"""
    run_inference(model, init_params;
                  n_chains=4, n_samples=1000, n_warmup=500,
                  target_accept=0.8, max_depth=10) -> Chains

Run NUTS MCMC with parallel chains and print convergence diagnostics.
Uses the MAP estimate for initialisation via Turing's initial_params.
"""
function run_inference(
    model, init_params::NamedTuple;
    n_chains::Int = 4,
    n_samples::Int = 1000,
    n_warmup::Int = 500,
    target_accept::Float64 = 0.8,
    max_depth::Int = 10,
)
    sampler = NUTS(n_warmup, target_accept; max_depth)

    println("\nRunning MCMC: $n_chains chains × $n_samples samples (+ $n_warmup warmup)")
    println("Target acceptance: $target_accept, max tree depth: $max_depth")

    # Turing expects init_params as a vector of values in parameter declaration order.
    base_vec = collect(values(init_params))

    chains = if n_chains > 1
        # Small multiplicative jitter to diversify chain starting points
        jitter_rng = Random.MersenneTwister(42)
        init_vecs = [base_vec .* (1.0 .+ 0.01 .* randn(jitter_rng, length(base_vec)))
                     for _ in 1:n_chains]
        sample(model, sampler, MCMCSerial(), n_samples, n_chains;
               initial_params=init_vecs, progress=true)
    else
        sample(model, sampler, n_samples;
               initial_params=base_vec, progress=true)
    end

    print_diagnostics(chains)
    return chains
end

"""
    print_diagnostics(chains)

Print convergence diagnostics: R-hat, ESS, and divergence summary.
"""
function print_diagnostics(chains::Chains)
    println("\n── Convergence diagnostics ──")
    println(chains)
end
