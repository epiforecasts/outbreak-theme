# ── 8. Turing model definitions ───────────────────────────────────────────────

using Turing
using Distributions

"""
Spillover-only model (Module A + C). 5 estimated parameters.
Captures early dynamics and HRZ contrast without farm-to-farm transmission.
"""
Turing.@model function spillover_model(data::ModelData, in_zone::BitMatrix)
    t₀ ~ truncated(Normal(15, 5), 1, 44)
    φ_hrz ~ LogNormal(log(1e-3), 1.0)
    φ_non ~ LogNormal(log(1e-4), 1.0)
    δ ~ Exponential(1/50)     # scale parameterisation: mean = 0.02/day
    β_duck ~ Beta(2, 8)

    Turing.@addlogprob! foi_loglik_spillover(data, in_zone, t₀, φ_hrz, φ_non, δ, β_duck)
end

"""
Full model (Module A + B + C). 8 estimated parameters.
Includes spillover, spatial transmission, and movement transmission.
"""
Turing.@model function full_model(data::ModelData, in_zone::BitMatrix)
    t₀ ~ truncated(Normal(15, 5), 1, 44)
    φ_hrz ~ LogNormal(log(1e-3), 1.0)
    φ_non ~ LogNormal(log(1e-4), 1.0)
    δ ~ Exponential(1/50)     # scale parameterisation: mean = 0.02/day
    β_duck ~ Beta(2, 8)
    β ~ LogNormal(log(1e-4), 1.5)
    α ~ LogNormal(log(3500), 0.5)
    p_mov ~ Beta(2, 20)

    Turing.@addlogprob! foi_loglik(data, in_zone, t₀, φ_hrz, φ_non, δ, β_duck, β, α, p_mov)
end

"""
Reparameterised full model. 8 estimated parameters.
Replaces (φ_non, β) with (ρ_non, κ) to reduce parameter competition:
- ρ_non = φ_non / φ_hrz  (non-HRZ spillover as fraction of HRZ)
- κ = β / φ_hrz  (spatial-to-spillover hazard ratio)
This ties spillover and spatial scales together, preventing the two from
competing freely during inference.
"""
Turing.@model function full_reparam_model(data::ModelData, in_zone::BitMatrix)
    t₀ ~ truncated(Normal(15, 5), 1, 44)
    φ_hrz ~ LogNormal(log(1e-3), 1.0)
    ρ_non ~ Beta(2, 10)                   # φ_non / φ_hrz ratio
    δ ~ Exponential(1/50)
    β_duck ~ Beta(2, 8)
    κ ~ LogNormal(log(0.1), 1.5)          # β / φ_hrz ratio
    α ~ LogNormal(log(3500), 0.5)
    p_mov ~ Beta(2, 20)

    # Derived parameters
    φ_non = ρ_non * φ_hrz
    β = κ * φ_hrz

    Turing.@addlogprob! foi_loglik(data, in_zone, t₀, φ_hrz, φ_non, δ, β_duck, β, α, p_mov)
end
