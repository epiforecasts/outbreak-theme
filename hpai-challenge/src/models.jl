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
