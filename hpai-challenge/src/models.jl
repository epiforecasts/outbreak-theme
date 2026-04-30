# ── 8. Turing model definitions ───────────────────────────────────────────────

using Turing
using Distributions

"""
Spillover-only model (Module A + C). 6 estimated parameters.
Captures early dynamics and HRZ contrast without farm-to-farm transmission.
Spillover profile is a Bateman function: rises from t₀ with rate σ, decays with rate δ.
"""
Turing.@model function spillover_model(data::ModelData, in_zone::BitMatrix)
    t₀ ~ truncated(Normal(15, 5), 1, 75)
    φ_hrz ~ LogNormal(log(1e-3), 1.0)
    φ_non ~ LogNormal(log(1e-4), 1.0)
    δ ~ Exponential(1/50)     # scale parameterisation: mean = 0.02/day
    β_duck ~ Beta(2, 8)
    σ ~ LogNormal(log(0.3), 1.0)  # spillover rise rate

    Turing.@addlogprob! foi_loglik_spillover(data, in_zone, t₀, φ_hrz, φ_non, δ, β_duck, σ)
end

"""
Full model (Module A + B + C). 9 estimated parameters.
Includes spillover, spatial transmission, and movement transmission.
Spillover profile is a Bateman function: rises from t₀ with rate σ, decays with rate δ.
"""
Turing.@model function full_model(data::ModelData, in_zone::BitMatrix)
    t₀ ~ truncated(Normal(15, 5), 1, 75)
    φ_hrz ~ LogNormal(log(1e-3), 1.0)
    φ_non ~ LogNormal(log(1e-4), 1.0)
    δ ~ Exponential(1/50)     # scale parameterisation: mean = 0.02/day
    β_duck ~ Beta(2, 8)
    σ ~ LogNormal(log(0.3), 1.0)  # spillover rise rate
    β ~ LogNormal(log(1e-4), 1.5)
    α ~ LogNormal(log(3500), 0.5)
    p_mov ~ Beta(2, 20)

    Turing.@addlogprob! foi_loglik(data, in_zone, t₀, φ_hrz, φ_non, δ, β_duck, σ, β, α, p_mov)
end
