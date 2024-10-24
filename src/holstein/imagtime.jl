
"""
    independentbosons_Gτ(spectrum::SpectrumFunction, τ::Real; β, ϵ_d, n₀, U, Δ)

Matsubara Green's function in the imaginary-time axis
"""
function independentbosons_Gτ(spectrum::SpectrumFunction, τ::Real; β::Real, ϵ_d::Real, U::Real=0, Δ::Real=_compute_Δ(spectrum))
    μ = -ϵ_d
    return (U == 0) ? -(1-freefermion_occupation(β, μ))*exp(τ*(μ+Δ))*_exponent_f(spectrum, τ, β) : -_interact(τ, β, μ, Δ, U)*_exponent_f(spectrum, τ, β)
end

