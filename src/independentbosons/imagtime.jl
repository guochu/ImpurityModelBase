
"""
    independentbosons_Gτ(spectrum::SpectrumFunction, τ::Real; β, ϵ_d, n₀, U, Δ)

Matsubara Green's function in the imaginary-time axis
"""
function independentbosons_Gτ(spectrum::SpectrumFunction, τ::Real; β::Real, ϵ_d::Real, U::Real=0, Δ::Real=_compute_Δ(spectrum))
    μ = -ϵ_d
    return (U == 0) ? (1-freefermion_occupation(β, μ))*exp(τ*(μ+Δ))*_exponent_f(spectrum, τ, β) : _interact(τ, β, μ, Δ, U)*_exponent_f(spectrum, τ, β)
end
function independentbosons_Gτ(spectrum::SpectrumFunction; β::Real, N::Int, ϵ_d::Real, U::Real=0)
    Δ = _compute_Δ(spectrum)
    δτ = β / N
    gτ = zeros(Float64, N+1)
    for i in 1:N
        τ = (i-1) * δτ
        tmp = independentbosons_Gτ(spectrum, τ, β=β, ϵ_d=ϵ_d, U=U, Δ=Δ)
        (abs(imag(tmp)) < 1.0e-8) || error("imaginary part of Gτ is too large")
        gτ[i] = real(tmp)
    end
    gτ[end] = 1 - gτ[1]
    return gτ
end
