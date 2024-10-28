
"""
    independentbosons_Gτ(spectrum::AbstractSpectrumFunction, τ::Real; β, ϵ_d, n₀, U, Δ)

Matsubara Green's function in the imaginary-time axis
"""
function independentbosons_Gτ(spectrum::AbstractSpectrumFunction, τ::Real; β::Real, ϵ_d::Real, U::Real=0, nbands::Int=1, Δ::Real=_compute_Δ(spectrum))
    (nbands in (1, 2)) || throw(ArgumentError("nbands must be 1 or 2"))
    μ′ = -ϵ_d + Δ
    if nbands == 1
        (U == 0) || println("nonzero U=$(U) ignored for nbands=1")
        return (1-freefermion_occupation(β, μ′))*exp(τ*μ′)*_exponent_f(spectrum, τ, β)
    else
        return _interact(τ, β, -ϵ_d, Δ, U)*_exponent_f(spectrum, τ, β)
    end
end
function independentbosons_Gτ(spectrum::AbstractSpectrumFunction; β::Real, N::Int, ϵ_d::Real, U::Real=0, nbands::Int=1)
    Δ = _compute_Δ(spectrum)
    δτ = β / N
    gτ = zeros(Float64, N+1)
    for i in 1:N
        τ = (i-1) * δτ
        tmp = independentbosons_Gτ(spectrum, τ, β=β, ϵ_d=ϵ_d, U=U, Δ=Δ, nbands=nbands)
        (abs(imag(tmp)) < 1.0e-8) || error("imaginary part of Gτ is too large")
        gτ[i] = real(tmp)
    end
    gτ[end] = 1 - gτ[1]
    return gτ
end
