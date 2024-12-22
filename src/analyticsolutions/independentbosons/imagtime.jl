
"""
    independentbosons_Gτ(spectrum::AbstractSpectrumFunction, τ::Real; β, ϵ_d, U, bands, Δ)

Matsubara Green's function in the imaginary-time axis for the independent bosons model

ϵ_d is the on-site energy of the localized electron
U is the interaction strength of the localized electron
bands is the total number of bands, which can be 1 (noninteracting) 
or 2 (interacting)
"""
function independentbosons_Gτ(spectrum::AbstractSpectrumFunction, τ::Real; β::Real, ϵ_d::Real, U::Real=0, bands::Int=1, Δ::Real=_compute_Δ(spectrum))
    (bands in (1, 2)) || throw(ArgumentError("bands must be 1 or 2"))
    μ′ = -ϵ_d + Δ
    if bands == 1
        (U == 0) || println("nonzero U=$(U) ignored for bands=1")
        return freefermion_Gτ(τ, β=β, μ=μ′)*_exponent_f(spectrum, τ, β)
    else
        U′ = U - 2Δ
        return fermion_Gτ(τ, β=β, μ=μ′, U=U′)*_exponent_f(spectrum, τ, β)
    end
end
function independentbosons_Gτ(spectrum::AbstractSpectrumFunction; β::Real, Nτ::Int, ϵ_d::Real, U::Real=0, bands::Int=1)
    Δ = _compute_Δ(spectrum)
    δτ = β / Nτ
    gτ = zeros(Float64, Nτ+1)
    for i in 1:Nτ
        τ = (i-1) * δτ
        tmp = independentbosons_Gτ(spectrum, τ, β=β, ϵ_d=ϵ_d, U=U, Δ=Δ, bands=bands)
        (abs(imag(tmp)) < 1.0e-8) || error("imaginary part of Gτ is too large")
        gτ[i] = real(tmp)
    end
    gτ[end] = 1 - gτ[1]
    return gτ
end
