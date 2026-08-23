
"""
    independentbosons_Gτ(spectrum::AbstractBoundedFunction, τ::Real; β, ϵ_d, U, bands, Δ)

Matsubara Green's function in the imaginary-time axis for the independent bosons model.

The independent bosons model describes a fermionic impurity (localized electrons)
coupled to a phonon heat reservoir (a bosonic bath of independent harmonic oscillators),
with `spectrum` the bath spectrum density characterizing the electron-phonon coupling
gₖ.

For `bands = 1` (a spinless, single-band fermion), the Hamiltonian reads

H = ϵ_d d†d + ∑ₖ ωₖ aₖ† aₖ + d†d ∑ₖ gₖ (aₖ + aₖ†),

where d† (d) creates (annihilates) the spinless electron on the impurity with on-site
energy ϵ_d, and aₖ† (aₖ) creates (annihilates) a phonon of frequency ωₖ.

For `bands = 2` (a spinful, two-band fermion), the Hamiltonian reads

H = ϵ_d ∑_σ d_σ† d_σ + U d_↑† d_↑ d_↓† d_↓ + ∑ₖ ωₖ aₖ† aₖ + ∑_σ d_σ† d_σ ∑ₖ gₖ (aₖ + aₖ†),

where σ ∈ {↑, ↓} labels the spin of the impurity fermion and U is the interaction
strength between the two spin species.

ϵ_d is the on-site energy of the localized electron
U is the interaction strength of the localized electron
bands is the number of bands (spin species) of the impurity fermion:
- `bands = 1` corresponds to the exact solution for a spinless (single-band) fermion
- `bands = 2` corresponds to the exact solution for a spinful (two-band) fermion,
  where U is the interaction strength between the two spin species
"""
function independentbosons_Gτ(spectrum::AbstractBoundedFunction, τ::Real; β::Real, ϵ_d::Real, U::Real=0, bands::Int=1, Δ::Real=_compute_Δ(spectrum))
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
function independentbosons_Gτ(spectrum::AbstractBoundedFunction; β::Real, Nτ::Int, ϵ_d::Real, δτ::Real=β/Nτ, U::Real=0, bands::Int=1)
    Δ = _compute_Δ(spectrum)
    # δτ = β / Nτ
    gτ = zeros(Float64, Nτ+1)
    for i in 1:Nτ+1
        τ = (i-1) * δτ
        tmp = independentbosons_Gτ(spectrum, τ, β=β, ϵ_d=ϵ_d, U=U, Δ=Δ, bands=bands)
        (abs(imag(tmp)) < 1.0e-8) || error("imaginary part of Gτ is too large")
        gτ[i] = real(tmp)
    end
    # gτ[end] = 1 - gτ[1]
    return gτ
end
