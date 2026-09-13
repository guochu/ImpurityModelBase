# real-time Green's function of the Toulouse model
# (a single fermionic/bosonic impurity coupled to a quadratic normal bath)

"""
    fermionic_toulouse_Gw(spectrum::AbstractBoundedFunction, ω::Real; ϵ_d::Real, μ::Real, δ::Real)

Retarded Green's function in the real-frequency axis for the fermionic Toulouse model
(a single localized electron coupled to a noninteracting fermionic bath). Its Hamiltonian
is the resonant-level (interacting-impurity-free) model

H = ϵ_d d†d + ∑ₖ εₖ cₖ† cₖ + ∑ₖ Vₖ (d† cₖ + cₖ† d),

for which the impurity retarded Green's function is exactly

G(ω) = 1 / (ω + iδ - ϵ_d - Δ(ω)),  Δ(ω) = ∫ dε J(ε)/(ω + μ - ε + iδ),

where `spectrum` provides the bath spectral density J(ε), `ϵ_d` is the impurity on-site
energy, `μ` the bath chemical potential and `δ` an infinitesimal broadening. `β` is not a
parameter as G(ω) is independent of β.
"""
function fermionic_toulouse_Gw(f::AbstractBoundedFunction, ω::Real; ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8)
    g(ϵ) = ω+μ-ϵ+im*δ
	return 1.0/(ω+im*δ-ϵ_d-quadgkwrapper(f/g))
end
"""
    toulouse_Gw(bath::AbstractFermionicNormalBath, ω::Real; ϵ_d::Real, kwargs...)

Retarded Green's function of the fermionic Toulouse model for the impurity embedded
in the thermal fermionic bath `bath`; the bath chemical potential is taken from
`bath`. Delegates to `fermionic_toulouse_Gw` with the bath's spectrum and `μ`.
"""
toulouse_Gw(bath::AbstractFermionicNormalBath, ω::Real; kwargs...) = fermionic_toulouse_Gw(bath.spectrum, ω; μ=bath.μ, kwargs...)

"""
    bosonic_toulouse_Gw(spectrum::AbstractBoundedFunction, ω::Real; ϵ_d::Real, μ::Real, δ::Real)

Retarded Green's function in the real-frequency axis for the bosonic Toulouse model
(a single bosonic mode coupled linearly to a noninteracting bosonic bath). For this
quadratic Hamiltonian the impurity retarded Green's function has the same single-particle
form as the fermionic one,

G(ω) = 1 / (ω + iδ - ϵ_d - Δ(ω)),  Δ(ω) = ∫ dε J(ε)/(ω + μ - ε + iδ),

the particle statistics only enter the temperature dependence of the thermal (greater/
lesser) Green's functions. See [`fermionic_toulouse_Gw`](@ref) for the parameter meanings.
"""
bosonic_toulouse_Gw(f::AbstractBoundedFunction, ω::Real; kwargs...) = fermionic_toulouse_Gw(f, ω; kwargs...)
toulouse_Gw(bath::AbstractBosonicNormalBath, ω::Real; kwargs...) = fermionic_toulouse_Gw(bath.spectrum, ω; μ=bath.μ, kwargs...)

"""
    fermionic_toulouse_Gt(spectrum::AbstractBoundedFunction, t::Real; ϵ_d, μ, wmax, wmin, δ)

Retarded Green's function in the real-time axis for the fermionic Toulouse model, obtained
by Fourier transforming `fermionic_toulouse_Gw` back to the time domain.

`spectrum` is the bath spectral density, `ϵ_d` the impurity on-site energy, `μ` the bath
chemical potential; `wmax`/`wmin` bound the frequency integration window.
"""
function fermionic_toulouse_Gt(spectrum::AbstractBoundedFunction, t::Real; ϵ_d::Real, μ::Real=0, wmax::Real=20., wmin::Real=-wmax, δ::Real=1.0e-8)
    A = quadgkwrapper(bounded(ω -> (fermionic_toulouse_Gw(spectrum, ω; ϵ_d=ϵ_d, μ=μ, δ=δ)-1.0/(ω+im*δ))*exp(-im*ω*t), wmin, wmax))
    return A/(2π)-im
end
"""
    toulouse_Gt(bath::AbstractFermionicNormalBath, t::Real; ϵ_d::Real, kwargs...)

Real-time retarded Green's function of the fermionic Toulouse model for the impurity
embedded in the thermal fermionic bath `bath`; the bath chemical potential is taken
from `bath`. Delegates to `fermionic_toulouse_Gt`.
"""
toulouse_Gt(bath::AbstractFermionicNormalBath, t::Real; kwargs...) = fermionic_toulouse_Gt(bath.spectrum, t; μ=bath.μ, kwargs...)

"""
    fermionic_toulouse_Gw_semicircular(ω::Real; ϵ_d::Real, μ::Real=0, t::Real=1.0, δ::Real=1.0e-8)

Closed-form retarded Green's function of the fermionic Toulouse model for a
semicircular bath spectrum J(ε) = (2/πt²)√(t²-ε²) (normalized to unity with
half-bandwidth `t`, i.e. exactly the spectrum returned by [`semicircular`](@ref)).

The Hilbert transform of the semicircle is available analytically,

Δ(z) = 2(w - √(w²-t²))/t²,  w = z,  branch √(w²-t²) ~ w at |w|→∞,

so the impurity Green's function

G(ω) = 1 / (ω + iδ - ϵ_d - Δ(ω + μ + iδ))

is evaluated in closed form without any numerical integration. Equivalent to
calling `fermionic_toulouse_Gw` with `semicircular(t)` but exact and much faster.
"""
function fermionic_toulouse_Gw_semicircular(ω::Real; ϵ_d::Real, μ::Real=0, t::Real=1.0, δ::Real=1.0e-8)
	w = ω + μ + im * δ
	# retarded branch s ~ w at |w|→∞: the product of principal square roots
	# (plain sqrt(w²-t²) picks the wrong sign for Re(w)<0)
	s = sqrt(w - t) * sqrt(w + t)
	Δ = 2 * (w - s) / t^2
	return 1 / (ω + im * δ - ϵ_d - Δ)
end

"""
    fermionic_toulouse_Gt_semicircular(τ::Real; ϵ_d::Real, μ::Real=0, t::Real=1.0)

Closed-form retarded Green's function in the real-time axis for the fermionic
Toulouse model with the semicircular bath spectrum `semicircular(t)` (normalized
to unity, half-bandwidth `t`); see [`fermionic_toulouse_Gw_semicircular`](@ref).

With ϵ̃ = ϵ_d + μ and γ(ω) = √(t²-ω²), the exact result for τ > 0 is the sum of the
bound-state poles and a branch-cut integral (the discontinuity of G across the bath band):

G(τ) = e^{iμτ} { -i Σⱼ Zⱼ e^{-irⱼτ} - (2i/π) ∫₋ₜᵗ e^{-iωτ} γ(ω)/|((ω-ϵ̃)(ω+iγ(ω))-2)|² dω }

The bound states are the real solutions r (|r|>t) of the implicit equation

(r-ϵ̃)(r + sign(r)√(r²-t²)) = 2,

each with residue Z = t²/((t²-2) + 2r/√(r²-t²)). There can be up to two bound states:
a single one when the level is far outside the band (|ϵ̃| > t), and — for sufficiently
narrow bands (t < √2) — a particle-hole pair even at ϵ̃ = 0. In the special case
ϵ̃ = 0, t = 2 (where 𝒢 coincides with Δ) the result reduces to the well-known
G(τ) = -iθ(τ) 2J₁(tτ)/(tτ), and the sum rule G(0⁺) = -i holds in general.
The function is causal: G(τ<0) = 0.
"""
function fermionic_toulouse_Gt_semicircular(τ::Real; ϵ_d::Real, μ::Real=0, t::Real=1.0)
	ϵ_eff = ϵ_d + μ
	τ < 0 && return zero(ComplexF64)            # retarded causality
	τ == 0 && return -im                        # sum rule: ∫A(ω)dω = 1

	# bound-state poles: real roots of (r-ϵ_eff)(r+s) = 2 outside the band.
	# They are found among the roots of the quadratic obtained by squaring
	# (spurious roots are removed by the unsquared-equation check).
	D2 = float(t)^2
	cands = Float64[]
	if isapprox(D2, 4.0; rtol=0.0, atol=1.0e-12)
		ϵ_eff == 0 || push!(cands, ϵ_eff + 1 / ϵ_eff)
	elseif ϵ_eff^2 + 4.0 - D2 >= 0
		E = sqrt(ϵ_eff^2 + 4.0 - D2)
		append!(cands, ((ϵ_eff * (D2 - 2) + 2E) / (D2 - 4), (ϵ_eff * (D2 - 2) - 2E) / (D2 - 4)))
	end
	boundstates = Tuple{Float64,Float64}[]      # (r, Z)
	for u in cands
		abs(u) > t || continue
		s = sign(u) * sqrt(u^2 - t^2)
		abs((u - ϵ_eff) * (u + s) - 2) < 1.0e-6 * (1 + abs(u)) && push!(boundstates, (u, D2 / ((D2 - 2) + 2 * u / s)))
	end

	# branch-cut contribution: -(2i/π)∫ e^{-iωτ} γ(ω)/|H₊(ω)|² dω,
	# H₊(ω) = (ω-ϵ_eff)(ω+iγ(ω)) - 2 (|H₊|² = ((ω-ϵ_eff)ω-2)² + ((ω-ϵ_eff)γ)²)
	cut, _ = quadgk(-t, t; rtol=1.0e-11, atol=1.0e-13) do ω
		γ = sqrt(max(t^2 - ω^2, 0.0))
		Hre = (ω - ϵ_eff) * ω - 2
		Him = (ω - ϵ_eff) * γ
		γ * exp(-im * ω * τ) / (Hre^2 + Him^2)
	end

	Gτ = -im * sum((Z * cis(-r * τ) for (r, Z) in boundstates); init=0.0im) - (2im / π) * cut
	return exp(im * μ * τ) * Gτ
end

"""
    bosonic_toulouse_Gt(spectrum::AbstractBoundedFunction, t::Real; ϵ_d, μ, wmax, wmin, δ)

Retarded Green's function in the real-time axis for the bosonic Toulouse model; see
[`bosonic_toulouse_Gw`](@ref) and [`fermionic_toulouse_Gt`](@ref).
"""
function bosonic_toulouse_Gt(spectrum::AbstractBoundedFunction, t::Real; ϵ_d::Real, μ::Real=0, wmax::Real=20., wmin::Real=-wmax, δ::Real=1.0e-8)
    A = quadgkwrapper(bounded(ω -> (bosonic_toulouse_Gw(spectrum, ω; ϵ_d=ϵ_d, μ=μ, δ=δ)-1.0/(ω+im*δ))*exp(-im*ω*t), wmin, wmax))
    return A/(2π)-im
end
toulouse_Gt(bath::AbstractBosonicNormalBath, t::Real; kwargs...) = bosonic_toulouse_Gt(bath.spectrum, t; μ=bath.μ, kwargs...)


"""
    toulouse_Δw(f::AbstractBoundedFunction, ω::Real; δ)

The hybridization function in the real-frequency axis for the Toulouse model

f is the bath spectrum density
"""
function toulouse_Δw(f::AbstractBoundedFunction, ω::Real; δ::Real=1.0e-8)
    g(ϵ) = ω - ϵ + im*δ
    return quadgkwrapper(f / g)
end

# the relation between Δw and Jw for the Toulouse model
# toulouse_Jw(spectrum::AbstractBoundedFunction, ω::Real; kwargs...) = -imag(toulouse_Δw(spectrum, ω; kwargs...)) / π