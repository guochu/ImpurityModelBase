# imaginary-time (Matsubara) Green's function of the Toulouse model
# (a single fermionic/bosonic impurity coupled to a quadratic normal bath)

"""
    fermionic_toulouse_Giw(spectrum::AbstractBoundedFunction, ω::Real; ϵ_d::Real, μ::Real)

Matsubara Green's function in the imaginary-frequency axis for the fermionic Toulouse model.
Its Hamiltonian reads

H = ϵ_d d†d + ∑ₖ εₖ cₖ† cₖ + ∑ₖ Vₖ (d† cₖ + cₖ† d),

where d† (d) creates (annihilates) an electron on the impurity with on-site energy ϵ_d,
cₖ† (cₖ) creates (annihilates) a bath electron of energy εₖ, and Vₖ is the impurity-bath
coupling strength, whose distribution is characterized by the bath spectral density
`spectrum`. The impurity Matsubara Green's function is exactly

G(iω) = 1 / (iω - ϵ_d - Δ(iω)),  Δ(iω) = ∫ dε J(ε)/(iω + μ - ε).

`ϵ_d` is the on-site energy of the localized electron and `μ` the chemical potential of the
bath; `β` is not a parameter as G(iω) is independent of β for the Toulouse model.
"""
function fermionic_toulouse_Giw(f::AbstractBoundedFunction, ω::Real; ϵ_d::Real, μ::Real=0)
    g(ε) = im*ω+μ-ε
    1.0/(im*ω-ϵ_d-quadgkwrapper(f / g))
end
function fermionic_toulouse_Giw(spectrum::AbstractBoundedFunction; β::Real, ϵ_d::Real, μ::Real=0, n::Int=1000)
    return [fermionic_toulouse_Giw(spectrum, x; ϵ_d=ϵ_d, μ=μ) for x in ifrequencies(β, n)]
end
"""
    toulouse_Giw(bath::AbstractFermionicNormalBath, ω::Real; ϵ_d::Real, kwargs...)

Matsubara Green's function of the fermionic Toulouse model for the impurity embedded
in the thermal fermionic bath `bath`; the bath chemical potential is taken from `bath`.
Delegates to `fermionic_toulouse_Giw`.
"""
toulouse_Giw(bath::AbstractFermionicNormalBath, ω::Real; kwargs...) = fermionic_toulouse_Giw(bath.spectrum, ω; μ=bath.μ, kwargs...)
function toulouse_Giw(bath::AbstractFermionicNormalBath; β::Real=bath.β, ϵ_d::Real, μ::Real=bath.μ, n::Int=1000)
    return fermionic_toulouse_Giw(bath.spectrum; β=β, ϵ_d=ϵ_d, μ=μ, n=n)
end


"""
    bosonic_toulouse_Giw(spectrum::AbstractBoundedFunction, ω::Real; ϵ_d::Real, μ::Real)

Matsubara Green's function in the imaginary-frequency axis for the bosonic Toulouse model
(a single bosonic mode coupled linearly to a noninteracting bosonic bath). For this
quadratic Hamiltonian the impurity Matsubara Green's function has the same single-particle
form as the fermionic one,

G(iω) = 1 / (iω - ϵ_d - Δ(iω)),  Δ(iω) = ∫ dε J(ε)/(iω + μ - ε),

with bosonic (even) Matsubara frequencies ω_n = 2πn/β. See [`fermionic_toulouse_Giw`](@ref)
for the parameter meanings.
"""
function bosonic_toulouse_Giw(f::AbstractBoundedFunction, ω::Real; ϵ_d::Real, μ::Real=0)
    g(ε) = im*ω+μ-ε
    1.0/(im*ω-ϵ_d-quadgkwrapper(f / g))
end
function bosonic_toulouse_Giw(spectrum::AbstractBoundedFunction; β::Real, ϵ_d::Real, μ::Real=0, n::Int=1000)
    return [bosonic_toulouse_Giw(spectrum, 2π*nk/β; ϵ_d=ϵ_d, μ=μ) for nk in -n:n]
end
toulouse_Giw(bath::AbstractBosonicNormalBath, ω::Real; kwargs...) = bosonic_toulouse_Giw(bath.spectrum, ω; μ=bath.μ, kwargs...)
function toulouse_Giw(bath::AbstractBosonicNormalBath; β::Real=bath.β, ϵ_d::Real, μ::Real=bath.μ, n::Int=1000)
    return bosonic_toulouse_Giw(bath.spectrum; β=β, ϵ_d=ϵ_d, μ=μ, n=n)
end


"""
    fermionic_toulouse_Gτ(spectrum::AbstractBoundedFunction, τ::Real; β::Real, ϵ_d::Real, μ::Real, n::Int)

Matsubara Green's function in the imaginary-time axis for the fermionic Toulouse model
(fermionic, anti-periodic, odd Matsubara frequencies ω_n = (2n+1)π/β). `n` controls the
number of Matsubara frequencies kept in the series.
"""
function fermionic_toulouse_Gτ(spectrum::AbstractBoundedFunction, τ::Real; β::Real, ϵ_d::Real, μ::Real=0., n::Int=1000)
    res = 0.0
    for ω in ifrequencies(β, n)
        res += (fermionic_toulouse_Giw(spectrum, ω; ϵ_d=ϵ_d, μ=μ)-1/(im*ω))*exp(-im*τ*ω)
    end
    res = -(res/β-0.5)
end
"""
    toulouse_Gτ(bath::AbstractFermionicNormalBath, τ::Real; ϵ_d::Real, kwargs...)

Imaginary-time Matsubara Green's function of the fermionic Toulouse model for the
impurity embedded in the thermal fermionic bath `bath`; the bath chemical potential
is taken from `bath`. Delegates to `fermionic_toulouse_Gτ`.
"""
toulouse_Gτ(bath::AbstractFermionicNormalBath, τ::Real; kwargs...) = fermionic_toulouse_Gτ(bath.spectrum, τ; β=bath.β, μ=bath.μ, kwargs...)
function fermionic_toulouse_Gτ(spectrum::AbstractBoundedFunction; β::Real, Nτ::Int, ϵ_d::Real, μ::Real=0., n::Int=1000)
    δτ = β / Nτ
    gτ = zeros(Float64, Nτ+1)
    for i in 1:Nτ
        τ = (i-1) * δτ
        tmp = fermionic_toulouse_Gτ(spectrum, τ, β=β, ϵ_d=ϵ_d, μ=μ, n=n)
        (abs(imag(tmp)) < 1.0e-8) || error("imaginary part of Gτ is too large")
        gτ[i] = real(tmp)
    end
    gτ[end] = 1 - gτ[1]
    return gτ
end
toulouse_Gτ(bath::AbstractFermionicNormalBath; Nτ::Int, kwargs...) = fermionic_toulouse_Gτ(bath.spectrum; β=bath.β, μ=bath.μ, Nτ=Nτ, kwargs...)

"""
    fermionic_toulouse_Giw_semicircular(ω::Real; ϵ_d::Real, μ::Real=0, t::Real=1.0)

Closed-form Matsubara Green's function of the fermionic Toulouse model for a
semicircular bath spectrum J(ε) = (2/πt²)√(t²-ε²) (normalized to unity with
half-bandwidth `t`, i.e. exactly the spectrum returned by [`semicircular`](@ref)).

The Hilbert transform of the semicircle is available analytically,

Δ(iω) = 2(w - √(w-t)√(w+t))/t²,  w = μ + iω,  branch √(w-t)√(w+t) ~ w at |w|→∞,

so the impurity Matsubara Green's function

G(iω) = 1 / (iω - ϵ_d - Δ(iω))

is evaluated in closed form without any numerical integration (and without any
broadening δ). Equivalent to calling `fermionic_toulouse_Giw` with
`semicircular(t)` but exact and much faster.
"""
function fermionic_toulouse_Giw_semicircular(ω::Real; ϵ_d::Real, μ::Real=0, t::Real=1.0)
    w = μ + im * ω
    # retarded branch s ~ w at |w|→∞: the product of principal square roots
    # (plain sqrt(w²-t²) picks the wrong sign for Re(w)<0)
    s = sqrt(w - t) * sqrt(w + t)
    Δ = 2 * (w - s) / t^2
    return 1 / (im * ω - ϵ_d - Δ)
end

function fermionic_toulouse_Giw_semicircular(; β::Real, ϵ_d::Real, μ::Real=0, t::Real=1.0, n::Int=1000)
    return [fermionic_toulouse_Giw_semicircular(x; ϵ_d=ϵ_d, μ=μ, t=t) for x in ifrequencies(β, n)]
end

"""
    fermionic_toulouse_Gτ_semicircular(τ::Real; β::Real, ϵ_d::Real, μ::Real=0, t::Real=1.0)

Closed-form Matsubara Green's function in the imaginary-time axis for the fermionic
Toulouse model with the semicircular bath spectrum `semicircular(t)` (normalized to
unity, half-bandwidth `t`); see [`fermionic_toulouse_Giw_semicircular`](@ref).

There is no elementary closed form at finite β (the Fermi factors prevent it), but the
spectral representation gives an exact single-quadrature expression with everything in
closed form. With ϵ̃ = ϵ_d + μ, γ(ω) = √(t²-(ω+μ)²) and the numerically stable kernel
k(x) = e^{-x(τ-β/2)}/(2cosh(xβ/2)) = e^{-xτ}/(1+e^{-βx}),

G(τ) = Σⱼ Zⱼ k(rⱼ-μ) + ∫_{-t-μ}^{t-μ} dω γ(ω+μ)/(π t²·(a(ω)² + b(ω)²)) · k(ω)

where a(ω) = ω - ϵ_d - 2(ω+μ)/t², b(ω) = 2γ(ω)/t², and the bound states (rⱼ, Zⱼ) are
exactly those of the retarded function (see [`fermionic_toulouse_Gt_semicircular`](@ref)):
up to two real poles outside the bath band solving (r-ϵ̃)(r+sign(r)√(r²-t²)) = 2.
The result is antiperiodically extended for τ outside [0, β].
"""
function fermionic_toulouse_Gτ_semicircular(τ::Real; β::Real, ϵ_d::Real, μ::Real=0, t::Real=1.0)
    # antiperiodic reduction: τ = s·β + τr with τr ∈ [0, β), G(τ) = (-1)^s G(τr)
    s, τr = divrem(τ, β)
    τr += β * (τr < 0)          # divrem rounds toward zero: fix sign of remainder

    # bound-state poles (identical to the retarded case), in the bath frame w = ω + μ
    ϵ_eff = ϵ_d + μ
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
        sv = sign(u) * sqrt(u^2 - t^2)
        abs((u - ϵ_eff) * (u + sv) - 2) < 1.0e-6 * (1 + abs(u)) && push!(boundstates, (u, D2 / ((D2 - 2) + 2 * u / sv)))
    end

    # numerically stable Fermi kernel k(x) = e^{-xτ}/(1+e^{-βx})
    kernel(x) = exp(-x * (τr - β / 2)) / (2 * cosh(x * β / 2))

    # bound-state contribution (package convention: G(τ) = +⟨T̂ d(τ)d†(0)⟩)
	Gτ = sum(Z * kernel(r - μ) for (r, Z) in boundstates; init=0.0)

	# branch-band contribution from the closed-form spectral density
	# A(ω) = b(ω)/(π(a(ω)²+b(ω)²)), a = ω-ϵ_d-2(ω+μ)/t², b = 2√(t²-(ω+μ)²)/t²
	cut, _ = quadgk(-t - μ, t - μ; rtol=1.0e-10, atol=1.0e-12) do ω
		γ = sqrt(max(t^2 - (ω + μ)^2, 0.0))
		b = 2 * γ / t^2
		a = ω - ϵ_d - 2 * (ω + μ) / t^2
		b / (π * (a^2 + b^2)) * kernel(ω)
	end
	Gτ += cut
	return isodd(s) ? -Gτ : Gτ
end

# function toulouse_Δiw(spectrum::AbstractBoundedFunction; β::Real, n::Int=1000)
#     f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
#     ff(ω) = quadgkwrapper(bounded(ϵ -> f(ϵ) / (im*ω - ϵ), lb, ub))

#     return [ff((2*n-1)*π/β) for n in -n:n+1]
# end

"""
    bosonic_toulouse_Gτ(spectrum::AbstractBoundedFunction, τ::Real; β::Real, ϵ_d::Real, μ::Real, n::Int)

Matsubara Green's function in the imaginary-time axis for the bosonic Toulouse model.
For a normal (quadratic) bosonic system only the *single-particle* Green's function is
non-singular; it is obtained from `bosonic_toulouse_Giw` by the bosonic Matsubara sum
(bosonic/even frequencies ω_n = 2πn/β),

G(τ) = -1/β [ G(0) + Σ_{n≠0} (G(iω_n) - 1/(iω_n)) e^{-iω_n τ} + (τ - β/2) ],

where the `G(0)` and `τ - β/2` terms are the regularized n = 0 (static) contribution of the
1/(iω) tail of the Green's function. `n` controls the number of Matsubara frequencies kept.
"""
function bosonic_toulouse_Gτ(spectrum::AbstractBoundedFunction, τ::Real; β::Real, ϵ_d::Real, μ::Real=0., n::Int=1000)
    G0 = bosonic_toulouse_Giw(spectrum, 0.0; ϵ_d=ϵ_d, μ=μ)
    res = 0.0
    for nk in -n:n
        nk == 0 && continue
        ω = 2π*nk/β
        res += (bosonic_toulouse_Giw(spectrum, ω; ϵ_d=ϵ_d, μ=μ) - 1/(im*ω)) * exp(-im*τ*ω)
    end
    return -(res + G0 + (τ - β/2))/β
end
toulouse_Gτ(bath::AbstractBosonicNormalBath, τ::Real; kwargs...) = bosonic_toulouse_Gτ(bath.spectrum, τ; β=bath.β, μ=bath.μ, kwargs...)
function bosonic_toulouse_Gτ(spectrum::AbstractBoundedFunction; β::Real, Nτ::Int, ϵ_d::Real, μ::Real=0., n::Int=1000)
    δτ = β / Nτ
    gτ = zeros(Float64, Nτ+1)
    for i in 1:Nτ+1
        τ = (i-1) * δτ
        tmp = bosonic_toulouse_Gτ(spectrum, τ, β=β, ϵ_d=ϵ_d, μ=μ, n=n)
        gτ[i] = real(tmp)
    end
    return gτ
end
toulouse_Gτ(bath::AbstractBosonicNormalBath; Nτ::Int, kwargs...) = bosonic_toulouse_Gτ(bath.spectrum; β=bath.β, μ=bath.μ, Nτ=Nτ, kwargs...)


"""
    toulouse_Δiw(f::AbstractBoundedFunction; β::Real, n::Int=1000)

The hybridization function in the imaginary-frequency axis for the Toulouse model
(independent of the fermionic/bosonic statistics of the impurity).

f is the bath spectrum density
"""
function toulouse_Δiw(f::AbstractBoundedFunction; β::Real, n::Int=1000)
    function ff(ω)
        g(ϵ) = im*ω - ϵ
        return quadgkwrapper(f / g)
    end 

    return [ff(x) for x in ifrequencies(β, n)]
end
toulouse_Δiw(bath::AbstractFermionicNormalBath; n::Int=1000) = toulouse_Δiw(bath.spectrum, β=bath.β, n=n)
toulouse_Δiw(bath::AbstractBosonicNormalBath; n::Int=1000) = toulouse_Δiw(bath.spectrum, β=bath.β, n=n)

# function toulouse_Δτ(spectrum::AbstractBoundedFunction; β::Real, N::Int)
#     f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
#     δτ = β / N
#     ff(τ) = quadgkwrapper(bounede(ϵ -> -f(ϵ) / (exp(-ϵ*τ) / (1+exp(-β*ϵ)) ), lb, ub))
#     return [ff(i*δτ) for i in 0:N]
# end

"""
    toulouse_Δτ(f::AbstractBoundedFunction; β::Real, Nτ::Int)

The hybridization function in the imaginary-time axis for the Toulouse model
(independent of the fermionic/bosonic statistics of the impurity).

f is the bath spectrum density
"""
function toulouse_Δτ(f::AbstractBoundedFunction; β::Real, Nτ::Int)
    δτ = β / Nτ
    function ff(τ)
        g(ϵ) = -exp(-ϵ*τ) / (1+exp(-β*ϵ)) 
        return quadgkwrapper(f / g)
    end 
    return [ff(i*δτ) for i in 0:Nτ]
end
toulouse_Δτ(bath::AbstractFermionicNormalBath; Nτ::Int) = toulouse_Δτ(bath.spectrum, β=bath.β, Nτ=Nτ)
toulouse_Δτ(bath::AbstractBosonicNormalBath; Nτ::Int) = toulouse_Δτ(bath.spectrum, β=bath.β, Nτ=Nτ)