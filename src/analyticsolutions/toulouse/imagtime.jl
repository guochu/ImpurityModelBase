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