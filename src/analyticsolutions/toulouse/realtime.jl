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
toulouse_Gt(bath::AbstractFermionicNormalBath, t::Real; kwargs...) = fermionic_toulouse_Gt(bath.spectrum, t; μ=bath.μ, kwargs...)

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