# real-time Green's function of Toulouse model


"""
    toulouse_Gw(spectrum::AbstractSpectrumFunction, ω::Real; ϵ_d::Real, μ::Real)

Retarded Green's function in the real-frequency axis for the Toulouse model

f is the bath spectrum density
ϵ_d is the on-site energy of the localized electron
μ is the chemical potental of the bath
β is not a parameter as Gw is independent of β for the Toulouse model
"""
function toulouse_Gw(f::AbstractSpectrumFunction, ω::Real; ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8)
    g(ϵ) = ω+μ-ϵ+im*δ
	return 1.0/(ω+im*δ-ϵ_d-quadgkwrapper(f/g))
end
toulouse_Gw(bath::AbstractFermionicBath, ω::Real; kwargs...) = toulouse_Gw(bath.spectrum, ω; μ=bath.μ, kwargs...)

"""
    toulouse_Gt(spectrum::AbstractSpectrumFunction, t::Real; ϵ_d, μ, wmax, wmin, δ)

Retarded Green's function in the real-time axis for the Toulouse model

f is the bath spectrum density
ϵ_d is the on-site energy of the localized electron
μ is the chemical potental of the bath

"""
function toulouse_Gt(spectrum::AbstractSpectrumFunction, t::Real; ϵ_d::Real, μ::Real=0, wmax::Real=20., wmin::Real=-wmax, δ::Real=1.0e-8)
    A = quadgkwrapper(bounded(ω -> (toulouse_Gw(spectrum, ω; ϵ_d=ϵ_d, μ=μ, δ=δ)-1.0/(ω+im*δ))*exp(-im*ω*t), wmin, wmax))
    return im*(A/(2π)-im)
end
toulouse_Gt(bath::AbstractFermionicBath, t::Real; kwargs...) = toulouse_Gt(bath.spectrum, t; μ=bath.μ, kwargs...)


"""
    toulouse_Δw(f::AbstractSpectrumFunction, ω::Real; δ)

The hybridization function in the real-frequency axis for the Toulouse model

f is the bath spectrum density
"""
function toulouse_Δw(spectrum::AbstractSpectrumFunction, ω::Real; δ::Real=1.0e-8)
    f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
    g(ϵ) = ω - ϵ + im*δ
    return quadgkwrapper(f / g)
end

# the relation between Δw and Jw for the Toulouse model
# toulouse_Jw(spectrum::AbstractSpectrumFunction, ω::Real; kwargs...) = -imag(toulouse_Δw(spectrum, ω; kwargs...)) / π