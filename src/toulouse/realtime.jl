# real-time Green's function of Toulouse model

"""
    toulouse_Gw(spectrum::AbstractSpectrumFunction, ω::Real; ϵ_d::Real, μ::Real)

Retarded Green's function in the frequency axis
"""
function toulouse_Gw(spectrum::AbstractSpectrumFunction, ω::Real; ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8)
	f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
	1.0/(ω-ϵ_d-quadgkwrapper(bounded(ε -> f(ε)/(ω+μ-ε+im*δ), lb, ub)))
end
toulouse_Gw(bath::AbstractFermionicBath, ω::Real; kwargs...) = toulouse_Gw(bath.spectrum, ω; μ=bath.μ, kwargs...)

"""
    toulouse_Gt(spectrum::AbstractSpectrumFunction, t::Real; ϵ_d::Real, μ::Real, wmax::Real)

Retarded Green's function in the real time axis
"""
function toulouse_Gt(spectrum::AbstractSpectrumFunction, t::Real; ϵ_d::Real, μ::Real=0, wmax::Real=20., δ::Real=1.0e-8)
    A = quadgkwrapper(bounded(ω -> (toulouse_Gw(spectrum, ω; ϵ_d=ϵ_d, μ=μ)-1.0/(ω+im*δ))*exp(-im*ω*t), -wmax, wmax))
    return im*(A/(2π)-im)
end
toulouse_Gt(bath::AbstractFermionicBath, t::Real; kwargs...) = toulouse_Gt(bath.spectrum, t; μ=bath.μ, kwargs...)


function toulouse_Δw(spectrum::AbstractSpectrumFunction, ω::Real; δ::Real=1.0e-8)
    f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
    return quadgkwrapper(bounded(ϵ -> f(ϵ) / (ω - ϵ + im*δ), lb, ub))
end
toulouse_Jw(spectrum::AbstractSpectrumFunction, ω::Real; kwargs...) = -imag(toulouse_Δw(spectrum, ω; kwargs...)) / π