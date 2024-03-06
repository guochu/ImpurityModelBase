# real-time Green's function of Toulouse model

"""
    toulouse_Gw(spectrum::SpectrumFunction, ω::Real; ϵ_d::Real, μ::Real)

Retarded Green's function in the frequency axis
"""
function toulouse_Gw(spectrum::SpectrumFunction, ω::Real; ϵ_d::Real, μ::Real=0)
	δ = 1e-8
	f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
	1.0/(ω-ϵ_d-quadgk(ε -> f(ε)/(ω+μ-ε+im*δ), lb, ub)[1])
end
toulouse_Gw(bath::AbstractFermionicBath, ω::Real; kwargs...) = toulouse_Gw(bath.spectrum, ω; μ=bath.μ, kwargs...)

"""
    toulouse_Gt(spectrum::SpectrumFunction, t::Real; ϵ_d::Real, μ::Real, wmax::Real)

Retarded Green's function in the real time axis
"""
function toulouse_Gt(spectrum::SpectrumFunction, t::Real; ϵ_d::Real, μ::Real=0, wmax::Real=20.)
    δ = 1e-8
    A = quadgk(ω -> (toulouse_Gw(spectrum, ω; ϵ_d=ϵ_d, μ=μ)-1.0/(ω+im*δ))*exp(-im*ω*t), -wmax, wmax)[1]
    im*(A/(2π)-im)
end
toulouse_Gt(bath::AbstractFermionicBath, t::Real; kwargs...) = toulouse_Gt(bath.spectrum, t; μ=bath.μ, kwargs...)

# imaginary-time Green's function of Toulouse model

"""
    toulouse_Giw(spectrum::SpectrumFunction, ω::Real; ϵ_d::Real, μ::Real)

Matsubara Green's function in the frequency axis
"""
function toulouse_Giw(spectrum::SpectrumFunction, ω::Real; ϵ_d::Real, μ::Real=0)
	f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
    1.0/(im*ω-ϵ_d-quadgk(ε -> f(ε)/(im*ω+μ-ε), lb, ub)[1])
end
toulouse_Giw(bath::AbstractFermionicBath, ω::Real; kwargs...) = toulouse_Giw(bath.spectrum, ω; μ=bath.μ, kwargs...)

"""
    toulouse_Gτ(spectrum::SpectrumFunction, τ::Real; β::Real, ϵ_d::Real, μ::Real, nmax::Int)

Matsubara Green's function in the imaginary-time axis
"""
function toulouse_Gτ(spectrum::SpectrumFunction, τ::Real; β::Real, ϵ_d::Real, μ::Real=0., nmax::Int=1000)
    res = 0.0
    for n = -nmax:nmax+1
        ω = (2n-1)*π/β
        res += (toulouse_Giw(spectrum, ω; ϵ_d=ϵ_d, μ=μ)-1/(im*ω))*exp(-im*τ*ω)
    end
    res = -(res/β-0.5)
end
toulouse_Gτ(bath::AbstractFermionicBath, τ::Real; kwargs...) = toulouse_Gτ(bath.spectrum, τ; β=bath.β, μ=bath.μ, kwargs...)