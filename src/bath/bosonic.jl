# bosonic baths
struct BosonicBath{F <: AbstractSpectrumFunction} <: AbstractBosonicBath
	f::F
	β::Float64
	μ::Float64
end
BosonicBath(f::AbstractSpectrumFunction; β::Real, μ::Real=0) = BosonicBath(f, convert(Float64, β), convert(Float64, μ))


struct BosonicVacuum{F <: AbstractSpectrumFunction} <: AbstractBosonicBath
	f::F
	μ::Float64	
end
BosonicVacuum(f::AbstractSpectrumFunction; μ::Real=0) = BosonicVacuum(f, convert(Float64, μ))


function boseeinstein(β::Real, μ::Real, ϵ::Real)
	n_k = 0.
	(ϵ > μ) || throw(ArgumentError("energy must be larger than μ"))
	if β == Inf
		return 0.
	else
		return 1 / (exp(β * (ϵ - μ)) - 1)
	end
end

"""
	thermal_occupation(bath::AbstractBosonicBath, energy::Real)
	how to define this function if μ > 0
"""
thermaloccupation(bath::AbstractBosonicBath, ϵ::Real) = boseeinstein(bath.β, bath.μ, ϵ)

"""
	bosonic_bath(f; β::Real=Inf) 
	f should support y = f(x) with y positive
"""
bosonicbath(f::AbstractSpectrumFunction; β::Real=Inf, μ::Real=0) = (β == Inf) ? BosonicVacuum(f, μ=μ) : BosonicBath(f, β=β, μ=μ)


# function bosonicbath(f; β::Real=Inf, μ::Real=0) 
# 	y = f(0.)
# 	(y >= 0.) || throw(ArgumentError("spectrum function should return positive numbers."))
	
# 	(β == Inf) && return BosonicVacuum(f, convert(Float64, μ))
# 	(β >= 0.) || throw(ArgumentError("temperature should not be negative."))
# 	return BosonicBath(f, convert(Float64, β), convert(Float64, μ))
# end