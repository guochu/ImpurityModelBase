"""
	struct BosonicBath{F <: AbstractSpectrumFunction}

Bosonic bath container, includes a bath spectrum density,
the inverse temperature β and the chemical potential μ
"""
struct BosonicBath{F <: AbstractSpectrumFunction} <: AbstractBosonicBath
	f::F
	β::Float64
	μ::Float64
end
BosonicBath(f::AbstractSpectrumFunction; β::Real, μ::Real=0) = BosonicBath(f, convert(Float64, β), convert(Float64, μ))

"""
	struct BosonicVacuum{F <: AbstractSpectrumFunction}

Bosonic bath container, includes a bath spectrum density,
the chemical potential μ
the inverse temperature β=Inf
"""
struct BosonicVacuum{F <: AbstractSpectrumFunction} <: AbstractBosonicBath
	f::F
	μ::Float64	
end
BosonicVacuum(f::AbstractSpectrumFunction; μ::Real=0) = BosonicVacuum(f, convert(Float64, μ))


"""
	boseeinstein(β, μ, ϵ)

Boson-einstein distribution for a bosonic bath 
with β, μ at energy ϵ
"""
function boseeinstein(β::Real, μ::Real, ϵ::Real)
	n_k = 0.
	(ϵ > μ) || throw(ArgumentError("energy must be larger than μ"))
	return 1 / (exp(safe_mult(β, ϵ - μ)) - 1)
	# if β == Inf
	# 	return 0.
	# else
	# 	return 1 / (exp(β * (ϵ - μ)) - 1)
	# end
end

"""
	thermaloccupation(bath::AbstractBosonicBath, energy::Real)

return n(ϵ)
"""
thermaloccupation(bath::AbstractBosonicBath, ϵ::Real) = boseeinstein(bath.β, bath.μ, ϵ)

"""
	bosonicbath(f; β, μ) 

Return a bosonic bath with β and μ
"""
bosonicbath(f::AbstractSpectrumFunction; β::Real=Inf, μ::Real=0) = (β == Inf) ? BosonicVacuum(f, μ=μ) : BosonicBath(f, β=β, μ=μ)


# function bosonicbath(f; β::Real=Inf, μ::Real=0) 
# 	y = f(0.)
# 	(y >= 0.) || throw(ArgumentError("spectrum function should return positive numbers."))
	
# 	(β == Inf) && return BosonicVacuum(f, convert(Float64, μ))
# 	(β >= 0.) || throw(ArgumentError("temperature should not be negative."))
# 	return BosonicBath(f, convert(Float64, β), convert(Float64, μ))
# end