# definition of particle type
abstract type AbstractParticle end
struct Boson <: AbstractParticle end
struct Fermion <: AbstractParticle end


"""
	boseeinstein(β, μ, ϵ)

Boson-einstein distribution for a bosonic bath 
with β, μ at energy ϵ
μ = 0 by default
"""
function boseeinstein(β::Real, ϵ::Real)
	(ϵ > 0) || throw(ArgumentError("energy must be larger than μ"))
	x = exp(-safe_mult(β, ϵ))
	return x / (1 - x)
end

boseeinstein(β::Real, μ::Real, ϵ::Real) = boseeinstein(β, ϵ - μ)

# function boseeinstein(β::Real, μ::Real, ϵ::Real)
# 	(ϵ > μ) || throw(ArgumentError("energy must be larger than μ"))
# 	return 1 / (exp(safe_mult(β, ϵ - μ)) - 1)
# end
thermaloccupation(::Type{Boson}, β::Real, μ::Real, ϵ::Real) = boseeinstein(β, μ, ϵ)
thermaloccupation(::Type{Boson}, β::Real, ϵ::Real) = boseeinstein(β, ϵ)

"""
	fermidirac(β, μ, ϵ)

Return Fermi-Dirac distribution for a fermionic bath 
with β, μ at energy ϵ
"""
function fermidirac(β::Real, ϵ::Real)
	if ϵ >= 0
		x = exp(-safe_mult(β, ϵ))
		return x/(1+x)
	else
		return 1.0/(1.0+exp(safe_mult(β, ϵ)))
	end
end
fermidirac(β::Real, μ::Real, ϵ::Real) = fermidirac(β, ϵ-μ)


# function fermidirac(β::Real, μ::Real, ϵ::Real)
# 	x = exp(-safe_mult(β, ϵ-μ))
# 	return x/(1+x)
# end
thermaloccupation(::Type{Fermion}, β::Real, μ::Real, ϵ::Real) = fermidirac(β, μ, ϵ)
thermaloccupation(::Type{Fermion}, β::Real, ϵ::Real) = fermidirac(β, ϵ)


function safe_mult(β::Real, ϵ::Real)
	if (β == Inf) && (ϵ == 0)
		return zero(ϵ)
	end
	return β * ϵ
end