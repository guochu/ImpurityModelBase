abstract type AbstractFermionicBath end

# fermionic baths
struct FermionicBath{F <: SpectrumFunction} <: AbstractFermionicBath
	f::F
	β::Float64
	μ::Float64
end
FermionicBath(f::SpectrumFunction; β::Real, μ::Real=0) = FermionicBath(f, convert(Float64, β), convert(Float64, μ))

function Base.getproperty(m::FermionicBath, s::Symbol)
	if s == :T
		return 1 / m.β
	elseif s == :spectrum
		return m.f
	else
		return getfield(m, s)
	end
end

struct FermionicVacuum{F <: SpectrumFunction} <: AbstractFermionicBath
	f::F
	μ::Float64	
end
FermionicVacuum(f::SpectrumFunction; μ::Real=0) = FermionicVacuum(f, convert(Float64, μ))

function Base.getproperty(m::FermionicVacuum, s::Symbol)
	if s == :β
		return Inf
	elseif s == :T
		return 0.
	elseif s == :spectrum
		return m.f
	else
		return getfield(m, s)
	end
end

function thermaloccupation(β::Real, μ::Real, ϵ::Real)
	n_k = 0.
	if β == Inf
		if ϵ > μ
			n_k = 0.
		else
			n_k = 1.
		end
	else
		n_k = 1 / (exp(β * (ϵ - μ)) + 1)
	end
	return n_k	
end
thermaloccupation(bath::AbstractFermionicBath, ϵ::Real) = thermaloccupation(bath.β, bath.μ, ϵ)

# function thermaloccupation(bath::AbstractFermionicBath, ϵ::Real)
# 	beta = bath.β
# 	mu = bath.μ
# 	n_k = 0.
# 	if beta == Inf
# 		if ϵ > mu
# 			n_k = 0.
# 		else
# 			n_k = 1.
# 		end
# 	else
# 		n_k = 1 / (exp(beta * (ϵ - mu)) + 1)
# 	end
# 	return n_k
# end


"""
	fermionicbath(f; β::Real=Inf, μ::Real=0.) 
	f should support y = f(x) with y positive 
"""
fermionicbath(f::SpectrumFunction; β::Real=Inf, μ::Real=0.) = (β == Inf) ? FermionicVacuum(f, μ=μ) : FermionicBath(f, β=β, μ=μ)






