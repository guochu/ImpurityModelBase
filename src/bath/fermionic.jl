# fermionic baths
struct FermionicBath{F <: AbstractSpectrumFunction} <: AbstractFermionicBath
	f::F
	β::Float64
	μ::Float64
end
FermionicBath(f::AbstractSpectrumFunction; β::Real, μ::Real=0) = FermionicBath(f, convert(Float64, β), convert(Float64, μ))

struct FermionicVacuum{F <: AbstractSpectrumFunction} <: AbstractFermionicBath
	f::F
	μ::Float64	
end
FermionicVacuum(f::AbstractSpectrumFunction; μ::Real=0) = FermionicVacuum(f, convert(Float64, μ))


function fermidirac(β::Real, μ::Real, ϵ::Real)
	if β == Inf
		if ϵ > μ
			return 0.0
		elseif ϵ < μ
			return 1.0
		else
			return 0.5
		end
	else
		1.0/(1.0+exp(β * (ϵ-μ) ))
	end
end
thermaloccupation(bath::AbstractFermionicBath, ϵ::Real) = fermidirac(bath.β, bath.μ, ϵ)

"""
	fermionicbath(f; β::Real=Inf, μ::Real=0.) 
	f should support y = f(x) with y positive 
"""
fermionicbath(f::AbstractSpectrumFunction; β::Real=Inf, μ::Real=0) = (β == Inf) ? FermionicVacuum(f, μ=μ) : FermionicBath(f, β=β, μ=μ)

