"""
	struct FermionicBath{F <: AbstractSpectrumFunction}

Fermionic bath container, includes a bath spectrum density,
the inverse temperature β and the chemical potential μ
"""
struct FermionicBath{F <: AbstractSpectrumFunction} <: AbstractFermionicBath
	f::F
	β::Float64
	μ::Float64
end
FermionicBath(f::AbstractSpectrumFunction; β::Real, μ::Real=0) = FermionicBath(f, convert(Float64, β), convert(Float64, μ))

"""
	struct FermionicVacuum{F <: AbstractSpectrumFunction}

Fermionic bath container, includes a bath spectrum density,
the chemical potential μ
the inverse temperature β=Inf
"""
struct FermionicVacuum{F <: AbstractSpectrumFunction} <: AbstractFermionicBath
	f::F
	μ::Float64	
end
FermionicVacuum(f::AbstractSpectrumFunction; μ::Real=0) = FermionicVacuum(f, convert(Float64, μ))

"""
	fermidirac(β, μ, ϵ)

Return Fermi-Dirac distribution for a fermionic bath 
with β, μ at energy ϵ
"""
function fermidirac(β::Real, μ::Real, ϵ::Real)
	return 1.0/(1.0+exp(safe_mult(β, ϵ-μ)))
	# if β == Inf
	# 	if ϵ > μ
	# 		return 0.0
	# 	elseif ϵ < μ
	# 		return 1.0
	# 	else
	# 		return 0.5
	# 	end
	# else
	# 	1.0/(1.0+exp(β * (ϵ-μ) ))
	# end
end
thermaloccupation(bath::AbstractFermionicBath, ϵ::Real) = fermidirac(bath.β, bath.μ, ϵ)

"""
	fermionicbath(f; β, μ) 

Return a fermionic bath with β and μ
"""
fermionicbath(f::AbstractSpectrumFunction; β::Real=Inf, μ::Real=0) = (β == Inf) ? FermionicVacuum(f, μ=μ) : FermionicBath(f, β=β, μ=μ)

