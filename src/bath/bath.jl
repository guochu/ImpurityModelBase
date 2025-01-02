# definition of particle type
abstract type AbstractParticle end
struct Boson <: AbstractParticle end
struct Fermion <: AbstractParticle end

abstract type AbstractBath{P<:AbstractParticle} end
particletype(::Type{<:AbstractBath{P}}) where {P<:AbstractParticle} = P
particletype(x::AbstractBath) = particletype(typeof(x))

"""
	struct FermionicBath{F <: AbstractSpectrumFunction}

Fermionic bath container, includes a bath spectrum density,
the inverse temperature β and the chemical potential μ
"""
struct Bath{P<:AbstractParticle, F <: AbstractSpectrumFunction} <: AbstractBath{P}
	f::F
	β::Float64
	μ::Float64
end
Bath(::Type{P}, f::F; β::Real, μ::Real=0) where {P<:AbstractParticle, F<:AbstractSpectrumFunction} = Bath{P, F}(f, convert(Float64, β), convert(Float64, μ))


const BosonicBath{F} = Bath{Boson, F} where {F<:AbstractSpectrumFunction}
const FermionicBath{F} = Bath{Fermion, F} where {F<:AbstractSpectrumFunction}

BosonicBath(f::AbstractSpectrumFunction; kwargs...) = Bath(Boson, f; kwargs...)
FermionicBath(f::AbstractSpectrumFunction; kwargs...) = Bath(Fermion, f; kwargs...)

"""
	struct FermionicVacuum{F <: AbstractSpectrumFunction}

Fermionic bath container, includes a bath spectrum density,
the chemical potential μ
the inverse temperature β=Inf
"""
struct Vacuum{P<:AbstractParticle, F <: AbstractSpectrumFunction} <: AbstractBath{P}
	f::F
	μ::Float64	
end
Vacuum(::Type{P}, f::F; μ::Real=0) where {P<:AbstractParticle, F<:AbstractSpectrumFunction} = Vacuum{P, F}(f, convert(Float64, μ))

const BosonicVacuum{F} = Vacuum{Boson, F} where {F<:AbstractSpectrumFunction}
const FermionicVacuum{F} = Vacuum{Fermion, F} where {F<:AbstractSpectrumFunction}

FermionicVacuum(f::AbstractSpectrumFunction; kwargs...) = Vacuum(Fermion, f; kwargs...)
BosonicVacuum(f::AbstractSpectrumFunction; kwargs...) = Vacuum(Boson, f; kwargs...)


const AbstractBosonicBath = Union{BosonicBath{F}, BosonicVacuum{F}} where {F<:AbstractSpectrumFunction}
const AbstractFermionicBath = Union{FermionicBath{F}, FermionicVacuum{F}} where {F<:AbstractSpectrumFunction}

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
thermaloccupation(::Type{Boson}, β::Real, μ::Real, ϵ::Real) = boseeinstein(β, μ, ϵ)

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
thermaloccupation(::Type{Fermion}, β::Real, μ::Real, ϵ::Real) = fermidirac(β, μ, ϵ)

"""
	thermaloccupation(bath::BosonicBath, energy::Real)

return n(ϵ)
"""
thermaloccupation(bath::AbstractBosonicBath, ϵ::Real) = boseeinstein(bath.β, bath.μ, ϵ)
thermaloccupation(bath::AbstractFermionicBath, ϵ::Real) = fermidirac(bath.β, bath.μ, ϵ)


"""
	bosonicbath(f; β, μ) 

Return a bosonic bath with β and μ
"""
bosonicbath(f::AbstractSpectrumFunction; β::Real=Inf, μ::Real=0) = (β == Inf) ? BosonicVacuum(f, μ=μ) : BosonicBath(f, β=β, μ=μ)


"""
	fermionicbath(f; β, μ) 

Return a fermionic bath with β and μ
"""
fermionicbath(f::AbstractSpectrumFunction; β::Real=Inf, μ::Real=0) = (β == Inf) ? FermionicVacuum(f, μ=μ) : FermionicBath(f, β=β, μ=μ)


bath(::Type{Boson}, f::AbstractSpectrumFunction; kwargs...) = bosonicbath(f; kwargs...)
bath(::Type{Fermion}, f::AbstractSpectrumFunction; kwargs...) = fermionicbath(f; kwargs...)
vacuum(::Type{Boson}, f::AbstractSpectrumFunction; kwargs...) = BosonicVacuum(f; kwargs...)
vacuum(::Type{Fermion}, f::AbstractSpectrumFunction; kwargs...) = FermionicVacuum(f; kwargs...)

function Base.getproperty(m::Bath, s::Symbol)
	if s == :T
		return 1 / m.β
	elseif s == :spectrum
		return m.f
	else
		return getfield(m, s)
	end
end

function Base.getproperty(m::Vacuum, s::Symbol)
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

function safe_mult(β::Real, ϵ::Real)
	if (β == Inf) && (ϵ == 0)
		return zero(ϵ)
	end
	return β * ϵ
end