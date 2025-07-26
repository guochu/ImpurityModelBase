# definition of particle type
abstract type AbstractParticle end
struct Boson <: AbstractParticle end
struct Fermion <: AbstractParticle end

abstract type AbstractBath{P<:AbstractParticle} end
particletype(::Type{<:AbstractBath{P}}) where {P<:AbstractParticle} = P
particletype(x::AbstractBath) = particletype(typeof(x))

"""
	struct Bath{F <: AbstractBoundedFunction}

Fermionic bath container, includes a bath spectrum density,
the inverse temperature β and the chemical potential μ
"""
struct Bath{P<:AbstractParticle, F <: AbstractBoundedFunction} <: AbstractBath{P}
	f::F
	β::Float64
	μ::Float64
end
Bath(::Type{P}, f::F; β::Real, μ::Real=0) where {P<:AbstractParticle, F<:AbstractBoundedFunction} = Bath{P, F}(f, convert(Float64, β), convert(Float64, μ))
Base.similar(x::Bath, ::Type{P}, f::AbstractBoundedFunction; β::Real=x.β, μ::Real=x.μ) where {P} = Bath(P, f, β=β, μ=μ)
Base.similar(x::Bath{P}, f::AbstractBoundedFunction; β::Real=x.β, μ::Real=x.μ) where {P} = Bath(P, f, β=β, μ=μ)
Base.similar(x::Bath{P}; f::AbstractBoundedFunction=x.f, β::Real=x.β, μ::Real=x.μ) where {P} = Bath(P, f, β=β, μ=μ)


const BosonicBath{F} = Bath{Boson, F} where {F<:AbstractBoundedFunction}
const FermionicBath{F} = Bath{Fermion, F} where {F<:AbstractBoundedFunction}

BosonicBath(f::AbstractBoundedFunction; kwargs...) = Bath(Boson, f; kwargs...)
"""
	bosonicbath(f; β, μ) 

Return a bosonic bath with β and μ
"""
bosonicbath(f::AbstractBoundedFunction; kwargs...) = BosonicBath(f; kwargs...)
FermionicBath(f::AbstractBoundedFunction; kwargs...) = Bath(Fermion, f; kwargs...)
"""
	fermionicbath(f; β, μ) 

Return a fermionic bath with β and μ
"""
fermionicbath(f::AbstractBoundedFunction; kwargs...) = FermionicBath(f; kwargs...)



"""
	struct FermionicVacuum{F <: AbstractBoundedFunction}

Fermionic bath container, includes a bath spectrum density,
the chemical potential μ
the inverse temperature β=Inf
"""
struct Vacuum{P<:AbstractParticle, F <: AbstractBoundedFunction} <: AbstractBath{P}
	f::F
	μ::Float64	
end
Vacuum(::Type{P}, f::F; μ::Real=0) where {P<:AbstractParticle, F<:AbstractBoundedFunction} = Vacuum{P, F}(f, convert(Float64, μ))
Base.similar(x::Vacuum, ::Type{P}, f::AbstractBoundedFunction; μ::Real=x.μ) where {P} = Vacuum(P, f, μ=μ)
Base.similar(x::Vacuum{P}, f::AbstractBoundedFunction; μ::Real=x.μ) where {P} = Vacuum(P, f, μ=μ)
Base.similar(x::Vacuum{P}; f::AbstractBoundedFunction=x.f, μ::Real=x.μ) where {P} = Vacuum(P, f, μ=μ)

const BosonicVacuum{F} = Vacuum{Boson, F} where {F<:AbstractBoundedFunction}
const FermionicVacuum{F} = Vacuum{Fermion, F} where {F<:AbstractBoundedFunction}

FermionicVacuum(f::AbstractBoundedFunction; kwargs...) = Vacuum(Fermion, f; kwargs...)
fermionicvacuum(f::AbstractBoundedFunction; kwargs...) = FermionicVacuum(f; kwargs...)
BosonicVacuum(f::AbstractBoundedFunction; kwargs...) = Vacuum(Boson, f; kwargs...)
bosonicvacuum(f::AbstractBoundedFunction; kwargs...) = BosonicVacuum(f; kwargs...)


const AbstractBosonicBath = Union{BosonicBath{F}, BosonicVacuum{F}} where {F<:AbstractBoundedFunction}
const AbstractFermionicBath = Union{FermionicBath{F}, FermionicVacuum{F}} where {F<:AbstractBoundedFunction}

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

"""
	thermaloccupation(bath::AbstractBosonicBath, energy::Real)
	thermaloccupation(bath::AbstractFermionicBath, energy::Real)

return n(ϵ)
"""
thermaloccupation(bath::AbstractBath, ϵ::Real) = thermaloccupation(particletype(bath), bath.β, bath.μ, ϵ)

bath(::Type{Boson}, f::AbstractBoundedFunction; kwargs...) = bosonicbath(f; kwargs...)
bath(::Type{Fermion}, f::AbstractBoundedFunction; kwargs...) = fermionicbath(f; kwargs...)
vacuum(::Type{Boson}, f::AbstractBoundedFunction; kwargs...) = BosonicVacuum(f; kwargs...)
vacuum(::Type{Fermion}, f::AbstractBoundedFunction; kwargs...) = FermionicVacuum(f; kwargs...)

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