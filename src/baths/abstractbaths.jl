abstract type AbstractBath{P<:AbstractParticle} end
particletype(::Type{<:AbstractBath{P}}) where {P<:AbstractParticle} = P
particletype(x::AbstractBath) = particletype(typeof(x))
Base.eltype(x::AbstractBath) = eltype(typeof(x))


abstract type AbstractContinuousBath{P} <: AbstractBath{P} end
abstract type AbstractDiscreteBath{P} <: AbstractBath{P} end

abstract type AbstractContinuousNormalBath{P} <: AbstractContinuousBath{P} end
abstract type AbstractContinuousBCSBath <: AbstractContinuousBath{Fermion} end
abstract type AbstractContinuousBECBath <: AbstractContinuousBath{Boson} end

abstract type AbstractDiscreteNormalBath{P} <: AbstractDiscreteBath{P} end
abstract type AbstractDiscreteBCSBath <: AbstractDiscreteBath{Fermion} end
abstract type AbstractDiscreteBECBath <: AbstractDiscreteBath{Boson} end

const AbstractNormalBath{P} = Union{AbstractContinuousNormalBath{P}, AbstractDiscreteNormalBath{P}} where {P<:AbstractParticle}
const AbstractFermionicNormalBath = AbstractNormalBath{Fermion}
const AbstractBosonicNormalBath = AbstractNormalBath{Boson}
const AbstractBCSBath = Union{AbstractContinuousBCSBath, AbstractDiscreteBCSBath}
const AbstractBECBath = Union{AbstractContinuousBECBath, AbstractDiscreteBECBath}


"""
	thermaloccupation(bath::AbstractBosonicBath, energy::Real)
	thermaloccupation(bath::AbstractFermionicBath, energy::Real)

return n(ϵ)
"""
thermaloccupation(bath::Union{AbstractContinuousNormalBath, AbstractDiscreteNormalBath}, ϵ::Real) = thermaloccupation(particletype(bath), bath.β, bath.μ, ϵ)
