"""
	AbstractBath{P}

Abstract type for a particle bath in a quantum impurity model; the type parameter `P`
is the particle species (`Boson`/`Fermion`).
"""
abstract type AbstractBath{P<:AbstractParticle} end
"""
	particletype(bath)
	particletype(::Type)

Return the particle species `P` (`Boson` or `Fermion`) of a bath.
"""
particletype(::Type{<:AbstractBath{P}}) where {P<:AbstractParticle} = P
particletype(x::AbstractBath) = particletype(typeof(x))
Base.eltype(x::AbstractBath) = eltype(typeof(x))

"""
	AbstractNormalBath{P}

Abstract type for a normal (unpaired) particle bath, characterized by a spectral density,
an inverse temperature `β` and a chemical potential `μ`.
"""
abstract type AbstractNormalBath{P} <: AbstractBath{P} end
"""
	AbstractBCSBath

Abstract type for a BCS (paired) fermionic bath, additionally characterized by a pairing
parameter `Δ`.
"""
abstract type AbstractBCSBath <: AbstractBath{Fermion} end
"""
	AbstractBECBath

Abstract type for a Bose-Einstein condensate (BEC) bosonic bath.
"""
abstract type AbstractBECBath <: AbstractBath{Boson} end

"""
	AbstractFermionicNormalBath

Alias for a normal fermionic bath, equivalent to `AbstractNormalBath{Fermion}`.
"""
const AbstractFermionicNormalBath = AbstractNormalBath{Fermion}
"""
	AbstractBosonicNormalBath

Alias for a normal bosonic bath, equivalent to `AbstractNormalBath{Boson}`.
"""
const AbstractBosonicNormalBath = AbstractNormalBath{Boson}

# abstract type AbstractContinuousBath{P} <: AbstractBath{P} end
# abstract type AbstractDiscreteBath{P} <: AbstractBath{P} end

# abstract type AbstractContinuousNormalBath{P} <: AbstractContinuousBath{P} end
# abstract type AbstractContinuousBCSBath <: AbstractContinuousBath{Fermion} end
# abstract type AbstractContinuousBECBath <: AbstractContinuousBath{Boson} end

# abstract type AbstractDiscreteNormalBath{P} <: AbstractDiscreteBath{P} end
# abstract type AbstractDiscreteBCSBath <: AbstractDiscreteBath{Fermion} end
# abstract type AbstractDiscreteBECBath <: AbstractDiscreteBath{Boson} end

# const AbstractNormalBath{P} = Union{AbstractContinuousNormalBath{P}, AbstractDiscreteNormalBath{P}} where {P<:AbstractParticle}
# const AbstractFermionicNormalBath = AbstractNormalBath{Fermion}
# const AbstractBosonicNormalBath = AbstractNormalBath{Boson}
# const AbstractBCSBath = Union{AbstractContinuousBCSBath, AbstractDiscreteBCSBath}
# const AbstractBECBath = Union{AbstractContinuousBECBath, AbstractDiscreteBECBath}


"""
	thermaloccupation(bath::AbstractNormalBath, energy::Real)

return n(ϵ)
"""
thermaloccupation(bath::AbstractNormalBath, ϵ::Real) = thermaloccupation(particletype(bath), bath.β, bath.μ, ϵ)
