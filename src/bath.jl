# definition of particle type
abstract type AbstractParticle end
struct Boson <: AbstractParticle end
struct Fermion <: AbstractParticle end

abstract type AbstractBath{P<:AbstractParticle} end
particletype(::Type{<:AbstractBath{P}}) where {P<:AbstractParticle} = P
particletype(x::AbstractBath) = particletype(typeof(x))

"""
	struct Bath{F <: AbstractSpectrumFunction}

Fermionic bath container, includes a bath spectrum density,
the inverse temperature β and the chemical potential μ
"""
struct Bath{P<:AbstractParticle, F <: AbstractSpectrumFunction} <: AbstractBath{P}
	f::F
	β::Float64
	μ::Float64
end
Bath(::Type{P}, f::F; β::Real, μ::Real=0) where {P<:AbstractParticle, F<:AbstractSpectrumFunction} = Bath{P, F}(f, convert(Float64, β), convert(Float64, μ))
Base.similar(x::Bath, ::Type{P}, f::AbstractSpectrumFunction; β::Real=x.β, μ=x.μ) where {P} = Bath(P, f, β=β, μ=μ)
Base.similar(x::Bath{P}, f::AbstractSpectrumFunction; β::Real=x.β, μ=x.μ) where {P} = Bath(P, f, β=β, μ=μ)
Base.similar(x::Bath{P}; f::AbstractSpectrumFunction=x.f, β::Real=x.β, μ=x.μ) where {P} = Bath(P, f, β=β, μ=μ)


const BosonicBath{F} = Bath{Boson, F} where {F<:AbstractSpectrumFunction}
const FermionicBath{F} = Bath{Fermion, F} where {F<:AbstractSpectrumFunction}

BosonicBath(f::AbstractSpectrumFunction; kwargs...) = Bath(Boson, f; kwargs...)
"""
	bosonicbath(f; β, μ) 

Return a bosonic bath with β and μ
"""
bosonicbath(f::AbstractSpectrumFunction; kwargs...) = BosonicBath(f; kwargs...)
FermionicBath(f::AbstractSpectrumFunction; kwargs...) = Bath(Fermion, f; kwargs...)
"""
	fermionicbath(f; β, μ) 

Return a fermionic bath with β and μ
"""
fermionicbath(f::AbstractSpectrumFunction; kwargs...) = FermionicBath(f; kwargs...)



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
Base.similar(x::Vacuum, ::Type{P}, f::AbstractSpectrumFunction; μ=x.μ) where {P} = Vacuum(P, f, μ=μ)
Base.similar(x::Vacuum{P}, f::AbstractSpectrumFunction; β::Real=x.β, μ=x.μ) where {P} = Vacuum(P, f, μ=μ)
Base.similar(x::Vacuum{P}; f::AbstractSpectrumFunction=x.f, β::Real=x.β, μ=x.μ) where {P} = Vacuum(P, f, μ=μ)

const BosonicVacuum{F} = Vacuum{Boson, F} where {F<:AbstractSpectrumFunction}
const FermionicVacuum{F} = Vacuum{Fermion, F} where {F<:AbstractSpectrumFunction}

FermionicVacuum(f::AbstractSpectrumFunction; kwargs...) = Vacuum(Fermion, f; kwargs...)
fermionicvacuum(f::AbstractSpectrumFunction; kwargs...) = FermionicVacuum(f; kwargs...)
BosonicVacuum(f::AbstractSpectrumFunction; kwargs...) = Vacuum(Boson, f; kwargs...)
bosonicvacuum(f::AbstractSpectrumFunction; kwargs...) = BosonicVacuum(f; kwargs...)


const AbstractBosonicBath = Union{BosonicBath{F}, BosonicVacuum{F}} where {F<:AbstractSpectrumFunction}
const AbstractFermionicBath = Union{FermionicBath{F}, FermionicVacuum{F}} where {F<:AbstractSpectrumFunction}

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