
"""
	struct Bath{F <: AbstractBoundedFunction}

Fermionic bath container, includes a bath spectrum density,
the inverse temperature β and the chemical potential μ
"""
struct Bath{P<:AbstractParticle, F <: AbstractBoundedFunction} <: AbstractContinuousNormalBath{P}
	f::F
	β::Float64
	μ::Float64
end
Bath(::Type{P}, f::F; β::Real, μ::Real=0) where {P<:AbstractParticle, F<:AbstractBoundedFunction} = Bath{P, F}(f, float(β), float(μ))
Base.similar(x::Bath, ::Type{P}, f::AbstractBoundedFunction; β::Real=x.β, μ::Real=x.μ) where {P} = Bath(P, f, β=β, μ=μ)
Base.similar(x::Bath{P}, f::AbstractBoundedFunction; β::Real=x.β, μ::Real=x.μ) where {P} = Bath(P, f, β=β, μ=μ)
Base.similar(x::Bath{P}; f::AbstractBoundedFunction=x.f, β::Real=x.β, μ::Real=x.μ) where {P} = Bath(P, f, β=β, μ=μ)
Base.eltype(::Type{Bath{P, F}}) where {P, F} = Float64

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
struct Vacuum{P<:AbstractParticle, F <: AbstractBoundedFunction} <: AbstractContinuousNormalBath{P}
	f::F
	μ::Float64	
end
Vacuum(::Type{P}, f::F; μ::Real=0) where {P<:AbstractParticle, F<:AbstractBoundedFunction} = Vacuum{P, F}(f, convert(Float64, μ))
Base.similar(x::Vacuum, ::Type{P}, f::AbstractBoundedFunction; μ::Real=x.μ) where {P} = Vacuum(P, f, μ=μ)
Base.similar(x::Vacuum{P}, f::AbstractBoundedFunction; μ::Real=x.μ) where {P} = Vacuum(P, f, μ=μ)
Base.similar(x::Vacuum{P}; f::AbstractBoundedFunction=x.f, μ::Real=x.μ) where {P} = Vacuum(P, f, μ=μ)
Base.eltype(::Type{Vacuum{P, F}}) where {P, F} = Float64

const BosonicVacuum{F} = Vacuum{Boson, F} where {F<:AbstractBoundedFunction}
const FermionicVacuum{F} = Vacuum{Fermion, F} where {F<:AbstractBoundedFunction}

FermionicVacuum(f::AbstractBoundedFunction; kwargs...) = Vacuum(Fermion, f; kwargs...)
fermionicvacuum(f::AbstractBoundedFunction; kwargs...) = FermionicVacuum(f; kwargs...)
BosonicVacuum(f::AbstractBoundedFunction; kwargs...) = Vacuum(Boson, f; kwargs...)
bosonicvacuum(f::AbstractBoundedFunction; kwargs...) = BosonicVacuum(f; kwargs...)


const AbstractBosonicBath = Union{BosonicBath{F}, BosonicVacuum{F}} where {F<:AbstractBoundedFunction}
const AbstractFermionicBath = Union{FermionicBath{F}, FermionicVacuum{F}} where {F<:AbstractBoundedFunction}


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