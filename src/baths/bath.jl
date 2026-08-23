
"""
	struct Bath{P<:AbstractParticle, F<:AbstractBoundedFunction}

Bath container holding a bath spectral density `f`, the inverse temperature `β` and the
chemical potential `μ`. The type parameter `P` is the particle species (`Boson`/`Fermion`)
and `F` the spectral function type. The spectral density and temperature are accessible
via the properties `spectrum` and `T` (= 1/β).

	Bath(::Type{P}, f; β, μ=0)

Construct a bath of particle species `P`.
"""
struct Bath{P<:AbstractParticle, F <: AbstractBoundedFunction} <: AbstractNormalBath{P}
	f::F
	β::Float64
	μ::Float64
end
Bath(::Type{P}, f::F; β::Real, μ::Real=0) where {P<:AbstractParticle, F<:AbstractBoundedFunction} = Bath{P, F}(f, float(β), float(μ))
Base.similar(x::Bath, ::Type{P}, f::AbstractBoundedFunction; β::Real=x.β, μ::Real=x.μ) where {P} = Bath(P, f, β=β, μ=μ)
Base.similar(x::Bath{P}, f::AbstractBoundedFunction; β::Real=x.β, μ::Real=x.μ) where {P} = Bath(P, f, β=β, μ=μ)
Base.similar(x::Bath{P}; f::AbstractBoundedFunction=x.f, β::Real=x.β, μ::Real=x.μ) where {P} = Bath(P, f, β=β, μ=μ)
Base.eltype(::Type{Bath{P, F}}) where {P, F} = Float64

"""
	BosonicBath{F}

Type alias for a bosonic bath, equivalent to `Bath{Boson, F}`.
"""
const BosonicBath{F} = Bath{Boson, F} where {F<:AbstractBoundedFunction}
"""
	FermionicBath{F}

Type alias for a fermionic bath, equivalent to `Bath{Fermion, F}`.
"""
const FermionicBath{F} = Bath{Fermion, F} where {F<:AbstractBoundedFunction}

BosonicBath(f::AbstractBoundedFunction; kwargs...) = Bath(Boson, f; kwargs...)
"""
	bosonicbath(f; β, μ)

Construct a bosonic bath with spectral density `f`, inverse temperature `β` and
chemical potential `μ`.
"""
bosonicbath(f::AbstractBoundedFunction; kwargs...) = BosonicBath(f; kwargs...)
FermionicBath(f::AbstractBoundedFunction; kwargs...) = Bath(Fermion, f; kwargs...)
"""
	fermionicbath(f; β, μ)

Construct a fermionic bath with spectral density `f`, inverse temperature `β` and
chemical potential `μ`.
"""
fermionicbath(f::AbstractBoundedFunction; kwargs...) = FermionicBath(f; kwargs...)



"""
	struct Vacuum{P<:AbstractParticle, F<:AbstractBoundedFunction}

Vacuum (zero-temperature) bath container holding a bath spectral density `f` and the
chemical potential `μ`, with inverse temperature `β = Inf` (`T = 0`). The type parameter
`P` is the particle species (`Boson`/`Fermion`).

	Vacuum(::Type{P}, f; μ=0)

Construct a vacuum bath of particle species `P`.
"""
struct Vacuum{P<:AbstractParticle, F <: AbstractBoundedFunction} <: AbstractNormalBath{P}
	f::F
	μ::Float64	
end
Vacuum(::Type{P}, f::F; μ::Real=0) where {P<:AbstractParticle, F<:AbstractBoundedFunction} = Vacuum{P, F}(f, convert(Float64, μ))
Base.similar(x::Vacuum, ::Type{P}, f::AbstractBoundedFunction; μ::Real=x.μ) where {P} = Vacuum(P, f, μ=μ)
Base.similar(x::Vacuum{P}, f::AbstractBoundedFunction; μ::Real=x.μ) where {P} = Vacuum(P, f, μ=μ)
Base.similar(x::Vacuum{P}; f::AbstractBoundedFunction=x.f, μ::Real=x.μ) where {P} = Vacuum(P, f, μ=μ)
Base.eltype(::Type{Vacuum{P, F}}) where {P, F} = Float64

"""
	BosonicVacuum{F}

Type alias for a bosonic vacuum bath, equivalent to `Vacuum{Boson, F}`.
"""
const BosonicVacuum{F} = Vacuum{Boson, F} where {F<:AbstractBoundedFunction}
"""
	FermionicVacuum{F}

Type alias for a fermionic vacuum bath, equivalent to `Vacuum{Fermion, F}`.
"""
const FermionicVacuum{F} = Vacuum{Fermion, F} where {F<:AbstractBoundedFunction}

FermionicVacuum(f::AbstractBoundedFunction; kwargs...) = Vacuum(Fermion, f; kwargs...)
"""
	fermionicvacuum(f; μ)

Construct a fermionic vacuum bath (zero temperature) with spectral density `f` and
chemical potential `μ`.
"""
fermionicvacuum(f::AbstractBoundedFunction; kwargs...) = FermionicVacuum(f; kwargs...)
BosonicVacuum(f::AbstractBoundedFunction; kwargs...) = Vacuum(Boson, f; kwargs...)
"""
	bosonicvacuum(f; μ)

Construct a bosonic vacuum bath (zero temperature) with spectral density `f` and
chemical potential `μ`.
"""
bosonicvacuum(f::AbstractBoundedFunction; kwargs...) = BosonicVacuum(f; kwargs...)


# const AbstractBosonicBath = Union{BosonicBath{F}, BosonicVacuum{F}} where {F<:AbstractBoundedFunction}
# const AbstractFermionicBath = Union{FermionicBath{F}, FermionicVacuum{F}} where {F<:AbstractBoundedFunction}


"""
	bath(::Type{P}, f; β, μ=0)

Construct the corresponding particle bath for particle species `P` (`Boson`/`Fermion`).
"""
bath(::Type{Boson}, f::AbstractBoundedFunction; kwargs...) = bosonicbath(f; kwargs...)
bath(::Type{Fermion}, f::AbstractBoundedFunction; kwargs...) = fermionicbath(f; kwargs...)
"""
	vacuum(::Type{P}, f; μ=0)

Construct the corresponding vacuum (zero-temperature) bath for particle species `P`
(`Boson`/`Fermion`).
"""
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