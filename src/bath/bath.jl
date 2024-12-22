abstract type AbstractBath end
abstract type AbstractFermionicBath <: AbstractBath end
abstract type AbstractBosonicBath <: AbstractBath end


include("fermionic.jl")
include("bosonic.jl")

const ThermalBath{F} = Union{FermionicBath{F}, BosonicBath{F}} where F
const Vacuum{F} = Union{FermionicVacuum{F}, BosonicVacuum{F}} where F

function Base.getproperty(m::ThermalBath, s::Symbol)
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