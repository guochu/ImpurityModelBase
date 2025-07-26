

struct BCSBath{F <: AbstractBoundedFunction} <: AbstractBath{Fermion}
	f::F
	β::Float64
	μ::Float64
	Δ::Float64	
end

BCSBath(f::F; β::Real, Δ::Real=0, μ::Real=0) where {F<:AbstractBoundedFunction} = BCSBath{F}(f, convert(Float64, β), convert(Float64, μ), convert(Float64, Δ))
Base.similar(x::BCSBath, f::AbstractBoundedFunction; β::Real=x.β, μ::Real=x.μ, Δ::Real=x.Δ) = BCSBath(f, β=β, μ=μ, Δ=Δ)
Base.similar(x::BCSBath; f::AbstractBoundedFunction=x.f, β::Real=x.β, μ::Real=x.μ, Δ::Real=x.Δ) = BCSBath(f, β=β, μ=μ, Δ=Δ)


bcsbath(f::AbstractBoundedFunction; kwargs...) = BCSBath(f; kwargs...)
bcsbath(bath::FermionicBath; Δ::Real=0) = bcsbath(bath.spectrum, β=bath.β, μ=bath.μ, Δ=Δ)


struct BCSVacuum{F <: AbstractBoundedFunction} <: AbstractBath{Fermion}
	f::F
	μ::Float64
	Δ::Float64	
end
BCSVacuum(f::F; Δ::Real=0, μ::Real=0) where {F<:AbstractBoundedFunction} = BCSVacuum{F}(f, convert(Float64, μ), convert(Float64, Δ))
Base.similar(x::BCSVacuum, f::AbstractBoundedFunction; μ::Real=x.μ, Δ::Real=x.Δ) = BCSVacuum(f, μ=μ, Δ=Δ)
Base.similar(x::BCSVacuum; f::AbstractBoundedFunction=x.f, μ::Real=x.μ, Δ::Real=x.Δ) = BCSVacuum(f, μ=μ, Δ=Δ)

bcsvacuum(f::AbstractBoundedFunction; kwargs...) = BCSVacuum(f; kwargs...)
bcsvacuum(bath::FermionicVacuum; Δ::Real=0) = bcsvacuum(bath.spectrum, μ=bath.μ, Δ=Δ)

const AbstractBCSBath = Union{BCSBath{F}, BCSVacuum{F}} where {F<:AbstractBoundedFunction}


function Base.getproperty(m::BCSBath, s::Symbol)
	if s == :T
		return 1 / m.β
	elseif s == :spectrum
		return m.f
	else
		return getfield(m, s)
	end
end

function Base.getproperty(m::BCSVacuum, s::Symbol)
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