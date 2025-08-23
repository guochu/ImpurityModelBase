"""
	struct DiscreteBath{P<:AbstractParticle}

Fermionic bath container, includes a bath spectrum density,
the inverse temperature β and the chemical potential μ
"""
struct DiscreteBCSBath{T<:Number} <: AbstractDiscreteBath{Fermion}
	ws::Vector{Float64}
	fs::Vector{Float64}
	β::Float64
	μ::Float64
	Δ::T
end
function DiscreteBCSBath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}, Δ::T; β::Real, μ::Real=0) where {T<:Number}
	(length(ws) == length(fs)) || throw(DimensionMismatch("num frequencies mismatch with num spectrum values"))
	all(x->x>=0, fs) || throw(ArgumentError("spectrum values can not be negative"))
	issorted(ws) || throw("frequencies should be sorted")
	DiscreteBCSBath{float(T)}(convert(Vector{Float64}, ws), convert(Vector{Float64}, fs), float(β), float(μ), float(Δ))
end 
DiscreteBCSBath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; β::Real, μ::Real=0, Δ::Number=0) = DiscreteBCSBath(ws, fs, Δ, β=β, μ=μ)
Base.eltype(::Type{DiscreteBCSBath{T}}) where {T} = T

"""
	DiscreteBCSBath(ws, fs; β, μ, Δ) 

Return a fermionic bath with β and μ
"""
discretebcsbath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBCSBath(ws, fs; kwargs...)

"""
	struct DiscreteVacuum{P<:AbstractParticle}

Fermionic bath container, includes a bath spectrum density,
the chemical potential μ
the inverse temperature β=Inf
"""
struct DiscreteBCSVacuum{T<:Number} <: AbstractDiscreteBath{Fermion}
	ws::Vector{Float64}
	fs::Vector{Float64}
	μ::Float64	
	Δ::T
end
function DiscreteBCSVacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}, Δ::T; μ::Real=0) where {T<:Number}
	(length(ws) == length(fs)) || throw(DimensionMismatch("num frequencies mismatch with num spectrum values"))
	all(x->x>=0, fs) || throw(ArgumentError("spectrum values can not be negative"))
	issorted(ws) || throw("frequencies should be sorted")
	DiscreteBCSVacuum{float(T)}(convert(Vector{Float64}, ws), convert(Vector{Float64}, fs), float(μ), float(Δ))
end 
DiscreteBCSVacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; μ::Real=0, Δ::Number=0) = DiscreteBCSVacuum(ws, fs, Δ, μ=μ)
discretebcsvacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBCSVacuum(ws, fs; kwargs...)
Base.eltype(::Type{DiscreteBCSVacuum{T}}) where {T} = T

num_sites(b::Union{DiscreteBCSBath, DiscreteBCSVacuum}) = 2*length(frequencies(b))

function Base.getproperty(m::DiscreteBCSBath, s::Symbol)
	if s == :T
		return 1 / m.β
	else
		return getfield(m, s)
	end
end

function Base.getproperty(m::DiscreteBCSVacuum, s::Symbol)
	if s == :β
		return Inf
	elseif s == :T
		return 0.
	else
		return getfield(m, s)
	end
end


const AbstractDiscreteBCSBath = Union{DiscreteBCSBath, DiscreteBCSVacuum}


function discretebcsbath(freqs::Union{Vector{<:Real}, AbstractRange}, f::Function; atol::Real=1.0e-6, β::Real, μ::Real=0, Δ::Number=0, kwargs...)
	omegas, couplings = spectrum_couplings(freqs, f; atol=atol, kwargs...)
	return discretebcsbath(omegas, couplings; β=β, μ=μ, Δ=Δ)
end
function discretebath(b::BCSBath; δw::Real=0.1, kwargs...)
	f = b.spectrum
	freqs = lowerbound(f):δw:upperbound(f)
	return discretebcsbath(freqs, f; β=b.β, μ=b.μ, Δ=b.Δ, kwargs...)
end
function discretebcsvacuum(freqs::Union{Vector{<:Real}, AbstractRange}, f::Function; atol::Real=1.0e-6, μ::Real=0, Δ::Number=0, kwargs...) 
	omegas, couplings = spectrum_couplings(freqs, f; atol=atol, kwargs...)
	return discretebcsvacuum(omegas, couplings; μ=μ, Δ=Δ) 
end
function discretevacuum(b::BCSVacuum; δw::Real=0.1, kwargs...)
	f = b.spectrum
	freqs = lowerbound(f):δw:upperbound(f)
	return discretebcsvacuum(particletype(b), freqs, f; μ=b.μ, Δ=b.Δ, kwargs...)
end