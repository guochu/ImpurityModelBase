
# """
# 	struct DiscreteBath{P<:AbstractParticle}

# Fermionic bath container, includes a bath spectrum density,
# the inverse temperature β and the chemical potential μ
# """
# struct DiscreteBCSBath{T<:Number} <: AbstractDiscreteBCSBath
# 	ws::Vector{Float64}
# 	fs::Vector{Float64}
# 	β::Float64
# 	μ::Float64
# 	Δ::T
# end
# function DiscreteBCSBath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}, Δ::T; β::Real, μ::Real=0) where {T<:Number}
# 	(length(ws) == length(fs)) || throw(DimensionMismatch("num frequencies mismatch with num spectrum values"))
# 	all(x->x>=0, fs) || throw(ArgumentError("spectrum values can not be negative"))
# 	issorted(ws) || throw("frequencies should be sorted")
# 	DiscreteBCSBath{float(T)}(convert(Vector{Float64}, ws), convert(Vector{Float64}, fs), float(β), float(μ), float(Δ))
# end 
# DiscreteBCSBath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; β::Real, μ::Real=0, Δ::Number=0) = DiscreteBCSBath(ws, fs, Δ, β=β, μ=μ)
# Base.eltype(::Type{DiscreteBCSBath{T}}) where {T} = T


"""
	DiscreteBCSBath{T}

Type alias for a discrete BCS fermionic bath, equivalent to `BCSBath{DiscreteSpectrum, T}`.

	DiscreteBCSBath(ws, fs; β, μ=0, Δ=0)
"""
const DiscreteBCSBath{T<:Number} = BCSBath{DiscreteSpectrum, T}
DiscreteBCSBath(f::DiscreteSpectrum, Δ::Number; kwargs...) = BCSBath(f, Δ; kwargs...)
DiscreteBCSBath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}, Δ::Number; kwargs...) = DiscreteBCSBath(DiscreteSpectrum(ws, fs), Δ; kwargs...)
DiscreteBCSBath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; β::Real, μ::Real=0, Δ::Number=0) = DiscreteBCSBath(ws, fs, Δ, β=β, μ=μ)


"""
	discretebcsbath(ws, fs; β, μ=0, Δ=0)
	discretebcsbath(freqs, f; β, μ=0, Δ=0, atol=1e-6)

Construct a discrete BCS fermionic bath, either from frequencies `ws` and spectral values
`fs`, or by discretizing a continuous spectral function `f` on the intervals defined by
`freqs`.
"""
discretebcsbath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBCSBath(ws, fs; kwargs...)

# """
# 	struct DiscreteVacuum{P<:AbstractParticle}

# Fermionic bath container, includes a bath spectrum density,
# the chemical potential μ
# the inverse temperature β=Inf
# """
# struct DiscreteBCSVacuum{T<:Number} <: AbstractDiscreteBCSBath
# 	ws::Vector{Float64}
# 	fs::Vector{Float64}
# 	μ::Float64	
# 	Δ::T
# end
# function DiscreteBCSVacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}, Δ::T; μ::Real=0) where {T<:Number}
# 	(length(ws) == length(fs)) || throw(DimensionMismatch("num frequencies mismatch with num spectrum values"))
# 	all(x->x>=0, fs) || throw(ArgumentError("spectrum values can not be negative"))
# 	issorted(ws) || throw("frequencies should be sorted")
# 	DiscreteBCSVacuum{float(T)}(convert(Vector{Float64}, ws), convert(Vector{Float64}, fs), float(μ), float(Δ))
# end 
# DiscreteBCSVacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; μ::Real=0, Δ::Number=0) = DiscreteBCSVacuum(ws, fs, Δ, μ=μ)
# discretebcsvacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBCSVacuum(ws, fs; kwargs...)
# Base.eltype(::Type{DiscreteBCSVacuum{T}}) where {T} = T

"""
	DiscreteBCSVacuum{T}

Type alias for a discrete BCS fermionic vacuum bath, equivalent to
`BCSVacuum{DiscreteSpectrum, T}`.

	DiscreteBCSVacuum(ws, fs; μ=0, Δ=0)
"""
const DiscreteBCSVacuum{T<:Number} = BCSVacuum{DiscreteSpectrum, T}

DiscreteBCSVacuum(f::DiscreteSpectrum, Δ::Number; kwargs...) = BCSVacuum(f, Δ; kwargs...)
DiscreteBCSVacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}, Δ::Number; kwargs...) = DiscreteBCSVacuum(DiscreteSpectrum(ws, fs), Δ; kwargs...)
DiscreteBCSVacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; μ::Real=0, Δ::Number=0) = DiscreteBCSVacuum(ws, fs, Δ, μ=μ)
"""
	discretebcsvacuum(ws, fs; μ=0, Δ=0)
	discretebcsvacuum(freqs, f; μ=0, Δ=0, atol=1e-6)

Construct a discrete BCS fermionic vacuum (zero-temperature) bath; the parameters have the
same meaning as in [`discretebcsbath`](@ref).
"""
discretebcsvacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBCSVacuum(ws, fs; kwargs...)


const AbstractDiscreteBCSBath{T<:Number} = Union{DiscreteBCSBath{T}, DiscreteBCSVacuum{T}}

"""
	frequencies(b::AbstractDiscreteBCSBath)

Return the list of discrete frequencies of the discrete BCS bath `b`.
"""
frequencies(b::AbstractDiscreteBCSBath) = frequencies(b.f)
"""
	spectrumvalues(b::AbstractDiscreteBCSBath)

Return the list of spectral values of the discrete BCS bath `b` at each frequency.
"""
spectrumvalues(b::AbstractDiscreteBCSBath) = spectrumvalues(b.f)
"""
	spectrumcouplings(b::AbstractDiscreteBCSBath)

Return the coupling strengths of the discrete BCS bath `b`, given by the square roots of
the spectral values.
"""
spectrumcouplings(b::AbstractDiscreteBCSBath) = spectrumcouplings(b.f)
"""
	num_sites(b::AbstractDiscreteBCSBath)

Return the number of modes corresponding to the discrete BCS bath `b`, twice the number of
frequencies (particle and hole branches).
"""
num_sites(x::AbstractDiscreteBCSBath) = 2 * length(frequencies(x))


# num_sites(b::Union{DiscreteBCSBath, DiscreteBCSVacuum}) = 2*length(frequencies(b))

# function Base.getproperty(m::DiscreteBCSBath, s::Symbol)
# 	if s == :T
# 		return 1 / m.β
# 	else
# 		return getfield(m, s)
# 	end
# end

# function Base.getproperty(m::DiscreteBCSVacuum, s::Symbol)
# 	if s == :β
# 		return Inf
# 	elseif s == :T
# 		return 0.
# 	else
# 		return getfield(m, s)
# 	end
# end


# const AbstractDiscreteBCSBath = Union{DiscreteBCSBath, DiscreteBCSVacuum}


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
	return discretebcsvacuum(freqs, f; μ=b.μ, Δ=b.Δ, kwargs...)
end