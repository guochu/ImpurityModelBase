abstract type AbstractDiscreteFunction <: AbstractBoundedFunction end

"""
	struct DiscreteSpectrum{T<:Number}

Wrapper for discrete spectrum function,
linear interpolation is used for interior points,
zero value is assumed for frequency outside the range.
ws are discretized frequencies, fs are the functions 
values at the middle points
"""
struct DiscreteSpectrum <: AbstractDiscreteFunction
	ws::Vector{Float64}
	fs::Vector{Float64}
	# The validation is defined as an inner constructor so that it is not
	# shadowed by the default inner constructor for `Vector{Float64}` inputs
	function DiscreteSpectrum(ws::Vector{<:Real}, fs::AbstractVector{<:Real})
		(length(ws) == length(fs)) || throw(DimensionMismatch("n frequencies mismatch with n spectrum values"))
		all(x->x>=0, fs) || throw(ArgumentError("spectrum values can not be negative"))
		issorted(ws) || throw(ArgumentError("frequencies should be sorted"))
		return new(convert(Vector{Float64}, ws), convert(Vector{Float64}, fs))
	end
end

# lowerbound(x::DiscreteSpectrum) = x.ws[1]
# upperbound(x::DiscreteSpectrum) = x.ws[end]
Base.similar(x::DiscreteSpectrum, ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}) = DiscreteSpectrum(ws, fs)



lowerbound(x::DiscreteSpectrum) = -Inf
upperbound(x::DiscreteSpectrum) = Inf
"""
	frequencies(x::DiscreteSpectrum)

Return the list of discrete frequencies `ws` of the discrete spectrum `x`.
"""
frequencies(x::DiscreteSpectrum) = x.ws
"""
	spectrumvalues(x::DiscreteSpectrum)

Return the list of spectral values `fs` of the discrete spectrum `x` at each frequency.
"""
spectrumvalues(x::DiscreteSpectrum) = x.fs
"""
	spectrumcouplings(x::DiscreteSpectrum)

Return the coupling strengths of the discrete spectrum `x`, given by the square roots
of the spectral values `√fs`.
"""
spectrumcouplings(x::DiscreteSpectrum) = sqrt.(spectrumvalues(x))

# (x::DiscreteSpectrum)(ω::Real) = x.interp(ω)
# (x::DiscreteSpectrum)(ω::Union{Vector, AbstractRange}) = x.interp.(ω)

(x::DiscreteSpectrum)(ϵ::Real) = error("can not call to a discrete spectrum function")

# quadgkwrapper(m::DiscreteSpectrum; kwargs...) = _quadgk(m, lowerbound(m), upperbound(m); kwargs...)
# function quadgkwrapper(m::DiscreteSpectrum)
# 	ws, fs = frequencies(m.data), spectrumvalues(m.data)
# 	r = zero(eltype(fs))
# 	for i in 1:length(ws)-1
# 		dw = ws[i+1] - ws[i]
# 		f = (fs[i+1] + fs[i]) / 2
# 		r += f * dw
# 	end
# 	return r
# end

spectrumshift(m::DiscreteSpectrum, μ::Real) = DiscreteSpectrum(frequencies(m) .+ μ, spectrumvalues(m))