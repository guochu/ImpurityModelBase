
"""
	struct DiscreteSpectrum{T<:Number}

Wrapper for discrete spectrum function,
linear interpolation is used for interior points,
zero value is assumed for frequency outside the range.
ws are discretized frequencies, fs are the functions 
values at the middle points
"""
struct DiscreteSpectrum{T<:Number, I} <: AbstractBoundedFunction
	ws::Vector{Float64}
	fs::Vector{T}
	interp::I
end

function DiscreteSpectrum(ws::Vector{<:Real}, fs::AbstractVector{T}) where {T<:Number}
	all(fs .>= 0) || throw(ArgumentError("spectrum values should be positive"))
	ws = convert(Vector{Float64}, ws)
	fs = convert(Vector{T}, fs)
	interp = linear_interpolation(ws, fs)
	return DiscreteSpectrum(ws, fs, interp)
end 

lowerbound(x::DiscreteSpectrum) = x.ws[1]
upperbound(x::DiscreteSpectrum) = x.ws[end]

(x::DiscreteSpectrum)(ω::Real) = x.interp(ω)
(x::DiscreteSpectrum)(ω::Union{Vector, AbstractRange}) = x.interp.(ω)

# quadgkwrapper(m::DiscreteSpectrum; kwargs...) = _quadgk(m, lowerbound(m), upperbound(m); kwargs...)
function quadgkwrapper(m::DiscreteSpectrum)
	ws, fs = xvalues(m.data), fvalues(m.data)
	r = zero(eltype(fs))
	for i in 1:length(ws)-1
		dw = ws[i+1] - ws[i]
		f = (fs[i+1] + fs[i]) / 2
		r += f * dw
	end
	return r
end
spectrumshift(m::DiscreteSpectrum, μ::Real) = DiscreteSpectrum(xvalues(m.data) .+ μ, fvalues(m.data))
