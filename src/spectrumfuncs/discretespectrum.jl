# wrapper a discrete array into a continuous function

"""
	struct DiscreteSpectrum{T<:Number}
	omegas are discretized frequencies, fs are the functions values at the middle points.
"""
struct DiscreteSpectrum{T<:Number} <: AbstractBoundedFunction
	ws::Vector{Float64}
	fs::Vector{T}

function DiscreteSpectrum{T}(ws::Vector{<:Real}, fs::AbstractVector{<:Number}) where {T <: Number}
	(length(ws) == length(fs)) || throw(ArgumentError("Number of frequencies must be equal to the number of function values"))
	(ws == sort(ws)) || throw(ArgumentError("frequencies must be sorted from small to large"))
	new{T}(convert(Vector{Float64}, ws), convert(Vector{T}, fs))
end

end

DiscreteSpectrum(ws::Vector{<:Real}, fs::AbstractVector{T}) where {T<:Number} = DiscreteSpectrum{T}(ws, fs)

function (x::DiscreteSpectrum)(ω::Real)
	((ω < lowerbound(x)) || (ω > upperbound(x))) && return zero(eltype(x.fs))
	(ω == x.ws[end]) && return x.fs[end]
	i = 0
	for pos in 1:(length(x.fs)-1)
		if (ω >= x.ws[pos]) && (ω < x.ws[pos+1])
			i = pos
			break;
		end
	end
	(i == 0) && error("something wrong")
	dif = (ω - x.ws[i]) / (x.ws[i+1] - x.ws[i])
	return x.fs[i] * (1 - dif) + x.fs[i+1] * dif
end


lowerbound(x::DiscreteSpectrum) = x.ws[1]
upperbound(x::DiscreteSpectrum) = x.ws[end]


quadgkwrapper(m::DiscreteSpectrum; kwargs...) = _quadgk(m, lowerbound(m), upperbound(m); kwargs...)
