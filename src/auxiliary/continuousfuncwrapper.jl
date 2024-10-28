# wrapper a discrete array into a continuous function

"""
	struct ContinuousFuncWrapper{T<:Number}
	omegas are discretized frequencies, fs are the functions values at the middle points.
"""
struct ContinuousFuncWrapper{T<:Number} <: Function
	ws::Vector{Float64}
	fs::Vector{T}

function ContinuousFuncWrapper{T}(ws::Vector{<:Real}, fs::AbstractVector{<:Number}) where {T <: Number}
	(length(ws) == length(fs)) || throw(ArgumentError("Number of frequencies must be equal to the number of function values"))
	(ws == sort(ws)) || throw(ArgumentError("frequencies must be sorted from small to large"))
	new{T}(convert(Vector{Float64}, ws), convert(Vector{T}, fs))
end

end

ContinuousFuncWrapper(ws::Vector{<:Real}, fs::AbstractVector{T}) where {T<:Number} = ContinuousFuncWrapper{T}(ws, fs)

function (x::ContinuousFuncWrapper)(ω::Real)
	((ω < freqs_min(x)) || (ω > freqs_max(x))) && return zero(eltype(x.fs))
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


freqs_min(x::ContinuousFuncWrapper) = x.ws[1]
freqs_max(x::ContinuousFuncWrapper) = x.ws[end]

SpectrumFunction(ws::Vector{<:Real}, fs::AbstractVector{<:Number}) = SpectrumFunction(ContinuousFuncWrapper(ws, fs), lb=ws[1], ub=ws[end])
