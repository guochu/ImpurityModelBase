"""
	struct GreenFunction{T<:Number, I}

A simple definition of real-frequency Green's Function
"""
struct GreenFunction{T<:Number, I} <: Function
	ws::Vector{Float64}
	fs::Vector{T}
	δ::Float64
	interp::I
end
function GreenFunction(ws::AbstractVector{<:Real}, fs::AbstractVector{T}; δ::Real=1.0e-8) where {T<:Number}
	ws = convert(Vector{Float64}, ws)
	fs = convert(Vector{T}, fs)
	interp = linear_interpolation(ws, fs)
	return GreenFunction(ws, fs, convert(Float64, δ), interp)
end

function (x::GreenFunction)(ω::Real)
	if (ω < first(x.ws)) || (ω > last(x.ws))
		return 1/(ω+im*x.δ)
	end
	return x.interp(ω)
end