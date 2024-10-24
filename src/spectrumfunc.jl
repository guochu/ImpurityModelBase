"""
	SpectrumFunction{F}

Wrapper for spectrum function, including a lower bound and upper bound
"""
struct SpectrumFunction{F}
	f::F
	lb::Float64
	ub::Float64
end
function SpectrumFunction(f; lb::Real=-Inf, ub::Real=Inf) 
	(lb < ub) || throw(ArgumentError("lb must be less than ub"))
	return SpectrumFunction(f, convert(Float64, lb), convert(Float64, ub))
end	
lowerbound(x::SpectrumFunction) = x.lb
upperbound(x::SpectrumFunction) = x.ub
# function (x::SpectrumFunction)(ϵ::Real) 
# 	if lowerbound(x) <= ϵ <= upperbound(x)
# 		return x.f(ϵ)
# 	else
# 		return 0.
# 	end
# end
function (x::SpectrumFunction)(ϵ::Real) 
	@assert lowerbound(x) <= ϵ <= upperbound(x)
	return x.f(ϵ)
end

"""
	semicircular(t::Real)

Often used for fermionic bath spectrum density
"""
function semicircular(t::Real)
	t = convert(Float64, t)
	D = 2*t
	return SpectrumFunction(ϵ -> sqrt(1-(ϵ/D)^2) * (D/π), lb = -D, ub = D)
end
semicircular(; t::Real=1) = semicircular(t)

"""
	Leggett(; α::Real, d::Real, ωc::Real)

J(ω) = (α/(2ωc))(ωᵈ/ωcᵈ)e^(-ω/ωc)	
"""
function Leggett(; α::Real=1, d::Real=1, ωc::Real=1)
	d = convert(Float64, d)
	return SpectrumFunction(ϵ -> (α/2)*(ϵ^d/ωc^(d-1))*exp(-ϵ/ωc), lb = 0, ub = ωc)
end


# wrapper a discrete array into a continuous function

"""
	struct ContinuousFuncWrapper{T<:Number}
	omegas are discretized frequencies, fs are the functions values at the middle points.
"""
struct ContinuousFuncWrapper{T<:Number}
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
	(i == 0) && error("something wrong.")
	dif = (ω - x.ws[i]) / (x.ws[i+1] - x.ws[i])
	return x.fs[i] * (1 - dif) + x.fs[i+1] * dif
end


freqs_min(x::ContinuousFuncWrapper) = x.ws[1]
freqs_max(x::ContinuousFuncWrapper) = x.ws[end]

SpectrumFunction(ws::Vector{<:Real}, fs::AbstractVector{<:Number}) = SpectrumFunction(ContinuousFuncWrapper(ws, fs), lb=ws[1], ub=ws[end])
