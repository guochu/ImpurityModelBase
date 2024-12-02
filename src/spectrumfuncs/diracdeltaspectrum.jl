"""
	struct DiracDelta

Dirac delta function
"""
struct DiracDelta <: AbstractSpectrumFunction
	ω::Float64
	α::Float64
	lb::Float64
	ub::Float64
end
function DiracDelta(;ω::Real=0, α::Real=1, lb::Real=-Inf, ub::Real=Inf)
	(lb < ub) || throw(ArgumentError("lb must be less than ub"))
	return DiracDelta(convert(Float64, ω), convert(Float64, α), convert(Float64, lb), convert(Float64, ub))
end
(x::DiracDelta)(ϵ::Real) = error("can not call to a delta function")

quadgkwrapper(f::DiracDelta) = ifelse(lowerbound(f) <= f.ω <= upperbound(f), f.α, 0.)

spectrumshift(m::DiracDelta, μ::Real) = DiracDelta(ω=m.ω+μ, α=m.α, lb=lowerbound(m)-μ, ub=upperbound(m)-μ)
