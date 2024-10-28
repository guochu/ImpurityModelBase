"""
	struct DiracDelta

Dirac delta function
"""
struct DiracDelta <: AbstractBoundedFunction
	ω₀::Float64
	α::Float64
	lb::Float64
	ub::Float64
end
function DiracDelta(;ω₀::Real=0, α::Real=1, lb::Real=-Inf, ub::Real=Inf)
	(lb < ub) || throw(ArgumentError("lb must be less than ub"))
	return DiracDelta(convert(Float64, ω₀), convert(Float64, α), convert(Float64, lb), convert(Float64, ub))
end
(x::DiracDelta)(ϵ::Real) = error("can not call to a delta function")

quadgkwrapper(f::DiracDelta) = ifelse(lowerbound(f) <= f.ω₀ <= upperbound(f), f.α, 0.)