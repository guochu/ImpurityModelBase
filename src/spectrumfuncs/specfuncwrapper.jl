# the x-axis of boundary function is real!!!
abstract type AbstractBoundedFunction <: Function end
lowerbound(x::AbstractBoundedFunction) = x.lb
upperbound(x::AbstractBoundedFunction) = x.ub
(x::AbstractBoundedFunction)(ϵ::Real) = ifelse(lowerbound(x) <= ϵ <= upperbound(x), x.f(ϵ), 0.)
"""
	quadgkwrapper(m::AbstractBoundedFunction; kwargs...)

A simple wrapper of the quadgk function, integrates a single-variate function,
with a lowerbound and upperbound
"""
quadgkwrapper(m::AbstractBoundedFunction; kwargs...) = _quadgk(m.f, lowerbound(m), upperbound(m); kwargs...)

Base.:(-)(x::AbstractBoundedFunction) = bounded(ϵ->-x.f(ϵ), lowerbound(x), upperbound(x))
Base.adjoint(x::AbstractBoundedFunction) = bounded(ϵ->conj(x.f(ϵ)), lowerbound(x), upperbound(x))

"""
	BoundedFunction{F}

Wrapper for a general function, including a lower bound and upper bound
"""
struct BoundedFunction{F} <: AbstractBoundedFunction
	f::F
	lb::Float64
	ub::Float64
end
function BoundedFunction(f; lb::Real=-Inf, ub::Real=Inf) 
	(lb < ub) || throw(ArgumentError("lb must be less than ub"))
	return BoundedFunction(f, convert(Float64, lb), convert(Float64, ub))
end	
bounded(f, lb::Real, ub::Real) = BoundedFunction(f, lb=lb, ub=ub)
bounded(f; kwargs...) = BoundedFunction(f; kwargs...)
# Base.similar(x::BoundedFunction, f; lb::Real=x.lb, ub::Real=x.ub) = BoundedFunction(f, lb=lb, ub=ub)
# Base.similar(x::BoundedFunction; f=x.f, lb::Real=x.lb, ub::Real=x.ub) = BoundedFunction(f, lb=lb, ub=ub)

function spectrum(f, lb::Real, ub::Real) 
	check_spectrumfunction(f, lb, ub)
	return BoundedFunction(f, convert(Float64, lb), convert(Float64, ub))
end
spectrum(f; lb::Real=-Inf, ub::Real=Inf) = spectrum(f, lb, ub)
Base.similar(x::BoundedFunction, f; lb::Real=x.lb, ub::Real=x.ub) = BoundedFunction(f, lb=lb, ub=ub)
Base.similar(x::BoundedFunction; f=x.f, lb::Real=x.lb, ub::Real=x.ub) = BoundedFunction(f, lb=lb, ub=ub)

spectrumshift(m::BoundedFunction, μ::Real) = bounded(ϵ->m.f(ϵ+μ), lowerbound(m)-μ, upperbound(m)-μ)

function check_spectrumfunction(f, lb, ub)
	(lb < ub) || throw(ArgumentError("lb must be less than ub"))
	if lb == -Inf
		lb = -10000
	end
	if ub == Inf
		ub = 10000
	end
	r = f(rand(lb:ub))
	isa(r, Real) || throw(ArgumentError("the output of spectrum function should be real"))
	(r >= 0) || throw(ArgumentError("the output of spectrum function should be positive"))
end