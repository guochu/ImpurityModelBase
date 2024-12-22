# the x-axis of boundary function is real!!!
abstract type AbstractBoundedFunction <: Function end
abstract type AbstractSpectrumFunction <: AbstractBoundedFunction end
lowerbound(x::AbstractBoundedFunction) = x.lb
upperbound(x::AbstractBoundedFunction) = x.ub
(x::AbstractBoundedFunction)(ϵ::Real) = ifelse(lowerbound(x) <= ϵ <= upperbound(x), x.f(ϵ), 0.)
"""
	quadgkwrapper(m::AbstractBoundedFunction; kwargs...)

A simple wrapper of the quadgk function, integrates a single-variate function,
with a lowerbound and upperbound
"""
quadgkwrapper(m::AbstractBoundedFunction; kwargs...) = _quadgk(m.f, lowerbound(m), upperbound(m); kwargs...)

"""
	spectrumshift(m::AbstractSpectrumFunction, μ::Real)

Return a new spectrum by shift the x-axis by μ
"""
spectrumshift(m::AbstractSpectrumFunction, μ::Real) = error("spectrumshift not implemented for spectrum function type $(typeof(m))")

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
	return SpectrumFunction(f, convert(Float64, lb), convert(Float64, ub))
end	
bounded(f, lb::Real, ub::Real) = BoundedFunction(f, lb=lb, ub=ub)
bounded(f; kwargs...) = BoundedFunction(f; kwargs...)


"""
	SpectrumFunction{F}

Wrapper for spectrum function, including a lower bound and upper bound
Compared to bounded function, there are other requirements for spectrum 
function, such as that it should be positive
"""
struct SpectrumFunction{F} <: AbstractSpectrumFunction
	f::F
	lb::Float64
	ub::Float64
end
function SpectrumFunction(f; lb::Real=-Inf, ub::Real=Inf) 
	check_spectrumfunction(f, lb, ub)
	return SpectrumFunction(f, convert(Float64, lb), convert(Float64, ub))
end	
spectrum(f, lb::Real, ub::Real) = SpectrumFunction(f, lb=lb, ub=ub)
spectrum(f; kwargs...) = SpectrumFunction(f; kwargs...)

spectrumshift(m::SpectrumFunction, μ::Real) = spectrum(ϵ->m.f(ϵ+μ), lowerbound(m)-μ, upperbound(m)-μ)

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