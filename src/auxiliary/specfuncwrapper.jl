# the x-axis of boundary function is real!!!
abstract type AbstractBoundedFunction <: Function end
lowerbound(x::AbstractBoundedFunction) = x.lb
upperbound(x::AbstractBoundedFunction) = x.ub
(x::AbstractBoundedFunction)(ϵ::Real) = ifelse(lowerbound(x) <= ϵ <= upperbound(x), x.f(ϵ), 0.)
quadgkwrapper(m::AbstractBoundedFunction; kwargs...) = _quadgk(m.f, lowerbound(m), upperbound(m); kwargs...)


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
"""
struct SpectrumFunction{F} <: AbstractBoundedFunction
	f::F
	lb::Float64
	ub::Float64
end
function SpectrumFunction(f; lb::Real=-Inf, ub::Real=Inf) 
	(lb < ub) || throw(ArgumentError("lb must be less than ub"))
	return SpectrumFunction(f, convert(Float64, lb), convert(Float64, ub))
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