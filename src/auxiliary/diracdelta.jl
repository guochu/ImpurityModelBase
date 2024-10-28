"""
	struct DiracDelta

Dirac delta function
"""
struct DiracDelta <: Function
	ω₀::Float64
	α::Float64
end
DiracDelta(;ω₀::Real=0, α::Real=1) = DiracDelta(convert(Float64, ω₀), convert(Float64, α))
(x::DiracDelta)(ϵ::Real) = ifelse(x.ω₀ == ϵ, x.α, 0.)

_quadgk(f::DiracDelta, lb::Real, ub::Real) = ifelse(lb <= f.ω₀ <= ub, f.α, 0.)

struct DiracDeltaMultF{F} 
	δ::DiracDelta
	f::F
end


Base.:*(δ::DiracDelta, f::Function) = DiracDeltaMultF(δ, f)
Base.:*(f::Function, δ::DiracDelta) = DiracDeltaMultF(δ, f)


function _quadgk(m::DiracDeltaMultF, lb::Real, ub::Real)
	ω₀, α = m.δ.ω₀, m.δ.α
	return ifelse(lb <= ω₀ <= ub, α * m.f(ω₀), 0.)
end 

quadgkwrapper(m::DiracDeltaMultF{<:BoundedFunction}) = _quadgk(m, lowerbound(m.f), upperbound(m.f))