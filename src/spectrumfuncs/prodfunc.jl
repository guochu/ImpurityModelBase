Base.:*(x::AbstractBoundedFunction, y::Function) = bounded(ϵ->x(ϵ)*y(ϵ), lowerbound(x), upperbound(x))
Base.:/(x::AbstractBoundedFunction, y::Function) = bounded(ϵ->x(ϵ)/y(ϵ), lowerbound(x), upperbound(x))


"""
	struct DeltaMultF{F<:Function}

Delta function mult a spectrum function
"""
struct DeltaMultF{F<:Function} <: AbstractBoundedFunction
	δ::DiracDelta
	f::F
end
lowerbound(f::DeltaMultF) = lowerbound(f.δ)
upperbound(f::DeltaMultF) = upperbound(f.δ)

Base.:(-)(x::DeltaMultF) = DeltaMultF(x.δ, ϵ->-x.f(ϵ))
Base.adjoint(x::DeltaMultF) = DeltaMultF(x.δ, ϵ->conj(x.f(ϵ)))

Base.:*(x::DiracDelta, y::Function) = DeltaMultF(x, y)
Base.:/(x::DiracDelta, y::Function) = DeltaMultF(x, ϵ->1/y(ϵ))
Base.:*(x::DeltaMultF, y::Function) = DeltaMultF(x.δ, ϵ->x.f(ϵ)*y(ϵ))
Base.:/(x::DeltaMultF, y::Function) = DeltaMultF(x.δ, ϵ->x.f(ϵ)/y(ϵ))


function quadgkwrapper(m::DeltaMultF; kwargs...)
	ω₀, α = m.δ.ω, m.δ.α
	return ifelse(lowerbound(m) <= ω₀ <= upperbound(m), α * m.f(ω₀), 0.)
end 
