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

function spectrumshift(m::DeltaMultF, μ::Real)
	return DeltaMultF(spectrumshift(m.δ, μ), ϵ->m.f(ϵ+μ))
end

function quadgkwrapper(m::DeltaMultF; kwargs...)
	ω₀, α = m.δ.ω, m.δ.α
	return ifelse(lowerbound(m) <= ω₀ <= upperbound(m), α * m.f(ω₀), 0.)
end 


struct DeltasMultF{F<:Function} <: AbstractBoundedFunction
	δ::DiscreteSpectrum
	f::F
end
lowerbound(f::DeltasMultF) = lowerbound(f.δ)
upperbound(f::DeltasMultF) = upperbound(f.δ)

Base.:(-)(x::DeltasMultF) = DeltasMultF(x.δ, ϵ->-x.f(ϵ))
Base.adjoint(x::DeltasMultF) = DeltasMultF(x.δ, ϵ->conj(x.f(ϵ)))

Base.:*(x::DiscreteSpectrum, y::Function) = DeltasMultF(x, y)
Base.:/(x::DiscreteSpectrum, y::Function) = DeltasMultF(x, ϵ->1/y(ϵ))
Base.:*(x::DeltasMultF, y::Function) = DeltasMultF(x.δ, ϵ->x.f(ϵ)*y(ϵ))
Base.:/(x::DeltasMultF, y::Function) = DeltasMultF(x.δ, ϵ->x.f(ϵ)/y(ϵ))

function spectrumshift(m::DeltasMultF, μ::Real)
	return DeltasMultF(spectrumshift(m.δ, μ), ϵ->m.f(ϵ+μ))
end

function quadgkwrapper(m::DeltasMultF; kwargs...)
	ws, αs = frequencies(m.δ), spectrumvalues(m.δ)
	r = 0.
	for (w, α) in zip(ws, αs)
		r += α * m.f(w)
	end
	return r
end 
