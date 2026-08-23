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

"""
	spectrumshift(m::DeltaMultF, μ)

Shift the frequency axis of the product `m` of a δ function and a spectral function by `μ`.
"""
function spectrumshift(m::DeltaMultF, μ::Real)
	return DeltaMultF(spectrumshift(m.δ, μ), ϵ->m.f(ϵ+μ))
end

"""
	quadgkwrapper(m::DeltaMultF; kwargs...)

"Integrate" the product of a single δ peak and a spectral function `f`, giving
`α * f(ω₀)` (the δ peak is located at `ω₀` with strength `α`).
"""
function quadgkwrapper(m::DeltaMultF; kwargs...)
	ω₀, α = m.δ.ω, m.δ.α
	return ifelse(lowerbound(m) <= ω₀ <= upperbound(m), α * m.f(ω₀), 0.)
end 


"""
	struct DeltasMultF

Product of multiple δ peaks (a `DiscreteSpectrum`) and a spectral function `f`.
"""
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

"""
	spectrumshift(m::DeltasMultF, μ)

Shift the frequency axis of the discrete spectrum product `m` by `μ`.
"""
function spectrumshift(m::DeltasMultF, μ::Real)
	return DeltasMultF(spectrumshift(m.δ, μ), ϵ->m.f(ϵ+μ))
end

"""
	quadgkwrapper(m::DeltasMultF; kwargs...)

"Integrate" the product of a discrete spectrum and a spectral function `f`, giving
∑ₙ αₙ f(ωₙ).
"""
function quadgkwrapper(m::DeltasMultF; kwargs...)
	ws, αs = frequencies(m.δ), spectrumvalues(m.δ)
	r = 0.
	for (w, α) in zip(ws, αs)
		r += α * m.f(w)
	end
	return r
end 
