
"""
	spinboson_dephasingdynamics(spectrum::AbstractBoundedFunction, t::Real; β::Real, Δ::Real)

H = Δσz + σz∑ₖVₖ(aₖ + aₖ+) + ∑ₖωₖaₖ+aₖ
"""
function spinboson_dephasingdynamics(spectrum::AbstractBoundedFunction, t::Real, ρ₀::AbstractMatrix; β::Real, Δ::Real=0)
	(size(ρ₀, 1) == size(ρ₀, 2) == 2) || throw(ArgumentError("initial state should be a 2×2 density matrix"))
	ρout = Matrix{ComplexF64}(ρ₀)
	c = _renormalized_phase(spectrum, t, β, Δ)
	# println("c is ", c)
	ρout[1,2] = ρout[1,2] * c
	ρout[2,1] = conj(ρout[1,2])
	return ρout
end



function _renormalized_phase(f::AbstractBoundedFunction, t, β, Δ)
    g(ω) = (coth(β*ω/2)/ω^2) * (1 - cos(ω*t))
    _e = quadgkwrapper(f * g)
    phase = -im*Δ*t - _e
    return exp(phase)
end