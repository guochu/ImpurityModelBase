
"""
	spinboson_dephasingdynamics(spectrum::AbstractBoundedFunction, t::Real, ρ₀::AbstractMatrix; β::Real, Δ::Real=0)

Exact dephasing dynamics of the spin-boson model in the absence of the tunneling term
(no σₓ coupling). The Hamiltonian reads

H = Δσz + σz ∑ₖ Vₖ (aₖ + aₖ†) + ∑ₖ ωₖ aₖ† aₖ

where `spectrum` is the bath spectrum density (an `AbstractBoundedFunction`), `β` the
inverse temperature and `Δ` the energy splitting of the spin.

Given the 2×2 reduced density matrix `ρ₀` of the spin at time t=0, return the evolved
2×2 density matrix at time `t`. Since only σz couples to the bosonic bath, the
populations (diagonal elements) are unchanged and only the coherences (off-diagonal
elements) acquire a renormalized phase factor, which gives rise to the dephasing
dynamics.
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

# DD sequence XX
"""
	ddxx_spinboson_dephasingdynamics(spectrum, N, ρ₀; β, δt)

Evolve the spin-boson dephasing model under an XX dynamical-decoupling (DD) pulse sequence,
returning the evolved 2×2 density matrix. `N` is an even number of DD steps and `δt` the
time interval per step.
"""
function ddxx_spinboson_dephasingdynamics(spectrum::AbstractBoundedFunction, N::Int, ρ₀::AbstractMatrix; β::Real, δt::Real)
	iseven(N) || throw("Even number of DD steps assumed")
	(size(ρ₀, 1) == size(ρ₀, 2) == 2) || throw(ArgumentError("initial state should be a 2×2 density matrix"))
	ρout = Matrix{ComplexF64}(ρ₀)
	c = _renormalized_dd_phase(spectrum, δt, N, β)
	# println("c is ", c)
	ρout[1,2] = ρout[1,2] * c
	ρout[2,1] = conj(ρout[1,2])
	return ρout
end


function _renormalized_dd_phase(f::AbstractBoundedFunction, δt, N, β)
    phase = 0.
    for j in 1:N
    	tmp = _jj(f, δt, β)
    	for k in 1:j-1
    		sgn = isodd(j+k) ? -1 : 1
    		tmp += _jk(f, δt, j-k, β) * sgn
    	end
    	phase += tmp
    end
    return exp(-phase)
end


function _jj(f::AbstractBoundedFunction, δt, β)
	g(ω) = (coth(β*ω/2)/ω^2) * (1-cos(ω*δt))
	return quadgkwrapper(f * g)
end

function _jk(f::AbstractBoundedFunction, δt, Δk, β)
	g(ω) = (coth(β*ω/2)/ω^2) * (1-cos(ω*δt)) * cos(ω*Δk*δt)
	return 2 * quadgkwrapper(f * g)
end