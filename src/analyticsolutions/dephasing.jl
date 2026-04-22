
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

# DD sequence XX
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