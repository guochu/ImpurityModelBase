const holstein_finiteT_rtol = 1.0e-8

"""
	holstein_G0w_to_Gw_finiteT(G0w::Function, ϵ::Real; β::Real, g::Real, ω::Real, maxiter::Int, order)

calculate G(ω) from G₀(ω) for the holstein model at finite T,
note that this analytical is only an approximate solution
"""
function holstein_G0w_to_Gw_finiteT(G0w::Function, ϵ::Real; β::Real, g::Real, ω::Real, maxiter::Int=10, rtol::Real=holstein_finiteT_rtol)
	(β == Inf) && error("use holstein_G0w_to_Gw_zeroT instead for T=0")
	c = 1 - exp(-β*ω)
	r = 0.
	order = 100
	for n in 0:order
		rj = exp(-β*n*ω) * holstein_G0w_to_Gw_finiteT_n(G0w, ϵ, g=g, ω=ω, n=n, maxiter=maxiter)
		r += rj
		if abs(rj / r) <= rtol
			# println("converged in $n iterations")
			break
		end
	end
	return c * r
end

function holstein_G0w_to_Gw_finiteT_n(G0w::Function, ϵ::Real; g::Real, ω::Real, n::Int, maxiter::Int=10)
	A = holstein_A(G0w, ϵ, g=g, ω=ω, n=n)
	B = holstein_B(G0w, ϵ, g=g, ω=ω, n=n, maxiter=maxiter)
	return 1 / (1/G0w(ϵ) - A - B)
end 

holstein_A(G0w::Function, ϵ::Real; kwargs...) = holstein_A_util(G0w, ϵ, 1; kwargs...)
function holstein_A_util(G0w::Function, ϵ::Real, iter::Int; g::Real, ω::Real, n::Int)
	(n == 0) && return 0.
	if iter == n
		return (n+1-iter) * g^2 * G0w(ϵ+iter*ω)
	else
		return (n+1-iter) * g^2 / (1/G0w(ϵ+iter*ω) - holstein_A_util(G0w, ϵ, iter+1, g=g, ω=ω, n=n) )
	end
end

holstein_B(G0w::Function, ϵ::Real; kwargs...) = holstein_B_util(G0w, ϵ, 1; kwargs...)
function holstein_B_util(G0w::Function, ϵ::Real, iter::Int; g::Real, ω::Real, n::Int, maxiter::Int=10)
	if iter >= maxiter
		return (n+iter) * g^2 * G0w(ϵ-iter*ω)
	else
		return (n+iter) * g^2 / (1/G0w(ϵ-iter*ω) - holstein_B_util(G0w, ϵ, iter+1, g=g, ω=ω, n=n, maxiter=maxiter) )
	end
end