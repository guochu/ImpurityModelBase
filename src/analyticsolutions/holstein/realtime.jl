include("zeroT.jl")
include("finiteT.jl")

function holstein_G0w_to_Σw(G0w::Function, ϵ::Real; g::Real, ω::Real, β::Real=Inf, maxiter::Int=10, rtol::Real=holstein_finiteT_rtol)
	if β == Inf
		return holstein_G0w_to_Σw_zeroT(G0w, ϵ, g=g, ω=ω, maxiter=maxiter)
	else
		error("holstein_G0w_to_Σw not implemented for zero temperature")
	end
end

function holstein_G0w_to_Gw(G0w::Function, ϵ::Real; g::Real, ω::Real, β::Real=Inf, maxiter::Int=10, rtol::Real=holstein_finiteT_rtol)
	if β == Inf
		return holstein_G0w_to_Gw_zeroT(G0w, ϵ, g=g, ω=ω, maxiter=maxiter)
	else
		return holstein_G0w_to_Gw_finiteT(G0w, ϵ, β=β, g=g, ω=ω, maxiter=maxiter, rtol=rtol)
	end
end


function bethe_holstein_dmft_iteration(G0w::Function, wk::Real; g::Real, ω::Real, β::Real=Inf, t::Real=1, maxiter::Int=10, rtol::Real=holstein_finiteT_rtol)
	Σk = holstein_G0w_to_Σw(G0w, wk, β=β, g=g, ω=ω, maxiter=maxiter, rtol=rtol) 
	a = -(Σk-wk)
	b = sqrt((Σk-wk)^2-t^2)
	Gk = (a + b) / (t^2/2)
	if imag(Gk) > 0
		Gk = (a - b) / (t^2/2)
	end
	return Gk
end

function holstein_Gt(f::AbstractSpectrumFunction, t::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, 
						β::Real=Inf, wmax::Real=20., wmin::Real=-wmax, δ::Real=1.0e-8, maxiter::Int=10, rtol::Real=holstein_finiteT_rtol)
    A = quadgkwrapper(bounded(ϵ -> (holstein_Gw(f, ϵ; β=β, ϵ_d=ϵ_d, g=g, ω=ω, μ=μ, δ=δ, maxiter=maxiter, rtol=rtol)-1.0/(ϵ+im*δ))*exp(-im*ϵ*t), wmin, wmax))
    return A/(2π)-im
end

function holstein_Gw(f::AbstractSpectrumFunction, ϵ::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, β::Real=Inf, 
						δ::Real=1.0e-8, maxiter::Int=10, rtol::Real=holstein_finiteT_rtol)
	G0w(y) = toulouse_Gw(f, y, ϵ_d=ϵ_d, μ=μ, δ=δ)
	return holstein_G0w_to_Gw(G0w, ϵ, g=g, ω=ω, β=β, maxiter=maxiter, rtol=rtol)
end