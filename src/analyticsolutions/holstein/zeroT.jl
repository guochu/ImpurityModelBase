# G(ω) for the Holstein model

function holstein_Σw_zeroT(f::AbstractSpectrumFunction, ϵ::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8, maxiter::Int=10)
	G0w(y) = toulouse_Gw(f, y, ϵ_d=ϵ_d, μ=μ, δ=δ)
	return holstein_G0w_to_Σw_zeroT(G0w, ϵ, g=g, ω=ω, maxiter)
end

"""
	holstein_G0w_to_Gw_zeroT(G0w::Function, ϵ::Real; β::Real, g::Real, ω::Real, maxiter::Int)

calculate G(ω) from G₀(ω) for the holstein model at T=0
"""
holstein_G0w_to_Gw_zeroT(G0w::Function, ϵ::Real; kwargs...) = 1 / (1/G0w(ϵ) - holstein_G0w_to_Σw_zeroT(G0w, ϵ; kwargs...))

"""
	holstein_G0w_to_Σw_zeroT(G0w::Function, ϵ::Real; β::Real, g::Real, ω::Real, maxiter::Int)

calculate Σ(ω) from G₀(ω) for the holstein model at T=0
"""
holstein_G0w_to_Σw_zeroT(G0w::Function, ϵ::Real; kwargs...) = holstein_G0w_to_Σw_zeroT_util(G0w, ϵ, 1; kwargs...)

function holstein_G0w_to_Σw_zeroT_util(G0w::Function, ϵ::Real, iter::Int; g::Real, ω::Real, maxiter::Int=10)
	if iter >= maxiter
		return iter * g^2 * G0w(ϵ - iter *ω)
	else
		return iter * g^2 / (1 / G0w(ϵ - iter *ω) - holstein_G0w_to_Σw_zeroT_util(G0w, ϵ, iter+1, g=g, ω=ω, maxiter=maxiter))
	end
end

# function holstein_G0w_t0_Gw_zeroT_util(G0w::Function, ϵ::Real, maxiter::Int, x; g::Real, ω::Real)
# 	deno = (maxiter == 0) ? 1 : maxiter * g^2
# 	return deno / (1/G0w(ϵ-maxiter*ω) - x)
# end


