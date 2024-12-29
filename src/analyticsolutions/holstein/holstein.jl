include("greenfunction.jl")
include("realtime.jl")

"""
	holstein_scaleless_parameters(; g::Real, ω::Real, t::Real=1)

Return scaleless parameters from bare parameters of the holstein model

t is the scale of the bath spectrum density
The Hamiltonian of the holstein model is
Ĥ = ϵ_d â†â + ωb̂†b̂ + g â†â (b̂† + b̂) + ∑ₖϵₖĉ†ĉ + ∑ₖVₖ(â†ĉ + ĉ†â)
where â is the fermionic impurity mode, b̂ is the bosonic mode
ĉ is the fermionic bath mode

The scaleless parameters are defined as
λ = g^2 / (ωt)
γ = ω / t
α = g / ω
"""
function holstein_scaleless_parameters(; g::Real, ω::Real, t::Real=1)
	λ = g^2 / (ω * t)
	γ = ω / t
	α = g / ω

	return (λ=λ, γ=γ, α=α)
end

"""
	holstein_bare_parameters(; λ::Real, γ::Real, t::Real=1)

Return bare parameters from scaleless parameters of the holstein model
"""
function holstein_bare_parameters(; λ::Real, γ::Real, t::Real=1)
	ω = γ * t
	g = sqrt(λ * ω * t)
	return (ω=ω, g=g)
end