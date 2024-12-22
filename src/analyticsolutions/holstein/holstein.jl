include("realtime.jl")

function holstein_scaleless_parameters(; g::Real, ω::Real, t::Real=1)
	λ = g^2 / (ω * t)
	γ = ω / t
	α = g / ω

	return (λ=λ, γ=γ, α=α)
end

function holstein_bare_parameters(; λ::Real, γ::Real, t::Real=1)
	ω = γ * t
	g = sqrt(λ * ω * t)
	return (ω=ω, g=g)
end