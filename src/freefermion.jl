"""
	free_greater(t::Real; β::Real, μ::Real)

Greater Green's function of a free particle
"""
free_greater(t::Real; β::Real, μ::Real) = 1.0/(1+exp(β*μ))*exp(im*μ*t)

"""
	free_lesser(t::Real; β::Real, μ::Real)

Lesser Green's function of a free particle
"""
free_lesser(t::Real; β::Real, μ::Real) = 1.0/(exp(-β*μ)+1)*exp(im*μ*t)

"""
	free_greater(t::Real; β::Real, μ::Real)

Retarded Green's function of a free particle
"""
free_Gt(t::Real; kwargs...) = free_greater(t; kwargs...) + free_lesser(t; kwargs...)

"""
	free_Gτ(τ::Float64; β::Real, μ::Real)

Matsubara Green's function of a free particle in the imaginary-time axis
"""
free_Gτ(τ::Float64; β::Real, μ::Real) = exp(τ*μ)/(1+exp(β*μ))
