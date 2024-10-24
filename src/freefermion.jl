"""
	freefermion_greater(t::Real; β::Real, μ::Real)

Greater Green's function of a free particle
"""
freefermion_greater(t::Real; β::Real, μ::Real) = 1.0/(1+exp(β*μ))*exp(im*μ*t)

"""
	freefermion_lesser(t::Real; β::Real, μ::Real)

Lesser Green's function of a free particle
"""
freefermion_lesser(t::Real; β::Real, μ::Real) = 1.0/(exp(-β*μ)+1)*exp(im*μ*t)

"""
	freefermion_Gt(t::Real; β::Real, μ::Real)

Retarded Green's function of a free particle
"""
freefermion_Gt(t::Real; kwargs...) = freefermion_greater(t; kwargs...) + freefermion_lesser(t; kwargs...)

"""
	freefermion_Gτ(τ::Float64; β::Real, μ::Real)

Matsubara Green's function of a free particle in the imaginary-time axis
"""
freefermion_Gτ(τ::Float64; β::Real, μ::Real) = exp(τ*μ)/(1+exp(β*μ))

freefermion_occupation(β::Real, μ::Real) = 1 / (1+exp(-β*μ))