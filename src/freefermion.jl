"""
	freefermion_greater(t::Real; β::Real, μ::Real)

Greater Green's function of a free particle in thermal state
μ = -ϵ_d
"""
freefermion_greater(t::Real; β::Real, μ::Real) = -im*exp(im*μ*t)/(1+exp(β*μ))

"""
	freefermion_lesser(t::Real; β::Real, μ::Real)

Lesser Green's function of a free particle
"""
freefermion_lesser(t::Real; β::Real, μ::Real) = im*exp(im*μ*t)/(exp(-β*μ)+1)

"""
	freefermion_Gt(t::Real; β::Real, μ::Real)

Retarded Green's function of a free particle
"""
freefermion_Gt(t::Real; kwargs...) = freefermion_greater(t; kwargs...) + freefermion_lesser(t; kwargs...)

"""
	freefermion_Gτ(τ::Float64; β::Real, μ::Real)

Matsubara Green's function of a free particle in the imaginary-time axis
"""
freefermion_Gτ(τ::Float64; β::Real, μ::Real) = exp(τ*μ)*(1 - freefermion_occupation(β, μ))

freefermion_occupation(β::Real, μ::Real) = 1 / (1+exp(-β*μ))


# interacting fermion site
"""
	fermion_greater(t::Real; β::Real, μ::Real, U::Real)

Greater Green's function of an interacting fermion in thermal state
μ = -ϵ_d
"""
fermion_greater(t::Real; β::Real, μ::Real, U::Real) = -im*_fermion_greater_util(im*t, β, μ, U)
fermion_lesser(t::Real; β::Real, μ::Real, U::Real) = im*_fermion_lesser_util(im*t, β, μ, U)
fermion_Gt(t::Real; kwargs...) = fermion_greater(t; kwargs...) + fermion_lesser(t; kwargs...)
fermion_Gτ(τ::Float64; β::Real, μ::Real, U::Real) = _fermion_greater_util(τ, β, μ, U)


function _fermion_greater_util(τ, β, μ, U)
    x = exp(-β * μ)
    return (exp(τ*μ) * x + exp(τ*(μ-U))) / (x + 2 + exp(β*(μ-U)))
end

function _fermion_lesser_util(τ, β, μ, U)
    x = exp(-β * μ)
    y = exp(β*(μ-U))
    return (exp(τ*μ) + y*exp(τ*(μ-U))) / (x + 2 + y)
end