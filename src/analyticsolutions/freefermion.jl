"""
	freefermion_greater(t::Real; β::Real, μ::Real)

Greater Green's function (G>) of a free particle in thermal state
μ = -ϵ_d

G>(t,t′) is defined as G>(t, t′) = -im⟨â(t)â†(t′)⟩
and G>(t) = G>(t, 0)
"""
freefermion_greater(t::Real; β::Real, μ::Real) = -im*exp(im*μ*t)/(1+exp(safe_mult(β, μ)))

"""
	freefermion_lesser(t::Real; β::Real, μ::Real)

Lesser Green's function (G<) of a free particle in thermal state
μ = -ϵ_d

G<(t,t′) is defined as G>(t, t′) = im⟨â†(t′)â(t)⟩
and G<(t) = G<(t, 0)

"""
freefermion_lesser(t::Real; β::Real, μ::Real) = im*exp(im*μ*t)/(exp(-safe_mult(β, μ))+1)

"""
	freefermion_Gt(t::Real; β::Real, μ::Real)

Retarded Green's (G) function of a free particle
μ = -ϵ_d

G(t,t′) is defined as G(t, t′) = G>(t,t′) - G<(t,t′) for t > t′
and G(t) = G(t, 0)
"""
freefermion_Gt(t::Real; kwargs...) = freefermion_greater(t; kwargs...) + freefermion_lesser(t; kwargs...)

"""
	freefermion_Gτ(τ::Float64; β::Real, μ::Real)

Matsubara Green's function of a free particle in the imaginary-time axis
μ = -ϵ_d

The Matsubara Green's function is defined as
g(τ, τ′) = -⟨â(τ)â†(τ′)⟩ for τ > τ′
and g(τ) = g(τ, 0)
"""
freefermion_Gτ(τ::Float64; β::Real, μ::Real) = exp(τ*μ)*(1 - freefermion_occupation(β, μ))

"""
	freefermion_occupation(β::Real, μ::Real) 

The average occupation n̄ = ⟨n̂⟩ of a free fermion in thermal state
"""
freefermion_occupation(β::Real, μ::Real) = 1 / (1+exp(-β*μ))


# interacting fermion site
"""
	fermion_greater(t::Real; β::Real, μ::Real, U::Real)

Greater Green's function of an interacting fermion in thermal state
μ = -ϵ_d
"""
fermion_greater(t::Real; β::Real, μ::Real, U::Real) = -im*_fermion_greater_util(im*t, β, μ, U)

"""
	fermion_lesser(t::Real; β::Real, μ::Real, U::Real)

Lesser Green's function of an interacting fermion in thermal state
μ = -ϵ_d
"""
fermion_lesser(t::Real; β::Real, μ::Real, U::Real) = im*_fermion_lesser_util(im*t, β, μ, U)

"""
	fermion_Gt(t::Real; β::Real, μ::Real, U::Real)

Retarded Green's function of an interacting fermion in thermal state
μ = -ϵ_d
"""
fermion_Gt(t::Real; kwargs...) = fermion_greater(t; kwargs...) + fermion_lesser(t; kwargs...)

"""
	fermion_Gτ(τ::Float64; β::Real, μ::Real, U::Real)

Matsubara Green's function of an interacting fermion in thermal state
μ = -ϵ_d
"""
fermion_Gτ(τ::Float64; β::Real, μ::Real, U::Real) = _fermion_greater_util(τ, β, μ, U)


function _fermion_greater_util(τ, β, μ, U)
    x = exp(-safe_mult(β, μ))
    return (exp(τ*μ) * x + exp(τ*(μ-U))) / (x + 2 + exp(safe_mult(β, μ-U)))
end

function _fermion_lesser_util(τ, β, μ, U)
    x = exp(-safe_mult(β, μ))
    y = exp(safe_mult(β, μ-U))
    return (exp(τ*μ) + y*exp(τ*(μ-U))) / (x + 2 + y)
end
