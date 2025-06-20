# utilities for exact time evolution

# do h^t ρ - ρ h^t
# time evolution for the coefficient matrix rho of free fermions is dρ/dt = -i [h^t, ρ]
# h is the coefficient matrix
# h = [h₁₁ c₁†c₁, h₁₂ c₁†c₂, h₁₃ c₁†c₃...; h₂₁ c₂†c₁, h₂₂ c₂†c₂, h₂₃ c₂†c₃; ...]
struct EigenCache{T<:Number, R<:Real}
	m::Matrix{T}
	U::Matrix{T}
	λs::Vector{R}
end

function EigenCache(h::AbstractMatrix{T}) where {T<:Number}
	ishermitian(h) || throw(ArgumentError("EigenCache requires Hermitian matrix"))
	h2 = convert(Matrix{T}, h)
	λs, U = eigen(Hermitian(h2))
	return EigenCache(h2, U, λs)
end
eigencache(h::AbstractMatrix) = EigenCache(h)


# initializers
# do transpose here

"""
	freefermions_cache(h::AbstractMatrix)

h is the coefficient matrix, namely, 
h = [h₁₁ c₁†c₁, h₁₂ c₁†c₂, h₁₃ c₁†c₃...; h₂₁ c₂†c₁, h₂₂ c₂†c₂, h₂₃ c₂†c₃; ...]

time evolution for the coefficient matrix h of free fermions is dρ/dt = -i [h^t, ρ],
where ρ is the quadratic observables
"""
freefermions_cache(h::AbstractMatrix) = eigencache(transpose(h)) 

"""
	freefermions_timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Real, cache=freefermions_cache(h))

Return quadratic observables at time t
"""
freefermions_timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Real, cache::EigenCache=freefermions_cache(h)) = _generic_ed_timeevo(ρ₀, h, -im*t, cache)

# function freefermions_timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Real, cache::EigenCache=freefermions_cache(h))
# 	t2 = -im*t
# 	λs = [exp(λ*t2) for λ in cache.λs]
# 	exp_h = cache.U * Diagonal(λs) * adjoint(cache.U)
# 	return exp_h * ρ₀ * exp_h'
# end

timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Real, cache::EigenCache=eigencache(h)) = _generic_ed_timeevo(ρ₀, h, -im*t, cache)
itimeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, τ::Real, cache::EigenCache=eigencache(h)) = _generic_ed_timeevo(ρ₀, h, -τ, cache)

operator_timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Real, cache::EigenCache=eigencache(h)) = _generic_ed_timeevo(ρ₀, h, im*t, cache)
operator_itimeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, τ::Real, cache::EigenCache=eigencache(h)) = _generic_ed_timeevo(ρ₀, h, τ, cache)


function _generic_ed_timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Number, cache::EigenCache)
	λs = [exp(λ*t) for λ in cache.λs]
	exp_h = cache.U * Diagonal(λs) * adjoint(cache.U)
	return exp_h * ρ₀ * exp_h'
end