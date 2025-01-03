# utilities for exact time evolution

# do h^t ρ - ρ h^t
# time evolution for the coefficient matrix rho of free fermions is dρ/dt = -i [h^t, ρ]
# h is the coefficient matrix
# h = [h₁₁ c₁†c₁, h₁₂ c₁†c₂, h₁₃ c₁†c₃...; h₂₁ c₂†c₁, h₂₂ c₂†c₂, h₂₃ c₂†c₃; ...]
struct QuadraticCMatrixCache
	h::Matrix{ComplexF64}
	U::Matrix{ComplexF64}
	λs::Vector{Float64}
end

function QuadraticCMatrixCache(h::AbstractMatrix)
	@assert ishermitian(h)
	h2 = convert(Matrix{ComplexF64}, h)
	λs, U = eigen(Hermitian(h2))
	return QuadraticCMatrixCache(h2, U, λs)
end


# initializers
# do transpose here

"""
	freefermions_cache(h::AbstractMatrix)

h is the coefficient matrix, namely, 
h = [h₁₁ c₁†c₁, h₁₂ c₁†c₂, h₁₃ c₁†c₃...; h₂₁ c₂†c₁, h₂₂ c₂†c₂, h₂₃ c₂†c₃; ...]

time evolution for the coefficient matrix h of free fermions is dρ/dt = -i [h^t, ρ],
where ρ is the quadratic observables
"""
freefermions_cache(h::AbstractMatrix) = QuadraticCMatrixCache(convert(Matrix{ComplexF64}, transpose(h)))

"""
	freefermions_timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Real, cache=freefermions_cache(h))

Return quadratic observables at time t
"""
function freefermions_timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Real, cache::QuadraticCMatrixCache=freefermions_cache(h))
	t2 = -im*t
	λs = [exp(λ*t2) for λ in cache.λs]
	exp_h = cache.U * Diagonal(λs) * adjoint(cache.U)
	return exp_h * ρ₀ * exp_h'
end