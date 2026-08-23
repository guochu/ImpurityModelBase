# utilities for exact time evolution

# do h^t ρ - ρ h^t
# time evolution for the coefficient matrix rho of free fermions is dρ/dt = -i [h^t, ρ]
# h is the coefficient matrix
# h = [h₁₁ c₁†c₁, h₁₂ c₁†c₂, h₁₃ c₁†c₃...; h₂₁ c₂†c₁, h₂₂ c₂†c₂, h₂₃ c₂†c₃; ...]
"""
	struct EigenCache{T, R}

Eigendecomposition cache for a Hermitian matrix, storing the matrix `m`, the eigenvector
matrix `U` and the eigenvalues `λs`.
"""
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
"""
	eigencache(h)

Compute the eigendecomposition of the Hermitian matrix `h` and return an `EigenCache`.
"""
eigencache(h::AbstractMatrix) = EigenCache(h)
Base.conj(x::EigenCache) = EigenCache(conj(x.m), conj(x.U), x.λs)

# initializers
# do transpose here

# freefermions_cache(h::AbstractMatrix) = eigencache(transpose(h)) 
# cdmcache(h::AbstractMatrix) = eigencache(transpose(h))

"""
	freefermions_timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Real, cache=freefermions_cache(h))

Return quadratic observables at time t
"""
# freefermions_timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Real, cache::EigenCache=freefermions_cache(h)) = _generic_ed_timeevo(ρ₀, h, -im*t, cache)

# function freefermions_timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Real, cache::EigenCache=freefermions_cache(h))
# 	t2 = -im*t
# 	λs = [exp(λ*t2) for λ in cache.λs]
# 	exp_h = cache.U * Diagonal(λs) * adjoint(cache.U)
# 	return exp_h * ρ₀ * exp_h'
# end

"""
	timeevo(ρ₀, h, t, cache=eigencache(h))

Evolve the initial matrix `ρ₀` under the Hamiltonian `h`, returning
`exp(t·h) · ρ₀ · exp(t·h)†`; for real-time evolution use `t = -im·time`.
"""
timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Number, cache::EigenCache=eigencache(h)) = _generic_ed_timeevo(ρ₀, h, t, cache)
# itimeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, τ::Real, cache::EigenCache=eigencache(h)) = _generic_ed_timeevo(ρ₀, h, -τ, cache)

# operator_timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Real, cache::EigenCache=eigencache(h)) = _generic_ed_timeevo(ρ₀, h, im*t, cache)
# operator_itimeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, τ::Real, cache::EigenCache=eigencache(h)) = _generic_ed_timeevo(ρ₀, h, τ, cache)


function _generic_ed_timeevo(ρ₀::AbstractMatrix, h::AbstractMatrix, t::Number, cache::EigenCache)
	λs = [exp(λ*t) for λ in cache.λs]
	exp_h = cache.U * Diagonal(λs) * adjoint(cache.U)
	return exp_h * ρ₀ * exp_h'
end

"""
	thermocdm(::Type{P}, h, cache=eigencache(h); β, μ=0)
	thermocdm(::Type{P}, cache; β, μ=0)

Construct the thermal equilibrium coefficient density matrix (cdm) of a free-particle
system (`P` is `Boson` or `Fermion`) at inverse temperature `β` and chemical potential
`μ`.

Concept: in contrast to the density matrix (dm), which is the genuine many-body state
ρ = exp(-β(Ĥ-μN̂))/Z living in the exponentially large Fock/Hilbert space, the
coefficient density matrix (cdm) is an L×L single-particle matrix (L being the number
of single-particle states) encoding the single-particle correlations ⟨c_j† c_i⟩ of the
thermal Gaussian state. Any quadratic observable A = ∑_ij A_ij c_i† c_j can be evaluated
directly from the cdm, e.g. ⟨A⟩ = tr(cdm * A) with the appropriate convention, without
ever constructing the many-body density matrix.
"""
thermocdm(::Type{P}, h::AbstractMatrix, cache::EigenCache=eigencache(h); kwargs...) where {P <: AbstractParticle} = thermocdm(P, cache; kwargs...)

function thermocdm(::Type{P}, cache::EigenCache; β::Real, μ::Real=0) where {P <: AbstractParticle}
	U, λs = cache.U, cache.λs
	# println("eigenvalues...")
	# println(cache.m)
	# println(λs)
	# λs2 = exp.(-β .* λs)
	n = [thermaloccupation(P, β, μ, ϵ) for ϵ in λs]
	# λs2 ./= sum(λs2)
	# println("occupations....")
	# println(n)
	return transpose(U * Diagonal(n) * U')

end
"""
	fermionicthermocdm(cache; β, μ=0)
	bosonicthermocdm(cache; β, μ=0)

Construct the thermal equilibrium coefficient density matrix (cdm) of a fermionic/bosonic
free system at inverse temperature `β` and chemical potential `μ`; see `thermocdm` for
the concept of the cdm (single-particle correlation matrix) versus the true density
matrix.
"""
fermionicthermocdm(cache::EigenCache; kwargs...) = thermocdm(Fermion, cache; kwargs...)
bosonicthermocdm(cache::EigenCache; kwargs...) = thermocdm(Boson, cache; kwargs...)

