"""
	Gτ(h::FermionicQuadraticHamiltonian, i::Int, j::Int; β::Real, μ::Real)

The Matsubara Green's function 
Gᵢⱼ(τ, τ′) = -Tr[âᵢ(τ)â†ⱼ(τ′)exp(-β(Ĥ-μN̂)⟩]
Gᵢⱼ(τ) = Gᵢⱼ(τ, 0)
"""
Gτ(m::FermionicQuadraticHamiltonian, i::Int, j::Int=i; kwargs...) = freefermions_Gτ(cmatrix(m), i, j; kwargs...)

function freefermions_Gτ(h::AbstractMatrix, i::Int, j::Int=i; β::Real, μ::Real=0)
	@assert ishermitian(h)
	ham = h - μ .* LinearAlgebra.I
	λs, U = eigen(Hermitian(ham))
	ns = [fermidirac(β, 0, λs[k]) for k in 1:length(λs)]
	return τ -> _fermionic_Gτ_util(U, λs, ns, i, j, τ)
end

function _fermionic_Gτ_util(U, λs, ns, i::Int, j::Int, t::Real)
	r_g = zero(eltype(U))
	for k in 1:size(U, 1)
		ss = U[i, k] * conj(U[j, k])
		exp_t = exp(-λs[k] * t)
		n_k = ns[k]
		r_g += ss * (1 - n_k) * exp_t
	end
	return r_g
end



"""
	greater_lesser(m::FermionicQuadraticHamiltonian, i::Int, j::Int; kwargs...) 
	greater_lesser(m::FreeFermionicHamiltonian, ρ₀::AbstractMatrix, i::Int, j::Int)
FermionicQuadraticHamiltonian
Real-time greater and lesser Green's functions
"""
greater_lesser(m::FermionicQuadraticHamiltonian, i::Int, j::Int=i; kwargs...) = freefermions_greater_lesser(cmatrix(m), i, j; kwargs...)
greater_lesser(m::FermionicQuadraticHamiltonian, ρ₀::AbstractMatrix, i::Int, j::Int=i) = freefermions_greater_lesser(cmatrix(m), ρ₀, i, j)


function freefermions_greater_lesser(h::AbstractMatrix, i::Int, j::Int=i; β::Real, μ::Real=0) 
	@assert ishermitian(h)
	return _fermionic_eq_gf_util(h, i, j, β, μ)
end

function _fermionic_eq_gf_util(h::AbstractMatrix, i::Int, j::Int, β::Real, μ::Real)
	ham = convert(Matrix{eltype(h)}, h)
	evals, evecs = eigen(Hermitian(ham))
	L = size(ham, 1)
	function f(t::Number)
		r_g = zero(ComplexF64)
		r_l = zero(ComplexF64)
		for k in 1:L
			ss = evecs[i, k] * conj(evecs[j, k])
			exp_t = exp(-im * evals[k] * t)
			n_k = fermidirac(β, μ, evals[k])
			r_g += ss * (1 - n_k) * exp_t
			r_l += ss * n_k * exp_t
		end
		return -im * r_g, im * r_l
	end
	return f	
end


function freefermions_greater_lesser(h::AbstractMatrix, ρ₀::AbstractMatrix, i::Int, j::Int=i)
	(size(h) == size(ρ₀)) || throw(DimensionMismatch("Hamiltonian size mismatch with density matrix"))
	@assert ishermitian(h)
	ham = convert(Matrix{eltype(h)}, h)
	λs, U = eigen(Hermitian(ham))
	# ns = [real(ρ₀[k, k]) for k in 1:size(U, 1)]
	ρ₀′ = U' * ρ₀ * U
	# ρ₁ = one(ρ₀) - transpose(ρ₀′) 
	# ρ₂ = transpose(ρ₀′)
	ρ₁ = one(ρ₀) - ρ₀′
	ρ₂ = ρ₀′
	return t -> -im*_fermionic_neq_gf_util(U, λs, ρ₁, i, j, t), t -> im*_fermionic_neq_gf_util(U, λs, ρ₂, i, j, t)
end

# HU = Uλs
# U brings H into diagonal
# ns are diagonal terms of ρ in the diagonal representation of H
function _fermionic_neq_gf_util(U, λs, ρ, i::Int, j::Int, t::Real)
	r_g = zero(eltype(U))
	L = size(U, 1)
	for k in 1:L
		for k′ in 1:L
			r_g += U[i, k] * conj(U[j, k′]) * ρ[k, k′] * exp(-im * λs[k] * t)
		end
	end
	return r_g
end