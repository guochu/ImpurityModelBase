function toulouse_greater_lesser(b::AbstractDiscreteFermionicBath; ϵ_d::Real)
	h = toulouse_cmatrix(b, ϵ_d=ϵ_d)
	return freefermions_greater_lesser(h, β=b.β, μ=b.μ, i=1)
end
function toulouse_greater_lesser(b::AbstractDiscreteFermionicBath, ts::AbstractVector{<:Real}; kwargs...)
	f = toulouse_greater_lesser(b; kwargs...)
	r = f.(ts)
	return map(x->x[1], r), map(x->x[2], r)
end

"""
	freefermions_greater_lesser(h::AbstractMatrix, i::Int, j::Int; kwargs...) 
	freefermions_greater_lesser(h::AbstractMatrix, ρ₀::AbstractMatrix, i::Int, j::Int=i)
Real-time greater and lesser Green's functions
"""
function freefermions_greater_lesser(h::AbstractMatrix; β::Real, i::Int=1, j::Int=i, μ::Real=0) 
	@assert ishermitian(h)
	return _fermionic_eq_gf_util(h, i, j, β, μ)
end
function freefermions_greater_lesser(h::AbstractMatrix, ts::AbstractVector{<:Real}; kwargs...) 
	f = freefermions_greater_lesser(h; kwargs...)
	r = f.(ts)
	return map(x->x[1], r), map(x->x[2], r)
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


function freefermions_greater_lesser(h::AbstractMatrix, ρ₀::AbstractMatrix; i::Int=1, j::Int=i)
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
function freefermions_greater_lesser(h::AbstractMatrix, ρ₀::AbstractMatrix, ts::AbstractVector{<:Real}; kwargs...) 
	f1, f2 = freefermions_greater_lesser(h, ρ₀; kwargs...)
	return f1.(ts), f2.(ts)
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