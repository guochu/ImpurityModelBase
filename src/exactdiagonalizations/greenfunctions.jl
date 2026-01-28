# fermions
function freefermions_Gt(h::AbstractMatrix, i::Int, j::Int=i, cache::EigenCache=eigencache(h); kwargs...)
	f = freefermions_greater_lesser(h, i, j, cache; kwargs...)
	function f′(t)
		x1, x2 = f(t)
		return x1 - x2
	end
	return f′
end
function freefermions_Gt(h::AbstractMatrix, ts::AbstractVector{<:Real}, i::Int, j::Int=i, cache::EigenCache=eigencache(h); kwargs...)
	f = freefermions_Gt(h, i, j, cache; kwargs...)
	return f.(ts)
end


"""
	freefermions_greater_lesser(h::AbstractMatrix, i::Int, j::Int; kwargs...) 
	freefermions_greater_lesser(h::AbstractMatrix, ρ₀::AbstractMatrix, i::Int, j::Int=i)
Real-time greater and lesser Green's functions
"""
function freefermions_greater_lesser(h::AbstractMatrix, i::Int, j::Int=i, cache::EigenCache=eigencache(h); β::Real, μ::Real=0) 
	# @assert ishermitian(h)
	return _fermionic_eq_gf_util(cache, i, j, β, μ)
end
function freefermions_greater_lesser(h::AbstractMatrix, ts::AbstractVector{<:Real}, i::Int, j::Int=i, cache::EigenCache=eigencache(h); kwargs...) 
	f = freefermions_greater_lesser(h, i, j, cache; kwargs...)
	r = f.(ts)
	return map(x->x[1], r), map(x->x[2], r)
end

function _fermionic_eq_gf_util(cache::EigenCache, i::Int, j::Int, β::Real, μ::Real)
	λs, U = cache.λs, cache.U
	L = size(cache.m, 1)
	function f(t::Number)
		r_g = zero(ComplexF64)
		r_l = zero(ComplexF64)
		for k in 1:L
			ss = U[i, k] * conj(U[j, k]) * exp(-im * λs[k] * t)
			n_k = fermidirac(β, μ, λs[k])
			r_g += ss * (1 - n_k) 
			r_l += ss * n_k 
		end
		return -im * r_g, im * r_l
	end
	return f	
end


# ### neq Green's functions
# function freefermions_greater_lesser(h::AbstractMatrix, ρ₀::AbstractMatrix; i::Int=1, j::Int=i)
# 	(size(h) == size(ρ₀)) || throw(DimensionMismatch("Hamiltonian size mismatch with density matrix"))
# 	cache = eigencache(h)
# 	U = cache.U
# 	# ns = [real(ρ₀[k, k]) for k in 1:size(U, 1)]
# 	ρ₀′ = U' * ρ₀ * U
# 	# ρ₁ = one(ρ₀) - transpose(ρ₀′) 
# 	# ρ₂ = transpose(ρ₀′)
# 	ρ₁ = one(ρ₀) - ρ₀′
# 	ρ₂ = ρ₀′
# 	return t -> -im*_fermionic_neq_gf_util(cache, ρ₁, i, j, t), t -> im*_fermionic_neq_gf_util(cache, ρ₂, i, j, t)
# end
# function freefermions_greater_lesser(h::AbstractMatrix, ρ₀::AbstractMatrix, ts::AbstractVector{<:Real}; kwargs...) 
# 	f1, f2 = freefermions_greater_lesser(h, ρ₀; kwargs...)
# 	return f1.(ts), f2.(ts)
# end

# # HU = Uλs
# # U brings H into diagonal
# # ns are diagonal terms of ρ in the diagonal representation of H
# function _fermionic_neq_gf_util(cache::EigenCache, ρ, i::Int, j::Int, t::Real)
# 	U, λs = cache.U, cache.λs
# 	r_g = zero(eltype(U))
# 	L = size(U, 1)
# 	for k in 1:L
# 		for k′ in 1:L
# 			r_g += U[i, k] * conj(U[j, k′]) * ρ[k, k′] * exp(-im * λs[k] * t)
# 		end
# 	end
# 	return r_g
# end


function freefermions_Gt(h::AbstractMatrix, ρ₀::AbstractMatrix, i::Int, j::Int=i, cache::EigenCache=eigencache(h))
	f = freefermions_greater_lesser(h, ρ₀, i, j, cache)
	function f′(t)
		x1, x2 = f(t)
		return x1 - x2
	end
	return f′	
end
function freefermions_Gt(h::AbstractMatrix, ρ₀::AbstractMatrix, ts::AbstractVector{<:Real}, 
									 i::Int=1, j::Int=i, cache::EigenCache=eigencache(h))
	f = freefermions_Gt(h, ρ₀, i, j, cache)
	return f.(ts)
end

"""
	freefermions_greater_lesser(h::AbstractMatrix, ρ₀::AbstractMatrix, i::Int, j::Int=i, cache::EigenCache=eigencache(h))

Applies for both normal and bcs baths
"""
function freefermions_greater_lesser(h::AbstractMatrix, ρ₀::AbstractMatrix, i::Int, j::Int=i, cache::EigenCache=eigencache(h))
	(size(h) == size(ρ₀)) || throw(DimensionMismatch("Hamiltonian size mismatch with density matrix"))
	U = cache.U
	ρ₁ = one(ρ₀) - transpose(ρ₀)
	ρ₂′ = U' * transpose(ρ₀ ) * U
	# ρ₁ = one(ρ₀) - transpose(ρ₀′) 
	# ρ₂ = transpose(ρ₀′)
	ρ₁′ = U' * ρ₁ * U
	return t -> _fermionic_neq_gf_util(cache, ρ₁′, ρ₂′, i, j, t)
end
function freefermions_greater_lesser(h::AbstractMatrix, ρ₀::AbstractMatrix, ts::AbstractVector{<:Real}, 
									 i::Int=1, j::Int=i, cache::EigenCache=eigencache(h)) 
	f = freefermions_greater_lesser(h, ρ₀, i, j, cache)
	r = f.(ts)
	return map(x->x[1], r), map(x->x[2], r)
end

function _fermionic_neq_gf_util(cache::EigenCache, ρ₁, ρ₂, i::Int, j::Int, t::Real)
	U, λs = cache.U, cache.λs
	r_g = zero(eltype(U))
	r_l = zero(eltype(U))
	L = size(U, 1)
	for k in 1:L
		exp_t = exp(-im * λs[k] * t)
		for k′ in 1:L
			ss = U[i, k] * conj(U[j, k′]) * exp_t
			r_g +=  ss * ρ₁[k, k′] 
			# ss = U[i, k] * conj(U[j, k′])
			r_l += ss * ρ₂[k, k′]
		end
	end
	return -im * r_g, im * r_l
end


# Matsubara Green's functions
function freefermions_Gτ(h::AbstractMatrix, i::Int, j::Int=i, cache::EigenCache=eigencache(h); β::Real)
	λs, U = cache.λs, cache.U
	function f(τ)
		r_g = zero(eltype(U))
		for k in 1:size(U, 1)
			ss = U[i, k] * conj(U[j, k]) 
			n_k = fermidirac(β, 0, -λs[k])
			r_g += ss * (n_k * exp(-λs[k] * τ))
		end
		return r_g
	end
	return f
end
function freefermions_Gτ(h::AbstractMatrix, τs::AbstractVector{<:Real}, i::Int, j::Int=i, cache::EigenCache=eigencache(h); β::Real)
	f = freefermions_Gτ(h, i, j, cache, β=β)
	return f.(τs)
end


### bosons

function freebosons_Gt(h::AbstractMatrix, i::Int, j::Int=i, cache::EigenCache=eigencache(h); kwargs...)
	f = freebosons_greater_lesser(h, i, j, cache; kwargs...)
	function f′(t)
		x1, x2 = f(t)
		return x1 - x2
	end
	return f′
end
function freebosons_Gt(h::AbstractMatrix, ts::AbstractVector{<:Real}, i::Int, j::Int=i, cache::EigenCache=eigencache(h); kwargs...)
	f = freebosons_Gt(h, i, j, cache; kwargs...)
	return f.(ts)
end


"""
	freefermions_greater_lesser(h::AbstractMatrix, i::Int, j::Int; kwargs...) 
	freefermions_greater_lesser(h::AbstractMatrix, ρ₀::AbstractMatrix, i::Int, j::Int=i)
Real-time greater and lesser Green's functions
"""
function freebosons_greater_lesser(h::AbstractMatrix, i::Int, j::Int=i, cache::EigenCache=eigencache(h); β::Real, μ::Real=0) 
	# @assert ishermitian(h)
	return _bosonic_eq_gf_util(cache, i, j, β, μ)
end
function freebosons_greater_lesser(h::AbstractMatrix, ts::AbstractVector{<:Real}, i::Int, j::Int=i, cache::EigenCache=eigencache(h); kwargs...) 
	f = freebosons_greater_lesser(h, i, j, cache; kwargs...)
	r = f.(ts)
	return map(x->x[1], r), map(x->x[2], r)
end

function _bosonic_eq_gf_util(cache::EigenCache, i::Int, j::Int, β::Real, μ::Real)
	λs, U = cache.λs, cache.U
	L = size(cache.m, 1)
	function f(t::Number)
		r_g = zero(ComplexF64)
		r_l = zero(ComplexF64)
		for k in 1:L
			ss = U[i, k] * conj(U[j, k]) * exp(-im * λs[k] * t)
			n_k = boseeinstein(β, μ, λs[k])
			r_g += ss * (1 + n_k) 
			r_l += ss * n_k 
		end
		return -im * r_g, im * r_l
	end
	return f	
end

function freebosons_Gt(h::AbstractMatrix, ρ₀::AbstractMatrix, i::Int, j::Int=i, cache::EigenCache=eigencache(h))
	f = freebosons_greater_lesser(h, ρ₀, i, j, cache)
	function f′(t)
		x1, x2 = f(t)
		return x1 - x2
	end
	return f′	
end
function freebosons_Gt(h::AbstractMatrix, ρ₀::AbstractMatrix, ts::AbstractVector{<:Real}, 
									 i::Int=1, j::Int=i, cache::EigenCache=eigencache(h))
	f = freebosons_Gt(h, ρ₀, i, j, cache)
	return f.(ts)
end


function freebosons_greater_lesser(h::AbstractMatrix, ρ₀::AbstractMatrix, i::Int, j::Int=i, cache::EigenCache=eigencache(h))
	(size(h) == size(ρ₀)) || throw(DimensionMismatch("Hamiltonian size mismatch with density matrix"))
	U = cache.U
	ρ₁ = one(ρ₀) + transpose(ρ₀)
	ρ₂′ = U' * transpose(ρ₀ ) * U
	# ρ₁ = one(ρ₀) - transpose(ρ₀′) 
	# ρ₂ = transpose(ρ₀′)
	ρ₁′ = U' * ρ₁ * U
	return t -> _fermionic_neq_gf_util(cache, ρ₁′, ρ₂′, i, j, t)
end
function freebosons_greater_lesser(h::AbstractMatrix, ρ₀::AbstractMatrix, ts::AbstractVector{<:Real}, 
									 i::Int=1, j::Int=i, cache::EigenCache=eigencache(h)) 
	f = freebosons_greater_lesser(h, ρ₀, i, j, cache)
	r = f.(ts)
	return map(x->x[1], r), map(x->x[2], r)
end


# Matsubara Green's functions
function freebosons_Gτ(h::AbstractMatrix, i::Int, j::Int=i, cache::EigenCache=eigencache(h); β::Real)
	λs, U = cache.λs, cache.U
	function f(τ)
		r_g = zero(eltype(U))
		for k in 1:size(U, 1)
			ss = U[i, k] * conj(U[j, k]) 
			n_k = 1 + boseeinstein(β, 0, λs[k])
			r_g += ss * (n_k * exp(-λs[k] * τ))
		end
		return r_g
	end
	return f
end
function freebosons_Gτ(h::AbstractMatrix, τs::AbstractVector{<:Real}, i::Int, j::Int=i, cache::EigenCache=eigencache(h); β::Real)
	f = freebosons_Gτ(h, i, j, cache, β=β)
	return f.(τs)
end