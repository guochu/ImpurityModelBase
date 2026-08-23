using LinearAlgebra

# utilities functions
function random_normalquadratichamiltonian(m::AbstractMatrix) 
	L = size(m, 1)
	ham = NormalQuadraticHamiltonian(eltype(m), L)
	for i in 1:L, j in 1:L
		t = adaga(i, j, coeff=m[i, j])
		push!(ham, t)
	end
	return ham
end

function random_genericquadratichamiltonian(m::AbstractMatrix, m2::AbstractMatrix) 
	L = size(m, 1)
	ham = GenericQuadraticHamiltonian(eltype(m), L)
	for i in 1:L, j in 1:L
		t = adaga(i, j, coeff=m[i, j])
		push!(ham, t)
		t = adagadag(i, j, coeff=m2[i, j])
		if m2[i, j] != 0
			push!(ham, t)
			push!(ham, t')	
		end
	end
	return ham
end

function random_dm(::Type{T}, L::Int) where {T<:Number}
	dm = randn(T, L, L)
	dm = dm * dm'

	dm ./= tr(dm)
	return dm
end

function random_hermitian(::Type{T}, L::Int) where {T<:Number}
	m = randn(T, L, L)
	return m + m'
end

function normal_quadratic_obs(dm::AbstractMatrix)
	L = round(Int, log2(size(dm, 1)))
	obs = zeros(eltype(dm), L, L)
	tr_dm = tr(dm)
	for i in 1:L, j in 1:L
		t = adaga(i, j)
		op = fermionoperator(L, t)
		obs[i, j] = tr(op * dm) / tr_dm
	end
	return obs
end

function prod_boson_dm(L::Int; d::Int) 
	ns = [1 for i in 1:L]
	for i in 2:2:L
		ns[i] = 0
	end
	return bosonoccupationoperator(ns, d=d)
end

function boson_normal_quadratic_obs(dm::AbstractMatrix; d::Int)
	L = round(Int, log(d, size(dm, 1)))
	obs = zeros(eltype(dm), L, L)
	tr_dm = tr(dm)
	for i in 1:L, j in 1:L
		t = adaga(i, j)
		op = bosonoperator(L, t, d=d)
		obs[i, j] = tr(op * dm) / tr_dm
	end
	return obs
end

function generic_quadratic_obs(dm::AbstractMatrix)
	L = round(Int, log2(size(dm, 1)))
	obs = zeros(eltype(dm), L, L)
	tr_dm = tr(dm)
	for i in 1:L, j in 1:L
		t = adaga(i, j)
		op = fermionoperator(L, t)
		obs[i, j] = tr(op * dm) / tr_dm
	end
	obs2 = zeros(eltype(dm), L, L)
	obs3 = zeros(eltype(dm), L, L)
	for i in 1:L, j in 1:L
		t = adagadag(i, j)
		op = fermionoperator(L, t)
		obs2[i, j] = tr(op * dm) / tr_dm

		t = aa(i, j)
		op = fermionoperator(L, t)
		obs3[i, j] = tr(op * dm) / tr_dm
	end	
	(obs2 ≈ obs3') || throw(ArgumentError("something wrong"))
	return bcs_cdm(obs, obs2)	
end

# --- full-Hamiltonian reference helpers (merged from ed.jl) ---
# These "brute force" full-Hamiltonian Green's functions are the second
# (inefficient) method used to cross-check the efficient coefficient-matrix
# method; they are only practical for small systems and are used for debugging.

function spin_half_matrices()
	s_SP = Array{Float64, 2}([0 0; 1 0])
	s_SM = Array{Float64, 2}([0 1; 0 0])
	s_Z = Array{Float64, 2}([-1 0; 0 1])
	s_x = s_SP+s_SM
	s_y = -im*(s_SP-s_SM)
	n = Array{Float64, 2}([0 0; 0 1])
	return Dict("x"=>s_x, "y"=>s_y, "z"=>s_Z, "+"=>s_SP, "-"=>s_SM, "n"=>n)
end

function Aop(d::Int)
	(d <= 1) && error("d must be larger than 1.")
	a = zeros(Float64, d, d)
	for i = 1:(d - 1)
		a[i, i+1] = sqrt(i)
	end
	return a
end

ADAGop(d::Int) = Array(transpose(Aop(d)))

Nop(d::Int) = ADAGop(d) * Aop(d)


function boson(;d::Int=5)
	_N = Nop(d)
	_N2 = _N * _N
	return Dict("a"=>Aop(d),"adag"=>ADAGop(d), "n"=>_N, "n2"=>_N2)
end

# <e^τH A e^-τH B>
function gf_imag(H, A, B, β::Real, n::Int)
	δτ = β / n
	τs = 0:δτ:β
	λs, U = eigen(Hermitian(H))
	ρ = U * Diagonal(exp.(-β .* λs)) * U'
	tr_ρ = tr(ρ)
	g(τ) = tr(U * Diagonal(exp.(τ .* λs)) * U' * A * U * Diagonal(exp.(-τ .* λs)) * U' * B * ρ) / tr_ρ
	return g.(τs)
end

# <e^iHt A e^-iHt B>
function gf_real(H, A, B, β::Real, t::Real, n::Int, ρ=exp(-β*H))
	δt = t / n
	ts = 0:δt:t
	λs, U = eigen(Hermitian(H))
	# ρ = U * Diagonal(exp.(-β .* λs)) * U'
	tr_ρ = tr(ρ)
	gt(tj) = -im*tr(U * Diagonal(exp.(im*tj .* λs)) * U' * A * U * Diagonal(exp.(-im*tj .* λs)) * U' * B * ρ) / tr_ρ
	ls(tj) = im*tr(B * U * Diagonal(exp.(im*tj .* λs)) * U' * A * U * Diagonal(exp.(-im*tj .* λs)) * U' * ρ) / tr_ρ
	return gt.(ts), ls.(ts)
end