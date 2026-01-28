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