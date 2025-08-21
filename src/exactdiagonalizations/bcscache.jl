function bcsthermocdm(cache::EigenCache; β::Real, μ::Real=0)
	(μ == 0) || throw(ArgumentError("BCS bath should have μ=0"))
	U, λs = cache.U, cache.λs
	# println("eigenvalues...")
	# println(cache.m)
	# println(λs)

	L0 = size(U, 1)
	L = div(L0, 2)
	(2L == L0) || throw(ArgumentError("not a BCS cmatrix"))
	n = zeros(eltype(λs), L0)
	for i in 1:L
		ϵ = λs[i] - λs[2L+1-i]
		# println("energy is ", ϵ)
		n[i] = thermaloccupation(Fermion, β, μ, ϵ)
	end
	for i in 1:L
		n[2L+1-i] = 1 - n[i]
	end
	return transpose(U * Diagonal(n) * U')
end

function bcs_cmatrix(h::AbstractMatrix{<:Number}, g::AbstractMatrix{<:Number})
	(size(h) == size(g)) || throw(DimensionMismatch())
	# ishermitian(h) || throw(ArgumentError("h matrix should be hermitian"))
	L = size(h, 1)
	T = promote_type(eltype(h), eltype(g))
	m = zeros(T, 2L, 2L)
	m[1:L, 1:L] = h
	m[1:L, L+1:2L] = g
	m[L+1:2L, 1:L] = g'
	return bcs_symmetrize!(m)
end

function bcs_cdm(h::AbstractMatrix, g::AbstractMatrix)
	(size(h) == size(g)) || throw(DimensionMismatch())
	L = size(h, 1)
	T = promote_type(eltype(h), eltype(g))
	m = zeros(T, 2L, 2L)
	m[1:L, 1:L] = h
	m[1:L, L+1:2L] = g
	m[L+1:2L, 1:L] = g'
	m[L+1:2L, L+1:2L] = one(h) - transpose(h)
	return m
end

"""
	bcs_cdmcache(h::AbstractMatrix)
"""
bcs_cdmcache(h::AbstractMatrix) = eigencache(2*transpose(h))

function bcs_symmetrize!(h::AbstractMatrix)
	L0 = size(h, 1)
	L = div(L0, 2)
	(2L == L0) || throw(ArgumentError("not a BCS cmatrix"))
	m = (h[1:L, 1:L] - transpose(h[L+1:2L, L+1:2L])) /2
	h[1:L, 1:L] = m
	h[L+1:2L, L+1:2L] = -transpose(m)
	m = h[1:L, L+1:2L]
	h[1:L, L+1:2L] = antisymmetrize!(m)
	h[L+1:2L, 1:L] = m'
	return h
end

function antisymmetrize!(m::AbstractMatrix)
	L = size(m, 1)
	for i in 1:L
		m[i, i] = 0
		for j in i+1:L
			tmp = (m[i, j] - m[j, i]) / 2
			m[i, j] = tmp
			m[j, i] = -tmp
		end
	end
	return m
end

function bcs_hmatrix(m::AbstractMatrix)
	L0 = size(h, 1)
	L = div(L0, 2)
	(2L == L0) || throw(ArgumentError("not a BCS cmatrix"))
	return m[1:L, 1:L]	
end
function bcs_gmatrix(m::AbstractMatrix)
	L0 = size(h, 1)
	L = div(L0, 2)
	(2L == L0) || throw(ArgumentError("not a BCS cmatrix"))
	return m[1:L, L+1:2L]		
end