# Brute-force exact diagonalization on the full Fock space: matrix
# representations of the fermionic/bosonic ladder, density and occupation
# operators (dimension 2^L / d^L) and the thermal many-body density matrices
# (dm) built from them — the benchmark reference for the fast
# coefficient-matrix (cdm) methods in src/quadraticsystems.

# function spin_half_matrices()
# 	s_SP = Array{Float64, 2}([0 0; 1 0])
# 	s_SM = Array{Float64, 2}([0 1; 0 0])
# 	s_Z = Array{Float64, 2}([-1 0; 0 1])
# 	s_x = s_SP+s_SM
# 	s_y = -im*(s_SP-s_SM)
# 	JW = -s_Z
# 	n = Array{Float64, 2}([0 0; 0 1])
# 	return s_SP, s_SM, s_Z, JW, n
# end

# const σ₊, σ₋, σz, JW, n̂ = spin_half_matrices()

"""
	fermionadagoperator()
	fermionadagoperator(L, pos)

Return the matrix representation of a single fermionic creation operator ĉ†;
`fermionadagoperator(L, pos)` returns the creation operator at site `pos` in the occupation
number basis of `L` sites (including the Jordan-Wigner string).
"""
fermionadagoperator() = Array{Float64, 2}([0 0; 1 0])
"""
	fermionaoperator()
	fermionaoperator(L, pos)

Return the matrix representation of a single fermionic annihilation operator ĉ;
`fermionaoperator(L, pos)` returns the annihilation operator at site `pos` in the occupation
number basis of `L` sites.
"""
fermionaoperator() = adjoint(fermionadagoperator())
JWoperator() = Array{Float64, 2}([1 0; 0 -1])
"""
	fermiondensityoperator()
	fermiondensityoperator(L, pos)
	fermiondensityoperator(L)

Return the matrix representation of the fermionic density operator n̂ = ĉ†ĉ: the
single-site version, the version at site `pos`, and the total particle-number operator
of `L` sites.
"""
function fermiondensityoperator()
	adag = fermionadagoperator()
	return adag * adag'
end
function fermionadagoperator(L::Int, pos::Int)
	(1 <= pos <= L) || throw(BoundsError(1:L, pos))
	σ₊ = fermionadagoperator()
	JW = JWoperator()
	I2 = one(JW)
	ops = Vector{Matrix{Float64}}(undef, L)
	for i in 1:pos-1
		ops[i] = JW
	end
	ops[pos] = σ₊
	for i in pos+1:L
		ops[i] = I2
	end
	(L == 1) && return ops[1]
	return kron(ops...)
end
fermionaoperator(L::Int, pos::Int) = adjoint(fermionadagoperator(L, pos))
function fermiondensityoperator(L::Int, pos::Int)
	adag = fermionadagoperator(L, pos)
	return adag * adag'
end
"""
	fermionoccupationoperator(n)
	fermionoccupationoperator(L, pos, n)
	fermionoccupationoperator(ns)

Return the projection operator onto the occupation number `n` (0 or 1); the multi-site
version is the tensor product of the single-site projectors.
"""
function fermionoccupationoperator(n::Int)
	(n in (0, 1)) || throw(ArgumentError("occupation must be 0 or 1"))
	nop = fermiondensityoperator()
	(n == 1) ? nop : one(nop) - nop
end
function fermionoccupationoperator(L::Int, pos::Int, n::Int)
	(n in (0, 1)) || throw(ArgumentError("occupation must be 0 or 1"))
	nop = fermiondensityoperator(L, pos)
	(n == 1) ? nop : one(nop) - nop
end
function fermionoccupationoperator(ns::AbstractVector{Int})
	ops = [fermionoccupationoperator(n) for n in ns]
	return kron(ops...)
end

"""
	fermionoperator(L, term)
	fermionoperator(h::AbstractHamiltonian)

Map a quadratic/quartic term `term` or a Hamiltonian `h` to a matrix operator on the
fermionic Fock space (dimension 2^L).
"""
function fermionoperator(L::Int, h::AdagATerm{T}) where {T<:Number} 
	i, j = positions(h)
	m = fermionadagoperator(L, i) * fermionaoperator(L, j)
	return h.coeff * m
end
function fermionoperator(L::Int, h::AdagAdagTerm{T}) where {T<:Number} 
	i, j = positions(h)
	m = fermionadagoperator(L, i) * fermionadagoperator(L, j)
	return h.coeff * m
end
function fermionoperator(L::Int, h::AATerm{T}) where {T<:Number} 
	i, j = positions(h)
	m = fermionaoperator(L, i) * fermionaoperator(L, j)
	return h.coeff * m
end
function fermionoperator(L::Int, h::QuarticTerm{T}) where {T<:Number} 
	i, j, k, _l = positions(h)
	m = fermionadagoperator(L, i) * fermionadagoperator(L, j) * fermionaoperator(L, k) * fermionaoperator(L, _l)
	return h.coeff * m
end
function fermionoperator(h::AbstractHamiltonian{T}) where {T<:Number}
	L = num_sites(h)
	m = zeros(T, 2^L, 2^L)
	for t in h.data
		m .+= fermionoperator(L, t)
	end
	return m
end
function fermiondensityoperator(L::Int)
	m = zeros(Float64, 2^L, 2^L)
	for i in 1:L
		m .+= fermiondensityoperator(L, i)
	end
	return m
end
"""
	fermionicthermodm(h; β, μ=0)

Construct the fermionic thermal equilibrium (true) density matrix (dm)
exp(-β(Ĥ-μN̂))/Z of the Hamiltonian `h` at inverse temperature `β` and chemical
potential `μ`; the dm lives in the Fock space and can be used to evaluate the
expectation value of any observable, in contrast to the cdm which only encodes
quadratic observables (see `thermocdm`).
"""
function fermionicthermodm(h::AbstractHamiltonian; β::Real, μ::Real=0)
	m = fermionoperator(h)
	# @assert m ≈ m' atol=1.0e-12
	if μ != zero(μ)
		m = m - μ * fermiondensityoperator(num_sites(h))
	end
	# rho = exp(-β .* m)
	# rho ./= tr(rho)
	# return rho
	return thermodm(m, β=β)
end
"""
	thermodm(h, cache=eigencache(h); β)
	thermodm(cache; β)

Construct the thermal equilibrium (true) density matrix (dm) exp(-βh)/Z of a matrix
(Hamiltonian) `h` at inverse temperature `β`; `cache` is an eigendecomposition cache.
"""
thermodm(h::AbstractMatrix, cache::EigenCache=eigencache(h); β::Real) = thermodm(cache, β=β)
function thermodm(cache::EigenCache; β::Real)
	U, λs = cache.U, cache.λs
	λs2 = exp.(-β .* λs)
	rho = U * Diagonal(λs2) * U'
	rho ./= tr(rho)
	return rho
end


### bosonic operators
"""
	bosonaoperator(L, pos; d)
	bosonaoperator(; d)

Return the matrix representation of the bosonic annihilation operator â with truncation
dimension `d`; `bosonaoperator(L, pos; d)` returns the annihilation operator of mode `pos`
among `L` modes with truncation dimension `d`.
"""
function bosonaoperator(L::Int, pos::Int; d::Int)
	a = bosonaoperator(d=d)
	Ia = one(a)
	ops = [Ia for i in 1:L]
	ops[pos] = a
	return kron(ops...)
end
bosonadagoperator(L::Int, pos::Int; d::Int) = adjoint(bosonaoperator(L, pos, d=d))
"""
	bosondensityoperator(L, pos; d)
	bosondensityoperator(L; d)
	bosondensityoperator(; d)

Return the matrix representation of the bosonic density operator n̂ = â†â: for a single
mode, for mode `pos`, or the total particle-number operator of `L` modes.
"""
function bosondensityoperator(L::Int, pos::Int; d::Int)
	adag = bosonadagoperator(L, pos, d=d)
	return adag * adag'
end
"""
	bosonoccupationoperator(n; d)
	bosonoccupationoperator(L, pos, n; d)
	bosonoccupationoperator(ns; d)

Return the bosonic projection operator onto the occupation number `n` (0 ≤ n < d); the
multi-mode version is the tensor product of the single-mode projectors.
"""
function bosonoccupationoperator(L::Int, pos::Int, n::Int; d::Int)
	(0 <= n < d ) || throw(BoundsError(0:d-1, n))
	a = bosonoccupationoperator(n, d=d)
	Ia = one(a)
	ops = [Ia for i in 1:L]
	ops[pos] = a
	return kron(ops...)
end
function bosonoccupationoperator(ns::AbstractVector{Int}; d::Int)
	ops = [bosonoccupationoperator(n, d=d) for n in ns]
	return kron(ops...)
end

"""
	bosonoperator(L, term; d)
	bosonoperator(h::AbstractHamiltonian; d)

Map a quadratic/quartic term `term` or a Hamiltonian `h` to a matrix operator on the
bosonic Fock space (dimension d^L).
"""
function bosonoperator(L::Int, h::AdagATerm{T}; d::Int) where {T<:Number} 
	i, j = positions(h)
	m = bosonadagoperator(L, i, d=d) * bosonaoperator(L, j, d=d)
	return h.coeff * m
end
function bosonoperator(L::Int, h::AdagAdagTerm{T}; d::Int) where {T<:Number} 
	i, j = positions(h)
	m = bosonadagoperator(L, i, d=d) * bosonadagoperator(L, j, d=d)
	return h.coeff * m
end
function bosonoperator(L::Int, h::AATerm{T}; d::Int) where {T<:Number} 
	i, j = positions(h)
	m = bosonaoperator(L, i, d=d) * bosonaoperator(L, j, d=d)
	return h.coeff * m
end
function bosonoperator(L::Int, h::QuarticTerm{T}; d::Int) where {T<:Number} 
	i, j, k, _l = positions(h)
	m = bosonadagoperator(L, i, d=d) * bosonadagoperator(L, j, d=d) * bosonaoperator(L, k, d=d) * bosonaoperator(L, _l, d=d)
	return h.coeff * m
end
function bosonoperator(h::AbstractHamiltonian{T}; d::Int) where {T<:Number}
	L = num_sites(h)
	m = zeros(T, d^L, d^L)
	for t in h.data
		m .+= bosonoperator(L, t, d=d)
	end
	return m
end
function bosondensityoperator(L::Int; d::Int)
	m = zeros(Float64, d^L, d^L)
	for i in 1:L
		m .+= bosondensityoperator(L, i, d=d)
	end
	return m
end
"""
	bosonicthermodm(h; d, β, μ=0)

Construct the bosonic thermal equilibrium (true) density matrix (dm)
exp(-β(Ĥ-μN̂))/Z of the Hamiltonian `h` at inverse temperature `β` and chemical
potential `μ`, with truncation dimension `d` per mode.
"""
function bosonicthermodm(h::AbstractHamiltonian; d::Int, β::Real, μ::Real=0)
	m = bosonoperator(h, d=d)
	# @assert m ≈ m' atol=1.0e-12
	if μ != zero(μ)
		m = m - μ * bosondensityoperator(num_sites(h), d=d)
	end
	# rho = exp(-β .* m)
	# rho ./= tr(rho)
	# return rho
	return thermodm(m, β=β)
end


function bosonaoperator(; d::Int)
	(d <= 1) && error("d must be larger than 1")
	a = zeros(Float64, d, d)
	for i = 1:(d - 1)
		a[i, i+1] = sqrt(i)
	end
	return a
end
"""
	bosonadagoperator(L, pos; d)
	bosonadagoperator(; d)

Return the matrix representation of the bosonic creation operator â† with truncation
dimension `d`; `bosonadagoperator(L, pos; d)` returns the creation operator of mode `pos`
among `L` modes with truncation dimension `d`.
"""
bosonadagoperator(; d::Int) = adjoint(bosonaoperator(d=d))
function bosondensityoperator(; d::Int) 
	a = bosonaoperator(d=d)
	return a' * a
end
function bosonoccupationoperator(n::Int; d::Int)
	(0 <= n <= d-1) || throw(BoundsError(0:d-1, n))
	r = zeros(d, d)
	r[n+1, n+1] = 1
	return r
end


# function get_paulistring(h::AdagATerm{T}) where {T}
# 	r = Dict{Int, Matrix{T}}()
# 	i, j = positions(h)
# 	if i < j
# 		r[i] = σ₊
# 		r[j] = σ₋
# 		for k in i+1:j-1
# 			r[k] = I2
# 		end
# 	end
# end

# function fermionoperator(L::Int, opstr::AbstractDict{Int, <:AbstractMatrix})
# 	all(x->1<=x<=L, keys(opstr)) || throw(ArgumentError("operator position out of range"))
# 	isempty(opstr) && return zeros(2^L, 2^L)
# 	I2 = one(JW)
# 	ops = []
# 	for i in 1:L
# 		op = get(opstr, i, nothing)
# 		if isnothing(op)
# 			push!(ops, I2)
# 		else
# 			(size(op, 1) == size(op, 2) == 2) || throw(DimensionMismatch())
# 			push!(ops, op)
# 		end
# 	end
# 	return kron(ops...)
# end

Base.kron(a::AbstractMatrix) = a
