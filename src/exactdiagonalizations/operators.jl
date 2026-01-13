abstract type AbstractTerm{T<:Number} end
abstract type QuadraticTerm{T} <: AbstractTerm{T} end

positions(x::AbstractTerm) = x.positions
Base.eltype(::Type{<:AbstractTerm{T}}) where {T<:Number} = T
Base.eltype(x::AbstractTerm) = eltype(typeof(x))

Base.:*(m::Number, s::AbstractTerm) = s * m
Base.:/(s::AbstractTerm, m::Number) = s * (1 / m)
Base.:+(s::AbstractTerm) = s
Base.:-(s::AbstractTerm) = (-1) * s

Base.convert(::Type{M}, x::M) where {M<:AbstractTerm} = x


"""
	struct AdagATerm <: QuadraticTerm

AdagA term
"""
struct AdagATerm{T <: Number} <: QuadraticTerm{T}
	positions::Tuple{Int, Int}
	coeff::T
end

AdagATerm(pos::Tuple{Int, Int}; coeff::Number=1) = AdagATerm(pos, float(coeff))
AdagATerm(i::Int, j::Int; kwargs...) = AdagATerm((i, j); kwargs...)
tunneling(i::Int, j::Int; kwargs...) = AdagATerm(i, j; kwargs...)
adaga(i::Int, j::Int; kwargs...) = AdagATerm(i, j; kwargs...)

function Base.adjoint(x::AdagATerm)
	i, j = positions(x)
	return AdagATerm((j, i), coeff=conj(x.coeff))
end

Base.copy(x::AdagATerm) = AdagATerm(positions(x), copy(x.coeff))
Base.:*(s::AdagATerm, m::Number) = AdagATerm(positions(s), coeff=s.coeff * m)
Base.convert(::Type{AdagATerm{T}}, x::AdagATerm) where {T<:Number} = AdagATerm(positions(x), convert(T, x.coeff))


"""
	struct AdagAdagTerm <: QuadraticTerm

AdagA term
"""
struct AdagAdagTerm{T <: Number} <: QuadraticTerm{T}
	positions::Tuple{Int, Int}
	coeff::T
end

AdagAdagTerm(pos::Tuple{Int, Int}; coeff::Number=1) = AdagAdagTerm(pos, float(coeff))
AdagAdagTerm(i::Int, j::Int; kwargs...) = AdagAdagTerm((i, j); kwargs...)
adagadag(i::Int, j::Int; kwargs...) = AdagAdagTerm(i, j; kwargs...)

Base.copy(x::AdagAdagTerm) = AdagATerm(positions(x), copy(x.coeff))
Base.:*(s::AdagAdagTerm, m::Number) = AdagATerm(positions(s), coeff=s.coeff * m)
Base.convert(::Type{AdagAdagTerm{T}}, x::AdagAdagTerm) where {T<:Number} = AdagAdagTerm(positions(x), convert(T, x.coeff))


"""
	struct AATerm <: QuadraticTerm

AdagA term
"""
struct AATerm{T <: Number} <: QuadraticTerm{T}
	positions::Tuple{Int, Int}
	coeff::T
end

AATerm(pos::Tuple{Int, Int}; coeff::Number=1) = AATerm(pos, float(coeff))
AATerm(i::Int, j::Int; kwargs...) = AATerm((i, j); kwargs...)
aa(i::Int, j::Int; kwargs...) = AATerm(i, j; kwargs...)

Base.copy(x::AATerm) = AATerm(positions(x), copy(x.coeff))
Base.:*(s::AATerm, m::Number) = AATerm(positions(s), coeff=s.coeff * m)
Base.convert(::Type{AATerm{T}}, x::AATerm) where {T<:Number} = AATerm(positions(x), convert(T, x.coeff))


function Base.adjoint(x::AdagAdagTerm)
	i, j = positions(x)
	return AATerm((j, i), coeff=conj(x.coeff))
end
function Base.adjoint(x::AATerm)
	i, j = positions(x)
	return AdagAdagTerm((j, i), coeff=conj(x.coeff))
end

Base.convert(::Type{QuadraticTerm{T}}, x::AdagATerm) where {T} = convert(AdagATerm{T}, x)
Base.convert(::Type{QuadraticTerm{T}}, x::AdagAdagTerm) where {T} = convert(AdagAdagTerm{T}, x)
Base.convert(::Type{QuadraticTerm{T}}, x::AATerm) where {T} = convert(AATerm{T}, x)

# interacting term
"""
	struct FourBodyTerm <: AbstractFTerm

Fermionic fourbody term, ĉ₁†ĉ₂†ĉ₃ĉ₄
"""
struct QuarticTerm{T<:Number} <: AbstractTerm{T}
	positions::NTuple{4, Int}
	coeff::T
end

QuarticTerm(pos::NTuple{4, Int}; coeff::Real=1) = QuarticTerm(pos, float(coeff))
QuarticTerm(i::Int, j::Int, k::Int, l::Int; kwargs...) = QuarticTerm((i, j, k, l); kwargs...)
interaction(i::Int, j::Int, k::Int, l::Int; kwargs...) = QuarticTerm(i, j, k, l; kwargs...)

function Base.adjoint(x::QuarticTerm)
	i, j, k, l = positions(x)
	return interaction((l,k,j,i), coeff=conj(x.coeff))
end

Base.copy(x::QuarticTerm) = QuarticTerm(positions(x), copy(x.coeff))
Base.:*(s::QuarticTerm, m::Number) = QuarticTerm(positions(s), coeff=s.coeff * m)
Base.convert(::Type{QuarticTerm{T}}, x::QuarticTerm) where {T<:Number} = QuarticTerm(positions(x), convert(T, x.coeff))


const NormalTerm{T<:Number} = Union{AdagATerm{T}, QuarticTerm{T}}

abstract type AbstractHamiltonian{T<:Number} end
Base.eltype(::Type{<:AbstractHamiltonian{T}}) where {T<:Number} = T
Base.eltype(x::AbstractHamiltonian) = eltype(typeof(x))
abstract type QuadraticHamiltonian{T<:Number} <: AbstractHamiltonian{T} end
num_sites(x::AbstractHamiltonian) = x.n
function Base.:+(x::M, y::M) where {M<:AbstractHamiltonian}
	data = vcat(x.data, y.data)
	n = max(num_sites(x), num_sites(y))
	return M(data, n)
end

struct NormalHamiltonian{T<:Number} <: AbstractHamiltonian{T}
	data::Vector{NormalTerm{T}}
	n::Int
end
NormalHamiltonian(::Type{T}, n::Int) where {T<:Number} = NormalHamiltonian(Vector{NormalTerm{T}}(), n)
NormalHamiltonian(n::Int, x::Vector{<:NormalTerm{T}}) where {T} = NormalHamiltonian(convert(Vector{NormalTerm{T}}, x), n)

struct NormalQuadraticHamiltonian{T<:Number} <: QuadraticHamiltonian{T}
	data::Vector{AdagATerm{T}}
	n::Int
end
NormalQuadraticHamiltonian(::Type{T}, n::Int) where {T<:Number} = NormalQuadraticHamiltonian(Vector{AdagATerm{T}}(), n)
NormalQuadraticHamiltonian(n::Int, x::Vector{AdagATerm{T}}) where {T} = NormalQuadraticHamiltonian(x, n)

struct GenericQuadraticHamiltonian{T<:Number} <: QuadraticHamiltonian{T}
	data::Vector{QuadraticTerm{T}}
	n::Int
end
GenericQuadraticHamiltonian(::Type{T}, n::Int) where {T<:Number} = GenericQuadraticHamiltonian(Vector{QuadraticTerm{T}}(), n)
GenericQuadraticHamiltonian(n::Int, x::Vector{<:QuadraticTerm{T}}) where {T} = GenericQuadraticHamiltonian(convert(Vector{QuadraticTerm{T}}, x), n)

function Base.push!(x::AbstractHamiltonian, f::AbstractTerm)
	all(y->1<=y<=num_sites(x), positions(f)) || throw(BoundsError(1:num_sites(x), positions(f))) 
	push!(x.data, f)
end 

function quadratichamiltonian(n::Int, x::Vector{<:QuadraticTerm{T}}) where {T<:Number}
	if all(y->isa(y, AdagATerm), x)
		return NormalQuadraticHamiltonian(n, convert(Vector{AdagATerm{T}}, x))
	else
		return GenericQuadraticHamiltonian(n, x)
	end
end

# coefficient matrix of Quadratic Hamiltonians
function cmatrix(L::Int, h::AdagATerm{T}; normal::Bool=true) where {T<:Number} 
	i, j = positions(h)
	if normal
		m = zeros(T, L, L)
	else
		m = zeros(T, 2L, 2L)
	end
	m[i, j] = h.coeff
	return m
end
function cmatrix(L::Int, h::AdagAdagTerm{T}) where {T<:Number} 
	i, j = positions(h)
	m = zeros(T, 2L, 2L)
	m[i, L+j] = h.coeff
	return m
end
function cmatrix(L::Int, h::AATerm{T}) where {T<:Number} 
	i, j = positions(h)
	m = zeros(T, 2L, 2L)
	m[L+i, j] = h.coeff
	return m
end

function cmatrix(h::NormalQuadraticHamiltonian{T}) where {T<:Number}
	L = num_sites(h)
	m = zeros(T, L, L)
	for item in h.data
		i, j = positions(item)
		m[i, j] += item.coeff
	end
	return m
end
function cmatrix(h::GenericQuadraticHamiltonian{T}) where {T<:Number}
	L = num_sites(h) 
	m = zeros(T, 2L, 2L)
	for item in h.data
		i, j = positions(item)
		coef = item.coeff
		if isa(item, AdagATerm)
			m[i, j] += coef
			# m[L+i, L+j] += 1-coef
		elseif isa(item, AdagAdagTerm)
			m[i, L+j] += coef
		elseif isa(item, AATerm)
			m[L+i, j] += coef
		else
			error("unknown type $(typeof(item))")
		end
	end
	return bcs_symmetrize!(m)
end

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

fermionadagoperator() = Array{Float64, 2}([0 0; 1 0])
fermionaoperator() = adjoint(fermionadagoperator())
JWoperator() = Array{Float64, 2}([1 0; 0 -1])
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
function fermionoccupationoperator(L::Int, pos::Int, n::Int)
	(n in (0, 1)) || throw(ArgumentError("occupation must be 0 or 1"))
	nop = fermiondensityoperator(L, pos)
	(n == 1) ? nop : one(nop) - nop
end

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
thermodm(h::AbstractMatrix, cache::EigenCache=eigencache(h); β::Real) = thermodm(cache, β=β)
function thermodm(cache::EigenCache; β::Real)
	U, λs = cache.U, cache.λs
	λs2 = exp.(-β .* λs)
	rho = U * Diagonal(λs2) * U'
	rho ./= tr(rho)
	return rho
end


# bosonic operators
function bosonaoperator(; d::Int)
	(d <= 1) && error("d must be larger than 1.")
	a = zeros(Float64, d, d)
	for i = 1:(d - 1)
		a[i, i+1] = sqrt(i)
	end
	return a
end
bosonadagoperator(; d::Int) = adjoint(bosonaoperator(d=d))
function bosondensityoperator(; d::Int) 
	a = bosonaoperator(d=d)
	return a' * a
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

