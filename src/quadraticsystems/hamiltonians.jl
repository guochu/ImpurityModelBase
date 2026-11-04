# Quadratic Hamiltonians: the term types (`AdagATerm`, `AdagAdagTerm`, `AATerm`,
# `QuarticTerm`), the Hamiltonian containers and their coefficient-matrix
# (`cmatrix`) representation — the single-particle basis on which all fast
# quadratic-systems methods in this folder are built.

"""
	AbstractTerm{T}

Abstract type for a term in a Hamiltonian; the type parameter `T` is the coefficient type.
"""
abstract type AbstractTerm{T<:Number} end
"""
	QuadraticTerm{T}

Abstract type for a quadratic term (`AdagATerm`, `AdagAdagTerm`, `AATerm`), corresponding
to a free-particle Hamiltonian.
"""
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
"""
	tunneling(i, j; coeff=1)
	adaga(i, j; coeff=1)

Construct an `AdagATerm` (âᵢ†âⱼ term); `tunneling` and `adaga` are equivalent.
"""
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
"""
	adagadag(i, j; coeff=1)

Construct an `AdagAdagTerm` (âᵢ†âⱼ† term).
"""
adagadag(i::Int, j::Int; kwargs...) = AdagAdagTerm(i, j; kwargs...)

Base.copy(x::AdagAdagTerm) = AdagAdagTerm(positions(x), copy(x.coeff))
Base.:*(s::AdagAdagTerm, m::Number) = AdagAdagTerm(positions(s), coeff=s.coeff * m)
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
"""
	aa(i, j; coeff=1)

Construct an `AATerm` (âᵢâⱼ term).
"""
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
"""
	interaction(i, j, k, l; coeff=1)

Construct a `QuarticTerm` (four-body interaction term, e.g. ĉᵢ†ĉⱼ†ĉₖĉₗ for fermions).
"""
interaction(i::Int, j::Int, k::Int, l::Int; kwargs...) = QuarticTerm(i, j, k, l; kwargs...)

function Base.adjoint(x::QuarticTerm)
	i, j, k, l = positions(x)
	return interaction(l,k,j,i, coeff=conj(x.coeff))
end

Base.copy(x::QuarticTerm) = QuarticTerm(positions(x), copy(x.coeff))
Base.:*(s::QuarticTerm, m::Number) = QuarticTerm(positions(s), coeff=s.coeff * m)
Base.convert(::Type{QuarticTerm{T}}, x::QuarticTerm) where {T<:Number} = QuarticTerm(positions(x), convert(T, x.coeff))


"""
	NormalTerm{T}

Union type for normal (â†â-type only) terms, equivalent to
`Union{AdagATerm{T}, QuarticTerm{T}}`.
"""
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

"""
	struct NormalHamiltonian{T}

General Hamiltonian (possibly interacting), composed of `NormalTerm`s (â†â terms and
quartic terms).

	NormalHamiltonian(n, data)

Construct a Hamiltonian with `n` sites and the term list `data`.
"""
struct NormalHamiltonian{T<:Number} <: AbstractHamiltonian{T}
	data::Vector{NormalTerm{T}}
	n::Int
end
NormalHamiltonian(::Type{T}, n::Int) where {T<:Number} = NormalHamiltonian(Vector{NormalTerm{T}}(), n)
NormalHamiltonian(n::Int, x::Vector{<:NormalTerm{T}}) where {T} = NormalHamiltonian(convert(Vector{NormalTerm{T}}, x), n)

"""
	struct NormalQuadraticHamiltonian{T}

Quadratic Hamiltonian containing only âᵢ†âⱼ terms (e.g. a normal free particle/fermion
system).

	NormalQuadraticHamiltonian(n, data)
"""
struct NormalQuadraticHamiltonian{T<:Number} <: QuadraticHamiltonian{T}
	data::Vector{AdagATerm{T}}
	n::Int
end
NormalQuadraticHamiltonian(::Type{T}, n::Int) where {T<:Number} = NormalQuadraticHamiltonian(Vector{AdagATerm{T}}(), n)
NormalQuadraticHamiltonian(n::Int, x::Vector{AdagATerm{T}}) where {T} = NormalQuadraticHamiltonian(x, n)

"""
	struct GenericQuadraticHamiltonian{T}

General quadratic Hamiltonian, possibly containing â†â, â†â† and ââ terms (e.g. BCS-type
Hamiltonians).

	GenericQuadraticHamiltonian(n, data)
"""
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

"""
	quadratichamiltonian(n, x)

Construct a Hamiltonian from a list of quadratic terms `x`: returns a
`NormalQuadraticHamiltonian` if all terms are `AdagATerm`, otherwise a
`GenericQuadraticHamiltonian`.
"""
function quadratichamiltonian(n::Int, x::Vector{<:QuadraticTerm{T}}) where {T<:Number}
	if all(y->isa(y, AdagATerm), x)
		return NormalQuadraticHamiltonian(n, convert(Vector{AdagATerm{T}}, x))
	else
		return GenericQuadraticHamiltonian(n, x)
	end
end

# coefficient matrix of Quadratic Hamiltonians
"""
	cmatrix(h)

Return the coefficient matrix of the quadratic Hamiltonian `h` (for BCS-type Hamiltonians
this is the 2L×2L Bogoliubov-de Gennes form).

The coefficient matrix is defined such that the quadratic Hamiltonian reads
h = [h₁₁ c₁†c₁, h₁₂ c₁†c₂, h₁₃ c₁†c₃...; h₂₁ c₂†c₁, h₂₂ c₂†c₂, h₂₃ c₂†c₃; ...],
and the quadratic observables are encoded in the coefficient density matrix (cdm) ρ.
The time evolution of a free-fermion system is then dρ/dt = -i [hᵀ, ρ].
"""
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
