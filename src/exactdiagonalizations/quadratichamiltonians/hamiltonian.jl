abstract type AbstractQuadraticHamiltonian{P<:AbstractParticle} end
particletype(::Type{<:AbstractQuadraticHamiltonian{P}}) where {P} = P
particletype(x::AbstractQuadraticHamiltonian) = particletype(typeof(x))
Base.isempty(x::AbstractQuadraticHamiltonian) = isempty(x.data)
function Base.length(x::AbstractQuadraticHamiltonian)
	@assert !isempty(x)
	m = 0
	for item in x.data
		m = max(m, maximum(positions(item)))
	end
	return m
end



struct QuadraticHamiltonian{P, T<:Number} <: AbstractQuadraticHamiltonian{P}
	data::Vector{Tunneling{P, T}}
end

QuadraticHamiltonian{P, T}() where {P<:AbstractParticle, T<:Number} = QuadraticHamiltonian(Vector{Tunneling{P, T}}())

Base.push!(x::QuadraticHamiltonian{P}, m::Tunneling{P}) where {P} = push!(x.data, m)
Base.copy(x::QuadraticHamiltonian) = QuadraticHamiltonian(copy(x.data))

changepositions(x::QuadraticHamiltonian, m::AbstractDict{Int, Int}) = QuadraticHamiltonian([changepositions(item, m) for item in x.data])

"""
	cmatrix(x::QuadraticHamiltonian)

Return the coefficient matrix for free fermion model
"""
function cmatrix(x::QuadraticHamiltonian{P, T}, L::Int) where {T<:Number}
	@assert length(x) <= L
	m = zeros(T, L, L)
	for item in x.data
		i, j = positions(item)
		m[i, j] = coeff(item)
	end
	return m
end
cmatrix(x::QuadraticHamiltonian) = cmatrix(x, length(x))


const FermionicQuadraticHamiltonian{T} = QuadraticHamiltonian{Fermion, T} where {T<:Number}
const BosonicQuadraticHamiltonian{T} = QuadraticHamiltonian{Boson, T} where {T<:Number}