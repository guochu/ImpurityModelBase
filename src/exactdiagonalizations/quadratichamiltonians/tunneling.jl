abstract type AbstractTunnelingTerm{P<:AbstractParticle} end

particletype(::Type{<:AbstractTunnelingTerm{P}}) where {P} = P
particletype(x::AbstractTunnelingTerm) = particletype(typeof(x))

positions(x::AbstractTunnelingTerm) = x.positions
coeff(x::AbstractTunnelingTerm) = x.coeff

Base.:*(m::Number, s::AbstractTunnelingTerm) = s * m
Base.:/(s::AbstractTunnelingTerm, m::Number) = s * (1 / m)
Base.:+(s::AbstractTunnelingTerm) = s
Base.:-(s::AbstractTunnelingTerm) = (-1) * s


"""
	struct Tunneling{P, T<Number}
"""
struct Tunneling{P, T<Number} <: AbstractTunnelingTerm{P}
	positions::Tuple{Int, Int}
	coeff::T
end

tunneling(::Type{P}, pos::Tuple{Int, Int}; coeff::Number=1) where {P<:AbstractParticle} = Tunneling(P, pos, coeff)
fermionictunneling(pos::Tuple{Int, Int}; kwargs...) = tunneling(Fermion(), pos; kwargs...)
bosonictunneling(pos::Tuple{Int, Int}; kwargs...) = tunneling(Boson(), pos; kwargs...)

function Base.adjoint(x::Tunneling{P}) where {P}
	i, j = positions(x)
	return Tunneling{P}((j, i), conj(coeff(x)))
end

Base.copy(x::Tunneling{P}) where P = Tunneling{P}(positions(x), coeff(x))
Base.:*(s::Tunneling{P}, m::Number) where P = Tunneling{P}(positions(s), coeff=coeff(s) * m)

changepositions(x::Tunneling{P}, m::AbstractDict{Int, Int}) where P = Tunneling{P}(ntuple(k=>m[positions(x)[k]], Val(2)), coeff(x))
Base.convert(::Type{Tunneling{P, T}}, x::Tunneling{P}) where {P, T} = Tunneling{P}(positions(x), convert(T, coeff(x)))