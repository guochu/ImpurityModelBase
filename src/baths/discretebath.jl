# abstract type AbstractDiscreteBath{P<:AbstractParticle} end
# particletype(::Type{<:AbstractDiscreteBath{P}}) where {P<:AbstractParticle} = P
# particletype(x::AbstractDiscreteBath) = particletype(typeof(x))



# """
# 	struct DiscreteBath{P<:AbstractParticle}

# Fermionic bath container, includes a bath spectrum density,
# the inverse temperature β and the chemical potential μ
# """
# struct DiscreteBath{P<:AbstractParticle} <: AbstractDiscreteNormalBath{P}
# 	ws::Vector{Float64}
# 	fs::Vector{Float64}
# 	β::Float64
# 	μ::Float64
# end
# function DiscreteBath(::Type{P}, ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; β::Real, μ::Real=0) where {P<:AbstractParticle}
# 	(length(ws) == length(fs)) || throw(DimensionMismatch("num frequencies mismatch with num spectrum values"))
# 	all(x->x>=0, fs) || throw(ArgumentError("spectrum values can not be negative"))
# 	issorted(ws) || throw("frequencies should be sorted")
# 	DiscreteBath{P}(convert(Vector{Float64}, ws), convert(Vector{Float64}, fs), float(β), float(μ))
# end 
# Base.eltype(::Type{DiscreteBath{P}}) where {P} = Float64

"""
	DiscreteBath{P}

Type alias for a discrete particle bath, equivalent to `Bath{P, DiscreteSpectrum}`, whose
spectral density is given by discrete frequencies `ws` and spectral values `fs`.

	DiscreteBath(::Type{P}, ws, fs; β, μ=0)
"""
const DiscreteBath{P<:AbstractParticle} = Bath{P, DiscreteSpectrum}

DiscreteBath(::Type{P}, f::DiscreteSpectrum; kwargs...) where {P<:AbstractParticle} = Bath(P, f; kwargs...)
DiscreteBath(::Type{P}, ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) where {P<:AbstractParticle} = DiscreteBath(P, DiscreteSpectrum(ws, fs); kwargs...)

DiscreteBosonicBath(f::DiscreteSpectrum; kwargs...) = DiscreteBath(Boson, f; kwargs...)
DiscreteBosonicBath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBath(Boson, ws, fs; kwargs...)
"""
	discretebosonicbath(ws, fs; β, μ) 

Return a bosonic bath with β and μ
"""
discretebosonicbath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBath(Boson, ws, fs; kwargs...)

DiscreteFermionicBath(f::DiscreteSpectrum; kwargs...) = DiscreteBath(Fermion, f; kwargs...)
DiscreteFermionicBath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBath(Fermion, ws, fs; kwargs...)
"""
	discretefermionicbath(ws, fs; β, μ) 

Return a fermionic bath with β and μ
"""
discretefermionicbath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBath(Fermion, ws, fs; kwargs...)

# """
# 	struct DiscreteVacuum{P<:AbstractParticle}

# Fermionic bath container, includes a bath spectrum density,
# the chemical potential μ
# the inverse temperature β=Inf
# """
# struct DiscreteVacuum{P<:AbstractParticle} <: AbstractDiscreteNormalBath{P}
# 	ws::Vector{Float64}
# 	fs::Vector{Float64}
# 	μ::Float64	
# end
# function DiscreteVacuum(::Type{P}, ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; μ::Real=0) where {P<:AbstractParticle}
# 	(length(ws) == length(fs)) || throw(DimensionMismatch("num frequencies mismatch with num spectrum values"))
# 	all(x->x>=0, fs) || throw(ArgumentError("spectrum values can not be negative"))
# 	issorted(ws) || throw("frequencies should be sorted")
# 	DiscreteVacuum{P}(convert(Vector{Float64}, ws), convert(Vector{Float64}, fs), float(μ))
# end 
# Base.eltype(::Type{DiscreteVacuum{P}}) where {P} = Float64

"""
	DiscreteVacuum{P}

Type alias for a discrete vacuum (zero-temperature) particle bath, equivalent to
`Vacuum{P, DiscreteSpectrum}`.

	DiscreteVacuum(::Type{P}, ws, fs; μ=0)
"""
const DiscreteVacuum{P<:AbstractParticle} = Vacuum{P, DiscreteSpectrum} 


DiscreteVacuum(::Type{P}, f::DiscreteSpectrum; kwargs...) where {P<:AbstractParticle} = Vacuum(P, f; kwargs...)
DiscreteVacuum(::Type{P}, ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) where {P<:AbstractParticle} = DiscreteVacuum(P, DiscreteSpectrum(ws, fs); kwargs...)

DiscreteBosonicVacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteVacuum(Boson, ws, fs; kwargs...)
"""
	discretebosonicvacuum(ws, fs; μ)

Construct a bosonic discrete vacuum bath (zero temperature) with frequencies `ws` and
spectral values `fs`.
"""
discretebosonicvacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBosonicVacuum(ws, fs; kwargs...)
DiscreteFermionicVacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteVacuum(Fermion, ws, fs; kwargs...)
"""
	discretefermionicvacuum(ws, fs; μ)

Construct a fermionic discrete vacuum bath (zero temperature) with frequencies `ws` and
spectral values `fs`.
"""
discretefermionicvacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteFermionicVacuum(ws, fs; kwargs...)


const AbstractDiscreteNormalBath{P<:AbstractParticle} = Union{DiscreteBath{P}, DiscreteVacuum{P}}

"""
	frequencies(b::AbstractDiscreteNormalBath)

Return the list of discrete frequencies of the discrete bath `b`.
"""
frequencies(b::AbstractDiscreteNormalBath) = frequencies(b.f)
"""
	spectrumvalues(b::AbstractDiscreteNormalBath)

Return the list of spectral values of the discrete bath `b` at each frequency.
"""
spectrumvalues(b::AbstractDiscreteNormalBath) = spectrumvalues(b.f)
"""
	spectrumcouplings(b::AbstractDiscreteNormalBath)

Return the coupling strengths of the discrete bath `b`, given by the square roots of the
spectral values.
"""
spectrumcouplings(b::AbstractDiscreteNormalBath) = spectrumcouplings(b.f)
"""
	num_sites(x)

Return the number of sites (modes) of a discrete bath, Hamiltonian or Toulouse model.
"""
num_sites(x::AbstractDiscreteNormalBath) = length(frequencies(x))

# const AbstractDiscreteBosonicBath = Union{DiscreteBath{Boson}, DiscreteVacuum{Boson}} 
# const AbstractDiscreteFermionicBath = Union{DiscreteBath{Fermion}, DiscreteVacuum{Fermion}}

# thermaloccupation(bath::AbstractDiscreteBath, ϵ::Real) = thermaloccupation(particletype(bath), bath.β, bath.μ, ϵ)
# num_sites(b::Union{DiscreteBath, DiscreteVacuum}) = length(frequencies(b))

# function Base.getproperty(m::DiscreteBath, s::Symbol)
# 	if s == :T
# 		return 1 / m.β
# 	else
# 		return getfield(m, s)
# 	end
# end

# function Base.getproperty(m::DiscreteVacuum, s::Symbol)
# 	if s == :β
# 		return Inf
# 	elseif s == :T
# 		return 0.
# 	else
# 		return getfield(m, s)
# 	end
# end

"""
	discretebath(::Type{P}, ws, fs; β, μ=0)
	discretebath(::Type{P}, freqs, f; β, μ=0, atol=1e-6)
	discretebath(b::AbstractBath; δw=0.1)

Construct a discrete bath of particle species `P`:
- directly from frequencies `ws` and spectral values `fs`;
- by discretizing a continuous spectral function `f` on the intervals defined by `freqs`
  (see [`spectrum_couplings`](@ref));
- or by uniformly discretizing a continuous bath `b` with step size `δw`.
"""
discretebath(::Type{Boson}, ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBath(Boson, ws, fs; kwargs...)
discretebath(::Type{Fermion}, ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBath(Fermion, ws, fs; kwargs...)
"""
	discretevacuum(::Type{P}, ws, fs; μ=0)
	discretevacuum(::Type{P}, freqs, f; μ=0, atol=1e-6)
	discretevacuum(b::AbstractBath; δw=0.1)

Construct a discrete vacuum (zero-temperature) bath of particle species `P`; the
parameters have the same meaning as in [`discretebath`](@ref).
"""
discretevacuum(::Type{Boson}, ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteVacuum(Boson, ws, fs; kwargs...)
discretevacuum(::Type{Fermion}, ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteVacuum(Fermion, ws, fs; kwargs...)

"""
	spectrum_couplings(freqs::Vector{<:Real}, f; kwargs...)

Arguments
* f: pectrum function
* freqs: L frequences {w1, w2, ... , wL}, sorted from small to large

Return discretized frequencies Xn and hybridizations Vn of size at most L-1,
where |Vn|^2 = ∫_{w_{n}}^{w_{n+1}} dx f^2(x), and 
Xn = ∫_{w_{n}}^{w_{n+1}} dx x f^2(x) / |Vn|^2.
If |Vn| = 0 for some n, the corresponding Vn and Xn are removed.


Reference "How to discretize a quantum bath for real-time evolution".
Eqs.(10a, 10b) 
"""
function spectrum_couplings(freqs::Union{Vector{<:Real}, AbstractRange}, f::Function; atol::Real=1.0e-6, kwargs...)
	@assert atol > 0
	L = length(freqs)
	omegas = Float64[]
	couplings = Float64[]
	for i in 1:L-1
		v1, err1 = quadgk(f, freqs[i], freqs[i+1]; kwargs...)
		v2, err2 = quadgk(x -> x*f(x), freqs[i], freqs[i+1]; kwargs...)
		if v1 <= zero(v1)
			println("couplings is nonpositive for frequency interval $((round(freqs[i], digits=6), round(freqs[i+1], digits=6)))")
			v1 = zero(v1)
		end
		if v1 > atol
			omega = v2 / v1
			@assert (freqs[i] <= omega <= freqs[i+1]) 
			push!(omegas, omega)
			# push!(couplings, sqrt(v1))
			push!(couplings, v1)
		end
	end
	return omegas, couplings
end

function discretebath(::Type{P}, freqs::Union{Vector{<:Real}, AbstractRange}, f::Function; atol::Real=1.0e-6, β::Real, μ::Real=0, kwargs...) where {P<:AbstractParticle}
	omegas, couplings = spectrum_couplings(freqs, f; atol=atol, kwargs...)
	return discretebath(P, omegas, couplings; β=β, μ=μ)
end
function discretebath(b::AbstractBath; δw::Real=0.1, kwargs...)
	f = b.spectrum
	freqs = lowerbound(f):δw:upperbound(f)
	return discretebath(particletype(b), freqs, f; β=b.β, μ=b.μ, kwargs...)
end
function discretevacuum(::Type{P}, freqs::Union{Vector{<:Real}, AbstractRange}, f::Function; atol::Real=1.0e-6, μ::Real=0, kwargs...) where {P<:AbstractParticle}
	omegas, couplings = spectrum_couplings(freqs, f; atol=atol, kwargs...)
	return discretevacuum(P, omegas, couplings; μ=μ) 
end
function discretevacuum(b::AbstractBath; δw::Real=0.1, kwargs...)
	f = b.spectrum
	freqs = lowerbound(f):δw:upperbound(f)
	return discretevacuum(particletype(b), freqs, f; μ=b.μ, kwargs...)
end