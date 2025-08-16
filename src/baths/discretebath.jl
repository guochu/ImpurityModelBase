abstract type AbstractDiscreteBath{P<:AbstractParticle} end
particletype(::Type{<:AbstractDiscreteBath{P}}) where {P<:AbstractParticle} = P
particletype(x::AbstractDiscreteBath) = particletype(typeof(x))

frequencies(b::AbstractDiscreteBath) = b.ws
spectrumvalues(b::AbstractDiscreteBath) = b.fs
num_sites(x::AbstractDiscreteBath) = length(frequencies(x))

"""
	struct DiscreteBath{P<:AbstractParticle}

Fermionic bath container, includes a bath spectrum density,
the inverse temperature β and the chemical potential μ
"""
struct DiscreteBath{P<:AbstractParticle} <: AbstractDiscreteBath{P}
	ws::Vector{Float64}
	fs::Vector{Float64}
	β::Float64
	μ::Float64
end
function DiscreteBath(::Type{P}, ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; β::Real, μ::Real=0) where {P<:AbstractParticle}
	(length(ws) == length(fs)) || throw(DimensionMismatch("num frequencies mismatch with num spectrum values"))
	all(x->x>=0, fs) || throw(ArgumentError("spectrum values can not be negative"))
	issorted(ws) || throw("frequencies should be sorted")
	DiscreteBath{P}(convert(Vector{Float64}, ws), convert(Vector{Float64}, fs), float(β), float(μ))
end 

DiscreteBosonicBath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBath(Boson, ws, fs; kwargs...)
"""
	discretebosonicbath(ws, fs; β, μ) 

Return a bosonic bath with β and μ
"""
discretebosonicbath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBath(Boson, ws, fs; kwargs...)

DiscreteFermionicBath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBath(Fermion, ws, fs; kwargs...)
"""
	discretefermionicbath(ws, fs; β, μ) 

Return a fermionic bath with β and μ
"""
discretefermionicbath(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBath(Fermion, ws, fs; kwargs...)

"""
	struct DiscreteVacuum{P<:AbstractParticle}

Fermionic bath container, includes a bath spectrum density,
the chemical potential μ
the inverse temperature β=Inf
"""
struct DiscreteVacuum{P<:AbstractParticle} <: AbstractDiscreteBath{P}
	ws::Vector{Float64}
	fs::Vector{Float64}
	μ::Float64	
end
function DiscreteVacuum(::Type{P}, ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; μ::Real=0) where {P<:AbstractParticle}
	(length(ws) == length(fs)) || throw(DimensionMismatch("num frequencies mismatch with num spectrum values"))
	all(x->x>=0, fs) || throw(ArgumentError("spectrum values can not be negative"))
	issorted(ws) || throw("frequencies should be sorted")
	DiscreteVacuum{P}(convert(Vector{Float64}, ws), convert(Vector{Float64}, fs), float(μ))
end 

DiscreteBosonicVacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteVacuum(Boson, ws, fs; kwargs...)
discretebosonicvacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBosonicVacuum(Boson, ws, fs; kwargs...)
DiscreteFermionicVacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteVacuum(Fermion, ws, fs; kwargs...)
discretefermionicvacuum(ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteFermionicVacuum(Fermion, ws, fs; kwargs...)


const AbstractDiscreteBosonicBath = Union{DiscreteBath{Boson}, DiscreteVacuum{Boson}} 
const AbstractDiscreteFermionicBath = Union{DiscreteBath{Fermion}, DiscreteVacuum{Fermion}}

thermaloccupation(bath::AbstractDiscreteBath, ϵ::Real) = thermaloccupation(particletype(bath), bath.β, bath.μ, ϵ)


function Base.getproperty(m::DiscreteBath, s::Symbol)
	if s == :T
		return 1 / m.β
	else
		return getfield(m, s)
	end
end

function Base.getproperty(m::DiscreteVacuum, s::Symbol)
	if s == :β
		return Inf
	elseif s == :T
		return 0.
	else
		return getfield(m, s)
	end
end

discretebath(::Type{Boson}, ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBath(Boson, ws, fs; kwargs...)
discretebath(::Type{Fermion}, ws::AbstractVector{<:Real}, fs::AbstractVector{<:Real}; kwargs...) = DiscreteBath(Fermion, ws, fs; kwargs...)
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
			push!(couplings, sqrt(v1))
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