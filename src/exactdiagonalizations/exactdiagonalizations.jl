include("cache.jl")
include("bcscache.jl")
include("operators.jl")
include("correlations.jl")

include("greenfunctions.jl")
include("toulouse.jl")


function hamiltonian(b::AbstractDiscreteBath; include_chemical=false)
	T = eltype(b)
	data = AdagATerm{T}[]
	for (i, ϵ) in enumerate(frequencies(b))
		if include_chemical
			ϵ = ϵ - b.μ
		end
		push!(data, adaga(i, j, coeff=ϵ))
	end
	return NormalQuadraticHamiltonian(num_sites(b), data)	
end

function hamiltonian(b::AbstractDiscreteBCSBath)
	T = eltype(b)
	data = QuadraticTerm{T}[]
	Δ = b.Δ
	for (i, ϵ) in enumerate(frequencies(b))
		push!(data, adaga(2i-1, 2i-1, coeff=ϵ))
		push!(data, adaga(2i, 2i, coeff=ϵ))

		t = adagadag(2i-1, 2i, coeff=-Δ)
		push!(data, t)
		push!(data, t')
	end	
	return GenericQuadraticHamiltonian(num_sites(b), data)	
end

function thermocdm(b::AbstractDiscreteBath)
	L = num_sites(b)
	ρ = zeros(Float64, L, L)
	for (i, ϵ) in enumerate(frequencies(b))
		ρ[i, i] = thermaloccupation(b, ϵ)
	end
	return ρ
end

function cmatrix(b::AbstractDiscreteBath)
	L = num_sites(b)
	h = zeros(Float64, L, L)
	for (i, ϵ) in enumerate(frequencies(b))
		h[i, i] = ϵ
	end
	return h
end

function thermocdm(b::AbstractDiscreteBCSBath)
	h = cmatrix(b)
	cache = eigencache(h)
	return fermionicthermocdm(cache, β=b.β, μ=b.μ)
end

function cmatrix(b::AbstractDiscreteBCSBath)
	L = div(num_sites(b), 2)
	Δ = b.Δ
	h = zeros(typeof(Δ), 2L, 2L)
	g = zero(h)
	for (i, ϵ) in enumerate(frequencies(b))
		h[2i-1, 2i-1] = ϵ
		h[2i, 2i] = ϵ

		g[2i-1, 2i] = -Δ
	end	

	return bcs_cmatrix(h, g)
end

include("boundarydriving.jl")