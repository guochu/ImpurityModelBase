include("cache.jl")
include("bcscache.jl")
include("operators.jl")
include("correlations.jl")

include("greenfunctions.jl")
include("toulouse.jl")

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