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
	return thermocdm(cache, β=b.β, μ=b.μ)
end

function cmatrix(b::AbstractDiscreteBCSBath)
	L = num_sites(b)
	Δ = b.Δ
	h = zeros(typeof(Δ), 4L, 4L)
	for (i, ϵ) in enumerate(frequencies(b))
		h[2i-1, 2i-1] = ϵ
		h[2i, 2i] = ϵ
		
		h[2L+2i-1, 2L+2i-1] = 1-ϵ
		h[2L+2i, 2L+2i] = 1-ϵ	

		h[2i-1, 2L+2i] = -Δ
		h[2L+2i, 2i-1] = -conj(Δ)
	end	
	return h
end

include("boundarydriving.jl")