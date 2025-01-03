include("cache.jl")
include("greenfunctions.jl")
include("toulouse.jl")

function thermalstate(b::AbstractDiscreteBath)
	L = num_sites(b)
	ρ = zeros(Float64, L, L)
	for (i, ϵ) in enumerate(frequencies(b))
		ρ[i, i] = thermaloccupation(b, ϵ)
	end
	return ρ
end


include("boundarydriving.jl")