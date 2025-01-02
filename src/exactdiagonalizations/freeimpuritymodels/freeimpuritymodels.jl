"""
	struct FreeImpurityModel
"""
struct FreeImpurityModel{H<:QuadraticHamiltonian, B<:AbstractBath}
	hsys::H
	baths::Vector{Union{Vector{B}, Missing}}
end

function FreeImpurityModel(Hsys::QuadraticHamiltonian, baths::Dict{Int, Vector{B}}) where {B<:AbstractBath}
	(particletype(Hsys) === particletype(B)) || throw(ArgumentError("particletype mismatch"))
	L = length(Hsys)
	all(x -> 1<=x<=L, keys(baths)) || throw(BoundsError(1:L, keys(baths)))
	baths_new = Vector{Union{Vector{B}, Missing}}(missing, L)
	for band in 1:L
		b = get(baths, band, nothing)
		if !isnothing(b)
			baths_new[band] = b
		end
	end
	FreeImpurityModel(Hsys, baths_new)
end

function num_sites(x::FreeImpurityModel)
	L = length(x.hsys)
	for b in values(x.baths)
		for bj in b
			L += num_sites(bj)
		end
	end
	return L
end


function cmatrix(m::FreeImpurityModel)
	hsys = cmatrix(m.hsys)
	L = size(Hsys, 1)
	N = num_sites(m)
	h = zeros(eltype(hsys), N, N)
	h[1:L, 1:L] = hsys
	pos = L
	for band in 1:L
		b = m.baths[band]
		if !ismissing(b)
			for bj in b
				Lj = num_sites(bj)
			end
		end
	end

end
