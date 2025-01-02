abstract type AbstractFreeImpurityModel{P<:AbstractParticle} end
particletype(::Type{<:AbstractFreeImpurityModel{P}}) where {P} = P
particletype(x::AbstractFreeImpurityModel) = particletype(typeof(x))

num_bands(x::AbstractFreeImpurityModel) = length(x.hsys)
Base.eltype(x::AbstractFreeImpurityModel) = eltype(typeof(x))
syssites(x::AbstractFreeImpurityModel) = 1:num_bands(x)

function thermalstate(b::AbstractDiscreteBath)
	L = num_sites(b)
	ρ = zeros(Float64, L, L)
	for (i, ϵ) in enumerate(frequencies(b))
		ρ[i, i] = thermaloccupation(particletype(b), b.β, b.μ, ϵ)
	end
	return ρ
end


include("onebandonebath.jl")
include("boundarydriving.jl")


# """
# 	struct FreeImpurityModel
# """
# struct FreeImpurityModel{P<:AbstractParticle, T<:Number, B<:AbstractDiscreteBath{P}}
# 	hsys::QuadraticHamiltonian{P, T}
# 	baths::Vector{Union{Vector{B}, Missing}}
# end

# function FreeImpurityModel(Hsys::QuadraticHamiltonian{P}, baths::Dict{Int, Vector{B}}) where {P<:AbstractParticle, B<:AbstractDiscreteBath{P}}
# 	L = length(Hsys)
# 	all(x -> 1<=x<=L, keys(baths)) || throw(BoundsError(1:L, keys(baths)))
# 	baths_new = Vector{Union{Vector{B}, Missing}}(missing, L)
# 	for band in 1:L
# 		b = get(baths, band, nothing)
# 		if !isnothing(b)
# 			baths_new[band] = b
# 		end
# 	end
# 	FreeImpurityModel(Hsys, baths_new)
# end

# function num_sites(x::FreeImpurityModel)
# 	L = length(x.hsys)
# 	for b in values(x.baths)
# 		for bj in b
# 			L += num_sites(bj)
# 		end
# 	end
# 	return L
# end


# function cmatrix(m::FreeImpurityModel)
# 	hsys = cmatrix(m.hsys)
# 	L = size(Hsys, 1)
# 	N = num_sites(m)
# 	h = zeros(eltype(hsys), N, N)
# 	h[1:L, 1:L] = hsys
# 	pos = L
# 	for band in 1:L
# 		b = m.baths[band]
# 		if !ismissing(b)
# 			for bj in b
# 				Lj = num_sites(bj)
# 				ws, fs = frequencies(bj), spectrumvalues(bj)
# 				for j in 1:Lj
# 					h[pos + j, pos + j] = ws[j]
# 					h[band, pos + j] = fs[j]
# 					h[pos + j, band] = fs[j]
# 				end
# 				pos += Lj
# 			end
# 		end
# 	end
# 	return h
# end	
