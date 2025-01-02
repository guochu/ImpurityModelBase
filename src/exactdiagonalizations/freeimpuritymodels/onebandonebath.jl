"""
	struct OneBandOneBath

Each band is coupled to a bath
"""
struct OneBandOneBath{P<:AbstractParticle, T<:Number, B<:AbstractDiscreteBath{P}} <: AbstractFreeImpurityModel{P}
	baths::Vector{B}
	hsys::QuadraticHamiltonian{P, T}
end

function OneBandOneBath(hsys::QuadraticHamiltonian{P}, baths::Vector{B}) where {P<:AbstractParticle, B<:AbstractDiscreteBath{P}}
	L = length(hsys)
	(L == length(baths)) || throw(DimensionMismatch("hamiltonian size mismatch with num baths"))
	return OneBandOneBath(hsys, baths)
end

Base.eltype(::Type{OneBandOneBath{P, T, B}}) where {P, T, B} = T

function num_sites(x::OneBandOneBath)
	L = num_bands(x)
	for b in x.baths
		L += num_sites(b)
	end
	return L
end

function cmatrix(m::OneBandOneBath)
	hsys = cmatrix(m.hsys)
	L = size(hsys, 1)
	N = num_sites(m)
	h = zeros(eltype(hsys), N, N)
	h[1:L, 1:L] = hsys
	pos = L
	for band in 1:L
		bj = m.baths[band]
		Lj = num_sites(bj)
		ws, fs = frequencies(bj), spectrumvalues(bj)
		for j in 1:Lj
			h[pos + j, pos + j] = ws[j]
			h[band, pos + j] = fs[j]
			h[pos + j, band] = fs[j]
		end
		pos += Lj
	end
	return h
end	

function separablestate(m::OneBandOneBath, ρ_sys::AbstractMatrix)
	L = num_bands(m)
	(size(ρ_sys, 1) == size(ρ_sys, 2) == L) || throw(DimensionMismatch("impurity density matrix size mismatch with impurity Hamiltonian size"))
	N = num_sites(m)
	ρ = zeros(eltype(m), N, N)
	ρ[1:L, 1:L] = ρ_sys

	pos = L
	for band in 1:L
		bj = m.baths[band]
		Lj = num_sites(bj)
		ρ[pos+1:pos+Lj, pos+1:pos+Lj] = thermalstate(bj)
		pos += Lj
	end
end

function thermalstate(m::OneBandOneBath)
	β, μ = m.baths[1].β, m.baths[1].μ
	all(x->(x.β==β) && (x.μ==μ), m.baths) || throw(ArgumentError("thermalstate requires all the baths to have the same β and μ"))
	h = cmatrix(m)
	evals, U = eigen(Hermitian(h))
	evals = [thermaloccupation(particletype(m), β, μ, item) for item in evals]
	return U * Diagonal(evals) * U'
end
