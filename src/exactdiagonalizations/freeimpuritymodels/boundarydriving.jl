"""
	struct BoundaryDriving

Each band is coupled to a bath
"""
struct BoundaryDriving{P<:AbstractParticle, T<:Number, B<:AbstractDiscreteBath{P}} <: AbstractFreeImpurityModel{P}
	hsys::QuadraticHamiltonian{P, T}
	leftbath::B
	rightbath::B
end

num_sites(x::BoundaryDriving) = num_bands(x) + num_sites(x.leftbath) + num_sites(x.rightbath)

Base.eltype(::Type{BoundaryDriving{P, T, B}}) where {P, T, B} = T

function cmatrix(m::BoundaryDriving)
	hsys = cmatrix(m.hsys)
	L = size(hsys, 1)
	N = num_sites(m)
	h = zeros(eltype(m), N, N)
	h[1:L, 1:L] = hsys
	pos = L

	for (band, bj) in ((1, m.leftbath), (L, m.rightbath))
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

function separablestate(m::BoundaryDriving, ρ_sys::AbstractMatrix)
	L = num_bands(m)
	(size(ρ_sys, 1) == size(ρ_sys, 2) == L) || throw(DimensionMismatch("impurity density matrix size mismatch with impurity Hamiltonian size"))
	N = num_sites(m)
	ρ = zeros(eltype(m), N, N)
	ρ[1:L, 1:L] = ρ_sys

	pos = L
	for (band, bj) in ((1, m.leftbath), (L, m.rightbath))
		bj = m.baths[band]
		Lj = num_sites(bj)
		ρ[pos+1:pos+Lj, pos+1:pos+Lj] = thermalstate(bj)
		pos += Lj
	end
	return ρ
end

function thermalstate(m::OneBandOneBath)
	β, μ = m.leftbath.β, m.leftbath.μ
	((β==m.rightbath.β) && (μ==m.rightbath.μ)) || throw(ArgumentError("thermalstate requires all the baths to have the same β and μ"))
	h = cmatrix(m)
	evals, U = eigen(Hermitian(h))
	evals = [thermaloccupation(particletype(m), β, μ, item) for item in evals]
	return U * Diagonal(evals) * U'
end

function leftparticlecurrent_cmatrix(m::BoundaryDriving)
	N = num_sites(m)
	return _particlecurrent_util!(zeros(eltype(m), N, N), m.leftbath, leftbathsites(m), 1)
end 
function rightparticlecurrent_cmatrix(m::BoundaryDriving)
	N = num_sites(m)
	return _particlecurrent_util!(zeros(eltype(m), N, N), m.rightbath, rightbathsites(m), 1)
end

function leftbathsites(m::BoundaryDriving) 
	L = num_bands(m)
	return L+1:L+num_sites(m.leftbath)
end
function rightbathsites(m::BoundaryDriving)
	L = num_bands(m) + num_sites(m.leftbath)
	return L+1:L+num_sites(m.rightbath)
end

function _particlecurrent_util!(h::AbstractMatrix, b::AbstractDiscreteBath, bsites, band::Int)
	fs = spectrumvalues(b)
	for (j, v) in zip(bsites, fs)
		h[band, j] = v
		h[j, band] = v
	end
	return h
end