struct BoundaryDriving{B<:AbstractDiscreteBath, M<:AbstractMatrix}
	hsys::M
	leftbath::B
	rightbath::B
end

num_bands(m::BoundaryDriving) = size(m.hsys, 1)
num_sites(m::BoundaryDriving) = num_bands(m) + num_sites(m.leftbath) + num_sites(m.rightbath)
Base.eltype(::Type{BoundaryDriving{B, M}}) where {B, M} = eltype(M)
Base.eltype(x::BoundaryDriving) = eltype(typeof(x))
particletype(::Type{BoundaryDriving{B, M}}) where {B, M} = particletype(B)
particletype(x::BoundaryDriving) = particletype(typeof(x))

# hamiltonian and cmatrix
function freehamiltonian(h::BoundaryDriving)
	T = eltype(h)
	data = AdagATerm{T}[]

	hsys, leftbath, rightbath = h.hsys, h.leftbath, h.rightbath
	L = size(hsys, 1)
	for i in 1:L, j in 1:L
		c = hsys[i, j]
		if c != zero(c)
			push!(data, adaga(i, j, coeff=c))
		end
	end

	pos = L
	for (band, bj) in ((1, leftbath), (L, rightbath))
		Lj = num_sites(bj)
		ws = frequencies(bj)
		for j in 1:Lj
			push!(data, adaga(pos+j, pos+j, coeff=ws[j]))
		end
		pos += Lj
	end
	return NormalQuadraticHamiltonian(num_sites(h), data)	
end
function hamiltonian(h::BoundaryDriving)
	leftbath, rightbath = h.leftbath, h.rightbath
	L = num_bands(h)
	ham = freehamiltonian(h)

	pos = L
	for (band, bj) in ((1, leftbath), (L, rightbath))
		Lj = num_sites(bj)
		fs = spectrumvalues(bj)
		for j in 1:Lj
			t = adaga(band, pos + j, coeff=fs[j])
			push!(ham, t)
			push!(ham, t')
		end
		pos += Lj
	end
	return ham
end


function cmatrix(m::BoundaryDriving)
	hsys, leftbath, rightbath = m.hsys, m.leftbath, m.rightbath
	L = size(hsys, 1)
	N = num_sites(leftbath) + L + num_sites(rightbath)
	h = zeros(eltype(hsys), N, N)
	h[1:L, 1:L] = hsys

	pos = L
	for (band, bj) in ((1, leftbath), (L, rightbath))
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

# thermal state
function thermocdm(m::BoundaryDriving)
	β, μ = m.leftbath.β, m.leftbath.μ
	((β==m.rightbath.β) && (μ==m.rightbath.μ)) || throw(ArgumentError("thermocdm requires all the baths to have the same β and μ"))
	h = cmatrix(m)
	return thermocdm(particletype(m), eigencache(h), β=β, μ=μ)
end

# separable state
function separablecdm(m::BoundaryDriving, ρ_sys::AbstractMatrix)
	(size(ρ_sys) == size(m.hsys)) || throw(DimensionMismatch("Hamiltonian size mismatch with density matrix size"))
	leftbath, rightbath = m.leftbath, m.rightbath
	L = num_bands(m)
	N = num_sites(m)
	ρ = zeros(eltype(ρ_sys), N, N)
	ρ[1:L, 1:L] = ρ_sys

	pos = L
	for (band, bj) in ((1, leftbath), (L, rightbath))
		Lj = num_sites(bj)
		ρ[pos+1:pos+Lj, pos+1:pos+Lj] = thermocdm(bj)
		pos += Lj
	end
	return ρ
end



# currents
function leftparticlecurrent_cmatrix(m::BoundaryDriving)
	N = num_sites(m)
	return _particlecurrent_util!(zeros(ComplexF64, N, N), m.leftbath, leftbathsites(m), 1)
end 
function rightparticlecurrent_cmatrix(m::BoundaryDriving)
	N = num_sites(m)
	return _particlecurrent_util!(zeros(ComplexF64, N, N), m.rightbath, rightbathsites(m), num_bands(m))
end
function leftheatcurrent_cmatrix(m::BoundaryDriving)
	N = num_sites(m)
	return _heatcurrent_util!(zeros(ComplexF64, N, N), m.leftbath, leftbathsites(m), 1)
end 
function rightheatcurrent_cmatrix(m::BoundaryDriving)
	N = num_sites(m)
	return _heatcurrent_util!(zeros(ComplexF64, N, N), m.rightbath, rightbathsites(m), num_bands(m))
end

function leftbathsites(m::BoundaryDriving) 
	L = num_bands(m)
	return L+1:L+num_sites(m.leftbath)
end
function rightbathsites(m::BoundaryDriving)
	L = num_bands(m) + num_sites(m.leftbath)
	return L+1:L+num_sites(m.rightbath)
end
