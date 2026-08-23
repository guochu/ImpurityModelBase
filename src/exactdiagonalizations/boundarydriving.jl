# BoundaryDriving does not support BCS or BEC type bath
"""
	struct BoundaryDriving

Boundary-driven model: a system with Hamiltonian `hsys`, coupled on the left and right to
discrete baths `leftbath` and `rightbath`; used for non-equilibrium transport simulations
(BCS/BEC baths are not supported).
"""
struct BoundaryDriving{B<:AbstractDiscreteNormalBath, M<:AbstractMatrix}
	hsys::M
	leftbath::B
	rightbath::B
end

num_bands(m::BoundaryDriving) = size(m.hsys, 1)
num_sites(m::BoundaryDriving) = num_bands(m) + num_sites(m.leftbath) + num_sites(m.rightbath)
Base.eltype(::Type{BoundaryDriving{B, M}}) where {B, M} = promote_type(eltype(B), eltype(M)) 
Base.eltype(x::BoundaryDriving) = eltype(typeof(x))
particletype(::Type{BoundaryDriving{B, M}}) where {B, M} = particletype(B)
particletype(x::BoundaryDriving) = particletype(typeof(x))

# hamiltonian and cmatrix
function freehamiltonian(h::BoundaryDriving; include_chemical::Bool=false)
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
			wj = include_chemical ? ws[j] - bj.μ : ws[j]
			push!(data, adaga(pos+j, pos+j, coeff=wj))
		end
		pos += Lj
	end
	return NormalQuadraticHamiltonian(num_sites(h), data)	
end
"""
	hamiltonian(m::BoundaryDriving; include_chemical=false)

Construct the full quadratic Hamiltonian (system + left/right baths + couplings) of the
boundary-driven model; `include_chemical=true` absorbs the chemical potential into the bath
energies.
"""
function hamiltonian(h::BoundaryDriving; include_chemical::Bool=false)
	leftbath, rightbath = h.leftbath, h.rightbath
	L = num_bands(h)
	ham = freehamiltonian(h, include_chemical=include_chemical)

	pos = L
	for (band, bj) in ((1, leftbath), (L, rightbath))
		Lj = num_sites(bj)
		fs = spectrumcouplings(bj)
		for j in 1:Lj
			t = adaga(band, pos + j, coeff=fs[j])
			push!(ham, t)
			push!(ham, t')
		end
		pos += Lj
	end
	return ham
end


"""
	cmatrix(m::BoundaryDriving)

Return the coefficient matrix of the boundary-driven model (system + left/right baths +
couplings).
"""
function cmatrix(m::BoundaryDriving)
	hsys, leftbath, rightbath = m.hsys, m.leftbath, m.rightbath
	L = size(hsys, 1)
	N = num_sites(leftbath) + L + num_sites(rightbath)
	h = zeros(eltype(hsys), N, N)
	h[1:L, 1:L] = hsys

	pos = L
	for (band, bj) in ((1, leftbath), (L, rightbath))
		Lj = num_sites(bj)
		ws, fs = frequencies(bj), spectrumcouplings(bj)
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
"""
	thermocdm(m::BoundaryDriving)

Construct the thermal equilibrium coefficient density matrix (cdm) of the boundary-driven
model; the cdm is a single-particle correlation matrix of the quadratic Hamiltonian, from
which all quadratic observables can be evaluated. The left and right baths must have the
same `β` and `μ`.
"""
function thermocdm(m::BoundaryDriving)
	β, μ = m.leftbath.β, m.leftbath.μ
	((β==m.rightbath.β) && (μ==m.rightbath.μ)) || throw(ArgumentError("thermocdm requires all the baths to have the same β and μ"))
	h = cmatrix(m)
	return thermocdm(particletype(m), eigencache(h), β=β, μ=μ)
end

"""
	fermionicthermodm(m::BoundaryDriving)

Construct the thermal equilibrium (true) density matrix (dm) ρ = exp(-β(Ĥ-μN̂))/Z of the
boundary-driven model in the many-body Fock space. Note that this is the genuine density
matrix, not the coefficient density matrix (cdm, see `thermocdm`). The left and right
baths must have the same `β` and `μ`.
"""
function fermionicthermodm(m::BoundaryDriving)
	(particletype(m) == Fermion) || throw(ArgumentError("Fermion particletype assumed"))
	β, μ = m.leftbath.β, m.leftbath.μ
	((β==m.rightbath.β) && (μ==m.rightbath.μ)) || throw(ArgumentError("thermodm requires all the baths to have the same β and μ"))
	h = hamiltonian(m, include_chemical=true)
	return fermionicthermodm(h, β=β)
end

# separable state
"""
	separablecdm(m::BoundaryDriving, ρ_sys)

Construct a separable coefficient density matrix (cdm) with the system initialized in
`ρ_sys` and the left/right baths in their own thermal equilibrium.
"""
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

"""
	fermionicseparabledm(m::BoundaryDriving, sysdm)

Construct the separable (true) fermionic density matrix (dm) with the system initialized
in `sysdm` and the left/right baths in their own thermal equilibrium. Note that this is
the genuine density matrix in the many-body Fock space, not the coefficient density
matrix (cdm, see `separablecdm`).
"""
function fermionicseparabledm(m::BoundaryDriving, sysdm::AbstractMatrix)
	(particletype(m) == Fermion) || throw(ArgumentError("Fermion particletype assumed"))
	(size(sysdm, 1) == 2^(size(m.hsys, 1))) || throw(DimensionMismatch("Hamiltonian size mismatch with density operator size"))
	leftbath, rightbath = m.leftbath, m.rightbath
	ρ_l = fermionicthermodm(hamiltonian(leftbath, include_chemical=true), β=leftbath.β)
	ρ_r = fermionicthermodm(hamiltonian(rightbath, include_chemical=true), β=rightbath.β)
	return kron(sysdm, ρ_l, ρ_r)
end

# currents
"""
	leftparticlecurrent_cmatrix(m::BoundaryDriving)
	rightparticlecurrent_cmatrix(m::BoundaryDriving)

Return the coefficient-matrix representation of the particle current operators at the
left/right contacts (used for averages with reduced density matrices).
"""
function leftparticlecurrent_cmatrix(m::BoundaryDriving)
	N = num_sites(m)
	return _particlecurrent_util!(zeros(ComplexF64, N, N), m.leftbath, leftbathsites(m), 1)
end 
function rightparticlecurrent_cmatrix(m::BoundaryDriving)
	N = num_sites(m)
	return _particlecurrent_util!(zeros(ComplexF64, N, N), m.rightbath, rightbathsites(m), num_bands(m))
end
"""
	leftheatcurrent_cmatrix(m::BoundaryDriving)
	rightheatcurrent_cmatrix(m::BoundaryDriving)

Return the coefficient-matrix representation of the heat current operators at the left/right
contacts (used for averages with reduced density matrices).
"""
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


"""
	leftparticlecurrent_hamiltonian(m::BoundaryDriving)
	rightparticlecurrent_hamiltonian(m::BoundaryDriving)

Return the Hamiltonian representation of the particle current operators at the left/right
contacts (used for averages with full Fock-space density matrices).
"""
function leftparticlecurrent_hamiltonian(m::BoundaryDriving)
	h = NormalQuadraticHamiltonian(ComplexF64, num_sites(m))
	return _particlecurrent_hamiltonian_util!(h, m.leftbath, leftbathsites(m), 1)
end 
function rightparticlecurrent_hamiltonian(m::BoundaryDriving)
	h = NormalQuadraticHamiltonian(ComplexF64, num_sites(m))
	return _particlecurrent_hamiltonian_util!(h, m.rightbath, rightbathsites(m), num_bands(m))
end
"""
	leftheatcurrent_hamiltonian(m::BoundaryDriving)
	rightheatcurrent_hamiltonian(m::BoundaryDriving)

Return the Hamiltonian representation of the heat current operators at the left/right
contacts (used for averages with full Fock-space density matrices).
"""
function leftheatcurrent_hamiltonian(m::BoundaryDriving)
	h = NormalQuadraticHamiltonian(ComplexF64, num_sites(m))
	return _heatcurrent_hamiltonian_util!(h, m.leftbath, leftbathsites(m), 1)
end 
function rightheatcurrent_hamiltonian(m::BoundaryDriving)
	h = NormalQuadraticHamiltonian(ComplexF64, num_sites(m))
	return _heatcurrent_hamiltonian_util!(h, m.rightbath, rightbathsites(m), num_bands(m))
end