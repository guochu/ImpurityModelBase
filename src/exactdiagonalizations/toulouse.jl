struct Toulouse{B<:AbstractDiscreteBath{Fermion}}
	bath::B
	ϵ_d::Float64
end
Toulouse(b::AbstractDiscreteNormalBath{Fermion}; ϵ_d::Real) = Toulouse(b, float(ϵ_d))
function Toulouse(b::AbstractDiscreteBCSBath; ϵ_d::Real)
	(b.μ == 0) || throw(ArgumentError("BCS bath should have μ=0"))
	return Toulouse(b, float(ϵ_d))
end
Base.eltype(::Type{Toulouse{B}}) where B = eltype(B)
Base.eltype(x::Toulouse) = eltype(typeof(x))

const NormalToulouse = Toulouse{T} where {T<:AbstractDiscreteNormalBath{Fermion}}
const BCSToulouse =Toulouse{T} where {T<:AbstractDiscreteBCSBath}

num_sites(m::NormalToulouse) = num_sites(m.bath) + 1
num_sites(m::BCSToulouse) = num_sites(m.bath) + 2


# hamiltonian and cmatrix
function freehamiltonian(h::NormalToulouse; include_chemical::Bool=false)
	T = eltype(h)
	data = AdagATerm{T}[]
	L = num_sites(h)
	ϵ_d = h.ϵ_d
	push!(data, adaga(1, 1, coeff=ϵ_d))
	bath = h.bath
	ws = frequencies(bath)
	for i in 1:L-1
		wj = include_chemical ? ws[i] - bath.μ : ws[i]
		push!(data, adaga(1+i, 1+i, coeff=wj))
	end
	return NormalQuadraticHamiltonian(L, data)
	
end
function hamiltonian(h::NormalToulouse; include_chemical::Bool=false)
	ham = freehamiltonian(h, include_chemical=include_chemical)
	fs = spectrumvalues(h.bath)
	for i in 1:num_sites(h)-1
		t = adaga(1, 1+i, coeff=fs[i])
		push!(ham, t)
		push!(ham, t')
	end
	return ham
end
function freehamiltonian(h::BCSToulouse)
	T = eltype(h)
	L = div(num_sites(h), 2)
	ham = GenericQuadraticHamiltonian(T, num_sites(h))
	bath = h.bath
	ws = frequencies(bath)
	Δ = bath.Δ
	ϵ_d = h.ϵ_d

	push!(ham, adaga(1,1, coeff=ϵ_d))
	push!(ham, adaga(2,2, coeff=ϵ_d))

	for i in 1:L-1
		push!(ham, adaga(2i+1,2i+1, coeff=ws[i]))
		push!(ham, adaga(2i+2, 2i+2, coeff=ws[i]))

		t = adagadag(2i+1, 2i+2, coeff=-Δ)
		push!(ham, t)
		push!(ham, t')
	end
	return ham
end
function hamiltonian(h::BCSToulouse)
	ham = freehamiltonian(h)
	fs =spectrumvalues(h.bath)
	L = div(num_sites(h), 2)
	for i in 1:L-1
		t = adaga(1, 2i+1, coeff=fs[i])
		push!(ham, t)
		push!(ham, t')
		t = adaga(2, 2i+2, coeff=fs[i])
		push!(ham, t)
		push!(ham, t')
	end
	return ham
end
# function hamiltonian(h::BCSToulouse)
# 	T = eltype(h)
# 	data = QuadraticTerm{T}[]

# 	L = num_sites(h)
# 	bath = h.bath
# 	ws, fs = frequencies(bath), spectrumvalues(bath)	
# 	Δ = bath.Δ
# 	ϵ_d = h.ϵ_d

# 	push!(data, adaga(1,1, coeff=ϵ_d))
# 	push!(data, adaga(2,2, coeff=ϵ_d))

# 	for i in 1:n-1
# 		push!(data, adaga(2i+1,2i+1, coeff=ws[i]))
# 		push!(data, adaga(2i+2, 2i+2, coeff=ws[i]))
# 		t = adaga(1, 2i+1, coeff=fs[i])
# 		push!(data, t)
# 		push!(data, t')
# 		t = adaga(2, 2i+2, coeff=fs[i])
# 		push!(data, t)
# 		push!(data, t')

# 		t = adagadag(2i+1, 2L+2i+2, coeff=Δ)
# 		push!(data, t)
# 		push!(data, t')
# 	end
# 	return GenericQuadraticHamiltonian(2L, data)
# end

cmatrix(m::Toulouse) = toulouse_cmatrix(m.bath, ϵ_d=m.ϵ_d)
function toulouse_cmatrix(b::AbstractDiscreteNormalBath{Fermion}; ϵ_d::Real)
	n = num_sites(b)
	m = zeros(Float64, n+1, n+1)
	m[1, 1] = ϵ_d
	ws, fs = frequencies(b), spectrumvalues(b)
	for i in 1:n
		m[1+i, 1+i] = ws[i]
		m[1, 1+i] = fs[i]
		m[1+i, 1] = fs[i]
	end
	return m
end

"""
	toulouse_cmatrix(b::AbstractDiscreteBCSBath; ϵ_d::Real)

The fermion modes are ordered as â₁†â-₁†â₂†â-₂† ⋯ âₖ†â-ₖ†  â₁â-₁â₂â-₂ ⋯ âₖâ-ₖ
"""
function toulouse_cmatrix(b::AbstractDiscreteBCSBath; ϵ_d::Real)
	n = div(num_sites(b), 2)
	L = n + 1
	Δ = b.Δ
	ws, fs = frequencies(b), spectrumvalues(b)

	h = zeros(typeof(Δ), 2L, 2L)
	h[1, 1] = ϵ_d
	h[2, 2] = ϵ_d
	for i in 1:n
		h[2i+1, 2i+1] = ws[i]
		h[2i+2, 2i+2] = ws[i]
		h[1, 2i+1] = fs[i]
		h[2i+1, 1] = fs[i]
		h[2, 2i+2] = fs[i]
		h[2i+2, 2] = fs[i]		
	end
	g = zeros(typeof(Δ), 2L, 2L)
	for i in 1:n
		g[2i+1, 2i+2] = -Δ
	end


	return bcs_cmatrix(h, g)
end


# thermal state
function thermocdm(m::NormalToulouse)
	h = cmatrix(hamiltonian(m, include_chemical=true))
	cache = eigencache(h)
	return fermionicthermocdm(cache, β=m.bath.β)
end
function thermocdm(m::BCSToulouse)
	h = cmatrix(hamiltonian(m))
	cache = eigencache(h)
	return fermionicthermocdm(cache, β=m.bath.β)
end
function thermodm(m::NormalToulouse)
	h = hamiltonian(m, include_chemical=true)
	return thermodm(h, β=m.bath.β)
end
function thermodm(m::BCSToulouse)
	h = hamiltonian(m)
	return thermodm(h, β=m.bath.β)
end

# separable state
separablecdm(m::Toulouse; nsys::Real=fermidirac(m.bath.β, 0, m.ϵ_d)) = separablecdm(m, nsys)
function separablecdm(m::NormalToulouse, nsys::Real)
	N = num_sites(m) 
	ρ = zeros(Float64, N, N)
	ρ[1, 1] = nsys

	Lj = num_sites(m.bath)
	ρ[2:N, 2:N] = thermocdm(m.bath)
	return ρ
end
function separabledm(m::NormalToulouse)
	h = freehamiltonian(m, include_chemical=true)
	return thermodm(h, β=m.bath.β)
end

function separablecdm(m::BCSToulouse, nsys::Real)
	L = div(num_sites(m), 2)
	ρ₀ = thermocdm(m.bath)

	ρ = zeros(eltype(ρ₀), 4L, 4L)
	ρ[1, 1] = nsys
	ρ[2, 2] = nsys
	ρ[2L+1, 2L+1] = 1-nsys
	ρ[2L+2, 2L+2] = 1-nsys	


	N = div(num_sites(m.bath), 2)
	ρ[3:2L, 3:2L] = ρ₀[1:2N, 1:2N]
	ρ[3:2L, 2L+3:4L] = ρ₀[1:2N, 2N+1:4N]
	ρ[2L+3:4L, 3:2L] = ρ₀[2N+1:4N, 1:2N]
	ρ[2L+3:4L, 2L+3:4L] = ρ₀[2N+1:4N, 2N+1:4N]
	return ρ
end
function separabledm(m::BCSToulouse)
	h = freehamiltonian(m)
	return thermodm(h, β=m.bath.β)
end

# green functions

"""
	toulouse_Gτ(b::AbstractDiscreteNormalBath{Fermion}; ϵ_d::Real)

The Matsubara Green's function for the Toulouse model
Gᵢⱼ(τ, τ′) = -Tr[âᵢ(τ)â†ⱼ(τ′)exp(-β(Ĥ-μN̂)⟩]
Gᵢⱼ(τ) = Gᵢⱼ(τ, 0)
"""
function toulouse_Gτ(b::AbstractDiscreteNormalBath{Fermion}; ϵ_d::Real)
	# @assert ishermitian(h)
	# ham = h - μ .* LinearAlgebra.I
	h = toulouse_cmatrix(b, ϵ_d=ϵ_d)
	β, μ = b.β, b.μ
	c = μ * one(h)
	c[1, 1] = 0
	ham = h - c
	# λs, U = eigen(Hermitian(ham))
	# cache = eigencache(ham)
	# ns = [fermidirac(β, 0, cache.λs[k]) for k in 1:length(cache.λs)]
	# return τ -> _fermionic_Gτ_util(cache, ns, 1, 1, τ)
	return freefermions_Gτ(ham, 1, 1, β=β)
end
# function toulouse_Gτ(b::AbstractDiscreteBCSBath; ϵ_d::Real)
# 	(b.μ == 0) || throw(ArgumentError("BCS bath should have μ=0"))
# 	# ham = h - μ .* LinearAlgebra.I
# 	h = toulouse_cmatrix(b, ϵ_d=ϵ_d)
# 	β, μ = b.β, b.μ
# 	# λs, U = eigen(Hermitian(ham))
# 	cache = eigencache(h)
# 	ns = [fermidirac(β, 0, cache.λs[k]) for k in 1:length(cache.λs)]
# 	return τ -> _fermionic_Gτ_util(cache, ns, 1, 1, τ)
# end
function toulouse_Gτ(b::AbstractDiscreteBCSBath; ϵ_d::Real)
	(b.μ == 0) || throw(ArgumentError("BCS bath should have μ=0"))
	# ham = h - μ .* LinearAlgebra.I
	h = toulouse_cmatrix(b, ϵ_d=ϵ_d)
	return freefermions_Gτ(h, 1, 1, β=b.β)
end
toulouse_Gτ(m::Toulouse) = toulouse_Gτ(m.bath; ϵ_d=m.ϵ_d)
toulouse_Gτ(m::Toulouse, τs::AbstractVector{<:Real}) = toulouse_Gτ(m.bath, τs, ϵ_d=m.ϵ_d)
function toulouse_Gτ(h::AbstractDiscreteBath{Fermion}, τs::AbstractVector{<:Real}; kwargs...)
	f = toulouse_Gτ(h; kwargs...)
	return f.(τs)
end


function _fermionic_Gτ_util(cache::EigenCache, ns, i::Int, j::Int, t::Real)
	λs, U = cache.λs, cache.U
	r_g = zero(eltype(U))
	for k in 1:size(U, 1)
		ss = U[i, k] * conj(U[j, k])
		exp_t = exp(-λs[k] * t)
		n_k = ns[k]
		r_g += ss * (1 - n_k) * exp_t
	end
	return r_g
end

function toulouse_Gt(m::Toulouse)
	f = toulouse_greater_lesser(m)
	return t -> begin
		gt, lt = f(t)
		return gt - lt
	end
end
function toulouse_Gt(m::Toulouse, ts::AbstractVector{<:Real})
	f = toulouse_Gt(m)
	return f.(ts)
end


function toulouse_greater_lesser(model::Toulouse)
	h = cmatrix(model)
	b = model.bath
	ρ = thermocdm(model)
	return freefermions_greater_lesser(h, ρ, 1)
end
function toulouse_greater_lesser(b::Toulouse, ts::AbstractVector{<:Real})
	f = toulouse_greater_lesser(b)
	r = f.(ts)
	return map(x->x[1], r), map(x->x[2], r)
end

function toulouse_neq_greater_lesser(b::Toulouse; kwargs...)
	h = cmatrix(b)
	ρ = separablecdm(b; kwargs...)
	return freefermions_greater_lesser(h, ρ, 1)
end 
function toulouse_neq_greater_lesser(b::Toulouse, ts::AbstractVector{<:Real}; kwargs...) 
	f = toulouse_neq_greater_lesser(b; kwargs...)
	r = f.(ts)
	return map(x->x[1], r), map(x->x[2], r)
end


# currents
function particlecurrent_cmatrix(m::NormalToulouse)
	N = num_sites(m)
	return _particlecurrent_util!(zeros(ComplexF64, N, N), m.bath, bathsites(m), 1)
end
function heatcurrent_cmatrix(m::NormalToulouse)
	N = num_sites(m)
	return _heatcurrent_util!(zeros(ComplexF64, N, N), m.bath, bathsites(m), 1)
end

bathsites(m::NormalToulouse) = 2:num_sites(m)

function _particlecurrent_util!(h::AbstractMatrix, b::AbstractDiscreteBath, bsites, band::Int)
	fs = spectrumvalues(b)
	for (j, v) in zip(bsites, fs)
		h[j, band] = -2*im*v
	end
	return h
end

function _heatcurrent_util!(h::AbstractMatrix, b::AbstractDiscreteBath, bsites, band::Int)
	ws, fs = frequencies(b), spectrumvalues(b)
	for (j, w, v) in zip(bsites, ws, fs)
		h[j, band] = -2*im*w*v
	end
	return h
end

function particlecurrent_hamiltonian(m::NormalToulouse)
	h = NormalQuadraticHamiltonian(ComplexF64, num_sites(m))
	return _particlecurrent_hamiltonian_util!(h, m.bath, bathsites(m), 1)
end
function heatcurrent_hamiltonian(m::NormalToulouse)
	h = NormalQuadraticHamiltonian(ComplexF64, num_sites(m))
	return _heatcurrent_hamiltonian_util!(h, m.bath, bathsites(m), 1)
end

function _particlecurrent_hamiltonian_util!(h::AbstractHamiltonian, b::AbstractDiscreteBath, bsites, band::Int)
	fs = spectrumvalues(b)
	for (j, v) in zip(bsites, fs)
		t = adaga(j, band, coeff=-2*im*v)
		push!(h, t)
	end
	return h	
end
function _heatcurrent_hamiltonian_util!(h::AbstractHamiltonian, b::AbstractDiscreteBath, bsites, band::Int)
	ws, fs = frequencies(b), spectrumvalues(b)
	for (j, w, v) in zip(bsites, ws, fs)
		t = adaga(j, band, coeff=-2*im*w*v)
		push!(h, t)
	end
	return h
end
