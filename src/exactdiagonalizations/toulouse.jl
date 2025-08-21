struct Toulouse{B<:AbstractDiscreteBath{Fermion}}
	bath::B
	ϵ_d::Float64
end
Toulouse(b::AbstractDiscreteBath{Fermion}; ϵ_d::Real) = Toulouse(b, convert(Float64, ϵ_d))
Base.eltype(::Type{Toulouse{B}}) where B = eltype(B)
Base.eltype(x::Toulouse) = eltype(typeof(x))

num_sites(m::Toulouse) = num_sites(m.bath) + 1

function hamiltonian(h::Toulouse{<:AbstractDiscreteFermionicBath})
	T = eltype(h)
	data = Vector{AdagATerm{T}}[]
	L = num_sites(h)
	ϵ_d = h.ϵ_d
	push!(data, adaga(1, 1, coeff=ϵ_d))
	bath = h.bath
	ws, fs = frequencies(bath), spectrumvalues(bath)
	for i in 1:L-1
		push!(data, adaga(1+i, 1+i, coeff=ws[i]))
		t = adaga(1, 1+i, coeff=fs[i])
		push!(data, t)
		push!(data, t')
	end
	return NormalQuadraticHamiltonian(L, data)
end
function hamiltonian(h::Toulouse{<:AbstractDiscreteBCSBath})
	T = eltype(h)
	data = Vector{QuadraticTerm{T}}[]

	L = num_sites(h)
	bath = h.bath
	ws, fs = frequencies(bath), spectrumvalues(bath)	
	Δ = bath.Δ
	ϵ_d = h.ϵ_d

	push!(data, adaga(1,1, coeff=ϵ_d))
	push!(data, adaga(2,2, coeff=ϵ_d))

	for i in 1:n-1
		push!(data, adaga(2i+1,2i+1, coeff=ws[i]))
		push!(data, adaga(2i+2, 2i+2, coeff=ws[i]))
		t = adaga(1, 2i+1, coeff=fs[i])
		push!(data, t)
		push!(data, t')
		t = adaga(2, 2i+2, coeff=fs[i])
		push!(data, t)
		push!(data, t')

		t = adagadag(2i+1, 2L+2i+2, coeff=Δ)
		push!(data, t)
		push!(data, t')
	end
	return GenericQuadraticHamiltonian(2L, data)
end

"""
	toulouse_Gτ(b::AbstractDiscreteFermionicBath; ϵ_d::Real)

The Matsubara Green's function for the Toulouse model
Gᵢⱼ(τ, τ′) = -Tr[âᵢ(τ)â†ⱼ(τ′)exp(-β(Ĥ-μN̂)⟩]
Gᵢⱼ(τ) = Gᵢⱼ(τ, 0)
"""
function toulouse_Gτ(b::AbstractDiscreteFermionicBath; ϵ_d::Real)
	# @assert ishermitian(h)
	# ham = h - μ .* LinearAlgebra.I
	h = toulouse_cmatrix(b, ϵ_d=ϵ_d)
	β, μ = b.β, b.μ
	c = μ * one(h)
	c[1, 1] = 0
	ham = h - c
	# λs, U = eigen(Hermitian(ham))
	cache = eigencache(ham)
	ns = [fermidirac(β, 0, cache.λs[k]) for k in 1:length(cache.λs)]
	return τ -> _fermionic_Gτ_util(cache, ns, 1, 1, τ)
end
function toulouse_Gτ(b::AbstractDiscreteBCSBath; ϵ_d::Real)
	(b.μ == 0) || throw(ArgumentError("BCS bath should have μ=0"))
	# ham = h - μ .* LinearAlgebra.I
	h = toulouse_cmatrix(b, ϵ_d=ϵ_d)
	β, μ = b.β, b.μ
	# λs, U = eigen(Hermitian(ham))
	cache = eigencache(h)
	ns = [fermidirac(β, 0, cache.λs[k]) for k in 1:length(cache.λs)]
	return τ -> _fermionic_Gτ_util(cache, ns, 1, 1, τ)
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


cmatrix(m::Toulouse) = toulouse_cmatrix(m.bath, ϵ_d=m.ϵ_d)
function toulouse_cmatrix(b::AbstractDiscreteFermionicBath; ϵ_d::Real)
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
	n = num_sites(b)
	L = n + 1
	Δ = b.Δ

	m = zeros(typeof(Δ), 4L, 4L)
	m[1, 1] = ϵ_d
	m[2, 2] = ϵ_d
	m[2L+1, 2L+1] = 1 - ϵ_d
	m[2L+2, 2L+2] = 1 - ϵ_d
	ws, fs = frequencies(b), spectrumvalues(b)

	for i in 1:n
		m[2i+1, 2i+1] = ws[i]
		m[2i+2, 2i+2] = ws[i]
		m[1, 2i+1] = fs[i]
		m[2i+1, 1] = fs[i]
		m[2, 2i+2] = fs[i]
		m[2i+2, 2] = fs[i]


		m[2L+2i+1, 2L+2i+1] = 1-ws[i]
		m[2L+2i+2, 2L+2i+2] = 1-ws[i]
		m[2L+1, 2L+2i+1] = 1-fs[i]
		m[2L+2i+1, 2L+1] = 1-fs[i]
		m[2L+2, 2L+2i+2] = 1-fs[i]
		m[2L+2i+2, 2L+2] = 1-fs[i]

		m[2i+1, 2L+2i+2] = -Δ
		m[2L+2i+2, 2i+1] = -conj(Δ)
	end

	return m
end

toulouse_Gt(m::Toulouse) = toulouse_Gt(m.bath, ϵ_d=m.ϵ_d)
toulouse_Gt(m::Toulouse, ts::AbstractVector{<:Real}) = toulouse_Gt(m.bath, ts, ϵ_d=m.ϵ_d)
function toulouse_Gt(b::AbstractDiscreteFermionicBath; kwargs...)
	f = toulouse_greater_lesser(b; kwargs...)
	return t -> begin
		gt, lt = f(t)
		return gt - lt
	end
end
function toulouse_Gt(b::AbstractDiscreteFermionicBath, ts::AbstractVector{<:Real}; kwargs...)
	f = toulouse_Gt(b; kwargs...)
	return f.(ts)
end

toulouse_greater_lesser(b::Toulouse) = toulouse_greater_lesser(b.bath, ϵ_d=b.ϵ_d)
toulouse_greater_lesser(b::Toulouse, ts::AbstractVector{<:Real}) = toulouse_greater_lesser(b.bath, ts, ϵ_d=b.ϵ_d)
function toulouse_greater_lesser(b::AbstractDiscreteFermionicBath; ϵ_d::Real)
	h = toulouse_cmatrix(b, ϵ_d=ϵ_d)
	return freefermions_greater_lesser(h, 1, β=b.β, μ=b.μ)
end
function toulouse_greater_lesser(b::AbstractDiscreteFermionicBath, ts::AbstractVector{<:Real}; kwargs...)
	f = toulouse_greater_lesser(b; kwargs...)
	r = f.(ts)
	return map(x->x[1], r), map(x->x[2], r)
end

function toulouse_neq_greater_lesser(b::Toulouse; nsys::Real=0)
	h = cmatrix(b)
	ρ = separablecdm(b, nsys)
	return freefermions_greater_lesser(h, ρ, 1)
end 
function toulouse_neq_greater_lesser(b::Toulouse, ts::AbstractVector{<:Real}; nsys::Real=0) 
	f1, f2 = toulouse_neq_greater_lesser(b, nsys=nsys)
	return f1.(ts), f2.(ts)
end

function separablecdm(m::Toulouse{<:AbstractDiscreteFermionicBath}, nsys::Real)
	N = num_sites(m) 
	ρ = zeros(Float64, N, N)
	ρ[1, 1] = nsys

	Lj = num_sites(m.bath)
	ρ[2:N, 2:N] = thermocdm(m.bath)
	return ρ
end

function separablecdm(m::Toulouse{<:AbstractDiscreteBCSBath}, nsys::Real)
	L = num_sites(m) 
	ρ₀ = thermocdm(m.bath)

	ρ = zeros(eltype(ρ₀), 4L, 4L)
	ρ[1, 1] = nsys
	ρ[2, 2] = nsys
	ρ[2L+1, 2L+1] = 1-nsys
	ρ[2L+2, 2L+2] = 1-nsys	


	N = num_sites(m.bath)
	ρ[3:2L, 3:2L] = ρ₀[1:2N, 1:2N]
	ρ[2L+3:4L, 2L+3:4L] = ρ₀[2N+1:4N, 2N+1:4N]
	return ρ
end

function thermocdm(m::Toulouse)
	# β, μ = m.bath.β, m.bath.μ
	h = cmatrix(m)
	cache = eigencache(h)
	# evals = [thermaloccupation(m.bath, item) for item in cache.λs]
	# return cache.U * Diagonal(evals) * cache.U'
	return thermocdm(cache, β=m.bath.β, μ=m.bath.μ)
end

function particlecurrent_cmatrix(m::Toulouse)
	N = num_sites(m)
	return _particlecurrent_util!(zeros(ComplexF64, N, N), m.bath, bathsites(m), 1)
end
function heatcurrent_cmatrix(m::Toulouse)
	N = num_sites(m)
	return _heatcurrent_util!(zeros(ComplexF64, N, N), m.bath, bathsites(m), 1)
end

bathsites(m::Toulouse) = 2:num_sites(m)

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