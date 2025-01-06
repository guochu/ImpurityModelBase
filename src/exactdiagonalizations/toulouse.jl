struct Toulouse{B<:AbstractDiscreteBath}
	bath::B
	ϵ_d::Float64
end
Toulouse(b::AbstractDiscreteBath; ϵ_d::Real) = Toulouse(b, convert(Float64, ϵ_d))
Base.eltype(::Type{Toulouse{B}}) where B = Float64
Base.eltype(x::Toulouse) = eltype(typeof(x))

num_sites(m::Toulouse) = num_sites(m.bath) + 1
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
	λs, U = eigen(Hermitian(ham))
	ns = [fermidirac(β, 0, λs[k]) for k in 1:length(λs)]
	return τ -> _fermionic_Gτ_util(U, λs, ns, 1, 1, τ)
end
toulouse_Gτ(m::Toulouse) = toulouse_Gτ(m.bath; ϵ_d=m.ϵ_d)
toulouse_Gτ(m::Toulouse, τs::AbstractVector{<:Real}) = toulouse_Gτ(m.bath, τs, ϵ_d=m.ϵ_d)
function toulouse_Gτ(h::AbstractDiscreteFermionicBath, τs::AbstractVector{<:Real}; kwargs...)
	f = toulouse_Gτ(h; kwargs...)
	return f.(τs)
end


function _fermionic_Gτ_util(U, λs, ns, i::Int, j::Int, t::Real)
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


function toulouse_greater_lesser(b::AbstractDiscreteFermionicBath; ϵ_d::Real)
	h = toulouse_cmatrix(b, ϵ_d=ϵ_d)
	return freefermions_greater_lesser(h, β=b.β, μ=b.μ, i=1)
end
function toulouse_greater_lesser(b::AbstractDiscreteFermionicBath, ts::AbstractVector{<:Real}; kwargs...)
	f = toulouse_greater_lesser(b; kwargs...)
	r = f.(ts)
	return map(x->x[1], r), map(x->x[2], r)
end

function separablestate(m::Toulouse, nsys::Real)
	N = num_sites(m) 
	ρ = zeros(Float64, N, N)
	ρ[1, 1] = nsys

	Lj = num_sites(m.bath)
	ρ[2:N, 2:N] = thermalstate(m.bath)
	return ρ
end

function thermalstate(m::Toulouse)
	# β, μ = m.bath.β, m.bath.μ
	h = cmatrix(m)
	evals, U = eigen(Hermitian(h))
	evals = [thermaloccupation(bath, item) for item in evals]
	return U * Diagonal(evals) * U'
end

function particlecurrent_cmatrix(m::Toulouse)
	N = num_sites(m)
	return _particlecurrent_util!(zeros(ComplexF64, N, N), m.bath, bathsites(m), 1)
end

bathsites(m::Toulouse) = 2:num_sites(m)

function _particlecurrent_util!(h::AbstractMatrix, b::AbstractDiscreteBath, bsites, band::Int)
	fs = spectrumvalues(b)
	for (j, v) in zip(bsites, fs)
		h[j, band] = -2*im*v
	end
	return h
end