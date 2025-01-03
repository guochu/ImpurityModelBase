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