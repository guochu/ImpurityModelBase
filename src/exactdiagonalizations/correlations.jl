
# function correlation_2op_1t(h::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix, times::AbstractVector{<:Real}, cache::EigenCache=eigencache(h); 
# 							β::Real, reverse::Bool=false) 
# 	return correlation_2op_1t(h, A, B, thermocdm(cache, β=β), times, cache, reverse=reverse)
# end

"""
	correlation_2op_1t
	compute <a(t)b> if revere=false and <a b(t)> if reverse=true
"""
function correlation_2op_1t(h::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix, ρ::AbstractMatrix, times::AbstractVector{<:Real}, 
							cache::EigenCache=eigencache(h); reverse::Bool=false)
	tr_ρ = tr(ρ)
	return reverse ? [tr(A * timeevo(B, h, im*t, cache) * ρ) / tr_ρ  for t in times] : [tr(timeevo(A, h, im*t, cache) * B * ρ) / tr_ρ  for t in times]	
end

"""
	correlation_2op_1τ
	compute <a(τ)b> if revere=false and <a b(τ)> if reverse=true
"""
function correlation_2op_1τ(h::AbstractMatrix, A::AbstractMatrix, B::AbstractMatrix, times::AbstractVector{<:Real}, 
							cache::EigenCache=eigencache(h); β::Real, reverse::Bool=false)
	U, λs = cache.U, cache.λs
	tr_ρ = sum(exp.(-β .* λs))
	A′, B′ = U' * A * U, U' * B * U
	if reverse
		r = [tr(A′ * Diagonal(exp.(τ .* λs)) * B′ *  Diagonal(exp.(-(τ +β) .* λs))) / tr_ρ for τ in times]
	else
		r = [tr(Diagonal(exp.((τ -β) .* λs)) * A′ * Diagonal(exp.(-τ .* λs)) * B′ ) / tr_ρ for τ in times]
	end
	return r
end


# temporary (inefficient) solution for lindblad dynamics
struct LindbladOperator
	m::Array{ComplexF64, 4}
end

function lrmult(a, b)
	d = size(a, 1)
	r = kron(b, a)
	r4 = reshape(r, d, d, d, d)
	return permutedims(r4, (1,4,3,2))
end

# function lindbladoperator(H::AbstractMatrix, jumpops::Vector)
# 	I2 = one(H)
# 	@tensor L[1,4,2,3] := -im * H[1,2] * I2[3,4]
# 	@tensor L[1,4,2,3] += im * I2[1,2] * H[3,4]
# 	for A in jumpops
# 		AdagA = A' * A
# 		@tensor L[1,4,2,3] -= AdagA[1,2] * I2[3,4]
# 		@tensor L[1,4,2,3] -= I2[1,2] * AdagA[3,4]
# 		@tensor L[1,4,2,3] += 2.0 * A[1,2] * A'[3,4]
# 	end
# 	return LindbladOperator(L)
# end

function lindbladoperator(H::AbstractMatrix, jumpops::Vector)
	I2 = one(H)
	L = -im * lrmult(H, I2)
	L += im * lrmult(I2, H)
	for A in jumpops
		AdagA = A' * A
		L -= lrmult(AdagA, I2)
		L -= lrmult(I2, AdagA)
		L += 2.0 * lrmult(A, A')
	end
	return LindbladOperator(L)
end

function (op::LindbladOperator)(rho::AbstractMatrix)
	m = op.m
	d = size(m, 1)
	d2 = d * d
	r = reshape(m, d2, d2) * reshape(rho, d2)
	return reshape(r, d, d)
end

function timeevo(rho::AbstractMatrix, L::LindbladOperator, t::Real)
	d = size(rho, 1)
	d2 = d * d
	L2 = exp(t .* reshape(L.m, d2, d2))
	return reshape(L2 * reshape(rho, d2), d, d)
end

function correlation_2op_1t(h::LindbladOperator, A::AbstractMatrix, B::AbstractMatrix, ρ::AbstractMatrix, times::AbstractVector{<:Real}; 
							reverse::Bool=false)
	tr_ρ = tr(ρ)
	state = ρ
	if reverse
		state_right = state * A
	else
		state_right = B * state
	end
	r = ComplexF64[]
	for t in times
		state = timeevo(state_right, h, t)
		if reverse
			v = B * state
		else
			v = A * state
		end
		push!(r, tr(v) / tr_ρ)
	end
	return r
end
