
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