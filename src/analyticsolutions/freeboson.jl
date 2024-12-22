"""
	freeboson_greater(t::Real; β::Real, ω::Real)

Greater Green's function of a free particle
"""
freeboson_greater(t::Real; β::Real, ω::Real) = -im*exp(-im*ω*t)/(1-exp(-safe_mult(β, ω)))

"""
	freeboson_lesser(t::Real; β::Real, μ::Real)

Lesser Green's function of a free particle
"""
function freeboson_lesser(t::Real; β::Real, ω::Real)
	x = exp(-safe_mult(β, ω))
	return im*exp(-im*ω*t) * (x / (1-x))
end

"""
	freeboson_Gt(t::Real; kwargs...)

Retarded Green's function of a free particle
"""
freeboson_Gt(t::Real; kwargs...) = freeboson_greater(t; kwargs...) + freeboson_lesser(t; kwargs...)

"""
	freeboson_Gτ(τ::Float64; β::Real, ω::Real)

Matsubara Green's function of a free particle in the imaginary-time axis
"""
freeboson_Gτ(τ::Float64; β::Real, ω::Real) = -exp(-τ*ω)/(1-exp(-β*ω))

"""
	freeboson_occupation(β::Real, ω::Real)

The average occupation n̄ = ⟨n̂⟩ of a free boson in thermal state
"""
function freeboson_occupation(β::Real, ω::Real)
	x = exp(-safe_mult(β, ω))
	return x / (1-x)
end