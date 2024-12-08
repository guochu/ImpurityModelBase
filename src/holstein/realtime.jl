# G(ω) for the Holstein model


function holstein_Gw(f::AbstractSpectrumFunction, ϵ::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8, order::Int=10)
	x = 0.
	for o in order:-1:0
		x = holstein_Gw_util(f, ϵ, o, x, g=g, ω=ω, ϵ_d=ϵ_d, μ=μ, δ=δ)
	end
	return x
end

function holstein_Σw(f::AbstractSpectrumFunction, ϵ::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8, order::Int=10)
	x = 0.
	for o in order:-1:1
		x = holstein_Gw_util(f, ϵ, o, x, g=g, ω=ω, ϵ_d=ϵ_d, μ=μ, δ=δ)
	end
	return x
end

function holstein_Gw_util(f::AbstractSpectrumFunction, ϵ::Real, order::Int, x; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8)
	G0w(y) = toulouse_Gw(f, y, ϵ_d=ϵ_d, μ=μ, δ=δ)
	deno = (order == 0) ? 1 : order * g^2
	return deno / (1/G0w(ϵ-order*ω) - x)
end