# G(ω) for the Holstein model


function holstein_Gw(f::AbstractSpectrumFunction, ϵ::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8, order::Int=10)
	G0w(y) = toulouse_Gw(f, y, ϵ_d=ϵ_d, μ=μ, δ=δ)
	return holstein_G0w_to_Gw(G0w, ϵ, g=g, ω=ω, δ=δ)
end

function holstein_Σw(f::AbstractSpectrumFunction, ϵ::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8, order::Int=10)
	G0w(y) = toulouse_Gw(f, y, ϵ_d=ϵ_d, μ=μ, δ=δ)
	return holstein_G0w_to_Σw(G0w, ϵ, g=g, ω=ω, δ=δ)
end

function holstein_G0w_to_Gw(G0w::Function, ϵ::Real; g::Real, ω::Real, δ::Real=1.0e-8, order::Int=10)
	x = 0.
	for o in order:-1:0
		x = holstein_G0w_t0_Gw_util(G0w, ϵ, o, x, g=g, ω=ω, δ=δ)
	end
	return x
end

function holstein_G0w_to_Σw(G0w::Function, ϵ::Real; g::Real, ω::Real, δ::Real=1.0e-8, order::Int=10)
	x = 0.
	for o in order:-1:1
		x = holstein_G0w_t0_Gw_util(G0w, ϵ, o, x, g=g, ω=ω, δ=δ)
	end
	return x
end

function holstein_Gw_util(f::AbstractSpectrumFunction, ϵ::Real, order::Int, x; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8)
	G0w(y) = toulouse_Gw(f, y, ϵ_d=ϵ_d, μ=μ, δ=δ)
	return holstein_G0w_t0_Gw_util(G0w, ϵ, order, x, g=g, ω=ω)
end

function holstein_G0w_t0_Gw_util(G0w, ϵ::Real, order::Int, x; g::Real, ω::Real, δ::Real=1.0e-8)
	deno = (order == 0) ? 1 : order * g^2
	return deno / (1/G0w(ϵ-order*ω) + im*δ - x)
end