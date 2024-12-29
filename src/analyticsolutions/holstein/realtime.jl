# G(ω) for the Holstein model

# function holstein_Gw(f::AbstractSpectrumFunction, ϵ::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8, order::Int=10)
# 	G0w(y) = toulouse_Gw(f, y, ϵ_d=ϵ_d, μ=μ, δ=δ)
# 	return holstein_G0w_to_Gw(G0w, ϵ, g=g, ω=ω, δ=δ)
# end

# function holstein_Σw(f::AbstractSpectrumFunction, ϵ::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8, order::Int=10)
# 	G0w(y) = toulouse_Gw(f, y, ϵ_d=ϵ_d, μ=μ, δ=δ)
# 	return holstein_G0w_to_Σw(G0w, ϵ, g=g, ω=ω, δ=δ)
# end

# function holstein_G0w_to_Gw(G0w::Function, ϵ::Real; g::Real, ω::Real, δ::Real=1.0e-8, order::Int=10)
# 	x = 0.
# 	for o in order:-1:0
# 		x = holstein_G0w_t0_Gw_util(G0w, ϵ, o, x, g=g, ω=ω, δ=δ)
# 	end
# 	return x
# end

# function holstein_G0w_to_Σw(G0w::Function, ϵ::Real; g::Real, ω::Real, δ::Real=1.0e-8, order::Int=10)
# 	x = 0.
# 	for o in order:-1:1
# 		x = holstein_G0w_t0_Gw_util(G0w, ϵ, o, x, g=g, ω=ω, δ=δ)
# 	end
# 	return x
# end

# function holstein_Gw_util(f::AbstractSpectrumFunction, ϵ::Real, order::Int, x; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8)
# 	G0w(y) = toulouse_Gw(f, y, ϵ_d=ϵ_d, μ=μ, δ=δ)
# 	return holstein_G0w_t0_Gw_util(G0w, ϵ, order, x, g=g, ω=ω)
# end

# function holstein_G0w_t0_Gw_util(G0w, ϵ::Real, order::Int, x; g::Real, ω::Real, δ::Real=1.0e-8)
# 	deno = (order == 0) ? 1 : order * g^2
# 	return deno / (1/G0w(ϵ-order*ω) + im*δ - x)
# end


struct GreenFunction{T<:Number, I} <: Function
	ws::Vector{Float64}
	fs::Vector{T}
	δ::Float64
	interp::I
end
function GreenFunction(ws::AbstractVector{<:Real}, fs::AbstractVector{T}; δ::Real=1.0e-8) where {T<:Number}
	ws = convert(Vector{Float64}, ws)
	fs = convert(Vector{T}, fs)
	interp = linear_interpolation(ws, fs)
	return GreenFunction(ws, fs, convert(Float64, δ), interp)
end

function (x::GreenFunction)(ω::Real)
	if (ω < first(x.ws)) || (ω > last(x.ws))
		return 1/(ω+im*x.δ)
	end
	return x.interp(ω)
end

function bethe_holstein_dmft_iteration(G0w::Function, wk::Real; g::Real, ω::Real, t::Real=1, δ::Real=1.0e-8, order::Int=10)
	Σk = holstein_G0w_to_Σw(G0w, wk, g=g, ω=ω, δ=δ, order=order) 
	a = -(Σk-wk)
	b = sqrt((Σk-wk)^2-t^2)
	Gk = (a + b) / (t^2/2)
	if imag(Gk) > 0
		Gk = (a - b) / (t^2/2)
	end
	return Gk
end

function holstein_Gt(f::AbstractSpectrumFunction, t::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, wmax::Real=20., wmin::Real=-wmax, δ::Real=1.0e-8, order::Int=10)
    A = quadgkwrapper(bounded(ϵ -> (holstein_Gw(f, ϵ; ϵ_d=ϵ_d, g=g, ω=ω, μ=μ, δ=δ, order=order)-1.0/(ϵ+im*δ))*exp(-im*ϵ*t), wmin, wmax))
    return im*(A/(2π)-im)
end

function holstein_Gw(f::AbstractSpectrumFunction, ϵ::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8, order::Int=10)
	G0w(y) = toulouse_Gw(f, y, ϵ_d=ϵ_d, μ=μ, δ=δ)
	return holstein_G0w_to_Gw(G0w, ϵ, g=g, ω=ω, δ=δ, order=order)
end

function holstein_Σw(f::AbstractSpectrumFunction, ϵ::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, δ::Real=1.0e-8, order::Int=10)
	G0w(y) = toulouse_Gw(f, y, ϵ_d=ϵ_d, μ=μ, δ=δ)
	return holstein_G0w_to_Σw(G0w, ϵ, g=g, ω=ω, δ=δ, order)
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

function holstein_G0w_t0_Gw_util(G0w::Function, ϵ::Real, order::Int, x; g::Real, ω::Real, δ::Real=1.0e-8)
	deno = (order == 0) ? 1 : order * g^2
	return deno / (1/G0w(ϵ-order*ω) + im*δ - x)
end


