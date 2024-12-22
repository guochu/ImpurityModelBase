"""
	Gτ_to_Giw(gτ; kwargs...)

For δτ = β/N, we assume Gτ is a vector of length N+1, with elements 
G(0), G(δτ), …, G(N*δτ) = G(β), with G(β) = 1-G(0)
The last element is thus redunant, and only the first N elementes are 
used for the Fourier transfermation

The imaginary-frequency points are 2(n-1)π/β for n in -nmax:nmax+1 
Therefore there are 2nmax+2 points in total
"""
Gτ_to_Giw(gτ::Vector{<:Real}; β::Real, n::Int=1000) = ifourier(gτ, β=β, n=n)

"""
	Δτ_to_Δiw(Δτ; kwargs...)

Similar to Gτ_to_Giw, convert the imaginary-time axis hybridization 
function Δτ to the imaginary-frequency domain
"""
Δτ_to_Δiw(gτ::AbstractVector{<:Real}; kwargs...) = Gτ_to_Giw(gτ, ; kwargs...)

"""
	Giw_to_Gτ(Giw; β, N)

Convert Giw to the imaginary-time domain with δτ=β/N
"""
function Giw_to_Gτ(Giw::AbstractVector{<:Number}; β::Real, Nτ::Int)
	iseven(length(Giw)) || throw(ArgumentError("even number of frequencies expected"))
	nmax = div(length(Giw), 2) - 1
	δτ = β / Nτ
	f(τ) = sum((Giw[i]-1/(im*(2n-1)*π/β)) * exp(-im*τ*(2n-1)*π/β) for (i, n) in enumerate(-nmax:nmax+1))
	gτ = zeros(Float64, Nτ+1)
	for i in 1:Nτ
		tmp = f((i-1) * δτ)
		(abs(imag(tmp)) < 1.0e-8) || error("imaginary part of Gτ is too large")
		gτ[i] = real(tmp)/β - 0.5
	end
	gτ[end] = 1 - gτ[1]
	return gτ
end

"""
	ifrequency(β, n)

Return imaginary-time frequency (2n-1)π/β
"""
ifrequency(β::Real, n::Int) = (2*n-1)*π/β
ifrequencies(β::Real, nmax::Int=1000) = [ifrequency(β, n) for n in -nmax:nmax+1]
ifrequencies(Giw::AbstractVector; β::Real) = ifrequencies(β, div(length(Giw), 2)-1)

function ifourier(gτ::Vector{<:Real}; β::Real, n::Int=1000)
	# (gτ[1] + gτ[end] ≈ 1) || throw(ArgumentError("sum of the first and last elementes should be 1"))
	Nτ = length(gτ)-1
	δτ = β / Nτ
	f(ω) = sum(gτ[k]*exp(im*(k-1)*δτ*ω) for k in 1:Nτ+1)*δτ
	return [f((2*nj-1)*π/β) for nj in -n:n+1]
end


# function _imagfourier(τs::Vector{<:Real}, gτ::Vector{<:Real}, w::Real)
# 	Nτ = length(gτ)-1
# 	f(ω) = begin
# 		r = 0.
# 		for k in 1:Nτ
# 			dτ = τs[k+1] - τs[k]
# 			r += gτ[k] * exp(im*(k-1)*τs[k]*ω) * dτ
# 		end
# 		return r
# 	end
# 	return f(w)
# end

# function _inverse_imagfourier(ws::Vector{<:Real}, fs::Vector{<:Number}, τ::Real)
# 	r = zero(eltype(fs))
# 	for i in 1:length(ws)-1
# 		dw = ws[i+1] - ws[i]
# 		r += (Giw[i] - 1/(im *ws[i])) * exp(-im*τ*ws[i]) * dw
# 	end
# 	(abs(imag(r)) < 1.0e-8) || error("imaginary part of Gτ is too large")
# 	return real(tmp) / (2π) - 0.5	
# end
