"""
	Gτ_to_Giw

For δτ = β/N, we assume Gτ is a vector of length N+1, with elements 
G(0), G(δτ), …, G(N*δτ) = G(β), with G(β) - 1-G(0)
The last element is thus redunant, and only the first N elementes are 
used for the Fourier transfermation

"""
function Gτ_to_Giw(gτ::AbstractVector{<:Real}; kwargs...)
	(gτ[1] + gτ[end] ≈ 1) || throw(ArgumentError("sum of the first and last elementes should be 1"))
	return ifourier(gτ; kwargs...)
end
Δτ_to_Δiw(gτ::AbstractVector{<:Real}; kwargs...) = ifourier(gτ; kwargs...)

function Giw_to_Gτ(Giw::AbstractVector{<:Number}; β::Real, N::Int)
	iseven(length(Giw)) || throw(ArgumentError("even number of frequencies expected"))
	nmax = div(length(Giw), 2) - 1
	δτ = β / N
	f(τ) = sum((Giw[i]-1/(im*(2n-1)*π/β)) * exp(-im*τ*(2n-1)*π/β) for (i, n) in enumerate(-nmax:nmax+1))
	gτ = zeros(Float64, N+1)
	for i in 1:N
		tmp = f((i-1) * δτ)
		(abs(imag(tmp)) < 1.0e-8) || error("imaginary part of Gτ is too large")
		gτ[i] = real(tmp)/β - 0.5
	end
	gτ[end] = 1 - gτ[1]
	return gτ
end

ifrequency(β::Real, n::Int) = (2*n-1)*π/β
ifrequencies(β::Real, nmax::Int=1000) = [ifrequency(β, n) for n in -nmax:nmax+1]
ifrequencies(Giw::AbstractVector; β::Real) = ifrequencies(β, div(length(Giw), 2)-1)

function ifourier(gτ::Vector{<:Real}; β::Real, nmax::Int=1000)
	# (gτ[1] + gτ[end] ≈ 1) || throw(ArgumentError("sum of the first and last elementes should be 1"))
	Nτ = length(gτ)-1
	δτ = β / Nτ
	f(ω) = sum(gτ[k]*exp(im*(k-1)*δτ*ω) for k in 1:Nτ+1)*δτ
	return [f((2*n-1)*π/β) for n in -nmax:nmax+1]
end