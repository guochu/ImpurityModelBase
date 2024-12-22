# imaginary-time Green's function of Toulouse model

"""
    toulouse_Giw(spectrum::AbstractSpectrumFunction, ω::Real; ϵ_d::Real, μ::Real)

Matsubara Green's function in the imaginary-frequency axis

ϵ_d is the on-site energy of the localized electron
μ is the chemical potental of the bath
β is not a parameter as Gw is independent of β for the Toulouse model
"""
function toulouse_Giw(f::AbstractSpectrumFunction, ω::Real; ϵ_d::Real, μ::Real=0)
    g(ε) = im*ω+μ-ε
    1.0/(im*ω-ϵ_d-quadgkwrapper(f / g))
end
function toulouse_Giw(spectrum::AbstractSpectrumFunction; β::Real, ϵ_d::Real, μ::Real=0, n::Int=1000)
    return [toulouse_Giw(spectrum, (2*nj-1)*π/β; ϵ_d=ϵ_d, μ=μ) for nj in ifrequencies(β, n)]
end
toulouse_Giw(bath::AbstractFermionicBath, ω::Real; kwargs...) = toulouse_Giw(bath.spectrum, ω; μ=bath.μ, kwargs...)


"""
    toulouse_Gτ(spectrum::AbstractSpectrumFunction, τ::Real; β::Real, ϵ_d::Real, μ::Real, n::Int)

Matsubara Green's function in the imaginary-time axis
n+1 is the half the number of the imaginary frequency points (-n:n+1),
for each n, the imaginary frequency point is (2n-1)π/β
"""
function toulouse_Gτ(spectrum::AbstractSpectrumFunction, τ::Real; β::Real, ϵ_d::Real, μ::Real=0., n::Int=1000)
    res = 0.0
    for nj in ifrequencies(β, n)
        ω = (2nj-1)*π/β
        res += (toulouse_Giw(spectrum, ω; ϵ_d=ϵ_d, μ=μ)-1/(im*ω))*exp(-im*τ*ω)
    end
    res = -(res/β-0.5)
end
toulouse_Gτ(bath::AbstractFermionicBath, τ::Real; kwargs...) = toulouse_Gτ(bath.spectrum, τ; β=bath.β, μ=bath.μ, kwargs...)
function toulouse_Gτ(spectrum::AbstractSpectrumFunction; β::Real, Nτ::Int, ϵ_d::Real, μ::Real=0., n::Int=1000)
    δτ = β / Nτ
    gτ = zeros(Float64, Nτ+1)
    for i in 1:Nτ
        τ = (i-1) * δτ
        tmp = toulouse_Gτ(spectrum, τ, β=β, ϵ_d=ϵ_d, μ=μ, n=n)
        (abs(imag(tmp)) < 1.0e-8) || error("imaginary part of Gτ is too large")
        gτ[i] = real(tmp)
    end
    gτ[end] = 1 - gτ[1]
    return gτ
end
toulouse_Gτ(bath::AbstractFermionicBath; Nτ::Int, kwargs...) = toulouse_Gτ(bath.spectrum; β=bath.β, μ=bath.μ, Nτ=Nτ, kwargs...)

# function toulouse_Δiw(spectrum::AbstractSpectrumFunction; β::Real, n::Int=1000)
#     f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
#     ff(ω) = quadgkwrapper(bounded(ϵ -> f(ϵ) / (im*ω - ϵ), lb, ub))

#     return [ff((2*n-1)*π/β) for n in -n:n+1]
# end

"""
    toulouse_Δiw(f::AbstractSpectrumFunction; β::Real, n::Int=1000)

The hybridization function in the imaginary-frequency axis for the Toulouse model

f is the bath spectrum density
"""
function toulouse_Δiw(f::AbstractSpectrumFunction; β::Real, n::Int=1000)
    function ff(ω)
        g(ϵ) = im*ω - ϵ
        return quadgkwrapper(f / g)
    end 

    return [ff((2*nj-1)*π/β) for nj in ifrequencies(β, n)]
end
toulouse_Δiw(bath::AbstractFermionicBath; n::Int=1000) = toulouse_Δiw(bath.spectrum, β=bath.β, n=n)

# function toulouse_Δτ(spectrum::AbstractSpectrumFunction; β::Real, N::Int)
#     f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
#     δτ = β / N
#     ff(τ) = quadgkwrapper(bounede(ϵ -> -f(ϵ) / (exp(-ϵ*τ) / (1+exp(-β*ϵ)) ), lb, ub))
#     return [ff(i*δτ) for i in 0:N]
# end

"""
    toulouse_Δτ(f::AbstractSpectrumFunction; β::Real, Nτ::Int)

The hybridization function in the imaginary-time axis for the Toulouse model

f is the bath spectrum density
"""
function toulouse_Δτ(f::AbstractSpectrumFunction; β::Real, Nτ::Int)
    δτ = β / Nτ
    function ff(τ)
        g(ϵ) = -exp(-ϵ*τ) / (1+exp(-β*ϵ)) 
        return quadgkwrapper(f / g)
    end 
    return [ff(i*δτ) for i in 0:Nτ]
end
toulouse_Δτ(bath::AbstractFermionicBath; Nτ::Int) = toulouse_Δτ(bath.spectrum, β=bath.β, Nτ=Nτ)


