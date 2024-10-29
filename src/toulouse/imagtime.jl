# imaginary-time Green's function of Toulouse model

# function toulouse_Giw(spectrum::AbstractSpectrumFunction, ω::Real; ϵ_d::Real, μ::Real=0)
#     f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
#     1.0/(im*ω-ϵ_d-quadgkwrapper(bounded(ε -> f(ε)/(im*ω+μ-ε), lb=lb, ub=ub)))
# end

"""
    toulouse_Giw(spectrum::AbstractSpectrumFunction, ω::Real; ϵ_d::Real, μ::Real)

Matsubara Green's function in the frequency axis
"""
function toulouse_Giw(f::AbstractSpectrumFunction, ω::Real; ϵ_d::Real, μ::Real=0)
    g(ε) = im*ω+μ-ε
    1.0/(im*ω-ϵ_d-quadgkwrapper(f / g))
end


function toulouse_Giw(spectrum::AbstractSpectrumFunction; β::Real, ϵ_d::Real, μ::Real=0, nmax::Int=1000)
    return [toulouse_Giw(spectrum, (2*n-1)*π/β; ϵ_d=ϵ_d, μ=μ) for n in -nmax:nmax+1]
end

toulouse_Giw(bath::AbstractFermionicBath, ω::Real; kwargs...) = toulouse_Giw(bath.spectrum, ω; μ=bath.μ, kwargs...)


"""
    toulouse_Gτ(spectrum::AbstractSpectrumFunction, τ::Real; β::Real, ϵ_d::Real, μ::Real, nmax::Int)

Matsubara Green's function in the imaginary-time axis
"""
function toulouse_Gτ(spectrum::AbstractSpectrumFunction, τ::Real; β::Real, ϵ_d::Real, μ::Real=0., nmax::Int=1000)
    res = 0.0
    for n = -nmax:nmax+1
        ω = (2n-1)*π/β
        res += (toulouse_Giw(spectrum, ω; ϵ_d=ϵ_d, μ=μ)-1/(im*ω))*exp(-im*τ*ω)
    end
    res = -(res/β-0.5)
end
toulouse_Gτ(bath::AbstractFermionicBath, τ::Real; kwargs...) = toulouse_Gτ(bath.spectrum, τ; β=bath.β, μ=bath.μ, kwargs...)
function toulouse_Gτ(spectrum::AbstractSpectrumFunction; β::Real, N::Int, ϵ_d::Real, μ::Real=0., nmax::Int=1000)
    δτ = β / N
    gτ = zeros(Float64, N+1)
    for i in 1:N
        τ = (i-1) * δτ
        tmp = toulouse_Gτ(spectrum, τ, β=β, ϵ_d=ϵ_d, μ=μ, nmax=nmax)
        (abs(imag(tmp)) < 1.0e-8) || error("imaginary part of Gτ is too large")
        gτ[i] = real(tmp)
    end
    gτ[end] = 1 - gτ[1]
    return gτ
end
toulouse_Gτ(bath::AbstractFermionicBath; N::Int, kwargs...) = toulouse_Gτ(bath.spectrum; β=bath.β, μ=bath.μ, N=N, kwargs...)

# function toulouse_Δiw(spectrum::AbstractSpectrumFunction; β::Real, nmax::Int=1000)
#     f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
#     ff(ω) = quadgkwrapper(bounded(ϵ -> f(ϵ) / (im*ω - ϵ), lb, ub))

#     return [ff((2*n-1)*π/β) for n in -nmax:nmax+1]
# end

function toulouse_Δiw(f::AbstractSpectrumFunction; β::Real, nmax::Int=1000)
    function ff(ω)
        g(ϵ) = im*ω - ϵ
        return quadgkwrapper(f / g)
    end 

    return [ff((2*n-1)*π/β) for n in -nmax:nmax+1]
end
toulouse_Δiw(bath::AbstractFermionicBath; nmax::Int=1000) = toulouse_Δiw(bath.spectrum, β=bath.β, nmax=nmax)

# function toulouse_Δτ(spectrum::AbstractSpectrumFunction; β::Real, N::Int)
#     f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
#     δτ = β / N
#     ff(τ) = quadgkwrapper(bounede(ϵ -> -f(ϵ) / (exp(-ϵ*τ) / (1+exp(-β*ϵ)) ), lb, ub))
#     return [ff(i*δτ) for i in 0:N]
# end

function toulouse_Δτ(f::AbstractSpectrumFunction; β::Real, N::Int)
    δτ = β / N
    function ff(τ)
        g(ϵ) = -exp(-ϵ*τ) / (1+exp(-β*ϵ)) 
        return quadgkwrapper(f / g)
    end 
    return [ff(i*δτ) for i in 0:N]
end
toulouse_Δτ(bath::AbstractFermionicBath; N::Int) = toulouse_Δτ(bath.spectrum, β=bath.β, N=N)
