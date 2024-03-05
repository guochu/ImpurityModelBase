# real-time Green's function of Toulouse model
function toulouse_Gw(spectrum::SpectrumFunction, ω::Float64; ϵ_d::Real, μ::Real=0)
	δ = 1e-8
	f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
	1.0/(ω-ϵ_d-quadgk(ε -> f(ε)/(ω+μ-ε+im*δ), lb, ub)[1])
end

function toulouse_Gt(spectrum::SpectrumFunction, t::Float64; ϵ_d::Real, μ::Real=0, wmax::Real=20.)
    δ = 1e-8
    A = quadgk(ω -> (toulouse_Gw(spectrum, ω; ϵ_d=ϵ_d, μ=μ)-1.0/(ω+im*δ))*exp(-im*ω*t), -wmax, wmax)[1]
    im*(A/(2π)-im)
end


# imaginary-time Green's function of Toulouse model
function toulouse_Giw(spectrum::SpectrumFunction, ω::Float64; ϵ_d::Real, μ::Real=0)
	f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
    1.0/(im*ω-ϵ_d-quadgk(ε -> f(ε)/(im*ω+μ-ε), lb, ub)[1])
end

function toulouse_Gτ(spectrum::SpectrumFunction, τ::Float64; β::Real, ϵ_d::Real, μ::Real=0., nmax::Int=1000)
    res = 0.0
    for n = -nmax:nmax+1
        ω = (2n-1)*π/β
        res += (toulouse_Giw(spectrum, ω; ϵ_d=ϵ_d, μ=μ)-1/(im*ω))*exp(-im*τ*ω)
    end
    res = -(res/β-0.5)
end