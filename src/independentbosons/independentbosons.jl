_compute_Δ(s::SpectrumFunction) = quadgkwrapper(bounded(ω->s.f(ω)/ω, lb=lowerbound(s), ub=upperbound(s)))

function _exponent_f(spectrum::SpectrumFunction, τ, β)
    f = spectrum.f
    function ff(ω) 
        x = exp(-β*ω)
        (f(ω)/ω^2) * (1-exp(-τ*ω)+x-exp(-(β-τ)*ω))/(1-x)
    end
    _e = quadgkwrapper(bounded(ω -> ff(ω), lb=lowerbound(spectrum), ub=upperbound(spectrum)))
    return exp(-_e)
end

# _interact(τ, β, μ, Δ, U) = (exp(τ*(μ+Δ)) + exp(β*(μ+Δ) + τ*(μ+Δ-U))) / (1+2*exp(β*(μ+Δ))+exp(β*(2μ+2Δ-U)))

function _interact(τ, β, μ, Δ, U)
    # we have used μ ← μ + Δ
    μ′ = μ + Δ
    x = exp(-β * μ′)
    return (exp(τ*μ′) * x + exp(τ*(μ+3Δ-U))) / (x + 2 + exp(β*(μ+3Δ-U)))
end

include("realtime.jl")
include("imagtime.jl")