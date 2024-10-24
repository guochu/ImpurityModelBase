function _compute_Δ(spectrum::SpectrumFunction)
    f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
    return quadgk(ω -> f(ω) / ω, lb, ub)[1]
end

function _exponent_f(spectrum::SpectrumFunction, τ, β)
    f, lb, ub = spectrum.f, lowerbound(spectrum), upperbound(spectrum)
    function ff(ω) 
        x = exp(-β*ω)
        (f(ω)/ω^2) * (1-exp(-τ*ω)+x-exp(-(β-τ)*ω))/(1-x)
    end
    _e = quadgk(ω -> ff(ω), lb, ub)[1]
    return exp(-_e)
end

# _interact(τ, β, μ, Δ, U) = (exp(τ*(μ+Δ)) + exp(β*(μ+Δ) + τ*(μ+Δ-U))) / (1+2*exp(β*(μ+Δ))+exp(β*(2μ+2Δ-U)))

function _interact(τ, β, μ, U)
    # we have used μ ← μ + Δ
    x = exp(-β * μ)
    return (exp(τ*μ) * x + exp(τ*(μ-U))) / (x + 2 + exp(β*(μ-U)))
end

include("realtime.jl")
include("imagtime.jl")