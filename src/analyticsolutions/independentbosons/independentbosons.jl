function _compute_Δ(s::AbstractBoundedFunction)
    g(ω) = 1/ω
    return quadgkwrapper(s * g)
end

function _exponent_f(f::AbstractBoundedFunction, τ, β)
    function g(ω)
        x = exp(-β*ω)
        return (1-exp(-τ*ω)+x-exp(-(β-τ)*ω))/((1-x)*ω^2)
    end
    _e = quadgkwrapper(f * g)
    return exp(-_e)
end

# Bath transition factor for entering/leaving the doubly-occupied impurity
# state (bands=2) with a product initial state:
# exp(-2i ∫ dω J(ω)/ω² sin(ωt)).
# It encodes the extra displacement phase picked up when the electron is added
# on top of an already-occupied (opposite-spin) impurity state.
function _double_transition(f::AbstractBoundedFunction, t::Real)
    g(ω) = sin(ω*t)/ω^2
    _e = quadgkwrapper(f * g)
    return exp(-2im * _e)
end

# _interact(τ, β, μ, Δ, U) = (exp(τ*(μ+Δ)) + exp(β*(μ+Δ) + τ*(μ+Δ-U))) / (1+2*exp(β*(μ+Δ))+exp(β*(2μ+2Δ-U)))

# function _interact(τ, β, μ, Δ, U)
#     # we have used μ ← μ + Δ
#     μ′ = μ + Δ
#     x = exp(-β * μ′)
#     return (exp(τ*μ′) * x + exp(τ*(μ+3Δ-U))) / (x + 2 + exp(β*(μ+3Δ-U)))
# end
# function _interact(τ, β, μ, U)
#     # we have used μ ← μ + Δ
#     x = exp(-β * μ)
#     return (exp(τ*μ) * x + exp(τ*(μ-U))) / (x + 2 + exp(β*(μ-U)))
# end

include("realtime.jl")
include("imagtime.jl")