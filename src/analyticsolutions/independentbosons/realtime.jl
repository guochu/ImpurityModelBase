"""
    independentbosons_greater(spectrum::AbstractBoundedFunction, t; β, ϵ_d, U, bands, Δ)

Greater Green's function in the real time axis for the independent bosons model
"""
function independentbosons_greater(spectrum::AbstractBoundedFunction, t::Real; β::Real, ϵ_d::Real, U::Real=0, 
                                    bands::Int=1, Δ::Real=_compute_Δ(spectrum))
    μ, U = _normalized_paras(-ϵ_d, U, Δ, bands)
    r::ComplexF64 = 0
    if bands == 1
        r = freefermion_greater(t, β=β, μ=μ) * _exponent_f(spectrum, im*t, β)
    else
        r = fermion_greater(t, β=β, μ=μ, U=U) * _exponent_f(spectrum, im*t, β)
    end
    return r
end

"""
    independentbosons_lesser(spectrum::AbstractBoundedFunction, t; β, ϵ_d, U, bands, Δ)

Lesser Green's function in the real time axis for the independent bosons model
"""
function independentbosons_lesser(spectrum::AbstractBoundedFunction, t::Real; β::Real, ϵ_d::Real, U::Real=0,
                                    bands::Int=1, Δ::Real=_compute_Δ(spectrum))
    μ, U = _normalized_paras(-ϵ_d, U, Δ, bands)
    r::ComplexF64 = 0
    if bands == 1
        r = freefermion_lesser(t, β=β, μ=μ) * _exponent_f(spectrum, -im*t, β)
    else               
        r = fermion_lesser(t, β=β, μ=μ, U=U) * _exponent_f(spectrum, -im*t, β)        
    end
    return r 
end

"""
    independentbosons_greater(spectrum::AbstractBoundedFunction, t, ρ_0; β, ϵ_d, U, bands, Δ)

Greater Green's function in the real time axis for the independent bosons model with an
arbitrary impurity initial density matrix `ρ_0` (the bosonic bath is still taken to be in
thermal equilibrium `ρ_bath = e^{-βH_b}/Z_b`).

`ρ_0` is a full impurity density matrix:
- `bands = 1`: 2×2 matrix in the `[|0⟩, |1⟩]` basis;
- `bands = 2`: 4×4 matrix in the `[|0⟩, |↑⟩, |↓⟩, |↑↓⟩]` basis.

Because `[H, n̂_d] = 0`, the impurity occupation is conserved and only the *diagonal*
elements of `ρ_0` enter the Green's function; coherences (off-diagonal elements) do not
contribute. For `bands = 2` the transition into the doubly-occupied state carries the
extra bath phase `exp(-2i∫dω J(ω)/ω² sin(ωt))` (`_double_transition`).
"""
function independentbosons_greater(spectrum::AbstractBoundedFunction, t::Real, ρ_0::AbstractMatrix{<:Number};
                                    β::Real, ϵ_d::Real, U::Real=0, bands::Int=1, Δ::Real=_compute_Δ(spectrum))
    d = _ρ₀_diag(ρ_0, bands)
    μ, U = _normalized_paras(-ϵ_d, U, Δ, bands)
    r::ComplexF64 = 0
    if bands == 1
        r = -im * d[1] * exp(im*μ*t) * _exponent_f(spectrum, im*t, β)
    else
        r = -im * (d[1] * exp(im*μ*t) + d[3] * exp(im*(μ-U)*t) * _double_transition(spectrum, t)) *
            _exponent_f(spectrum, im*t, β)
    end
    return r
end

"""
    independentbosons_lesser(spectrum::AbstractBoundedFunction, t, ρ_0; β, ϵ_d, U, bands, Δ)

Lesser Green's function in the real time axis for the independent bosons model with an
arbitrary impurity initial density matrix `ρ_0` (the bosonic bath is still taken to be in
thermal equilibrium `ρ_bath = e^{-βH_b}/Z_b`).

`ρ_0` is a full impurity density matrix:
- `bands = 1`: 2×2 matrix in the `[|0⟩, |1⟩]` basis;
- `bands = 2`: 4×4 matrix in the `[|0⟩, |↑⟩, |↓⟩, |↑↓⟩]` basis.

Only the diagonal elements of `ρ_0` enter (impurity occupation is conserved). For a
product initial state the lesser Green's function carries the same polaron factor
`exp(-Φ(it))` as the greater one. For `bands = 2` the transition leaving the
doubly-occupied state carries the extra bath phase `exp(-2i∫dω J(ω)/ω² sin(ωt))`
(`_double_transition`).
"""
function independentbosons_lesser(spectrum::AbstractBoundedFunction, t::Real, ρ_0::AbstractMatrix{<:Number};
                                   β::Real, ϵ_d::Real, U::Real=0, bands::Int=1, Δ::Real=_compute_Δ(spectrum))
    d = _ρ₀_diag(ρ_0, bands)
    μ, U = _normalized_paras(-ϵ_d, U, Δ, bands)
    r::ComplexF64 = 0
    if bands == 1
        r = im * d[2] * exp(im*μ*t) * _exponent_f(spectrum, im*t, β)
    else
        r = im * (d[2] * exp(im*μ*t) + d[4] * exp(im*(μ-U)*t) * _double_transition(spectrum, t)) *
            _exponent_f(spectrum, im*t, β)
    end
    return r
end

# extract the diagonal probabilities of a full impurity density matrix, validating
# the dimension against `bands`
function _ρ₀_diag(ρ_0::AbstractMatrix{<:Number}, bands::Int)
    n = size(ρ_0, 1)
    expected = bands == 1 ? 2 : 4
    (n == expected && size(ρ_0, 2) == expected) ||
        throw(ArgumentError("ρ_0 must be $(expected)x$(expected) for bands=$bands"))
    return real.(diag(ρ_0))
end



# function _interact_lesser(τ, β, μ, Δ, U)
#     # we have used μ ← μ + Δ
#     μ′ = μ + Δ
#     x = exp(-β * μ′)
#     y = exp(β*(μ+3Δ-U))
#     return (exp(τ*μ′) + y*exp(τ*(μ+3Δ-U))) / (x + 2 + y)
# end

function _normalized_paras(μ, U, Δ, bands::Int)
    (bands in (1, 2)) || throw(ArgumentError("bands must be 1 or 2"))
    μ = μ + Δ
    if bands == 2
        U = U - 2Δ
    else
        (U == 0) || println("nonzero U=$(U) ignored for bands=1")
    end
    return μ, U   
end