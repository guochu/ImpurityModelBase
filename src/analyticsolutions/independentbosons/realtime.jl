"""
    independentbosons_greater(spectrum::AbstractSpectrumFunction, t; β, ϵ_d, U, bands, Δ)

Greater Green's function in the real time axis for the independent bosons model
"""
function independentbosons_greater(spectrum::AbstractSpectrumFunction, t::Real; β::Real, ϵ_d::Real, U::Real=0, 
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
    independentbosons_lesser(spectrum::AbstractSpectrumFunction, t; β, ϵ_d, U, bands, Δ)

Lesser Green's function in the real time axis for the independent bosons model
"""
function independentbosons_lesser(spectrum::AbstractSpectrumFunction, t::Real; β::Real, ϵ_d::Real, U::Real=0,
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