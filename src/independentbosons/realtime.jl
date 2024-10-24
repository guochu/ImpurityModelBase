"""
    free_holstein_greater(spectrum::SpectrumFunction, t::Real; ϵ_d::Real, n₀, Δ)

Retarded Green's function in the real time axis
"""
function independentbosons_greater(spectrum::SpectrumFunction, t::Real; β::Real, ϵ_d::Real, U::Real=0, 
                                    n₀::Union{Real, Nothing}=nothing, Δ::Real=_compute_Δ(spectrum))
    μ′ = -ϵ_d + Δ
    r::ComplexF64 = 0
    if U == 0
        if isnothing(n₀)
            n₀ = freefermion_occupation(β, μ′)
        end
        r = -im*(1-n₀)*exp(im*t*μ′)*_exponent_f(spectrum, im*t, β)
    else
        isnothing(n₀) || println("only local thermal initial state is supported for interacting case, n₀ is ignored")
        r = -im*_interact(im*t, β, μ′, U)*_exponent_f(spectrum, im*t, β)
    end
    return r
end


function independentbosons_lesser(spectrum::SpectrumFunction, t::Real; β::Real, ϵ_d::Real, U::Real=0,
                                   n₀::Union{Real, Nothing}=nothing, Δ::Real=_compute_Δ(spectrum))
    μ′ = -ϵ_d + Δ
    r::ComplexF64 = 0
    if U == 0
        if isnothing(n₀)
            n₀ = freefermion_occupation(β, μ′)
        end
        r = im*n₀*exp(im*t*μ′)*_exponent_f(spectrum, im*t, β)
    else  
        isnothing(n₀) || println("only local thermal initial state is supported for interacting case, n₀ is ignored")
        r = im*_interact(im*t, β, μ′, U)*_exponent_f(spectrum, im*t, β)
    end
    return r 
end
