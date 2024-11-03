"""
    free_holstein_greater(spectrum::AbstractSpectrumFunction, t::Real; ϵ_d::Real, n₀, Δ)

Retarded Green's function in the real time axis
"""
function independentbosons_greater(spectrum::AbstractSpectrumFunction, t::Real; β::Real, ϵ_d::Real, U::Real=0, 
                                    n₀::Union{Real, Nothing}=nothing, nbands::Int=1, Δ::Real=_compute_Δ(spectrum))
    (nbands in (1, 2)) || throw(ArgumentError("nbands must be 1 or 2"))
    μ′ = -ϵ_d + Δ
    r::ComplexF64 = 0
    if nbands == 1
        (U == 0) || println("nonzero U=$(U) ignored for nbands=1")
        if isnothing(n₀)
            n₀ = freefermion_occupation(β, μ′)
        end
        r = -im*(1-n₀)*exp(im*t*μ′)*_exponent_f(spectrum, im*t, β)
    else
        isnothing(n₀) || println("only local thermal initial state is supported for interacting case, n₀ is ignored")
        r = -im*_interact(im*t, β, -ϵ_d, Δ, U)*_exponent_f(spectrum, im*t, β)
    end
    return r
end


function independentbosons_lesser(spectrum::AbstractSpectrumFunction, t::Real; β::Real, ϵ_d::Real, U::Real=0,
                                   n₀::Union{Real, Nothing}=nothing, nbands::Int=1, Δ::Real=_compute_Δ(spectrum))
    (nbands in (1, 2)) || throw(ArgumentError("nbands must be 1 or 2"))
    μ′ = -ϵ_d + Δ
    r::ComplexF64 = 0
    if nbands == 1
        (U == 0) || println("nonzero U=$(U) ignored for nbands=1")
        if isnothing(n₀)
            n₀ = freefermion_occupation(β, μ′)
        end
        r = im*n₀*exp(im*t*μ′)*_exponent_f(spectrum, -im*t, β)
    else               
        isnothing(n₀) || println("only local thermal initial state is supported for interacting case, n₀ is ignored")
        r = im*_interact_lesser(im*t, β, -ϵ_d, Δ, U)*_exponent_f(spectrum, -im*t, β)        
    end
    return r 
end



function _interact_lesser(τ, β, μ, Δ, U)
    # we have used μ ← μ + Δ
    μ′ = μ + Δ
    x = exp(-β * μ′)
    y = exp(β*(μ+3Δ-U))
    return (exp(τ*μ′) + y*exp(τ*(μ+3Δ-U))) / (x + 2 + y)
end