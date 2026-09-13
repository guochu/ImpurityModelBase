include("zeroT.jl")
include("finiteT.jl")

"""
	holstein_G0w_to_Gw(G0w::Function, ϵ::Real; g::Real, ω::Real, β::Real=Inf, maxiter::Int=10, rtol=holstein_finiteT_rtol)

Impurity Green's function of the Holstein model in the continued-fraction expansion
(CFE), given the free (bare) propagator `G0w` as a function of the complex frequency.
At zero temperature (`β = Inf`, default) the CFE is exact; at finite temperature the
implementation treats the phonon-induced spectral redistribution perturbatively.
`g` is the electron-phonon coupling and `ω` the Einstein phonon frequency.
"""
function holstein_G0w_to_Gw(G0w::Function, ϵ::Real; g::Real, ω::Real, β::Real=Inf, maxiter::Int=10, rtol::Real=holstein_finiteT_rtol)
	if β == Inf
		return holstein_G0w_to_Gw_zeroT(G0w, ϵ, g=g, ω=ω, maxiter=maxiter)
	else
		return holstein_G0w_to_Gw_finiteT(G0w, ϵ, β=β, g=g, ω=ω, maxiter=maxiter, rtol=rtol)
	end
end


"""
	holstein_G0w_to_Σw(G0w::Function, ϵ::Real; g::Real, ω::Real, β::Real=Inf, maxiter::Int=10)

Impurity self-energy of the Holstein model from the CFE, given the free (bare)
propagator `G0w`. Only implemented at zero temperature (`β = Inf`, default).
"""
function holstein_G0w_to_Σw(G0w::Function, ϵ::Real; g::Real, ω::Real, β::Real=Inf, maxiter::Int=10)
	β == Inf || error("holstein_G0w_to_Σw is only implemented at zero temperature; use holstein_G0w_to_Gw_finiteT for β < Inf")
	return holstein_G0w_to_Σw_zeroT(G0w, ϵ, g=g, ω=ω, maxiter=maxiter)
end

"""
	holstein_Gt(f::AbstractBoundedFunction, t::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, β::Real=Inf, wmax::Real=20, wmin::Real=-wmax, δ::Real=1e-8, maxiter::Int=10, rtol=holstein_finiteT_rtol)

Real-time retarded impurity Green's function of the Holstein model for the bounded
bath spectrum `f`, obtained by Fourier transforming `holstein_Gw` over
`[wmin, wmax]` (the free-particle tail `1/(ω+iδ)` is subtracted and its `-i`
contribution added back analytically). `g` is the electron-phonon coupling,
`ω` the Einstein phonon frequency, `ϵ_d`/`μ` the impurity level and chemical potential.
"""
function holstein_Gt(f::AbstractBoundedFunction, t::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0,
						β::Real=Inf, wmax::Real=20., wmin::Real=-wmax, δ::Real=1.0e-8, maxiter::Int=10, rtol::Real=holstein_finiteT_rtol)
    A = quadgkwrapper(bounded(ϵ -> (holstein_Gw(f, ϵ; β=β, ϵ_d=ϵ_d, g=g, ω=ω, μ=μ, δ=δ, maxiter=maxiter, rtol=rtol)-1.0/(ϵ+im*δ))*exp(-im*ϵ*t), wmin, wmax))
    return A/(2π)-im
end

"""
	holstein_Gw(f::AbstractBoundedFunction, ϵ::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, β::Real=Inf, δ::Real=1e-8, maxiter::Int=10, rtol=holstein_finiteT_rtol)

Retarded impurity Green's function of the Holstein model at real frequency `ϵ` for
the bounded bath spectrum `f`: the free propagator is
``G_0(y) = \\int \\mathrm{d}\\epsilon'\\, f(\\epsilon')/(y + \\mu - \\epsilon' + i\\delta)``
(computed with the Toulouse-model formula) and the CFE maps `G_0` to the dressed
impurity Green's function. See [`holstein_G0w_to_Gw`](@ref).
"""
function holstein_Gw(f::AbstractBoundedFunction, ϵ::Real; g::Real, ω::Real, ϵ_d::Real, μ::Real=0, β::Real=Inf,
						δ::Real=1.0e-8, maxiter::Int=10, rtol::Real=holstein_finiteT_rtol)
	G0w(y) = fermionic_toulouse_Gw(f, y, ϵ_d=ϵ_d, μ=μ, δ=δ)
	return holstein_G0w_to_Gw(G0w, ϵ, g=g, ω=ω, β=β, maxiter=maxiter, rtol=rtol)
end