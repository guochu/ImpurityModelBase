# DMFT solution of the Holstein model on the Bethe lattice
# following Ciuchi, de Pasquale, Fratini & Feinberg, PRB 56, 4494 (1997)

"""
    bethe_holstein_G0inv_from_G(z, G; t::Real=1)

Bethe-lattice self-consistency condition (paper Eq. 37):
G0⁻¹(ω) = ω - (t²/4) G(ω).

`z = ω + iδ` is the retarded complex frequency, `G` the local (site) Green's
function and `t` the half bandwidth of the Bethe lattice.
"""
bethe_holstein_G0inv_from_G(z, G; t::Real=1) = z - (t^2/4) * G

"""
    bethe_holstein_G_from_G0inv(g0invfun::Function, ω::Real; g, ω0, β, maxiter)

Dyson equation for the local Green's function:
G(ω) = 1 / (G0⁻¹(ω) - Σ(ω)),
with the self-energy Σ(ω) obtained from the continued-fraction expansion
(CFE) of the small polaron impurity model (paper Eqs. 35/36 at T=0 and
Eqs. 39-42 at finite T), expressed as a functional of `g0invfun`.
"""
function bethe_holstein_G_from_G0inv(g0invfun::Function, ω::Real; g::Real, ω0::Real, β::Real=Inf, maxiter::Int=10)
    G0w(y) = 1 / g0invfun(y)
    return holstein_G0w_to_Gw(G0w, ω; g=g, ω=ω0, β=β, maxiter=maxiter)
end

"""
    holstein_dmft_bethe(ws; g, ω, [λ, γ], t=1, β=Inf, δ, maxiter, mixing, tol, maxit, nw, init)

Solve the DMFT self-consistency of the Holstein model on an infinite-coordination
Bethe lattice (half bandwidth `t`), at zero (`β=Inf`) or finite temperature.

The impurity problem (single site coupled to Einstein phonons and embedded in a
non-interacting effective medium) is solved by the continued-fraction expansion
(CFE) of Ref. [Ciuchi et al., PRB 56, 4494 (1997)]. The effective medium is
fixed self-consistently by the Bethe-lattice condition
    G0⁻¹(ω) = ω - (t²/4) G(ω).                                        (Eq. 37)
Both the free propagator G0 and the local propagator G are retarded: the
calculation is performed on the grid `z = ω + iδ`.

The electron-phonon coupling can be specified either with the bare parameters
(`g`, `ω`), or with the scaleless ones (`λ`, `γ`) where λ = g²/(ω t) and
γ = ω/t (see `holstein_scaleless_parameters`).

# Arguments
- `ws`: frequency window (in units of `t`) on which the spectrum is requested.
- `g, ω`: bare electron-phonon coupling and Einstein phonon frequency.
- `λ, γ`: scaleless coupling λ = g²/(ω t) and adiabatic ratio γ = ω/t.
- `t=1`: half bandwidth of the Bethe lattice.
- `β=Inf`: inverse temperature (Inf = zero temperature).
- `δ=1e-3`: broadening of the retarded Green's function.
- `maxiter`: truncation order of the CFE; default ≈ max(10, 5α²+5) with α²=λ/γ.
- `mixing=0.5`: linear mixing parameter for the self-consistency iteration.
- `tol=1e-6`: convergence threshold on max |ΔA(ω)|.
- `maxit=1000`: maximum number of DMFT iterations.
- `nw=2000`: number of grid points used internally (the requested window is
  padded to cover all frequency shifts ω ± nω0 needed by the CFE).
- `init=:bethe`: initial guess for G0⁻¹, `:atomic` (G0⁻¹=ω+iδ) or `:bethe`
  (non-interacting Bethe-lattice solution).

# Returns
A `NamedTuple` with fields `ws`, `A` (spectral density -Im G/π), `G`, `ImΣ`,
`ReΣ`, `converged` and `iterations`.
"""
function holstein_dmft_bethe(ws::AbstractVector{<:Real};
        g::Union{Nothing,Real}=nothing, ω::Union{Nothing,Real}=nothing,
        λ::Union{Nothing,Real}=nothing, γ::Union{Nothing,Real}=nothing,
        t::Real=1, β::Real=Inf,
        δ::Real=1.0e-3, maxiter::Union{Nothing,Int}=nothing,
        mixing::Real=0.5, tol::Real=1.0e-6, maxit::Int=1000,
        nw::Int=2000, init::Symbol=:bethe, verbose::Bool=false)

    # ---- parameter conversion (bare or scaleless) ----
    if (g === nothing) != (ω === nothing)
        error("give both bare parameters g and ω, or both scaleless λ and γ")
    end
    if g === nothing
        (λ === nothing || γ === nothing) && error("give λ and γ, or g and ω")
        ω, g = holstein_bare_parameters(λ=λ, γ=γ, t=t)
    end
    λ = g^2 / (ω * t)
    γ = ω / t
    α2 = λ / γ                       # α² = g²/ω²

    maxiter === nothing && (maxiter = max(10, ceil(Int, 5α2) + 5))

    # ---- frequency grid, padded to cover the CFE shifts ω ± n ω0 ----
    wmin, wmax = extrema(ws)
    ntherm = (β == Inf) ? 0 : max(5, ceil(Int, 10 / (β * ω)))
    pad = (maxiter + ntherm) * ω + 2δ
    wgrid = range(wmin - pad, wmax + pad; length=nw)
    z = collect(wgrid) .+ im * δ

    # ---- initial guess for G0⁻¹ ----
    g0inv = if init === :atomic
        copy(z)
    else
        # non-interacting Bethe-lattice solution: G0⁻¹ = (z + sqrt(z²-t²))/2,
        # keeping the retarded branch (Im G0⁻¹ > 0)
        s = sqrt.(z .^ 2 .- t^2)
        g0inv0 = (z .+ s) ./ 2
        @. ifelse(imag(g0inv0) < 0, (z - s) / 2, g0inv0)
    end

    interp = linear_interpolation(wgrid, g0inv; extrapolation_bc=Line())
    g0invfun(y::Real) = (wmin - pad ≤ y ≤ wmax + pad) ? interp(y) : (y + im * δ)

    # ---- DMFT self-consistency loop ----
    A_prev = nothing
    converged = false
    iterations = 0
    G = Vector{ComplexF64}(undef, nw)
    for it in 1:maxit
        iterations = it
        for k in eachindex(wgrid)
            G[k] = bethe_holstein_G_from_G0inv(g0invfun, wgrid[k]; g=g, ω0=ω, β=β, maxiter=maxiter)
        end
        A = -imag.(G) ./ π
        g0inv_new = bethe_holstein_G0inv_from_G(z, G; t=t)
        g0inv .= (1 - mixing) .* g0inv .+ mixing .* g0inv_new
        interp = linear_interpolation(wgrid, g0inv; extrapolation_bc=Line())
        if A_prev !== nothing && maximum(abs.(A .- A_prev)) < tol
            converged = true
            verbose && println("DMFT converged in $it iterations")
            break
        end
        A_prev = copy(A)
    end
    converged || verbose && println("DMFT did not converge within $maxit iterations")

    # ---- final evaluation at the fixed-point G0⁻¹ ----
    for k in eachindex(wgrid)
        G[k] = bethe_holstein_G_from_G0inv(g0invfun, wgrid[k]; g=g, ω0=ω, β=β, maxiter=maxiter)
    end
    # self-energy at the fixed point: Σ = G0⁻¹ - G⁻¹, with G0⁻¹ = z - (t²/4)G
    Σ = z .- (t^2 / 4) .* G .- 1 ./ G

    # ---- interpolate the padded-grid result back onto the requested `ws` ----
    A_fine = -imag.(G) ./ π
    mask = (wmin .≤ wgrid .≤ wmax)
    if all(w -> wmin ≤ w ≤ wmax, ws)
        ia = linear_interpolation(wgrid, A_fine; extrapolation_bc=Line())
        iG = linear_interpolation(wgrid, G; extrapolation_bc=Line())
        iΣ = linear_interpolation(wgrid, Σ; extrapolation_bc=Line())
        return (ws=collect(ws),
                A=ia.(ws),
                G=iG.(ws),
                ImΣ=imag.(iΣ.(ws)),
                ReΣ=real.(iΣ.(ws)),
                converged=converged,
                iterations=iterations)
    end
    return (ws=collect(wgrid[mask]),
            A=A_fine[mask],
            G=G[mask],
            ImΣ=imag.(Σ)[mask],
            ReΣ=real.(Σ)[mask],
            converged=converged,
            iterations=iterations)
end
