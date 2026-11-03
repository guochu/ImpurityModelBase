# Contour-ordered single-particle Green's functions of free (quadratic) systems,
# discretized on the closed Keldysh contour (RealGreenFunctions) and on the
# L-shaped Kadanoff-Baym contour (MixedGreenFunctions), following the block
# structure of QuAPI's RealCorrelationFunction / MixedCorrelationFunction:
# only the independent cross blocks are stored and the diagonal blocks are
# derived on demand through `branch`.


"""
    RealGreenFunctions

Contour-ordered single-particle Green's functions of a free-particle system on
the closed real-time (Keldysh) contour, between a fixed pair of sites (i, j).
The statistics enters only through the constructions
`freefermions_keldysh_gfs` / `freebosons_keldysh_gfs`, which fill the blocks
with the fermionic/bosonic conventions described below; the fields are plain
dense matrices.

Following the minimal-block storage of QuAPI's `RealCorrelationFunction`, only
the two independent CROSS blocks are stored (both `Matrix{ComplexF64}` of size
(Nt+1)×(Nt+1) over the real-time nodes t₀ = 0, …, t_{Nt}):

- `G₊₋`:  the (+,−) block, the lesser Green's function,
          `G₊₋[a, b] = G^<(tₐ − t_b)`;
- `G₋₊`:  the (−,+) block, the greater Green's function,
          `G₋₊[a, b] = G^>(tₐ − t_b)`.

With the Green's-function convention `G^<(t) = i⟨x†(0)x(t)⟩`-type kernel of the
lesser weight and `G^>(t) = -i⟨x(t)x†(0)⟩`-type kernel of the greater weight
(the statistics-dependent weights being n̄ for the lesser and 1−n̄ (fermions) /
1+n̄ (bosons) for the greater kernel), the diagonal blocks are derived exactly
by `branch(x; b1, b2)` through the contour-ordering sector rule:

    (+,+)[a, b] = a > b ? G₋₊[a, b] : G₊₋[a, b]   (time ordered)
    (−,−)[a, b] = a > b ? G₊₋[a, b] : G₋₊[a, b]   (anti-time ordered)

so the equal-time contacts take the lesser partner for (+,+) and the greater
partner for (−,−), the same contact convention as QuAPI's derived diagonal
blocks.

Construct with `freefermions_keldysh_gfs` / `freebosons_keldysh_gfs`.
"""
struct RealGreenFunctions
    G₊₋::Matrix{ComplexF64}
    G₋₊::Matrix{ComplexF64}
end


"""
    MixedGreenFunctions

Contour-ordered Green's functions of a free-particle system in thermal
equilibrium on the L-shaped (Kadanoff-Baym) contour between a fixed pair of
sites (i, j): an imaginary-time branch `:τ` of nodes 0 ≤ τ < β followed by the
closed real-time Keldysh part with branches `:+`/`:-`.  The statistics enters
only through the constructions `freefermions_kadanoff_gfs` /
`freebosons_kadanoff_gfs`; the fields are plain dense matrices.

Following the minimal-block storage of QuAPI's `MixedCorrelationFunction`, the
fields are the two real-time cross blocks (`G₊₋`, `G₋₊`, of size
(Nt+1)×(Nt+1), see `RealGreenFunctions`) plus the three mixed blocks

- `G₊τ`:  the (+,τ) block, `G₊τ[a, m] = G(tₐ, τₘ)`, (Nt+1)×Nτ;
- `Gτ₊`:  the (τ,+) block, `Gτ₊[m, a] = G(τₘ, tₐ)`, Nτ×(Nt+1);
- `Gττ`:  the (τ,τ) block, the Matsubara Green's function `+i·Gᴹ(τₘ − τₙ)`
  with `Gᴹ(τ) = -⟨T̂ x(τ) x†(0)⟩`, Nτ×Nτ, equal-time contact taken in the
  τ→0⁺ (greater-weight) convention.

`branch(x; b1, b2)` with `b1, b2 ∈ (:+, :-, :τ)` returns the nine contour
sub-blocks; the backward-branch combinations are exact copies of the forward
ones, `branch(:−, :τ) == branch(:+, :τ)` and `branch(:τ, :−) == branch(:τ, :+)`,
since the relative contour ordering of the imaginary-time and real-time
operators does not depend on the Keldysh branch.

Construct with `freefermions_kadanoff_gfs` / `freebosons_kadanoff_gfs`.
"""
struct MixedGreenFunctions
    G₊₋::Matrix{ComplexF64}
    G₋₊::Matrix{ComplexF64}
    G₊τ::Matrix{ComplexF64}
    Gτ₊::Matrix{ComplexF64}
    Gττ::Matrix{ComplexF64}
end

"""
    ImagGreenFunctions{T<:Number}

Matsubara (imaginary-time) single-particle Green's functions of a free-particle
system in thermal equilibrium between a fixed pair of sites (i, j), evaluated on
an imaginary-time grid 0 ≤ τ ≤ β.  The single field `data::Matrix{T}` stores

    data[m, n] = -i⟨T̂ x(τₘ) x†(τₙ)⟩ = +i·Gᴹ(τₘ − τₙ),

the (τ,τ) block convention of `MixedGreenFunctions` (so on matching grids
`branch(free*_kadanoff_gfs(h, ts, τs, ...), :τ, :τ) == free*_matsubara_gfs(h, τs, ...).data`),
with `Gᴹ(τ) = -⟨T̂ x(τ) x†(0)⟩` the standard Matsubara Green's function and the
equal-time contact taken in the τ→0⁺ (greater-weight) convention.
`branch(x; b1, b2)` accesses the data with the only valid labels `b1 = b2 = :τ`.

Construct with `freefermions_matsubara_gfs` / `freebosons_matsubara_gfs`.
"""
struct ImagGreenFunctions{T<:Number}
    data::Matrix{T}
end
ImagGreenFunctions(data::AbstractMatrix) = ImagGreenFunctions(Matrix{ComplexF64}(data))


Base.show(io::IO, ::MIME"text/plain", x::RealGreenFunctions) =
    print(io, "real-time Green's functions [$(size(x.G₊₋, 1))] (2 cross blocks)")
Base.show(io::IO, ::MIME"text/plain", x::MixedGreenFunctions) =
    print(io, "Kadanoff Green's functions [$(size(x.Gττ, 1))+$(size(x.G₊₋, 1))]")
Base.show(io::IO, ::MIME"text/plain", x::ImagGreenFunctions) =
    print(io, "Matsubara Green's functions [$(size(x.data, 1))]")

Base.size(x::RealGreenFunctions, i::Int...) = size(x.G₊₋, i...)
Base.size(x::MixedGreenFunctions, i::Int...) = i === () ? (size(x.Gττ, 1), size(x.G₊₋, 1)) : size(x.G₊₋, i...)
Base.size(x::ImagGreenFunctions, i::Int...) = size(x.data, i...)

Base.:(==)(x::RealGreenFunctions, y::RealGreenFunctions) = (x.G₊₋ == y.G₊₋) && (x.G₋₊ == y.G₋₊)
Base.:(==)(x::MixedGreenFunctions, y::MixedGreenFunctions) =
    (x.G₊₋ == y.G₊₋) && (x.G₋₊ == y.G₋₊) && (x.G₊τ == y.G₊τ) && (x.Gτ₊ == y.Gτ₊) && (x.Gττ == y.Gττ)
Base.:(==)(x::ImagGreenFunctions, y::ImagGreenFunctions) = x.data == y.data
Base.isapprox(x::RealGreenFunctions, y::RealGreenFunctions; kwargs...) =
    isapprox(x.G₊₋, y.G₊₋; kwargs...) && isapprox(x.G₋₊, y.G₋₊; kwargs...)
Base.isapprox(x::MixedGreenFunctions, y::MixedGreenFunctions; kwargs...) =
    isapprox(x.G₊₋, y.G₊₋; kwargs...) && isapprox(x.G₋₊, y.G₋₊; kwargs...) &&
    isapprox(x.G₊τ, y.G₊τ; kwargs...) && isapprox(x.Gτ₊, y.Gτ₊; kwargs...) &&
    isapprox(x.Gττ, y.Gττ; kwargs...)
Base.isapprox(x::ImagGreenFunctions, y::ImagGreenFunctions; kwargs...) = isapprox(x.data, y.data; kwargs...)


"""
    branches(x::RealGreenFunctions) -> (:+, :-)
    branches(x::MixedGreenFunctions) -> (:+, :-, :τ)
    branches(x::ImagGreenFunctions) -> (:τ,)

Contour branches carried by the Green's functions, in contour order.
"""
branches(::RealGreenFunctions) = (:+, :-)
branches(::MixedGreenFunctions) = (:+, :-, :τ)
branches(::ImagGreenFunctions) = (:τ,)


# —— diagonal real-time blocks derived from the cross pair (branch 与 Mixed 共用) ——

# time-ordered (+,+) block: contour-ordering sector rule — the strict lower
# triangle (tₐ > t_b) from the greater partner, the rest (including the
# equal-time contact) from the lesser partner.
function _G₊₊(G₊₋::AbstractMatrix, G₋₊::AbstractMatrix)
    n = size(G₊₋, 1)
    M = similar(G₊₋)
    for a in 1:n, b in 1:n
        M[a, b] = a > b ? G₋₊[a, b] : G₊₋[a, b]
    end
    return M
end

# anti-time-ordered (−,−) block: contour-ordering sector rule — the strict
# lower triangle (tₐ > t_b) from the lesser partner, the rest (including the
# equal-time contact) from the greater partner.
function _G₋₋(G₊₋::AbstractMatrix, G₋₊::AbstractMatrix)
    n = size(G₊₋, 1)
    M = similar(G₊₋)
    for a in 1:n, b in 1:n
        M[a, b] = a > b ? G₊₋[a, b] : G₋₊[a, b]
    end
    return M
end


"""
    branch(x::RealGreenFunctions, b1::Symbol, b2::Symbol)
    branch(x::RealGreenFunctions; b1, b2)

Return the real-time contour sub-block `G[b1, b2]` as a `Matrix{ComplexF64}`,
with `b1, b2 ∈ (:+, :-)` labelling the forward/backward branch.  The stored
cross blocks are returned for the mixed branch pairs, the diagonal blocks are
derived by the contour-ordering sector rule (see `RealGreenFunctions`).
"""
function branch(x::RealGreenFunctions, b1::Symbol, b2::Symbol)
    (b1 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    (b2 in (:+, :-)) || throw(ArgumentError("branch must be :+ or :-"))
    if b1 == :+
        return b2 == :+ ? _G₊₊(x.G₊₋, x.G₋₊) : x.G₊₋
    else
        return b2 == :+ ? x.G₋₊ : _G₋₋(x.G₊₋, x.G₋₊)
    end
end

"""
    branch(x::MixedGreenFunctions, b1::Symbol, b2::Symbol)
    branch(x::MixedGreenFunctions; b1, b2)

Return the Kadanoff contour sub-block `G[b1, b2]` as a `Matrix{ComplexF64}`,
with `b1, b2 ∈ (:+, :-, :τ)`.  The real-time diagonal blocks are derived by the
contour-ordering sector rule; the imaginary-time combinations
`(−,τ) = (+,τ)` and `(τ,−) = (τ,+)` are exact copies of the forward ones.
"""
function branch(x::MixedGreenFunctions, b1::Symbol, b2::Symbol)
    (b1 in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
    (b2 in (:+, :-, :τ)) || throw(ArgumentError("branch must be one of :+, :- or :τ"))
    if b1 == :τ
        b2 == :τ && return x.Gττ
        return x.Gτ₊                       # (τ,±) share one block
    elseif b2 == :τ
        return x.G₊τ                       # (±,τ) share one block
    end
    if b1 == :+
        return b2 == :+ ? _G₊₊(x.G₊₋, x.G₋₊) : x.G₊₋
    else
        return b2 == :+ ? x.G₋₊ : _G₋₋(x.G₊₋, x.G₋₊)
    end
end

branch(x::RealGreenFunctions; b1::Symbol, b2::Symbol) = branch(x, b1, b2)
branch(x::MixedGreenFunctions; b1::Symbol, b2::Symbol) = branch(x, b1, b2)

"""
    branch(x::ImagGreenFunctions, :τ, :τ)

Return the stored Matsubara data as a `Matrix`; the only valid branch labels
are `b1 = b2 = :τ`.
"""
function branch(x::ImagGreenFunctions, b1::Symbol, b2::Symbol)
    (b1 == :τ && b2 == :τ) || throw(ArgumentError("branch must be :τ for ImagGreenFunctions"))
    return x.data
end
branch(x::ImagGreenFunctions; b1::Symbol, b2::Symbol) = branch(x, b1, b2)


# —— equilibrium spectral kernels shared by the constructors ——————————————————
# Single-particle basis of mode k: weight wₖ = U[i, k]·conj(U[j, k]); the lesser
# kernel carries the occupation n̄ₖ, the greater kernel the weight 1−n̄ₖ (fermions)
# or 1+n̄ₖ (bosons).  All blocks are -i⟨T_c x x†⟩-type contour-ordered functions;
# the statistics sign s = +1 (fermions) / -1 (bosons) multiplies the lesser-type
# kernels (whose raw expectation enters with opposite swap signs in T_c).

function _lesser_greater_blocks(::Type{P}, cache::EigenCache, ts::AbstractVector{<:Real},
                                i::Int, j::Int; β::Real, μ::Real) where {P<:AbstractParticle}
    U, λs = cache.U, cache.λs
    s = (P === Fermion) ? 1.0 : -1.0
    ns = [thermaloccupation(P, β, μ, λ) for λ in λs]
    ws = [U[i, k] * conj(U[j, k]) for k in eachindex(λs)]
    Nt = length(ts)
    G₊₋ = zeros(ComplexF64, Nt, Nt)
    G₋₊ = zeros(ComplexF64, Nt, Nt)
    for a in 1:Nt, b in 1:Nt
        Δt = ts[a] - ts[b]
        G₊₋[a, b] = s * im * sum(ws[k] * ns[k] * exp(-im * λs[k] * Δt) for k in eachindex(λs))
        G₋₊[a, b] = -im * sum(ws[k] * ((P === Fermion ? 1.0 - ns[k] : 1.0 + ns[k])) * exp(-im * λs[k] * Δt) for k in eachindex(λs))
    end
    return G₊₋, G₋₊
end

# (τ,τ) block shared by the Kadanoff contour and the standalone Matsubara
# container: +i·Gᴹ(Δτ), contact in the τ→0⁺ (greater weight) convention;
# x†(|Δτ|) = e^{+τε}x† grows in imaginary time (KMS-bounded for |Δτ| ≤ β)
function _matsubara_block(::Type{P}, cache::EigenCache, τs::AbstractVector{<:Real},
                          i::Int, j::Int; β::Real, μ::Real) where {P<:AbstractParticle}
    U, λs = cache.U, cache.λs
    s = (P === Fermion) ? 1.0 : -1.0
    ns = [thermaloccupation(P, β, μ, λ) for λ in λs]
    ws = [U[i, k] * conj(U[j, k]) for k in eachindex(λs)]
    gk(k) = (P === Fermion) ? 1.0 - ns[k] : 1.0 + ns[k]
    Nτ = length(τs)
    Gττ = zeros(ComplexF64, Nτ, Nτ)
    for m in 1:Nτ, n in 1:Nτ
        Δτ = τs[m] - τs[n]
        Gττ[m, n] = Δτ >= 0 ? -im * sum(ws[k] * gk(k) * exp(-λs[k] * Δτ) for k in eachindex(λs)) :
                              s * im * sum(ws[k] * ns[k] * exp(-λs[k] * Δτ) for k in eachindex(λs))
    end
    return Gττ
end

function _mixed_blocks(::Type{P}, cache::EigenCache, ts::AbstractVector{<:Real},
                       τs::AbstractVector{<:Real}, i::Int, j::Int; β::Real, μ::Real) where {P<:AbstractParticle}
    U, λs = cache.U, cache.λs
    s = (P === Fermion) ? 1.0 : -1.0
    ns = [thermaloccupation(P, β, μ, λ) for λ in λs]
    ws = [U[i, k] * conj(U[j, k]) for k in eachindex(λs)]
    gk(k) = (P === Fermion) ? 1.0 - ns[k] : 1.0 + ns[k]
    Nt, Nτ = length(ts), length(τs)
    G₊₋, G₋₊ = _lesser_greater_blocks(P, cache, ts, i, j; β=β, μ=μ)
    # (+,τ): no contour swap, -i⟨x(t)x†(τ)⟩ — greater weight
    G₊τ = [-im * sum(ws[k] * gk(k) * exp(λs[k] * τs[m] - im * λs[k] * ts[a]) for k in eachindex(λs))
           for a in 1:Nt, m in 1:Nτ]
    # (τ,+): one contour swap, +i⟨x†(t)x(τ)⟩ (fermions) / -i⟨x†(t)x(τ)⟩ (bosons) — lesser weight;
    # x(τ) = e^{-τε}x decays in imaginary time
    Gτ₊ = [s * im * sum(ws[k] * ns[k] * exp(-λs[k] * τs[m] + im * λs[k] * ts[a]) for k in eachindex(λs))
           for m in 1:Nτ, a in 1:Nt]
    Gττ = _matsubara_block(P, cache, τs, i, j; β=β, μ=μ)
    return G₊₋, G₋₊, G₊τ, Gτ₊, Gττ
end

# —— single-(site-)pair helpers shared by the normal constructors ————
function _keldysh_pair(::Type{P}, cache::EigenCache, ts::AbstractVector{<:Real},
                       i::Int, j::Int; β::Real, μ::Real) where {P<:AbstractParticle}
    G₊₋, G₋₊ = _lesser_greater_blocks(P, cache, ts, i, j; β=β, μ=μ)
    return RealGreenFunctions(G₊₋, G₋₊)
end
function _mixed_pair(::Type{P}, cache::EigenCache, ts::AbstractVector{<:Real}, τs::AbstractVector{<:Real},
                     i::Int, j::Int; β::Real, μ::Real) where {P<:AbstractParticle}
    G₊₋, G₋₊, G₊τ, Gτ₊, Gττ = _mixed_blocks(P, cache, ts, τs, i, j; β=β, μ=μ)
    return MixedGreenFunctions(G₊₋, G₋₊, G₊τ, Gτ₊, Gττ)
end
function _matsubara_pair(::Type{P}, cache::EigenCache, τs::AbstractVector{<:Real},
                         i::Int, j::Int; β::Real, μ::Real) where {P<:AbstractParticle}
    return ImagGreenFunctions(_matsubara_block(P, cache, τs, i, j; β=β, μ=μ))
end


# —— equilibrium Keldysh contour ——————————————————————————————————————————————

"""
    freefermions_keldysh_gfs(h, ts, i, j=i, cache=eigencache(h); β, μ=0) -> RealGreenFunctions

Contour-ordered Green's functions of free fermions in thermal equilibrium
(coefficient density matrix `thermocdm` at inverse temperature `β` and chemical
potential `μ`) on the closed real-time Keldysh contour, between the sites `i`
and `j` of the single-particle Hamiltonian `h`, evaluated on the real-time
nodes `ts` (t₁ = 0, …).  The blocks are accessed through `branch(x; b1, b2)`;
see `RealGreenFunctions` for the conventions.
"""
function freefermions_keldysh_gfs(h::AbstractMatrix, ts::AbstractVector{<:Real}, i::Int, j::Int=i,
                                  cache::EigenCache=eigencache(h); β::Real, μ::Real=0)
    return _keldysh_pair(Fermion, cache, ts, i, j; β=β, μ=μ)
end

"""
    freebosons_keldysh_gfs(h, ts, i, j=i, cache=eigencache(h); β, μ=0) -> RealGreenFunctions

Same as `freefermions_keldysh_gfs` for free bosons; see `RealGreenFunctions`
for the block conventions.
"""
function freebosons_keldysh_gfs(h::AbstractMatrix, ts::AbstractVector{<:Real}, i::Int, j::Int=i,
                                cache::EigenCache=eigencache(h); β::Real, μ::Real=0)
    return _keldysh_pair(Boson, cache, ts, i, j; β=β, μ=μ)
end


# —— nonequilibrium Keldysh contour (initial cdm ρ₀) ——————————————————————————

"""
    freefermions_keldysh_gfs(h, ρ₀, ts, i, j=i, cache=eigencache(h)) -> RealGreenFunctions
    freebosons_keldysh_gfs(h, ρ₀, ts, i, j=i, cache=eigencache(h)) -> RealGreenFunctions

Nonequilibrium contour-ordered Green's functions starting from the initial
coefficient density matrix `ρ₀` (see `freefermions_greater_lesser` for the
semantics of `ρ₀`), on the closed real-time Keldysh contour between sites `i`
and `j`, evaluated on the real-time nodes `ts`.  The mixed (imaginary-time)
blocks of the Kadanoff contour have no nonequilibrium generalization, so no
`ρ₀` variant of `free*_kadanoff_gfs` exists.
"""
function freefermions_keldysh_gfs(h::AbstractMatrix, ρ₀::AbstractMatrix, ts::AbstractVector{<:Real},
                                  i::Int, j::Int=i, cache::EigenCache=eigencache(h))
    f = freefermions_greater_lesser(h, ρ₀, i, j, cache)
    return _keldysh_from_greater_lesser(Fermion, f, ts)
end
function freebosons_keldysh_gfs(h::AbstractMatrix, ρ₀::AbstractMatrix, ts::AbstractVector{<:Real},
                                i::Int, j::Int=i, cache::EigenCache=eigencache(h))
    f = freebosons_greater_lesser(h, ρ₀, i, j, cache)
    return _keldysh_from_greater_lesser(Boson, f, ts)
end

# f(t) -> (G>(t), G<(t)) in the greater/lesser convention of
# free*_greater_lesser; the bosonic physical lesser block is the negative of
# the codebase lesser kernel (contour value -i⟨b†b⟩).
function _keldysh_from_greater_lesser(::Type{P}, f::Function, ts::AbstractVector{<:Real}) where {P<:AbstractParticle}
    Nt = length(ts)
    G₊₋ = zeros(ComplexF64, Nt, Nt)
    G₋₊ = zeros(ComplexF64, Nt, Nt)
    for a in 1:Nt, b in 1:Nt
        g, l = f(ts[a] - ts[b])
        G₊₋[a, b] = (P === Fermion) ? l : -l
        G₋₊[a, b] = g
    end
    return RealGreenFunctions(G₊₋, G₋₊)
end


# —— Kadanoff (L-shaped) contour ———————————————————————————————————————————————

"""
    freefermions_kadanoff_gfs(h, ts, τs, i, j=i, cache=eigencache(h); β, μ=0) -> MixedGreenFunctions

Contour-ordered Green's functions of free fermions in thermal equilibrium on
the L-shaped (Kadanoff-Baym) contour between sites `i` and `j`: an
imaginary-time branch with nodes `τs` (0 ≤ τ < β) followed by the closed
real-time part on the nodes `ts`.  The blocks are accessed through
`branch(x; b1, b2)`; see `MixedGreenFunctions` for the conventions.
"""
function freefermions_kadanoff_gfs(h::AbstractMatrix, ts::AbstractVector{<:Real}, τs::AbstractVector{<:Real},
                                   i::Int, j::Int=i, cache::EigenCache=eigencache(h); β::Real, μ::Real=0)
    G₊₋, G₋₊, G₊τ, Gτ₊, Gττ = _mixed_blocks(Fermion, cache, ts, τs, i, j; β=β, μ=μ)
    return MixedGreenFunctions(G₊₋, G₋₊, G₊τ, Gτ₊, Gττ)
end

"""
    freebosons_kadanoff_gfs(h, ts, τs, i, j=i, cache=eigencache(h); β, μ=0) -> MixedGreenFunctions

Same as `freefermions_kadanoff_gfs` for free bosons; see `MixedGreenFunctions`
for the block conventions.
"""
function freebosons_kadanoff_gfs(h::AbstractMatrix, ts::AbstractVector{<:Real}, τs::AbstractVector{<:Real},
                                 i::Int, j::Int=i, cache::EigenCache=eigencache(h); β::Real, μ::Real=0)
    G₊₋, G₋₊, G₊τ, Gτ₊, Gττ = _mixed_blocks(Boson, cache, ts, τs, i, j; β=β, μ=μ)
    return MixedGreenFunctions(G₊₋, G₋₊, G₊τ, Gτ₊, Gττ)
end


# —— Matsubara (imaginary-time) axis ——————————————————————————————————————————

"""
    freefermions_matsubara_gfs(h, τs, i, j=i, cache=eigencache(h); β, μ=0) -> ImagGreenFunctions

Matsubara Green's functions of free fermions in thermal equilibrium
(coefficient density matrix `thermocdm` at inverse temperature `β` and
chemical potential `μ`) between the sites `i` and `j` of the single-particle
Hamiltonian `h`, evaluated on the imaginary-time nodes `τs` (0 ≤ τ ≤ β).  The
field `data[m, n] = -i⟨T̂ x(τₘ) x†(τₙ)⟩` follows the (τ,τ) block convention of
`MixedGreenFunctions`, so on matching grids
`branch(freefermions_kadanoff_gfs(h, ts, τs, ...), :τ, :τ) == freefermions_matsubara_gfs(h, τs, ...).data`.
"""
function freefermions_matsubara_gfs(h::AbstractMatrix, τs::AbstractVector{<:Real}, i::Int, j::Int=i,
                                    cache::EigenCache=eigencache(h); β::Real, μ::Real=0)
    return ImagGreenFunctions(_matsubara_block(Fermion, cache, τs, i, j; β=β, μ=μ))
end

"""
    freebosons_matsubara_gfs(h, τs, i, j=i, cache=eigencache(h); β, μ=0) -> ImagGreenFunctions

Same as `freefermions_matsubara_gfs` for free bosons; see `ImagGreenFunctions`
for the data conventions.
"""
function freebosons_matsubara_gfs(h::AbstractMatrix, τs::AbstractVector{<:Real}, i::Int, j::Int=i,
                                  cache::EigenCache=eigencache(h); β::Real, μ::Real=0)
    return ImagGreenFunctions(_matsubara_block(Boson, cache, τs, i, j; β=β, μ=μ))
end
