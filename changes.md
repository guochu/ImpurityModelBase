# Changes

## Unreleased

### New: contour-ordered Green's functions (Keldysh / Kadanoff / Matsubara containers)

New container types for the contour-ordered single-particle Green's functions of
free-particle baths, with uniform `branch` accessors and Fock-benchmarked
constructors, intended as the ED reference for downstream TEMPO / GTEMPO
benchmarks (all in `src/exactdiagonalizations/contourgreenfunctions.jl`):

- `RealGreenFunctions` — closed real-time (Keldysh) contour; stores the lesser
  (`G₊₋`) and greater (`G₋₊`) cross blocks over the real-time nodes;
  `branch(x; b1, b2)` with `b1, b2 ∈ (:+, :-)` derives all four contour blocks
  (the diagonal blocks through the contour-ordering sector rule).
- `MixedGreenFunctions` — L-shaped (Kadanoff–Baym) contour; additionally stores
  the mixed (`G₊τ`, `Gτ₊`) and Matsubara (`Gττ`) blocks; `branch` accepts
  `b1, b2 ∈ (:+, :-, :τ)`.
- `ImagGreenFunctions{T}` — Matsubara-only container with a single field
  `data::Matrix{T}`, `data[m, n] = −i⟨T̂ x(τₘ) x†(τₙ)⟩`; the data equals the
  `(τ,τ)` Kadanoff block, so
  `branch(kadanoff, :τ, :τ) == matsubara.data` on matching grids.
- `branches(x)` — the contour branches carried by the container
  (`(:+, :-)` / `(:+, :-, :τ)` / `(:τ,)`).

Constructors (thermal equilibrium; the lesser/greater weights and statistics
signs follow the existing `free*_greater_lesser` conventions):

- `freefermions_keldysh_gfs(h, ts, i, j=i, cache; β, μ=0)` /
  `freebosons_keldysh_gfs(h, ts, i, j=i, cache; β, μ=0)` → `RealGreenFunctions`
- `freefermions_kadanoff_gfs(h, ts, τs, i, j=i, cache; β, μ=0)` /
  `freebosons_kadanoff_gfs(h, ts, τs, i, j=i, cache; β, μ=0)` → `MixedGreenFunctions`
- `freefermions_matsubara_gfs(h, τs, i, j=i, cache; β, μ=0)` /
  `freebosons_matsubara_gfs(h, τs, i, j=i, cache; β, μ=0)` → `ImagGreenFunctions`

All benchmarks (single- and two-mode baths; Keldysh, Kadanoff and Matsubara
blocks; KMS identities; retarded-function consistency with `freefermions_Gt` /
`freebosons_Gt`) live in `test/exactdiagonalizations/greenfunctions.jl`.

### `freebosons_Gt` now returns the physical retarded Green's function (interface change)

`freebosons_Gt` (both the equilibrium and the `ρ₀` variants) returned the wrong object:
it applied the fermionic formula `G>(t) − G<(t)` verbatim, but the codebase bosonic
lesser carries the correlation convention

    G<(t) = +i⟨b†(0)b(t)⟩ = −(physical, Haug–Jauho lesser)

so the old return value was `G>(t) − G<(t) = −i(1 + 2n̄ₖ)e^{−iεₖt}` per mode — *not* a
Green's function of the commutator. It is now the physical retarded function

    Gᴿ(t) = −iθ(t)⟨[b(t), b†(0)]⟩ = G>(t) + G<(t) = −i e^{−iεₖt} per mode,

computed as the **sum** of the two outputs of `freebosons_greater_lesser`, consistent
with the analytic solutions, which already used
`freeboson_Gt = freeboson_greater + freeboson_lesser`.

**Migration.** Code that relied on the old return value should replace it with
`G>(t) − G<(t)` built from `freebosons_greater_lesser` directly:

```julia
# old (buggy) behaviour, recoverable as:
gl = freebosons_greater_lesser(h, i, j; β=β, μ=μ)
old_gt(t) = gl(t)[1] - gl(t)[2]
# new behaviour:
gt = freebosons_Gt(h, i, j; β=β, μ=μ)   # gt(t) == gl(t)[1] + gl(t)[2] == Gᴿ(t)
```

i.e. `new = old + 2·(second output of freebosons_greater_lesser)`.

**Unaffected:** `freefermions_Gt` (its `G> − G<` was already the physical retarded
function), the outputs of `freebosons_greater_lesser` / `freefermions_greater_lesser`
(convention unchanged, now documented), `free*_Gτ`, and all analytic solutions.
The contour Green's functions (`free*_keldysh_gfs`) were already built on the physical
convention and are unaffected.

### BCS contour Green's functions: attempted and withdrawn

An implementation of `BCSGreenFunctions{C}` (normal uu / pair-creation ud /
pair-annihilation du families over the particle-hole doubled BdG space, with
`bcs_keldysh_gfs` / `bcs_kadanoff_gfs` / `bcs_matsubara_gfs`) was attempted and
**removed again**: the normal family (uu) and the Δ = 0 reduction validate against
a Fock-space benchmark, but the conventions of the anomalous families (the ud/du
Wick contractions and their Kadanoff τ-sectors) could not be pinned down — the
directly computed slot-slot covariance disagrees with every assembly formula
derivable from the extended cdm's particle-hole blocks. The benchmark (kept as a
reference for a future attempt) can be recovered from the git history of
`test/exactdiagonalizations/greenfunctions.jl`. Completing this requires
confirming the conventions against the Fock-tested BCS machinery of the QuAPI
package.

### Test environment: `Pkg.test()` fixed

`[extras]`/`[targets]` sections were added to `Project.toml` declaring the `Test` and
`Random` stdlibs, so `Pkg.test()` no longer fails with "Package Test not found".
