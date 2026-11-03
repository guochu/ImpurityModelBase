println("------------------------------------")
println("|     Contour Green's functions    |")
println("------------------------------------")

using LinearAlgebra

# The contour Green's functions (RealGreenFunctions / MixedGreenFunctions) are
# benchmarked against the exact many-body solution on the (truncated) Fock
# space: single-mode (L=1) and two-mode (L=2) baths, with every contour block
# (including the derived diagonal blocks of `branch`) compared against
# independently assembled contour-ordered correlators.  The bosonic Fock spaces
# are truncated with cutoffs chosen so that the occupation tails stay far below
# the test tolerance.

# —— full-Fock-space machinery ————————————————————————————————————————————————

function _thermal_fock_dm(H::AbstractMatrix, β::Real)
	ρ = exp(-β .* H)
	ρ ./= tr(ρ)
	return ρ
end

_corr(ρ, A, B, trρ) = tr(ρ * A * B) / trρ

# diagonal product-state (geometric) bosonic Fock dm realizing the diagonal cdm n̄
function _boson_product_dm(n̄::AbstractVector{<:Real}, d::Int)
	L = length(n̄)
	dim = d^L
	p = zeros(Float64, dim)
	for idx0 in 0:dim-1
		tmp, w = idx0, 1.0
		for k in L:-1:1        # mode k is the (L-k)-th least significant digit
			n, tmp = tmp % d, tmp ÷ d
			r = n̄[k] / (1 + n̄[k])
			w *= r^n * (1 - r)
		end
		p[idx0+1] = w
	end
	return Diagonal(p ./ sum(p))
end

# one-particle (Slater determinant) Fock state vector; its cdm is w·w' with w = conj(v)
function _slater_state(v::AbstractVector{<:Number})
	L = length(v)
	ψ = zeros(ComplexF64, 2^L)
	for p in 1:L
		ψ[1 + 2^(L - p)] = v[p]
	end
	return ψ
end

function _fermion_fock_ctx(h::AbstractMatrix; β::Real=0.0, μ::Real=0.0, ρ=nothing)
	L = size(h, 1)
	cs = [ComplexF64.(fermionaoperator(L, p)) for p in 1:L]
	cds = [adjoint(c) for c in cs]
	H = zeros(ComplexF64, 2^L, 2^L)
	for a in 1:L, b in 1:L
		H .+= h[a, b] .* ComplexF64.(fermionoperator(L, adaga(a, b, coeff=1.0)))
	end
	ρv = ρ === nothing ? _thermal_fock_dm(H - μ * ComplexF64.(fermiondensityoperator(L)), β) : ρ
	λs, U = eigen(Hermitian(H))
	Ut(t) = U * Diagonal(exp.(im .* λs .* t)) * U'     # e^{+iHt}
	Udt(t) = U * Diagonal(exp.(-im .* λs .* t)) * U'   # e^{-iHt}
	Tm(τ) = U * Diagonal(exp.(λs .* τ)) * U'           # e^{+τH}
	Tmi(τ) = U * Diagonal(exp.(-λs .* τ)) * U'         # e^{-τH}
	return (ρ=ρv, trρ=tr(ρv), cs=cs, cds=cds, Ut=Ut, Udt=Udt, Tm=Tm, Tmi=Tmi)
end

function _boson_fock_ctx(h::AbstractMatrix, d::Int; β::Real=0.0, μ::Real=0.0, ρ=nothing)
	L = size(h, 1)
	cs = [ComplexF64.(bosonaoperator(L, p, d=d)) for p in 1:L]
	cds = [adjoint(c) for c in cs]
	H = zeros(ComplexF64, d^L, d^L)
	for a in 1:L, b in 1:L
		H .+= h[a, b] .* ComplexF64.(bosonoperator(L, adaga(a, b, coeff=1.0), d=d))
	end
	ρv = ρ === nothing ? _thermal_fock_dm(H - μ * ComplexF64.(bosondensityoperator(L, d=d)), β) : ρ
	λs, U = eigen(Hermitian(H))
	Ut(t) = U * Diagonal(exp.(im .* λs .* t)) * U'
	Udt(t) = U * Diagonal(exp.(-im .* λs .* t)) * U'
	Tm(τ) = U * Diagonal(exp.(λs .* τ)) * U'
	Tmi(τ) = U * Diagonal(exp.(-λs .* τ)) * U'
	return (ρ=ρv, trρ=tr(ρv), cs=cs, cds=cds, Ut=Ut, Udt=Udt, Tm=Tm, Tmi=Tmi)
end

# real-time Wightman data: ⟨xᵢ(t)xⱼ†(0)⟩ (greater) and ⟨xⱼ†(0)xᵢ(t)⟩ (lesser)
_wightman(ctx, i, j, t) =
	(_corr(ctx.ρ, ctx.Ut(t) * ctx.cs[i] * ctx.Udt(t), ctx.cds[j], ctx.trρ),
	 _corr(ctx.ρ, ctx.cds[j], ctx.Ut(t) * ctx.cs[i] * ctx.Udt(t), ctx.trρ))

@testset "fermionic keldysh gfs: single & two modes vs full Fock space" begin
	tol = 1.0e-10
	ts = collect(0.0:0.25:1.5)
	Nt = length(ts)
	for h in (ComplexF64[1.0;;], ComplexF64[1.0 0.3; 0.3 2.0])
		L = size(h, 1)
		for (β, μ) in ((1.3, 0.0), (0.8, 0.45))
			ctx = _fermion_fock_ctx(h; β=β, μ=μ)
			for i in 1:min(L, 2), j in unique((i, min(L, i % L + 1)))
				gfs = freefermions_keldysh_gfs(h, ts, i, j; β=β, μ=μ)
				Wg(t) = _wightman(ctx, i, j, t)[1]
				Wl(t) = _wightman(ctx, i, j, t)[2]
				# stored cross blocks
				G₊₋ = [im * Wl(ts[a] - ts[b]) for a in 1:Nt, b in 1:Nt]
				G₋₊ = [-im * Wg(ts[a] - ts[b]) for a in 1:Nt, b in 1:Nt]
				@test branch(gfs, :+, :-) ≈ G₊₋ atol = tol
				@test branch(gfs, :-, :+) ≈ G₋₊ atol = tol
				# derived diagonal blocks (contour-ordering sector rule)
				G₊₊ = [a > b ? G₋₊[a, b] : G₊₋[a, b] for a in 1:Nt, b in 1:Nt]
				G₋₋ = [a > b ? G₊₋[a, b] : G₋₊[a, b] for a in 1:Nt, b in 1:Nt]
				@test branch(gfs, :+, :+) ≈ G₊₊ atol = tol
				@test branch(gfs, :-, :-) ≈ G₋₋ atol = tol
				# equal-time canonical anticommutator contact: G>(0) - G<(0) = -i δᵢⱼ
				@test maximum(abs.(G₋₊[a, a] - G₊₋[a, a] .+ (i == j ? im : 0.0im)) for a in 1:Nt) < tol
			end
		end
	end
end

@testset "bosonic keldysh gfs: single & two modes vs truncated Fock space" begin
	tol = 1.0e-9
	ts = collect(0.0:0.2:1.2)
	Nt = length(ts)
	cases = ((ComplexF64[1.2;;], 32, [(2.0, 0.0), (2.0, 0.3)]),
			 (ComplexF64[1.1 0.25; 0.25 1.9], 16, [(3.0, 0.0)]))
	for (h, d, βμs) in cases
		L = size(h, 1)
		for (β, μ) in βμs
			ctx = _boson_fock_ctx(h, d; β=β, μ=μ)
			for i in 1:min(L, 2), j in unique((i, min(L, i % L + 1)))
				gfs = freebosons_keldysh_gfs(h, ts, i, j; β=β, μ=μ)
				Wg(t) = _wightman(ctx, i, j, t)[1]
				Wl(t) = _wightman(ctx, i, j, t)[2]
				# stored cross blocks (bosonic lesser contour value is -i⟨b†b⟩)
				G₊₋ = [-im * Wl(ts[a] - ts[b]) for a in 1:Nt, b in 1:Nt]
				G₋₊ = [-im * Wg(ts[a] - ts[b]) for a in 1:Nt, b in 1:Nt]
				@test branch(gfs, :+, :-) ≈ G₊₋ atol = tol
				@test branch(gfs, :-, :+) ≈ G₋₊ atol = tol
				# derived diagonal blocks (same sector rule, bosonic swap signs)
				G₊₊ = [a > b ? G₋₊[a, b] : G₊₋[a, b] for a in 1:Nt, b in 1:Nt]
				G₋₋ = [a > b ? G₊₋[a, b] : G₋₊[a, b] for a in 1:Nt, b in 1:Nt]
				@test branch(gfs, :+, :+) ≈ G₊₊ atol = tol
				@test branch(gfs, :-, :-) ≈ G₋₋ atol = tol
				# equal-time canonical commutator contact: G>(0) - G<(0) = -i δᵢⱼ
				@test maximum(abs.(G₋₊[a, a] - G₊₋[a, a] .+ (i == j ? im : 0.0im)) for a in 1:Nt) < tol
			end
		end
	end
end

@testset "nonequilibrium keldysh gfs (ρ₀) vs Fock space" begin
	tol = 1.0e-10
	ts = collect(0.0:0.25:1.5)
	Nt = length(ts)

	# fermionic: Slater determinant with cdm ρ₀ = w w′, w = conj(v)
	h = ComplexF64[1.0 0.35; 0.35 2.0]
	ctx = _fermion_fock_ctx(h)
	v = normalize(ComplexF64[0.7, 0.3 + 0.4im])
	ψ = _slater_state(v)
	w = conj(v)
	gfs = freefermions_keldysh_gfs(h, w * w', ts, 1, 2)
	G₊₋ = zeros(ComplexF64, Nt, Nt)
	G₋₊ = zeros(ComplexF64, Nt, Nt)
	for a in 1:Nt, b in 1:Nt
		t = ts[a] - ts[b]
		G₊₋[a, b] = im * dot(ψ, ctx.cds[2] * (ctx.Ut(t) * ctx.cs[1] * ctx.Udt(t) * ψ))
		G₋₊[a, b] = -im * dot(ψ, ctx.Ut(t) * ctx.cs[1] * ctx.Udt(t) * (ctx.cds[2] * ψ))
	end
	@test branch(gfs, :+, :-) ≈ G₊₋ atol = tol
	@test branch(gfs, :-, :+) ≈ G₋₊ atol = tol
	@test branch(gfs, :+, :+) ≈ [a > b ? G₋₊[a, b] : G₊₋[a, b] for a in 1:Nt, b in 1:Nt] atol = tol
	@test branch(gfs, :-, :-) ≈ [a > b ? G₊₋[a, b] : G₋₊[a, b] for a in 1:Nt, b in 1:Nt] atol = tol

	# bosonic: product state of modes with occupations n̄ (diagonal cdm)
	hb = ComplexF64[1.1 0.25; 0.25 1.9]
	d = 20
	n̄ = [0.15, 0.25]
	ρmb = _boson_product_dm(n̄, d)
	ctxb = _boson_fock_ctx(hb, d; ρ=ρmb)
	gfsb = freebosons_keldysh_gfs(hb, Diagonal(n̄), ts, 1, 2)
	G₊₋ = zeros(ComplexF64, Nt, Nt)
	G₋₊ = zeros(ComplexF64, Nt, Nt)
	for a in 1:Nt, b in 1:Nt
		t = ts[a] - ts[b]
		Wg = _corr(ρmb, ctxb.Ut(t) * ctxb.cs[1] * ctxb.Udt(t), ctxb.cds[2], ctxb.trρ)
		Wl = _corr(ρmb, ctxb.cds[2], ctxb.Ut(t) * ctxb.cs[1] * ctxb.Udt(t), ctxb.trρ)
		G₊₋[a, b] = -im * Wl
		G₋₊[a, b] = -im * Wg
	end
	@test branch(gfsb, :+, :-) ≈ G₊₋ atol = tol
	@test branch(gfsb, :-, :+) ≈ G₋₊ atol = tol
	@test branch(gfsb, :+, :+) ≈ [a > b ? G₋₊[a, b] : G₊₋[a, b] for a in 1:Nt, b in 1:Nt] atol = tol
	@test branch(gfsb, :-, :-) ≈ [a > b ? G₊₋[a, b] : G₋₊[a, b] for a in 1:Nt, b in 1:Nt] atol = tol
end

@testset "fermionic kadanoff gfs: single & two modes vs full Fock space" begin
	tol = 1.0e-10
	ts = collect(0.0:0.25:1.5)
	Nt = length(ts)
	for h in (ComplexF64[1.0;;], ComplexF64[1.0 0.3; 0.3 2.0])
		L = size(h, 1)
		for (β, μ) in ((1.3, 0.0), (0.8, 0.45))
			Nτ = 5
			τs = (0:Nτ-1) .* (β / Nτ)
			ctx = _fermion_fock_ctx(h; β=β, μ=μ)
			for i in 1:min(L, 2), j in unique((i, min(L, i % L + 1)))
				gfs = freefermions_kadanoff_gfs(h, ts, τs, i, j; β=β, μ=μ)
				Wg(t) = _wightman(ctx, i, j, t)[1]
				Wl(t) = _wightman(ctx, i, j, t)[2]
				G₊₋ = [im * Wl(ts[a] - ts[b]) for a in 1:Nt, b in 1:Nt]
				G₋₊ = [-im * Wg(ts[a] - ts[b]) for a in 1:Nt, b in 1:Nt]
				# mixed blocks: (±,τ) greater weight, (τ,±) lesser weight
				G₊τ = [-im * _corr(ctx.ρ, ctx.Ut(ts[a]) * ctx.cs[i] * ctx.Udt(ts[a]),
								   ctx.Tm(τs[m]) * ctx.cds[j] * ctx.Tmi(τs[m]), ctx.trρ)
					   for a in 1:Nt, m in 1:Nτ]
				Gτ₊ = [im * _corr(ctx.ρ, ctx.Ut(ts[a]) * ctx.cds[j] * ctx.Udt(ts[a]),
								  ctx.Tm(τs[m]) * ctx.cs[i] * ctx.Tmi(τs[m]), ctx.trρ)
					   for m in 1:Nτ, a in 1:Nt]
				Gττ = zeros(ComplexF64, Nτ, Nτ)
				for m in 1:Nτ, n in 1:Nτ
					Δ = τs[m] - τs[n]
					# +i·Gᴹ(Δ), Gᴹ = -⟨T̂ cᵢ(Δ)cⱼ†⟩ with fermionic swap sign
					Gττ[m, n] = Δ >= 0 ? -im * _corr(ctx.ρ, ctx.Tm(Δ) * ctx.cs[i] * ctx.Tmi(Δ), ctx.cds[j], ctx.trρ) :
										 im * _corr(ctx.ρ, ctx.Tm(-Δ) * ctx.cds[j] * ctx.Tmi(-Δ), ctx.cs[i], ctx.trρ)
				end
				# all nine contour blocks
				@test branch(gfs, :+, :-) ≈ G₊₋ atol = tol
				@test branch(gfs, :-, :+) ≈ G₋₊ atol = tol
				@test branch(gfs, :+, :+) ≈ [a > b ? G₋₊[a, b] : G₊₋[a, b] for a in 1:Nt, b in 1:Nt] atol = tol
				@test branch(gfs, :-, :-) ≈ [a > b ? G₊₋[a, b] : G₋₊[a, b] for a in 1:Nt, b in 1:Nt] atol = tol
				@test branch(gfs, :+, :τ) ≈ G₊τ atol = tol
				@test branch(gfs, :-, :τ) ≈ G₊τ atol = tol
				@test branch(gfs, :τ, :+) ≈ Gτ₊ atol = tol
				@test branch(gfs, :τ, :-) ≈ Gτ₊ atol = tol
				@test branch(gfs, :τ, :τ) ≈ Gττ atol = tol
			end
		end
	end
end

@testset "bosonic kadanoff gfs: single & two modes vs truncated Fock space" begin
	tol = 1.0e-9
	ts = collect(0.0:0.2:1.2)
	Nt = length(ts)
	# τ grid kept at ≤ β/2 so that the e^{+τε} growth of the (±,τ) blocks does
	# not amplify the Fock truncation tails beyond the tolerance
	cases = ((ComplexF64[1.2;;], 40, [(2.0, 0.0), (2.0, 0.3)]),
			 (ComplexF64[1.1 0.25; 0.25 1.9], 18, [(3.0, 0.0)]))
	for (h, d, βμs) in cases
		L = size(h, 1)
		for (β, μ) in βμs
			Nτ = 5
			τs = (0:Nτ-1) .* (0.5 * β / Nτ)
			ctx = _boson_fock_ctx(h, d; β=β, μ=μ)
			for i in 1:min(L, 2), j in unique((i, min(L, i % L + 1)))
				gfs = freebosons_kadanoff_gfs(h, ts, τs, i, j; β=β, μ=μ)
				Wg(t) = _wightman(ctx, i, j, t)[1]
				Wl(t) = _wightman(ctx, i, j, t)[2]
				G₊₋ = [-im * Wl(ts[a] - ts[b]) for a in 1:Nt, b in 1:Nt]
				G₋₊ = [-im * Wg(ts[a] - ts[b]) for a in 1:Nt, b in 1:Nt]
				G₊τ = [-im * _corr(ctx.ρ, ctx.Ut(ts[a]) * ctx.cs[i] * ctx.Udt(ts[a]),
								   ctx.Tm(τs[m]) * ctx.cds[j] * ctx.Tmi(τs[m]), ctx.trρ)
					   for a in 1:Nt, m in 1:Nτ]
				Gτ₊ = [-im * _corr(ctx.ρ, ctx.Ut(ts[a]) * ctx.cds[j] * ctx.Udt(ts[a]),
								   ctx.Tm(τs[m]) * ctx.cs[i] * ctx.Tmi(τs[m]), ctx.trρ)
					   for m in 1:Nτ, a in 1:Nt]
				Gττ = zeros(ComplexF64, Nτ, Nτ)
				for m in 1:Nτ, n in 1:Nτ
					Δ = τs[m] - τs[n]
					# +i·Gᴹ(Δ), bosonic Matsubara T̂ without swap signs
					Gττ[m, n] = Δ >= 0 ? -im * _corr(ctx.ρ, ctx.Tm(Δ) * ctx.cs[i] * ctx.Tmi(Δ), ctx.cds[j], ctx.trρ) :
										 -im * _corr(ctx.ρ, ctx.Tm(-Δ) * ctx.cds[j] * ctx.Tmi(-Δ), ctx.cs[i], ctx.trρ)
				end
				@test branch(gfs, :+, :-) ≈ G₊₋ atol = tol
				@test branch(gfs, :-, :+) ≈ G₋₊ atol = tol
				@test branch(gfs, :+, :+) ≈ [a > b ? G₋₊[a, b] : G₊₋[a, b] for a in 1:Nt, b in 1:Nt] atol = tol
				@test branch(gfs, :-, :-) ≈ [a > b ? G₊₋[a, b] : G₋₊[a, b] for a in 1:Nt, b in 1:Nt] atol = tol
				@test branch(gfs, :+, :τ) ≈ G₊τ atol = tol
				@test branch(gfs, :-, :τ) ≈ G₊τ atol = tol
				@test branch(gfs, :τ, :+) ≈ Gτ₊ atol = tol
				@test branch(gfs, :τ, :-) ≈ Gτ₊ atol = tol
				@test branch(gfs, :τ, :τ) ≈ Gττ atol = tol
			end
		end
	end
end

@testset "fermionic matsubara gfs: single & two modes vs full Fock space" begin
	tol = 1.0e-10
	for h in (ComplexF64[1.0;;], ComplexF64[1.0 0.3; 0.3 2.0])
		L = size(h, 1)
		for (β, μ) in ((1.3, 0.0), (0.8, 0.45))
			Nτ = 6
			τs = (0:Nτ-1) .* (β / Nτ)
			ctx = _fermion_fock_ctx(h; β=β, μ=μ)
			for i in 1:min(L, 2), j in unique((i, min(L, i % L + 1)))
				gfs = freefermions_matsubara_gfs(h, τs, i, j; β=β, μ=μ)
				G = zeros(ComplexF64, Nτ, Nτ)
				for m in 1:Nτ, n in 1:Nτ
					Δ = τs[m] - τs[n]
					# +i·Gᴹ(Δ), Gᴹ = -⟨T̂ cᵢ(Δ)cⱼ†⟩ with fermionic swap sign
					G[m, n] = Δ >= 0 ? -im * _corr(ctx.ρ, ctx.Tm(Δ) * ctx.cs[i] * ctx.Tmi(Δ), ctx.cds[j], ctx.trρ) :
									   im * _corr(ctx.ρ, ctx.Tm(-Δ) * ctx.cds[j] * ctx.Tmi(-Δ), ctx.cs[i], ctx.trρ)
				end
				@test gfs.data ≈ G atol = tol
				# fermionic KMS Gᴹ(τ − β) = -e^{βμ}·Gᴹ(τ) on a grid through β
				τsβ = (0:Nτ) .* (β / Nτ)
				gβ = freefermions_matsubara_gfs(h, τsβ, i, j; β=β, μ=μ).data
				for k in 2:Nτ
					@test gβ[k, Nτ+1] ≈ -exp(β * μ) * gβ[k, 1] atol = tol
				end
			end
		end
	end
end

@testset "bosonic matsubara gfs: single & two modes vs truncated Fock space" begin
	tol = 1.0e-9
	# cutoffs keep the e^{±λΔ} dynamic range of the many-body reference
	# (λmax·Δmax) in the numerically safe regime while the thermal tails
	# stay far below the tolerance
	cases = ((ComplexF64[1.2;;], 24, [(2.0, 0.0), (2.0, 0.3)], 6),
			 (ComplexF64[1.1 0.25; 0.25 1.9], 10, [(3.0, 0.0)], 3))
	for (h, d, βμs, Nτβ) in cases
		L = size(h, 1)
		for (β, μ) in βμs
			Nτ = 6
			τs = (0:Nτ-1) .* (β / Nτ)
			ctx = _boson_fock_ctx(h, d; β=β, μ=μ)
			for i in 1:min(L, 2), j in unique((i, min(L, i % L + 1)))
				gfs = freebosons_matsubara_gfs(h, τs, i, j; β=β, μ=μ)
				G = zeros(ComplexF64, Nτ, Nτ)
				for m in 1:Nτ, n in 1:Nτ
					Δ = τs[m] - τs[n]
					# +i·Gᴹ(Δ), bosonic Matsubara T̂ without swap signs
					G[m, n] = Δ >= 0 ? -im * _corr(ctx.ρ, ctx.Tm(Δ) * ctx.cs[i] * ctx.Tmi(Δ), ctx.cds[j], ctx.trρ) :
									   -im * _corr(ctx.ρ, ctx.Tm(-Δ) * ctx.cds[j] * ctx.Tmi(-Δ), ctx.cs[i], ctx.trρ)
				end
				@test gfs.data ≈ G atol = tol
				# bosonic KMS Gᴹ(τ − β) = +e^{βμ}·Gᴹ(τ) on a grid through β
				τsβ = (0:Nτβ) .* (β / Nτβ)
				gβ = freebosons_matsubara_gfs(h, τsβ, i, j; β=β, μ=μ).data
				for k in 2:Nτβ
					@test gβ[k, Nτβ+1] ≈ exp(β * μ) * gβ[k, 1] atol = tol
				end
			end
		end
	end
end

@testset "contour identities and API" begin
	h = ComplexF64[1.0 0.3; 0.3 2.0]
	β, μ = 1.3, 0.2
	ts = collect(0.0:0.25:1.5)
	Nt = length(ts)
	ts6 = (0:6) .* (1.5 / 6)
	τs6 = (0:5) .* (β / 6)
	for (gfs, Gt, gfs11) in
		((freefermions_keldysh_gfs(h, ts, 1, 2; β=β, μ=μ), freefermions_Gt(h, ts, 1, 2; β=β, μ=μ),
		  freefermions_keldysh_gfs(h, ts, 1, 1; β=β, μ=μ)),
		 (freebosons_keldysh_gfs(h, ts, 1, 2; β=β, μ=μ), freebosons_Gt(h, ts, 1, 2; β=β, μ=μ),
		  freebosons_keldysh_gfs(h, ts, 1, 1; β=β, μ=μ)))
		# canonical (anti)commutator contact: G>(0) - G<(0) = -i δᵢⱼ
		@test maximum(abs.(branch(gfs, :-, :+)[a, a] - branch(gfs, :+, :-)[a, a]) for a in 1:Nt) < 1.0e-12
		@test maximum(abs.(branch(gfs11, :-, :+)[a, a] - branch(gfs11, :+, :-)[a, a] + im) for a in 1:Nt) < 1.0e-12
		# retarded function from the diagonal blocks: for t > 0,
		# Gᵀ(t) - Gᵀ̃(t) = G>(t) - G<(t) = Gᴿ(t), and the sector rule makes both
		# sides exact rearrangements of the stored cross blocks
		@test maximum(abs.(branch(gfs, :+, :+)[a, 1] - branch(gfs, :-, :-)[a, 1] -
						   (branch(gfs, :-, :+)[a, 1] - branch(gfs, :+, :-)[a, 1])) for a in 2:Nt) < 1.0e-12
		# both statistics' Gt references are the physical retarded function
		@test maximum(abs.(branch(gfs, :+, :+)[a, 1] - branch(gfs, :-, :-)[a, 1] - Gt[a]) for a in 2:Nt) < 1.0e-12
		# keyword form of branch
		@test branch(gfs; b1=:+, b2=:+) == branch(gfs, :+, :+)
		# sizes and equality
		@test size(gfs) == (Nt, Nt)
		@test gfs == deepcopy(gfs)
		@test branches(gfs) == (:+, :-)
		@test_throws ArgumentError branch(gfs, :τ, :+)
		@test_throws ArgumentError branch(gfs, :+, :x)
	end

	# Kadanoff: mixed blocks are shared by the backward branch; contour
	# continuity ties the τ=0 mixed blocks to the real-time Wightman functions
	for gfs in (freefermions_kadanoff_gfs(h, ts6, τs6, 1, 1; β=β, μ=μ),
				freebosons_kadanoff_gfs(h, ts6, τs6, 1, 1; β=β, μ=μ))
		Nτ, Nt2 = size(gfs)
		@test branch(gfs, :-, :τ) == branch(gfs, :+, :τ)
		@test branch(gfs, :τ, :-) == branch(gfs, :τ, :+)
		@test branches(gfs) == (:+, :-, :τ)
		@test branch(gfs, :τ, :τ) == gfs.Gττ
		# G₊τ(t, τ=0) is the greater function, Gτ₊(τ=0, t) the lesser at -t
		@test branch(gfs, :+, :τ)[:, 1] ≈ branch(gfs, :-, :+)[:, 1] atol = 1.0e-13
		@test branch(gfs, :τ, :+)[1, :] ≈ branch(gfs, :+, :-)[1, :] atol = 1.0e-13
	end

	# the Imag container shares the (τ,τ) block convention of the Kadanoff one
	@test freefermions_matsubara_gfs(h, τs6, 1, 2; β=β, μ=μ).data ==
		  branch(freefermions_kadanoff_gfs(h, ts, τs6, 1, 2; β=β, μ=μ), :τ, :τ)
	@test freebosons_matsubara_gfs(h, τs6, 1, 2; β=β, μ=μ).data ==
		  branch(freebosons_kadanoff_gfs(h, ts, τs6, 1, 2; β=β, μ=μ), :τ, :τ)
	for gm in (freefermions_matsubara_gfs(h, τs6, 1; β=β, μ=μ),
			   freebosons_matsubara_gfs(h, τs6, 1; β=β, μ=μ))
		@test branch(gm, :τ, :τ) === gm.data
		@test size(gm) == (length(τs6), length(τs6))
		@test branches(gm) == (:τ,)
		@test gm == deepcopy(gm)
		@test_throws ArgumentError branch(gm, :+, :τ)
		@test sprint(show, MIME"text/plain"(), gm) isa String
	end

	@test sprint(show, MIME"text/plain"(), freefermions_keldysh_gfs(h, ts, 1; β=β, μ=μ)) isa String
	@test sprint(show, MIME"text/plain"(), freebosons_kadanoff_gfs(h, ts, τs6, 1; β=β, μ=μ)) isa String
end

