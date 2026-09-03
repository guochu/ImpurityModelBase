println("------------------------------------")
println("|        Holstein Model            |")
println("------------------------------------")

# ---------------------------------------------------------------------------
# helper for exact-diagonalization benchmark (impurity real-time G(t))
# ---------------------------------------------------------------------------
function holstein_operators(ϵ_d; ω₀=1, α₀=0.5, ω₁=1, α₁=1, d=100)
	p1 = spin_half_matrices()
	n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
	Is = one(n̂)
	p2 = boson(d=d)
	b̂, b̂′, n̂b = p2["a"], p2["adag"], p2["n"]
	Ib = one(b̂)
	# total Hamiltonian
	Himp = -ϵ_d * kron(n̂, kron(Is, Ib))
	Hbath0 = ω₀ * kron(Is, kron(Is, n̂b))
	Hbath1 = ω₁ * kron(Is, kron(n̂, Ib))
	Hhyb0 = sqrt(α₀) * kron(n̂, kron(Is, b̂′ + b̂))
	tmp = kron(σ₊, kron(JW*σ₋, Ib))
	Hhyb1 = sqrt(α₁) * (tmp + tmp')
	H = Himp + Hhyb0 + Hhyb1 + Hbath0 + Hbath1
	A, B = kron(σ₋, kron(Is, Ib)), kron(σ₊, kron(Is, Ib))

	# Himp = kron(Is - n̂, kron(Is, Ib))
	Hbath = ω₀ * kron(Is, n̂b) + ω₁ * kron(n̂, Ib)
	rhosys = Is - n̂

	return H, A, B,  Hbath, rhosys
end

function holstein_neq(ϵ_d; β=1, t=1, N=100, ω₀=1, α₀=0.5, ω₁=1, α₁=1, d=100)
	δt=t/N

	H, a, adag, Hbath, rhosys = holstein_operators(ϵ_d, ω₀=ω₀, α₀=α₀, ω₁=ω₁, α₁=α₁, d=d)
	if β == Inf
		λs, U = eigen(Hermitian(Hbath))
		gs = U[:, 1:1]
		rhobath = gs * gs'
	else
		rhobath = exp(-β*Hbath)
	end
	ρ = kron(rhosys, rhobath)
	g1, g2 = gf_real(H, a, adag, β, t, N, ρ)

	return g1, g2
end

@testset "Holstein model: real time" begin
	rtol = 5.0e-2

	δt=0.1
	N = 5
	t = N * δt

	ω₀=0.8
	α₀=0.5
	ω₁=1.2
	α₁=0.7

	spec = DiracDelta(ω=ω₁, α=α₁)

	# lightened ED cross-check: small boson truncation, a reduced integration
	# window and a shallow CFE keep the reference calculation fast
	for β in (10., Inf)
		for ϵ_d in (-0.5, 0., 0.7)

			g1 = [holstein_Gt(spec, tj, β=β, ϵ_d=-ϵ_d, g=sqrt(α₀), ω=ω₀, wmax=4.0, maxiter=6) for tj in 0:δt:t]

			g1′, g2′ = holstein_neq(ϵ_d, β=β, t=t, N=N, ω₀=ω₀, α₀=α₀, ω₁=ω₁, α₁=α₁, d=8)

			@test norm(g1 - g1′) / norm(g1′) < rtol

		end
	end

end

# ---------------------------------------------------------------------------
# analytical solution: atomic limit (CFE must reduce to the exact atomic
# propagator). At T=0 the spectral density of the atomic model has peaks at
# ω_n = ω0 (n - α²) with Poisson weights e^{-α²} α^{2n}/n!, α = g/ω0.
# ---------------------------------------------------------------------------
@testset "Holstein model: atomic limit (T=0)" begin
	g = 1.0
	ω0 = 1.0
	α2 = g^2 / ω0^2
	δ = 1.0e-3
	G0w(y) = 1 / (y + im * δ)          # atomic free propagator

	# with a finite δ each δ-function peak is broadened to height w/(πδ)
	for n in 0:4
		wp = ω0 * (n - α2)
		A = -imag(holstein_G0w_to_Gw(G0w, wp; β=Inf, g=g, ω=ω0, maxiter=60)) / π
		exact = exp(-α2) * α2^n / factorial(n) / (π * δ)
		@test A ≈ exact rtol = 1.0e-2
	end
end

# ---------------------------------------------------------------------------
# analytical solution: atomic limit at finite temperature. For large β the
# spectrum must reduce to the T=0 Poisson result; a small-β test checks the
# thermal redistribution of the zero-phonon peak.
# ---------------------------------------------------------------------------
@testset "Holstein model: atomic limit (finite T)" begin
	g = 1.0
	ω0 = 1.0
	α2 = g^2 / ω0^2
	δ = 1.0e-3
	G0w(y) = 1 / (y + im * δ)

	# large β → converges to the T=0 spectrum
	β = 100.0
	peaks = [ω0 * (n - α2) for n in 0:3]
	A_largeβ = [-imag(holstein_G0w_to_Gw(G0w, wp; β=β, g=g, ω=ω0, maxiter=40)) / π for wp in peaks]
	for (n, A) in enumerate(A_largeβ)
		exact = exp(-α2) * α2^(n - 1) / factorial(n - 1) / (π * δ)
		@test A ≈ exact rtol = 1.0e-2
	end

	# small β → thermal phonons deplete the zero-phonon peak (its integrated
	# weight W₀ is reduced compared to the low-temperature limit)
	β = 0.5
	ws = collect(range(-1.5, 3.0; length=1501))
	dw = ws[2] - ws[1]
	A_smallβ = [-imag(holstein_G0w_to_Gw(G0w, w; β=β, g=g, ω=ω0, maxiter=60)) / π for w in ws]
	W0_small = sum(A_smallβ[(peaks[1] - 0.5 .<= ws) .& (ws .<= peaks[1] + 0.5)]) * dw
	A_largeβ2 = [-imag(holstein_G0w_to_Gw(G0w, w; β=100.0, g=g, ω=ω0, maxiter=60)) / π for w in ws]
	W0_large = sum(A_largeβ2[(peaks[1] - 0.5 .<= ws) .& (ws .<= peaks[1] + 0.5)]) * dw
	@test W0_small < W0_large * 0.85
end

# ---------------------------------------------------------------------------
# DMFT on the Bethe lattice
# ---------------------------------------------------------------------------
@testset "Holstein model: DMFT on the Bethe lattice" begin

	# non-interacting limit (λ=0) must reproduce the Bethe semicircle. The
	# finite broadening δ smears the square-root singularities at the band edges
	# ω=±1, so the pointwise check is restricted to the smooth interior; the
	# normalization is verified over the whole band.
	ws = collect(range(-1.5, 1.5; length=301))
	r0 = holstein_dmft_bethe(ws; λ=0.0, γ=2.0, δ=5.0e-3, nw=1200, maxit=200, tol=1.0e-5)
	@test r0.converged
	interior = abs.(r0.ws) .<= 0.95
	Aexact = 2 .* sqrt.(max.(1 .- r0.ws[interior] .^ 2, 0.0)) ./ π
	@test maximum(abs.(r0.A[interior] .- Aexact)) < 0.02
	dw0 = r0.ws[2] - r0.ws[1]
	@test sum(r0.A) * dw0 ≈ 1 rtol = 0.05

	# finite coupling: converged, normalized, with a low-energy coherent peak
	ws = collect(range(-4.0, 4.0; length=401))
	r = holstein_dmft_bethe(ws; λ=0.75, γ=2.0, δ=5.0e-3, nw=1400, maxit=300, tol=1.0e-5)
	@test r.converged
	dw = r.ws[2] - r.ws[1]
	@test sum(r.A) * dw ≈ 1 rtol = 0.05
	# coherent polaron peak in the low-energy window
	low = r.ws .< 0.5
	@test maximum(r.A[low]) > 0.3
	# spectral density is non-negative
	@test minimum(r.A) > -0.02
end
