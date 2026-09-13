println("------------------------------------")
println("|      Analytical Solutions        |")
println("------------------------------------")



@testset "GF-imaginary time: benchmarking ED with Analytic solutions" begin
	N = 25
	δτ = 0.01
	ϵ_d = 1.25*pi
	dw = 0.1
	β = N * δτ
	τs = collect(0:δτ:β)
	rtol = 1.0e-2

	spec = spectrum_func()

	for μ in (-5, 0, 5)

		b1 = bath(Fermion, spec, β=β, μ=μ)
		b2 = discretebath(b1, δw=dw)

		g₁ = toulouse_Gτ(Toulouse(b2, ϵ_d=ϵ_d), τs)
		g₂ = [real(toulouse_Gτ(b1, τ, ϵ_d = ϵ_d)) for τ in τs]

		@test norm(g₁ - g₂) / norm(g₁) < rtol
	end

end


@testset "GF-imaginary time (bosons): benchmarking ED with Analytic solutions" begin
	N = 25
	δτ = 0.01
	ϵ_d = 1.25*pi
	dw = 0.1
	β = N * δτ
	τs = collect(0:δτ:β)
	rtol = 1.0e-2

	spec = spectrum_func()

	# for a (normal) bosonic Toulouse model the chemical potential must be negative to
	# keep the single-particle spectrum positive (required for a stable thermal state)
	for μ in (-5, -2)

		b1 = bath(Boson, spec, β=β, μ=μ)
		b2 = discretebath(b1, δw=dw)

		g₁ = toulouse_Gτ(Toulouse(b2, ϵ_d=ϵ_d), τs)
		g₂ = [real(toulouse_Gτ(b1, τ, ϵ_d = ϵ_d)) for τ in τs]

		@test norm(g₁ - g₂) / norm(g₁) < rtol
	end

end


@testset "GF-real time: benchmarking ED with Analytic solutions" begin
	N = 10
	δt = 0.01
	t = N * δt
	ϵ_d = 1.25*pi
	dw = 0.1
	β = 1.
	ts = [i*δt for i in 0:N]
	rtol = 1.0e-2
		

	# println("μ = ", μ)
	for spec in (spectrum_func(),)
		b1 = bath(Fermion, spec, β=β, μ=0.)
		b2 = discretebath(b1, δw=dw)
		gf1 = [toulouse_Gt(b1, tj, ϵ_d = ϵ_d, wmax=100) for tj in ts]
		gf2 = toulouse_Gt(Toulouse(b2, ϵ_d = ϵ_d), ts)
		@test norm(gf1 - gf2) / norm(gf1) < rtol

	end
end

@testset "Analytical solutions: closed forms" begin
	β = 2.0
	# free fermion
	@test freefermion_occupation(β, 0.5) ≈ 1 / (1 + exp(-β * 0.5))
	@test freefermion_Gτ(0.0; β=β, μ=0.5) ≈ 1 - freefermion_occupation(β, 0.5)
	t = 0.7
	@test freefermion_Gt(t; β=β, μ=0.5) ≈ freefermion_greater(t; β=β, μ=0.5) + freefermion_lesser(t; β=β, μ=0.5)

	# interacting fermion site
	U = 0.8
	@test fermion_Gt(t; β=β, μ=-0.3, U=U) ≈ fermion_greater(t; β=β, μ=-0.3, U=U) + fermion_lesser(t; β=β, μ=-0.3, U=U)

	# free boson
	ω = 1.0
	@test freeboson_occupation(β, ω) ≈ boseeinstein(β, ω)
	@test freeboson_Gτ(0.3; β=β, ω=ω) ≈ -exp(-ω * 0.3) / (1 - exp(-β * ω))
	@test freeboson_Gt(t; β=β, ω=ω) ≈ freeboson_greater(t; β=β, ω=ω) + freeboson_lesser(t; β=β, ω=ω)

	# Toulouse model with a single-mode (delta) bath: analytic checks
	# Δ(iω) = α/(iω - ω0),  Δ(τ) = -α (1+e^{-βω0}) e^{ω0 τ},
	# G(ω) = 1/(ω+iδ-ϵ_d-Δ(ω))
	ω0 = 1.0
	α = 1.0
	δd = DiracDelta(ω0, α=α)
	n = 200
	Δiw = toulouse_Δiw(δd; β=β, n=n)
	@test Δiw ≈ [α / (im * w - ω0) for w in ifrequencies(β, n)]

	Nτ = 100
	δτ = β / Nτ
	Δτ = toulouse_Δτ(δd; β=β, Nτ=Nτ)
	@test length(Δτ) == Nτ + 1
	for i in 0:2
		@test Δτ[i + 1] ≈ -α * (1 + exp(-β * ω0)) * exp(ω0 * i * δτ)
	end

	# toulouse_Gw: decoupled bath -> free-impurity pole,
	# single delta mode -> closed analytic form
	ϵ_d = 0.2
	δ = 1.0e-8
	gspec = spectrum(ϵ -> 0.0, -1, 1)
	@test toulouse_Gw(fermionicbath(gspec, β=β, μ=0.0), 0.3; ϵ_d=ϵ_d, δ=δ) ≈ 1 / (0.3 + im * δ - ϵ_d)
	@test toulouse_Gw(fermionicbath(δd, β=β), 0.3; ϵ_d=ϵ_d, δ=δ) ≈ 1 / (0.3 + im * δ - ϵ_d - α / (0.3 - ω0 + im * δ))
end

@testset "Toulouse model: semicircular bath closed forms" begin
	# --- G(ω): closed form vs numerical Hilbert transform of semicircular(t) ---
	for (ϵ_d, μ, t) in ((0.3, 0.0, 1.0), (-0.7, 0.4, 1.0), (0.5, -0.3, 2.0), (2.3, 0.0, 1.0), (-2.5, 0.6, 2.0))
		spec = semicircular(t)
		bath = fermionicbath(spec, β=1.0, μ=μ)
		for ω in (-3.0, -1.5, -0.7, 0.0, 0.5, 1.2, 2.5)
			@test fermionic_toulouse_Gw_semicircular(ω; ϵ_d=ϵ_d, μ=μ, t=t) ≈
				toulouse_Gw(bath, ω; ϵ_d=ϵ_d) rtol = 1.0e-8
		end
	end

	# --- G(τ): closed form vs numerical Fourier transform of the generic G(ω) ---
	τs = collect(0.1:0.2:2.0)
	for (ϵ_d, μ, t) in ((0.3, 0.0, 1.0), (-0.7, 0.4, 1.0), (2.3, 0.0, 1.0), (-2.5, 0.6, 2.0), (0.0, 0.0, 1.6), (0.0, 0.0, 2.5))
		spec = semicircular(t)
		bath = fermionicbath(spec, β=1.0, μ=μ)
		g1 = [fermionic_toulouse_Gt_semicircular(τ; ϵ_d=ϵ_d, μ=μ, t=t) for τ in τs]
		g2 = [toulouse_Gt(bath, τ; ϵ_d=ϵ_d) for τ in τs]
		# the generic reference carries finite-window/δ/quadrature errors (its
		# quadgk must resolve the near-δ bound-state Lorentzians), hence 2e-2
		@test norm(g1 - g2) / norm(g2) < 2.0e-2
	end

	# --- exact known limit: ϵ_d = μ = 0, t = 2 (where 𝒢 coincides with Δ) →
	#     G(τ) = -iθ(τ) 2J₁(tτ)/(tτ); for t ≠ 2 no such simple form exists ---
	besselj1(x) = sum((-1.0)^k * (x / 2)^(2k + 1) / (factorial(big(k)) * factorial(big(k + 1))) for k in 0:30)
	for τ in (0.1, 0.7, 1.5, 2.5)
		@test fermionic_toulouse_Gt_semicircular(τ; ϵ_d=0.0, t=2.0) ≈
			-im * 2 * besselj1(2.0 * τ) / (2.0 * τ) rtol = 1.0e-8
	end

	# --- sum rule G(0⁺) = -i, causality, and the bound-state poles ---
	@test fermionic_toulouse_Gt_semicircular(1.0e-12; ϵ_d=0.7, μ=0.3, t=1.0) ≈ -im rtol = 1.0e-6
	@test fermionic_toulouse_Gt_semicircular(-0.4; ϵ_d=0.7, t=1.0) == 0.0im
	# single bound state far outside the band: for t=2, r = ϵ_d + 1/ϵ_d, Z = t²/((t²-2)+2r/√(r²-t²));
	# at large τ the pole dominates the vanishing cut contribution
	t, ϵ_d = 2.0, 3.0
	r = ϵ_d + 1 / ϵ_d
	s = sqrt(r^2 - t^2)
	Z = t^2 / ((t^2 - 2) + 2 * r / s)
	τ = 15.0
	@test fermionic_toulouse_Gt_semicircular(τ; ϵ_d=ϵ_d, t=t) ≈ -im * Z * cis(-r * τ) atol = 3.0e-2 * Z
	# narrow band (t < √2): a particle-hole pair of bound states even at ϵ_d = μ = 0,
	# at r = ±2/√(4-t²) with Z = t²/((t²-2)+2r/√(r²-t²)) each; the cut is negligible at large τ
	t = 1.0
	r = 2 / sqrt(4 - t^2)
	s = sqrt(r^2 - t^2)
	Z = t^2 / ((t^2 - 2) + 2 * r / s)
	τ = 15.0
	@test fermionic_toulouse_Gt_semicircular(τ; ϵ_d=0.0, t=t) ≈ -im * Z * (cis(-r * τ) + cis(r * τ)) atol = 3.0e-2 * 2Z

	# --- Matsubara G(iω): closed form vs numerical Hilbert transform of semicircular(t) ---
	for (ϵ_d, μ, t) in ((0.3, 0.0, 1.0), (-0.7, 0.4, 1.0), (2.3, 0.0, 1.0), (-2.5, 0.6, 2.0), (0.0, 0.0, 1.6))
		spec = semicircular(t)
		bath = fermionicbath(spec, β=2.0, μ=μ)
		for ω in (0.4, 1.7, -3.1, 9.3)  # ω=0 avoided: the generic integral is PV-divergent for |μ|<t
			@test fermionic_toulouse_Giw_semicircular(ω; ϵ_d=ϵ_d, μ=μ, t=t) ≈
				toulouse_Giw(bath, ω; ϵ_d=ϵ_d) rtol = 1.0e-8
		end
	end

	# --- imaginary-time G(τ): closed form vs generic Matsubara sum (interior points) ---
	β = 2.0
	τs = (0.2, 0.6, 1.0, 1.4, 1.8)
	for (ϵ_d, μ, t) in ((0.3, 0.0, 1.0), (-0.7, 0.4, 1.0), (2.3, 0.0, 1.0), (-2.5, 0.6, 2.0), (0.0, 0.0, 1.6))
		spec = semicircular(t)
		bath = fermionicbath(spec, β=β, μ=μ)
		g1 = [fermionic_toulouse_Gτ_semicircular(τ; β=β, ϵ_d=ϵ_d, μ=μ, t=t) for τ in τs]
		g2 = [toulouse_Gτ(bath, τ; ϵ_d=ϵ_d, n=2000) for τ in τs]
		# the generic Matsubara sum carries a slow O(1/n) tail error, hence 5e-3
		@test norm(g1 - g2) / norm(g2) < 5.0e-3
	end

	# --- particle-hole symmetry: at the PH-invariant point ϵ_d = μ = 0 → G(τ) = G(β-τ) ---
	for t in (1.0, 2.0)
		for τ in (0.3, 0.9, 1.4)
			@test fermionic_toulouse_Gτ_semicircular(τ; β=β, ϵ_d=0.0, μ=0.0, t=t) ≈
				fermionic_toulouse_Gτ_semicircular(β - τ; β=β, ϵ_d=0.0, μ=0.0, t=t) rtol = 1.0e-8
		end
	end
end