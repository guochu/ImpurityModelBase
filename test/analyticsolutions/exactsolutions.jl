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