println("------------------------------------")
println("|       noninteracting fermions    |")
println("------------------------------------")

spectrum_func() = DiracDelta(ω=1, α=0.5)


@testset "Noninteracting fermions: imaginary time" begin
	tol = 1.0e-4
	δτ=0.01
	N = 10
	β = N * δτ
	for ϵ_d in (-0.5, 0., 0.7)
		g1 = independentbosons_Gτ(spectrum_func(), β=β, ϵ_d=-ϵ_d, Nτ=N)
		H, a, adag = noninteracting_operators(ϵ_d, ω₀=1, α=0.5, d=100)
		g2 = gf_imag(H, a, adag, β, N)
		@test norm(g1 - g2) / norm(g1) < tol
		g2′ = correlation_2op_1τ(H, a, adag, 0:δτ:β, β=β)
		@test norm(g1 - g2′) / norm(g1) < tol
	end
end


@testset "Noninteracting fermions: real time" begin
	tol = 1.0e-4

	β = 0.2
	δt=0.01
	N = 10
	t = N * δt
	for ϵ_d in (-0.5, 0., 0.7)

		g1 = [independentbosons_greater(spectrum_func(), tj, β=β, ϵ_d=-ϵ_d) for tj in 0:δt:t]
		g2 = [independentbosons_lesser(spectrum_func(), tj, β=β, ϵ_d=-ϵ_d) for tj in 0:δt:t]

		H, a, adag = noninteracting_operators(ϵ_d, ω₀=1, α=0.5, d=100)
		g1′, g2′ = gf_real(H, a, adag, β, t, N)
		@test norm(g1 - g1′) / norm(g1) < tol
		@test norm(g2 - g2′) / norm(g2) < tol

		d1 = -im .* correlation_2op_1t(H, a, adag, exp(-β * H), 0:δt:t, reverse = false)
		d2 = im .* correlation_2op_1t(H, adag, a, exp(-β * H), 0:δt:t, reverse = true)

		@test norm(g1 - d1) / norm(g1) < tol
		@test norm(g2 - d2) / norm(g2) < tol		
	end
	
end


@testset "Noninteracting fermions with ρ_0: real time" begin
	tol = 1.0e-4
	β = 0.2
	δt = 0.01
	N = 10
	t = N * δt
	d = 100
	# ρ_0 = Diagonal([P(|0⟩), P(|1⟩)])
	for ϵ_d in (-0.5, 0., 0.7)
		for ρ_0 in (Diagonal([1.0, 0.0]), Diagonal([0.0, 1.0]), Diagonal([0.4, 0.6]))
			g1 = [independentbosons_greater(spectrum_func(), tj, ρ_0, β=β, ϵ_d=-ϵ_d) for tj in 0:δt:t]
			g2 = [independentbosons_lesser(spectrum_func(), tj, ρ_0, β=β, ϵ_d=-ϵ_d) for tj in 0:δt:t]

			H, a, adag = noninteracting_operators(ϵ_d, ω₀=1, α=0.5, d=d)
			ρ = prod_state(ρ_0, 1, β, d)
			d1 = -im .* correlation_2op_1t(H, a, adag, ρ, 0:δt:t, reverse = false)
			d2 = im .* correlation_2op_1t(H, adag, a, ρ, 0:δt:t, reverse = true)

			@test norm(g1 - d1) / max(norm(g1), eps()) < tol
			@test norm(g2 - d2) / max(norm(g2), eps()) < tol
		end
	end
end


@testset "Noninteracting fermions with full ρ_0 matrix: real time" begin
	tol = 1.0e-4
	β = 0.2
	δt = 0.01
	N = 10
	t = N * δt
	d = 100
	# arbitrary (non-diagonal) impurity density matrix; only its diagonal elements
	# can enter G</G> because the impurity occupation is conserved
	ρ_0 = [0.4  0.2+0.1im;
	       0.2-0.1im  0.6]
	for ϵ_d in (-0.5, 0., 0.7)
		g1 = [independentbosons_greater(spectrum_func(), tj, ρ_0, β=β, ϵ_d=-ϵ_d) for tj in 0:δt:t]
		g2 = [independentbosons_lesser(spectrum_func(), tj, ρ_0, β=β, ϵ_d=-ϵ_d) for tj in 0:δt:t]

		H, a, adag = noninteracting_operators(ϵ_d, ω₀=1, α=0.5, d=d)
		ρ = prod_state(ρ_0, 1, β, d)
		d1 = -im .* correlation_2op_1t(H, a, adag, ρ, 0:δt:t, reverse = false)
		d2 = im .* correlation_2op_1t(H, adag, a, ρ, 0:δt:t, reverse = true)

		@test norm(g1 - d1) / max(norm(g1), eps()) < tol
		@test norm(g2 - d2) / max(norm(g2), eps()) < tol
	end
end