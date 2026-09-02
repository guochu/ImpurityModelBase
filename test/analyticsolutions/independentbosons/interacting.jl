println("------------------------------------")
println("|         interacting fermions     |")
println("------------------------------------")


@testset "Interacting fermions: imaginary time" begin
	tol = 1.0e-4
	δτ=0.01
	N = 10
	β = N * δτ
	for U in (0, 1)
		for ϵ_d in (-0.5, 0., 0.7)
			g1 = independentbosons_Gτ(spectrum_func(), β=β, ϵ_d=-ϵ_d, Nτ=N, U=U, bands=2)
			H, a, adag = interacting_operators(U, ϵ_d, ω₀=1, α=0.5, d=100)
			g2 = gf_imag(H, a, adag, β, N)
			@test norm(g1 - g2) / norm(g1) < tol

			g2′ = correlation_2op_1τ(H, a, adag, 0:δτ:β, β=β)
			@test norm(g1 - g2′) / norm(g1) < tol
		end
	end
end


@testset "Interacting fermions: real time" begin
	tol = 1.0e-4

	β = 0.2
	δt=0.01
	N = 10
	t = N * δt
	for U in (0, 1)
		for ϵ_d in (-0.5, 0., 0.7)

			g1 = [independentbosons_greater(spectrum_func(), tj, β=β, ϵ_d=-ϵ_d, U=U, bands=2) for tj in 0:δt:t]
			g2 = [independentbosons_lesser(spectrum_func(), tj, β=β, ϵ_d=-ϵ_d, U=U, bands=2) for tj in 0:δt:t]

			H, a, adag = interacting_operators(U, ϵ_d, ω₀=1, α=0.5, d=100)
			g1′, g2′ = gf_real(H, a, adag, β, t, N)
			@test norm(g1 - g1′) / norm(g1) < tol
			@test norm(g2 - g2′) / norm(g2) < tol

			ρ = exp(-β * H)
			d1 = -im .* correlation_2op_1t(H, a, adag, ρ, 0:δt:t, reverse = false)
			d2 = im .* correlation_2op_1t(H, adag, a, ρ, 0:δt:t, reverse = true)

			@test norm(g1 - d1) / norm(g1) < tol
			@test norm(g2 - d2) / norm(g2) < tol

		end
	end
end


@testset "Interacting fermions with ρ_0: real time" begin
	tol = 1.0e-4
	β = 0.2
	δt = 0.01
	N = 10
	t = N * δt
	d = 100
	# ρ_0 = Diagonal([P(|0⟩), P(|↑⟩), P(|↓⟩), P(|↑↓⟩)])
	ρ_0s = (Diagonal([1.0, 0.0, 0.0, 0.0]), Diagonal([0.0, 1.0, 0.0, 0.0]), Diagonal([0.0, 0.0, 1.0, 0.0]),
			Diagonal([0.0, 0.0, 0.0, 1.0]), Diagonal([0.2, 0.3, 0.1, 0.4]))
	for U in (0, 1)
		for ϵ_d in (-0.5, 0., 0.7)
			for ρ_0 in ρ_0s
				g1 = [independentbosons_greater(spectrum_func(), tj, ρ_0, β=β, ϵ_d=-ϵ_d, U=U, bands=2) for tj in 0:δt:t]
				g2 = [independentbosons_lesser(spectrum_func(), tj, ρ_0, β=β, ϵ_d=-ϵ_d, U=U, bands=2) for tj in 0:δt:t]

				H, a, adag = interacting_operators(U, ϵ_d, ω₀=1, α=0.5, d=d)
				ρ = prod_state(ρ_0, 1, β, d)
				d1 = -im .* correlation_2op_1t(H, a, adag, ρ, 0:δt:t, reverse = false)
				d2 = im .* correlation_2op_1t(H, adag, a, ρ, 0:δt:t, reverse = true)

				@test norm(g1 - d1) / max(norm(g1), eps()) < tol
				@test norm(g2 - d2) / max(norm(g2), eps()) < tol
			end
		end
	end
end


@testset "Interacting fermions with full ρ_0 matrix: real time" begin
	tol = 1.0e-4
	β = 0.2
	δt = 0.01
	N = 10
	t = N * δt
	d = 100
	# arbitrary (non-diagonal) 4x4 impurity density matrix in the
	# [|0⟩, |↑⟩, |↓⟩, |↑↓⟩] basis; only diagonal elements can enter G</G>
	ρ_0 = [0.2   0.05+0.02im  0.03       0.01;
	       0.05-0.02im  0.3  0.04-0.01im  0.02;
	       0.03  0.04+0.01im  0.25  0.02+0.03im;
	       0.01  0.02  0.02-0.03im  0.25]
	for U in (0, 1)
		for ϵ_d in (-0.5, 0., 0.7)
			g1 = [independentbosons_greater(spectrum_func(), tj, ρ_0, β=β, ϵ_d=-ϵ_d, U=U, bands=2) for tj in 0:δt:t]
			g2 = [independentbosons_lesser(spectrum_func(), tj, ρ_0, β=β, ϵ_d=-ϵ_d, U=U, bands=2) for tj in 0:δt:t]

			H, a, adag = interacting_operators(U, ϵ_d, ω₀=1, α=0.5, d=d)
			ρ = prod_state(ρ_0, 1, β, d)
			d1 = -im .* correlation_2op_1t(H, a, adag, ρ, 0:δt:t, reverse = false)
			d2 = im .* correlation_2op_1t(H, adag, a, ρ, 0:δt:t, reverse = true)

			@test norm(g1 - d1) / max(norm(g1), eps()) < tol
			@test norm(g2 - d2) / max(norm(g2), eps()) < tol
		end
	end
end