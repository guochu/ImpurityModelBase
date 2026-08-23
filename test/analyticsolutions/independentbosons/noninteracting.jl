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