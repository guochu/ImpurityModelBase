println("------------------------------------")
println("|       Lindblad Correlation       |")
println("------------------------------------")

@testset "Lindblad steady state" begin
	tol = 1.0e-10
	for d in [2, 3]
		H = randn(ComplexF64, d, d)
		H = H + H'
		jumps = [randn(ComplexF64, d, d) for i in 1:d]
		L = lindbladoperator(H, jumps)
		rho = steady_state(L)
		@test maximum(abs.(rho - rho')) < tol
		rho′ = L(rho)
		@test maximum(abs.(rho′)) < tol
	end
end

@testset "Correlations" begin
	tol = 1.0e-6

	ts = [0, 0.3, 0.7, 1]
	for d in [2,3]
		H = randn(ComplexF64, d, d)
		H = H + H'
		L = lindbladoperator(H, [])
		rho = random_dm(ComplexF64, d)

		rho1 = timeevo(rho, H, -im)
		rho2 = timeevo(rho, L, 1)
		@test norm(rho1 - rho2) / norm(rho1) < tol

		A = randn(ComplexF64, d, d)
		B = randn(ComplexF64, d, d)
		for reverse in (true, false)
			corr1 = correlation_2op_1t(H, A, B, rho, ts, reverse=reverse)
			corr2 = correlation_2op_1t(L, A, B, rho, ts, reverse=reverse)

			@test norm(corr1 - corr2) / norm(corr1) < tol
		end
	end

end
