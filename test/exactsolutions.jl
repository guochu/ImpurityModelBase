println("------------------------------------")
println("|        Exact Solutions           |")
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

		g₁ = toulouse_Gτ(b2, τs, ϵ_d=ϵ_d)
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
		gf1 = [toulouse_Gt(spec, tj, ϵ_d = ϵ_d, wmax=100) for tj in ts]
		gf2 = toulouse_Gt(b2, ts, ϵ_d = ϵ_d)
		@test norm(gf1 - gf2) / norm(gf1) < rtol

	end
end