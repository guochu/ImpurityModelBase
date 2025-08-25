println("------------------------------------")
println("|     Bogoliubov transformation    |")
println("------------------------------------")


@testset "Bogoliubov transformation: Toulouse Green's functions" begin
	atol=1.0e-8

	β = 10
	δτ = 1
	τs = collect(0:δτ:β)
	δt = 0.7
	t = 3.5
	ts = collect(0:δt:t)

	dw = 1

	ϵ_d = 0.5

	Δ = 0.2+0.3*im
	spec = semicircular(t = 1)
	# spec = Leggett()
	for Δ in (0.4, 0.2+0.3*im)
		# println("Δ = ", Δ)
		b = discretebath(bcsbath(fermionicbath(spec, β=β), Δ=Δ), δw=dw)
		model = Toulouse(b, ϵ_d = ϵ_d)
		# @test num_sites(model) == 6
			
		# Matsubara Green's function
		g1 = toulouse_Gτ(model, τs)

		bg_h = bogoliubov_cmatrix(model)
		g2 = freefermions_Gτ(bg_h, τs, 1, 1, β=β)
		@test g1 ≈ g2 atol=atol

		# equilibrium greater and lesser
		g1, l1 = toulouse_greater_lesser(model, ts)
		bg_cdm = bogoliubov_thermocdm(model)
		g2, l2 = freefermions_greater_lesser(bg_h, bg_cdm, ts, 1, 1)

		@test g1 ≈ g2 atol=atol
		@test l1 ≈ l2 atol=atol

		cdm = thermocdm(model)
		g4, l4 = freefermions_greater_lesser(cmatrix(model), cdm, ts, num_sites(model)+1, 2)

		g5, l5 = freefermions_greater_lesser(bg_h, bg_cdm, ts, num_sites(model)+1, 2)

		@test g4 ≈ g5 atol=atol
		@test l4 ≈ l5 atol=atol

		# nonequilibrium greater and lesser
		
		g1, l1 = toulouse_neq_greater_lesser(model, ts)	

		bg_cdm = bogoliubov_separablecdm(model)		
		g2, l2 = freefermions_greater_lesser(bg_h, bg_cdm, ts, 1, 1)
		
		@test g1 ≈ g2 atol=atol
		@test l1 ≈ l2 atol=atol
	end

end