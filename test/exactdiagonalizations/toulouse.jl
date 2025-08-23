println("------------------------------------")
println("|         ED Toulouse Model        |")
println("------------------------------------")



@testset "Toulouse Green's functions, normal bath" begin
	# normal
	atol=1.0e-8

	β = 10
	δτ = 1
	τs = collect(0:δτ:β)
	δt = 0.7
	t = 3.5
	ts = collect(0:δt:t)

	# μ = 0
	dw = 0.5

	ϵ_d = -0.5

	for μ in (0.5, 0, -0.5)
		b = discretebath(fermionicbath(semicircular(t=1), β=β, μ=μ), δw=dw)
		model = Toulouse(b, ϵ_d = ϵ_d)
		@test num_sites(model) == 5

		ham = hamiltonian(model)
		# @test cmatrix(ham) ≈ cmatrix(model) atol=atol

		h = fermionoperator(ham)
		chemical = zero(h)
		for i in 2:num_sites(model)
			chemical .+= fermiondensityoperator(num_sites(model), i)
		end

		a = fermionaoperator(num_sites(model), 1)
		adag = a'

		h1 = h - μ * chemical
		cache1 = eigencache(h1)
			
		# Matsubara Green's function
		g1 = correlation_2op_1τ(h1, a, adag, τs, cache1, β=β)
		g2 = toulouse_Gτ(model, τs)

		@test g1 ≈ g2 atol=atol

		# equilibrium greater and lesser
		hh = hamiltonian(model, include_chemical=true)
		dm = thermodm(hh, β=β)
		@test dm ≈ thermodm(model) atol=atol
		g1 = -im .* correlation_2op_1t(h, a, adag, dm, ts, reverse = false)
		l1 = im .* correlation_2op_1t(h, adag, a, dm, ts, reverse = true)

		cdm = fermionicthermocdm(eigencache(cmatrix(hh)), β=β)
		@test cdm ≈ thermocdm(model) atol=atol
		@test normal_quadratic_obs(dm) ≈ cdm atol=atol
		g2, l2 = freefermions_greater_lesser(cmatrix(ham), cdm, ts, 1, 1)

		@test g1 ≈ g2 atol=atol
		@test l1 ≈ l2 atol=atol

		g3, l3 = toulouse_greater_lesser(model, ts)

		@test g1 ≈ g3 atol=atol
		@test l1 ≈ l3 atol=atol

		# nonequilibrium greater and lesser
		dm = separabledm(model)
		cdm = separablecdm(model)
		@test normal_quadratic_obs(dm) ≈ cdm atol=atol

		g1 = -im .* correlation_2op_1t(h, a, adag, dm, ts, reverse = false)
		l1 = im .* correlation_2op_1t(h, adag, a, dm, ts, reverse = true)
		
		g2, l2 = toulouse_neq_greater_lesser(model, ts)	

		@test g1 ≈ g2 atol=atol
		@test l1 ≈ l2 atol=atol
	end
end


@testset "Toulouse Green's functions, BCS bath" begin
	atol=1.0e-8

	β = 10
	δτ = 1
	τs = collect(0:δτ:β)
	δt = 0.7
	t = 3.5
	ts = collect(0:δt:t)

	dw = 1

	ϵ_d = -0.5

	for Δ in (0.4, 0.2+0.3*im)
		b = discretebath(bcsbath(fermionicbath(semicircular(t=1), β=β), Δ=Δ), δw=dw)
		model = Toulouse(b, ϵ_d = ϵ_d)
		@test num_sites(model) == 6

		ham = hamiltonian(model)
		@test cmatrix(ham) ≈ cmatrix(model) atol=atol

		h = fermionoperator(ham)

		a = fermionaoperator(num_sites(ham), 1)
		adag = a'

		cache = eigencache(h)
			
		# Matsubara Green's function
		g1 = correlation_2op_1τ(h, a, adag, τs, cache, β=β)
		g2 = toulouse_Gτ(model, τs)

		@test g1 ≈ g2 atol=atol

		# equilibrium greater and lesser
		hh = hamiltonian(model)
		dm = thermodm(hh, β=β)
		@test dm ≈ thermodm(model) atol=atol
		g1 = -im .* correlation_2op_1t(h, a, adag, dm, ts, reverse = false)
		l1 = im .* correlation_2op_1t(h, adag, a, dm, ts, reverse = true)

		cdm = fermionicthermocdm(eigencache(cmatrix(hh)), β=β)
		@test cdm ≈ thermocdm(model) atol=atol
		@test generic_quadratic_obs(dm) ≈ cdm atol=atol
		g2, l2 = freefermions_greater_lesser(cmatrix(ham), cdm, ts, 1, 1)

		@test g1 ≈ g2 atol=atol
		@test l1 ≈ l2 atol=atol

		g3, l3 = toulouse_greater_lesser(model, ts)

		@test g1 ≈ g3 atol=atol
		@test l1 ≈ l3 atol=atol

		# nonequilibrium greater and lesser
		dm = separabledm(model)
		cdm = separablecdm(model)
		@test generic_quadratic_obs(dm) ≈ cdm atol=atol

		g1 = -im .* correlation_2op_1t(h, a, adag, dm, ts, reverse = false)
		l1 = im .* correlation_2op_1t(h, adag, a, dm, ts, reverse = true)
		
		g2, l2 = toulouse_neq_greater_lesser(model, ts)	

		@test g1 ≈ g2 atol=atol
		@test l1 ≈ l2 atol=atol
	end
end
