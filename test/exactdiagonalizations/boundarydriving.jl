println("------------------------------------")
println("|        ED Boundary Driving       |")
println("------------------------------------")


@testset "BoundaryDriving Green's functions and currents" begin

	atol=1.0e-8

	βl = 10
	βr = 1
	μl = 0.5
	μr = -0.5

	δt = 0.6
	t = 1.8
	ts = collect(0:δt:t)

	# μ = 0
	dw = 1


	leftbath = discretebath(fermionicbath(semicircular(t=1), β=βl, μ=μl), δw=dw)
	rightbath = discretebath(fermionicbath(semicircular(t=1), β=βr, μ=μr), δw=dw)

	L = 2
	for T in (Float64, ComplexF64)

		hsys = random_hermitian(T, L)

		model = BoundaryDriving(hsys, leftbath, rightbath)
		@test num_sites(model) == 6

		ham = hamiltonian(model)
		h2 = cmatrix(ham)
		@test h2 ≈ cmatrix(model) atol=atol


		i = 1

		h1 = fermionoperator(ham)
		a = fermionaoperator(num_sites(model), i)
		adag = a'


		# nonequilibrium greater and lesser

		rho_sys = random_dm(T, 2^L)
		dm = fermionicseparabledm(model, rho_sys)
		cdm = separablecdm(model, normal_quadratic_obs(rho_sys))
		@test normal_quadratic_obs(dm) ≈ cdm atol=atol

		g1 = -im .* correlation_2op_1t(h1, a, adag, dm, ts, reverse = false)
		l1 = im .* correlation_2op_1t(h1, adag, a, dm, ts, reverse = true)
		
		g2, l2 = freefermions_greater_lesser(h2, cdm, ts, i, i)

		@test g1 ≈ g2 atol=atol
		@test l1 ≈ l2 atol=atol


		# currents
		cache1 = eigencache(h1)
		cache2 = eigencache(transpose(h2))

		rho1 = timeevo(dm, h1, -im*t, cache1)
		rho2_cdm = timeevo(cdm, cache2.m, im*t, cache2)
		rho1_cdm = normal_quadratic_obs(rho1)

		@test rho1_cdm ≈ rho2_cdm atol=atol


		ob_ham = leftparticlecurrent_hamiltonian(model)
		ob_cdm = leftparticlecurrent_cmatrix(model)
		@test cmatrix(ob_ham) ≈ ob_cdm atol=atol
		ob_m = fermionoperator(ob_ham)

		x = tr(ob_m * rho1)
		y = sum(ob_cdm .* rho2_cdm)
		@test x ≈ y atol=atol

		ob_ham = rightparticlecurrent_hamiltonian(model)
		ob_cdm = rightparticlecurrent_cmatrix(model)
		@test cmatrix(ob_ham) ≈ ob_cdm atol=atol
		ob_m = fermionoperator(ob_ham)

		x = tr(ob_m * rho1)
		y = sum(ob_cdm .* rho2_cdm)
		@test x ≈ y atol=atol

		ob_ham = leftheatcurrent_hamiltonian(model)
		ob_cdm = leftheatcurrent_cmatrix(model)
		@test cmatrix(ob_ham) ≈ ob_cdm atol=atol
		ob_m = fermionoperator(ob_ham)

		x = tr(ob_m * rho1)
		y = sum(ob_cdm .* rho2_cdm)
		@test x ≈ y atol=atol

		ob_ham = rightheatcurrent_hamiltonian(model)
		ob_cdm = rightheatcurrent_cmatrix(model)
		@test cmatrix(ob_ham) ≈ ob_cdm atol=atol
		ob_m = fermionoperator(ob_ham)

		x = tr(ob_m * rho1)
		y = sum(ob_cdm .* rho2_cdm)
		@test x ≈ y atol=atol
	end


end