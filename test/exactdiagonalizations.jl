println("------------------------------------")
println("|      Exact Diagonalizations      |")
println("------------------------------------")


# utilities functions
function random_normalquadratichamiltonian(m::AbstractMatrix) 
	L = size(m, 1)
	ham = NormalQuadraticHamiltonian(eltype(m), L)
	for i in 1:L, j in 1:L
		t = adaga(i, j, coeff=m[i, j])
		push!(ham, t)
	end
	return ham
end

function random_genericquadratichamiltonian(m::AbstractMatrix, m2::AbstractMatrix) 
	L = size(m, 1)
	ham = GenericQuadraticHamiltonian(eltype(m), L)
	for i in 1:L, j in 1:L
		t = adaga(i, j, coeff=m[i, j])
		push!(ham, t)
		t = adagadag(i, j, coeff=m2[i, j])
		if m2[i, j] != 0
			push!(ham, t)
			push!(ham, t')	
		end
	end
	return ham
end

function random_dm(::Type{T}, L::Int) where {T<:Number}
	dm = randn(T, L, L)
	dm = dm * dm'

	dm ./= tr(dm)
	return dm
end

function random_hermitian(::Type{T}, L::Int) where {T<:Number}
	m = randn(T, L, L)
	return m + m'
end

function normal_quadratic_obs(dm::AbstractMatrix)
	L = round(Int, log2(size(dm, 1)))
	obs = zeros(eltype(dm), L, L)
	tr_dm = tr(dm)
	for i in 1:L, j in 1:L
		t = adaga(i, j)
		op = fermionoperator(L, t)
		obs[i, j] = tr(op * dm) / tr_dm
	end
	return obs
end

function generic_quadratic_obs(dm::AbstractMatrix)
	L = round(Int, log2(size(dm, 1)))
	obs = zeros(eltype(dm), L, L)
	tr_dm = tr(dm)
	for i in 1:L, j in 1:L
		t = adaga(i, j)
		op = fermionoperator(L, t)
		obs[i, j] = tr(op * dm) / tr_dm
	end
	obs2 = zeros(eltype(dm), L, L)
	obs3 = zeros(eltype(dm), L, L)
	for i in 1:L, j in 1:L
		t = adagadag(i, j)
		op = fermionoperator(L, t)
		obs2[i, j] = tr(op * dm) / tr_dm

		t = aa(i, j)
		op = fermionoperator(L, t)
		obs3[i, j] = tr(op * dm) / tr_dm
	end	
	(obs2 ≈ obs3') || throw(ArgumentError("something wrong"))
	return bcs_cdm(obs, obs2)	
end

@testset "thermal state" begin
	atol = 1.0e-8
	β = 1.4

	# normal
	for T in (Float64, ComplexF64)
		for μ in (0.5, 0, -0.5)
			for L in 1:4
				m = random_hermitian(T, L)

				h = random_normalquadratichamiltonian(m)
				@test m ≈ cmatrix(h) atol = atol
				dm = thermodm(h, β=β, μ=μ)
				cdm = fermionicthermocdm(eigencache(m), β=β, μ=μ)
				# println("cdm ", cdm)
				tr_dm = tr(dm)
				ob1 = zeros(T, L, L)
				ob2 = zeros(T, L, L)
				for i in 1:L, j in 1:L
					t = adaga(i, j)
					x1 = cmatrix(L, t)
					ob1[i, j] = sum(x1 .*  cdm)
					x2 = fermionoperator(L, t)
					ob2[i, j] = tr(x2 * dm) / tr_dm
				end

				@test ob1 ≈ ob2 atol=atol
			end
		end
	end

	# bcs
	for T in (Float64, ComplexF64)
		for L in 1:4
			m = random_hermitian(T, L)
			m2 = rand(T, L, L)
			h = random_genericquadratichamiltonian(m, m2)
			dm = thermodm(h, β=β)

			mm = bcs_cmatrix(m, m2)
			@test mm ≈ cmatrix(h) atol=atol
			cdm = fermionicthermocdm(eigencache(mm), β=β)

			tr_dm = tr(dm)

			ob1 = zeros(T, L, L)
			ob2 = zeros(T, L, L)
			for i in 1:L, j in 1:L
				t = adaga(i, j)
				x1 = cmatrix(L, t, normal=false)
				ob1[i, j] = sum(x1 .*  cdm)
				x2 = fermionoperator(L, t)
				ob2[i, j] = tr(x2 * dm) / tr_dm
			end
			@test ob1 ≈ ob2 atol=atol
			# fill!(ob1, 0)
			# fill!(ob2, 0)

			for i in 1:L, j in 1:L
				t = adagadag(i, j)
				x1 = cmatrix(L, t)
				ob1[i, j] = sum(x1 .*  cdm)
				x2 = fermionoperator(L, t)
				ob2[i, j] = tr(x2 * dm) / tr_dm
			end		
			@test ob1 ≈ ob2 atol=atol

			for i in 1:L, j in 1:L
				t = aa(i, j)
				x1 = cmatrix(L, t)
				ob1[i, j] = sum(x1 .*  cdm)
				x2 = fermionoperator(L, t)
				ob2[i, j] = tr(x2 * dm) / tr_dm
			end		
			@test ob1 ≈ ob2 atol=atol
		end
	end
end

@testset "real time evolution" begin
	atol=1.0e-8

	β = 10
	dt = 0.7
	n = 5

	# normal
	for T in (Float64, ComplexF64)
		for L in 1:4
			dm = random_dm(T, 2^L)

			cdm = normal_quadratic_obs(dm)

			# hamiltonian
			m = random_hermitian(T, L)

			ham = random_normalquadratichamiltonian(m)
			h = fermionoperator(ham)
			# time evolution
			cache1 = eigencache(h)
			cache2 = eigencache(transpose(m))
			for k in 1:n
				t = k * dt
				rho1 = timeevo(dm, h, -im*t, cache1)
				rho2_cdm = timeevo(cdm, cache2.m, im*t, cache2)
				rho1_cdm = normal_quadratic_obs(rho1)

				@test rho1_cdm ≈ rho2_cdm atol=atol
			end

		end
	end

	# bcs
	for T in (Float64, ComplexF64)
		for L in 1:4
			dm = random_dm(T, 2^L)

			cdm = generic_quadratic_obs(dm)

			# hamiltonian
			m = random_hermitian(T, L)
			m2 = rand(T, L, L)
			mm = bcs_cmatrix(m, m2)

			ham = random_genericquadratichamiltonian(m, m2)
			h = fermionoperator(ham)

			# time evolution
			cache1 = eigencache(h)
			cache2 = eigencache(transpose(mm))
			for k in 1:n
				t = k * dt
				rho1 = timeevo(dm, h, -im*t, cache1)
				rho2_cdm = timeevo(cdm, cache2.m, im*t, cache2)
				rho1_cdm = generic_quadratic_obs(rho1)

				@test rho1_cdm ≈ rho2_cdm atol=atol
			end

		end
	end
end

@testset "real-time Green's functions" begin
	# normal
	atol=1.0e-8

	β = 10
	t = 0.7
	δt = 0.1
	ts = collect(0:δt:t)

	for T in (Float64, ComplexF64)
		for L in 1:4
			# hamiltonian
			m = random_hermitian(T, L)
			ham = random_normalquadratichamiltonian(m)

			h = fermionoperator(ham)
			# time evolution
			cache1 = eigencache(h)
			cache2 = eigencache(m)


			dm = random_dm(T, 2^L)
			cdm = normal_quadratic_obs(dm)

			for i in 1:L, j in 1:L
				a_i = fermionaoperator(L, i)
				adag_j = fermionadagoperator(L, j)
				g1 = -im .* correlation_2op_1t(h, a_i, adag_j, dm, ts, cache1, reverse = false)
				l1 = im .* correlation_2op_1t(h, adag_j, a_i, dm, ts, cache1, reverse = true)

				g2, l2 = freefermions_greater_lesser(m, cdm, ts, i, j, cache2)

				@test g1 ≈ g2 atol=atol
				@test l1 ≈ l2 atol=atol
			end

			# equilibrium green's function
			for μ in (0.5, 0, -0.5)
				dm = thermodm(ham, β=β, μ=μ)
				cdm = normal_quadratic_obs(dm)

				for i in 1:L, j in 1:L
					a_i = fermionaoperator(L, i)
					adag_j = fermionadagoperator(L, j)
					g1 = -im .* correlation_2op_1t(h, a_i, adag_j, dm, ts, cache1, reverse = false)
					l1 = im .* correlation_2op_1t(h, adag_j, a_i, dm, ts, cache1, reverse = true)

					g2, l2 = freefermions_greater_lesser(m, cdm, ts, i, j, cache2)

					@test g1 ≈ g2 atol=atol
					@test l1 ≈ l2 atol=atol

					g3, l3 = freefermions_greater_lesser(m, ts, i, j, cache2, β=β, μ=μ)

					@test g1 ≈ g3 atol=atol
					@test l1 ≈ l3 atol=atol

				end
			end
		end
	end	

	# bcs 
	μ = 0
	for T in (Float64, ComplexF64)
		for L in 1:4
			# hamiltonian
			m = random_hermitian(T, L)
			m2 = rand(T, L, L)
			mm = bcs_cmatrix(m, m2)
			ham = random_genericquadratichamiltonian(m, m2)

			h = fermionoperator(ham)
			# time evolution
			cache1 = eigencache(h)
			cache2 = eigencache(mm)


			dm = random_dm(T, 2^L)
			cdm = generic_quadratic_obs(dm)

			for i in 1:L, j in 1:L
				a_i = fermionaoperator(L, i)
				adag_j = fermionadagoperator(L, j)

				g1 = -im .* correlation_2op_1t(h, a_i, adag_j, dm, ts, cache1, reverse = false)
				l1 = im .* correlation_2op_1t(h, adag_j, a_i, dm, ts, cache1, reverse = true)

				g2, l2 = freefermions_greater_lesser(mm, cdm, ts, i, j, cache2)

				@test g1 ≈ g2 atol=atol
				@test l1 ≈ l2 atol=atol

				adag_i = a_i'
				g3 = -im .* correlation_2op_1t(h, adag_i, adag_j, dm, ts, cache1, reverse = false)
				l3 = im .* correlation_2op_1t(h, adag_j, adag_i, dm, ts, cache1, reverse = true)

				g4, l4 = freefermions_greater_lesser(mm, cdm, ts, L+i, j, cache2)
				@test g3 ≈ g4 atol=atol
				@test l3 ≈ l4 atol=atol				
			end

			# equilibrium green's function
			dm = thermodm(ham, β=β, μ=μ)
			cdm = generic_quadratic_obs(dm)
			@test cdm ≈ fermionicthermocdm(cache2, β=β, μ=μ)

			for i in 1:L, j in 1:L
				a_i = fermionaoperator(L, i)
				adag_j = fermionadagoperator(L, j)
				g1 = -im .* correlation_2op_1t(h, a_i, adag_j, dm, ts, cache1, reverse = false)
				l1 = im .* correlation_2op_1t(h, adag_j, a_i, dm, ts, cache1, reverse = true)

				g2, l2 = freefermions_greater_lesser(mm, cdm, ts, i, j, cache2)

				@test g1 ≈ g2 atol=atol
				@test l1 ≈ l2 atol=atol

				g3, l3 = freefermions_greater_lesser(mm, ts, i, j, cache2, β=β, μ=μ)

				@test g1 ≈ g3 atol=atol
				@test l1 ≈ l3 atol=atol
			end

		end
	end	
end

@testset "Matsubara Green's functions" begin
	# normal
	atol=1.0e-8

	β = 10
	δτ = 1
	τs = collect(0:δτ:β)

	for T in (Float64, ComplexF64)
		for L in 1:4
			# hamiltonian
			m = random_hermitian(T, L)
			ham = random_normalquadratichamiltonian(m)

			h = fermionoperator(ham)
			# time evolution
			cache1 = eigencache(h)
			cache2 = eigencache(m)

			for i in 1:L, j in 1:L
				a_i = fermionaoperator(L, i)
				adag_j = fermionadagoperator(L, j)
				g1 = correlation_2op_1τ(h, a_i, adag_j, τs, cache1, β=β)

				g2 = freefermions_Gτ(m, τs, i, j, cache2, β=β)

				@test g1 ≈ g2 atol=atol
			end
		end
	end	

	# bcs 
	for T in (Float64, ComplexF64)
		for L in 1:4
			# hamiltonian
			m = random_hermitian(T, L)
			m2 = rand(T, L, L)
			mm = bcs_cmatrix(m, m2)
			ham = random_genericquadratichamiltonian(m, m2)

			h = fermionoperator(ham)
			# time evolution
			cache1 = eigencache(h)
			cache2 = eigencache(mm)
			for i in 1:L, j in 1:L
				a_i = fermionaoperator(L, i)
				adag_j = fermionadagoperator(L, j)

				g1 = correlation_2op_1τ(h, a_i, adag_j, τs, cache1, β=β)

				g2 = freefermions_Gτ(mm, τs, i, j, cache2, β=β)

				@test g1 ≈ g2 atol=atol
			end
		end
	end	

end
