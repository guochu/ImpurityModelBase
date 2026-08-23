# Tests for the exported exact-diagonalization operators and terms
# (src/exactdiagonalizations/operators.jl) and the Lindblad super-operator,
# not exercised elsewhere in the test suite. The (efficient) coefficient-matrix
# method is cross-checked against the (inefficient, debug-only)
# full-Hamiltonian method where applicable.

println("------------------------------------")
println("|     ED operators and terms       |")
println("------------------------------------")

@testset "Hamiltonian terms" begin
	# tunneling / adaga
	t = tunneling(1, 2, coeff=2.0)
	@test t isa AdagATerm
	@test positions(t) == (1, 2)
	@test adaga(1, 2, coeff=2.0) == t
	@test t' == AdagATerm(2, 1, coeff=2.0)
	@test (2.0 * t) isa AdagATerm

	# adagadag: scalar operations keep the term type (regression)
	a = adagadag(1, 2, coeff=1.0)
	@test a isa AdagAdagTerm
	@test (2.0 * a) isa AdagAdagTerm
	@test (-a) isa AdagAdagTerm
	@test copy(a) isa AdagAdagTerm
	@test (a / 2.0) isa AdagAdagTerm
	@test a' isa AATerm

	# aa
	@test aa(1, 2) isa AATerm
	@test aa(1, 2, coeff=1.5)' isa AdagAdagTerm

	# interaction / QuarticTerm
	q = interaction(1, 2, 3, 4, coeff=1.5)
	@test q isa QuarticTerm
	@test positions(q) == (1, 2, 3, 4)
	@test q' isa QuarticTerm

	# abstract / union type relations (concrete instantiations, since bare
	# UnionAlls with different parameter bounds do not compare with <:)
	@test AdagATerm{Float64} <: QuadraticTerm{Float64}
	@test AdagAdagTerm{Float64} <: QuadraticTerm{Float64}
	@test AATerm{Float64} <: QuadraticTerm{Float64}
	@test QuadraticTerm{Float64} <: AbstractTerm{Float64}
	@test QuarticTerm{Float64} <: AbstractTerm{Float64}
	@test AdagATerm{Float64} <: NormalTerm{Float64}
	@test QuarticTerm{Float64} <: NormalTerm{Float64}

	# quadratichamiltonian selects the Hamiltonian type
	h1 = quadratichamiltonian(2, [adaga(1, 2, coeff=1.0)])
	@test h1 isa NormalQuadraticHamiltonian
	h2 = quadratichamiltonian(2, [adagadag(1, 2, coeff=1.0)])
	@test h2 isa GenericQuadraticHamiltonian

	# NormalHamiltonian with bounds checking on push!
	nh = NormalHamiltonian(Float64, 2)
	@test length(nh.data) == 0
	push!(nh, adaga(1, 2, coeff=1.0))
	@test length(nh.data) == 1
	@test num_sites(nh) == 2
	@test_throws BoundsError push!(nh, adaga(3, 1, coeff=1.0))
end

@testset "Occupation operators and free Gt consistency" begin
	# fermion occupation projectors
	p1 = fermionoccupationoperator(1)
	@test p1 == [0 0; 0 1]
	p0 = fermionoccupationoperator(0)
	@test p0 == [1 0; 0 0]
	@test p1 + p0 == [1 0; 0 1]
	@test_throws ArgumentError fermionoccupationoperator(2)

	p = fermionoccupationoperator(2, 1, 1)
	@test size(p) == (4, 4)
	pmulti = fermionoccupationoperator([1, 0])
	@test pmulti == kron([0 0; 0 1], [1 0; 0 0])

	# Gt is G> - G< for both fermions and bosons
	L = 3
	h = quadratichamiltonian(L, [adaga(i, i, coeff=0.3 * (i - 1)) for i in 1:L])
	hc = cmatrix(h)
	β = 10.0
	t = 0.7

	f = freefermions_Gt(hc, 1, 2; β=β)
	gl = freefermions_greater_lesser(hc, 1, 2; β=β)   # single function gl(t) -> (G>, G<)
	g, l = gl(t)
	@test f(t) ≈ g - l

	fb = freebosons_Gt(hc, 1, 1; β=β, μ=-1.0)
	glb = freebosons_greater_lesser(hc, 1, 1; β=β, μ=-1.0)
	gb, lb = glb(t)
	@test fb(t) ≈ gb - lb
end

@testset "Lindblad operator" begin
	# fixed small example
	H = [1.0 0.2; 0.2 0.5]
	J = [0.0 0.0; 0.1 0.0]
	L = lindbladoperator(H, [J])
	@test L isa LindbladOperator
	@test size(L.m) == (2, 2, 2, 2)

	ρ = [0.7 0.1; 0.1 0.3]
	rhs = -im * (H * ρ - ρ * H) + 2 * (J * ρ * J') - (J' * J * ρ + ρ * J' * J)

	Lρ = zeros(ComplexF64, 2, 2)
	for i in 1:2, j in 1:2, k in 1:2, l in 1:2
		Lρ[i, j] += L.m[i, j, k, l] * ρ[k, l]
	end
	@test Lρ ≈ rhs atol=1.0e-12
end
