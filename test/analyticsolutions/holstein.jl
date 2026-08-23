println("------------------------------------")
println("|        Holstein Model            |")
println("------------------------------------")

function holstein_operators(ϵ_d; ω₀=1, α₀=0.5, ω₁=1, α₁=1, d=100)
	p1 = spin_half_matrices()
	n̂, σ₊, σ₋, JW = p1["n"], p1["+"], p1["-"], -p1["z"]
	Is = one(n̂)
	p2 = boson(d=d)
	b̂, b̂′, n̂b = p2["a"], p2["adag"], p2["n"]
	Ib = one(b̂)
	# total Hamiltonian
	Himp = -ϵ_d * kron(n̂, kron(Is, Ib))
	Hbath0 = ω₀ * kron(Is, kron(Is, n̂b))
	Hbath1 = ω₁ * kron(Is, kron(n̂, Ib))
	Hhyb0 = sqrt(α₀) * kron(n̂, kron(Is, b̂′ + b̂))
	tmp = kron(σ₊, kron(JW*σ₋, Ib))
	Hhyb1 = sqrt(α₁) * (tmp + tmp')
	H = Himp + Hhyb0 + Hhyb1 + Hbath0 + Hbath1
	A, B = kron(σ₋, kron(Is, Ib)), kron(σ₊, kron(Is, Ib))

	# Himp = kron(Is - n̂, kron(Is, Ib))
	Hbath = ω₀ * kron(Is, n̂b) + ω₁ * kron(n̂, Ib)
	rhosys = Is - n̂

	return H, A, B,  Hbath, rhosys
end

function holstein_neq(ϵ_d; β=1, t=1, N=100, ω₀=1, α₀=0.5, ω₁=1, α₁=1, d=100)
	δt=t/N

	H, a, adag, Hbath, rhosys = holstein_operators(ϵ_d, ω₀=ω₀, α₀=α₀, ω₁=ω₁, α₁=α₁, d=d)
	if β == Inf
		λs, U = eigen(Hermitian(Hbath))
		gs = U[:, 1:1]
		rhobath = gs * gs'
	else
		rhobath = exp(-β*Hbath)
	end
	ρ = kron(rhosys, rhobath)
	g1, g2 = gf_real(H, a, adag, β, t, N, ρ)

	return g1, g2
end

@testset "Holstein model: real time" begin
	rtol = 1.0e-2

	δt=0.1
	N = 10
	t = N * δt

	ω₀=0.8
	α₀=0.5
	ω₁=1.2
	α₁=0.7

	spec = DiracDelta(ω=ω₁, α=α₁)

	for β in (10., Inf)
		for ϵ_d in (-0.5, 0., 0.7)

			g1 = [holstein_Gt(spec, tj, β=β, ϵ_d=-ϵ_d, g=sqrt(α₀), ω=ω₀) for tj in 0:δt:t]
			# println(g1)

			g1′, g2′ = holstein_neq(ϵ_d, β=β, t=t, N=N, ω₀=ω₀, α₀=α₀, ω₁=ω₁, α₁=α₁)
			# println(g1′)
			# println(g2′)

			# println("---------------")
			# println(g1 - g1′)

			@test norm(g1 - g1′) / norm(g1′) < rtol

		end
	end

	# for ϵ_d in (-0.5, 0., 0.7)

	# 	g1 = [holstein_Gt(spec, tj, β=β, ϵ_d=-ϵ_d, g=sqrt(α₀), ω=ω₀) for tj in 0:δt:t]

	# 	g1′, g2′ = holstein_neq(-ϵ_d, β=β, t=t, N=N, ω₀=ω₀, α₀=α₀, ω₁=ω₁, α₁=α₁)

	# 	@test norm(g1 - g1′) / norm(g1) < tol

	# end
	
end