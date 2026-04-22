println("------------------------------------")
println("|        dephasing dynamics        |")
println("------------------------------------")


function toy_dephasing_model(Δ; d=100, α=0.5, ω₀=1)
	z = [0.5 0; 0 -0.5]
	a = bosonaoperator(d=d)
	n = bosondensityoperator(d=d)

	H = Δ .* kron(z, one(n)) + ω₀ .* kron(one(z), n) + sqrt(α) .* kron(z, a+a')

	return H, ω₀ .* n
end

@testset "Real time dephasing dynamics" begin
	tol = 1.0e-5

	Δ = 0.25
	α = 0.5
	ω₀ = 0.7
	β = 1.3
	d = 50

	ts = [0.5, 1]

	m = randn(ComplexF64, 2, 2)
	m = m' * m
	ρ₀ = m ./ tr(m)

	H, Hbath = toy_dephasing_model(Δ, α=α, ω₀=ω₀, d=d)

	ρbath = exp(-β .* Hbath)
	ρbath ./= tr(ρbath)

	ρ = kron(ρ₀, ρbath)

	cache = eigencache(H)

	spec = DiracDelta(ω=ω₀, α=α)

	for t in ts
		ρout = timeevo(ρ, H, -im*t, cache)

		tmp = reshape(ρout, d, 2, d, 2)
		ρout1 = zeros(eltype(tmp), 2, 2)
		for i in 1:2, j in 1:2
			ρout1[i, j] = tr(tmp[:, i, :, j])
		end

		ρout2 = spinboson_dephasingdynamics(spec, t, ρ₀, β=β, Δ=Δ)

		@test norm(ρout1 - ρout2) / norm(ρout1) < tol
	end
end

@testset "Real time dephasing dynamics after DD (XX sequence)" begin
	tol = 1.0e-5

	Δ = 0.25
	α = 0.5
	ω₀ = 0.7
	β = 1.3
	d = 50

	δt = 0.1
	N = 10

	function ptrace(rho)
		tmp = reshape(rho, d, 2, d, 2)
		ρout1 = zeros(eltype(tmp), 2, 2)
		for i in 1:2, j in 1:2
			ρout1[i, j] = tr(tmp[:, i, :, j])
		end	
		return ρout1
	end

	m = randn(ComplexF64, 2, 2)
	m = m' * m
	ρ₀ = m ./ tr(m)

	H, Hbath = toy_dephasing_model(Δ, α=α, ω₀=ω₀, d=d)

	Ib = one(Hbath)
	xop = kron([0 1.; 1. 0], Ib)

	ρbath = exp(-β .* Hbath)
	ρbath ./= tr(ρbath)

	ρ = kron(ρ₀, ρbath)

	cache = eigencache(H)

	spec = DiracDelta(ω=ω₀, α=α)

	ρout = ρ
	for i in 1:N
		ρout = timeevo(ρout, H, -im*δt, cache)
		ρout = xop * ρout * xop'
	end

	ρout1 = ptrace(ρout)
	
	# analytic solution
	ρout2 = ddxx_spinboson_dephasingdynamics(spec, N, ρ₀, β=β, δt=δt)
	@test norm(ρout1 - ρout2) / norm(ρout1) < tol
end