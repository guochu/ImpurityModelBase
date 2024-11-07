println("------------------------------------")
println("|       Free fermionic GF          |")
println("------------------------------------")

function freefermion_operators(ϵ_d)
	p1 = spin_half_matrices()
	n̂, σ₊, σ₋ = p1["n"], p1["+"], p1["-"]

	H = -ϵ_d * n̂

	return H, σ₋, σ₊
end


function fermion_operators(U, ϵ_d=U/2)
	p1 = spin_half_matrices()
	n̂, σ₊, σ₋ = p1["n"], p1["+"], p1["-"]
	Is = one(n̂)
	n_ud = kron(n̂, Is) + kron(Is, n̂)
	nn = kron(n̂,n̂)

	H = -ϵ_d*n_ud + U * nn
	A, B = kron(σ₋, Is), kron(σ₊, Is)

	return H, A, B
end



@testset "Noninteracting fermions: imaginary time" begin
	tol = 1.0e-4
	δτ=0.01
	N = 10
	β = N * δτ
	for ϵ_d in (-0.5, 0., 0.7)
		g1 = [freefermion_Gτ(τ, β=β, μ=-ϵ_d) for τ in 0:δτ:β]
		H, a, adag = freefermion_operators(-ϵ_d)
		g2 = gf_imag(H, a, adag, β, N)
		@test norm(g1 - g2) / norm(g1) < tol
	end
end


@testset "Noninteracting fermions: real time" begin
	tol = 1.0e-4
	δt=0.01
	N = 10
	t = N * δt
	β = 1.
	for ϵ_d in (-0.5, 0., 0.7)
		g1 = [freefermion_greater(τ, β=β, μ=-ϵ_d) for τ in 0:δt:t]
		g2 = [freefermion_lesser(τ, β=β, μ=-ϵ_d) for τ in 0:δt:t]
		H, a, adag = freefermion_operators(-ϵ_d)
		g1′, g2′ = gf_real(H, a, adag, β, t, N)
		@test norm(g1 - g1′) / norm(g1) < tol
		@test norm(g2 - g2′) / norm(g2) < tol
	end
end


@testset "Interacting fermions: imaginary time" begin
	tol = 1.0e-4
	δτ=0.01
	N = 10
	β = N * δτ
	for ϵ_d in (-0.5, 0., 0.7)
		for U in (0, 1)
			g1 = [fermion_Gτ(τ, β=β, μ=-ϵ_d, U=U) for τ in 0:δτ:β]
			H, a, adag = fermion_operators(U, -ϵ_d)
			g2 = gf_imag(H, a, adag, β, N)
			@test norm(g1 - g2) / norm(g1) < tol
		end
	end
end


@testset "Interacting fermions: real time" begin
	tol = 1.0e-4
	δt=0.01
	N = 10
	t = N * δt
	β = 1.
	for ϵ_d in (-0.5, 0., 0.7)
		for U in (0, 1)
			g1 = [fermion_greater(τ, β=β, μ=-ϵ_d, U=U) for τ in 0:δt:t]
			g2 = [fermion_lesser(τ, β=β, μ=-ϵ_d, U=U) for τ in 0:δt:t]
			H, a, adag = fermion_operators(U, -ϵ_d)
			g1′, g2′ = gf_real(H, a, adag, β, t, N)
			@test norm(g1 - g1′) / norm(g1) < tol
			@test norm(g2 - g2′) / norm(g2) < tol
		end
	end
end
