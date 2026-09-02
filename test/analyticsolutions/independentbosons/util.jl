
# ϵ_d n̂ + α n̂(b̂ + b̂†) + ω₀b̂†b̂ 
function noninteracting_operators(ϵ_d; ω₀=1, α=0.5, d=100)
	p1 = spin_half_matrices()
	n̂, σ₊, σ₋ = p1["n"], p1["+"], p1["-"]
	Is = one(n̂)
	p2 = boson(d=d)
	b̂, b̂′, n̂b = p2["a"], p2["adag"], p2["n"]
	Ib = one(b̂)
	# total Hamiltonian
	Himp = -ϵ_d * kron(n̂, Ib)
	Hbath = ω₀ * kron(Is, n̂b)
	Hhyb = sqrt(α) * kron(n̂, b̂′ + b̂)
	H = Himp + Hhyb + Hbath
	A, B = kron(σ₋, Ib), kron(σ₊, Ib)

	return H, A, B
end

# ϵ_d(n̂↑ + n̂↓) + U n̂↑n̂↓ + α (n̂↑ + n̂↓)(b̂ + b̂†) + ω₀b̂†b̂ 
function interacting_operators(U, ϵ_d=U/2; ω₀=1, α=0.5, d=100)
	p1 = spin_half_matrices()
	n̂, σ₊, σ₋ = p1["n"], p1["+"], p1["-"]
	Is = one(n̂)
	n_ud = kron(n̂, Is) + kron(Is, n̂)
	nn = kron(n̂,n̂)
	p2 = boson(d=d)
	b̂, b̂′, n̂b = p2["a"], p2["adag"], p2["n"]
	Ib = one(b̂)
	# total Hamiltonian
	Himp = kron(-ϵ_d*n_ud + U * nn, Ib)
	Hbath = ω₀ * kron(kron(Is, Is), n̂b)
	Hhyb = sqrt(α) * kron(n_ud, b̂′ + b̂)
	H =  Himp + Hhyb + Hbath

	A, B = kron(kron(σ₋, Is), Ib), kron(kron(σ₊, Is), Ib)

	return H, A, B
end

# function gen_initstate(H, Himp, Hbath, β, init_state::Symbol)
# 	if init_state == :globalthermal
# 		return exp(-β*H)
# 	else
# 		# return exp(-β*Himp) * exp(-β*Hbath)
# 		return exp(-β*(Himp + Hbath)) 
# 	end
# end

# product initial state: impurity density matrix ρ_0 ⊗ thermal equilibrium of a
# single bosonic mode of frequency ω₀, truncated to d levels
#
# `ρ_0` is given in the analytical-solution convention:
#   bands=1: 2×2 matrix in the |0⟩,|1⟩ basis
#   bands=2: 4×4 matrix in the |0⟩,|↑⟩,|↓⟩,|↑↓⟩ basis, where |↑⟩ = first spin
#            occupied, |↓⟩ = second spin occupied.
# The ED basis is the natural kron basis |00⟩,|0↓⟩,|↑0⟩,|↑↓⟩, so the rows/columns
# of the singly-occupied states (indices 2 and 3) must be swapped for bands=2.
function prod_state(ρ_0::AbstractMatrix{<:Number}, ω₀::Real, β::Real, d::Int)
	ρ_bath = Diagonal(exp.(-β .* (0:d-1) .* ω₀))
	if size(ρ_0, 1) == 4
		# code convention [|0⟩,|↑⟩,|↓⟩,|↑↓⟩] -> ED kron basis [|0⟩,|0↓⟩,|↑0⟩,|↑↓⟩]
		ρ_0 = ρ_0[[1, 3, 2, 4], [1, 3, 2, 4]]
	end
	return kron(ρ_0, ρ_bath)
end