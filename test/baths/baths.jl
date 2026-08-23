println("------------------------------------")
println("|              baths              |")
println("------------------------------------")

f(D, ϵ) = sqrt(1-(ϵ/D)^2) / π
spectrum_func(D) = spectrum(ϵ->f(D, ϵ), lb=0, ub=D)
spectrum_func() = spectrum_func(10)

@testset "Fermionic Bath" begin
	bath = fermionicvacuum(spectrum_func(), μ=1)
	@test isa(bath, FermionicVacuum)
	@test bath.β == Inf
	@test bath.T == 0
	@test bath.μ == 1
	@test thermaloccupation(bath, 0.99) == 1
	@test thermaloccupation(bath, 1.01) == 0

	β = Inf
	bath = fermionicbath(spectrum_func(), β=β, μ=1)
	@test isa(bath, FermionicBath)
	@test bath.β == β
	@test bath.T == 0
	@test bath.μ == 1
	@test thermaloccupation(bath, 0.99) == 1
	@test thermaloccupation(bath, 1.01) == 0
	β = 0.25
	bath = fermionicbath(spectrum_func(), β=β)
	@test isa(bath, FermionicBath)
	@test bath.β == β
	@test bath.T == 4
	@test bath.μ == 0
	@test thermaloccupation(bath, 0.3) ≈ 1 / (exp(bath.β * (0.3 - bath.μ)) + 1)

end


@testset "Bosonic Bath" begin
	bath = bosonicvacuum(spectrum_func(), μ=-1)
	@test isa(bath, BosonicVacuum)
	@test bath.β == Inf
	@test bath.T == 0
	@test bath.μ == -1
	@test thermaloccupation(bath, -0.9) == 0
	β = Inf
	bath = bosonicbath(spectrum_func(), β=β, μ=-1)
	@test isa(bath, BosonicBath)
	@test bath.β == β
	@test bath.T == 0
	@test bath.μ == -1
	@test thermaloccupation(bath, -0.9) == 0
	β = 0.25
	bath = bosonicbath(spectrum_func(), β=β)
	@test isa(bath, BosonicBath)
	@test bath.β == β
	@test bath.T == 4
	@test bath.μ == 0
	@test thermaloccupation(bath, 0.3) ≈ 1 / (exp(bath.β * (0.3 - bath.μ)) - 1)

end

@testset "Particle types and baths" begin
	@test Boson <: AbstractParticle
	@test Fermion <: AbstractParticle

	spec = spectrum(ϵ -> 1 - ϵ^2, -1, 1)

	# generic bath / vacuum constructors
	b1 = bath(Boson, spec; β=2.0, μ=-0.5)
	@test b1 isa BosonicBath
	@test b1 isa AbstractBosonicNormalBath
	@test b1 isa AbstractNormalBath{Boson}
	@test b1 isa AbstractBath{Boson}
	@test particletype(b1) == Boson
	@test b1.β == 2.0 && b1.μ == -0.5 && b1.T == 0.5

	b2 = bath(Fermion, spec; β=2.0, μ=0.3)
	@test b2 isa FermionicBath
	@test b2 isa AbstractFermionicNormalBath
	@test particletype(b2) == Fermion

	v1 = vacuum(Boson, spec; μ=0.1)
	@test v1 isa BosonicVacuum
	@test v1 isa AbstractBath{Boson}
	@test v1.β == Inf && v1.T == 0

	v2 = vacuum(Fermion, spec; μ=0.1)
	@test v2 isa FermionicVacuum

	# abstract bath type relations
	bc = BCSBath(spec; β=2.0, Δ=0.3)
	@test bc isa AbstractBCSBath
	@test bc isa AbstractBath{Fermion}
	@test AbstractBECBath <: AbstractBath{Boson}

	# thermal distributions
	@test fermidirac(1.0, 0.5) ≈ 1 / (1 + exp(0.5))
	@test fermidirac(1.0, 0.5, 0.2) ≈ 1 / (1 + exp(-0.3))   # fermidirac(β, μ, ϵ) = f(ϵ-μ)
	@test thermaloccupation(Fermion, 1.0, 0.5, 0.2) ≈ fermidirac(1.0, 0.5, 0.2)
	@test boseeinstein(1.0, 0.0, 1.0) ≈ 1 / (exp(1.0) - 1)
	@test thermaloccupation(Boson, 1.0, 0.0, 1.0) ≈ boseeinstein(1.0, 0.0, 1.0)
	@test thermaloccupation(b2, 0.2) ≈ fermidirac(2.0, 0.3, 0.2)
end

@testset "Discrete baths" begin
	ws = [0.5, 1.0, 1.5]
	fs = [0.3, 0.5, 0.2]
	β = 5.0

	bf = discretefermionicbath(ws, fs; β=β, μ=0.2)
	@test bf isa DiscreteBath{Fermion}
	@test bf isa FermionicBath
	@test frequencies(bf) == ws
	@test spectrumvalues(bf) == fs
	@test spectrumcouplings(bf) ≈ sqrt.(fs)
	@test num_sites(bf) == 3

	bb = discretebosonicbath(ws, fs; β=β)
	@test bb isa DiscreteBath{Boson}
	@test particletype(bb) == Boson

	vf = discretefermionicvacuum(ws, fs; μ=0.2)
	@test vf isa DiscreteVacuum{Fermion}
	@test vf.β == Inf
	vb = discretebosonicvacuum(ws, fs)
	@test vb isa DiscreteVacuum{Boson}

	# generic discrete constructors
	bf2 = discretebath(Fermion, ws, fs; β=β, μ=0.2)
	@test bf2 == bf
	@test discretebath(Boson, ws, fs; β=β) isa DiscreteBath{Boson}
	vf2 = discretevacuum(Fermion, ws, fs; μ=0.2)
	@test vf2 isa DiscreteVacuum{Fermion}
	@test discretevacuum(Boson, ws, fs) isa DiscreteVacuum{Boson}
end

@testset "BCS baths" begin
	ws = [0.5, 1.0, 1.5]
	fs = [0.3, 0.5, 0.2]
	spec = spectrum(ϵ -> 1 - ϵ^2, -1, 1)

	bc = bcsbath(spec; β=3.0, μ=0.1, Δ=0.4)
	@test bc isa BCSBath
	@test bc isa AbstractBCSBath
	@test bc.Δ == 0.4

	bv = bcsvacuum(spec; μ=0.1, Δ=0.4)
	@test bv isa BCSVacuum
	@test bv.β == Inf && bv.T == 0

	# construct from a normal fermionic bath
	fb = fermionicbath(spec; β=3.0, μ=0.1)
	bcf = bcsbath(fb; Δ=0.4)
	@test bcf.Δ == 0.4 && bcf.β == fb.β
	fv = fermionicvacuum(spec; μ=0.1)
	@test bcsvacuum(fv; Δ=0.4).μ == fv.μ

	# discrete BCS baths
	db = discretebcsbath(ws, fs; β=3.0, μ=0.1, Δ=0.4)
	@test db isa DiscreteBCSBath
	@test frequencies(db) == ws
	@test spectrumcouplings(db) ≈ sqrt.(fs)
	@test num_sites(db) == 2 * length(ws)

	dv = discretebcsvacuum(ws, fs; μ=0.1, Δ=0.4)
	@test dv isa DiscreteBCSVacuum
	@test dv.β == Inf

	# discretize a continuous BCS bath (regression: discretevacuum(BCSVacuum))
	dbc = discretebath(bc; δw=0.1)
	@test dbc isa DiscreteBCSBath
	@test num_sites(dbc) > 0
	dvc = discretevacuum(bv; δw=0.1)
	@test dvc isa DiscreteBCSVacuum
	@test num_sites(dvc) > 0
end