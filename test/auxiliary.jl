println("------------------------------------")
println("|            auxiliary             |")
println("------------------------------------")

f(D, ϵ) = sqrt(1-(ϵ/D)^2) / π
spectrum_func(D) = SpectrumFunction(ϵ->f(D, ϵ), lb=-D, ub=D)
spectrum_func() = spectrum_func(10)

@testset "Fermionic Bath" begin
	β = Inf
	bath = fermionicbath(spectrum_func(), β=β, μ=1)
	@test isa(bath, FermionicVacuum)
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
	@test thermaloccupation(bath, 0.3) == 1 / (exp(bath.β * (0.3 - bath.μ)) + 1)

end