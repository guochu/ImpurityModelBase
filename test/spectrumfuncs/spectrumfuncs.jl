# Tests for the exported spectrum functions (src/spectrumfuncs: bounded spectrum
# wrappers and discrete spectra) not exercised elsewhere in the test suite.

println("------------------------------------")
println("|      spectrum functions          |")
println("------------------------------------")

@testset "Bounded spectrum functions" begin
	# bounded / BoundedFunction / spectrum
	f = bounded(ϵ -> 1 - ϵ^2, -1, 1)
	@test f isa BoundedFunction
	@test f isa AbstractBoundedFunction
	@test AbstractBoundedFunction <: Function
	@test lowerbound(f) == -1
	@test upperbound(f) == 1
	@test f(0.5) ≈ 1 - 0.5^2
	@test f(1.5) == 0.0          # outside the interval -> 0
	@test f(-1.5) == 0.0

	g = spectrum(ϵ -> 0.5, 0, 2)
	@test lowerbound(g) == 0 && upperbound(g) == 2

	# quadgkwrapper: semi-circular spectrum integrates to 1
	semi = semicircular(1)
	@test quadgkwrapper(semi) ≈ 1.0 atol=1.0e-8

	# negation and adjoint of a bounded function
	neg = -semi
	@test neg(0.5) ≈ -semi(0.5)
	adj = semi'
	@test adj(0.5) ≈ conj(semi(0.5))

	# spectrumshift moves the frequency axis (use a spectrum whose defining
	# function is safe outside its interval, since the call evaluates f first)
	sp = spectrum(ϵ -> 1 - ϵ^2, -1, 1)
	sh = spectrumshift(sp, 0.3)
	@test lowerbound(sh) ≈ -1.3
	@test upperbound(sh) ≈ 0.7
	@test sh(-1.0) ≈ 1 - 0.7^2      # f(ϵ+μ) = f(-0.7)
	@test sh(0.8) == 0.0            # outside the shifted interval
end

@testset "Discrete spectrum" begin
	ws = [0.1, 0.5, 1.0, 2.0]
	fs = [0.4, 0.8, 0.6, 0.2]
	d = DiscreteSpectrum(ws, fs)
	@test frequencies(d) == ws
	@test spectrumvalues(d) == fs
	@test spectrumcouplings(d) ≈ sqrt.(fs)
	@test lowerbound(d) == -Inf
	@test upperbound(d) == Inf

	ds = spectrumshift(d, 0.2)
	@test frequencies(ds) ≈ ws .+ 0.2
	@test spectrumvalues(ds) ≈ fs

	@test_throws DimensionMismatch DiscreteSpectrum([0.1, 0.2], [1.0])
	@test_throws ArgumentError DiscreteSpectrum([0.1, 0.2], [-1.0, 1.0])
end
