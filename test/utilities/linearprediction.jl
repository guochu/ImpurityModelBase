# Tests for the linear prediction utilities (src/utilities/linearprediction.jl).

println("------------------------------------")
println("|       Linear Prediction          |")
println("------------------------------------")

@testset "Linear prediction" begin
	# linear sequence -> exact linear extrapolation
	obs = collect(0.0:1.0:9.0)
	lp = LinearPrediction(obs, stepsize=1.0, nfit=8, p=2)
	@test lp isa LinearPrediction
	@test lp isa AbstractPredictionScheme
	@test lp[11] ≈ 10.0 atol=1.0e-4
	@test lp[12] ≈ 11.0 atol=1.0e-4

	# exponential decay sequence (error accumulates over extrapolation steps)
	obs2 = exp.(-0.1 .* (0:9))
	lp2 = LinearPrediction(obs2, stepsize=1.0, nfit=10, p=2)
	@test lp2[12] ≈ exp(-0.1 * 11) atol=5.0e-3     # 2 steps beyond the data
	@test lp2[15] ≈ exp(-0.1 * 14) atol=5.0e-2     # 5 steps beyond the data

	# linear_predict end-to-end
	ext = linear_predict(obs, 1.0; δt=1.0, maxiter=20, tol=1.0e-12, verbosity=0)
	@test length(ext) == 30
	@test ext[end] ≈ 29.0 atol=1.0e-4
end
