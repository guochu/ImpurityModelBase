println("------------------------------------")
println("|       Fourier Transform          |")
println("------------------------------------")



@testset "Real time and real frequency" begin

	spec = spectrum(ϵ->(2/π)*sqrt(1-ϵ^2), -1, 1)

	wmin = -20.
	wmax = 20.
	dw = 1.0e-2
	ws = collect(wmin:dw:wmax)

	δt = 1.0e-2
	ts = 0:δt:100

	gw = [toulouse_Δw(spec, w) for w in ws]
	# println("gw[1]=", gw[1], " gw[end]=", gw[end])

	gt = Gw_to_Gt(gw, ts; wmin=wmin, δw=dw)
	# println("gt[1]=", gt[1], " gt[end]=", gt[end])

	gw′ = Gt_to_Gw(gt, ws, δt=δt)

	@test norm(gw - gw′) / norm(gw) < 2.0e-2
end


@testset "Imaginary time and Imaginary frequency" begin

	spec = spectrum(ϵ->(2/π)*sqrt(1-ϵ^2), -1, 1)

	n = 1000

	β = 10
	Nτ = 100000
	ϵ_d = 0.


	giw = toulouse_Giw(spec, β=β, ϵ_d=ϵ_d, n=n) 
	# println("giw[1]=", giw[1], " giw[end]=", giw[end])

	gτ = Giw_to_Gτ(giw; β=β, Nτ=Nτ)
	# println("gτ[1]=", gτ[1], " gτ[end]=", gτ[end])

	giw′ = Gτ_to_Giw(gτ, β=β, n=n)

	@test norm(giw - giw′) / norm(giw) < 1.0e-2
end


