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


	giw = toulouse_Giw(fermionicbath(spec, β=β); ϵ_d=ϵ_d, n=n) 
	# println("giw[1]=", giw[1], " giw[end]=", giw[end])

	gτ = Giw_to_Gτ(giw; β=β, Nτ=Nτ)
	# println("gτ[1]=", gτ[1], " gτ[end]=", gτ[end])

	giw′ = Gτ_to_Giw(gτ, β=β, n=n)

	@test norm(giw - giw′) / norm(giw) < 1.0e-2
end

@testset "Fourier utilities" begin
	β = 10.0
	@test ifrequency(β, 1) ≈ π / β
	@test ifrequency(β, 2) ≈ 3π / β
	fs = ifrequencies(β, 2)
	@test length(fs) == 6
	@test fs ≈ [(2n - 1) * π / β for n in -2:3]

	# Gw_to_Aw: A(ω) = -Im G(ω)/π, negative parts truncated
	gw = [1.0 - 2.0im, 0.5 + 0.0im]
	@test Gw_to_Aw(gw) ≈ [2 / π, 0.0]
	@test Gw_to_Aw([0.5 + 0.02im]; verbosity=0) == [0.0]

	# Δw_to_Jw is the same transformation as Gw_to_Aw (regression)
	@test Δw_to_Jw(gw) ≈ Gw_to_Aw(gw)

	# Δτ_to_Δiw delegates to the fermionic imaginary-time transform
	gτ = [0.6, 0.55, 0.5, 0.45, 0.4]
	@test Δτ_to_Δiw(gτ; β=β, n=50) ≈ Gτ_to_Giw(gτ; β=β, n=50)

	# Aw_to_Gτ for a narrow gaussian peak at ω0
	ω0 = 0.5
	σ = 0.05
	dw = 0.005
	wmin = -2.0
	ws = collect(wmin:dw:3.0)
	Aw = exp.(-(ws .- ω0) .^ 2 / (2σ^2)) / (σ * sqrt(2π))
	Aw ./= sum(Aw) * dw
	Nτ = 40
	Gτ = Aw_to_Gτ(Aw; β=β, Nτ=Nτ, wmin=wmin, δw=dw)
	@test length(Gτ) == Nτ + 1
	τs = range(0, β, length=Nτ + 1)
	for i in 1:3
		@test Gτ[i] ≈ -exp(-ω0 * τs[i]) / (1 + exp(-β * ω0)) atol=5.0e-3
	end
end


