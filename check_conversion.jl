
# using ImpurityModelBase
include("src/includes.jl")
using LinearAlgebra
using Interpolations
using DelimitedFiles

f(D, ϵ) = sqrt(1-(ϵ/D)^2) / π
spectrum_func(D=1) = SpectrumFunction(ϵ->f(D, ϵ), lb=0, ub=D)


# function check_linearpredict()
# 	stepsize = 0.1
# 	ts = 0:stepsize:50.
# 	f = spectrum_func()
# 	Gt = [toulouse_Gt(f, t, ϵ_d=0.1) for t in ts]

# 	Gt_predicted = linear_predict(Gt, stepsize, maxiter=200)
# 	Gt′ = [toulouse_Gt(f, (i-1)*stepsize, ϵ_d=0.1) for i in 1:length(Gt_predicted)]

# 	return Gt′, Gt_predicted
# end


function main()
	β = 10.
	ϵ_d = 0.
	spec = spectrum_func()
	gτ = toulouse_Gτ(spec, β=β, ϵ_d=ϵ_d, N=100)

	giw = Gτ_to_Giw(gτ, β=β)

	giw′ = toulouse_Giw(spec, β=β, ϵ_d=ϵ_d) 

	# gτ′ = Giw_to_Gτ(giw, β=β, N=100)
	# return norm(gτ - gτ′) / norm(gτ)

	# return norm(giw - giw′) / norm(giw)
    return giw, giw′
end

β = 10.
ϵ_d = 0.
spec = spectrum_func()
gτ = toulouse_Gτ(spec, β=β, ϵ_d=ϵ_d, N=100)

giw1 = Gτ_to_Giw(gτ, β=β)

giw2 = toulouse_Giw(spec, β=β, ϵ_d=ϵ_d) 


interp=linear_interpolation(0:0.1:10.0, gτ)
gτ1 = interp.(0:0.01:10.0)
giw3 = Gτ_to_Giw(gτ1, β=10.0)

gτ2 = Giw_to_Gτ(giw2, β=10, N=100)

writedlm("a.dat", [1:2002 real.(giw1) imag.(giw1)])
writedlm("b.dat", [1:2002 real.(giw2) imag.(giw2)])
writedlm("c.dat", [1:2002 real.(giw3) imag.(giw3)])
writedlm("d.dat", [1:101 gτ gτ2])
