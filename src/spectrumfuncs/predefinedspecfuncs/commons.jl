# function semicircular(t::Real)
# 	t = convert(Float64, t)
# 	D = 2*t
# 	return spectrum(ϵ -> sqrt(1-(ϵ/D)^2) * (D/π), lb = -D, ub = D)
# end
# semicircular(; t::Real=1) = semicircular(t)

"""
	semicircular(t::Real)

Semi-circular bath spectrum density, often 
used for fermions
"""
semicircular(t) = spectrum(ϵ->(2/(pi*t^2)) * sqrt(t^2 - ϵ^2), lb=-t, ub=t)
semicircular(; t::Real=1) = semicircular(t)

"""
	Leggett(; α::Real, d::Real, ωc::Real)

phenomelogical bath spectrum density,
J(ω) = (α ωc/2)(ωᵈ/ωcᵈ)e^(-ω/ωc),	 
often used for bosons
"""
function Leggett(; α::Real=1, d::Real=1, ωc::Real=1)
	d = convert(Float64, d)
	return spectrum(ϵ -> (α/2)*(ϵ^d/ωc^(d-1))*exp(-ϵ/ωc), lb = 0, ub = 100*ωc)
end