# # simple fourier transform
# abstract type FourierTransformScheme end

# struct FourierTransform <: FourierTransformScheme
# 	Gt::Vector{ComplexF64}
# 	δt::Float64
# 	δ::Float64
# end
# function FourierTransform(Gt::Vector{<:Number}; δt::Real, δ::Real=0.) 
# 	(δ >= 0.) || throw(ArgumentError("δ should not be negative"))
# 	return FourierTransform(convert(Vector{ComplexF64}, Gt), convert(Float64, δt), convert(Float64, δ))
# end 


# # function (x::FourierTransform)(ω::Real)
# # 	# exp_t = exp(im * ω * x.δt)
# # 	r = zero(ComplexF64)
# # 	Gt = x.Gt

# # 	# coef = 1.

# # 	for i in 1:length(Gt)
# # 		# r += coef * Gt[i] * x.δt
# # 		# coef *= exp_t
# # 		r += Gt[i] * x.δt * exp(im * ω * i * x.δt - x.δ * i * x.δt)
# # 	end

# # 	return r
# # end

# function (x::FourierTransform)(ω::Real)
# 	# exp_t = exp(im * ω * x.δt)
# 	r = zero(ComplexF64)
# 	Gt = x.Gt

# 	# coef = 1.
# 	ω′ = im * ω - x.δ 

# 	for i in 1:length(Gt)
# 		# r += coef * Gt[i] * x.δt
# 		# coef *= exp_t
# 		tj = i * x.δt 
# 		r += Gt[i] * exp(ω′ * tj)
# 	end

# 	return r * x.δt
# end

# function gf_retarded_ω(gt::Vector{<:Number}, ws::Vector{<:Real}; δt::Real)
# 	x = FourierTransform(gt, δt=δt)
# 	return [x(w) for w in ws]
# end 

# function Gt_to_Gw(gt::Vector{<:Number}, stepsize::Real, ws::Vector{<:Real})
# 	x = FourierTransform(gt, δt=stepsize)
# 	return [x(w) for w in ws]
# end 

"""
	Gw_to_Aw(Gw::Vector{<:Number}; verbosity)

Obtaining spectrum density A(ω) from frequency-domain Green's function G(ω),
using the relation A(ω) = -imag(G(ω)) / π
"""
function Gw_to_Aw(Gw::Vector{<:Number}; verbosity::Int=1)
	Aw = zeros(real(eltype(Gw)), length(Gw))
	for i in 1:length(Gw)
		Gwi = -imag(Gw[i])/π
		if (verbosity > 0) && (Gwi < 0) && (abs(Gwi) > 1.0e-3)
			println("negative Aw[$(i)] = ", Gwi)
		end
		Aw[i] = (Gwi >= 0) ? Gwi : zero(Gwi)
	end
	return Aw
end

"""
	Δw_to_Jw(Δw::Vector{<:Number}; verbosity)

Obtaining the bath spectrum density J(ω) from the frequency-domain hybridization function Δ(ω),
using the relation J(ω) = -imag(Δ(ω)) / π (the relation is the same as Gw_to_Aw)
"""
Δw_to_Jw(Δw::Vector{<:Number}; verbosity::Int=1) = Gw_to_Aw(Δω, verbosity=verbosity)

"""
	Gt_to_Gw(gt::Vector{<:Number}, ws; δt, δ::Real=1.0e-8)

Convert the real-time Green's function G(t) to the real-frequency Green's function G(ω)
using Fourier transformation

G(t) is assumed to be equally spaced with stepsize δt, and should have converged to 0
"""
Gt_to_Gw(gt::Vector{<:Number}, ws::Union{Vector{<:Real}, AbstractRange}; δt::Real, δ::Real=1.0e-8) = [_fourier(gt, w, δt=δt, δ=δ) for w in ws]

"""
	Gw_to_Gt(gw::Vector{<:Number}, ts; wmin, δw, δ::Real=1.0e-8)

Convert the real-frequency Green's function G(ω) tp the real-time Green's function G(t)
using inverse Fourier frequency

G(ω) is assumed to be equally spaced with stepsize δw,
wmin is the smallest frequency

In the implemented the assympototic behavior 1/ω is used for G(ω) with ω < wmin
or w > wmax (wmax is calculated based on the size of G(ω), wmin and δw)
"""
Gw_to_Gt(gw::Vector{<:Number}, ts::Union{Vector{<:Real}, AbstractRange}; wmin::Real, δw::Real=1.0e-6, δ::Real=1.0e-8) = [_inverse_fourier(gw, t, wmin=wmin, δw=δw, δ=δ) for t in ts]



"""
	Aw_to_Gτ(Aw::Vector{<:Real}; β::Real, wmin::Real, δw::Real=1.0e-4, Nτ::Int)

Convert spectrum density A(ω) to imaginary-time Green's function G(τ), using 
inverse analytical continuation

wmin is the smallest real frequency for A(ω), δw is the real-frequency 
stepsize of A(ω)

β is the inverse temperature, Nτ is the number of discrete points 
on the imaginary-time axis
"""
function Aw_to_Gτ(Aw::Vector{<:Real}; β::Real, wmin::Real, δw::Real=1.0e-4, Nτ::Int)
	# check the summation of Aw
	nrm = sum(Aw) * δw
	if abs(nrm-1) > 1.0e-2
		println("sum(A(ω)) = ", nrm)
	end
	τs = range(0, β, length=Nτ+1)
	Gτ = zeros(eltype(Aw), length(τs))
	for (j, Awi) in enumerate(Aw)
		w = wmin + (j-1) * δw
		tmp = -Awi * δw / (1+exp(-β*w))
		for (i, τ) in enumerate(τs)
			# Gτ[i] += -exp(-w*τ) * Awi * dw / (1+exp(-β*w))
			Gτ[i] += exp(-w*τ) * tmp
		end
	end
	return Gτ
end

# frequencies(;lb::Real, ub::Real, dw::Real=1.0e-4) = lb:dw:ub


function _fourier(gt::Vector{<:Number}, ω::Real; δt::Real, δ::Real)
	# exp_t = exp(im * ω * x.δt)
	r = zero(ComplexF64)
	ω′ = im * ω - δ 

	for i in 1:length(gt)-1
		tj = (i-1) * δt
		r += gt[i] * exp(ω′ * tj)
	end

	return r * δt
end



function _inverse_fourier(gw::Vector{<:Number}, t::Real; wmin::Real, δw::Real, δ::Real)
	r = zero(eltype(gw))
	for i in 1:length(gw)-1
		wj = wmin + δw*(i-1)
		# r += (gw[i] - 1.0/(wj+im*δ)) * exp(-im*wj*t) 
		# r += (gw[i] - 0.5/(wj+im) - 0.5/(wj-im)) * exp(-im*wj*t) 
		r += (gw[i] - 1.0/(wj+im)) * exp(-im*wj*t) 
	end
	return (r*δw) / (2π) - im*exp(-t)
end
