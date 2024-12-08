# simple fourier transform
abstract type FourierTransformScheme end

struct FourierTransform <: FourierTransformScheme
	Gt::Vector{ComplexF64}
	δt::Float64
	δ::Float64
end
function FourierTransform(Gt::Vector{<:Number}; δt::Real, δ::Real=0.) 
	(δ >= 0.) || throw(ArgumentError("δ should not be negative"))
	return FourierTransform(convert(Vector{ComplexF64}, Gt), convert(Float64, δt), convert(Float64, δ))
end 


# function (x::FourierTransform)(ω::Real)
# 	# exp_t = exp(im * ω * x.δt)
# 	r = zero(ComplexF64)
# 	Gt = x.Gt

# 	# coef = 1.

# 	for i in 1:length(Gt)
# 		# r += coef * Gt[i] * x.δt
# 		# coef *= exp_t
# 		r += Gt[i] * x.δt * exp(im * ω * i * x.δt - x.δ * i * x.δt)
# 	end

# 	return r
# end

function (x::FourierTransform)(ω::Real)
	# exp_t = exp(im * ω * x.δt)
	r = zero(ComplexF64)
	Gt = x.Gt

	# coef = 1.
	ω′ = im * ω - x.δ 

	for i in 1:length(Gt)
		# r += coef * Gt[i] * x.δt
		# coef *= exp_t
		tj = i * x.δt 
		r += Gt[i] * exp(ω′ * tj)
	end

	return r * x.δt
end

function gf_retarded_ω(gt::Vector{<:Number}, ws::Vector{<:Real}; δt::Real)
	x = FourierTransform(gt, δt=δt)
	return [x(w) for w in ws]
end 

# function Gt_to_Gw(gt::Vector{<:Number}, stepsize::Real, ws::Vector{<:Real})
# 	x = FourierTransform(gt, δt=stepsize)
# 	return [x(w) for w in ws]
# end 

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
Δw_to_Jw(Δw::Vector{<:Number}; verbosity::Int=1) = Gw_to_Aw(Δω, verbosity=verbosity)

function Gt_to_Gw(gt::Vector{<:Number}, stepsize::Real; lb::Real, ub::Real, dw::Real=1.0e-4, δ::Real=1.0e-8)
	x = FourierTransform(gt, δt=stepsize, δ=δ)
	return [x(w) for w in lb:dw:ub]
end 
function Gw_to_Gt(Gw::Vector{<:Number}, dw::Real; lb::Real, ub::Real, δt::Real=1.0e-4, δ::Real=1.0e-8)
	iseven(length(Gw)) && throw(ArgumentError("odd number of frequencies expected"))
	n = div(length(Gw), 2)
	ts = lb:δt:ub
	Gt = zeros(eltype(Gw), length(ts))
	for (i, tj) in enumerate(ts)
		r = zero(eltype(Gw))
		for nj in -n:n
			wj = nj * dw
			if nj != 0
				r += (Gw[nj+n+1] - 1.0/(wj+im*δ)) * exp(-im*wj*tj)
			end
			# r += (Gw[nj+n+1] - 1.0/(wj+im*δ)) * exp(-im*wj*tj)
		end
		Gt[i] = r * dw / (2*pi) - im
	end
	return Gt
end

function Gt_to_Aw(gt::Vector{<:Number}, stepsize::Real; lb::Real, ub::Real, dw::Real=1.0e-4, normalize::Bool=true, δ::Real=1.0e-8, verbosity::Int=1)
	if (verbosity > 0) && (abs(gt[end]) > 1.0e-6)
		println("last element of input GF has abs value ", abs(gt[end]))
	end
	Gw = Gt_to_Gw(gt, stepsize, lb=lb, ub=ub, dw=dw, δ=δ)
	Aw = Gw_to_Aw(Gw, verbosity=verbosity)
	nrm = sum(Aw) * dw
	if (verbosity > 0) && (abs(nrm-1) > 1.0e-2)
		println("sum(A(ω)) = ", nrm)
	end	
	if normalize
		Aw ./= nrm
	end
	return Aw
end

function Aw_to_Gτ(Aw::Vector{<:Real}; β::Real, lb::Real, ub::Real, dw::Real=1.0e-4, δτ::Real=0.1)
	# check the summation of Aw
	nrm = sum(Aw) * dw
	if abs(nrm-1) > 1.0e-2
		println("sum(A(ω)) = ", nrm)
	end
	τs = 0:δτ:β
	Gτ = zeros(eltype(Aw), length(τs))
	for (w, Awi) in zip(lb:dw:(ub-dw), Aw)
		tmp = -Awi * dw / (1+exp(-β*w))
		for (i, τ) in enumerate(τs)
			# Gτ[i] += -exp(-w*τ) * Awi * dw / (1+exp(-β*w))
			Gτ[i] += exp(-w*τ) * tmp
		end
	end
	return Gτ
end

frequencies(;lb::Real, ub::Real, dw::Real=1.0e-4) = lb:dw:ub

# function Aw_to_Git_single(Aw::Vector{<:Real}, τ::Real; β::Real, lb::Real, ub::Real, dw::Real=1.0e-4)
# 	gτ = zero(eltype(Aw))
# 	for w, Awi in zip(lb:dw:(ub-dw), Aw)
# 		gτ += -exp(w*τ) * Awi / (1+exp(-β*w))
# 	end
# 	return gτ
# end

