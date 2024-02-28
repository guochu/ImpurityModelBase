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

Gt_to_Gw(gt::Vector{<:Number}, ws::Vector{<:Real}; δt::Real) = gf_retarded_ω(gt, ws, δt=δt)

