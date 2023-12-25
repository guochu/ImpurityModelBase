"""
	SpectrumFunction{F}

Wrapper for spectrum function, including a lower bound and upper bound
"""
struct SpectrumFunction{F}
	f::F
	lb::Float64
	ub::Float64
end
function SpectrumFunction(f; lb::Real=-Inf, ub::Real=Inf) 
	(lb < ub) || throw(ArgumentError("lb must be less than ub"))
	return SpectrumFunction(f, convert(Float64, lb), convert(Float64, ub))
end	
lowerbound(x::SpectrumFunction) = x.lb
upperbound(x::SpectrumFunction) = x.ub
# function (x::SpectrumFunction)(ϵ::Real) 
# 	if lowerbound(x) <= ϵ <= upperbound(x)
# 		return x.f(ϵ)
# 	else
# 		return 0.
# 	end
# end
function (x::SpectrumFunction)(ϵ::Real) 
	@assert lowerbound(x) <= ϵ <= upperbound(x)
	return x.f(ϵ)
end
