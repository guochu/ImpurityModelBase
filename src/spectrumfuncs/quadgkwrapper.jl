# simple wrapper for quadgk function
# add support for delta function

const QuadGKTolerance = 1.0e-6

"""
	quadgkwrapper(f; kwargs...)

A wrapper for quadgk function, where f is a single-variate function of real numbers
It is assumed that the lower and upper bound information is contained in f, so f
is in general a customised function
"""
quadgkwrapper(f; kwargs...) = error("quadgkwrapper not implemented for function type $(typeof(f))")
function _quadgk(f, lb::Number, ub::Number; kwargs...)
	r, err = quadgk(f, lb, ub; kwargs...)
	(err <= QuadGKTolerance) || println("quadgk error $(err) larger than tol $(QuadGKTolerance)")
	return r
end