# simple wrapper for quadgk function
# add support for delta function

const QuadGKTolerance = 1.0e-6

function _quadgk(f, lb::Number, ub::Number; kwargs...)
	r, err = quadgk(f, lb, ub; kwargs...)
	(err <= QuadGKTolerance) || println("quadgk error $(err) larger than tol $(QuadGKTolerance)")
	return r
end