using QuadGK, Interpolations, LinearAlgebra

# spectrum functions
include("spectrumfuncs/spectrumfuncs.jl")

# baths
include("bath.jl")
include("discretebath.jl")

# bcs baths
include("bcsbath.jl")

# exact diagonalizations
include("exactdiagonalizations/exactdiagonalizations.jl")

# fourier transformations
include("fouriertransforms/fouriertransforms.jl")


# some analytical solutions
include("analyticsolutions/analyticsolutions.jl")


# utilities for real spectrum functions
include("utilities/utilities.jl")
