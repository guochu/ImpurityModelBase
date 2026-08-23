using QuadGK, Interpolations, LinearAlgebra

# spectrum functions: wrappers and predefined spectrum densities
include("spectrumfuncs/spectrumfuncs.jl")

# baths: basic definitions of particle types, spectrum functions and temperature
include("baths/baths.jl")

# exact diagonalizations: exact numerical solutions for noninteracting cases
include("exactdiagonalizations/exactdiagonalizations.jl")


# analytical solutions for specific models
include("analyticsolutions/analyticsolutions.jl")


# utilities: fourier transforms between time and frequency domains
# (in subfolder fouriertransforms) and other helper functions
include("utilities/utilities.jl")