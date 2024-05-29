using QuadGK

include("spectrumfunc.jl")
include("bath.jl")


# fourier transformations
include("fourier/fourier.jl")

# collections of some analytical solutions
include("freefermion.jl")
include("toulouse/toulouse.jl")


# utilities for real spectrum functions
include("utilities/utilities.jl")