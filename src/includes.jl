using QuadGK

include("spectrumfuncs/spectrumfuncs.jl")
include("bath/bath.jl")


# fourier transformations
include("fourier/fourier.jl")

# collections of some analytical solutions
include("freefermion.jl")
include("toulouse/toulouse.jl")

include("freeboson.jl")
include("independentbosons/independentbosons.jl")


# utilities for real spectrum functions
include("utilities/utilities.jl")