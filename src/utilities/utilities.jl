module Utilities

# utilities for real spectrum functions
export ExtrapolationScheme, LinearExtrapolation
export FourierTransformScheme, FourierTransform, gf_retarded_ω

include("extrapolation.jl")
include("fouriertransform.jl")


end