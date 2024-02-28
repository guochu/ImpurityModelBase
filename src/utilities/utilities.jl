module Utilities

# utilities for real spectrum functions
export AbstractPredictionScheme, LinearPrediction
# export FourierTransformScheme, FourierTransform, gf_retarded_ω
export Gt_to_Gw

include("linearprediction.jl")
include("Gt2Gw.jl")


end