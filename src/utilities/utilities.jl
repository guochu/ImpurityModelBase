module Utilities

using QuadGK
using ImpurityModelBase

# utilities for real spectrum functions
export AbstractPredictionScheme, LinearPrediction
export FourierTransformScheme, FourierTransform, gf_retarded_ω
export Gt_to_Gw

# analytical solutions in the noninteracting case
export toulouse_Gw, toulouse_Gt, toulouse_Giw, toulouse_Gτ

include("linearprediction.jl")
include("Gt2Gw.jl")
include("toulouse.jl")


end