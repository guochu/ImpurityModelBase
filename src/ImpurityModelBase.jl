module ImpurityModelBase

# definition of bath, spectrum density
export AbstractFermionicBath, FermionicBath, FermionicVacuum, fermionicbath, thermaloccupation
export SpectrumFunction, lowerbound, upperbound

# analytic
export free_greater, free_lesser, free_Gt, free_Gτ
export toulouse_Gw, toulouse_Gt, toulouse_Giw, toulouse_Gτ
export toulouse_Δiw, toulouse_Δτ

# utilities
export AbstractPredictionScheme, LinearPrediction, linear_predict
export FourierTransformScheme, FourierTransform, Gt_to_Gw
export Gt_to_Aw, Aw_to_Gτ, Gτ_to_Giw, Giw_to_Gτ


using QuadGK

include("spectrumfunc.jl")
include("bath.jl")

# collections of some analytical solutions
include("analytic/analytic.jl")


# utilities for real spectrum functions
include("utilities/utilities.jl")


end