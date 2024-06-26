module ImpurityModelBase

# definition of bath, spectrum density
export AbstractFermionicBath, FermionicBath, FermionicVacuum, fermionicbath, thermaloccupation
export SpectrumFunction, lowerbound, upperbound, semicircular

# Fourier transformations
export Gt_to_Gw, Gt_to_Δw, Aw_to_Gτ, Gw_to_Aw, frequencies
export Gτ_to_Giw, Giw_to_Gτ, Δτ_to_Δiw, ifrequency, ifrequencies

# analytical solutions 
export free_greater, free_lesser, free_Gt, free_Gτ
export toulouse_Gw, toulouse_Gt, toulouse_Δw, toulouse_Jw
export toulouse_Giw, toulouse_Gτ, toulouse_Δiw, toulouse_Δτ

# utilities
export AbstractPredictionScheme, LinearPrediction, linear_predict
export bethe_Gw_to_Δw, bethe_Giw_to_Δiw

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


end