module ImpurityModelBase

# definition of bath, spectrum density
export SpectrumFunction, lowerbound, upperbound, semicircular, Leggett

export AbstractBath, AbstractFermionicBath, AbstractBosonicBath
export FermionicBath, FermionicVacuum, fermionicbath, thermaloccupation, fermidirac
export BosonicBath, BosonicVacuum, bosonicbath, boseeinstein


# Fourier transformations
export Gt_to_Gw, Gt_to_Δw, Aw_to_Gτ, Gw_to_Aw, frequencies
export Gτ_to_Giw, Giw_to_Gτ, Δτ_to_Δiw, ifrequency, ifrequencies

# analytical solutions 
export freefermion_greater, freefermion_lesser, freefermion_Gt, freefermion_Gτ, freefermion_occupation
export toulouse_Gw, toulouse_Gt, toulouse_Δw, toulouse_Jw
export toulouse_Giw, toulouse_Gτ, toulouse_Δiw, toulouse_Δτ

export freeboson_greater, freeboson_lesser, freeboson_Gt, freeboson_Gτ, freeboson_occupation
export independentbosons_Gτ, independentbosons_greater, independentbosons_lesser

# utilities
export AbstractPredictionScheme, LinearPrediction, linear_predict
export bethe_Gw_to_Δw, bethe_Giw_to_Δiw

using QuadGK

include("spectrumfunc.jl")
include("bath/bath.jl")


# fourier transformations
include("fourier/fourier.jl")

# collections of some analytical solutions
include("freefermion.jl")
include("toulouse/toulouse.jl")

include("freeboson.jl")
include("holstein/holstein.jl")


# utilities for real spectrum functions
include("utilities/utilities.jl")


end