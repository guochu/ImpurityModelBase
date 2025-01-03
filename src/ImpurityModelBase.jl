module ImpurityModelBase

# spectrum density wrapper and some predefined spectrum functions
export AbstractBoundedFunction, AbstractSpectrumFunction
export BoundedFunction, SpectrumFunction, bounded, spectrum
export lowerbound, upperbound, semicircular, Leggett, DiscreteSpectrum
export DiracDelta, quadgkwrapper, spectrumshift

# definition of bosonic and fermionic bath, thermal distributions
export AbstractParticle, Boson, Fermion, particletype
export AbstractBath, AbstractBosonicBath, AbstractFermionicBath, Bath, Vacuum
export FermionicBath, FermionicVacuum, fermionicbath, thermaloccupation, fermidirac
export BosonicBath, BosonicVacuum, bosonicbath, boseeinstein, bath, vacuum

# discrete bath
export DiscreteBath, DiscreteVacuum, AbstractDiscreteBosonicBath, AbstractDiscreteFermionicBath
export discretebosonicbath, discretefermionicbath, discretebath, discretevacuum
export frequencies, spectrumvalues, num_sites


# exact diagonalizations
export freefermions_cache, freefermions_timeevo, freefermions_greater_lesser
export Toulouse, toulouse_greater_lesser, toulouse_cmatrix
export BoundaryDriving, cmatrix, separablestate, thermalstate, leftparticlecurrent_cmatrix, rightparticlecurrent_cmatrix

# Fourier transformations
export Gt_to_Gw, Gw_to_Gt, Aw_to_Gτ, Gw_to_Aw, Δw_to_Jw
export Gτ_to_Giw, Giw_to_Gτ, Δτ_to_Δiw, ifrequency, ifrequencies

# analytical solutions 
export freefermion_greater, freefermion_lesser, freefermion_Gt, freefermion_Gτ, freefermion_occupation
export fermion_greater, fermion_lesser, fermion_Gt, fermion_Gτ
export toulouse_Gw, toulouse_Gt, toulouse_Δw
export toulouse_Giw, toulouse_Gτ, toulouse_Δiw, toulouse_Δτ

export freeboson_greater, freeboson_lesser, freeboson_Gt, freeboson_Gτ, freeboson_occupation
export independentbosons_Gτ, independentbosons_greater, independentbosons_lesser

export holstein_scaleless_parameters, holstein_bare_parameters, GreenFunction
export holstein_G0w_to_Gw, holstein_G0w_to_Σw, bethe_holstein_dmft_iteration
export holstein_Gt, holstein_Gw, holstein_Σw

# utilities
export AbstractPredictionScheme, LinearPrediction, linear_predict
# bethe lattice
# export bethe_Gw_to_Δw, bethe_Giw_to_Δiw

using QuadGK, Interpolations, LinearAlgebra

# spectrum functions
include("spectrumfuncs/spectrumfuncs.jl")

# baths
include("bath.jl")
include("discretebath.jl")

# exact diagonalizations
include("exactdiagonalizations/exactdiagonalizations.jl")

# fourier transformations
include("fouriertransforms/fouriertransforms.jl")


# some analytical solutions
include("analyticsolutions/analyticsolutions.jl")


# utilities for real spectrum functions
include("utilities/utilities.jl")


end