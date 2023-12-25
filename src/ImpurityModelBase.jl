module ImpurityModelBase

export AbstractFermionicBath, FermionicBath, FermionicVacuum, fermionicbath, thermaloccupation
export SpectrumFunction, lowerbound, upperbound

include("spectrumfunc.jl")
include("bath.jl")

end