using Test, Random
using ImpurityModelBase

Random.seed!(12354)

# The test suite is organized into 5 parts, mirroring the source code layout
# (src/spectrumfuncs, src/baths, src/exactdiagonalizations,
# src/analyticsolutions, src/utilities):
#   1. spectrumfuncs          — spectrum density wrappers and predefined spectra
#   2. baths                  — particle types, thermal distributions, bath containers
#   3. exactdiagonalizations  — exact numerical solutions; the (efficient)
#                               coefficient-matrix method is cross-checked against
#                               the (debug-only) full-Hamiltonian method
#   4. analyticsolutions      — analytical solutions, benchmarked against the
#                               full-Hamiltonian ED reference
#   5. utilities              — Fourier transforms and linear prediction

# 1. spectrum functions
include("spectrumfuncs/spectrumfuncs.jl")

# 2. baths
include("baths/baths.jl")

# 3. exact diagonalizations
# util.jl provides the shared reference helpers (full-Hamiltonian Green's
# functions) used both here and by the analytical-solution benchmarks below
include("exactdiagonalizations/exactdiagonalizations.jl")

# 4. analytical solutions
include("analyticsolutions/freefermion.jl")
include("analyticsolutions/exactsolutions.jl")
include("analyticsolutions/independentbosons/independentbosons.jl")
include("analyticsolutions/dephasing.jl")
# include("analyticsolutions/holstein.jl")

# 5. utilities
include("utilities/fourier.jl")
include("utilities/linearprediction.jl")
