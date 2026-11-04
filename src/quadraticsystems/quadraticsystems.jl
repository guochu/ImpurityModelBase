# Fast solutions of quadratic (free-particle) Hamiltonians: everything based on
# the single-particle coefficient matrix (cmatrix) and coefficient density
# matrix (cdm) — Hamiltonian containers, thermal cdm states, time evolution,
# Green's functions and the quadratic Toulouse / boundary-driven models. In
# contrast to the brute-force full-Fock-space exact diagonalization in
# src/exactdiagonalizations, which serves as the benchmark reference.

include("hamiltonians.jl")

include("cache.jl")

include("bcscache.jl")

# discrete baths as quadratic systems: hamiltonian / cmatrix / thermocdm
include("baths.jl")

include("greenfunctions.jl")

include("contourgreenfunctions.jl")

include("toulouse.jl")

include("bogoliubov.jl")

include("boundarydriving.jl")
