# Brute-force exact diagonalization on the full Fock space: matrix
# representations of the fermionic/bosonic operators and the many-body density
# matrices (dimension 2^L / d^L), plus full-space correlation functions and
# Lindblad dynamics. Only for benchmarking the fast quadratic-systems methods
# in src/quadraticsystems.

include("operators.jl")

include("correlations.jl")
