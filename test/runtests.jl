using Test, Random

# push!(LOAD_PATH, dirname(Base.@__DIR__) * "/src")
# using ImpurityModelBase

include("../src/includes.jl")

Random.seed!(12354)

# util
include("ed.jl")


include("auxiliary.jl")
include("fourier.jl")
include("freefermion.jl")
include("exactsolutions.jl")

include("exactdiagonalizations.jl")

include("interactingbosons/util.jl")
include("interactingbosons/noninteracting.jl")
include("interactingbosons/interacting.jl")



### include("holstein.jl")