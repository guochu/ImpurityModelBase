using Test, Random

push!(LOAD_PATH, dirname(Base.@__DIR__) * "/src")
using ImpurityModelBase


Random.seed!(12354)

### auxiliary
include("auxiliary.jl")
