using BenchmarkTools
using Random
include("src/DataTransformation.jl")
include("src/DataGenerate.jl")
include("src/ArrivableGraph.jl")
include("src/SPN.jl")
include("src/Utils.jl")
using .SPNGenerator # Wait, SPNGenerator is not imported

# Let's test the above optimization!
