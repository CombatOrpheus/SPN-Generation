using BenchmarkTools
using Random
include("src/DataTransformation.jl")
include("src/ArrivableGraph.jl")
include("src/SPN.jl")

# From _generate_candidate_matrices_core optimization we see that it allocates:
#     candidate_matrices = Vector{typeof(base_petri_matrix)}(undef, total_matrices)
# Wait, total_matrices is correct and we don't push!
# Then in _generate_candidate_matrices we do:
#    if convert(Bool, get(config, "enable_add_place", false)) && num_transitions > 0
#        ...
#        push!(candidate_matrices, modified_matrix)
