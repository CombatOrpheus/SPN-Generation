using BenchmarkTools
using Random
include("src/DataTransformation.jl")

function generate_petri_net_variations_opt(petri_matrix, config)
    base_petri_matrix = petri_matrix
    candidate_matrices = _generate_candidate_matrices(base_petri_matrix, config)

    max_candidates = convert(Int, get(config, "max_candidates_per_structure", 50))
    if length(candidate_matrices) > max_candidates
        len = length(candidate_matrices)
        # Random subset selection via Fisher-Yates shuffle
        # We only need to pick max_candidates.
        @inbounds for i in 1:max_candidates
            j = rand(i:len)
            tmp = candidate_matrices[i]
            candidate_matrices[i] = candidate_matrices[j]
            candidate_matrices[j] = tmp
        end
        resize!(candidate_matrices, max_candidates)
    end
end
