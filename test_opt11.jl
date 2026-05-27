using BenchmarkTools
using Random
include("src/DataGenerate.jl")

function _generate_candidate_matrices_core_opt2(
    base_petri_matrix,
    enable_delete_edge,
    enable_add_edge,
    enable_add_token,
    enable_delete_token,
)
    num_places, num_cols = size(base_petri_matrix)

    # Optimization: Calculate exact number of matrices to pre-allocate
    # Avoid count and views to prevent small array allocations in count(view(...))
    num_delete_edge = 0
    num_add_edge = 0

    # We can do this in one pass over the matrix without allocating views.
    if enable_delete_edge || enable_add_edge
        # We only look at columns 1 to num_cols-1
        @inbounds for j in 1:(num_cols-1)
            for i in 1:num_places
                if base_petri_matrix[i, j] == 1
                    num_delete_edge += 1
                else
                    num_add_edge += 1
                end
            end
        end
    end

    if !enable_delete_edge
        num_delete_edge = 0
    end
    if !enable_add_edge
        num_add_edge = 0
    end

    num_add_token = enable_add_token ? num_places : 0

    token_sum = 0
    num_delete_token = 0
    if enable_delete_token
        @inbounds for i in 1:num_places
            val = base_petri_matrix[i, end]
            token_sum += val
            if val >= 1
                num_delete_token += 1
            end
        end
        if token_sum <= 1
            num_delete_token = 0
        end
    end

    total_matrices = num_delete_edge + num_add_edge + num_add_token + num_delete_token

    # Typed array avoids Any allocations
    candidate_matrices = Vector{typeof(base_petri_matrix)}(undef, total_matrices)
    idx = 1

    # Optimization: Loop column-major without allocating views or using findall
    if enable_delete_edge
        @inbounds for j in 1:num_cols-1
            for i in 1:num_places
                if base_petri_matrix[i, j] == 1
                    modified_matrix = copy(base_petri_matrix)
                    modified_matrix[i, j] = 0
                    candidate_matrices[idx] = modified_matrix
                    idx += 1
                end
            end
        end
    end

    if enable_add_edge
        @inbounds for j in 1:num_cols-1
            for i in 1:num_places
                if base_petri_matrix[i, j] == 0
                    modified_matrix = copy(base_petri_matrix)
                    modified_matrix[i, j] = 1
                    candidate_matrices[idx] = modified_matrix
                    idx += 1
                end
            end
        end
    end

    if enable_add_token
        @inbounds for i in 1:num_places
            modified_matrix = copy(base_petri_matrix)
            modified_matrix[i, end] += 1
            candidate_matrices[idx] = modified_matrix
            idx += 1
        end
    end

    if enable_delete_token && token_sum > 1
        @inbounds for i in 1:num_places
            if base_petri_matrix[i, end] >= 1
                modified_matrix = copy(base_petri_matrix)
                modified_matrix[i, end] -= 1
                candidate_matrices[idx] = modified_matrix
                idx += 1
            end
        end
    end

    return candidate_matrices
end
