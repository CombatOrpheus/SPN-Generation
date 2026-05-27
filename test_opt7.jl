using BenchmarkTools
using Random
include("src/DataGenerate.jl")

function delete_excess_edges_opt(petri_matrix, num_transitions)
    num_places = size(petri_matrix, 1)
    num_cols = 2 * num_transitions

    edge_indices = Vector{Int}(undef, max(num_places, num_cols))

    @inbounds for i in 1:num_places
        edge_count = 0
        for j in 1:num_cols
            if petri_matrix[i, j] == 1
                edge_count += 1
                edge_indices[edge_count] = j
            end
        end

        if edge_count >= 3
            num_to_remove = edge_count - 2
            # Use partial Fisher-Yates shuffle directly in the loop to avoid allocating intermediate subsets
            # or picking randomly via repeatedly checking items.
            for k in 1:num_to_remove
                idx_to_swap = rand(k:edge_count)
                tmp = edge_indices[k]
                edge_indices[k] = edge_indices[idx_to_swap]
                edge_indices[idx_to_swap] = tmp
                petri_matrix[i, edge_indices[k]] = 0
            end
        end
    end

    @inbounds for j in 1:num_cols
        edge_count = 0
        for i in 1:num_places
            if petri_matrix[i, j] == 1
                edge_count += 1
                edge_indices[edge_count] = i
            end
        end

        if edge_count >= 3
            num_to_remove = edge_count - 2
            for k in 1:num_to_remove
                idx_to_swap = rand(k:edge_count)
                tmp = edge_indices[k]
                edge_indices[k] = edge_indices[idx_to_swap]
                edge_indices[idx_to_swap] = tmp
                petri_matrix[edge_indices[k], j] = 0
            end
        end
    end

    return petri_matrix
end

pn = generate_random_petri_net(10, 10)
println("Original:")
b = @benchmarkable delete_excess_edges(copy(pn), 10) setup=(pn=$pn)
display(run(b))
println("Optimized:")
b2 = @benchmarkable delete_excess_edges_opt(copy(pn), 10) setup=(pn=$pn)
display(run(b2))
