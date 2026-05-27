using BenchmarkTools
using Random
include("src/DataGenerate.jl")

function delete_excess_edges_opt(petri_matrix, num_transitions)
    num_places = size(petri_matrix, 1)
    num_cols = 2 * num_transitions

    # Avoid Vector{Int} allocation for edge indices
    # We can pre-allocate it outside or just use a small tuple/static array for indices since we only keep a few? No, a node might have all edges connected.
    # We can do this in-place! No, we don't need a separate array. We can just keep the indices, but we need to pick random elements.
    # In Julia, pre-allocating an array per function call is still an allocation.

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
            # Instead of a loop of rands, do a partial Fisher-Yates
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

pn = generate_random_petri_net(100, 100)
println("Original:")
b = @benchmarkable delete_excess_edges(copy(pn), 100) setup=(pn=$pn)
display(run(b))
println("Optimized:")
b2 = @benchmarkable delete_excess_edges_opt(copy(pn), 100) setup=(pn=$pn)
display(run(b2))
