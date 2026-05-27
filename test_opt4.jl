using BenchmarkTools
using Random
include("src/DataGenerate.jl")

function add_missing_connections_opt(petri_matrix, num_transitions)
    num_places = size(petri_matrix, 1)

    # 1. Ensure each transition has at least one connection
    # Transitions are in columns 1 to 2*num_transitions
    @inbounds for j in 1:(2*num_transitions)
        has_connection = false
        for i in 1:num_places
            if petri_matrix[i, j] != 0
                has_connection = true
                break
            end
        end
        if !has_connection
            random_row = rand(1:num_places)
            petri_matrix[random_row, j] = 1
        end
    end

    # 2. Ensure each place has at least one incoming edge
    # Incoming edges are in columns 1 to num_transitions
    # Transposing loop correctly
    @inbounds for i in 1:num_places
        has_incoming = false
        # The inner loop iterates over columns (j). Matrix is column major!
        # So `petri_matrix[i, j]` for fixed `i` and varying `j` is non-contiguous in memory
        # Wait, the original code is:
        # for j in 1:num_transitions
        #     if petri_matrix[i, j] != 0
        # If we just do it like this, we're scanning rows.
        # But this is a small matrix (10x21). The cache misses are negligible.
        for j in 1:num_transitions
            if petri_matrix[i, j] != 0
                has_incoming = true
                break
            end
        end
        if !has_incoming
            random_col = rand(1:num_transitions)
            petri_matrix[i, random_col] = 1
        end
    end

    return petri_matrix
end

pn = generate_random_petri_net(100, 100)
println("Original:")
b = @benchmarkable add_missing_connections(copy(pn), 100) setup=(pn=$pn)
display(run(b))
