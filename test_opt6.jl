using BenchmarkTools
using Random
include("src/DataGenerate.jl")

function prune_petri_net_opt(petri_matrix)
    num_transitions = (size(petri_matrix, 2) - 1) ÷ 2
    petri_matrix = delete_excess_edges(petri_matrix, num_transitions)

    # Original logic of add_missing_connections optimized directly here
    num_places = size(petri_matrix, 1)

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

    has_incoming = falses(num_places)
    @inbounds for j in 1:num_transitions
        for i in 1:num_places
            if petri_matrix[i, j] != 0
                has_incoming[i] = true
            end
        end
    end

    @inbounds for i in 1:num_places
        if !has_incoming[i]
            random_col = rand(1:num_transitions)
            petri_matrix[i, random_col] = 1
        end
    end

    has_outgoing = falses(num_places)
    @inbounds for j in (num_transitions + 1):(2*num_transitions)
        for i in 1:num_places
            if petri_matrix[i, j] != 0
                has_outgoing[i] = true
            end
        end
    end

    @inbounds for i in 1:num_places
        if !has_outgoing[i]
            random_col = rand(1:num_transitions) + num_transitions
            petri_matrix[i, random_col] = 1
        end
    end

    return petri_matrix
end

pn = generate_random_petri_net(10, 10)
println("Original:")
b = @benchmarkable prune_petri_net(copy(pn)) setup=(pn=$pn)
display(run(b))
println("Optimized:")
b2 = @benchmarkable prune_petri_net_opt(copy(pn)) setup=(pn=$pn)
display(run(b2))
