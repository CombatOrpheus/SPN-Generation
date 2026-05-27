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
    # Cache locality: this is iterating rows then columns, which is cache unfriendly
    # Let's write an optimized version that iterates column-major
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

    # 3. Ensure each place has at least one outgoing edge
    # Outgoing edges are in columns (num_transitions + 1) to 2*num_transitions
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
b = @benchmarkable add_missing_connections(copy(pn), 10) setup=(pn=$pn)
display(run(b))
println("Optimized:")
b2 = @benchmarkable add_missing_connections_opt(copy(pn), 10) setup=(pn=$pn)
display(run(b2))
