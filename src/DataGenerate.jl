# Contents of DataGenerate.jl
"""
Initializes the Petri net matrix and selects the first connection.
"""
function initialize_petri_net(num_places::Int, num_transitions::Int)
    remaining_nodes = collect(1:(num_places + num_transitions))
    petri_matrix = zeros(Int32, num_places, 2 * num_transitions + 1)

    first_place = rand(1:num_places)
    first_transition = rand((num_places + 1):(num_places + num_transitions))

    filter!(x -> x != first_place, remaining_nodes)
    filter!(x -> x != first_transition, remaining_nodes)

    if rand() <= 0.5
        petri_matrix[first_place, first_transition - num_places] = 1
    else
        petri_matrix[first_place, first_transition - num_places + num_transitions] = 1
    end

    shuffle!(remaining_nodes)
    sub_graph = Vector{Int}([first_place, first_transition])
    return petri_matrix, remaining_nodes, sub_graph
end

"""
Connects the remaining nodes to the sub-graph.
"""
function connect_remaining_nodes(petri_matrix, remaining_nodes, sub_graph, num_places, num_transitions)
    for node in shuffle(remaining_nodes)
        sub_places = filter(x -> x <= num_places, sub_graph)
        sub_transitions = filter(x -> x > num_places, sub_graph)

        place, transition = if node <= num_places
            (node, rand(sub_transitions))
        else
            (rand(sub_places), node)
        end

        if rand() <= 0.5
            petri_matrix[place, transition - num_places] = 1
        else
            petri_matrix[place, transition - num_places + num_transitions] = 1
        end

        push!(sub_graph, node)
    end
    return petri_matrix
end

"""
Generates a random Petri net matrix.
"""
function generate_random_petri_net(num_places::Int, num_transitions::Int)
    petri_matrix, remaining_nodes, sub_graph = initialize_petri_net(num_places, num_transitions)
    petri_matrix = connect_remaining_nodes(petri_matrix, remaining_nodes, sub_graph, num_places, num_transitions)

    # Add an initial marking
    random_place = rand(1:num_places)
    petri_matrix[random_place, end] = 1

    return petri_matrix
end

"""
Deletes excess edges from the Petri net.
"""
function delete_excess_edges(petri_matrix, num_transitions)
    num_places = size(petri_matrix, 1)

    # Places
    place_sums = sum(petri_matrix[:, 1:end-1], dims=2)
    for i in findall(s -> s >= 3, vec(place_sums))
        edge_indices = findall(isone, petri_matrix[i, 1:end-1])
        if length(edge_indices) > 2
            indices_to_remove = shuffle(edge_indices)[1:(length(edge_indices) - 2)]
            petri_matrix[i, indices_to_remove] .= 0
        end
    end

    # Transitions
    transition_sums = sum(petri_matrix, dims=1)
    for i in findall(s -> s >= 3, transition_sums[1, 1:(2*num_transitions)])
        edge_indices = findall(isone, petri_matrix[:, i])
        if length(edge_indices) > 2
            indices_to_remove = shuffle(edge_indices)[1:(length(edge_indices) - 2)]
            petri_matrix[indices_to_remove, i] .= 0
        end
    end

    return petri_matrix
end

"""
Adds connections to ensure the Petri net is valid.
"""
function add_missing_connections(petri_matrix, num_transitions)
    num_places = size(petri_matrix, 1)

    # Ensure each transition has at least one connection
    zero_sum_cols = findall(iszero, vec(sum(view(petri_matrix, :, 1:(2*num_transitions)), dims=1)))
    if !isempty(zero_sum_cols)
        random_rows = rand(1:num_places, length(zero_sum_cols))
        petri_matrix[random_rows, zero_sum_cols] .= 1
    end

    pre_matrix = view(petri_matrix, :, 1:num_transitions)
    post_matrix = view(petri_matrix, :, (num_transitions + 1):(2*num_transitions))

    # Ensure each place has at least one incoming edge
    rows_with_zero_pre_sum = findall(iszero, vec(sum(pre_matrix, dims=2)))
    if !isempty(rows_with_zero_pre_sum)
        random_cols_pre = rand(1:num_transitions, length(rows_with_zero_pre_sum))
        petri_matrix[rows_with_zero_pre_sum, random_cols_pre] .= 1
    end

    # Ensure each place has at least one outgoing edge
    rows_with_zero_post_sum = findall(iszero, vec(sum(post_matrix, dims=2)))
    if !isempty(rows_with_zero_post_sum)
        random_cols_post = rand(1:num_transitions, length(rows_with_zero_post_sum))
        petri_matrix[rows_with_zero_post_sum, random_cols_post .+ num_transitions] .= 1
    end

    return petri_matrix
end

"""
Prunes a Petri net by deleting edges and adding nodes.
"""
function prune_petri_net(petri_matrix)
    num_transitions = (size(petri_matrix, 2) - 1) ÷ 2
    petri_matrix = delete_excess_edges(petri_matrix, num_transitions)
    petri_matrix = add_missing_connections(petri_matrix, num_transitions)
    return petri_matrix
end

"""
Adds tokens to random places in the Petri net.
"""
function add_tokens_randomly(petri_matrix)
    num_places = size(petri_matrix, 1)
    random_values = rand(0:9, num_places)
    petri_matrix[:, end] .+= (random_values .<= 2)
    return petri_matrix
end
