# Contents of DataGenerate.jl
"""
Initializes the Petri net matrix and selects the first connection.
"""
function _initialize_petri_net(num_places::Int, num_transitions::Int)
    remaining_nodes = collect(1:(num_places + num_transitions))
    petri_matrix = zeros(Int32, num_places, 2 * num_transitions + 1)

    first_place = rand(1:num_places)
    first_transition = rand((num_places + 1):(num_places + num_transitions))

    filter!(x -> x != first_place, remaining_nodes)
    filter!(x -> x != first_transition, remaining_nodes)

    if rand(Bool)
        petri_matrix[first_place, first_transition - num_places] = 1
    else
        petri_matrix[first_place, first_transition - num_places + num_transitions] = 1
    end

    shuffle!(remaining_nodes)
    # Optimize: Maintain explicit typed collections to avoid filter! inside loops
    sub_places = Int[first_place]
    sub_transitions = Int[first_transition]
    return petri_matrix, remaining_nodes, sub_places, sub_transitions
end

"""
Connects the remaining nodes to the sub-graph.
"""
function _connect_remaining_nodes(petri_matrix, remaining_nodes, sub_places, sub_transitions, num_places, num_transitions)
    # Optimize: remaining_nodes is already shuffled in _initialize_petri_net, avoid shuffle allocation here
    for node in remaining_nodes
        # Optimize: eliminate intermediate sub_places/sub_transitions filter array allocations

        place, transition = if node <= num_places
            (node, rand(sub_transitions))
        else
            (rand(sub_places), node)
        end

        if rand(Bool)
            petri_matrix[place, transition - num_places] = 1
        else
            petri_matrix[place, transition - num_places + num_transitions] = 1
        end

        # Dynamically push to typed collection instead of single sub_graph
        if node <= num_places
            push!(sub_places, node)
        else
            push!(sub_transitions, node)
        end
    end
    return petri_matrix
end

"""
Generates a random Petri net matrix.
"""
function generate_random_petri_net(num_places::Int, num_transitions::Int)
    petri_matrix, remaining_nodes, sub_places, sub_transitions = _initialize_petri_net(num_places, num_transitions)
    petri_matrix = _connect_remaining_nodes(petri_matrix, remaining_nodes, sub_places, sub_transitions, num_places, num_transitions)

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
    num_cols = 2 * num_transitions

    # Pre-allocate array for indices
    edge_indices = Vector{Int}(undef, max(num_places, num_cols))

    # Places
    @inbounds for i in 1:num_places
        edge_count = 0
        for j in 1:num_cols
            if petri_matrix[i, j] == 1
                edge_count += 1
                edge_indices[edge_count] = j
            end
        end

        if edge_count >= 3
            # We need to remove edge_count - 2 edges randomly
            num_to_remove = edge_count - 2
            # Shuffle only the part we care about to pick `num_to_remove` items
            for k in 1:num_to_remove
                idx_to_swap = rand(k:edge_count)
                # Swap
                tmp = edge_indices[k]
                edge_indices[k] = edge_indices[idx_to_swap]
                edge_indices[idx_to_swap] = tmp

                # Remove edge
                petri_matrix[i, edge_indices[k]] = 0
            end
        end
    end

    # Transitions
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

"""
Adds connections to ensure the Petri net is valid.
"""
function add_missing_connections(petri_matrix, num_transitions)
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
    @inbounds for i in 1:num_places
        has_incoming = false
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

    # 3. Ensure each place has at least one outgoing edge
    # Outgoing edges are in columns (num_transitions + 1) to 2*num_transitions
    @inbounds for i in 1:num_places
        has_outgoing = false
        for j in (num_transitions + 1):(2*num_transitions)
            if petri_matrix[i, j] != 0
                has_outgoing = true
                break
            end
        end
        if !has_outgoing
            random_col = rand(1:num_transitions) + num_transitions
            petri_matrix[i, random_col] = 1
        end
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

    @inbounds for i in 1:num_places
        if rand() <= 0.3
            petri_matrix[i, end] += 1
        end
    end
    return petri_matrix
end
