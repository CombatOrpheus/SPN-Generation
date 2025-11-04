# Contents of ArrivableGraph submodule
"""
Identifies enabled transitions and calculates the resulting markings.
"""
function get_enabled_transitions(pre_condition_matrix, change_matrix, current_marking_vector)
    enabled_mask = all(current_marking_vector .>= pre_condition_matrix, dims=1)
    enabled_transitions = findall(vec(enabled_mask))

    if isempty(enabled_transitions)
        num_places = size(pre_condition_matrix, 1)
        return Matrix{Int64}(undef, 0, num_places), Vector{Int64}(undef, 0)
    end

    new_markings = current_marking_vector .+ change_matrix[:, enabled_transitions]
    return new_markings', enabled_transitions
end

function initialize_bfs(initial_marking)
    marking_index_counter = 1 # Julia is 1-indexed
    visited_markings_list = Vector{Vector{Int}}([initial_marking])
    explored_markings_dict = Dict{Vector{Int}, Int}(initial_marking => marking_index_counter)
    processing_queue = Queue{Int}()
    push!(processing_queue, marking_index_counter)
    return marking_index_counter, visited_markings_list, explored_markings_dict, processing_queue
end

function process_marking(
    current_marking_index,
    visited_markings_list,
    pre_matrix,
    change_matrix,
    place_upper_limit,
    max_markings_to_explore,
)
    current_marking = visited_markings_list[current_marking_index]

    if length(visited_markings_list) >= max_markings_to_explore
        return nothing, nothing, true
    end

    enabled_next_markings, enabled_transition_indices = get_enabled_transitions(
        pre_matrix, change_matrix, current_marking
    )

    if !isempty(enabled_next_markings) && any(enabled_next_markings .> place_upper_limit)
        return nothing, nothing, true
    end

    return enabled_next_markings, enabled_transition_indices, false
end

function update_graph!(
    new_marking,
    enabled_transition_index,
    current_marking_index,
    marking_index_counter,
    visited_markings_list,
    explored_markings_dict,
    processing_queue,
    reachability_edges,
    edge_transition_indices,
    max_markings_to_explore,
)
    if !haskey(explored_markings_dict, new_marking)
        marking_index_counter += 1
        if marking_index_counter >= max_markings_to_explore
            push!(reachability_edges, [current_marking_index, marking_index_counter])
            push!(edge_transition_indices, enabled_transition_index)
            return marking_index_counter, true
        end

        push!(visited_markings_list, new_marking)
        explored_markings_dict[new_marking] = marking_index_counter
        push!(processing_queue, marking_index_counter)
        push!(reachability_edges, [current_marking_index, marking_index_counter])
    else
        existing_index = explored_markings_dict[new_marking]
        push!(reachability_edges, [current_marking_index, existing_index])
    end

    push!(edge_transition_indices, enabled_transition_index)
    return marking_index_counter, false
end

function generate_reachability_graph(incidence_matrix_with_initial; place_upper_limit=10, max_markings_to_explore=500)
    incidence_matrix = incidence_matrix_with_initial
    num_transitions = size(incidence_matrix, 2) ÷ 2
    pre_matrix = view(incidence_matrix, :, 1:num_transitions)
    post_matrix = view(incidence_matrix, :, (num_transitions + 1):(2*num_transitions))
    initial_marking = Vector{Int}(incidence_matrix[:, end])
    change_matrix = post_matrix - pre_matrix

    (
        marking_index_counter,
        visited_markings_list,
        explored_markings_dict,
        processing_queue,
    ) = initialize_bfs(initial_marking)

    reachability_edges = Vector{Vector{Int}}()
    edge_transition_indices = Vector{Int}()
    is_bounded = true

    while !isempty(processing_queue)
        current_marking_index = popfirst!(processing_queue)

        enabled_next_markings, enabled_transition_indices, stop_exploration = process_marking(
            current_marking_index,
            visited_markings_list,
            pre_matrix,
            change_matrix,
            place_upper_limit,
            max_markings_to_explore,
        )

        if stop_exploration
            is_bounded = false
            break
        end

        if isnothing(enabled_next_markings)
            continue
        end

        for i in 1:size(enabled_next_markings, 1)
            new_marking = enabled_next_markings[i, :]
            enabled_transition_index = enabled_transition_indices[i]

            marking_index_counter, stop = update_graph!(
                new_marking,
                enabled_transition_index,
                current_marking_index,
                marking_index_counter,
                visited_markings_list,
                explored_markings_dict,
                processing_queue,
                reachability_edges,
                edge_transition_indices,
                max_markings_to_explore,
            )
            if stop
                is_bounded = false
                break
            end
        end
        if !is_bounded
            break
        end
    end

    return (
        visited_markings_list,
        reachability_edges,
        edge_transition_indices,
        num_transitions,
        is_bounded,
    )
end
