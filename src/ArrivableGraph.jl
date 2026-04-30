# Contents of ArrivableGraph submodule
"""
Identifies enabled transitions and calculates the resulting markings.
"""
function get_enabled_transitions!(pre_condition_matrix, change_matrix, current_marking_vector, enabled_transitions_buffer, new_markings_buffer)
    # Optimize hot loop: use nested loops with early break instead of allocating vectorized operations
    num_places, num_transitions = size(pre_condition_matrix)

    num_enabled = 0
    @inbounds for j in 1:num_transitions
        enabled = true
        for i in 1:num_places
            if current_marking_vector[i] < pre_condition_matrix[i, j]
                enabled = false
                break
            end
        end
        if enabled
            num_enabled += 1
            enabled_transitions_buffer[num_enabled] = j
        end
    end

    if num_enabled == 0
        return view(new_markings_buffer, 1:0, 1:num_places), view(enabled_transitions_buffer, 1:0)
    end

    @inbounds for j in 1:num_enabled
        trans_idx = enabled_transitions_buffer[j]
        for i in 1:num_places
            new_markings_buffer[j, i] = current_marking_vector[i] + change_matrix[i, trans_idx]
        end
    end

    return view(new_markings_buffer, 1:num_enabled, :), view(enabled_transitions_buffer, 1:num_enabled)
end

function _initialize_bfs(initial_marking)
    marking_index_counter = 1 # Julia is 1-indexed
    visited_markings_list = Vector{Vector{Int}}([initial_marking])
    explored_markings_dict = Dict{Vector{Int}, Int}(initial_marking => marking_index_counter)
    processing_queue = Int[]
    # ⚡ Bolt Optimization: Replace DataStructures.Queue with a simple Vector
    # to significantly reduce memory allocations and queue operations overhead in BFS.
    push!(processing_queue, marking_index_counter)
    return marking_index_counter, visited_markings_list, explored_markings_dict, processing_queue
end

function _process_marking(
    current_marking_index,
    visited_markings_list,
    pre_matrix,
    change_matrix,
    place_upper_limit,
    max_markings_to_explore,
    enabled_transitions_buffer,
    new_markings_buffer
)
    current_marking = visited_markings_list[current_marking_index]

    if length(visited_markings_list) >= max_markings_to_explore
        return nothing, nothing, true
    end

    enabled_next_markings, enabled_transition_indices = get_enabled_transitions!(
        pre_matrix, change_matrix, current_marking, enabled_transitions_buffer, new_markings_buffer
    )

    exceeds_limit = false
    if !isempty(enabled_next_markings)
        num_enabled = size(enabled_next_markings, 1)
        num_places = size(enabled_next_markings, 2)
        @inbounds for j in 1:num_places
            for i in 1:num_enabled
                if enabled_next_markings[i, j] > place_upper_limit
                    exceeds_limit = true
                    break
                end
            end
            if exceeds_limit
                break
            end
        end
    end

    if exceeds_limit
        return nothing, nothing, true
    end

    return enabled_next_markings, enabled_transition_indices, false
end

function _update_graph!(
    new_marking_view,
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
    if !haskey(explored_markings_dict, new_marking_view)
        marking_index_counter += 1
        if marking_index_counter >= max_markings_to_explore
            push!(reachability_edges, (current_marking_index, marking_index_counter))
            push!(edge_transition_indices, enabled_transition_index)
            return marking_index_counter, true
        end

        new_marking_vector = Vector{Int}(new_marking_view)
        push!(visited_markings_list, new_marking_vector)
        explored_markings_dict[new_marking_vector] = marking_index_counter
        push!(processing_queue, marking_index_counter)
        push!(reachability_edges, (current_marking_index, marking_index_counter))
    else
        existing_index = explored_markings_dict[new_marking_view]
        push!(reachability_edges, (current_marking_index, existing_index))
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
    ) = _initialize_bfs(initial_marking)

    sizehint!(visited_markings_list, max_markings_to_explore)
    sizehint!(explored_markings_dict, max_markings_to_explore)
    sizehint!(processing_queue, max_markings_to_explore)

    reachability_edges = Vector{Tuple{Int, Int}}()
    edge_transition_indices = Vector{Int}()

    # We estimate edges to be about 2-3x the number of vertices for sparse graphs
    sizehint!(reachability_edges, max_markings_to_explore * 2)
    sizehint!(edge_transition_indices, max_markings_to_explore * 2)

    is_bounded = true

    num_places = size(incidence_matrix, 1)
    enabled_transitions_buffer = Vector{Int}(undef, num_transitions)
    new_markings_buffer = Matrix{Int}(undef, num_transitions, num_places)

    queue_head = 1
    while queue_head <= length(processing_queue)
        current_marking_index = processing_queue[queue_head]
        queue_head += 1

        enabled_next_markings, enabled_transition_indices, stop_exploration = _process_marking(
            current_marking_index,
            visited_markings_list,
            pre_matrix,
            change_matrix,
            place_upper_limit,
            max_markings_to_explore,
            enabled_transitions_buffer,
            new_markings_buffer
        )

        if stop_exploration
            is_bounded = false
            break
        end

        if isnothing(enabled_next_markings)
            continue
        end

        for i in 1:size(enabled_next_markings, 1)
            new_marking_view = view(enabled_next_markings, i, :)
            enabled_transition_index = enabled_transition_indices[i]

            marking_index_counter, stop = _update_graph!(
                new_marking_view,
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
