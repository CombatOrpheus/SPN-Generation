module SPNGenerator

using Random
using DataStructures
using SparseArrays
using HDF5
using IterativeSolvers
using LinearAlgebra
using TOML
using JSON3
using Base.Threads
using ProgressMeter

# Exports from DataGenerate and its submodules
export generate_random_petri_net, prune_petri_net, add_tokens_randomly
export generate_reachability_graph
export filter_spn, get_spn_info, generate_single_spn, augment_single_spn
export generate_petri_net_variations, generate_lambda_variations

# Exports from Utils
export create_directory, load_toml_file, save_data_to_json_file, load_json_file, count_json_files, sample_json_files_from_directory
export write_to_hdf5, write_to_jsonl, run_generation_from_config

function generate_single_spn(config)
    max_attempts = 100
    for _ in 1:max_attempts
        place_num = rand(config["minimum_number_of_places"]::Int : config["maximum_number_of_places"]::Int)
        trans_num = rand(config["minimum_number_of_transitions"]::Int : config["maximum_number_of_transitions"]::Int)

        petri_matrix = generate_random_petri_net(place_num, trans_num)
        if get(config, "enable_pruning", false)::Bool
            petri_matrix = prune_petri_net(petri_matrix)
        end
        if get(config, "enable_token_addition", false)::Bool
            petri_matrix = add_tokens_randomly(petri_matrix)
        end

        # Convert to sparse if it's large enough
        if place_num >= 500
            petri_matrix = sparse(petri_matrix)
        end

        results, success = filter_spn(
            petri_matrix,
            place_upper_bound=config["place_upper_bound"]::Int,
            marks_lower_limit=config["marks_lower_limit"]::Int,
            marks_upper_limit=config["marks_upper_limit"]::Int,
        )
        if success
            return results
        end
    end
    return nothing
end

function _augment_single_spn_impl(sample, petri_net::AbstractMatrix{Int32}, config)
    augmented_data::Vector{Dict{String, Any}} = generate_lambda_variations(
        sample,
        petri_net,
        convert(Int, config["lambda_variations_per_sample"]),
    )

    if isempty(augmented_data)
        return Dict{String, Any}[]
    end

    max_transforms = convert(Int, get(config, "maximum_transformations_per_sample", length(augmented_data)))
    if length(augmented_data) > max_transforms
        return sample(augmented_data, convert(Int, max_transforms), replace=false)
    end
    return augmented_data
end

function augment_single_spn(sample, config)
    if isnothing(sample) || !haskey(sample, "petri_net")
        return Dict{String, Any}[]
    end

    petri_net = sample["petri_net"]::AbstractMatrix{Int32}
    return _augment_single_spn_impl(sample, petri_net, config)
end

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

    if rand() <= 0.5
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

        if rand() <= 0.5
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
    processing_queue = Queue{Int}()
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
    ) = _initialize_bfs(initial_marking)

    reachability_edges = Vector{Vector{Int}}()
    edge_transition_indices = Vector{Int}()
    is_bounded = true

    num_places = size(incidence_matrix, 1)
    enabled_transitions_buffer = Vector{Int}(undef, num_transitions)
    new_markings_buffer = Matrix{Int}(undef, num_transitions, num_places)

    while !isempty(processing_queue)
        current_marking_index = popfirst!(processing_queue)

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
            new_marking = enabled_next_markings[i, :]
            enabled_transition_index = enabled_transition_indices[i]

            marking_index_counter, stop = _update_graph!(
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


# Contents of SPN submodule
function _compute_state_equation_core(num_vertices, edges, arc_transitions, lambda_values)
    num_edges = length(edges)
    I = Vector{Int}(undef, 2 * num_edges + num_vertices)
    J = Vector{Int}(undef, 2 * num_edges + num_vertices)
    V = Vector{Float64}(undef, 2 * num_edges + num_vertices)

    idx = 1
    for i in 1:num_edges
        edge = edges[i]
        trans_idx = arc_transitions[i]
        src_idx, dest_idx = edge[1], edge[2]
        rate = lambda_values[trans_idx]

        I[idx] = src_idx
        J[idx] = src_idx
        V[idx] = -rate
        idx += 1

        I[idx] = dest_idx
        J[idx] = src_idx
        V[idx] = rate
        idx += 1
    end

    for j in 1:num_vertices
        I[idx] = num_vertices + 1
        J[idx] = j
        V[idx] = 1.0
        idx += 1
    end

    return sparse(I, J, V, num_vertices + 1, num_vertices)
end

function compute_state_equation(vertices, edges, arc_transitions, lambda_values)
    num_vertices = length(vertices)
    state_matrix = _compute_state_equation_core(num_vertices, edges, arc_transitions, lambda_values)
    target_vector = zeros(Float64, num_vertices + 1)
    target_vector[end] = 1.0
    return state_matrix, target_vector
end

function compute_average_markings(vertices::Matrix{Int}, steady_state_probs::Vector{Float64})
    num_states, num_places = size(vertices)

    unique_token_values = sort(unique(vertices))
    num_unique = length(unique_token_values)

    token_to_idx = Dict{Int, Int}()
    for (idx, val) in enumerate(unique_token_values)
        token_to_idx[val] = idx
    end

    marking_density_matrix = zeros(Float64, num_places, num_unique)
    avg_tokens_per_place = zeros(Float64, num_places)

    @inbounds for j in 1:num_places
        for i in 1:num_states
            prob = steady_state_probs[i]
            val = vertices[i, j]
            idx = token_to_idx[val]
            marking_density_matrix[j, idx] += prob
            avg_tokens_per_place[j] += val * prob
        end
    end

    return marking_density_matrix, avg_tokens_per_place
end

function solve_for_steady_state(state_matrix, target_vector)
    try
        probs, history = lsmr(state_matrix, target_vector, atol=1e-6, btol=1e-6, conlim=1e7, maxiter=100 * size(state_matrix, 2), log=true)
        if history.isconverged
            probs[probs .< 0] .= 0
            prob_sum = sum(probs)
            if prob_sum > 1e-9
                return probs ./ prob_sum
            end
        end
    catch e
        # Handle potential numerical issues
    end
    return nothing
end

function _vertices_to_matrix(vertices::AbstractVector{<:AbstractVector{Int}})
    num_states = length(vertices)
    num_places = length(vertices[1])
    mat = Matrix{Int}(undef, num_states, num_places)
    @inbounds for j in 1:num_places
        for i in 1:num_states
            mat[i, j] = vertices[i][j]
        end
    end
    return mat
end

function _edges_to_matrix(edges::AbstractVector{<:AbstractVector{Int}})
    if isempty(edges)
        return Matrix{Int}(undef, 0, 2)
    end
    n = length(edges)
    mat = Matrix{Int}(undef, n, 2)
    @inbounds for i in 1:n
        mat[i, 1] = edges[i][1]
        mat[i, 2] = edges[i][2]
    end
    return mat
end

function _run_spn_task(vertices, edges, arc_transitions, transition_rates)
    if isempty(vertices)
        return nothing, nothing, nothing, false
    end

    vertices_np = _vertices_to_matrix(vertices)
    state_matrix, target_vector = compute_state_equation(vertices, edges, arc_transitions, transition_rates)
    steady_state_probs = solve_for_steady_state(state_matrix, target_vector)

    if isnothing(steady_state_probs)
        return nothing, nothing, nothing, false
    end

    marking_density, avg_markings = compute_average_markings(vertices_np, steady_state_probs)
    return steady_state_probs, marking_density, avg_markings, true
end

function generate_stochastic_net_task(vertices, edges, arc_transitions, num_transitions)
    transition_rates = rand(1:10, num_transitions)
    probs, density, markings, success = _run_spn_task(vertices, edges, arc_transitions, transition_rates)
    return probs, density, markings, transition_rates, success
end

function generate_stochastic_net_task_with_rates(vertices, edges, arc_transitions, transition_rates)
    return _run_spn_task(vertices, edges, arc_transitions, transition_rates)
end

function is_connected(petri_net_matrix)
    if isempty(petri_net_matrix) || ndims(petri_net_matrix) != 2
        return false
    end
    num_places, num_cols = size(petri_net_matrix)
    if num_places == 0 || num_cols < 3
        return false
    end
    num_transitions = (num_cols - 1) ÷ 2
    if num_transitions == 0
        return false
    end

    # Check if any place has 0 connections
    @inbounds for i in 1:num_places
        has_connection = false
        for j in 1:(2*num_transitions)
            if petri_net_matrix[i, j] != 0
                has_connection = true
                break
            end
        end
        if !has_connection
            return false
        end
    end

    # Check if any transition has 0 connections (both pre and post)
    # Using column-major traversal for better cache locality
    @inbounds for j in 1:num_transitions
        has_connection = false
        for i in 1:num_places
            if petri_net_matrix[i, j] != 0 || petri_net_matrix[i, j + num_transitions] != 0
                has_connection = true
                break
            end
        end
        if !has_connection
            return false
        end
    end

    return true
end

function _create_spn_result_dict(petri_net_matrix, vertices, edges, arc_transitions, firing_rates, steady_state_probs, marking_densities, average_markings)
    return Dict{String, Any}(
        "petri_net" => petri_net_matrix,
        "arr_vlist" => _vertices_to_matrix(vertices),
        "arr_edge" => _edges_to_matrix(edges),
        "arr_tranidx" => isempty(arc_transitions) ? Vector{Int}() : arc_transitions,
        "spn_labda" => firing_rates,
        "spn_steadypro" => steady_state_probs,
        "spn_markdens" => marking_densities,
        "spn_allmus" => average_markings,
        "spn_mu" => sum(average_markings),
    )
end

function filter_spn(petri_net_matrix; place_upper_bound=10, marks_lower_limit=4, marks_upper_limit=500)
    if !is_connected(petri_net_matrix)
        return Dict{String, Any}(), false
    end

    vertices, edges, arc_transitions, num_transitions, is_bounded = generate_reachability_graph(
        petri_net_matrix,
        place_upper_limit=place_upper_bound,
        max_markings_to_explore=marks_upper_limit,
    )

    if !is_bounded || isempty(vertices) || length(vertices) < marks_lower_limit
        return Dict{String, Any}(), false
    end

    probs, density, markings, rates, success = generate_stochastic_net_task(vertices, edges, arc_transitions, num_transitions)

    if !success || sum(markings) > 1000 || sum(markings) < -1000
        return Dict{String, Any}(), false
    end

    return _create_spn_result_dict(petri_net_matrix, vertices, edges, arc_transitions, rates, probs, density, markings), true
end

function get_spn_info(petri_net_matrix, vertices, edges, arc_transitions, transition_rates)
    if !is_connected(petri_net_matrix) || isempty(vertices)
        return Dict{String, Any}(), false
    end

    probs, density, markings, success = generate_stochastic_net_task_with_rates(vertices, edges, arc_transitions, transition_rates)

    if !success
        return Dict{String, Any}(), false
    end

    return _create_spn_result_dict(petri_net_matrix, vertices, edges, arc_transitions, transition_rates, probs, density, markings), true
end

# Contents of DataTransformation submodule
function _generate_candidate_matrices_core(
    base_petri_matrix,
    enable_delete_edge,
    enable_add_edge,
    enable_add_token,
    enable_delete_token,
)
    num_places, num_cols = size(base_petri_matrix)

    # Optimization: Calculate exact number of matrices to pre-allocate,
    # avoiding dynamically growing `Any` arrays with `push!`.
    num_delete_edge = enable_delete_edge ? count(x -> x == 1, view(base_petri_matrix, :, 1:num_cols-1)) : 0
    num_add_edge = enable_add_edge ? count(x -> x == 0, view(base_petri_matrix, :, 1:num_cols-1)) : 0
    num_add_token = enable_add_token ? num_places : 0

    token_sum = 0
    num_delete_token = 0
    if enable_delete_token
        @inbounds for i in 1:num_places
            token_sum += base_petri_matrix[i, end]
        end
        if token_sum > 1
            @inbounds for i in 1:num_places
                if base_petri_matrix[i, end] >= 1
                    num_delete_token += 1
                end
            end
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

function _generate_candidate_matrices(base_petri_matrix, config)
    candidate_matrices = _generate_candidate_matrices_core(
        base_petri_matrix,
        convert(Bool, get(config, "enable_delete_edge", false)),
        convert(Bool, get(config, "enable_add_edge", false)),
        convert(Bool, get(config, "enable_add_token", false)),
        convert(Bool, get(config, "enable_delete_token", false)),
    )

    num_places, num_cols = size(base_petri_matrix)
    num_transitions = (num_cols - 1) ÷ 2
    if convert(Bool, get(config, "enable_add_place", false)) && num_transitions > 0
        new_place_row = zeros(Int32, 1, num_cols)
        t_idx_to_connect = rand(1:(num_transitions * 2))
        new_place_row[1, t_idx_to_connect] = 1
        modified_matrix = vcat(base_petri_matrix, new_place_row)
        push!(candidate_matrices, modified_matrix)
    end

    return candidate_matrices
end

function _generate_rate_variations_impl(base_variation, p_net::AbstractMatrix{Int32}, num_variations)
    num_trans = (size(p_net, 2) - 1) ÷ 2
    if num_trans == 0
        return Dict{String, Any}[]
    end

    vlist_as_vecs = [v for v in eachrow(base_variation["arr_vlist"]::AbstractMatrix{Int})]

    rate_variations = Vector{Dict{String, Any}}()
    for _ in 1:num_variations
        new_rates = rand(1:10, num_trans)
        s_probs, m_dens, avg_marks, success = generate_stochastic_net_task_with_rates(
            vlist_as_vecs,
            [e for e in eachrow(base_variation["arr_edge"]::AbstractMatrix{Int})],
            base_variation["arr_tranidx"]::Vector{Int},
            new_rates,
        )

        if success
            push!(rate_variations, Dict{String, Any}(
                "petri_net" => p_net,
                "arr_vlist" => base_variation["arr_vlist"],
                "arr_edge" => base_variation["arr_edge"],
                "arr_tranidx" => base_variation["arr_tranidx"]::Vector{Int},
                "spn_labda" => new_rates,
                "spn_steadypro" => s_probs,
                "spn_markdens" => m_dens,
                "spn_allmus" => avg_marks,
                "spn_mu" => sum(avg_marks),
            ))
        end
    end

    return rate_variations
end

function _generate_petri_net_variations_impl(candidate_matrices, place_bound, marks_lower, marks_upper, config)
    results = [
        filter_spn(
            matrix,
            place_upper_bound=place_bound,
            marks_lower_limit=marks_lower,
            marks_upper_limit=marks_upper
        ) for matrix in candidate_matrices
    ]

    structural_variations = [res for (res, success) in results if success]

    all_augmented_data = Vector{Dict{String, Any}}()
    append!(all_augmented_data, structural_variations)

    if convert(Bool, get(config, "enable_rate_variations", false))
        num_rate_variations = convert(Int, get(config, "num_rate_variations_per_structure", 5))
        for base_variation in structural_variations
            p_net = base_variation["petri_net"]::AbstractMatrix{Int32}
            rate_variations = _generate_rate_variations_impl(base_variation, p_net, num_rate_variations)
            append!(all_augmented_data, rate_variations)
        end
    end

    return all_augmented_data
end

function generate_petri_net_variations(petri_matrix, config)
    base_petri_matrix = petri_matrix
    candidate_matrices = _generate_candidate_matrices(base_petri_matrix, config)

    max_candidates = convert(Int, get(config, "max_candidates_per_structure", 50))
    if length(candidate_matrices) > max_candidates
        candidate_matrices = sample(candidate_matrices, max_candidates, replace=false)
    end

    place_bound = convert(Int, get(config, "place_upper_bound", 10))
    marks_lower = convert(Int, get(config, "marks_lower_limit", 4))
    marks_upper = convert(Int, get(config, "marks_upper_limit", 500))

    return _generate_petri_net_variations_impl(candidate_matrices, place_bound, marks_lower, marks_upper, config)
end

function _generate_lambda_variations_impl(petri_dict, petri_net::AbstractMatrix{Int32}, vlist::AbstractMatrix{Int}, edge_list::AbstractMatrix{Int}, tranidx::Vector{Int}, num_lambda_variations)
    num_transitions = (size(petri_net, 2) - 1) ÷ 2
    vlist_as_vecs = [v for v in eachrow(vlist)]
    edge_as_vecs = [e for e in eachrow(edge_list)]

    lambda_variations = Vector{Dict{String, Any}}()
    for _ in 1:num_lambda_variations
        lambda_values = rand(1:10, num_transitions)
        results_dict, success = get_spn_info(
            petri_net,
            vlist_as_vecs,
            edge_as_vecs,
            tranidx,
            lambda_values,
        )
        if success
            push!(lambda_variations, results_dict)
        end
    end

    return lambda_variations
end

function generate_lambda_variations(petri_dict, petri_net::AbstractMatrix{Int32}, num_lambda_variations)
    vlist = petri_dict["arr_vlist"]::AbstractMatrix{Int}
    edge_list = petri_dict["arr_edge"]::AbstractMatrix{Int}
    tranidx = petri_dict["arr_tranidx"]::Vector{Int}
    return _generate_lambda_variations_impl(petri_dict, petri_net, vlist, edge_list, tranidx, num_lambda_variations)
end

# Contents of Utils.jl
"""
Creates a directory if it does not already exist.
"""
function create_directory(path::String)
    if !isdir(path)
        mkpath(path)
        println("Directory created: $path")
    end
end

"""
Loads data from a TOML file.
"""
function load_toml_file(file_path::String)
    return TOML.parsefile(file_path)
end

"""
Saves data to a JSON file.
"""
function save_data_to_json_file(file_path::String, data)
    open(file_path, "w") do f
        JSON3.write(f, data)
    end
end

"""
Loads data from a JSON file.
"""
function load_json_file(file_path::String)
    open(file_path, "r") do f
        return JSON3.read(f)
    end
end

"""
Counts the number of JSON files in a directory.
"""
function count_json_files(directory_path::String)
    files = readdir(directory_path)
    json_files = filter(x -> endswith(x, ".json"), files)
    # The python version sorts the files by a number in the filename.
    # We will replicate that logic here.
    sort!(json_files, by = x -> parse(Int, x[5:end-5]))
    return length(json_files), json_files
end

"""
Samples a specified number of JSON files from a directory.
"""
function sample_json_files_from_directory(num_samples::Int, directory_path::String)
    _, json_files = count_json_files(directory_path)
    if isempty(json_files)
        return Any[]
    end

    num_to_sample = min(num_samples, length(json_files))
    sampled_files = sample(json_files, num_to_sample, replace=false)

    sampled_data = Vector{Any}()
    for file_name in sampled_files
        file_path = joinpath(directory_path, file_name)
        data = load_json_file(file_path)
        # Filter out noisy data
        if -100 <= data[:spn_mu] <= 100
            push!(sampled_data, data)
        else
            # Replace noisy data with a new random sample
            while true
                new_file_name = rand(json_files)
                if new_file_name ∉ sampled_files
                    new_file_path = joinpath(directory_path, new_file_name)
                    new_data = load_json_file(new_file_path)
                    if -100 <= new_data[:spn_mu] <= 100
                        push!(sampled_data, new_data)
                        push!(sampled_files, new_file_name)
                        break
                    end
                end
            end
        end
    end
    return sampled_data
end


"""
Writes a sample to an HDF5 group.
"""
function write_to_hdf5(group, data; compression="gzip", compression_opts=4)
    for (key, value) in data
        try
            if value isa AbstractArray && ndims(value) > 0
                HDF5.create_dataset(group, key, value, ((true, true), (compression, compression_opts)))
            else
                HDF5.create_dataset(group, key, value)
            end
        catch e
            println("Warning: Could not save key '$key' for sample $(HDF5.name(group)). Error: $e")
        end
    end
end

"""
Appends a sample to a JSONL file.
"""
function write_to_jsonl(file_handler, data)
    JSON3.write(file_handler, data)
    write(file_handler, "\n")
end

function run_generation_from_config(config)
    output_format = convert(String, get(config, "output_format", "hdf5"))
    output_dir = joinpath(convert(String, config["output_data_location"]), "data_$(output_format)")
    create_directory(output_dir)
    output_path = joinpath(output_dir, convert(String, config["output_file"]))

    println("Generating $(convert(Int, config["number_of_samples_to_generate"])) initial SPN samples...")
    initial_samples = Vector{Union{Dict{String, Any}, Nothing}}(undef, convert(Int, config["number_of_samples_to_generate"]))
    p = Progress(convert(Int, config["number_of_samples_to_generate"]))
    Threads.@threads for i in 1:convert(Int, config["number_of_samples_to_generate"])
        initial_samples[i] = generate_single_spn(config)
        next!(p)
    end

    valid_samples = filter(x -> !isnothing(x), initial_samples)
    println("Generated $(length(valid_samples)) valid initial samples.")

    all_samples = Vector{Dict{String, Any}}()
    if convert(Bool, get(config, "enable_transformations", false))
        println("Augmenting samples...")
        augmented_lists = Vector{Vector{Dict{String, Any}}}(undef, length(valid_samples))
        p = Progress(length(valid_samples))
        Threads.@threads for i in 1:length(valid_samples)
            augmented_lists[i] = augment_single_spn(valid_samples[i], config)
            next!(p)
        end
        all_samples = vcat(augmented_lists...)
    else
        all_samples = valid_samples
    end

    if isempty(all_samples)
        println("No samples were generated. Skipping file writing and reporting.")
        return
    end

    if output_format == "hdf5"
        h5open(output_path, "w") do hf
            attrs(hf)["generation_config"] = JSON3.write(config)
            dataset_group = create_group(hf, "dataset_samples")
            println("Writing $(length(all_samples)) samples to HDF5...")
            @showprogress for (i, sample) in enumerate(all_samples)
                sample_group = create_group(dataset_group, "sample_$(lpad(i, 7, '0'))")
                write_to_hdf5(sample_group, sample)
            end
            attrs(hf)["total_samples_written"] = length(all_samples)
        end
        println("HDF5 file '$output_path' created successfully.")
    elseif output_format == "jsonl"
        open(output_path, "w") do f
            write(f, JSON3.write(config) * "\n")
            @showprogress for sample in all_samples
                write_to_jsonl(f, sample)
            end
        end
        println("JSONL file '$output_path' created successfully.")
    end
end

end # module SPNGenerator