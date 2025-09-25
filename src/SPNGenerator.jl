module SPNGenerator

using Random
using DataStructures
using SparseArrays
using IterativeSolvers
using LinearAlgebra
using TOML
using JSON3
using Base.Threads

# Exports from DataGenerate and its submodules
export generate_random_petri_net, prune_petri_net, add_tokens_randomly
export generate_reachability_graph
export filter_spn, get_spn_info
export generate_petri_net_variations, generate_lambda_variations

# Exports from Utils
export create_directory, load_toml_file, save_data_to_json_file, load_json_file, count_json_files, sample_json_files_from_directory
export write_to_hdf5, write_to_jsonl

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
    sub_graph = [first_place, first_transition]
    return petri_matrix, remaining_nodes, sub_graph
end

"""
Connects the remaining nodes to the sub-graph.
"""
function _connect_remaining_nodes(petri_matrix, remaining_nodes, sub_graph, num_places, num_transitions)
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
    petri_matrix, remaining_nodes, sub_graph = _initialize_petri_net(num_places, num_transitions)
    petri_matrix = _connect_remaining_nodes(petri_matrix, remaining_nodes, sub_graph, num_places, num_transitions)

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

function _initialize_bfs(initial_marking)
    marking_index_counter = 1 # Julia is 1-indexed
    visited_markings_list = [initial_marking]
    explored_markings_dict = Dict(initial_marking => marking_index_counter)
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

    reachability_edges = []
    edge_transition_indices = []
    is_bounded = true

    while !isempty(processing_queue)
        current_marking_index = popfirst!(processing_queue)

        enabled_next_markings, enabled_transition_indices, stop_exploration = _process_marking(
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
function _compute_state_equation_numba(num_vertices, edges, arc_transitions, lambda_values)
    state_matrix = spzeros(Float64, num_vertices + 1, num_vertices)
    for i in 1:length(edges)
        edge = edges[i]
        trans_idx = arc_transitions[i]
        src_idx, dest_idx = edge[1], edge[2]
        rate = lambda_values[trans_idx]
        state_matrix[src_idx, src_idx] -= rate
        state_matrix[dest_idx, src_idx] += rate
    end
    state_matrix[num_vertices + 1, :] .= 1.0
    return state_matrix
end

function compute_state_equation(vertices, edges, arc_transitions, lambda_values)
    num_vertices = length(vertices)
    state_matrix = _compute_state_equation_numba(num_vertices, edges, arc_transitions, lambda_values)
    target_vector = zeros(Float64, num_vertices + 1)
    target_vector[end] = 1.0
    return state_matrix, target_vector
end

function compute_average_markings(vertices::Matrix{Int}, steady_state_probs::Vector{Float64})
    avg_tokens_per_place = sum(vertices .* steady_state_probs, dims=1)

    unique_token_values = sort(unique(vertices))
    num_places = size(vertices, 2)
    marking_density_matrix = zeros(Float64, num_places, length(unique_token_values))

    for (token_idx, token_val) in enumerate(unique_token_values)
        states_with_token = vertices .== token_val
        marking_density_matrix[:, token_idx] = sum(states_with_token .* steady_state_probs, dims=1)'
    end

    return marking_density_matrix, vec(avg_tokens_per_place)
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

function _run_spn_task(vertices, edges, arc_transitions, transition_rates)
    if isempty(vertices)
        return nothing, nothing, nothing, false
    end

    vertices_np = Matrix(hcat(vertices...)')
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

    if any(sum(petri_net_matrix[:, 1:(2*num_transitions)], dims=2) .== 0)
        return false
    end

    pre_sum = sum(petri_net_matrix[:, 1:num_transitions], dims=1)
    post_sum = sum(petri_net_matrix[:, (num_transitions + 1):(2*num_transitions)], dims=1)
    if any(pre_sum + post_sum .== 0)
        return false
    end

    return true
end

function _create_spn_result_dict(petri_net_matrix, vertices, edges, arc_transitions, firing_rates, steady_state_probs, marking_densities, average_markings)
    return Dict(
        "petri_net" => petri_net_matrix,
        "arr_vlist" => hcat(vertices...)',
        "arr_edge" => isempty(edges) ? Matrix{Int}(undef, 0, 2) : reduce(vcat, permutedims.(edges)),
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
        return Dict(), false
    end

    vertices, edges, arc_transitions, num_transitions, is_bounded = generate_reachability_graph(
        petri_net_matrix,
        place_upper_limit=place_upper_bound,
        max_markings_to_explore=marks_upper_limit,
    )

    if !is_bounded || isempty(vertices) || length(vertices) < marks_lower_limit
        return Dict(), false
    end

    probs, density, markings, rates, success = generate_stochastic_net_task(vertices, edges, arc_transitions, num_transitions)

    if !success || sum(markings) > 1000 || sum(markings) < -1000
        return Dict(), false
    end

    return _create_spn_result_dict(petri_net_matrix, vertices, edges, arc_transitions, rates, probs, density, markings), true
end

function get_spn_info(petri_net_matrix, vertices, edges, arc_transitions, transition_rates)
    if !is_connected(petri_net_matrix) || isempty(vertices)
        return Dict(), false
    end

    probs, density, markings, success = generate_stochastic_net_task_with_rates(vertices, edges, arc_transitions, transition_rates)

    if !success
        return Dict(), false
    end

    return _create_spn_result_dict(petri_net_matrix, vertices, edges, arc_transitions, transition_rates, probs, density, markings), true
end

# Contents of DataTransformation submodule
function _generate_candidate_matrices_numba(
    base_petri_matrix,
    enable_delete_edge,
    enable_add_edge,
    enable_add_token,
    enable_delete_token,
)
    candidate_matrices = []
    num_places, num_cols = size(base_petri_matrix)

    if enable_delete_edge
        for idx in findall(isone, base_petri_matrix[:, 1:end-1])
            modified_matrix = copy(base_petri_matrix)
            modified_matrix[idx] = 0
            push!(candidate_matrices, modified_matrix)
        end
    end

    if enable_add_edge
        for idx in findall(iszero, base_petri_matrix[:, 1:end-1])
            modified_matrix = copy(base_petri_matrix)
            modified_matrix[idx] = 1
            push!(candidate_matrices, modified_matrix)
        end
    end

    if enable_add_token
        for r in 1:num_places
            modified_matrix = copy(base_petri_matrix)
            modified_matrix[r, end] += 1
            push!(candidate_matrices, modified_matrix)
        end
    end

    if enable_delete_token && sum(base_petri_matrix[:, end]) > 1
        for r in findall(x -> x >= 1, base_petri_matrix[:, end])
            modified_matrix = copy(base_petri_matrix)
            modified_matrix[r, end] -= 1
            push!(candidate_matrices, modified_matrix)
        end
    end

    return candidate_matrices
end

function _generate_candidate_matrices(base_petri_matrix, config)
    candidate_matrices = _generate_candidate_matrices_numba(
        base_petri_matrix,
        get(config, "enable_delete_edge", false),
        get(config, "enable_add_edge", false),
        get(config, "enable_add_token", false),
        get(config, "enable_delete_token", false),
    )

    num_places, num_cols = size(base_petri_matrix)
    num_transitions = (num_cols - 1) ÷ 2
    if get(config, "enable_add_place", false) && num_transitions > 0
        new_place_row = zeros(Int32, 1, num_cols)
        t_idx_to_connect = rand(1:(num_transitions * 2))
        new_place_row[1, t_idx_to_connect] = 1
        modified_matrix = vcat(base_petri_matrix, new_place_row)
        push!(candidate_matrices, modified_matrix)
    end

    return candidate_matrices
end

function _generate_rate_variations(base_variation, num_variations)
    p_net = base_variation["petri_net"]
    num_trans = (size(p_net, 2) - 1) ÷ 2
    if num_trans == 0
        return []
    end

    vlist_as_vecs = [v for v in eachrow(base_variation["arr_vlist"])]

    rate_variations = []
    for _ in 1:num_variations
        new_rates = rand(1:10, num_trans)
        s_probs, m_dens, avg_marks, success = generate_stochastic_net_task_with_rates(
            vlist_as_vecs,
            [e for e in eachrow(base_variation["arr_edge"])],
            base_variation["arr_tranidx"],
            new_rates,
        )

        if success
            push!(rate_variations, Dict(
                "petri_net" => p_net,
                "arr_vlist" => base_variation["arr_vlist"],
                "arr_edge" => base_variation["arr_edge"],
                "arr_tranidx" => base_variation["arr_tranidx"],
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

function generate_petri_net_variations(petri_matrix, config)
    base_petri_matrix = petri_matrix
    candidate_matrices = _generate_candidate_matrices(base_petri_matrix, config)

    max_candidates = get(config, "max_candidates_per_structure", 50)
    if length(candidate_matrices) > max_candidates
        candidate_matrices = sample(candidate_matrices, max_candidates, replace=false)
    end

    place_bound = get(config, "place_upper_bound", 10)
    marks_lower = get(config, "marks_lower_limit", 4)
    marks_upper = get(config, "marks_upper_limit", 500)

    results = [
        filter_spn(
            matrix,
            place_upper_bound=place_bound,
            marks_lower_limit=marks_lower,
            marks_upper_limit=marks_upper
        ) for matrix in candidate_matrices
    ]

    structural_variations = [res for (res, success) in results if success]

    all_augmented_data = []
    append!(all_augmented_data, structural_variations)

    if get(config, "enable_rate_variations", false)
        num_rate_variations = get(config, "num_rate_variations_per_structure", 5)
        for base_variation in structural_variations
            rate_variations = _generate_rate_variations(base_variation, num_rate_variations)
            append!(all_augmented_data, rate_variations)
        end
    end

    return all_augmented_data
end

function generate_lambda_variations(petri_dict, num_lambda_variations)
    petri_net = petri_dict["petri_net"]
    num_transitions = (size(petri_net, 2) - 1) ÷ 2
    vlist_as_vecs = [v for v in eachrow(petri_dict["arr_vlist"])]

    lambda_variations = []
    for _ in 1:num_lambda_variations
        lambda_values = rand(1:10, num_transitions)
        results_dict, success = get_spn_info(
            petri_net,
            vlist_as_vecs,
            [e for e in eachrow(petri_dict["arr_edge"])],
            petri_dict["arr_tranidx"],
            lambda_values,
        )
        if success
            push!(lambda_variations, results_dict)
        end
    end

    return lambda_variations
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
        return []
    end

    num_to_sample = min(num_samples, length(json_files))
    sampled_files = sample(json_files, num_to_sample, replace=false)

    sampled_data = []
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

end # module SPNGenerator