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

    # Optimization: Replaced dynamically allocating `sort(unique(vertices))` with
    # a single pass token tracking array, eliminating intermediate array allocations
    # and avoiding expensive `unique` checks.
    if isempty(vertices)
        min_val, max_val = 0, 0
    else
        min_val, max_val = extrema(vertices)
    end

    val_range = max_val - min_val + 1
    present = falses(val_range)
    @inbounds for x in vertices
        present[x - min_val + 1] = true
    end

    num_unique = sum(present)

    token_to_idx = zeros(Int, val_range)
    idx = 1
    for i in 1:val_range
        if present[i]
            token_to_idx[i] = idx
            idx += 1
        end
    end

    offset = 1 - min_val

    marking_density_matrix = zeros(Float64, num_places, num_unique)
    avg_tokens_per_place = zeros(Float64, num_places)

    @inbounds for j in 1:num_places
        for i in 1:num_states
            prob = steady_state_probs[i]
            val = vertices[i, j]
            idx = token_to_idx[val + offset]
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
