# Contents of SPN submodule
function compute_state_equation_numba(num_vertices, edges, arc_transitions, lambda_values)
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
    state_matrix = compute_state_equation_numba(num_vertices, edges, arc_transitions, lambda_values)
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

function run_spn_task(vertices, edges, arc_transitions, transition_rates)
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
    probs, density, markings, success = run_spn_task(vertices, edges, arc_transitions, transition_rates)
    return probs, density, markings, transition_rates, success
end

function generate_stochastic_net_task_with_rates(vertices, edges, arc_transitions, transition_rates)
    return run_spn_task(vertices, edges, arc_transitions, transition_rates)
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

function create_spn_result_dict(petri_net_matrix, vertices, edges, arc_transitions, firing_rates, steady_state_probs, marking_densities, average_markings)
    return Dict{String, Any}(
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

    return create_spn_result_dict(petri_net_matrix, vertices, edges, arc_transitions, rates, probs, density, markings), true
end

function get_spn_info(petri_net_matrix, vertices, edges, arc_transitions, transition_rates)
    if !is_connected(petri_net_matrix) || isempty(vertices)
        return Dict(), false
    end

    probs, density, markings, success = generate_stochastic_net_task_with_rates(vertices, edges, arc_transitions, transition_rates)

    if !success
        return Dict(), false
    end

    return create_spn_result_dict(petri_net_matrix, vertices, edges, arc_transitions, transition_rates, probs, density, markings), true
end
