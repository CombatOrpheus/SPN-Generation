using BenchmarkTools
using Random
include("src/DataTransformation.jl")
include("src/ArrivableGraph.jl")
include("src/SPN.jl")

function filter_spn_opt(petri_net_matrix; place_upper_bound=10, marks_lower_limit=4, marks_upper_limit=500)
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

# Check SPN Generator benchmark
