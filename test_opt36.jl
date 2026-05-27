using BenchmarkTools
using Random
include("src/DataTransformation.jl")
include("src/ArrivableGraph.jl")
include("src/SPN.jl")
include("src/DataGenerate.jl")

function _connect_remaining_nodes_opt(petri_matrix, remaining_nodes, sub_places, sub_transitions, num_places, num_transitions)
    places_count = 1
    transitions_count = 1
    @inbounds for node in remaining_nodes
        place, transition = if node <= num_places
            (node, sub_transitions[rand(1:transitions_count)])
        else
            (sub_places[rand(1:places_count)], node)
        end

        if rand() <= 0.5
            petri_matrix[place, transition - num_places] = 1
        else
            petri_matrix[place, transition - num_places + num_transitions] = 1
        end

        if node <= num_places
            places_count += 1
            push!(sub_places, node)
        else
            transitions_count += 1
            push!(sub_transitions, node)
        end
    end
    return petri_matrix
end

function _connect_remaining_nodes_opt2(petri_matrix, remaining_nodes, sub_places, sub_transitions, num_places, num_transitions)
    # Using the pre-allocated sub_places and sub_transitions and mutating an index instead of push!
    places_count = 1
    transitions_count = 1
    @inbounds for node in remaining_nodes
        place, transition = if node <= num_places
            (node, sub_transitions[rand(1:transitions_count)])
        else
            (sub_places[rand(1:places_count)], node)
        end

        if rand() <= 0.5
            petri_matrix[place, transition - num_places] = 1
        else
            petri_matrix[place, transition - num_places + num_transitions] = 1
        end

        if node <= num_places
            places_count += 1
            sub_places[places_count] = node
        else
            transitions_count += 1
            sub_transitions[transitions_count] = node
        end
    end
    return petri_matrix
end
