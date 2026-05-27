using BenchmarkTools
using Random
include("src/DataTransformation.jl")
include("src/ArrivableGraph.jl")
include("src/SPN.jl")

function is_connected_opt(petri_net_matrix)
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
