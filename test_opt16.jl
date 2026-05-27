using BenchmarkTools
using Random
include("src/ArrivableGraph.jl")

function get_enabled_transitions_opt!(pre_condition_matrix, change_matrix, current_marking_vector, enabled_transitions_buffer, new_markings_buffer, place_upper_limit)
    num_places, num_transitions = size(pre_condition_matrix)

    num_enabled = 0
    # pre_condition_matrix is num_places x num_transitions
    # column-major order is fast for j then i. But current_marking_vector[i] is read inside.
    # What if we use an optimized loop?
    # @simd doesn't help because of early break.
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
        return view(new_markings_buffer, :, 1:0), view(enabled_transitions_buffer, 1:0), false
    end

    @inbounds for j in 1:num_enabled
        trans_idx = enabled_transitions_buffer[j]
        for i in 1:num_places
            new_val = current_marking_vector[i] + change_matrix[i, trans_idx]
            if new_val > place_upper_limit
                return view(new_markings_buffer, :, 1:0), view(enabled_transitions_buffer, 1:0), true
            end
            new_markings_buffer[i, j] = new_val
        end
    end

    return view(new_markings_buffer, :, 1:num_enabled), view(enabled_transitions_buffer, 1:num_enabled), false
end
