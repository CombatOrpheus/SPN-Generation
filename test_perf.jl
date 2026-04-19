using BenchmarkTools
using Random
using InteractiveUtils

function add_missing_connections_orig(petri_matrix_in, num_transitions)
    petri_matrix = copy(petri_matrix_in)
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

function add_missing_connections_opt(petri_matrix_in, num_transitions)
    petri_matrix = copy(petri_matrix_in)
    num_places = size(petri_matrix, 1)

    # 1. Ensure each transition has at least one connection
    # Transitions are in columns 1 to 2*num_transitions
    @inbounds for j in 1:(2*num_transitions)
        has_connection = false
        for i in 1:num_places
            if petri_matrix[i, j] != 0
                has_connection = true
                break
            end
        end
        if !has_connection
            random_row = rand(1:num_places)
            petri_matrix[random_row, j] = 1
        end
    end

    # 2. Ensure each place has at least one incoming edge
    # Incoming edges are in columns 1 to num_transitions
    @inbounds for i in 1:num_places
        has_incoming = false
        for j in 1:num_transitions
            if petri_matrix[i, j] != 0
                has_incoming = true
                break
            end
        end
        if !has_incoming
            random_col = rand(1:num_transitions)
            petri_matrix[i, random_col] = 1
        end
    end

    # 3. Ensure each place has at least one outgoing edge
    # Outgoing edges are in columns (num_transitions + 1) to 2*num_transitions
    @inbounds for i in 1:num_places
        has_outgoing = false
        for j in (num_transitions + 1):(2*num_transitions)
            if petri_matrix[i, j] != 0
                has_outgoing = true
                break
            end
        end
        if !has_outgoing
            random_col = rand(1:num_transitions) + num_transitions
            petri_matrix[i, random_col] = 1
        end
    end

    return petri_matrix
end

Random.seed!(42)
num_places = 10
num_transitions = 5
petri_matrix = rand(0:1, num_places, 2 * num_transitions + 1)
# Create some zeros for the missing connections
petri_matrix[:, 2] .= 0
petri_matrix[:, 4] .= 0
petri_matrix[3, 1:num_transitions] .= 0
petri_matrix[4, (num_transitions+1):(2*num_transitions)] .= 0
petri_matrix[:, end] .= 0

println("Original (add_missing_connections):")
@btime add_missing_connections_orig(petri_matrix, num_transitions)

println("Optimized (add_missing_connections):")
@btime add_missing_connections_opt(petri_matrix, num_transitions)
