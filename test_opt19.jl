using BenchmarkTools
using Random
include("src/SPN.jl")

function compute_average_markings_opt(vertices::Matrix{Int}, steady_state_probs::Vector{Float64})
    num_states, num_places = size(vertices)

    # Instead of extrema which does 2 passes or allocates (wait, extrema is single pass, but extrema(vertices) allocates if we do it on a matrix?)
    # Wait, extrema on a matrix returns a single min and max over all elements, single pass.
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

    # Optimization: Transpose iteration for marking_density_matrix writes
    # marking_density_matrix is num_places x num_unique (column major!)
    # `marking_density_matrix[j, idx] += prob` accesses varying columns (idx) for a fixed row (j)
    # in the inner loop (i over states). So it's random access across columns.
    # What if marking_density_matrix is num_unique x num_places, and we transpose at the end?
    # Or just iterate over num_places (j) on the outside, which we already do!
    # Wait, the current code is:
    # for j in 1:num_places
    #     sum_tokens = 0.0
    #     for i in 1:num_states
    #         val = vertices[i, j]
    #         idx = token_to_idx[val + offset]
    #         marking_density_matrix[j, idx] += prob
    #     end
    #     avg_tokens_per_place[j] = sum_tokens
    # end
    # `vertices[i, j]` is accessed sequentially since i is inner loop.
    # `marking_density_matrix[j, idx]` is row j, column idx.
    # So we are writing to the SAME row `j` repeatedly, but different columns `idx`.
    # In column major, a row is NOT contiguous. This is a cache miss.

    @inbounds for j in 1:num_places
        sum_tokens = 0.0
        for i in 1:num_states
            prob = steady_state_probs[i]
            val = vertices[i, j]
            idx = token_to_idx[val + offset]
            marking_density_matrix[j, idx] += prob
            sum_tokens += val * prob
        end
        avg_tokens_per_place[j] = sum_tokens
    end

    return marking_density_matrix, avg_tokens_per_place
end

function compute_average_markings_opt2(vertices::Matrix{Int}, steady_state_probs::Vector{Float64})
    num_states, num_places = size(vertices)

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

    # Transposed density matrix to be column major for writes
    marking_density_matrix_t = zeros(Float64, num_unique, num_places)
    avg_tokens_per_place = zeros(Float64, num_places)

    @inbounds for j in 1:num_places
        sum_tokens = 0.0
        for i in 1:num_states
            prob = steady_state_probs[i]
            val = vertices[i, j]
            idx = token_to_idx[val + offset]
            # Write is contiguous in memory since j is fixed and we vary idx?
            # Wait, no. For fixed j (outer loop), idx changes.
            # We are writing to column j, row idx!
            # Since marking_density_matrix_t is num_unique x num_places, column j is contiguous memory!
            marking_density_matrix_t[idx, j] += prob
            sum_tokens += val * prob
        end
        avg_tokens_per_place[j] = sum_tokens
    end

    return collect(marking_density_matrix_t'), avg_tokens_per_place
end

v = rand(0:10, 500, 10)
p = rand(500)
p ./= sum(p)
println("Original:")
b = @benchmarkable compute_average_markings(v, p) setup=(v=$v; p=$p)
display(run(b))
println("Optimized:")
b2 = @benchmarkable compute_average_markings_opt2(v, p) setup=(v=$v; p=$p)
display(run(b2))
