using BenchmarkTools
using SparseArrays
using IterativeSolvers

function _compute_state_equation_core_old(num_vertices, edges::AbstractMatrix, arc_transitions, lambda_values)
    num_edges = size(edges, 1)
    I = Vector{Int}(undef, 2 * num_edges + num_vertices)
    J = Vector{Int}(undef, 2 * num_edges + num_vertices)
    V = Vector{Float64}(undef, 2 * num_edges + num_vertices)

    idx = 1
    @inbounds for i in 1:num_edges
        trans_idx = arc_transitions[i]
        src_idx = edges[i, 1]
        dest_idx = edges[i, 2]
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

    @inbounds for j in 1:num_vertices
        I[idx] = num_vertices + 1
        J[idx] = j
        V[idx] = 1.0
        idx += 1
    end

    return sparse(I, J, V, num_vertices + 1, num_vertices)
end

function _compute_state_equation_core_new(num_vertices, edges::AbstractMatrix, arc_transitions, lambda_values)
    num_edges = size(edges, 1)

    # We can use a more compact allocation here.
    I = Vector{Int}(undef, 2 * num_edges + num_vertices)
    J = Vector{Int}(undef, 2 * num_edges + num_vertices)
    V = Vector{Float64}(undef, 2 * num_edges + num_vertices)

    idx = 1
    @inbounds for i in 1:num_edges
        trans_idx = arc_transitions[i]
        src_idx = edges[i, 1]
        dest_idx = edges[i, 2]
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

    @inbounds for j in 1:num_vertices
        I[idx] = num_vertices + 1
        J[idx] = j
        V[idx] = 1.0
        idx += 1
    end

    return sparse(I, J, V, num_vertices + 1, num_vertices)
end

num_vertices = 100
edges = rand(1:num_vertices, 200, 2)
arc_transitions = rand(1:10, 200)
lambda_values = rand(1.0:10.0, 10)

println("Old:")
@btime _compute_state_equation_core_old($num_vertices, $edges, $arc_transitions, $lambda_values)
println("New:")
@btime _compute_state_equation_core_new($num_vertices, $edges, $arc_transitions, $lambda_values)
