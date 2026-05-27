using BenchmarkTools
using SparseArrays
using LinearAlgebra
using IterativeSolvers
include("src/SPN.jl")

# From profile, IterativeSolvers log=false is already used, but sparse solves take a lot of time (which is normal)
# Also SPN.jl compute_average_markings uses extrema.
# Let's write an optimized version of _edges_to_matrix and _vertices_to_matrix. Wait, are they allocating heavily?
function _vertices_to_matrix_opt(vertices::AbstractVector{<:AbstractVector{Int}})
    num_states = length(vertices)
    num_places = length(vertices[1])
    mat = Matrix{Int}(undef, num_states, num_places)
    # Cache locality: matrices are column-major, so loop over j then i
    @inbounds for j in 1:num_places
        for i in 1:num_states
            mat[i, j] = vertices[i][j]
        end
    end
    return mat
end
