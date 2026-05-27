using BenchmarkTools
using Random
include("src/SPN.jl")

function _edges_to_matrix_opt(edges::Union{AbstractVector{<:AbstractVector{Int}}, AbstractVector{<:Tuple{Int, Int}}})
    if isempty(edges)
        return Matrix{Int}(undef, 0, 2)
    end
    n = length(edges)
    mat = Matrix{Int}(undef, n, 2)
    @inbounds for i in 1:n
        mat[i, 1] = edges[i][1]
    end
    @inbounds for i in 1:n
        mat[i, 2] = edges[i][2]
    end
    return mat
end

e = [(rand(1:10), rand(1:10)) for _ in 1:1000]
println("Original:")
b = @benchmarkable _edges_to_matrix(e) setup=(e=$e)
display(run(b))
println("Optimized:")
b2 = @benchmarkable _edges_to_matrix_opt(e) setup=(e=$e)
display(run(b2))
