using BenchmarkTools

function test_haskey_vector(dict, mat, col)
    vec = mat[:, col]
    return haskey(dict, vec)
end

function test_haskey_view(dict, mat, col)
    v = view(mat, :, col)
    return haskey(dict, v)
end

function test_haskey_view_optimized(dict, mat, col)
    v = view(mat, :, col)
    return haskey(dict, v)
end

n = 10
d = Dict{Vector{Int}, Int}()
for i in 1:1000
    v = rand(1:10, n)
    d[v] = i
end

mat = rand(1:10, n, 1000)

println("Vector:")
@btime test_haskey_vector($d, $mat, 1)

println("View:")
@btime test_haskey_view($d, $mat, 1)
