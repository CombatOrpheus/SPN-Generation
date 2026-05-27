using BenchmarkTools
using Random
include("src/DataGenerate.jl")

function add_tokens_randomly_opt(petri_matrix)
    num_places = size(petri_matrix, 1)

    @inbounds for i in 1:num_places
        if rand() <= 0.3
            petri_matrix[i, end] += 1
        end
    end
    return petri_matrix
end

pn = generate_random_petri_net(10, 10)
println("Original:")
b = @benchmarkable add_tokens_randomly(copy(pn)) setup=(pn=$pn)
display(run(b))
println("Optimized:")
b2 = @benchmarkable add_tokens_randomly_opt(copy(pn)) setup=(pn=$pn)
display(run(b2))
