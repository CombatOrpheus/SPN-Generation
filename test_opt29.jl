using BenchmarkTools
using Random
include("src/DataGenerate.jl")

function _initialize_petri_net_opt(num_places::Int, num_transitions::Int)
    remaining_nodes = collect(1:(num_places + num_transitions))
    petri_matrix = zeros(Int32, num_places, 2 * num_transitions + 1)

    first_place = rand(1:num_places)
    first_transition = rand((num_places + 1):(num_places + num_transitions))

    # In place partial removal without filter!
    idx = 1
    @inbounds for i in 1:(num_places + num_transitions)
        if i != first_place && i != first_transition
            remaining_nodes[idx] = i
            idx += 1
        end
    end
    resize!(remaining_nodes, idx - 1)

    if rand() <= 0.5
        petri_matrix[first_place, first_transition - num_places] = 1
    else
        petri_matrix[first_place, first_transition - num_places + num_transitions] = 1
    end

    shuffle!(remaining_nodes)
    sub_places = Int[first_place]
    sub_transitions = Int[first_transition]

    sizehint!(sub_places, num_places)
    sizehint!(sub_transitions, num_transitions)

    return petri_matrix, remaining_nodes, sub_places, sub_transitions
end

function generate_random_petri_net_opt(num_places::Int, num_transitions::Int)
    petri_matrix, remaining_nodes, sub_places, sub_transitions = _initialize_petri_net_opt(num_places, num_transitions)
    petri_matrix = _connect_remaining_nodes(petri_matrix, remaining_nodes, sub_places, sub_transitions, num_places, num_transitions)

    # Add an initial marking
    random_place = rand(1:num_places)
    petri_matrix[random_place, end] = 1

    return petri_matrix
end

pn = generate_random_petri_net(10, 10)
println("Original:")
b = @benchmarkable generate_random_petri_net(10, 10)
display(run(b))
println("Optimized:")
b2 = @benchmarkable generate_random_petri_net_opt(10, 10)
display(run(b2))
