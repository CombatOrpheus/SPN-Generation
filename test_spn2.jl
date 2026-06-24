using BenchmarkTools

function test_sum1(vertices, probs)
    num_states, num_places = size(vertices)
    avg_tokens_per_place = zeros(Float64, num_places)
    for j in 1:num_places
        sum_tokens = 0.0
        for i in 1:num_states
            prob = probs[i]
            val = vertices[i, j]
            sum_tokens += val * prob
        end
        avg_tokens_per_place[j] = sum_tokens
    end
    return avg_tokens_per_place
end

function test_sum2(vertices, probs)
    num_states, num_places = size(vertices)
    avg_tokens_per_place = zeros(Float64, num_places)
    for i in 1:num_states
        prob = probs[i]
        for j in 1:num_places
            val = vertices[i, j]
            avg_tokens_per_place[j] += val * prob
        end
    end
    return avg_tokens_per_place
end

vertices = rand(1:10, 1000, 10)
probs = rand(1000)

println("test_sum1:")
@btime test_sum1($vertices, $probs)

println("test_sum2:")
@btime test_sum2($vertices, $probs)
