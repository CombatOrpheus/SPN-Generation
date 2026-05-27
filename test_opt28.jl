using BenchmarkTools
using Random

function map_filter_opt()
    results = [ (rand(), rand()>0.5) for _ in 1:100 ]
    return [r[1] for r in results if r[2]]
end

function explicit_loop_opt()
    results = Float64[]
    for _ in 1:100
        val = rand()
        success = rand()>0.5
        if success
            push!(results, val)
        end
    end
    return results
end

println("Original")
@btime map_filter_opt()
println("Optimized")
@btime explicit_loop_opt()
