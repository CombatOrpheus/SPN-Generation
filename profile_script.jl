using Profile
using SPNGenerator
using Random

Random.seed!(42)

println("Warming up...")
for _ in 1:10000
    SPNGenerator.generate_random_petri_net(10, 10)
end

println("Profiling generation...")
Profile.clear()
@profile for _ in 1:500000
    SPNGenerator.generate_random_petri_net(10, 10)
end

Profile.print(format=:flat, sortedby=:count)
