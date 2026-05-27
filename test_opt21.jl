using BenchmarkTools
using Random
include("src/DataGenerate.jl")

function prune_petri_net_opt(petri_matrix)
    num_transitions = (size(petri_matrix, 2) - 1) ÷ 2

    # We inline delete_excess_edges but use a local tuple/array for edges.
    # Actually wait, allocating `edge_indices = Vector{Int}(undef, max(num_places, num_cols))`
    # takes a tiny amount of time, but can be avoided if we do it differently.
    # Or we can just use delete_excess_edges since it takes 450ns.
    # We found `add_missing_connections_opt` makes a 2-3x difference (948ns to 1.2us? Wait, Original was 948ns, Optimized was 1.291us. NO, Optimized was SLOWER!)
    # Let me re-check test_opt6.jl:
    # Original median: 948ns.
    # Optimized median: 1.291us.
    # Ah! The optimized `add_missing_connections_opt` was slower because iterating rows then columns on small matrices (10x10) is perfectly fine and fits in L1 cache. The "optimization" of keeping a `falses` array actually adds overhead for small sizes!

end
