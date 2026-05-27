using BenchmarkTools
using Random
include("src/DataGenerate.jl")

# test replacing count with a simple nested loop on base_petri_matrix
# it looks like allocating view + running count is taking significant allocations
# wait, view doesn't allocate much, it creates a SubArray structure on stack.
# But wait, count with a closure might allocate if it's not well inferred or if view allocations escape.
# The original code:
# num_delete_edge = enable_delete_edge ? count(x -> x == 1, view(base_petri_matrix, :, 1:num_cols-1)) : 0
# The optimized one didn't allocate less because the main allocation is candidate_matrices.
# Let's check candidate_matrices allocations.
