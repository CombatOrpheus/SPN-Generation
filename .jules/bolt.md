
## 2024-04-11 - Vectorized vs Loops in Hot Paths
**Learning:** In Julia, relying on vectorized broadcasts like `all(A .>= B, dims=1)` inside hot loops (such as state-space exploration like `get_enabled_transitions`) allocates intermediate boolean vectors and views which drastically degrades performance.
**Action:** Replace these operations with explicit `@inbounds` nested loops and early return/break mechanisms to avoid memory allocations and achieve ~10-15x speedups in those functions.

## 2026-04-12 - Array Allocation via Vectorization in State Summaries
**Learning:** In Julia, vectorized equality checks like `vertices .== token_val` inside a loop over unique values create many intermediate arrays that heavily penalize performance during data aggregations.
**Action:** Replaced these allocations with a single pass using explicit `@inbounds` nested loops, accumulating counts/values into a pre-allocated matrix. This single-pass O(N) strategy drastically reduced allocations and provided a ~3x speedup.

## 2026-04-12 - Array Layout Matters
**Learning:** In Julia, multidimensional arrays are stored in column-major order. When writing explicit nested loops to replace vectorized code, iterating with the column index (`j`) in the outer loop and the row index (`i`) in the inner loop avoids memory jumps (improving cache locality).
**Action:** Always traverse multidimensional arrays column-by-column rather than row-by-row in hot loops.
## 2024-04-14 - [Vectorized vs Loop Allocations in Julia `is_connected`]
**Learning:** Using full vector operations like `any(sum(...) .== 0)` inside hot loops in Julia creates intermediate arrays that require heap allocation, severely impacting performance. The overhead from garbage collection and allocation often outweighs the benefits of vectorization for small matrices.
**Action:** When performing boolean checks on matrices (e.g., checking if a node has any connections), write explicit `@inbounds` nested loops with early `break` or `return` conditions instead of fully vectorized evaluations. This short-circuits the check and completely eliminates temporary array allocations.
## 2024-04-15 - [Vectorized vs Loop Allocations in Julia BFS]
**Learning:** Using full vector operations like `any(enabled_next_markings .> place_upper_limit)` inside a BFS hot loop (`_process_marking`) in Julia creates temporary arrays (like the boolean `.>`) that require heap allocation, severely impacting performance across thousands of iterations.
**Action:** When performing boolean limit checks on matrices in hot loops, write explicit `@inbounds` nested loops with early `break` conditions instead of fully vectorized evaluations. This short-circuits the check and completely eliminates temporary array allocations.
