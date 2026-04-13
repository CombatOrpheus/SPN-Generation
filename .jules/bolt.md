
## 2024-04-11 - Vectorized vs Loops in Hot Paths
**Learning:** In Julia, relying on vectorized broadcasts like `all(A .>= B, dims=1)` inside hot loops (such as state-space exploration like `get_enabled_transitions`) allocates intermediate boolean vectors and views which drastically degrades performance.
**Action:** Replace these operations with explicit `@inbounds` nested loops and early return/break mechanisms to avoid memory allocations and achieve ~10-15x speedups in those functions.

## 2026-04-12 - Array Allocation via Vectorization in State Summaries
**Learning:** In Julia, vectorized equality checks like `vertices .== token_val` inside a loop over unique values create many intermediate arrays that heavily penalize performance during data aggregations.
**Action:** Replaced these allocations with a single pass using explicit `@inbounds` nested loops, accumulating counts/values into a pre-allocated matrix. This single-pass O(N) strategy drastically reduced allocations and provided a ~3x speedup.

## 2026-04-12 - Array Layout Matters
**Learning:** In Julia, multidimensional arrays are stored in column-major order. When writing explicit nested loops to replace vectorized code, iterating with the column index (`j`) in the outer loop and the row index (`i`) in the inner loop avoids memory jumps (improving cache locality).
**Action:** Always traverse multidimensional arrays column-by-column rather than row-by-row in hot loops.

## 2024-04-13 - Avoiding allocations in Graph Connectivity Checks
**Learning:** In Julia, vectorized array reductions like `any(sum(matrix, dims=2) .== 0)` allocate intermediate vectors for the sums and boolean conditions. This becomes a major overhead (allocations) in graph algorithms repeatedly checking structural properties like `is_connected`.
**Action:** Replace these vectorized aggregations with explicit `@inbounds` nested loops and early returns. This completely eliminated allocations and resulted in an ~18x speedup.
