
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
## 2024-04-16 - Pre-allocation and Explicit Loops over `findall`
**Learning:** In Julia, dynamically growing an untyped array (like `[]` with `push!`) and using allocating functions like `findall(isone, matrix)` or `sum` inside candidate generation functions significantly degrades performance. It causes heavy memory allocation and garbage collection.
**Action:** Pre-calculate the exact required size by counting conditions, pre-allocate a typed `Vector{T}(undef, size)`, and use explicit `@inbounds` column-major nested loops to populate the array. This single-pass strategy without intermediate views reduces latency by ~20% and avoids growing `Any` arrays.
## 2024-04-18 - [Eliminating filter allocations in hot graph building loops]
**Learning:** In Julia, dynamically filtering arrays like `filter(x -> x <= num_places, sub_graph)` inside graph building hot loops allocates many intermediate arrays, severely degrading performance.
**Action:** Replace single `sub_graph` collections that need filtering with explicit, separated typed collections (e.g., `sub_places` and `sub_transitions`) that are pushed to dynamically, effectively eliminating redundant view/allocation creation during iteration.
## 2024-04-18 - Replacing `hcat` and `reduce(vcat)` with explicit loops for matrix building
**Learning:** In Julia, operations like `hcat(arrays...)` and `reduce(vcat, permutedims.(arrays))` on lists of small arrays (like vectors of state or edges) create many intermediate allocations and perform poorly. This creates a large performance bottleneck when generating SPNs.
**Action:** Pre-allocate the resulting matrix using `Matrix{T}(undef, rows, cols)` and populate it using explicit `@inbounds` nested loops. Always use abstract types like `AbstractVector{<:AbstractVector{T}}` in the function signatures instead of concrete `Vector{Vector{T}}` to gracefully handle `SubArray`s or views, avoiding `MethodError`s.
## 2024-04-19 - Replace `findall` and array subsetting with pre-allocated explicit loops
**Learning:** In Julia's hot paths (especially matrix parsing/generating), vectorized operations like `sum(matrix, dims=...)`, `findall`, and `shuffle(array)[1:N]` allocate heavy intermediate arrays.
**Action:** Replace reductions and subsetting with pre-allocated arrays (`Vector{Int}(undef, N)`) populated using explicit, nested `@inbounds` loops. For partial shuffling, manual element swapping is significantly faster and uses less memory.
