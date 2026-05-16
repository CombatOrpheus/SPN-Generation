# Julia Performance Guidelines for SPN Simulator

This guide captures important performance optimizations and anti-patterns specific to Julia, focusing on memory allocation, cache locality, and matrix operations.

## Matrix Operations and Cache Locality

*   **Column-Major Order:** Julia stores matrices in column-major order. When iterating over a matrix, the inner loop must iterate over the first index (rows) and the outer loop over the second index (columns) to maximize cache locality.
*   **Contiguous Matrix Updates:** In hot loops where matrices are repeatedly updated, transposing or reshaping buffer matrices (e.g., from `num_transitions x num_places` to `num_places x num_transitions`) can align both reads and writes. Looping with the row index on the inside ensures all accesses are contiguous, which significantly reduces cache misses.
*   **Avoid Vectorized Reductions in Hot Paths:** Operations like `sum(matrix, dims=...)`, `findall`, or boolean subsetting (`shuffle(array)[1:N]`) allocate intermediate arrays. Replace them with explicit `@inbounds` nested loops and pre-allocated output arrays.
*   **In-Place Array Operations:** Avoid chained vectorized operations like `probs[probs .< 0] .= 0` followed by `probs ./ sum(probs)`. These create intermediate arrays. Use a fused `@simd @inbounds` loop to clamp values, accumulate the sum, and scale in-place.
*   **Scalar Accumulation:** Updating an array inside a hot nested loop (`array[j] += val`) incurs overhead. Use a local scalar variable (e.g., `sum_tokens = 0.0`) inside the inner loop and assign it to the array only once after the inner loop finishes.

## Avoiding Intermediate Allocations

*   **Explicit Matrix Building:** Avoid using `hcat(arrays...)` or `reduce(vcat, permutedims.(arrays))` on lists of small arrays, as this creates many allocations. Instead, pre-allocate the final matrix (`Matrix{T}(undef, rows, cols)`) and populate it using explicit nested loops.
*   **Array Splatting:** Never use `vcat(lists...)` for large lists, as it splats elements into a massive tuple, destroying performance and increasing GC pressure. Use `reduce(vcat, lists)` or pre-allocate and `append!`.
*   **Zero-Allocation Sampling:** Avoid `sample(array, n, replace=false)` (which also requires importing `StatsBase`). Use an in-place partial Fisher-Yates shuffle followed by `resize!` or a `view` to eliminate allocations.
*   **Dense Token Hashing:** Calling `sort(unique(vertices))` inside hot loops allocates dictionaries and arrays. For dense integer matrices, use a single-pass boolean tracking array bounded by `extrema()`, creating an $O(1)$ mapping with zero hashing allocations.

## Optimizing Dictionaries and Collections

*   **Array-Based Lookups:** Using a `Dict` to map densely packed integers to indices in a hot loop introduces significant hashing overhead. Use a pre-allocated vector with offset indexing (`array[val + offset]`) for $O(1)$ lookups.
*   **Lazy Array Allocation for Dict Keys:** `Dict{Vector{Int}, Int}` can index keys using views (e.g., `view(Matrix, i, :)`). Only convert the view to a dense array (`Vector{Int}(new_marking_view)`) if the key doesn't exist and needs to be added. This saves extensive allocations during BFS explorations.
*   **Pre-sizing Collections:** In graph algorithms or BFS traversals where a maximum limit is known (e.g., `max_markings_to_explore`), call `sizehint!(collection, limit)` immediately after initialization to prevent costly dynamic array resizing.
*   **Vector vs. Queue in BFS:** `DataStructures.Queue` introduces memory allocation overhead from linked list or block updates. Use a pre-allocated `Vector{Int}` and a `queue_head` counter when the maximum number of nodes is known.

## Algorithm and Signature Design

*   **Fused Bounds Checking:** In state generation loops (like BFS), fuse bounds checking directly into the inner generation loop to enable early termination. Do not generate all valid states first and run a secondary $O(N \times M)$ pass to check limits.
*   **Abstract Array Types:** Function signatures should use abstract types like `AbstractVector{<:AbstractVector{T}}` rather than concrete types like `Vector{Vector{T}}` to gracefully handle `SubArray`s or views, avoiding `MethodError`s.
