## 2024-04-14 - [Vectorized vs Loop Allocations in Julia `is_connected`]
**Learning:** Using full vector operations like `any(sum(...) .== 0)` inside hot loops in Julia creates intermediate arrays that require heap allocation, severely impacting performance. The overhead from garbage collection and allocation often outweighs the benefits of vectorization for small matrices.
**Action:** When performing boolean checks on matrices (e.g., checking if a node has any connections), write explicit `@inbounds` nested loops with early `break` or `return` conditions instead of fully vectorized evaluations. This short-circuits the check and completely eliminates temporary array allocations.

## 2024-04-15 - [Vectorized vs Loop Allocations in Julia BFS]
**Learning:** Using full vector operations like `any(enabled_next_markings .> place_upper_limit)` inside a BFS hot loop (`_process_marking`) in Julia creates temporary arrays (like the boolean `.>`) that require heap allocation, severely impacting performance across thousands of iterations.
**Action:** When performing boolean limit checks on matrices in hot loops, write explicit `@inbounds` nested loops with early `break` conditions instead of fully vectorized evaluations. This short-circuits the check and completely eliminates temporary array allocations.

## 2024-04-18 - [Eliminating filter allocations in hot graph building loops]
**Learning:** In Julia, dynamically filtering arrays like `filter(x -> x <= num_places, sub_graph)` inside graph building hot loops allocates many intermediate arrays, severely degrading performance.
**Action:** Replace single `sub_graph` collections that need filtering with explicit, separated typed collections (e.g., `sub_places` and `sub_transitions`) that are pushed to dynamically, effectively eliminating redundant view/allocation creation during iteration.

## 2025-04-25 - Avoid eagerly allocating state arrays when checking BFS dictionaries
**Learning:** In Julia, `Dict{Vector{Int}, Int}` can correctly index and check keys using an `AbstractVector` view like `view(Matrix, i, :)`. Previously, checking if a state has been explored in `_update_graph!` eagerness allocated a full `Vector{Int}` for every potential state transition, even if the state was already visited.
**Action:** Always wrap the potential next state from a pre-allocated matrix as a `view` (`new_marking_view = view(enabled_next_markings, i, :)`), pass the view to `haskey(explored_markings_dict, new_marking_view)`, and ONLY convert it to a dense array (`Vector{Int}(new_marking_view)`) if the state is truly new and must be added to the queue/visited list. This saves `O(transitions * max_explored)` array allocations and significantly reduces garbage collection pressure.

## 2024-05-18 - Replacing Queue with a pre-allocated Vector in BFS
**Learning:** In Julia, dynamically allocating and tracking nodes in a BFS algorithm using `DataStructures.Queue` introduces measurable memory allocation overhead due to internal linked list updates (or deque block allocations).
**Action:** Replace `DataStructures.Queue` with a simple pre-allocated `Vector{Int}` and a `queue_head` index counter. This drastically reduces allocations since we already know the maximum number of nodes we will explore (`max_markings_to_explore`).

## 2024-05-20 - [Fused bounds checking in state generation]
**Learning:** In algorithmic state generation (like BFS graph generation), extracting bounds checking into a secondary pass iterating over all generated states allocates extra execution time.
**Action:** Fuse the bounds check directly into the inner state-generation loop (e.g., inside `get_enabled_transitions!`) to enable early termination. This avoids computing all valid markings first and running a secondary O(N*M) pass to check limits.

## 2024-05-23 - Contiguous Matrix Updates Avoid Cache Misses in BFS Loops
**Learning:** In hot loops where matrices are repeatedly updated (like in `get_enabled_transitions!`), writing to a matrix with standard column-major inner loop iteration over columns is cache aligned for the write, but reading an input matrix via the varying column index forces non-contiguous reads, causing many cache misses. Furthermore, hashing a sliced row view later is slightly slower than hashing a sliced column view.
**Action:** Transpose or reshape the buffer matrix (e.g. from `num_transitions x num_places` to `num_places x num_transitions`). Loop with the row index on the inside to make ALL matrix accesses completely contiguous across both read (`change_matrix`) and write (`new_markings_buffer`) arrays. Slicing column views for downstream hashing is also faster.

## 2024-04-24 - [Massive memory overhead from ConvergenceHistory in IterativeSolvers.jl]
**Learning:** In Julia's `IterativeSolvers.jl`, calling `lsmr` (or similar solvers) with `log=true` forces the internal allocation of a `ConvergenceHistory` object that tracks residuals up to `maxiter`. For large matrices where `maxiter` is set high (e.g., $100 \times N$), this causes a massive memory allocation on every call (e.g., ~1.3MB overhead per solve), severely degrading GC performance.
**Action:** When using `IterativeSolvers.jl` in hot paths, set `log=false` to skip the history allocation and manually compute a residual norm (e.g., `norm(A * x - b) < tol`) after the solve if a convergence check is required.
