# Contents of DataTransformation submodule
function _generate_candidate_matrices_core(
    base_petri_matrix,
    enable_delete_edge,
    enable_add_edge,
    enable_add_token,
    enable_delete_token,
)
    num_places, num_cols = size(base_petri_matrix)

    # Optimization: Calculate exact number of matrices to pre-allocate,
    # avoiding dynamically growing `Any` arrays with `push!`.
    num_delete_edge = enable_delete_edge ? count(x -> x == 1, view(base_petri_matrix, :, 1:num_cols-1)) : 0
    num_add_edge = enable_add_edge ? count(x -> x == 0, view(base_petri_matrix, :, 1:num_cols-1)) : 0
    num_add_token = enable_add_token ? num_places : 0

    token_sum = 0
    num_delete_token = 0
    if enable_delete_token
        @inbounds for i in 1:num_places
            token_sum += base_petri_matrix[i, end]
        end
        if token_sum > 1
            @inbounds for i in 1:num_places
                if base_petri_matrix[i, end] >= 1
                    num_delete_token += 1
                end
            end
        end
    end

    total_matrices = num_delete_edge + num_add_edge + num_add_token + num_delete_token

    # Typed array avoids Any allocations
    candidate_matrices = Vector{typeof(base_petri_matrix)}(undef, total_matrices)
    idx = 1

    # Optimization: Loop column-major without allocating views or using findall
    if enable_delete_edge
        @inbounds for j in 1:num_cols-1
            for i in 1:num_places
                if base_petri_matrix[i, j] == 1
                    modified_matrix = copy(base_petri_matrix)
                    modified_matrix[i, j] = 0
                    candidate_matrices[idx] = modified_matrix
                    idx += 1
                end
            end
        end
    end

    if enable_add_edge
        @inbounds for j in 1:num_cols-1
            for i in 1:num_places
                if base_petri_matrix[i, j] == 0
                    modified_matrix = copy(base_petri_matrix)
                    modified_matrix[i, j] = 1
                    candidate_matrices[idx] = modified_matrix
                    idx += 1
                end
            end
        end
    end

    if enable_add_token
        @inbounds for i in 1:num_places
            modified_matrix = copy(base_petri_matrix)
            modified_matrix[i, end] += 1
            candidate_matrices[idx] = modified_matrix
            idx += 1
        end
    end

    if enable_delete_token && token_sum > 1
        @inbounds for i in 1:num_places
            if base_petri_matrix[i, end] >= 1
                modified_matrix = copy(base_petri_matrix)
                modified_matrix[i, end] -= 1
                candidate_matrices[idx] = modified_matrix
                idx += 1
            end
        end
    end

    return candidate_matrices
end

function _generate_candidate_matrices(base_petri_matrix, config)
    candidate_matrices = _generate_candidate_matrices_core(
        base_petri_matrix,
        convert(Bool, get(config, "enable_delete_edge", false)),
        convert(Bool, get(config, "enable_add_edge", false)),
        convert(Bool, get(config, "enable_add_token", false)),
        convert(Bool, get(config, "enable_delete_token", false)),
    )

    num_places, num_cols = size(base_petri_matrix)
    num_transitions = (num_cols - 1) ÷ 2
    if convert(Bool, get(config, "enable_add_place", false)) && num_transitions > 0
        new_place_row = zeros(Int32, 1, num_cols)
        t_idx_to_connect = rand(1:(num_transitions * 2))
        new_place_row[1, t_idx_to_connect] = 1
        modified_matrix = vcat(base_petri_matrix, new_place_row)
        push!(candidate_matrices, modified_matrix)
    end

    return candidate_matrices
end

function _generate_rate_variations_impl(base_variation, p_net::AbstractMatrix{Int32}, num_variations)
    num_trans = (size(p_net, 2) - 1) ÷ 2
    if num_trans == 0
        return Dict{String, Any}[]
    end

    vlist_as_vecs = [v for v in eachrow(base_variation["arr_vlist"]::AbstractMatrix{Int})]

    rate_variations = Vector{Dict{String, Any}}()
    for _ in 1:num_variations
        new_rates = rand(1:10, num_trans)
        s_probs, m_dens, avg_marks, success = generate_stochastic_net_task_with_rates(
            vlist_as_vecs,
            [e for e in eachrow(base_variation["arr_edge"]::AbstractMatrix{Int})],
            base_variation["arr_tranidx"]::Vector{Int},
            new_rates,
        )

        if success
            push!(rate_variations, Dict{String, Any}(
                "petri_net" => p_net,
                "arr_vlist" => base_variation["arr_vlist"],
                "arr_edge" => base_variation["arr_edge"],
                "arr_tranidx" => base_variation["arr_tranidx"]::Vector{Int},
                "spn_labda" => new_rates,
                "spn_steadypro" => s_probs,
                "spn_markdens" => m_dens,
                "spn_allmus" => avg_marks,
                "spn_mu" => sum(avg_marks),
            ))
        end
    end

    return rate_variations
end

function _generate_petri_net_variations_impl(candidate_matrices, place_bound, marks_lower, marks_upper, config)
    results = [
        filter_spn(
            matrix,
            place_upper_bound=place_bound,
            marks_lower_limit=marks_lower,
            marks_upper_limit=marks_upper
        ) for matrix in candidate_matrices
    ]

    structural_variations = [res for (res, success) in results if success]

    all_augmented_data = Vector{Dict{String, Any}}()
    append!(all_augmented_data, structural_variations)

    if convert(Bool, get(config, "enable_rate_variations", false))
        num_rate_variations = convert(Int, get(config, "num_rate_variations_per_structure", 5))
        for base_variation in structural_variations
            p_net = base_variation["petri_net"]::AbstractMatrix{Int32}
            rate_variations = _generate_rate_variations_impl(base_variation, p_net, num_rate_variations)
            append!(all_augmented_data, rate_variations)
        end
    end

    return all_augmented_data
end

function generate_petri_net_variations(petri_matrix, config)
    base_petri_matrix = petri_matrix
    candidate_matrices = _generate_candidate_matrices(base_petri_matrix, config)

    max_candidates = convert(Int, get(config, "max_candidates_per_structure", 50))
    if length(candidate_matrices) > max_candidates
        # ⚡ Bolt Optimization: Replace undefined `sample` with an in-place partial Fisher-Yates
        # shuffle. This avoids allocating a completely new array (as `shuffle(array)[1:N]` would)
        # and directly modifies the array in place to save memory allocations.
        len = length(candidate_matrices)
        @inbounds for i in 1:max_candidates
            j = rand(i:len)
            tmp = candidate_matrices[i]
            candidate_matrices[i] = candidate_matrices[j]
            candidate_matrices[j] = tmp
        end
        resize!(candidate_matrices, max_candidates)
    end

    place_bound = convert(Int, get(config, "place_upper_bound", 10))
    marks_lower = convert(Int, get(config, "marks_lower_limit", 4))
    marks_upper = convert(Int, get(config, "marks_upper_limit", 500))

    return _generate_petri_net_variations_impl(candidate_matrices, place_bound, marks_lower, marks_upper, config)
end

function _generate_lambda_variations_impl(petri_dict, petri_net::AbstractMatrix{Int32}, vlist::AbstractMatrix{Int}, edge_list::AbstractMatrix{Int}, tranidx::Vector{Int}, num_lambda_variations)
    num_transitions = (size(petri_net, 2) - 1) ÷ 2
    vlist_as_vecs = [v for v in eachrow(vlist)]
    edge_as_vecs = [e for e in eachrow(edge_list)]

    lambda_variations = Vector{Dict{String, Any}}()
    for _ in 1:num_lambda_variations
        lambda_values = rand(1:10, num_transitions)
        results_dict, success = get_spn_info(
            petri_net,
            vlist_as_vecs,
            edge_as_vecs,
            tranidx,
            lambda_values,
        )
        if success
            push!(lambda_variations, results_dict)
        end
    end

    return lambda_variations
end

function generate_lambda_variations(petri_dict, petri_net::AbstractMatrix{Int32}, num_lambda_variations)
    vlist = petri_dict["arr_vlist"]::AbstractMatrix{Int}
    edge_list = petri_dict["arr_edge"]::AbstractMatrix{Int}
    tranidx = petri_dict["arr_tranidx"]::Vector{Int}
    return _generate_lambda_variations_impl(petri_dict, petri_net, vlist, edge_list, tranidx, num_lambda_variations)
end
