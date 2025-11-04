# Contents of DataTransformation submodule
function generate_candidate_matrices_numba(
    base_petri_matrix,
    enable_delete_edge,
    enable_add_edge,
    enable_add_token,
    enable_delete_token,
)
    candidate_matrices = Vector{Matrix{Int32}}()
    num_places, num_cols = size(base_petri_matrix)

    # Pre-allocate a single matrix for modifications
    modified_matrix = copy(base_petri_matrix)

    if enable_delete_edge
        for idx in findall(isone, base_petri_matrix[:, 1:end-1])
            modified_matrix[idx] = 0
            push!(candidate_matrices, copy(modified_matrix))
            modified_matrix[idx] = 1 # Restore original value
        end
    end

    if enable_add_edge
        for idx in findall(iszero, base_petri_matrix[:, 1:end-1])
            modified_matrix[idx] = 1
            push!(candidate_matrices, copy(modified_matrix))
            modified_matrix[idx] = 0 # Restore original value
        end
    end

    if enable_add_token
        for r in 1:num_places
            modified_matrix[r, end] += 1
            push!(candidate_matrices, copy(modified_matrix))
            modified_matrix[r, end] -= 1 # Restore original value
        end
    end

    if enable_delete_token && sum(base_petri_matrix[:, end]) > 1
        for r in findall(x -> x >= 1, base_petri_matrix[:, end])
            modified_matrix[r, end] -= 1
            push!(candidate_matrices, copy(modified_matrix))
            modified_matrix[r, end] += 1 # Restore original value
        end
    end

    return candidate_matrices
end

function generate_candidate_matrices(base_petri_matrix, config)
    candidate_matrices = generate_candidate_matrices_numba(
        base_petri_matrix,
        get(config, "enable_delete_edge", false),
        get(config, "enable_add_edge", false),
        get(config, "enable_add_token", false),
        get(config, "enable_delete_token", false),
    )

    num_places, num_cols = size(base_petri_matrix)
    num_transitions = (num_cols - 1) ÷ 2
    if get(config, "enable_add_place", false) && num_transitions > 0
        new_place_row = zeros(Int32, 1, num_cols)
        t_idx_to_connect = rand(1:(num_transitions * 2))
        new_place_row[1, t_idx_to_connect] = 1
        modified_matrix = vcat(base_petri_matrix, new_place_row)
        push!(candidate_matrices, modified_matrix)
    end

    return candidate_matrices
end

function generate_rate_variations(base_variation, num_variations)
    p_net = base_variation["petri_net"]
    num_trans = (size(p_net, 2) - 1) ÷ 2
    if num_trans == 0
        return Vector{Dict}()
    end

    vlist_as_vecs = [v for v in eachrow(base_variation["arr_vlist"])]

    rate_variations = Vector{Dict}()
    for _ in 1:num_variations
        new_rates = rand(1:10, num_trans)
        s_probs, m_dens, avg_marks, success = generate_stochastic_net_task_with_rates(
            vlist_as_vecs,
            [e for e in eachrow(base_variation["arr_edge"])],
            base_variation["arr_tranidx"],
            new_rates,
        )

        if success
            push!(rate_variations, Dict(
                "petri_net" => p_net,
                "arr_vlist" => base_variation["arr_vlist"],
                "arr_edge" => base_variation["arr_edge"],
                "arr_tranidx" => base_variation["arr_tranidx"],
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

function generate_petri_net_variations(petri_matrix, config)
    base_petri_matrix = petri_matrix
    candidate_matrices = generate_candidate_matrices(base_petri_matrix, config)

    max_candidates = get(config, "max_candidates_per_structure", 50)
    if length(candidate_matrices) > max_candidates
        candidate_matrices = sample(candidate_matrices, max_candidates, replace=false)
    end

    place_bound = get(config, "place_upper_bound", 10)
    marks_lower = get(config, "marks_lower_limit", 4)
    marks_upper = get(config, "marks_upper_limit", 500)

    results = [
        filter_spn(
            matrix,
            place_upper_bound=place_bound,
            marks_lower_limit=marks_lower,
            marks_upper_limit=marks_upper
        ) for matrix in candidate_matrices
    ]

    structural_variations = [res for (res, success) in results if success]

    all_augmented_data = []
    append!(all_augmented_data, structural_variations)

    if get(config, "enable_rate_variations", false)
        num_rate_variations = get(config, "num_rate_variations_per_structure", 5)
        for base_variation in structural_variations
            rate_variations = generate_rate_variations(base_variation, num_rate_variations)
            append!(all_augmented_data, rate_variations)
        end
    end

    return all_augmented_data
end

function generate_lambda_variations(petri_dict, num_lambda_variations)
    petri_net = petri_dict["petri_net"]
    num_transitions = (size(petri_net, 2) - 1) ÷ 2
    vlist_as_vecs = [v for v in eachrow(petri_dict["arr_vlist"])]

    lambda_variations = []
    for _ in 1:num_lambda_variations
        lambda_values = rand(1:10, num_transitions)
        results_dict, success = get_spn_info(
            petri_net,
            vlist_as_vecs,
            [e for e in eachrow(petri_dict["arr_edge"])],
            petri_dict["arr_tranidx"],
            lambda_values,
        )
        if success
            push!(lambda_variations, results_dict)
        end
    end

    return lambda_variations
end
