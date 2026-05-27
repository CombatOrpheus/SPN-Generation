using BenchmarkTools
using Random
include("src/DataTransformation.jl")

function _generate_petri_net_variations_impl_opt(candidate_matrices, place_bound, marks_lower, marks_upper, config)
    # Original:
    # results = [
    #     filter_spn(
    #         matrix,
    #         place_upper_bound=place_bound,
    #         marks_lower_limit=marks_lower,
    #         marks_upper_limit=marks_upper
    #     ) for matrix in candidate_matrices
    # ]
    # structural_variations = [res for (res, success) in results if success]

    # Optimized: Do it in one pass to avoid intermediate array allocations
    structural_variations = Vector{Dict{String, Any}}()

    for matrix in candidate_matrices
        res, success = filter_spn(
            matrix,
            place_upper_bound=place_bound,
            marks_lower_limit=marks_lower,
            marks_upper_limit=marks_upper
        )
        if success
            push!(structural_variations, res)
        end
    end

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
