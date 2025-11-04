module SPNGenerator

using Random
using DataStructures
using SparseArrays
using HDF5
using IterativeSolvers
using LinearAlgebra
using TOML
using JSON3
using Base.Threads
using ProgressMeter
using Logging

# Include the refactored files
include("DataGenerate.jl")
include("ArrivableGraph.jl")
include("SPN.jl")
include("DataTransformation.jl")
include("Utils.jl")

# Exports from DataGenerate and its submodules
export generate_random_petri_net, prune_petri_net, add_tokens_randomly
export generate_reachability_graph
export filter_spn, get_spn_info, generate_single_spn, augment_single_spn
export generate_petri_net_variations, generate_lambda_variations

# Exports from Utils
export create_directory, load_toml_file, save_data_to_json_file, load_json_file, count_json_files, sample_json_files_from_directory
export write_to_hdf5, write_to_jsonl, run_generation_from_config

function generate_single_spn(config)
    max_attempts = 100
    for _ in 1:max_attempts
        place_num = rand(config["minimum_number_of_places"]:config["maximum_number_of_places"])
        trans_num = rand(config["minimum_number_of_transitions"]:config["maximum_number_of_transitions"])

        petri_matrix = generate_random_petri_net(place_num, trans_num)
        if get(config, "enable_pruning", false)
            petri_matrix = prune_petri_net(petri_matrix)
        end
        if get(config, "enable_token_addition", false)
            petri_matrix = add_tokens_randomly(petri_matrix)
        end

        results, success = filter_spn(
            petri_matrix,
            place_upper_bound=config["place_upper_bound"],
            marks_lower_limit=config["marks_lower_limit"],
            marks_upper_limit=config["marks_upper_limit"],
        )
        if success
            return results
        end
    end
    return nothing
end

function augment_single_spn(sample, config)
    if isnothing(sample) || !haskey(sample, "petri_net")
        return Vector{Dict}()
    end

    petri_net = sample["petri_net"]
    augmented_data = generate_lambda_variations(
        sample,
        config["lambda_variations_per_sample"],
    )

    if isempty(augmented_data)
        return []
    end

    max_transforms = get(config, "maximum_transformations_per_sample", length(augmented_data))
    if length(augmented_data) > max_transforms
        return sample(augmented_data, max_transforms, replace=false)
    end
    return augmented_data
end

end # module SPNGenerator
