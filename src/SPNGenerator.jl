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
        place_num = rand(config["minimum_number_of_places"]::Int : config["maximum_number_of_places"]::Int)
        trans_num = rand(config["minimum_number_of_transitions"]::Int : config["maximum_number_of_transitions"]::Int)

        petri_matrix = generate_random_petri_net(place_num, trans_num)
        if get(config, "enable_pruning", false)::Bool
            petri_matrix = prune_petri_net(petri_matrix)
        end
        if get(config, "enable_token_addition", false)::Bool
            petri_matrix = add_tokens_randomly(petri_matrix)
        end

        # Convert to sparse if it's large enough
        if place_num >= 500
            petri_matrix = sparse(petri_matrix)
        end

        results, success = filter_spn(
            petri_matrix,
            place_upper_bound=config["place_upper_bound"]::Int,
            marks_lower_limit=config["marks_lower_limit"]::Int,
            marks_upper_limit=config["marks_upper_limit"]::Int,
        )
        if success
            return results
        end
    end
    return nothing
end

function _augment_single_spn_impl(sample, petri_net::AbstractMatrix{Int32}, config)
    augmented_data::Vector{Dict{String, Any}} = generate_lambda_variations(
        sample,
        petri_net,
        convert(Int, config["lambda_variations_per_sample"]),
    )

    if isempty(augmented_data)
        return Dict{String, Any}[]
    end

    max_transforms = convert(Int, get(config, "maximum_transformations_per_sample", length(augmented_data)))
    if length(augmented_data) > max_transforms
        len = length(augmented_data)
        n = convert(Int, max_transforms)
        @inbounds for i in 1:n
            j = rand(i:len)
            tmp = augmented_data[i]
            augmented_data[i] = augmented_data[j]
            augmented_data[j] = tmp
        end
        resize!(augmented_data, n)
        return augmented_data
    end
    return augmented_data
end

function augment_single_spn(sample, config)
    if isnothing(sample) || !haskey(sample, "petri_net")
        return Dict{String, Any}[]
    end

    petri_net = sample["petri_net"]::AbstractMatrix{Int32}
    return _augment_single_spn_impl(sample, petri_net, config)
end

end # module SPNGenerator