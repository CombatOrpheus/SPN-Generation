# Add src to the LOAD_PATH for the main process
push!(LOAD_PATH, abspath("src"))

using Distributed
using ArgParse
using ProgressMeter
using HDF5
using JSON3
using Random
using DataGenerate
using Utils

# Make the functions available on all workers
@everywhere push!(LOAD_PATH, abspath("src"))
@everywhere using DataGenerate, Utils
@everywhere using Random

@everywhere function generate_single_spn(config)
    max_attempts = 100
    for _ in 1:max_attempts
        place_num = rand(config["minimum_number_of_places"]:config["maximum_number_of_places"])
        trans_num = rand(config["minimum_number_of_transitions"]:config["maximum_number_of_transitions"])

        petri_matrix = DataGenerate.generate_random_petri_net(place_num, trans_num)
        if get(config, "enable_pruning", false)
            petri_matrix = DataGenerate.prune_petri_net(petri_matrix)
        end
        if get(config, "enable_token_addition", false)
            petri_matrix = DataGenerate.add_tokens_randomly(petri_matrix)
        end

        results, success = DataGenerate.SPN.filter_spn(
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

@everywhere function augment_single_spn(sample, config)
    if isnothing(sample) || !haskey(sample, "petri_net")
        return []
    end

    petri_net = sample["petri_net"]
    augmented_data = DataGenerate.DataTransformation.generate_petri_net_variations(
        petri_net,
        config,
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

function setup_arg_parser()
    s = ArgParseSettings(description="Generate Stochastic Petri Net (SPN) datasets.")
    @add_arg_table! s begin
        "--config"
            help = "Path to config TOML file."
            default = "config/DataConfig/SPNGenerate.toml"
        "--output_data_location"
            help = "Save directory."
        "--output_file"
            help = "Output filename."
        "--output_format"
            help = "Output file format (hdf5 or jsonl)."
        "--number_of_samples_to_generate"
            help = "Number of samples."
            arg_type = Int
        "--number_of_parallel_jobs"
            help = "Number of parallel jobs."
            arg_type = Int
        "--minimum_number_of_places"
            help = "Min number of places."
            arg_type = Int
        "--maximum_number_of_places"
            help = "Max number of places."
            arg_type = Int
        "--minimum_number_of_transitions"
            help = "Min number of transitions."
            arg_type = Int
        "--maximum_number_of_transitions"
            help = "Max number of transitions."
            arg_type = Int
        "--place_upper_bound"
            help = "Upper bound for places."
            arg_type = Int
        "--marks_lower_limit"
            help = "Lower limit for markings."
            arg_type = Int
        "--marks_upper_limit"
            help = "Upper limit for markings."
            arg_type = Int
        "--enable_pruning"
            action = :store_true
            help = "Enable pruning."
        "--enable_token_addition"
            action = :store_true
            help = "Enable adding tokens."
        "--enable_transformations"
            action = :store_true
            help = "Enable augmentation."
        "--maximum_transformations_per_sample"
            help = "Max number of transformations."
            arg_type = Int
    end
    return s
end

function load_config(args)
    config = isfile(args["config"]) ? Utils.load_toml_file(args["config"]) : Dict()
    for (key, value) in args
        if !isnothing(value)
            config[key] = value
        end
    end
    return config
end

function run_generation_from_config(config)
    output_format = get(config, "output_format", "hdf5")
    output_dir = joinpath(config["output_data_location"], "data_$(output_format)")
    Utils.create_directory(output_dir)
    output_path = joinpath(output_dir, config["output_file"])

    println("Generating $(config["number_of_samples_to_generate"]) initial SPN samples...")
    initial_samples = @showprogress pmap(1:config["number_of_samples_to_generate"]) do _
        generate_single_spn(config)
    end

    valid_samples = filter(x -> !isnothing(x), initial_samples)
    println("Generated $(length(valid_samples)) valid initial samples.")

    all_samples = []
    if get(config, "enable_transformations", false)
        println("Augmenting samples...")
        augmented_lists = @showprogress pmap(valid_samples) do sample
            augment_single_spn(sample, config)
        end
        all_samples = vcat(augmented_lists...)
    else
        all_samples = valid_samples
    end

    if isempty(all_samples)
        println("No samples were generated. Skipping file writing and reporting.")
        return
    end

    if output_format == "hdf5"
        h5open(output_path, "w") do hf
            attrs(hf)["generation_config"] = JSON3.write(config)
            dataset_group = create_group(hf, "dataset_samples")
            println("Writing $(length(all_samples)) samples to HDF5...")
            @showprogress for (i, sample) in enumerate(all_samples)
                sample_group = create_group(dataset_group, "sample_$(lpad(i, 7, '0'))")
                Utils.write_to_hdf5(sample_group, sample)
            end
            attrs(hf)["total_samples_written"] = length(all_samples)
        end
        println("HDF5 file '$output_path' created successfully.")
    elseif output_format == "jsonl"
        open(output_path, "w") do f
            write(f, JSON3.write(config) * "\n")
            @showprogress for sample in all_samples
                Utils.write_to_jsonl(f, sample)
            end
        end
        println("JSONL file '$output_path' created successfully.")
    end
end

function main()
    s = setup_arg_parser()
    args = parse_args(s)
    config = load_config(args)

    num_jobs = get(config, "number_of_parallel_jobs", 1)
    if num_jobs > nworkers()
        addprocs(num_jobs - nworkers())
    end

    # The @everywhere directives at the top of the file will handle loading on all workers.
    run_generation_from_config(config)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end