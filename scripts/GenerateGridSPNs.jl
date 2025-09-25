push!(LOAD_PATH, "src")

using SPNGenerator
using ArgParse
using ProgressMeter
using HDF5
using JSON3
using Random
using TOML

function get_grid_index(value, grid_boundaries)
    for (i, boundary) in enumerate(grid_boundaries)
        if value < boundary
            return i
        end
    end
    return length(grid_boundaries) + 1
end


function sample_and_transform_data(config, initial_samples)
    # Partition the generated samples into the grid
    grid = Dict()
    row_p = config["places_grid_boundaries"]
    col_m = config["markings_grid_boundaries"]

    for sample in initial_samples
        p_idx = get_grid_index(size(sample["petri_net"], 1), row_p)
        m_idx = get_grid_index(length(sample["arr_vlist"]), col_m)

        if !haskey(grid, (p_idx, m_idx))
            grid[(p_idx, m_idx)] = []
        end
        push!(grid[(p_idx, m_idx)], sample)
    end

    # Sample from each grid cell and then apply transformations
    all_data = []
    for (key, samples_in_cell) in grid
        num_to_sample = min(length(samples_in_cell), config["samples_per_grid"])
        sampled_list = rand(samples_in_cell, num_to_sample)
        append!(all_data, sampled_list)
    end

    transformed_data = Vector{Any}(undef, length(all_data))
    p = Progress(length(all_data), "Transforming samples: ")
    Threads.@threads for i in 1:length(all_data)
        new_data = SPNGenerator.generate_lambda_variations(all_data[i], config["lambda_variations_per_sample"])
        transformed_data[i] = new_data
        next!(p)
    end
    transformed_data = vcat(filter(x -> !isnothing(x), transformed_data)...)

    return transformed_data
end

function package_dataset(config, data)
    save_dir = config["output_grid_location"]
    output_format = get(config, "output_format", "hdf5")
    output_file = get(config, "output_file", "grid_dataset.$(output_format)")
    output_path = joinpath(save_dir, output_file)

    SPNGenerator.create_directory(save_dir)

    if output_format == "hdf5"
        h5open(output_path, "w") do hf
            attrs(hf)["generation_config"] = JSON3.write(config)
            dataset_group = create_group(hf, "dataset_samples")

            println("Writing $(length(data)) samples to HDF5...")
            @showprogress for (i, sample) in enumerate(data)
                sample_group = create_group(dataset_group, "sample_$(lpad(i, 7, '0'))")
                SPNGenerator.write_to_hdf5(sample_group, sample)
            end
            attrs(hf)["total_samples_written"] = length(data)
        end
        println("HDF5 file '$output_path' created successfully.")
    elseif output_format == "jsonl"
        open(output_path, "w") do f
            write(f, JSON3.write(config) * "\n")
            println("Writing $(length(data)) samples to JSONL...")
            @showprogress for sample in data
                SPNGenerator.write_to_jsonl(f, sample)
            end
        end
        println("JSONL file '$output_path' created successfully.")
    end
end

function main()
    s = ArgParseSettings(description="Generate a grid-based SPN dataset.")
    @add_arg_table! s begin
        "--config"
            help = "Path to config TOML file."
            required = true
    end
    args = parse_args(s)
    config = SPNGenerator.load_toml_file(args["config"])

    # 1. Generate initial SPN samples
    println("Generating $(config["number_of_samples_to_generate"]) initial SPN samples...")
    initial_samples = Vector{Any}(undef, config["number_of_samples_to_generate"])
    p = Progress(config["number_of_samples_to_generate"], "Generating initial samples: ")
    Threads.@threads for i in 1:config["number_of_samples_to_generate"]
        initial_samples[i] = SPNGenerator.generate_single_spn(config)
        next!(p)
    end
    valid_samples = filter(x -> !isnothing(x), initial_samples)
    println("Generated $(length(valid_samples)) valid initial samples.")

    # 2. Sample from grid and apply transformations
    processed_data = sample_and_transform_data(config, valid_samples)

    # 3. Package the final dataset
    package_dataset(config, processed_data)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end