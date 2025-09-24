using ..DataGenerate
using ..Utils
using ArgParse
using ProgressMeter
using HDF5
using JSON3
using Random

function get_grid_index(value, grid_boundaries)
    for (i, boundary) in enumerate(grid_boundaries)
        if value < boundary
            return i
        end
    end
    return length(grid_boundaries) + 1
end

function _initialize_grid(grid_dir, accumulate_data, config)
    config_path = joinpath(grid_dir, "config.json")
    if isfile(config_path) && accumulate_data
        grid_config = Utils.load_json_file(config_path)
    else
        row_p = get(config, "places_grid_boundaries", [5 + 2 * (i + 1) for i in 0:4])
        col_m = get(config, "markings_grid_boundaries", [4 + 4 * (i + 1) for i in 0:9])
        grid_config = Dict(
            "row_p" => row_p,
            "col_m" => col_m,
            "json_count" => zeros(Int, length(row_p) + 1, length(col_m) + 1),
        )
    end

    for i in 1:(length(grid_config["row_p"]) + 1)
        for j in 1:(length(grid_config["col_m"]) + 1)
            Utils.create_directory(joinpath(grid_dir, "p$(i)", "m$(j)"))
        end
    end

    return grid_config
end

function partition_data_into_grid(grid_dir, accumulate_data, raw_data_path, config)
    grid_config = _initialize_grid(grid_dir, accumulate_data, config)
    row_p = grid_config["row_p"]
    col_m = grid_config["col_m"]
    dir_counts = grid_config["json_count"]

    all_data = Utils.load_json_file(raw_data_path)

    for data in values(all_data)
        p_idx = get_grid_index(size(data["petri_net"], 1), row_p)
        m_idx = get_grid_index(length(data["arr_vlist"]), col_m)
        dir_counts[p_idx, m_idx] += 1

        save_path = joinpath(
            grid_dir,
            "p$(p_idx)",
            "m$(m_idx)",
            "data$(dir_counts[p_idx, m_idx]).json",
        )
        Utils.save_data_to_json_file(save_path, data)
    end

    grid_config["json_count"] = dir_counts
    Utils.save_data_to_json_file(joinpath(grid_dir, "config.json"), grid_config)
end

function sample_and_transform_data(config)
    grid_data_loc_template = joinpath(config["temporary_grid_location"], "p%s", "m%s")
    all_data = []

    num_place_bins = length(get(config, "places_grid_boundaries", [5 + 2 * (i + 1) for i in 0:4])) + 1
    num_marking_bins = length(get(config, "markings_grid_boundaries", [4 + 4 * (i + 1) for i in 0:9])) + 1

    @showprogress for i in 1:num_place_bins
        for j in 1:num_marking_bins
            dir_path = Printf.format(Printf.Format(grid_data_loc_template), i, j)
            if isdir(dir_path)
                sampled_list = Utils.sample_json_files_from_directory(
                    config["samples_per_grid"], dir_path
                )
                append!(all_data, sampled_list)
            end
        end
    end

    transformed_data = []
    @showprogress for data in all_data
        new_data = DataGenerate.DataTransformation.generate_lambda_variations(data, config["lambda_variations_per_sample"])
        append!(transformed_data, new_data)
    end

    return transformed_data
end

function package_dataset(config, data)
    save_dir = config["output_grid_location"]
    output_format = get(config, "output_format", "hdf5")
    output_file = get(config, "output_file", "grid_dataset.$(output_format)")
    output_path = joinpath(save_dir, output_file)

    Utils.create_directory(save_dir)

    if output_format == "hdf5"
        h5open(output_path, "w") do hf
            attrs(hf)["generation_config"] = JSON3.write(config)
            dataset_group = create_group(hf, "dataset_samples")

            println("Writing $(length(data)) samples to HDF5...")
            @showprogress for (i, sample) in enumerate(data)
                sample_group = create_group(dataset_group, "sample_$(lpad(i, 7, '0'))")
                Utils.write_to_hdf5(sample_group, sample)
            end
            attrs(hf)["total_samples_written"] = length(data)
        end
        println("HDF5 file '$output_path' created successfully.")
    elseif output_format == "jsonl"
        open(output_path, "w") do f
            write(f, JSON3.write(config) * "\n")
            println("Writing $(length(data)) samples to JSONL...")
            @showprogress for sample in data
                Utils.write_to_jsonl(f, sample)
            end
        end
        println("JSONL file '$output_path' created successfully.")
    end
end

function main()
    s = ArgParseSettings(description="Process raw data to generate a grid-based dataset for GNN training.")
    @add_arg_table! s begin
        "--config"
            help = "Path to config TOML file."
            default = "config/DataConfig/PartitionGrid.toml"
    end
    args = parse_args(s)
    config = Utils.load_toml_file(args["config"])

    partition_data_into_grid(
        config["temporary_grid_location"],
        config["accumulation_data"],
        config["raw_data_location"],
        config,
    )

    processed_data = sample_and_transform_data(config)

    package_dataset(config, processed_data)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end