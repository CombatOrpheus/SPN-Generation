# This script provides an interactive command-line interface to set up and run
# a series of SPN (Stochastic Petri Net) generation tasks in parallel.

# Import necessary packages
using TOML
using ProgressMeter
using Dates

# Add the src directory to the load path to be able to load the SPNGenerator module
push!(LOAD_PATH, "src")
using SPNGenerator

const TEMP_CONFIG_DIR = "temp_configs"
const COMPLETED_CONFIG_DIR = "completed_configs"

# --- Configuration Loading ---

function get_default_configs()
    spn_generate_defaults = Dict(
        "number_of_samples_to_generate" => 100,
        "number_of_parallel_jobs" => 4,
        "minimum_number_of_places" => 5,
        "maximum_number_of_places" => 10,
        "minimum_number_of_transitions" => 5,
        "maximum_number_of_transitions" => 10,
        "place_upper_bound" => 10,
        "marks_lower_limit" => 1,
        "marks_upper_limit" => 10,
        "enable_pruning" => false,
        "enable_token_addition" => false,
        "enable_transformations" => false,
        "maximum_transformations_per_sample" => 3,
    )
    partition_grid_defaults = Dict(
        "accumulation_data" => false,
        "samples_per_grid" => 100,
        "lambda_variations_per_sample" => 1,
        "output_format" => "jsonl",
        "output_file" => "grid_dataset",
    )
    return spn_generate_defaults, partition_grid_defaults
end

# --- User Input Functions ---

function get_user_input(prompt::String, default, type_cast::Type=String; help_text::String="")
    if !isempty(help_text)
        println("  ", help_text)
    end

    if type_cast == Bool
        print("$prompt (y/n) [default: $(default ? 'y' : 'n')]: ")
        input = lowercase(strip(readline()))
        if isempty(input)
            return default
        end
        return input == "y"
    end

    print("$prompt [default: $default]: ")
    input = strip(readline())
    if isempty(input)
        return default
    end

    while true
        try
            return parse(type_cast, input)
        catch e
            if e isa ArgumentError
                print("Invalid input. Please enter a value of type $type_cast: ")
                input = strip(readline())
                if isempty(input)
                    return default
                end
            else
                rethrow(e)
            end
        end
    end
end

function get_user_input_for_list(prompt::String, default::String; help_text::String="")
    if !isempty(help_text)
        println("  ", help_text)
    end

    while true
        print("$prompt [default: $default]: ")
        user_input_str = strip(readline())
        if isempty(user_input_str)
            user_input_str = default
        end

        try
            # Meta.parse can parse Julia expressions. We evaluate it in a sandbox module.
            expr = Meta.parse(user_input_str)
            result = @eval(Main, $expr)

            if !(result isa AbstractVector)
                println("Invalid input. Input must evaluate to a vector (e.g., [1, 2, 3] or 1:5).")
                continue
            end

            # Check if all elements are integers
            if all(x -> isa(x, Integer), result)
                return collect(Int, result) # Return as a Vector{Int}
            else
                println("Invalid input. All elements in the list must be integers.")
            end
        catch e
            println("Error evaluating input: $e")
            println("Please enter a valid Julia expression (e.g., [1, 2, 3] or 1:5).")
        end
    end
end


# --- Configuration Setup ---

function get_generation_mode()
    println("\n--- Select Generation Mode ---")
    while true
        mode = get_user_input("Select generation mode (random/grid)", "grid", String)
        if mode in ["random", "grid"]
            return mode
        end
        println("Invalid mode. Please choose 'random' or 'grid'.")
    end
end

function get_common_data_folder()
    println("\n--- Common Data Folder ---")
    return get_user_input(
        "Enter the common data folder path",
        "data",
        String,
        help_text="This folder will be used for all inputs, outputs, and temporary files.",
    )
end

function get_spn_generate_config(defaults, common_data_folder, generation_mode)
    println("\n--- Configuring SPNGenerate.jl ---")
    config = Dict{String, Any}()
    config["output_data_location"] = joinpath(common_data_folder, "raw")
    println("  Output data location is set to: ", config["output_data_location"])

    config["output_file"] = get_user_input(
        "Output filename",
        "spn_dataset.jsonl",
        String,
        help_text="Name of the output file. Must be .jsonl to be used by the next script.",
    )
    config["output_format"] = "jsonl"
    println("  Output format is set to 'jsonl' to be compatible with the next script.")

    if generation_mode == "random"
        config["dataset_sizes"] = get_user_input_for_list(
            "List of dataset sizes (e.g., [1000, 5000] or 1000:1000:5000)",
            "[1000, 5000]",
            help_text="A Julia vector or range for dataset sizes.",
        )
    else # grid mode
        config["number_of_samples_to_generate"] = get_user_input(
            "Number of samples to generate",
            defaults["number_of_samples_to_generate"],
            Int,
            help_text="The number of SPN samples to generate for the grid.",
        )
    end

    config["number_of_parallel_jobs"] = get_user_input(
        "Number of parallel jobs for SPNGenerate",
        defaults["number_of_parallel_jobs"],
        Int,
        help_text="The number of threads to use for data generation within SPNGenerate.jl.",
    )
    config["minimum_number_of_places"] = get_user_input("Minimum number of places", defaults["minimum_number_of_places"], Int)
    config["maximum_number_of_places"] = get_user_input("Maximum number of places", defaults["maximum_number_of_places"], Int)
    config["minimum_number_of_transitions"] = get_user_input("Minimum number of transitions", defaults["minimum_number_of_transitions"], Int)
    config["maximum_number_of_transitions"] = get_user_input("Maximum number of transitions", defaults["maximum_number_of_transitions"], Int)
    config["place_upper_bound"] = get_user_input("Place upper bound", defaults["place_upper_bound"], Int)
    config["marks_lower_limit"] = get_user_input("Markings lower limit", defaults["marks_lower_limit"], Int)
    config["marks_upper_limit"] = get_user_input("Markings upper limit", defaults["marks_upper_limit"], Int)
    config["enable_pruning"] = get_user_input("Enable pruning", defaults["enable_pruning"], Bool)
    config["enable_token_addition"] = get_user_input("Enable token addition", defaults["enable_token_addition"], Bool)
    config["enable_transformations"] = get_user_input("Enable transformations", defaults["enable_transformations"], Bool)
    config["maximum_transformations_per_sample"] = get_user_input("Maximum transformations per sample", defaults["maximum_transformations_per_sample"], Int)

    return config
end

function get_partition_grid_config(defaults, spn_config, common_data_folder)
    println("\n--- Configuring ObtainGridDS.jl ---")
    config = Dict{String, Any}()

    spn_output_dir = joinpath(spn_config["output_data_location"], "data_$(spn_config["output_format"])")
    config["raw_data_location"] = joinpath(spn_output_dir, spn_config["output_file"])
    println("  Raw data location is set to: ", config["raw_data_location"])

    config["temporary_grid_location"] = joinpath(common_data_folder, "temp_grid")
    println("  Temporary grid location is set to: ", config["temporary_grid_location"])

    config["output_grid_location"] = joinpath(common_data_folder, "grid")
    println("  Output grid location is set to: ", config["output_grid_location"])

    config["accumulation_data"] = get_user_input("Accumulate data", defaults["accumulation_data"], Bool)

    # Automatically calculate grid boundaries
    min_places = spn_config["minimum_number_of_places"]
    max_places = spn_config["maximum_number_of_places"]
    num_place_partitions = get_user_input("Number of partitions for places", 5, Int, help_text="How many segments to divide the 'places' range into.")
    config["places_grid_boundaries"] = calculate_boundaries(min_places, max_places, num_place_partitions)
    println("  Calculated places boundaries: ", config["places_grid_boundaries"])

    min_markings = spn_config["marks_lower_limit"]
    max_markings = spn_config["marks_upper_limit"]
    num_marking_partitions = get_user_input("Number of partitions for markings", 10, Int, help_text="How many segments to divide the 'markings' range into.")
    config["markings_grid_boundaries"] = calculate_boundaries(min_markings, max_markings, num_marking_partitions)
    println("  Calculated markings boundaries: ", config["markings_grid_boundaries"])

    config["samples_per_grid"] = get_user_input("Samples per grid", defaults["samples_per_grid"], Int)
    config["lambda_variations_per_sample"] = get_user_input("Lambda variations per sample", defaults["lambda_variations_per_sample"], Int)
    config["output_format"] = get_user_input("Output format (hdf5 or jsonl)", defaults["output_format"], String)
    config["output_file"] = get_user_input("Output filename", defaults["output_file"], String)

    return config
end

function calculate_boundaries(min_val, max_val, num_partitions)
    if num_partitions <= 1
        return [min_val, max_val]
    end
    # Generate `num_partitions + 1` points, and ensure they are integers
    return round.(Int, range(start=min_val, stop=max_val, length=num_partitions + 1))
end

function configure_random_scenarios(spn_defaults, common_data_folder, generation_mode, scenario_count_offset)
    base_spn_config = get_spn_generate_config(spn_defaults, common_data_folder, generation_mode)
    dataset_sizes = pop!(base_spn_config, "dataset_sizes")

    println("\n--- Creating $(length(dataset_sizes)) scenarios based on dataset sizes ---")
    for (i, size) in enumerate(dataset_sizes)
        scenario_name = "scenario_$(scenario_count_offset + i - 1)"
        scenario_dir = joinpath(TEMP_CONFIG_DIR, scenario_name)
        mkpath(scenario_dir)

        spn_config = deepcopy(base_spn_config)
        spn_config["number_of_samples_to_generate"] = size

        original_filename = get(spn_config, "output_file", "spn_dataset.jsonl")
        name, ext = splitext(original_filename)
        spn_config["output_file"] = "$(name)_$(size)_samples$(ext)"

        open(joinpath(scenario_dir, "SPNGenerate.toml"), "w") do f
            TOML.print(f, spn_config)
        end
        println("Configuration for scenario $scenario_name (size: $size) saved in $scenario_dir")
    end
    return length(dataset_sizes)
end

function configure_grid_scenario(spn_defaults, grid_defaults, common_data_folder, generation_mode, scenario_name)
    println("\n--- Configuring Scenario $scenario_name ($generation_mode mode) ---")
    scenario_dir = joinpath(TEMP_CONFIG_DIR, scenario_name)
    mkpath(scenario_dir)

    spn_config = get_spn_generate_config(spn_defaults, common_data_folder, generation_mode)
    open(joinpath(scenario_dir, "SPNGenerate.toml"), "w") do f
        TOML.print(f, spn_config)
    end

    grid_config = get_partition_grid_config(grid_defaults, spn_config, common_data_folder)
    temp_grid_folder = grid_config["temporary_grid_location"]
    open(joinpath(scenario_dir, "PartitionGrid.toml"), "w") do f
        TOML.print(f, grid_config)
    end

    println("\nConfiguration for $scenario_name saved in $scenario_dir")
    return temp_grid_folder
end

# --- Scenario Execution ---

function run_scenario(scenario_dir, scenario_name, generation_mode)
    println("\n--- Running Scenario: $scenario_name ---")

    # Create a log file for the scenario's output
    log_dir = joinpath("logs", scenario_name)
    mkpath(log_dir)
    log_file = joinpath(log_dir, "output.log")

    # Run SPNGenerate.jl
    spn_config_path = joinpath(scenario_dir, "SPNGenerate.toml")
    println("Running SPNGenerate.jl for $scenario_name... Log: $log_file")

    try
        # Redirect output to the log file
        open(log_file, "w") do log
            cmd = `julia scripts/SPNGenerate.jl --config $spn_config_path`
            process = run(pipeline(cmd, stdout=log, stderr=log))
            if process.exitcode != 0
                println("Error running SPNGenerate.jl for $scenario_name. See log for details.")
                return
            end
        end
        println("SPNGenerate.jl completed for $scenario_name.")
    catch e
        println("An error occurred while running SPNGenerate.jl for $scenario_name: $e")
        return
    end

    # Move SPN config file to completed_configs
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    new_spn_config_name = "$(scenario_name)-$(timestamp)-SPNGenerate.toml"
    mv(spn_config_path, joinpath(COMPLETED_CONFIG_DIR, new_spn_config_name))

    if generation_mode == "grid"
        grid_config_path = joinpath(scenario_dir, "PartitionGrid.toml")
        println("Running ObtainGridDS.jl for $scenario_name... Log: $log_file")

        try
            # Append output to the same log file
            open(log_file, "a") do log
                cmd = `julia scripts/ObtainGridDS.jl --config $grid_config_path`
                process = run(pipeline(cmd, stdout=log, stderr=log))
                if process.exitcode != 0
                    println("Error running ObtainGridDS.jl for $scenario_name. See log for details.")
                    return
                end
            end
            println("ObtainGridDS.jl completed for $scenario_name.")
        catch e
            println("An error occurred while running ObtainGridDS.jl for $scenario_name: $e")
            return
        end

        new_grid_config_name = "$(scenario_name)-$(timestamp)-PartitionGrid.toml"
        mv(grid_config_path, joinpath(COMPLETED_CONFIG_DIR, new_grid_config_name))
    end

    println("--- Scenario $scenario_name completed successfully! ---")
end


# Main function to run the interactive generator
function main()
    println("Welcome to the interactive SPN generation setup for Julia.")
    println("This script will guide you through setting up and running multiple scenarios in parallel.")

    common_data_folder = get_common_data_folder()
    temp_grid_folder = ""

    if isdir(TEMP_CONFIG_DIR)
        rm(TEMP_CONFIG_DIR, recursive=true)
    end
    mkpath(TEMP_CONFIG_DIR)

    if !isdir(COMPLETED_CONFIG_DIR)
        mkpath(COMPLETED_CONFIG_DIR)
    end

    spn_defaults, grid_defaults = get_default_configs()
    scenario_count = 1
    generation_modes = Dict{String, String}()

    while true
        add_scenario_input = get_user_input("\nAdd a new scenario?", "y", String)
        if lowercase(add_scenario_input) == "n"
            if scenario_count == 1
                println("No scenarios configured. Exiting.")
                return
            end
            break
        end

        generation_mode = get_generation_mode()

        if generation_mode == "random"
            num_created = configure_random_scenarios(spn_defaults, common_data_folder, generation_mode, scenario_count)
            # Store the mode for each created scenario
            for i in 1:num_created
                generation_modes["scenario_$(scenario_count + i - 1)"] = generation_mode
            end
            scenario_count += num_created
        else # grid mode
            scenario_name = "scenario_$(scenario_count)"
            temp_grid_folder = configure_grid_scenario(
                spn_defaults,
                grid_defaults,
                common_data_folder,
                generation_mode,
                scenario_name,
            )
            generation_modes[scenario_name] = generation_mode
            scenario_count += 1
        end
    end

    println("\nAll scenarios configured.")

    # --- Execution Phase ---
    println("\n--- Starting Execution Phase ---")
    scenarios_to_run = readdir(TEMP_CONFIG_DIR)

    # Ask user for the number of parallel scenarios
    max_parallel = get_user_input(
        "Enter the maximum number of scenarios to run in parallel",
        length(scenarios_to_run), # Default to the total number of scenarios
        Int,
        help_text="This controls how many `julia` processes are spawned at once."
    )

    p = Progress(length(scenarios_to_run), desc="Running Scenarios: ")

    @sync begin
        for scenario_name in scenarios_to_run
            scenario_dir = joinpath(TEMP_CONFIG_DIR, scenario_name)
            if !isdir(scenario_dir)
                continue
            end

            # Limit the number of concurrent tasks
            while length(Base.Workqueues.wq) >= max_parallel
                sleep(1)
            end

            @async begin
                run_scenario(scenario_dir, scenario_name, generation_modes[scenario_name])
                next!(p)
            end
        end
    end

    println("\nAll scenarios have been processed.")
    # Clean up temp directories
    rm(TEMP_CONFIG_DIR, recursive=true)
    if temp_grid_folder != "" && isdir(temp_grid_folder)
        println("Cleaning up temporary grid folder: $temp_grid_folder")
        rm(temp_grid_folder, recursive=true)
    end
end


# Entry point of the script
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end