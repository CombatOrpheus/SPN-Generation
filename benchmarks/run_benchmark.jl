using Pkg
Pkg.instantiate()

push!(LOAD_PATH, "src")
using SPNGenerator
using BenchmarkTools
using Printf
using JSON3

const CONFIG_FILE_JL = "config/DataConfig/test_config.toml"

function run_julia_benchmark()
    config = SPNGenerator.load_toml_file(CONFIG_FILE_JL)

    # Create a temporary directory for the output
    output_dir = mktempdir()
    config["output_data_location"] = output_dir
    config["output_file"] = "benchmark_data"

    println("Running benchmark with the following configuration:")
    display(config)

    b = @benchmark begin
        initial_samples = Vector{Any}(undef, $config["number_of_samples_to_generate"])
        Threads.@threads for i in 1:$config["number_of_samples_to_generate"]
            max_attempts = 100
            for _ in 1:max_attempts
                place_num = rand($config["minimum_number_of_places"]:$config["maximum_number_of_places"])
                trans_num = rand($config["minimum_number_of_transitions"]:$config["maximum_number_of_transitions"])

                petri_matrix = SPNGenerator.generate_random_petri_net(place_num, trans_num)
                if get($config, "enable_pruning", false)
                    petri_matrix = SPNGenerator.prune_petri_net(petri_matrix)
                end
                if get($config, "enable_token_addition", false)
                    petri_matrix = SPNGenerator.add_tokens_randomly(petri_matrix)
                end

                results, success = SPNGenerator.filter_spn(
                    petri_matrix,
                    place_upper_bound=$config["place_upper_bound"],
                    marks_lower_limit=$config["marks_lower_limit"],
                    marks_upper_limit=$config["marks_upper_limit"],
                )
                if success
                    initial_samples[i] = results
                    break
                end
            end
        end
    end

    println("\n--- Benchmark Results ---")
    display(b)

    # Clean up the temporary directory
    rm(output_dir, recursive=true)
end

run_julia_benchmark()