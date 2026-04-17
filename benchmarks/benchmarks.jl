using BenchmarkTools
using SPNGenerator
using Random

const SUITE = BenchmarkGroup()

SUITE["generation"] = BenchmarkGroup()
SUITE["generation"]["generate_random_petri_net"] = @benchmarkable SPNGenerator.generate_random_petri_net(10, 10)

# We need a reproducible Petri net that successfully passes the filter to benchmark filtering and augmentation
Random.seed!(42)
pn_valid = nothing
sample_valid = nothing
while true
    global pn_valid, sample_valid
    pn_valid = SPNGenerator.generate_random_petri_net(10, 10)
    res, success = SPNGenerator.filter_spn(pn_valid)
    if success
        sample_valid = res
        break
    end
end

SUITE["filtering"] = BenchmarkGroup()
SUITE["filtering"]["filter_spn"] = @benchmarkable SPNGenerator.filter_spn($pn_valid)

augment_config = Dict{String, Any}(
    "lambda_variations_per_sample" => 3,
    "maximum_transformations_per_sample" => 5
)

SUITE["augmentation"] = BenchmarkGroup()
SUITE["augmentation"]["augment_single_spn"] = @benchmarkable SPNGenerator.augment_single_spn($sample_valid, $augment_config)

whole_program_config = Dict{String, Any}(
    "output_data_location" => "benchmarks_output_tmp",
    "output_file" => "benchmark_data",
    "output_format" => "jsonl",
    "number_of_samples_to_generate" => 10,
    "number_of_parallel_jobs" => 1,
    "minimum_number_of_places" => 5,
    "maximum_number_of_places" => 10,
    "minimum_number_of_transitions" => 5,
    "maximum_number_of_transitions" => 10,
    "place_upper_bound" => 10,
    "marks_lower_limit" => 4,
    "marks_upper_limit" => 500,
    "enable_pruning" => true,
    "enable_token_addition" => true,
    "enable_transformations" => false
)

# Function to clean up after benchmark run
function run_and_clean(config)
    SPNGenerator.run_generation_from_config(config)
    rm(config["output_data_location"], recursive=true, force=true)
end

SUITE["whole_program"] = BenchmarkGroup()
SUITE["whole_program"]["run_generation_from_config"] = @benchmarkable run_and_clean($whole_program_config)

if abspath(PROGRAM_FILE) == @__FILE__
    println("Tuning benchmarks...")
    tune!(SUITE)
    println("Running benchmarks...")
    results = run(SUITE, verbose=true)
    println("Benchmark Results:")
    display(median(results))
end
