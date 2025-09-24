using Pkg
Pkg.instantiate()

using Printf

# Define paths and commands for the Julia script
const JULIA_SCRIPT = "scripts/SPNGenerate.jl"
const CONFIG_FILE_JL = "config/DataConfig/test_config.toml"
const JULIA_OUTPUT_DIR = "benchmarks/data/julia_output"

# Function to run a command and time it
function time_command(cmd; dir=".")
    # The backticks create a Cmd object, which we can time with @elapsed
    return @elapsed run(Cmd(cmd, dir=dir))
end

function run_julia_benchmark()
    # 1. Set up the output directory
    mkpath(JULIA_OUTPUT_DIR)

    # 2. Benchmark the Julia script
    @info "Benchmarking Julia script..."
    julia_cmd = `julia --project=. $JULIA_SCRIPT --config $CONFIG_FILE_JL --output_data_location $JULIA_OUTPUT_DIR`
    julia_time = time_command(julia_cmd)

    @printf "\n--- Julia Benchmark Result ---\n"
    @printf "Execution time: %.2f seconds\n" julia_time

    # 3. Clean up the output directory
    rm("benchmarks/data", recursive=true, force=true)
end

# Run the benchmark
run_julia_benchmark()