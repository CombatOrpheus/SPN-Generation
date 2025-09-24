using Pkg
Pkg.instantiate()

using Printf

# Define paths and commands
const JULIA_SCRIPT = "scripts/SPNGenerate.jl"
const CONFIG_FILE_JL = "config/DataConfig/benchmark_config.toml"
const JULIA_OUTPUT_DIR = "benchmarks/data/julia_output"
const TIME_FILE = "benchmarks/python_time.txt"

# Function to run a command and time it
function time_command(cmd; dir=".")
    return @elapsed run(Cmd(cmd, dir=dir))
end

function run_julia_benchmark()
    # 1. Set up directories
    mkpath(JULIA_OUTPUT_DIR)

    # 2. Benchmark Julia script
    @info "Benchmarking Julia script..."
    julia_cmd = `julia --project=. --threads auto $JULIA_SCRIPT --config $CONFIG_FILE_JL --output_data_location $JULIA_OUTPUT_DIR`
    julia_time = time_command(julia_cmd)

    python_time_str = read(TIME_FILE, String)
    python_time = parse(Float64, python_time_str)

    @printf "Julia execution time: %.2f seconds\n" julia_time

    # 3. Compare results
    @printf "\n--- Benchmark Results ---\n"
    @printf "Python: %.2f s\n" python_time
    @printf "Julia:  %.2f s\n" julia_time
    if julia_time < python_time
        @printf "Julia is %.2fx faster than Python.\n" python_time / julia_time
    else
        @printf "Python is %.2fx faster than Julia.\n" julia_time / python_time
    end

    # 4. Clean up
    rm("benchmarks/data", recursive=true, force=true)
    rm(TIME_FILE)
end

run_julia_benchmark()