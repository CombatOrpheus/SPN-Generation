using Pkg
Pkg.instantiate()

using Printf

# Define paths and commands
const PYTHON_SCRIPT = "SPNGenerate.py" # Relative to python/ dir
const JULIA_SCRIPT = "scripts/SPNGenerate.jl"
const CONFIG_FILE_PY = "../config/DataConfig/test_config.toml" # Relative to python/ dir
const CONFIG_FILE_JL = "config/DataConfig/test_config.toml"
const PYTHON_OUTPUT_DIR_REL = "../benchmarks/data/python_output" # Relative to python/ dir
const JULIA_OUTPUT_DIR = "benchmarks/data/julia_output"

# Function to run a command and time it
function time_command(cmd; dir=".")
    return @elapsed run(Cmd(cmd, dir=dir))
end

function run_benchmarks()
    # 1. Set up directories
    mkpath(JULIA_OUTPUT_DIR)

    # 2. Benchmark Python script
    @info "Benchmarking Python script..."
    python_cmd = `uv run python $PYTHON_SCRIPT --config $CONFIG_FILE_PY --output_data_location $PYTHON_OUTPUT_DIR_REL`
    python_time = time_command(python_cmd, dir="python")
    @printf "Python execution time: %.2f seconds\n" python_time

    # 3. Benchmark Julia script
    @info "Benchmarking Julia script..."
    julia_cmd = `julia --project=. $JULIA_SCRIPT --config $CONFIG_FILE_JL --output_data_location $JULIA_OUTPUT_DIR`
    julia_time = time_command(julia_cmd)
    @printf "Julia execution time: %.2f seconds\n" julia_time

    # 4. Compare results
    @printf "\n--- Benchmark Results ---\n"
    @printf "Python: %.2f s\n" python_time
    @printf "Julia:  %.2f s\n" julia_time
    if julia_time < python_time
        @printf "Julia is %.2fx faster than Python.\n" python_time / julia_time
    else
        @printf "Python is %.2fx faster than Julia.\n" julia_time / python_time
    end

    # 5. Clean up
    rm("benchmarks/data", recursive=true, force=true)
end

run_benchmarks()