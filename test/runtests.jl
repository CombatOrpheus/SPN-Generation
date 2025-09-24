using Pkg
Pkg.instantiate()

using Test
using JSON3
using Random
using TOML
using SPNBenchmarks

# Define paths and commands
const PYTHON_SCRIPT = "SPNGenerate.py" # Relative to python/ dir
const JULIA_SCRIPT = "scripts/SPNGenerate.jl"
const CONFIG_FILE_PATH = "config/DataConfig/test_config.toml"
const PYTHON_OUTPUT_DIR_REL = "../test/data/python_output" # Relative to python/ dir
const JULIA_OUTPUT_DIR = "test/data/julia_output"
const PYTHON_OUTPUT_FILE = "test/data/python_output/data_jsonl/test_data"
const JULIA_OUTPUT_FILE = "test/data/julia_output/data_jsonl/test_data"

# Function to run a command and check for success
function run_command(cmd; dir=".")
    full_cmd = Cmd(cmd, dir=dir)
    try
        run(full_cmd)
        return true
    catch e
        println("Error running command: $full_cmd")
        println(e)
        return false
    end
end

# Function to compare the structure of two JSON objects
function compare_json_objects(obj1, obj2)
    @test Set(keys(obj1)) == Set(keys(obj2))
    for key in keys(obj1)
        if obj1[key] isa AbstractDict && haskey(obj2, key) && obj2[key] isa AbstractDict
            compare_json_objects(obj1[key], obj2[key])
        elseif obj1[key] isa AbstractArray && haskey(obj2, key) && obj2[key] isa AbstractArray
            @test obj1[key] isa AbstractArray
            @test obj2[key] isa AbstractArray
        elseif haskey(obj2, key)
            @test typeof(obj1[key]) == typeof(obj2[key]) ||
                  (obj1[key] isa Number && obj2[key] isa Number)
        end
    end
end


@testset "SPN-Benchmarks.jl" begin
    # 1. Set up directories
    mkpath(JULIA_OUTPUT_DIR) # Python output dir will be created by the script from within python/

    # 2. Run Python script to generate reference data
    @info "Running Python script to generate reference data..."
    # Note: uv must be run from the python directory to find the venv
    python_cmd = `uv run python $PYTHON_SCRIPT --config ../$CONFIG_FILE_PATH --output_data_location $PYTHON_OUTPUT_DIR_REL`
    @test run_command(python_cmd, dir="python")

    # 3. Run Julia script with the correct number of processes
    @info "Running Julia script..."
    config = TOML.parsefile(CONFIG_FILE_PATH)
    num_procs = get(config, "number_of_parallel_jobs", 1)

    # The --project=@. flag ensures the workers inherit the project environment
    julia_cmd = `julia -p $num_procs --project=@. $JULIA_SCRIPT --config $CONFIG_FILE_PATH --output_data_location $JULIA_OUTPUT_DIR`
    @test run_command(julia_cmd)

    # 4. Compare the output files
    @info "Comparing output files..."
    @test isfile(PYTHON_OUTPUT_FILE)
    @test isfile(JULIA_OUTPUT_FILE)

    python_lines = readlines(PYTHON_OUTPUT_FILE)
    julia_lines = readlines(JULIA_OUTPUT_FILE)

    # First line is config
    @test length(python_lines) > 1
    @test length(julia_lines) > 1

    # Compare the structure of each JSON object up to the minimum number of lines
    num_to_compare = min(length(python_lines), length(julia_lines))
    @info "Comparing the structure of the first $(num_to_compare - 1) generated samples..."
    for i in 2:num_to_compare
        py_obj = JSON3.read(python_lines[i])
        jl_obj = JSON3.read(julia_lines[i])
        compare_json_objects(py_obj, jl_obj)
    end

    # 5. Clean up
    rm("test/data", recursive=true, force=true)
end