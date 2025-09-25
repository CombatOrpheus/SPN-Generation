using ArgParse
using TOML
using ProgressMeter
using Base.Threads

function main()
    s = ArgParseSettings(description="Run multiple SPN generation scenarios concurrently.")
    @add_arg_table! s begin
        "scenarios_dir"
            help = "Path to the directory containing scenario TOML files."
            required = true
        "--max_parallel"
            help = "Maximum number of scenarios to run in parallel."
            arg_type = Int
            default = 4
    end
    args = parse_args(s)

    scenarios_dir = args["scenarios_dir"]
    if !isdir(scenarios_dir)
        println("Error: Directory not found at ", scenarios_dir)
        return
    end

    config_files = filter(f -> endswith(f, ".toml"), readdir(scenarios_dir, join=true))
    if isempty(config_files)
        println("No .toml configuration files found in ", scenarios_dir)
        return
    end

    println("Found $(length(config_files)) scenarios to run.")
    p = Progress(length(config_files), desc="Running Scenarios: ")

    tasks = []
    for config_path in config_files
        task = @spawn begin
            config = TOML.parsefile(config_path)
            task_type = get(config, "task_type", "")

            script_to_run = ""
            if task_type == "random"
                script_to_run = "scripts/GenerateRandomSPNs.jl"
            elseif task_type == "grid"
                script_to_run = "scripts/GenerateGridSPNs.jl"
            else
                println("Warning: Skipping unknown task_type '$task_type' in $config_path")
                return
            end

            scenario_name = first(splitext(basename(config_path)))
            log_dir = joinpath("logs", scenario_name)
            mkpath(log_dir)
            log_file = joinpath(log_dir, "output.log")

            println("Starting $scenario_name ($task_type)... Log: $log_file")

            try
                open(log_file, "w") do log
                    cmd = `julia --project=. $script_to_run --config $config_path`
                    process = run(pipeline(cmd, stdout=log, stderr=log))
                    if process.exitcode == 0
                        println("$scenario_name completed successfully.")
                    else
                        println("Error in $scenario_name. See log for details.")
                    end
                end
            catch e
                println("Failed to run $scenario_name: $e")
            end
            next!(p)
        end
        push!(tasks, task)
    end

    # Wait for all tasks to complete
    for task in tasks
        wait(task)
    end

    println("\nAll scenarios have been processed.")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end