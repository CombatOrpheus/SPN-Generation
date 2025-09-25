push!(LOAD_PATH, "src")

using SPNGenerator
using ArgParse
using TOML

function main()
    s = ArgParseSettings(description="Generate a random SPN dataset.")
    @add_arg_table! s begin
        "--config"
            help = "Path to config TOML file."
            required = true
    end
    args = parse_args(s)
    config = SPNGenerator.load_toml_file(args["config"])

    # The core logic is in the SPNGenerator module, so we just call it
    SPNGenerator.run_generation_from_config(config)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end