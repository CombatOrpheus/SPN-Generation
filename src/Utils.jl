module Utils

using TOML
using JSON3
using Random

export create_directory, load_toml_file, save_data_to_json_file, load_json_file, count_json_files, sample_json_files_from_directory

"""
Creates a directory if it does not already exist.
"""
function create_directory(path::String)
    if !isdir(path)
        mkpath(path)
        println("Directory created: $path")
    end
end

"""
Loads data from a TOML file.
"""
function load_toml_file(file_path::String)
    return TOML.parsefile(file_path)
end

"""
Saves data to a JSON file.
"""
function save_data_to_json_file(file_path::String, data)
    open(file_path, "w") do f
        JSON3.write(f, data)
    end
end

"""
Loads data from a JSON file.
"""
function load_json_file(file_path::String)
    open(file_path, "r") do f
        return JSON3.read(f)
    end
end

"""
Counts the number of JSON files in a directory.
"""
function count_json_files(directory_path::String)
    files = readdir(directory_path)
    json_files = filter(x -> endswith(x, ".json"), files)
    # The python version sorts the files by a number in the filename.
    # We will replicate that logic here.
    sort!(json_files, by = x -> parse(Int, x[5:end-5]))
    return length(json_files), json_files
end

"""
Samples a specified number of JSON files from a directory.
"""
function sample_json_files_from_directory(num_samples::Int, directory_path::String)
    _, json_files = count_json_files(directory_path)
    if isempty(json_files)
        return []
    end

    num_to_sample = min(num_samples, length(json_files))
    sampled_files = sample(json_files, num_to_sample, replace=false)

    sampled_data = []
    for file_name in sampled_files
        file_path = joinpath(directory_path, file_name)
        data = load_json_file(file_path)
        # Filter out noisy data
        if -100 <= data[:spn_mu] <= 100
            push!(sampled_data, data)
        else
            # Replace noisy data with a new random sample
            while true
                new_file_name = rand(json_files)
                if new_file_name ∉ sampled_files
                    new_file_path = joinpath(directory_path, new_file_name)
                    new_data = load_json_file(new_file_path)
                    if -100 <= new_data[:spn_mu] <= 100
                        push!(sampled_data, new_data)
                        push!(sampled_files, new_file_name)
                        break
                    end
                end
            end
        end
    end
    return sampled_data
end


"""
Writes a sample to an HDF5 group.
"""
function write_to_hdf5(group, data; compression="gzip", compression_opts=4)
    for (key, value) in data
        try
            if value isa AbstractArray && ndims(value) > 0
                HDF5.create_dataset(group, key, value, ((true, true), (compression, compression_opts)))
            else
                HDF5.create_dataset(group, key, value)
            end
        catch e
            println("Warning: Could not save key '$key' for sample $(HDF5.name(group)). Error: $e")
        end
    end
end

"""
Appends a sample to a JSONL file.
"""
function write_to_jsonl(file_handler, data)
    JSON3.write(file_handler, data)
    write(file_handler, "\n")
end

end # module Utils