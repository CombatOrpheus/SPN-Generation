# Scenario Configuration Documentation

This document describes the configuration parameters for the data generation scenarios. Each scenario is defined by a `.toml` file and is executed by the `scripts/run_scenarios.jl` script.

## Common Parameters

These parameters are required for all scenarios.

- `task_type`: (String) The type of generation task to perform. Must be either `"random"` or `"grid"`.
- `output_data_location`: (String) The directory where the output data will be saved. For random scenarios, this will be the root directory for the output. For grid scenarios, this is the root directory for the grid output.
- `output_file`: (String) The name of the output file.
- `output_format`: (String) The format of the output file. Can be `"jsonl"` or `"hdf5"`.

## "random" Scenario Parameters

These parameters are used when `task_type = "random"`.

- `number_of_samples_to_generate`: (Integer) The number of SPN samples to generate.
- `number_of_parallel_jobs`: (Integer) The number of threads to use for data generation.
- `minimum_number_of_places`: (Integer) The minimum number of places in the SPN.
- `maximum_number_of_places`: (Integer) The maximum number of places in the SPN.
- `minimum_number_of_transitions`: (Integer) The minimum number of transitions in the SPN.
- `maximum_number_of_transitions`: (Integer) The maximum number of transitions in the SPN.
- `place_upper_bound`: (Integer) The upper bound for the number of places in the SPN.
- `marks_lower_limit`: (Integer) The lower limit for the number of markings in the SPN.
- `marks_upper_limit`: (Integer) The upper limit for the number of markings in the SPN.
- `enable_pruning`: (Boolean) A boolean indicating whether to prune the SPN.
- `enable_token_addition`: (Boolean) A boolean indicating whether to add tokens to the SPN.
- `enable_transformations`: (Boolean) A boolean indicating whether to apply transformations to the SPNs.
- `maximum_transformations_per_sample`: (Integer) The maximum number of transformations to apply to each SPN.
- `lambda_variations_per_sample`: (Integer) The number of lambda variations to generate for each sample (if `enable_transformations` is true).

## "grid" Scenario Parameters

These parameters are used when `task_type = "grid"`.

### SPN Generation Parameters

These parameters are used for the initial generation of SPNs before they are partitioned into a grid.

- `number_of_samples_to_generate`: (Integer) The number of initial SPN samples to generate.
- `number_of_parallel_jobs`: (Integer) The number of threads to use for data generation.
- `minimum_number_of_places`: (Integer) The minimum number of places in the SPN.
- `maximum_number_of_places`: (Integer) The maximum number of places in the SPN.
- `minimum_number_of_transitions`: (Integer) The minimum number of transitions in the SPN.
- `maximum_number_of_transitions`: (Integer) The maximum number of transitions in the SPN.
- `place_upper_bound`: (Integer) The upper bound for the number of places in the SPN.
- `marks_lower_limit`: (Integer) The lower limit for the number of markings in the SPN.
- `marks_upper_limit`: (Integer) The upper limit for the number of markings in the SPN.
- `enable_pruning`: (Boolean) A boolean indicating whether to prune the SPN.
- `enable_token_addition`: (Boolean) A boolean indicating whether to add tokens to the SPN.
- `enable_transformations`: (Boolean) A boolean indicating whether to apply transformations to the SPNs.

### Grid Partitioning Parameters

- `places_grid_boundaries`: (Array of Integers) The boundaries for the grid rows (number of places).
- `markings_grid_boundaries`: (Array of Integers) The boundaries for the grid columns (number of markings).
- `samples_per_grid`: (Integer) The number of samples to take from each grid cell.
- `lambda_variations_per_sample`: (Integer) The number of lambda variations to generate for each sample.