# SPNGenerator

SPNGenerator is a Julia project for generating, analyzing, and transforming Stochastic Petri Nets (SPNs).

## Architecture

The project's architecture is modular and logically divided into several core components:
- **Data Generation** (`DataGenerate.jl`): Functions to initialize, connect, prune, and add tokens to random Petri nets.
- **Graph Analysis** (`ArrivableGraph.jl`): Identification of enabled transitions, marking calculations, and reachability graph generation using BFS.
- **SPN Calculations** (`SPN.jl`): State equation formulation, average marking computations, and steady-state probability solving for stochastic nets.
- **Data Transformation** (`DataTransformation.jl`): Creation of variations of Petri net structures (e.g., adding/deleting edges or tokens) and lambda (transition rate) variations.
- **Utilities** (`Utils.jl`): Functions to load/save JSON and TOML files, write data to HDF5 and JSONL formats, and manage directories.

## Setup

1. **Install Julia:**
   Julia must be installed to run this project. You can install it using:
   ```bash
   curl -fsSL https://install.julialang.org | sh -s -- -y
   ```

2. **Update Shell Environment:**
   After installing Julia, ensure it is available in your PATH:
   ```bash
   source ~/.bashrc
   ```

3. **Instantiate the Project:**
   Navigate to the project root and instantiate the package environment:
   ```bash
   julia --project -e 'using Pkg; Pkg.instantiate()'
   ```

## Usage

You can generate SPN datasets by running the provided scripts with a configuration TOML file.

**Generate Random SPNs:**
```bash
julia --project scripts/GenerateRandomSPNs.jl --config path/to/your/config.toml
```

**Run Scenarios:**
There is also a script to run scenarios:
```bash
julia --project scripts/run_scenarios.jl
```

## Testing

Tests are run using the standard Julia package manager command. To execute the full test suite, run:
```bash
julia --project -e 'using Pkg; Pkg.test()'
```
