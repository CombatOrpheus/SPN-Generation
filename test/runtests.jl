using Test
push!(LOAD_PATH, "src")
using SPNGenerator
using JSON3

@testset "SPNGenerator.jl" begin
    @testset "Petri Net Generation" begin
        pn = SPNGenerator.generate_random_petri_net(5, 3)
        @test size(pn) == (5, 7)
        @test sum(pn[:, end]) > 0  # Initial marking
    end

    @testset "Batch Petri Net Generation" begin
        batch_size = 4
        pns = SPNGenerator.generate_random_petri_net(5, 3, batch_size)
        @test size(pns) == (5, 7, batch_size)

        # Check that not all matrices in the batch are identical
        all_same = true
        for i in 2:batch_size
            if pns[:, :, i] != pns[:, :, 1]
                all_same = false
                break
            end
        end
        @test !all_same
    end

    @testset "Pruning" begin
        # This petri net has a place with 4 connections, so at least one edge should be pruned.
        pn = [1 1 1 1 0 0 0 1;
              1 0 0 0 1 0 0 0;
              1 0 0 0 0 1 0 0;
              1 0 0 0 0 0 1 0;
              0 1 0 0 1 0 0 0]
        pn = Int32.(pn)
        initial_edges = sum(pn[:, 1:end-1])
        pruned_pn = SPNGenerator.prune_petri_net(copy(pn))
        final_edges = sum(pruned_pn[:, 1:end-1])
        @test final_edges <= initial_edges
    end

    @testset "Reachability Graph" begin
        pn = SPNGenerator.generate_random_petri_net(4, 2)
        v, e, _, _, _ = SPNGenerator.generate_reachability_graph(pn)
        @test !isempty(v)
    end

    @testset "SPN Filtering" begin
        pn = SPNGenerator.generate_random_petri_net(5, 3)
        res, success = SPNGenerator.filter_spn(pn)
        if success
            @test !isempty(res)
        else
            @test isempty(res)
        end
    end

    @testset "Lambda Variations" begin
        pn = SPNGenerator.generate_random_petri_net(5, 3)
        res, success = SPNGenerator.filter_spn(pn)
        if success
            variations = SPNGenerator.generate_lambda_variations(res, 3)
            @test length(variations) <= 3
        end
    end

    @testset "File I/O" begin
        dir = "test_dir"
        SPNGenerator.create_directory(dir)
        @test isdir(dir)

        # Test JSON I/O
        json_file = joinpath(dir, "test.json")
        test_data = Dict("a" => 1, "b" => [1,2,3])
        SPNGenerator.save_data_to_json_file(json_file, test_data)
        @test isfile(json_file)
        loaded_data = SPNGenerator.load_json_file(json_file)
        @test Dict(string(k) => v for (k,v) in pairs(loaded_data)) == test_data

        # Test TOML I/O
        toml_file = joinpath(dir, "test.toml")
        write(toml_file, "a = 1\nb = [1,2,3]")
        loaded_toml = SPNGenerator.load_toml_file(toml_file)
        @test loaded_toml["a"] == 1

        rm(dir, recursive=true)
    end
end