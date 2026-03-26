using ERGMEgo
using Test

@testset "ERGMEgo.jl" begin
    @testset "Module loading" begin
        @test @isdefined(ERGMEgo)
    end

    @testset "EgoNetwork construction" begin
        alters = [1, 2, 3]
        alter_ties = zeros(Bool, 3, 3)
        alter_ties[1, 2] = true
        alter_ties[2, 1] = true

        ego_net = EgoNetwork(0, alters, alter_ties)
        @test ego_net isa EgoNetwork{Int}
        @test length(ego_net.alters) == 3

        # Mismatched dimensions should error
        @test_throws ArgumentError EgoNetwork(0, alters, zeros(Bool, 2, 2))
    end

    @testset "EgoData construction" begin
        ego1 = EgoNetwork(1, [1, 2], zeros(Bool, 2, 2))
        ego2 = EgoNetwork(2, [1, 2, 3], zeros(Bool, 3, 3))
        ed = EgoData([ego1, ego2])
        @test ed isa EgoData
        @test length(ed) == 2
    end

    @testset "Ego ERGM terms" begin
        @test EgoEdges() isa EgoEdges
        @test EgoNodeMatch(:race) isa EgoNodeMatch
        @test EgoDegree() isa EgoDegree
        @test EgoDegree(3) isa EgoDegree
        @test EgoGWDegree() isa EgoGWDegree
        @test EgoTriangle() isa EgoTriangle
        @test EgoMixingMatrix(:gender) isa EgoMixingMatrix
    end

    @testset "Estimation API" begin
        @test ergm_ego === fit_ego_ergm
    end

    @testset "Population size estimation" begin
        @test isdefined(ERGMEgo, :estimate_popsize)
    end

    @testset "Simulation" begin
        @test isdefined(ERGMEgo, :simulate_ego_sample)
    end

    @testset "Diagnostics" begin
        @test isdefined(ERGMEgo, :ego_gof)
    end
end
