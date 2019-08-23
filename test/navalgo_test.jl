include("../src/NaviSimu_adder.jl")
using Main.NaviSimu, Test

@testset "Position Algorithms" begin
    @testset "Reader Pieter Visser" begin
        # First, some test from Reader_AE2223-II_Analysis_4April2016.pdf (by Prof. Visser)
        gps_positions = [(2.0e7, 0.0, 0.0), (-2.0e7, 0.0, 0.0), (0.0, 2.0e7, 0.0), (0.0, -2.0e7, 0.0), (0.0, 0.0, 2.0e7)]
        observations = [19.999e6, 20.001e6, 20.000000025e6, 20.000000025e6, 20.000000025e6]
        # Tests for point position solution
        @test NaviSimu.pointPositionIteration((0.0, 0.0, 0.0, 0.0), observations, gps_positions, 0.0) ≈
            [1e3, 0.0, -0.0125, 0.0125] atol = 1e-8
        @test [i for i in pointPosition(observations, gps_positions, 0.0; maxIter = 5, correctionLimit = 1e-20, lighttimeCorrection=false).estimation] ≈
            [1000.0, 0.0, 0.0, 0.0] atol = 1e-8

        #Tests example calcultion DOP
        @test findNavGDOP((1000.0, 0.0, 0.0), gps_positions; checkLineOfSight = false) ≈ 1.58 atol = 1e-2
        @test findNavPDOP((1000.0, 0.0, 0.0), gps_positions; checkLineOfSight = false) ≈ 1.5 atol = 1e-8

        # Tests for DOPS in another reference system

        # Tests TDOP
    end

    @testset "Line of slight" begin
        # Tests for single satellites with line of sight

        # Tests for multiple satellites without line of sight

        # Test for multiple satellties along body with and without line of sight

        #Test for line of sight with a constellation
    end

    @testset "Perfect measurements" begin
        # Test with generating perfect measurements

        # Test with the point position solution

        # Test with kinematic algorithm
    end

    @testset "Nonperfect measurements" begin
        # Tests with generating simualted measureemnts

        # Tests with point position solutions from this

        # Tests with kinematic algorithm
    end

    @testset "Lighttime correction" begin
        #Use SPACE to verify the function findTransmitter for lighttime effect
    end

end
