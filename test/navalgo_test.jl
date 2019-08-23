include("../src/NaviSimu_adder.jl")
using Main.NaviSimu, Test, DoubleFloats

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

    @testset "Lighttime correction SPICE test" begin
        # Use SPACE to verify the function findTransmitter for lighttime effect
        using SPICE

        # Load SPICE Kernels
        furnsh("test/de435.bsp")    # SPK
        furnsh("test/naif0012.tls") # Leap seconds kernel

        # Set an epoch
        t0 = utc2et("January 01, 2000 12:00:00")

        # Obtain uncorrected and corrected state
        earthState = spkezr("EARTH_BARYCENTER", t0, "J2000", "None", "SOLAR SYSTEM BARYCENTER")
        earthState_lt = spkezr("EARTH_BARYCENTER", t0, "J2000", "CN", "SOLAR SYSTEM BARYCENTER")

        # helicentric gravitational constant, body for the solar system and earth orbit from earth state
        hgc = df64"1.32712440042e20"
        statefun_solar(t) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        solarSystemBody = Body("solar system", hgc, 0.0, statefun_solar)
        eS = [Double64(i)*1000 for i in earthState[1][1:6]]
        keplerEarthOrbit = cartesianToKepler(CartesianOrbit(eS[1], eS[2], eS[3], eS[4], eS[5], eS[6], solarSystemBody))

        # Solved lighttime correction from the library
        transmitterEarth= NaviSimu.transmitterFinder(0.0, (0.0, 0.0, 0.0), keplerEarthOrbit; maximumCorrection=1e-20)

        # Compare results from SPICE and NaviSimu library
        @test transmitterEarth.travelTime ≈ earthState_lt[2] atol =1e-8
        @test (earthState_lt[1]*1e3 ≈ NaviSimu.state(keplerEarthOrbit, transmitterEarth.transmissionTime)) rtol =1e-9
        @test (earthState_lt[1][1:3]*1e3 ≈ [i for i in transmitterEarth.pos]) rtol = 1e-9
    end

end
