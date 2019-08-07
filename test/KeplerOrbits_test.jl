include("../src/NaviSimu_adder.jl")
using Main.NaviSimu, Test, DoubleFloats

@testset "KeplerOrbits" begin

    # Definition of Earth as per Orbit and constellation design and management book
    earth_ocdm_statefun(t::Number) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    earth_ocdm = NaviSimu.Body(name="earth_ocdm", gravitationalParameter=df64"398600.441e9", radius=6378.136e3, stateFunction = earth_ocdm_statefun)

    # ISS orbit (conversion) tests
    orbit_iss_cartesian = CartesianOrbit(-2700816.14, -3314092.80, 5266346.42, 5168.606550, -5597.546618, -868.878445, earth_ocdm)
    orbit_iss_kepler = KeplerOrbit(6787746.891, 0.000731104, deg2rad(51.68714486), deg2rad(127.5486706),
        deg2rad(74.21987137), deg2rad(24.1002767), earth_ocdm)

    @test keplerToCartesian(orbit_iss_kepler).x ≈ orbit_iss_cartesian.x
    @test keplerToCartesian(orbit_iss_kepler).y ≈ orbit_iss_cartesian.y
    @test keplerToCartesian(orbit_iss_kepler).z ≈ orbit_iss_cartesian.z
    @test keplerToCartesian(orbit_iss_kepler).vx ≈ orbit_iss_cartesian.vx
    @test keplerToCartesian(orbit_iss_kepler).vy ≈ orbit_iss_cartesian.vy
    @test keplerToCartesian(orbit_iss_kepler).vz ≈ orbit_iss_cartesian.vz
    @test rad2deg(NaviSimu.findEccentricAnomaly(orbit_iss_kepler)) ≈ 24.08317766
    @test rad2deg(NaviSimu.findMeanAnomaly(orbit_iss_kepler)) ≈ 24.06608426

    # Cryosat orbit tests
    orbit_cryosat_cartesian = CartesianOrbit(3126974.99, -6374445.74, 28673.59, -254.91197, -83.30107, 7485.70674, earth_ocdm)
    orbit_cryosat_ecc = df64"0.0011219"
    # The following kepler orbit is defined from the mean anomaly to not lose accuracy
    orbit_cryosat_kepler = KeplerOrbit(df64"7096137.00", orbit_cryosat_ecc, deg2rad(df64"92.0316"), deg2rad(df64"296.1384"), deg2rad(df64"120.6878"),
        NaviSimu.meanToTrueAnomaly(deg2rad(df64"239.6546"), orbit_cryosat_ecc), earth_ocdm)

    @test globalState(orbit_cryosat_cartesian) ≈ globalState(orbit_cryosat_kepler)
    b= copy(orbit_cryosat_cartesian)
    @test b === orbit_cryosat_cartesian
    @test copy(orbit_cryosat_kepler) === orbit_cryosat_kepler

    # Constellation of ISS and crysat tests
    test_constellation = KeplerConstellation(orbit_iss_kepler)
    @test size(test_constellation) == 1
    @test push!(test_constellation, orbit_cryosat_kepler)[2] == KeplerConstellation([orbit_iss_kepler, orbit_cryosat_kepler])[2]
    @test size(test_constellation) == 2
    @test test_constellation[1] == orbit_iss_kepler
    @test localPosition(test_constellation)[1] == localPosition(orbit_iss_kepler)
    @test globalState(test_constellation)[2] ≈ globalState(orbit_cryosat_cartesian)

    test_constellation_propagated = propagateKeplerOrbit(test_constellation, df64"378")
    @test test_constellation_propagated[1] == propagateKeplerOrbit(orbit_iss_kepler, df64"378")
    @test test_constellation_propagated[2] == propagateKeplerOrbit(orbit_cryosat_kepler, df64"378")
    @test globalState(test_constellation_propagated, df64"999")[1] ≈
        [4.534898326006999e6, -5.005955629866256e6, -689296.9878195832, 3133.2226716658233, 3652.2124952973027, -5961.040511467195]

    # GEO tests, including propagation
    orbit_geo = KeplerOrbit(df64"42164172.83", df64"0", df64"0", df64"0", df64"0", df64"0", earth_ocdm)
    @test orbitalPeriod(orbit_geo) ≈ 86164.1004
    @test Float64(propagateKeplerOrbit(orbit_geo, 43082.0502).tanom) ≈ pi   #DoubleFloat accuracy is too high to be approximately pi for half an orbit
    @test Float64(propagateKeplerOrbit(orbit_geo, 21541.0251).tanom) ≈ 0.5*pi
end
