include("keplerOrbits.jl")


#Example orbit / data from mission geometry and orbit design
mu = 398600.441e9

iss = KeplerOrbit(6787746.891, 0.000731104,
    deg2rad(51.68714486), deg2rad(127.5486706),
    deg2rad(74.21987137), deg2rad(24.10027677), mu)

iss_new = propagateKeplerOrbit(iss, 2000)
