using Plots, LinearAlgebra

include("bodies.jl")
include("keplerOrbits.jl")
include("mathTools.jl")
include("navGeom.jl")
include("plotTools.jl")
include("gpsData.jl")

n = 3
# pecmeo = createCircPecmeo(26400e3, (2, 2, 2), earth, ((2/3)*pi, (2/3)*pi, (2/3)*pi); initialOrbitShift=(0.0, 0.0, 0.75*pi))
pecmeo = createCircPecmeo(26.4e6, (3, 3, 3), earth;
    initialOrbitShift=(0.0, 0.0, (2/6)*pi),
    equatorialRotation = 0*pi/8,
    inclination = 0.0)

#Simple Lunar GDOP / pecmeo calculation without shadowing
timevec_pecmeo = 0:100:orbitalPeriod(pecmeo[1])
@time gdop_pecmeo = map(x -> findNavGDOP(bodyPosition("Moon", 0), propagateKeplerOrbit(pecmeo, x)), timevec_pecmeo)

#plot(timevec_pecmeo, gdop_pecmeo)



#Lunar gps GDOP calculation without shadowing
# gpsfile = "lunarGDOP/cod20000.eph"
# @time gpsdata = openOrbitData(gpsfile)
# gpsPositionData(i) = [(gpsdata.x[i, prn], gpsdata.y[i, prn], gpsdata.z[i, prn]) for prn in 1:32]
# timevec_gps = 1:96
# @time gdop_gps = map( x -> findNavGDOP(bodyPosition("Moon", 0), gpsPositionData(x)), timevec_gps)
