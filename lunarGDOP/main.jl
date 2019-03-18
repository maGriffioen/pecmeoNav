using Plots, LinearAlgebra

include("bodies.jl")
include("keplerOrbits.jl")
include("mathTools.jl")
include("navGeom.jl")
include("plotTools.jl")
include("gpsData.jl")

n = 3
pecmeo = createCircPecmeo(26400e3, (3, 3, 3), earth)

#Simple Lunar GDOP / pecmeo calculation without shadowing
timevec_pecmeo = 0:100:28000
@time gdop_pecmeo = map(x -> findNavGDOP(bodyPosition("Moon", x), propagateKeplerOrbit(pecmeo, x)), timevec_pecmeo)
# gdop = map(x -> findNavGDOP(position(propagateKeplerOrbit(iss, x)), propagateKeplerOrbit(pecmeo, x)), timevec)

#plot(timevec_pecmeo, gdop_pecmeo)


#Pecmeo lunar GDOP calculation with shadowing
gdop_pecmeo_shadow = []
nsat_pecmeo = []

for time in timevec_pecmeo
    pecmeopos = position(propagateKeplerOrbit(pecmeo, time))
    moonPos = bodyPosition("Moon", time)
    los = hasLineOfSight(moonPos, pecmeopos, bodyPosition("Earth", 0), earth.radius)
    append!(gdop_pecmeo_shadow, findNavGDOP((384e6, 0, 0), pecmeopos[los]))
    append!(nsat_pecmeo, sum(los))
end


#Lunar gps GDOP calculation without shadowing
gpsfile = "lunarGDOP/cod20000.eph"
@time gpsdata = openOrbitData(gpsfile)
gpsPositionData(i) = [(gpsdata.x[i, prn], gpsdata.y[i, prn], gpsdata.z[i, prn]) for prn in 1:32]
timevec_gps = 1:96
@time gdop_gps = map( x -> findNavGDOP(bodyPosition("Moon", 0), gpsPositionData(x)), timevec_gps)
