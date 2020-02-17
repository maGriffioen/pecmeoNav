# Script used to compare GDOP plots of various PECMEO constellations
#
#
using Plots, LinearAlgebra
pyplot()
include("../src/NaviSimu_adder.jl")
using Main.NaviSimu



pecmeo_nav = lunarnavPECMEODesigner( [0.431459575985078,
 0.999729596168082,
 0.962413062424563,
 0.999922037056928])
pecmeo_level = lunarnavPECMEODesigner( [0.431205002151572,   0.511002410029895,   0.285692220934846,   0.639705340836302])



dt = 10
opc = orbitalPeriod(pecmeo_nav[1])
oc = 1:trunc(Int, opc/dt)
opm = orbitalPeriod(lunarOrbit)
om = 1:trunc(Int, opm/dt)

gdop_moon(navcon::Orbit, timevec) = map(x -> findNavGDOP(bodyPosition("Moon", x), propagateKeplerOrbit(navcon, x)), timevec)

timevec_pecmeo = 0:dt:orbitalPeriod(lunarOrbit)
@time gdop_pecmeo = map(x -> findNavGDOP(bodyPosition("Moon", x), propagateKeplerOrbit(pecmeo_nav, x)), timevec_pecmeo)
@time gdop_pecmeo_level = map(x -> findNavGDOP(bodyPosition("Moon", x), propagateKeplerOrbit(pecmeo_level, x)), timevec_pecmeo)


p1 = plot(timevec_pecmeo, [gdop_pecmeo], label=["Skewed"], xlabel="Time [s]", ylabel="GDOP", title="T of moon, comparison of constellations")
p2 = plot(timevec_pecmeo, [gdop_pecmeo_level], label=["Level"], xlabel="Time [s]", ylabel="GDOP", title="T of moon, comparison of constellations")
p3 = plot(timevec_pecmeo, [gdop_pecmeo gdop_pecmeo_level], label=["Skewed", "Level"], xlabel="Time [s]", ylabel="GDOP", title="T of moon, comparison of constellations")
