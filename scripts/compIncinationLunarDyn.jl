# Script used to compare GDOP plots of various PECMEO constellations
#
#
using Plots, LinearAlgebra

include("../src/NaviSimu.jl")
using Main.NaviSimu

n = 3
# pecmeo = createCircPecmeo(26400e3, (2, 2, 2), earth, ((2/3)*pi, (2/3)*pi, (2/3)*pi); initialOrbitShift=(0.0, 0.0, 0.75*pi))
#Pecmeo in equatorial plane, no shift
pecmeo1 = createCircPecmeo(26.4e6, (3, 3, 3), earth;
    initialOrbitShift=(0.0, 0.0, 0*(2/6)*pi),
    equatorialRotation = 0*pi/8,
    inclination = 0.0)
#Pecmeo in equatorial plane, relative shift
pecmeo2 = createCircPecmeo(26.4e6, (3, 3, 3), earth;
    initialOrbitShift=(0.0, 0.0, (2/6)*pi),
    equatorialRotation = 0*pi/8,
    inclination = 0.0)
#Pecmeo in lunar plane
pecmeo3 = createCircPecmeo(26.4e6, (3, 3, 3), earth;
    initialOrbitShift=(0.0, 0.0, (2/6)*pi),
    equatorialRotation = 0*pi/8,
    inclination = deg2rad(23.4 + 5.14))
#Pecmeo after optimization run
pecmeo4 = createCircPecmeo(26.4e6, (3, 3, 3), earth,
    (5.372554642808982, 2.5348806135682063, 0.7211160718288704);
    initialOrbitShift=(5.372554642808982, 2.5348806135682063, 0.7211160718288704),
    equatorialRotation = 0*pi/8,
    inclination = 0.0651384236528878)
#Pecmeo of 3x2 satellites
pecmeo5 = createCircPecmeo(26.4e6, (2, 2, 2), earth,
    (2.5891201128936614, 2.496350075955589, 0.6181176114714648);
    initialOrbitShift = (3.2100811591091225, 4.798965601931746, 0.5712986167177687),
    equatorialRotation = 0.059619646441608075,
    inclination = 0.16112265271507054)


dt = 100
opc = orbitalPeriod(pecmeo1[1])
oc = 1:trunc(Int, opc/dt)
opm = orbitalPeriod(lunarOrbit)
om = 1:trunc(Int, opm/dt)

timevec_pecmeo = 0:dt:orbitalPeriod(lunarOrbit)
@time gdop_pecmeo1_nd = map(x -> findNavGDOP(bodyPosition("Moon", 0), propagateKeplerOrbit(pecmeo1, x)), timevec_pecmeo)
@time gdop_pecmeo2_nd = map(x -> findNavGDOP(bodyPosition("Moon", 0), propagateKeplerOrbit(pecmeo2, x)), timevec_pecmeo)
@time gdop_pecmeo3_nd = map(x -> findNavGDOP(bodyPosition("Moon", 0), propagateKeplerOrbit(pecmeo3, x)), timevec_pecmeo)
@time gdop_pecmeo1_id = map(x -> findNavGDOP(bodyPosition("Moon", x), propagateKeplerOrbit(pecmeo1, x)), timevec_pecmeo)
@time gdop_pecmeo2_id = map(x -> findNavGDOP(bodyPosition("Moon", x), propagateKeplerOrbit(pecmeo2, x)), timevec_pecmeo)
@time gdop_pecmeo3_id = map(x -> findNavGDOP(bodyPosition("Moon", x), propagateKeplerOrbit(pecmeo3, x)), timevec_pecmeo)
@time gdop_pecmeo4_id = map(x -> findNavGDOP(bodyPosition("Moon", x), propagateKeplerOrbit(pecmeo4, x)), timevec_pecmeo)
@time gdop_pecmeo5_id = map(x -> findNavGDOP(bodyPosition("Moon", x), propagateKeplerOrbit(pecmeo5, x)), timevec_pecmeo)


p1 = plot(timevec_pecmeo[oc], [gdop_pecmeo1_nd[oc] gdop_pecmeo2_nd[oc]], label=["No shift", "Shifted"], xlabel="Time [s]", ylabel="GDOP", title="T of constellation, static moon")
p2 = plot(timevec_pecmeo, [gdop_pecmeo2_nd gdop_pecmeo2_id], label=["Static Moon", "Moving Moon"], xlabel="Time [s]", ylabel="GDOP", title="T of moon, comparison of moving vs static moon")
p3 = plot(timevec_pecmeo, [gdop_pecmeo2_nd gdop_pecmeo2_id gdop_pecmeo3_id], label=["Static Moon", "Moving Moon, equatorial", "Moving Moon, lunar"], xlabel="Time [s]", ylabel="GDOP", title="T of moon, comparison of orbital planes")
p4 = plot(timevec_pecmeo, [gdop_pecmeo2_nd gdop_pecmeo2_id gdop_pecmeo4_id], label=["Static Moon", "Moving Moon, equatorial", "Moving Moon, optimized"], xlabel="Time [s]", ylabel="GDOP", title="T of moon, comparison of constellations")
p5 = plot(timevec_pecmeo, [gdop_pecmeo4_id gdop_pecmeo5_id], label=["Optimized, 3x3", "Optimized, 3x2"], xlabel="Time [s]", ylabel="GDOP", title="T of moon, comparison of constellations")
