using Plots, LinearAlgebra, BlackBoxOptim

include("bodies.jl")
include("keplerOrbits.jl")
include("mathTools.jl")
include("navGeom.jl")
include("plotTools.jl")
include("gpsData.jl")

function evalLunarPecmeo(x; gdop = true, los = true)
    nsat = (3, 3, 3)
    maxspacing = (2 * pi) ./ nsat
    spacing = maxspacing .* (x[1], x[2], x[3])
    shift =  (2* pi) .*(x[4], x[5], x[6])
    r = 26400e3

    # pecmeo = createCircPecmeo(r, nsat, earth, spacing; initialOrbitShift=shift)
    pecmeo = createCircPecmeo(r, nsat, earth, shift; initialOrbitShift=shift)
    moonpos = bodyPosition("Moon", 0)

    timevec_pecmeo = 0:100:42700
    if gdop
        gdop_pecmeo = map(x -> findNavGDOP(moonpos, propagateKeplerOrbit(pecmeo, x)), timevec_pecmeo)
    else
        gdop = []
    end
    if los
        los = map(x-> sum(hasLineOfSightEarth(moonpos, position(propagateKeplerOrbit(pecmeo, x)))), timevec_pecmeo)
    else
        los = []
    end
    return (t = timevec_pecmeo, gdop = gdop_pecmeo, los = los)
end

function evalLunarPecmeoFitness1d(x)
    res = evalLunarPecmeo(x; los = false)

    return maximum(res.gdop)
end
function evalLunarPecmeoFitness2d(x)
    res = evalLunarPecmeo(x; los = false)

    return (maximum(res.gdop), sum(res.gdop) ./ length(res.gdop))
end

f = []
c = []
p = plot(layout = (2, 1), legend = false)

for i in 1:2
    use2d = false

    if use2d
    global res = bboptimize(evalLunarPecmeoFitness2d; Method=:borg_moea,
                FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                SearchRange=(0.0, 1.0), NumDimensions=6, ϵ=0.05,
                MaxSteps=10000, TraceInterval=5.0, TraceMode=:verbose);
    else
    global res = bboptimize(evalLunarPecmeoFitness1d;
                SearchRange=(0.0, 1.0), NumDimensions=6, ϵ=0.05,
                MaxSteps=10000, TraceInterval=5.0, TraceMode=:compact);
    end

    global f = vcat(f, best_fitness(res))
    global c = vcat(c, Tuple(best_candidate(res)))
    q = evalLunarPecmeo(best_candidate(res))
    plot!(q.t, [q.gdop q.los] ,layout = (2, 1), legend = false)
end
