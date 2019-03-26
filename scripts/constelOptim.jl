using Plots, LinearAlgebra, BlackBoxOptim, ProgressMeter

include("../src/NaviSimu.jl")
using Main.NaviSimu

# Best constellation (3, 3, 3) (for now) @ max_gdop = 1104.8
# 0.804656413854996
# 0.2625454045892302
# 0.29464753687011463
# 0.8550686284343616
# 0.4034387797971966
# 0.11476918737457498
# 0.5114186508013766
# 0.04146840843828451
# spacing = (1.6852684522871757, 0.5498738095275236, 0.6171083581529865)
# shift = (5.372554642808982, 2.5348806135682063, 0.7211160718288704)
# rota = 0.8033345381332042
# incl = 0.0651384236528878

function optimDesignToPecmeo(x; rawparam = false)
    #Transform normalized design vector to constellation parameters
    nsat = (2, 2, 2)
    maxspacing = (2 * pi) ./ nsat
    spacing = maxspacing .* (x[1], x[2], x[3])
    shift =  (2* pi) .*(x[4], x[5], x[6])
    rota = x[7] * pi/2
    incl = x[8] * pi/2
    r = 26400e3
    if rawparam
        return (spacing = spacing, shift = shift,
            rotation = rota, inclination = incl)
    end

    # Create PECMEO constellation
    pecmeo = createCircPecmeo(r, nsat, earth, spacing; initialOrbitShift=shift,
        inclination = incl, equatorialRotation = rota)
    return pecmeo
end

function evalLunarPecmeo(x; gdop = true, los = true, fast_max = false)
    #Create pecmeo from normalized design vector
    pecmeo = optimDesignToPecmeo(x)

    dt = 1000   #Simulation timestep
    timevec_pecmeo = 0:dt:orbitalPeriod(lunarOrbit)
    if fast_max
        max_gdop = 0
        for time in timevec_pecmeo
            gdop = findNavGDOP(bodyPosition(moon, time), propagateKeplerOrbit(pecmeo, time))
            if gdop > max_gdop
                max_gdop = gdop
            end
        end
        gdop_pecmeo = []
        los = []
    else
    if gdop
        gdop_pecmeo = map(x -> findNavGDOP(bodyPosition(moon, x), propagateKeplerOrbit(pecmeo, x)), timevec_pecmeo)
    else
        gdop = []
    end
    if los
        los = map(x-> sum(hasLineOfSightEarth(bodyPosition(moon, x), position(propagateKeplerOrbit(pecmeo, x)))), timevec_pecmeo)
    else
        los = []
    end
    max_gdop = maximum(gdop_pecmeo)
    end
    return (t = timevec_pecmeo, gdop = gdop_pecmeo, los = los, max_gdop = max_gdop)
end

function evalLunarPecmeoFitness1d(x)
    res = evalLunarPecmeo(x; fast_max = true)

    return res.max_gdop
end
function evalLunarPecmeoFitness2d(x)
    res = evalLunarPecmeo(x; los = false)

    return (maximum(res.gdop), sum(res.gdop) ./ length(res.gdop))
end

f = []
c = []
p = plot(layout = (2, 1), legend = false)

@showprogress for i in 1:10
    use2d = true

    if use2d
        #Multi objective constellation optimization
        global res = bboptimize(evalLunarPecmeoFitness2d; Method=:borg_moea,
                    FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                    SearchRange=(0.0, 1.0), NumDimensions=8, ϵ=0.05,
                    MaxSteps=50000, TraceInterval=5.0, TraceMode=:compact);
    else
        #Single objective constellation optimization
        global res = bboptimize(evalLunarPecmeoFitness1d;
                    SearchRange=(0.0, 1.0), NumDimensions=8, ϵ=0.05,
                    MaxSteps=50000, TraceInterval=5.0, TraceMode=:compact, populationSize = 200);
    end

    global f = vcat(f, best_fitness(res))
    global c = vcat(c, Tuple(best_candidate(res)))
    q = evalLunarPecmeo(best_candidate(res))
    plot!(q.t, [q.gdop q.los] ,layout = (2, 1), legend = false)
end
