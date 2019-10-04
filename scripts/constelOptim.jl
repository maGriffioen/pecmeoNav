using Plots, LinearAlgebra, BlackBoxOptim, ProgressMeter

include("../src/NaviSimu_lessModule.jl")
# using Main.NaviSimu


function optimDesignToPecmeo(x; rawparam = false)
    #Transform normalized design vector to constellation parameters
    nsat = (3, 3, 3)
    maxspacing = (2 * pi) ./ nsat
    spacing = 2*pi/3 .* (1.0, 1.0, 1.0)
    shift =  (2* pi) .*(x[1], x[2], x[3])
    rota = x[4] * pi/2
    incl = x[5] * pi/2
    r = 14000e3
    if rawparam
        return (shift = shift,
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

    dt = 50   #Simulation timestep
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
ress = []
p = plot(layout = (2, 1), legend = false)

@showprogress for i in 1:10
    use2d = true

    if use2d
        #Multi objective constellation optimization
        global res = bboptimize(evalLunarPecmeoFitness2d; Method=:borg_moea,
                    FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true),
                    SearchRange=(0.0, 1.0), NumDimensions=5, ϵ=0.05,
                    MaxSteps=2500, TraceInterval=5.0, TraceMode=:compact);
    else
        #Single objective constellation optimization
        global res = bboptimize(evalLunarPecmeoFitness1d;
                    SearchRange=(0.0, 1.0), NumDimensions=8, ϵ=0.05,
                    MaxSteps=50000, TraceInterval=5.0, TraceMode=:compact, populationSize = 200);
    end

    global f = vcat(f, best_fitness(res))
    global c = vcat(c, Tuple(best_candidate(res)))
    push!(ress, res)
    # q = evalLunarPecmeo(best_candidate(res))
    # plot!(q.t, [q.gdop q.los] ,layout = (2, 1), legend = false)
end
