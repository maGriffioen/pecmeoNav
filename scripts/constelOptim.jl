using Plots, LinearAlgebra, BlackBoxOptim, ProgressMeter, Statistics, Distributed, SharedArrays
# addprocs(5);
include("E:/Julia/pecmeoNav/src/NaviSimu_lessModule.jl")
# using Main.NaviSimu

global r = 14000e3

function H2inc(H)
    i = acos(H[3] / norm(H))
    return i
end
function H2raan(H)
    N = cross([0.0, 0.0, 1.0], H)
    Nxy = norm(N[1:2])
    raan = atan(N[2] / Nxy, N[1] /Nxy)
    return raan
end
function lunarnavPECMEODesigner(x; radius = 14e6)
    nsat = (3, 3, 3)
    spacing = (2*pi) ./ nsat
    moon_raan = lunarOrbit.raan
    moon_inc = lunarOrbit.i
    H = sqrt(radius * earth.gravitationalParameter)
    pecmeo = KeplerConstellation()

    # Find satellites in a plane polar w.r.t the moon.
    raan1 = moon_raan
    inc1 = moon_inc + (pi/2) + (x[1]-0.5)*(pi/2)
    H1 = [H * sin(inc1) * sin(raan1),
        -H * sin(inc1) * cos(raan1),
        H*cos(inc1)]

    # Add satellites for plane 1
    for i in 1:nsat[1]
        push!(pecmeo, KeplerOrbit(radius, 0.0, inc1, raan1, 0.0, 0*spacing[1] + (i-1) * spacing[1], earth))
    end

    # Find second polar plane w.r.t. moon & orthogonal to H1
    raan2prime = raan1 + (pi/2)
    inc2prime = (pi/2)
    H2prime = [H * sin(inc2prime) * sin(raan2prime),
        -H * sin(inc2prime) * cos(raan2prime),
        H*cos(inc2prime)]
    H2 = vector_rotation(H2prime, H1, (x[2]-0.5)*(pi/2))
    raan2 = H2raan(H2)
    inc2 = H2inc(H2)

    # Add satellites for plane 2
    for i in 1:nsat[2]
        push!(pecmeo, KeplerOrbit(radius, 0.0, inc2, raan2, 0.0, x[3]*spacing[2] + (i-1) * spacing[2], earth))
    end

    # Find last orbital plane, equatorial for moon

    H3 = vector_rotation(H2, H1, (pi/2))
    H32 = H * cross(normalize(H1), normalize(H2))
    raan3 = H2raan(H3)
    inc3 = H2inc(H3)

    # Add satellites for plane 3
    for i in 1:nsat[3]
        push!(pecmeo, KeplerOrbit(radius, 0.0, inc3, raan3, 0.0, x[4]*spacing[3] + (i-1) * spacing[3], earth))
    end

    return pecmeo
end

function newPECMEODesigner_old(x)
    nsat = (3, 3, 3)
    spacing = (2*pi) ./ nsat
    moon_raan = lunarOrbit.raan
    moon_inc = lunarOrbit.i
    radius = 14e6

    pecmeo = KeplerConstellation()

    raan1 = moon_raan
    inc1 = moon_inc + (x[1]-0.5)*(pi/2)
    o1 = KeplerOrbit(radius, 0.0, inc1, raan1, 0.0, x[1]*spacing[1], earth)
    state1 = localState(o1)
    H1 = cross(state1[1:3], state1[4:6])
    H = norm(H1)

    # Add stats for plane 1
    for i in 1:nsat[1]
        push!(pecmeo, KeplerOrbit(radius, 0.0, inc1, raan1, 0.0, 0*spacing[1] + (i-1) * spacing[1], earth))
    end

    raan2prime = raan1 + (pi/2)
    inc2prime = (pi/2)
    H2prime = [H * sin(inc2prime) * sin(raan2prime),
        -H * sin(inc2prime) * cos(raan2prime),
        H*cos(inc2prime)]
    H2 = vector_rotation(H2prime, H1, (x[2]-0.5)*(pi/2))
    raan2 = H2raan(H2)
    inc2 = H2inc(H2)

    # Add stats for plane 2
    for i in 1:nsat[2]
        push!(pecmeo, KeplerOrbit(radius, 0.0, inc2, raan2, 0.0, x[3]*spacing[2] + (i-1) * spacing[2], earth))
    end

    H3 = vector_rotation(H2, H1, (pi/2))
    raan3 = H2raan(H3)
    inc3 = H2inc(H3)

    # Add stats for plane 3
    for i in 1:nsat[3]
        push!(pecmeo, KeplerOrbit(radius, 0.0, inc3, raan3, 0.0, x[4]*spacing[3] + (i-1) * spacing[3], earth))
    end

    return pecmeo
end



function optimDesignToPecmeo(x; rawparam = false)
    #Transform normalized design vector to constellation parameters
    nsat = (3, 3, 3)
    maxspacing = (2 * pi) ./ nsat
    spacing = 2*pi/3 .* (1.0, 1.0, 1.0)
    shift =  (2* pi) .*(x[1], x[2], x[3])
    rota = x[4] * pi/2
    incl = x[5] * pi/2


    if rawparam
        return (shift = shift,
            rotation = rota, inclination = incl)
    end

    # Create PECMEO constellation
    pecmeo = createCircPecmeo(r, nsat, earth, spacing; initialOrbitShift=shift,
        inclination = incl, equatorialRotation = rota)
    return pecmeo
end

function saveParetoToFile(filepathname, optimresult)
    open(filepathname*"_fitness.txt", "w") do f
        for ln in pareto_frontier(optimresult)
            n1, n2 = fitness(ln)
            write(f, "$n1, $n2 \n")
        end
    end
    open(filepathname*"_params.txt", "w") do f
        for ln in pareto_frontier(optimresult)
            n1, n2, n3, n4 = params(ln)
            write(f, "$n1, $n2, $n3, $n4 \n")
        end
    end

end

function evalLunarPecmeo(x; gdop = true, los = true, fast_max = false, dt=200)
    #Create pecmeo from normalized design vector
    pecmeo = lunarnavPECMEODesigner(x)

    # dt = 200   #Simulation timestep
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
        gdop_pecmeo = Array{Float64}(undef, length(timevec_pecmeo))
        for i in 1:length(timevec_pecmeo)
            t = timevec_pecmeo[i]
            gdop_pecmeo[i] = findNavGDOP(bodyPosition(moon, t), propagateKeplerOrbit(pecmeo, t))
        end
        # gdop_pecmeo = map(x -> findNavGDOP(bodyPosition(moon, x), propagateKeplerOrbit(pecmeo, x)), timevec_pecmeo)
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

function evalLunarPecmeo_dist(x; gdop = true, los = true, fast_max = false, dt=200)
    #Create pecmeo from normalized design vector
    pecmeo = newPECMEODesigner(x)

    # dt = 200   #Simulation timestep
    timevec_pecmeo = 0:dt:orbitalPeriod(lunarOrbit)
    if gdop
        # gdop = []
        gdop_pecmeo = SharedArray{Float64}(length(timevec_pecmeo))
        @distributed for i in 1:length(timevec_pecmeo)
            t = timevec_pecmeo[i]
            gdop_pecmeo[i] =  findNavGDOP(bodyPosition(moon, t), propagateKeplerOrbit(pecmeo, t))
        end
        # gdop_pecmeo = pmap(x -> findNavGDOP(bodyPosition(moon, x), propagateKeplerOrbit(pecmeo, x)), timevec_pecmeo)

    else
        gdop_pecmeo = []
    end
    if los
        los = map(x-> sum(hasLineOfSightEarth(bodyPosition(moon, x), position(propagateKeplerOrbit(pecmeo, x)))), timevec_pecmeo)
    else
        los = []
    end
    max_gdop = maximum(gdop_pecmeo)

    return (t = timevec_pecmeo, gdop = gdop_pecmeo, los = los, max_gdop = max_gdop)
end

function evalLunarPecmeoFitness1d(x)
    res = evalLunarPecmeo_dist(x; los = false, dt=200)
    return 0.5*mean(res.gdop) + 0.5*maximum(res.gdop)
end
function evalLunarPecmeoFitness2d(x)
    res = evalLunarPecmeo(x; los = false, dt=200)
    # println( length(res.gdop))
    # fit1 = mean(gdop[gdop.>quantile(gdop, 0.99)])
    fit2 = mean(res.gdop)
    fit3 = maximum(res.gdop)
    return(fit2, fit3)
end

f = []
c = []
ress = []
p = plot(layout = (2, 1), legend = false)

# for q in 10000e3:2000e3:26000e3
#     global r = q
if false
    use2d = true
    if use2d
        #Multi objective constellation optimization
        weightedfitness(f) = f[1] * 0.5 + f[2]*0.5
        global res2 =bboptimize(evalLunarPecmeoFitness2d; Method=:borg_moea,
                    FitnessScheme=ParetoFitnessScheme{2}(is_minimizing=true, aggregator = weightedfitness),
                    SearchRange=(0.0, 1.0), NumDimensions=4, ϵ=0.05,
                    MaxSteps=10, TraceInterval=5.0, TraceMode=:compact,PopulationSize=50,
                    maxTime = 900);
    else
        #Single objective constellation optimization
        global res1 = bboptimize(evalLunarPecmeoFitness1d;
                    SearchRange=(0.0, 1.0), NumDimensions=4, ϵ=0.05,
                    MaxSteps=50000, TraceInterval=5.0, TraceMode=:compact, populationSize = 50,
                    maxTime = 900);
    end
    #
    # global f = vcat(f, best_fitness(res))
    # global c = vcat(c, Tuple(best_candidate(res)))
    # push!(ress, res)
    # q = evalLunarPecmeo(best_candidate(res))
    # plot!(q.t, [q.gdop q.los] ,layout = (2, 1), legend = false)
# end
end
