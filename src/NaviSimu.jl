module NaviSimu

using Plots
using LinearAlgebra

export
    Tup3d,
    #keplerOrbits.jl
    Body,
    KeplerOrbit,
    CartesianOrbit,
    CartesianState,
    KeplerConstellation,
    propagateKeplerOrbit,
    keplerToCartesian,
    orbitalPeriod,
    createCircPecmeo,

    #bodies.jl
    earth,
    moon,
    bodyPosition,
    lunarOrbit,
    lightConst,

    #naviTools.jl
    hasLineOfSightEarth,
    hasLineOfSight,
    findNavGDOP,
    findNavPDOP,
    pointPosition,

    #naviSimu.jl
    instaSignal,

    #io/*
    openOrbitData,
    gpsKepler,
    plotConstellation,
    plotConstellationConnections,
    animateConstellation

Tup3d = Tuple{Float64, Float64, Float64}

include("keplerOrbits.jl")
include("bodies.jl")
include("naviTools.jl")
include("naviSignals.jl")
include("io/plotTools.jl")
include("io/gpsData.jl")

end
