module NaviSimu

using Plots
using LinearAlgebra

export
    Tup3d,
    #keplerOrbits.jl
    Body,
    Orbit,
    KeplerOrbit,
    CartesianOrbit,
    CartesianState,
    KeplerConstellation,
    propagateKeplerOrbit,
    keplerToCartesian,
    orbitalPeriod,
    createCircPecmeo,
    localPosition,
    globalPosition,
    localState,
    globalState,

    #bodies.jl
    earth,
    moon,
    bodyPosition,
    lunarOrbit,
    lightConst,

    # ephemeris.jl
    Ephemeris,
    KeplerEphemeris,
    trueKeplerEphemeris,
    noisyKeplerEphemeris,

    #naviTools.jl
    hasLineOfSightEarth,
    hasLineOfSight,
    findNavGDOP,
    findNavPDOP,
    pointPosition,
    sequentialPointPosition,
    kinematicEstimation,


    #naviSimu.jl
    instaSignal,



    #io/*
    openOrbitData,
    # plotConstellation,
    # plotConstellationConnections,
    # animateConstellation


    Tup3d

Tup3d = Tuple{Float64, Float64, Float64}

include("keplerOrbits.jl")
include("bodies.jl")
include("ephemeris.jl")
include("naviTools.jl")
include("naviSignals.jl")
include("io/plotTools.jl")
include("io/gpsData.jl")


end
