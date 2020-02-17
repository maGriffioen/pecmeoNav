module NaviSimu

using Plots
using LinearAlgebra
using SparseArrays
using DoubleFloats

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
    cartesianToKepler,
    orbitalPeriod,
    createCircPecmeo,
    localPosition,
    globalPosition,
    localState,
    globalState,
    lunarnavPECMEODesigner,

    #bodies.jl
    earth,
    moon,
    bodyPosition,
    lunarOrbit,
    lightConst,

    # ephemeris.jl
    Ephemeris,
    KeplerEphemeris,
    KeplerEphemerisSD,
    trueKeplerEphemeris,
    noisyKeplerEphemeris,

    #navLink.jl
    RangeTransmitterSettings,
    RangeReceiverSettings,

    #naviTools.jl
    hasLineOfSightEarth,
    hasLineOfSight,
    findNavGDOP,
    findNavPDOP,
    pointPosition,
    sequentialPointPosition,
    kinematicEstimation,

    #naviSignals.jl
    # transmitterFinder,
    instantMeasurements,
    simulateMeasurements,

    #integration.jl
    integratorRK4,
    integration_interpolationRK4,

    #mathTools.jl
    movingAverage,
    mean,
    lin2log,
    log2lin,


    #io/*
    openOrbitData,
    # plotConstellation,
    # plotConstellationConnections,
    # animateConstellation


    Tup3d

Tup3d = Tuple{<:Number, <:Number, <:Number}
NaviSimuAsModule = true

include("keplerOrbits.jl")
include("bodies.jl")
include("ephemeris.jl")
include("navLink.jl")
include("naviTools.jl")
include("naviSignals.jl")
include("integration.jl")
include("mathTools.jl")
include("io/plotTools.jl")
include("io/gpsData.jl")


end

using Main.NaviSimu
