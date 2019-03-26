module NaviSimu

using Plots
using LinearAlgebra

export
    Body,
    KeplerOrbit,
    CartesianOrbit,
    CartesianState,
    KeplerConstellation,
    propagateKeplerOrbit,
    keplerToCartesian,
    orbitalPeriod,
    createCircPecmeo,

    earth,
    moon,
    bodyPosition,
    lunarOrbit,

    hasLineOfSightEarth,
    hasLineOfSight,
    findNavGDOP,

    openOrbitData,

    plotConstellation,
    plotConstellationConnections,
    animateConstellation


include("keplerOrbits.jl")
include("bodies.jl")
include("naviTools.jl")
include("io/plotTools.jl")
include("io/gpsData.jl")

end
