using Plots, LinearAlgebra, Random, Statistics, DoubleFloats
# plotly()
include("../src/NaviSimu_lessModule.jl")
# using Main.NaviSimu

gpsSat = KeplerOrbit(26559.8e3, 0.00, deg2rad(55.0), deg2rad(272.85 ), 0.0, deg2rad(11.68 ), earth)
# pecmeo_level = lunarnavPECMEODesigner( [0.431205002151572,   0.511002410029895,   0.285692220934846,   0.639705340836302])

satellite = gpsSat

# ODEs to describe time dilation of satellite with respect to Earth-fixed-rotating clock
function receiverClockODE(x::Array{<:Number, 1}, t::Number)
    dxdy = zeros(length(x))
    state = globalState(satellite, t)
    velocity = norm(state[4:6])
    c = lightConst

    dxdy[1] = sqrt(1- (velocity^2 / c^2)) #Local time rate
    dxdy[2] = dxdy[1] * (1/sqrt(1- (3873.96^2 / c^2)))  #Clock time rate
    #x[1] = t',
    #x[2] = t'_c
    return dxdy
end



# (Receiver clock) times at which to perform meaasurements
measurements = 3600*0:30:3600*24
# Find the inertial times for these measurement times
@time clockIntegratorResults = integration_interpolationRK4(receiverClockODE, measurements, 0.0, 10.0; interParam=2)

# Rename these results into usable variable names
trueMeasurementTimes = clockIntegratorResults.t  # Inertial times during measurements
localReceiverMeasurementTimes = clockIntegratorResults.y[1, :]   # local time progression during measurements
receiverClockMeasurementTimes = clockIntegratorResults.y[2, :]   # Time on receiver clock during measurements. Includes clock freq correction (and clock errors).
actualClockMeasurementError = (receiverClockMeasurementTimes - trueMeasurementTimes) * lightConst
