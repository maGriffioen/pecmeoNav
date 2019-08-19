using Plots, LinearAlgebra, Random
plotly()
include("../src/NaviSimu_lessModule.jl")
# using Main.NaviSimu

moonSat = KeplerOrbit(moon.radius + 100e3, 0.0,
    deg2rad(0), 0.0, 0.0, 0.0, moon)
pecmeo32 = createCircPecmeo(26.4e6, (3, 3, 3), earth;
    initialOrbitShift=(0.0, (2/6)*pi, (2/6)*pi),
    equatorialRotation = 0*pi/8,
    inclination = 0.0)

#Define constellations
receiverOrbit = moonSat
navcon = pecmeo32


# ODEs to describe time dilation of satellite with respect to Earth-fixed-rotating clock
function receiverClockODE(x::Array{<:Number, 1}, t::Number)
    dxdy = zeros(length(x))
    state = globalState(receiverOrbit, t)
    velocity = norm(state[4:6])
    c = lightConst

    dxdy[1] = sqrt(1- (velocity^2 / c^2)) #Local time rate
    dxdy[2] = dxdy[1] * (1/sqrt(1- (1925.63^2 / c^2)))  #Clock time rate
    #x[1] = t',
    #x[2] = t'_c
    return dxdy
end


# Block for setting up system properties for LuNNaC and Lunar satellite
begin
    recGain(theta, phi) = 0.0                       #db
    txGain(theta, phi) = (rad2deg(theta) < 10.5) * 15.0 #db
    recPointing(pos, t) = (0.0, 0.0, 0.0)
    function txPointing_PECMEO(satellitePosition, t)
        targetPosition = bodyPosition(moon, t)
        pointingVector = (targetPosition .- satellitePosition)
        pointingVector = pointingVector ./ norm(pointingVector)
        return pointingVector
    end

    Random.seed!(0)
    operatingFrequency = 1574.42e6
    recSettings = RangeReceiverSettings( 0.0, 513.0, 2e6, 15, 20e-3, operatingFrequency, recGain,
        recPointing, 0.0)
    nSats = size(navcon)
    txSettings = [RangeTransmitterSettings( rand(Float64), 300.0, operatingFrequency, txGain,
        txPointing_PECMEO, 0.0) for prn in 1:nSats]

    lte = true
    noise = true
end

# (Receiver clock) times at which to perform meaasurements
measurements = 0.0:30:7200
# Find the inertial times for these measurement times
clockIntegratorResults = integration_interpolationRK4(receiverClockODE, measurements, 0.0, 1.0; interParam=2)
# Rename these results into usable variable names
trueMeasurementTimes = clockIntegratorResults.t  # Inertial times during measurements
localReceiverMeasurementTimes = clockIntegratorResults.y[1, :]   # local time progression during measurements
receiverClockMeasurementTimes = clockIntegratorResults.y[2, :]   # Time on receiver clock during measurements. Includes clock freq correction (and clock errors).

#Lets make some plots for the integrated times
if false
    trueOffset = results.y[1,:]- results.t
    clockOffset= results.y[2,:]- results.t

    myplot_receivertime = plot(xlabel="Time [hours]", ylabel="Time offset [s]")
    plot!(results.t / 3600, trueOffset, label="Local time offset")
    plot!(results.t / 3600, clockOffset, label="Clock offset")

    myplot_timeInterpolationError = plot(xlabel="Time [hours]", ylabel="Time offset [s]")
    plot!(inter.t / 3600, inter.y[2, :]-inter.t, label="Time offset from inertial time during measurement")
    plot!(inter.t / 3600, inter.y[2,:]-measurements, label="Numerical measurement time error")
end

# Generate code range and phase measurements
Random.seed!(2) #Seed for measurement noise
simulatedMeasurements = simulateMeasurements(trueMeasurementTimes, receiverClockMeasurementTimes, moonSat, navcon, recSettings, txSettings;
    lighttimeEffect = lte, addNoise=noise)


# Find true satellite positions
correct_positions = [globalPosition(moonSat, trueMeasurementTimes[i]) for i in 1:length(measurements)]

# Calculate delution of precission based on true locations
gdop = [findNavGDOP(correct_positions[i], globalPosition(navcon, trueMeasurementTimes[i])) for i in 1:length(measurements)]
pdop = [findNavPDOP(correct_positions[i], globalPosition(navcon, trueMeasurementTimes[i])) for i in 1:length(measurements)]
myplot_dop = plot(trueMeasurementTimes/3600, gdop, label="GDOP", yaxis="DOP [-]")
plot!(trueMeasurementTimes/3600, pdop, label="PDOP")


# Create navigation satellite ephemeres. These are the broadcast ephemeres
Random.seed!(1) #Seed for ephemeres
input_navconEphemeres = [noisyKeplerEphemeris([0, 1800, 3600, 5400], navcon[i], KeplerEphemerisSD(0.03, 0.0, 0.0, 0.0, 0.0, 0.0)) for i in 1:size(navcon)]

# Define input parameters for the position estimation
input_measurementEpochs = simulatedMeasurements.timeStamp
input_pseudoRanges = simulatedMeasurements.codeObs
input_phases = simulatedMeasurements.phaseObs
input_availability = simulatedMeasurements.availability

#Point position estimation
pointPosition_apriori = vcat([vcat([j for j in bodyPosition(moon, simulatedMeasurements.timeStamp[i])], 0.0)' for i in 1:length(measurements)]...)
pointPosition_estimation = sequentialPointPosition(input_measurementEpochs, input_navconEphemeres, input_pseudoRanges, input_availability;
    aprioriEstimations = pointPosition_apriori, maxIter = 5, lighttimeCorrection = lte)

# Process point position results
pointPosition_xyzErrors = [correct_positions[i] .- pointPosition_estimation.estimation[i][1:3] for i in 1:length(measurements)]
pointPosition_rssErrors = [norm(correct_positions[i] .- pointPosition_estimation.estimation[i][1:3]) for i in 1:length(measurements)]
pointPosition_meanRssError = sum(pointPosition_rssErrors[pointPosition_estimation.resultValidity]) / length(pointPosition_rssErrors[pointPosition_estimation.resultValidity])
println("\n PointPosi error: \t", pointPosition_meanRssError)
myplot_ppErrors = plot(input_measurementEpochs/3600, [pointPosition_rssErrors movingAverage(pointPosition_rssErrors; n=5)]./1000.0, yaxis=("Point Position Error [km]", (0, 2.7e1)), label=["Direct" "Moving Average"], w=[1.0 2.0])
plot!(input_navconEphemeres[1].timeReferences/3600, linetype=:vline, c=:black, w=0.5, label="Ephemeris update")

#Kinematic positioning. First horrible line formats point positioning results into usable apriori
kinematicApriori = reshape(collect(Iterators.flatten(pointPosition_estimation.estimation[pointPosition_estimation.resultValidity])), (4, sum(pointPosition_estimation.resultValidity)))'
kinEstim = kinematicEstimation(input_navconEphemeres, trueMeasurementTimes[pointPosition_estimation.resultValidity], input_pseudoRanges[pointPosition_estimation.resultValidity, :], input_phases[pointPosition_estimation.resultValidity, :], input_availability[pointPosition_estimation.resultValidity, :];
    ppApriori = kinematicApriori, maxIter_pp=0, maxIter_kin = 4,
    codeWeight = 1.0, phaseWeight = 1.0e6, correctionLimit_kin = 1e-3,
    lighttimeCorrection = lte)

# Process kinematic estimation results
kinematic_positionTime = kinEstim.positionTimeEstimation
kinematic_phaseBiases   = kinEstim.biasEstimation
kinematic_rssErrors = [norm(correct_positions[pointPosition_estimation.resultValidity][e] .- kinematic_positionTime[e][1:3]) for e in 1:length(kinematic_positionTime)]
kinematic_meanAccuracy = sum(kinematic_rssErrors) /  length(kinematic_rssErrors)
println("\n Kinematic Position error: \t", kinematic_meanAccuracy)
myplot_kinErrors = plot(kinEstim.kinTimes/3600, [kinematic_rssErrors movingAverage(kinematic_rssErrors; n=5)], yaxis=("Kinematic Error [m]", (0, 12)), label=["Direct" "Moving Average"], w=[1.0 2.0])
plot!(input_navconEphemeres[1].timeReferences/3600, linetype=:vline, c=:black, w=0.5, label="Ephemeris update")

# 3d plot of the errors. Dont use it. Really, it doesn't show anything interesting...
if false
    moon_pos = [bodyPosition(moon, t) for t in trueMeasurementTimes]
    x1 = moon_pos[1][1] *0
    y1 = moon_pos[2][2] *0
    z1 = moon_pos[3][3] *0
    xm = [moon_pos[i][1] for i in 1:length(moon_pos)] .- x1
    ym = [moon_pos[i][2] for i in 1:length(moon_pos)] .- y1
    zm = [moon_pos[i][3] for i in 1:length(moon_pos)] .- z1
    x = [correct_positions[i][1] for i in 1:length(correct_positions)] .-x1
    y = [correct_positions[i][2] for i in 1:length(correct_positions)] .-y1
    z = [correct_positions[i][3] for i in 1:length(correct_positions)] .-z1

    x_k = [kinematic_positionTime[i][1] for i in 1:length(kinematic_positionTime)] .-x1
    y_k = [kinematic_positionTime[i][2] for i in 1:length(kinematic_positionTime)] .-y1
    z_k = [kinematic_positionTime[i][3] for i in 1:length(kinematic_positionTime)] .-z1
    x_p = [pointPosition_estimation.estimation[i][1] for i in 1:length(pointPosition_estimation.estimation)] .-x1
    y_p = [pointPosition_estimation.estimation[i][2] for i in 1:length(pointPosition_estimation.estimation)] .-y1
    z_p = [pointPosition_estimation.estimation[i][3] for i in 1:length(pointPosition_estimation.estimation)] .-z1


    p5 = plot3d(x_k, y_k, z_k, w=1.0, label="Kinematic")
    plot3d!(x_p, y_p, z_p, w=1.0, label="Point pos")
    plot3d!(x, y, z, line=(1.0, :dot), label ="True trajectory")
    plot3d!(xm, ym, zm, label="Moon trajectory", w=5.0)
    # plotLineSphere(earth.radius*1e3, 24; center=(-x1, -y1, -z1))
    # plot3d!(xlims = (-4e8, 4e8), ylims = (-4e8, 4e8), zlims = (-4e8, 4e8))
end
