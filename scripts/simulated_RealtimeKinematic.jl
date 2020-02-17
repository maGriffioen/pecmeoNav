using Plots, LinearAlgebra, Random, Statistics, DoubleFloats, ProgressMeter
# plotly()
include("../src/NaviSimu_lessModule.jl")
# using Main.NaviSimu

moonSat = KeplerOrbit(moon.radius + 100e3, 0.0,
    deg2rad(0), 0.0, 0.0, 0.0, moon)
pecmeo_nav = lunarnavPECMEODesigner( [0.431459575985078,
 0.999729596168082,
 0.962413062424563,
 0.999922037056928])

pecmeo_level = lunarnavPECMEODesigner( [0.431205002151572,   0.511002410029895,   0.285692220934846,   0.639705340836302])

#Define constellations
receiverOrbit = moonSat
navcon = pecmeo_nav


# ODEs to describe time dilation of satellite with respect to Earth-fixed-rotating clock
function receiverClockODE(x::Array{<:Number, 1}, t::Number)
    dxdy = zeros(length(x))
    state = globalState(receiverOrbit, t)
    velocity = norm(state[4:6])
    c = lightConst

    dxdy[1] = sqrt(1- (velocity^2 / c^2)) #Local time rate
    dxdy[2] = dxdy[1] * (1/sqrt(1- (1924.2^2 / c^2)))  #Clock time rate
    #x[1] = t',
    #x[2] = t'_c
    return dxdy
end


# Block for setting up system properties for LuNNaC and Lunar satellite
begin
    recGain(theta, phi) = 0.0                       #db
    function txGain(theta, phi)
        if rad2deg(theta) < 10.5
            return 16.5 #db
        else
            return -100 #db
        end
    end
    recPointing(pos, t) = (0.0, 0.0, 0.0)
    function txPointing_PECMEO(satellitePosition, t)
        targetPosition = bodyPosition(moon, t)
        pointingVector = (targetPosition .- satellitePosition)
        pointingVector = pointingVector ./ norm(pointingVector)
        return pointingVector
    end

    Random.seed!(0)
    operatingFrequency = 1575.42e6
    receptionNoiseTemp = 290.0
    receptionBandwidth = 2e6
    trackingLoopBandwidth = 5
    predetectionIntegrationInterval = 20e-3
    # operatingFrequency = 1e9
    recSettings = RangeReceiverSettings( 0.0, receptionNoiseTemp, receptionBandwidth, trackingLoopBandwidth, predetectionIntegrationInterval, operatingFrequency, recGain,
        recPointing, 0.0)
    nSats = size(navcon)
    transmitPower = 119
    txSettings = [RangeTransmitterSettings( rand(Float64), transmitPower, operatingFrequency, txGain,
        txPointing_PECMEO, 0.0) for prn in 1:nSats]

    lte = true
    noise = true
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

# trueMeasurementTimes = collect(measurements)
# receiverclockMeasurementTimes = collect(measurements)
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

    plot(trueMeasurementTimes/3600, hcat((localReceiverMeasurementTimes-trueMeasurementTimes)*1e6, (receiverClockMeasurementTimes-trueMeasurementTimes)*1e6), xaxis="Time since epoch [h]", yaxis="Observed clock offset [ms]", label=["Uncorrected clock", "Corrected clock"], title="Time velocity dillation", legend=:bottomleft, ylims=(-0.5, 0.1))
end

# Generate code range and phase measurements
Random.seed!(2) #Seed for measurement noise
@time simulatedMeasurements = simulateMeasurements(trueMeasurementTimes, receiverClockMeasurementTimes, receiverOrbit, navcon, recSettings, txSettings;
    lighttimeEffect = lte, addNoise=noise)


# Find true satellite positions
correct_positions = [globalPosition(receiverOrbit, trueMeasurementTimes[i]) for i in 1:length(measurements)]

# Calculate delution of precission based on true locations
gdop = [findNavGDOP(correct_positions[i], globalPosition(navcon, trueMeasurementTimes[i]); checkLOSMoon = true, time = trueMeasurementTimes[i]) for i in 1:length(measurements)]
pdop = [findNavPDOP(correct_positions[i], globalPosition(navcon, trueMeasurementTimes[i]); checkLOSMoon = true, time = trueMeasurementTimes[i]) for i in 1:length(measurements)]
tdop = [findNavTDOP(correct_positions[i], globalPosition(navcon, trueMeasurementTimes[i]); checkLOSMoon = true, time = trueMeasurementTimes[i]) for i in 1:length(measurements)]
myplot_dop = plot(trueMeasurementTimes/3600, gdop, label="GDOP", yaxis="DOP [-]", ylims=(0, 2000))
plot!(trueMeasurementTimes/3600, pdop, label="PDOP")


# Create navigation satellite ephemeres. These are the broadcast ephemeres
Random.seed!(1) #Seed for ephemeres
input_navconEphemeres = [noisyKeplerEphemeris(collect(0:1800:maximum(measurements)), navcon[i], 0.01) for i in 1:size(navcon)]
if true
    ephemerisErrors = []
    r_ephemerisErrors = []
    for t in trueMeasurementTimes
        ephErrorsThisEpoch = []
        repherrorThisepoch = []
        for prn in 1:9
            push!(ephErrorsThisEpoch, norm(globalPosition(input_navconEphemeres[prn], t) .- globalPosition(navcon[prn], t)))
            r1 = norm(globalPosition(input_navconEphemeres[prn], t) .- globalPosition(receiverOrbit, t))
            r2 = norm(globalPosition(navcon[prn], t) .- globalPosition(receiverOrbit, t))
            push!(repherrorThisepoch, r1-r2)
        end
        push!(ephemerisErrors, ephErrorsThisEpoch)
        push!(r_ephemerisErrors, repherrorThisepoch)
    end
    ephemerisErrors = [y[i] for y in ephemerisErrors, i in 1:9]
    r_ephemerisErrors = [y[i] for y in r_ephemerisErrors, i in 1:9]
    myplot_ephemerisErrors = plot(trueMeasurementTimes/3600, ephemerisErrors, yaxis="Ephemeris Error [m]", xaxis="Time [h]")
end

Random.seed!(1000)

# Define input parameters for the position estimation
input_measurementEpochs = simulatedMeasurements.timeStamp
# input_measurementEpochs = trueMeasurementTimes
input_pseudoRanges = simulatedMeasurements.codeObs
input_phases = simulatedMeasurements.phaseObs
input_availability = simulatedMeasurements.availability


Random.seed!(1000)
realtimeKinSolutions = []

arcs = findMeasurementArchs(input_availability)
boolarcs = [false for i in 1:length(measurements)]
for arc in arcs
    for i in arc
        boolarcs[i] = true
    end
end
n = 0
ns = []
for i in boolarcs
    if i == 1
        global n+=1
        push!(ns, n)
    else
        global n = 0
    end
end
kinematic_times_true = trueMeasurementTimes[boolarcs]

@showprogress for arc in arcs
    arcstart = arc.start
for epoch_now in arc
    # println(epoch_now)
    passed_epochs = arcstart:epoch_now
    # Define input parameters for the position estimation
    input_measurementEpochs = simulatedMeasurements.timeStamp[passed_epochs]
    # input_measurementEpochs = trueMeasurementTimes
    input_pseudoRanges = simulatedMeasurements.codeObs[passed_epochs, :]
    input_phases = simulatedMeasurements.phaseObs[passed_epochs, :]
    input_availability = simulatedMeasurements.availability[passed_epochs, :]

    #Point position estimation
    pointPosition_apriori = vcat([vcat([j for j in bodyPosition(moon, input_measurementEpochs[i])], 0.0)' for i in 1:length(input_measurementEpochs)]...)
    pointPosition_estimation = sequentialPointPosition(input_measurementEpochs, input_navconEphemeres, input_pseudoRanges, input_availability;
        aprioriEstimations = pointPosition_apriori, correctionLimit = 1e-4, maxIter = 10, lighttimeCorrection = lte)

    #Kinematic positioning. First horrible line formats point positioning results into usable apriori
    kinematicApriori = reshape(reinterpret(Float64, pointPosition_estimation.estimation), (4, length(input_measurementEpochs)))'
    kinEstim = kinematicEstimation(input_navconEphemeres, input_measurementEpochs,
        input_pseudoRanges, input_phases,
        input_availability;
        ppApriori = kinematicApriori, maxIter_pp=0, maxIter_kin = 8,
        codeWeight = 1.0, phaseWeight = 1.0e6, correctionLimit_kin = 2e-4,
        lighttimeCorrection = lte)

    # Process kinematic estimation results
    kinematic_epochs = vcat([collect(x) for x in kinEstim.archs]...)
    kinematic_positionTime = kinEstim.positionTimeEstimation
    push!(realtimeKinSolutions, kinematic_positionTime[length(kinematic_positionTime)])
end
end

kin_cor_solution = correct_positions[boolarcs]
rtk_rssErrors = [norm(realtimeKinSolutions[i][1:3] .-kin_cor_solution[i] ) for i in 1:sum(boolarcs)]

if true
    # Reshape errors into different reference frame
    # trueMeasurementTimes_reduced = trueMeasurementTimes[pointPosition_estimation.resultValidity]
    kinematic_rvnErrors = []
    for (i, epoch) in enumerate(kinematic_epochs)
        referenceTraj = globalState(moon, trueMeasurementTimes[epoch])
        r = referenceTraj[1:3]
        v = referenceTraj[4:6]

        rp_norm = r ./ norm(r)
        n = cross(rp_norm, v)
        np_norm = n ./ norm(n)
        vp_norm = cross(np_norm, rp_norm)
        # The relative reference frame is (rp_norm, vp_norm, np_norm) -> Unit vectors:
        #       - Along r
        #       - In orbit plane (normal to r)
        #       - Normal to orbit plane
        xyzError_cur = kinematic_xyzErrors[i]
        rvnError = [dot(xyzError_cur, rp_norm), dot(xyzError_cur, vp_norm), dot(xyzError_cur, np_norm)]
        global kinematic_rvnErrors = vcat(kinematic_rvnErrors, rvnError)
    end
    kinematic_rvnErrors = reshape(kinematic_rvnErrors, (3, Int64(length(kinematic_rvnErrors)/3)))'

    kin_rvn_meanErrors = [mean(kinematic_rvnErrors[:,i]) for i in 1:3]
end


plot(trueMeasurementTimes, hcat(gdop, tdop, pdop), label=["GDOP", "TDOP", "PDOP"], xaxis="Time [h]", yaxis = "DOP [-]", legend=:right, ylims=(0, 850))
# scatter(kinematic_times_true/3600, hcat(kinematic_rvnErrors), lab=["r error", "v error", "n error"], xaxis="Time [s]", yaxis="Navigation error [m]", legend=:bottomleft)
scatter(kinematic_times_true/3600, hcat(kinematic_rvnErrors, -kinematic_timeError), lab=["r error", "v error", "n error", "-dt error"], xaxis="Time [h]", yaxis="Navigation error [m]")
scatter(kinematic_times_true/3600, kinematic_rssErrors, xaxis="Time [h]", yaxis="Navigation error [m]", legend=false)

# Calcualte pearson correlation matrix between errors)
splitCorrelation = cor(hcat(kinematic_rvnErrors, kinematic_timeError))
rssErrorCorrelation = cor(kinematic_rssErrors, kinematic_timeError)


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
