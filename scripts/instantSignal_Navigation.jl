using Plots, ProgressMeter, LinearAlgebra, Random, DoubleFloats
include("../src/NaviSimu_lessModule.jl")
# using Main.NaviSimu

include("gpsKepler.jl")
# include("../src/keplerOrbits.jl")
# include("../src/bodies.jl")
# include("../src/naviTools.jl")
# include("../src/naviSignals.jl")
# include("../src/io/plotTools.jl")
# include("../src/io/gpsData.jl")


#Example orbit / data from mission geometfry and orbit design
#   Intl Space Station
mu = 398600.441e9
iss = KeplerOrbit(6787746.891, 0.000731104,
    deg2rad(51.68714486), deg2rad(127.5486706),
    deg2rad(74.21987137), deg2rad(24.10027677), earth)
moonSat = KeplerOrbit(moon.radius + 100e3, 0.0,
    deg2rad(0), 0.0, 0.0, 0.0, moon)
moonSat_p = KeplerOrbit(moon.radius + df64"100e3", df64"0.0", df64"0.0", df64"0.0", df64"0.0", df64"0.0", moon)
KeplerOrbit(1.8381e6, 0.0, 0.0, 0.0, 0.0, 0.0, moon)
#Pecmeo after optimization run
pecmeo_333 = createCircPecmeo(26.4e6, (3, 3, 3), earth;
    initialOrbitShift=(0.0, (2/6)*pi, (2/6)*pi),
    equatorialRotation = 0*pi/8,
    inclination = 0.0)
pecmeo_333_p = createCircPecmeo(df64"26.4e6", (3, 3, 3), earth;
    initialOrbitShift=(df64"0.0", 2 * Double64(pi) /6, 2 * Double64(pi) /6),
    equatorialRotation = df64"0.0",
    inclination = df64"0.0")
#Pecmeo of 3x2 satellites
pecmeo_222 = createCircPecmeo(26.4e6, (2, 2, 2), earth,
    (2.5891201128936614, 2.496350075955589, 0.6181176114714648);
    initialOrbitShift = (3.2100811591091225, 4.798965601931746, 0.5712986167177687),
    equatorialRotation = 0.059619646441608075,
    inclination = -0.24798)


# gpsfile = "data/cod20000.eph"
# if !(:gpsdata in names(Main))
#    @time gpsdata = openOrbitData(gpsfile)
# end
# gpsPositionData(i) = [(gpsdata.x[i, prn], gpsdata.y[i, prn], gpsdata.z[i, prn]) for prn in 1:32]
# gps1 = gpsPositionData(1)

### Set up epochs and time vector ###
nepochs = 2880      #Number of time steps
timestep = 30    #seconds
t0 = 0
timevec = t0 .+ ((1:nepochs).-1) * timestep
timevec_org = collect(timevec)
### Set up constellation and measurement storage ###
navcon = pecmeo_333
receiverOrbit = moonSat

# Generate true ephemeris for each navigation satelltie
navconEphemeres = [trueKeplerEphemeris([0, 3600], navcon[i]) for i in 1:size(navcon)]

# Generate distubed ephemeres
# navconEphemeres = [noisyKeplerEphemeris([0, 3600], navcon[i], KeplerEphemerisSD(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)) for i in 1:size(navcon)]

nsats = size(navcon) #Total number of satellites in the navigation constellation
codeSig = Array{Float64}(undef, nepochs, nsats)    #Per-epoch, per-prn code measurements
phaseSig = Array{Float64}(undef, nepochs, nsats)   #Per-epoch, per-prn phase measurements
availability = BitArray(undef, nepochs, nsats)     #Per-epoch, per-prn availability boolean

### Real satellite positon ###
truePositions = [globalPosition(receiverOrbit, t) for t in timevec]
truePositions_p = [globalPosition(moonSat_p, Double64(t)) for t in timevec]
truePosError = [norm(truePositions[i] .- truePositions_p[i]) for i in 1:length(truePositions)]
pdop = [findNavPDOP(truePositions[epoch], globalPosition(navcon, timevec[epoch])) for epoch in 1:nepochs]
gdop = [findNavGDOP(truePositions[epoch], globalPosition(navcon, timevec[epoch])) for epoch in 1:nepochs]

conPosErrors = []
conPosErrorRSS = []
### Generation of measurements ###
for epoch in 1:nepochs
   conPos = globalPosition(navcon, timevec[epoch])         #Constellation positions
   conPos_p = globalPosition(pecmeo_333_p, Double64(timevec[epoch]))
   append!(conPosErrors, [norm(conPos[i] .- conPos_p[i]) for i in 1:length(conPos)])
   append!(conPosErrorRSS, norm([norm(conPos[i] .- conPos_p[i]) for i in 1:length(conPos)]))

   curSignals = instantMeasurements(truePositions[epoch], conPos, timevec[epoch])

   #Store data
   codeSig[epoch, :] = curSignals.code
   phaseSig[epoch, :] = curSignals.phase
   availability[epoch, :] = curSignals.avail
end
conPosErrors = reshape(conPosErrors, (size(navcon), nepochs))

# Remove epochs with too few satellites --- TODO: Improve
satCount = [sum(availability[e, :]) for e in 1:size(availability)[1]]
usableEpoch = satCount .>4
nepochs = sum(usableEpoch)
timevec = collect(timevec)[usableEpoch]
codeSig = codeSig[usableEpoch, :]
phaseSig = phaseSig[usableEpoch,:]
availability = availability[usableEpoch,:]
truePositions = truePositions[usableEpoch]
pdop = pdop[usableEpoch]
gdop = gdop[usableEpoch]


Random.seed!(1)

#Convert code and phase measurements to ranges and phase lengths (convert to meters)
freq = 1575.42e6
waveLen = lightConst / freq
pseudoRanges = codeSig * lightConst
pseudoRanges += randn(size(pseudoRanges)) .* (pseudoRanges.!=0.0) * 1   #1m normal errors
phaseLens = phaseSig * waveLen
phaseLens += randn(size(phaseLens)) .* (phaseLens.!=0.0) *1e-3        #1mm normal errors

### Point position estimation ###
ppApriori = vcat([vcat([j for j in bodyPosition(moon, timevec[i])], 0.0)' for i in 1:nepochs]...)
ppesti = sequentialPointPosition(timevec, navconEphemeres, pseudoRanges, availability;
aprioriEstimations = ppApriori, maxIter = 100, lighttimeCorrection = false).estimation
ppPosErrors = [norm(truePositions[i] .- ppesti[i][1:3]) for i in 1:nepochs]

### Kinematic estimation
@time kinEstim = kinematicEstimation(navconEphemeres, timevec, pseudoRanges, phaseLens, availability;
    ppApriori = ppApriori, maxIter_pp=100, maxIter_kin = 5, codeWeight = 1.0, phaseWeight = 1.0e6,
    lighttimeCorrection = false, solver = 2)
kinPosTime = kinEstim.positionTimeEstimation
kinBias = kinEstim.biasEstimation

kinPosErrors = [norm(truePositions[e] .- kinPosTime[e][1:3]) for e in 1:nepochs]


### Average 3d position accuracies for estimation methods
acc1 = sum(ppPosErrors) / length(ppPosErrors)      #pdop * normal error
acc2 = sum(kinPosErrors) / length(kinPosErrors)
print("\n PointPosi error: \t", acc1)
print("\n Kinematic error: \t", acc2)


# scatter!(timevec/3600, kinPosErrors, marker=(:dot, 3), lab="Navigation Error, no time dependent error")
# # plot!(timevec/3600, pdop/1000000, lab="PDOP/1000")
# xaxis!("Time since epoch [h]")
# yaxis!("Navigation error [m] AND pdop/1E7 [-]")
# plot!(legend=:topleft)
