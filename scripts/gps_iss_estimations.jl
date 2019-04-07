using Plots, ProgressMeter, LinearAlgebra, Random
include("../src/NaviSimu_adder.jl")
using Main.NaviSimu

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


# gpsfile = "data/cod20000.eph"
# if !(:gpsdata in names(Main))
#    @time gpsdata = openOrbitData(gpsfile)
# end
# gpsPositionData(i) = [(gpsdata.x[i, prn], gpsdata.y[i, prn], gpsdata.z[i, prn]) for prn in 1:32]
# gps1 = gpsPositionData(1)

### Set up epochs and time vector ###
nepochs = 50   #Number of time steps
timestep = 20    #seconds
timevec = ((1:nepochs).-1) * timestep

### Set up constellation and measurement storage ###
navcon = gpsKepler
receiverOrbit = iss
nsats = size(navcon) #Total number of satellites in the navigation constellation
codeSig = Array{Float64}(undef, nepochs, nsats)    #Per-epoch, per-prn code measurements
phaseSig = Array{Float64}(undef, nepochs, nsats)   #Per-epoch, per-prn phase measurements
availability = BitArray(undef, nepochs, nsats)     #Per-epoch, per-prn availability boolean

### Real satellite positon ###
truePositions = [globalPosition(propagateKeplerOrbit(receiverOrbit, t)) for t in timevec]
pdop = [findNavPDOP(truePositions[epoch], globalPosition(propagateKeplerOrbit(navcon, timevec[epoch]))) for epoch in 1:nepochs]

### Generation of measurements ###
for epoch in 1:nepochs
   conPos = globalPosition(propagateKeplerOrbit(navcon, timevec[epoch]))        #Constellation positions
   curSignals = instaSignal(truePositions[epoch], conPos, timevec[epoch])

   #Store data
   codeSig[epoch, :] = curSignals.code
   phaseSig[epoch, :] = curSignals.phase
   availability[epoch, :] = curSignals.avail
end

Random.seed!(1)

#Convert code and phase measurements to ranges and phase lengths (convert to meters)
freq = 1575.42e6
waveLen = lightConst / freq
pseudoRanges = codeSig * lightConst
pseudoRanges += randn(size(pseudoRanges)) .* (pseudoRanges.!=0.0) * 1   #1m normal errors
phaseLens = phaseSig * waveLen
phaseLens += randn(size(phaseLens)) .* (phaseLens.!=0.0) *1e-3       #1mm normal errors

### Point position estimation ###
ppesti = sequentialPointPosition(timevec, navcon, pseudoRanges, availability; maxIter = 5)
ppPosErrors = [norm(truePositions[i] .- ppesti[i][1:3]) for i in 1:nepochs]

### Kinematic estimation
kinEstim = kinematicEstimation(navcon, timevec, pseudoRanges, phaseLens, availability)
kinPosTime = kinEstim.positionTimeEstimation
kinBias = kinEstim.biasEstimation

kinPosErrors = [norm(truePositions[e] .- kinPosTime[e][1:3]) for e in 1:nepochs]


### Average 3d position accuracies for estimation methods
acc1 = sum(ppPosErrors) / length(ppPosErrors)      #pdop * normal error
acc2 = sum(kinPosErrors) / length(kinPosErrors)
