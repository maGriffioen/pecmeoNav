using Plots, ProgressMeter, LinearAlgebra, Random
include("../src/NaviSimu_adder.jl")
using Main.NaviSimu

# Example orbit / data from mission geometfry and orbit design
#   Intl Space Station
iss = KeplerOrbit(6787746.891, 0.000731104,
    deg2rad(51.68714486), deg2rad(127.5486706),
    deg2rad(74.21987137), deg2rad(24.10027677), earth)

### Set up epochs and time vector ###
nepochs = 50   #Number of time steps
timestep = 20    #seconds
timevec = ((1:nepochs).-1) * timestep

### Set up constellation ###
navcon = gpsKepler
include("gpsKepler.jl")

### Allocate empty data Matrices ###
nsats = size(gpsKepler) #Total number of satellites in the navigation constellation
codeSig = Array{Float64}(undef, nepochs, nsats)    #Per-epoch, per-prn code measurements
phaseSig = Array{Float64}(undef, nepochs, nsats)   #Per-epoch, per-prn phase measurements
availability = BitArray(undef, nepochs, nsats)     #Per-epoch, per-prn availability boolean

### 'True' satellite positon ###
trueIss = [position(propagateKeplerOrbit(iss, t)) for t in timevec]
pdop = [findNavPDOP(trueIss[epoch], position(propagateKeplerOrbit(gpsKepler, timevec[epoch]))) for epoch in 1:nepochs]

### Generation of measurements ###
for epoch in 1:nepochs
   gpsPos = position(propagateKeplerOrbit(gpsKepler, timevec[epoch]))        #Constellation positions
   curSignals = instaSignal(trueIss[epoch], gpsPos, timevec[epoch])

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
ppPosErrors = [norm(trueIss[i] .- ppesti[i][1:3]) for i in 1:nepochs]

### Kinematic estimation
kinEstim = kinematicEstimation(navcon, timevec, pseudoRanges, phaseLens, availability)
kinPosTime = kinEstim.positionTimeEstimation
kinBias = kinEstim.biasEstimation

kinPosErrors = [norm(trueIss[e] .- kinPosTime[e][1:3]) for e in 1:nepochs]


### Average 3d position accuracies for estimation methods
acc1 = sum(ppPosErrors) / length(ppPosErrors)      #pdop * normal error
acc2 = sum(kinPosErrors) / length(kinPosErrors)
