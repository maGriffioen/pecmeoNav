using Plots, ProgressMeter, LinearAlgebra, Random
include("../src/NaviSimu_adder.jl")
using Main.NaviSimu

# include("../src/keplerOrbits.jl")
# include("../src/bodies.jl")
# include("../src/naviTools.jl")
# include("../src/naviSignals.jl")
# include("../src/io/plotTools.jl")
# include("../src/io/gpsData.jl")



function distances(refPoint, pointVector)
   distVec = map(x -> norm(pointVector[x] .- refPoint), 1:length(pointVector))
   return distVec
end

function kinematicIter(aPriEst_c, rangeData, phaseData, availability)

   avail = availability
   n_epochs = size(avail, 1)
   n_sats = size(avail, 2)
   n_measurements = sum(avail) *2
   prns = collect(1:size(avail)[2])

   #bias filter
   global biasNumber = zeros(Int16, n_epochs, maximum(prns))
   for epoch in 1:n_epochs
      for prn in prns[avail[epoch, :]]
         if epoch > 1 && biasNumber[epoch-1, prn] > 0
            #Continue using old bias if a number has been assigned in previous epoch
            biasNumber[epoch, prn] = biasNumber[epoch-1, prn]
         else
            #Create new bias Number
            biasNumber[epoch, prn] = maximum(biasNumber) + 1
         end
      end
   end

   n_bias = maximum(biasNumber)
   #global designMat = hcat(zeros(epochs * sats*2, epochs * 4), zeros(epochs*sats*2, nbias))

   global designMat = zeros(Float64,  n_measurements, n_epochs *4 + n_bias)
   measurements = zeros(Float64, n_measurements)
   model = zeros(Float64, n_measurements)
   weightMatrix = Matrix{Float64}(I, n_measurements, n_measurements)

   aPriEst = aPriEst_c
   while length(apriEst) < 4*n_epochs
      apriEst.append(0.0)
   end
   aPrioriLength = 4*n_epochs + n_bias
   while length(apriEst) < aPrioriLength
      iAddedBias = length(apriEst) - 4*n_epochs +1
      dataFilter = biasNumber.==iAddedBias
      aPriBias = sum((rangeData[dataFilter] - phaseData[dataFilter] )) / sum(dataFilter)
      append!(aPriEst, aPriBias)
   end

   rowIter = 1
   #Loop over all epochs
   for epoch in 1:n_epochs
      gpspos_e = gpsPositionData(epoch) #navigation constellation positions
      apri_e = aPriEst_c[4*epoch-3:4*epoch]   #apriori pos en t estimation for this epoch

      #Loop over available satellites
      for prn in prns[avail[epoch, :]]
         navsatpos = gpspos_e[prn]
         posdiff = apri_e[1:3].- collect(navsatpos)
         r = norm(posdiff) #Geometric range
         global sats_uvec = posdiff/r

         #Add phase measurement to design matrix and weight matrix
         designMat[rowIter, 4*epoch-3:4*epoch] = vcat(sats_uvec, 1.0)
         designMat[rowIter, n_epochs*4 + biasNumber[epoch, prn]] = -1   #bias
         weightMatrix[rowIter, rowIter] = 1e3      #Weight for phase measurement

         #Add code measurement to design matrix and weight matrix
         designMat[rowIter+1, 4*epoch-3:4*epoch] = vcat(sats_uvec, 1.0)
         weightMatrix[rowIter+1, rowIter+1] = 1    #Weight for code measurement

         measurements[rowIter] = phaseData[epoch, prn]
         measurements[rowIter+1] = rangeData[epoch, prn]

         model_phase = r + apri_e[4] - aPriEst[n_epochs*4 + biasNumber[epoch, prn]] #model for phase
         model[rowIter] = model_phase
         model_pseudorange = r + apri_e[4]               #model for pseudo range
         model[rowIter+1] = model_pseudorange

         rowIter += 2
      end
   end

   # global curr_error = (reshape(phaseLens', iter[end] * sats) - modelRange)
   # correction = pinv(designMat' * designMat) * designMat' * (reshape(phaseLens', iter[end] * sats) - modelRange)
   # correction = pinv(designMat' * designMat) * designMat' * curr_error
   global curr_error = measurements-model
   correction = inv(designMat' *weightMatrix* designMat) * designMat' *weightMatrix* (measurements - model)
return correction

end

#Example orbit / data from mission geometfry and orbit design
#   Intl Space Station
mu = 398600.441e9
iss = KeplerOrbit(6787746.891, 0.000731104,
    deg2rad(51.68714486), deg2rad(127.5486706),
    deg2rad(74.21987137), deg2rad(24.10027677), earth)


# Lunar gps GDOP calculation without shadowing
gpsfile = "data/cod20000.eph"
if !(:gpsdata in names(Main))
   @time gpsdata = openOrbitData(gpsfile)
end
gpsPositionData(i) = [(gpsdata.x[i, prn], gpsdata.y[i, prn], gpsdata.z[i, prn]) for prn in 1:32]
gps1 = gpsPositionData(1)

nepochs = 96
nsats = 32
codeSig = Array{Float64}(undef, nepochs, nsats)
phaseSig = Array{Float64}(undef, nepochs, nsats)
availability = BitArray(undef, nepochs, nsats)

timevec = ((1:nepochs).-1) * 15 * 60
trueIss = [position(propagateKeplerOrbit(iss, t)) for t in timevec]
pdop = [findNavPDOP(trueIss[epoch], gpsPositionData(epoch)) for epoch in 1:nepochs]

for epoch in 1:nepochs
   gpsPos = gpsPositionData(epoch)        #Constellation positions
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

#Perform point position estimation and compare position error (magnitude seems to be 1e-9)
ppesti = [pointPosition(pseudoRanges[epoch, availability[epoch, :]],
            gpsPositionData(epoch)[availability[epoch, :]]; niter=5) for epoch in 1:nepochs]
ppPosErrors = [norm(trueIss[i] .- ppesti[i][1:3]) for i in 1:nepochs]

#Perform kinematic position estimation
apriPosTime = vcat(map(x-> collect(ppesti[x]), 1:nepochs)...)
apriBias = (pseudoRanges[1,:] - phaseLens[1,:])
apriEst = vcat(apriPosTime, apriBias)

kinEstim = apriEst
corrections = []

@progress for i in 1:3
      #weights, singulariteit
   corr = kinematicIter(kinEstim, pseudoRanges,  phaseLens, availability)
   append!(corrections, norm(corr))
   global kinEstim += corr
end
kinPosErrors = [norm(trueIss[i] .- kinEstim[4*i-3:4*i-1]) for i in 1:nepochs]


acc1 = sum(ppPosErrors) / length(ppPosErrors)      #pdop * normal error
acc2 = sum(kinPosErrors) / length(kinPosErrors)
