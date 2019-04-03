using Plots, ProgressMeter, LinearAlgebra, Random
include("../src/NaviSimu_adder.jl")
using Main.NaviSimu

function pointPosition(ranges, navSats; niter = 10)
   esti = [0, 0, 0, 0]
   for i in 1:niter
      esti = pointPositionIteration(esti, ranges, navSats)
   end
   return Tuple(esti)
end

function distances(refPoint, pointVector)
   distVec = map(x -> norm(pointVector[x] .- refPoint), 1:length(pointVector))
   return distVec
end

function kinematicIter_old(aPriEst_c, rangeData, phaseData)

   epochs = size(phaseLens, 1)
   sats = size(phaseLens, 2)
   global designMat = hcat(zeros(epochs * sats, epochs * 4), -ones(epochs*sats, sats))

   for e in iter
      designMat[(e-1)*32+1:e*32, (e-1)*4+1:e*4] = NaviSimu.findPPDesignMatrix(Tuple(aPriEst_c[4*e-3:4*e-1]), gpsPositionData(e))
   end

   global modelRange = []
   for e in 1:epochs
      for isat in 1:sats
         recPos = aPriEst_c[4*e-3:4*e-1]
         traPos = gpsPositionData(e)[isat]
         recClockError = aPriEst_c[4*e]
         bias = aPriEst_c[4*96+isat]
         geomrange = norm(recPos .- traPos)
         pseudorange = geomrange + recClockError - bias
         append!(modelRange, pseudorange)
      end
   end

   global curr_error = (reshape(phaseLens', iter[end] * sats) - modelRange)
   correction = pinv(designMat' * designMat) * designMat' * (reshape(phaseLens', iter[end] * sats) - modelRange)
   correction = pinv(designMat' * designMat) * designMat' * curr_error
return correction

end
function kinematicIter(aPriEst_c, rangeData, phaseData)

   epochs = size(phaseLens, 1)
   sats = size(phaseLens, 2)
   nbias = sats
   global designMat = hcat(zeros(epochs * sats*2, epochs * 4), zeros(epochs*sats*2, nbias))
   measurements = []
   model = []
   weightMatrix = zeros(epochs*sats*2, epochs*sats*2)
   rowIter = 1
   aPriBias = aPriEst_c[epochs*4+1:end]

   #Loop over all epochs
   for e in 1:epochs
      designMat[rowIter:rowIter + 2*sats-1, 4*e] = ones(2*sats)   #1s for time
      gpspos_e = gpsPositionData(e) #navigation constellation positions
      apri = aPriEst_c[4*e-3:4*e]   #per-epoch apriori estimation
      #Loop over measurements / satellites
      for s in 1:sats
         navsatpos = gpspos_e[s]
         posdiff = apri[1:3].- collect(navsatpos)
         r = norm(posdiff) #Geometric range
         designMat[rowIter, 4*e-3:4*e-1] = posdiff/r
         #designMat[rowIter, epochs*4+1:end] = -ones(32)
         designMat[rowIter, epochs*4+s] = -1
         designMat[rowIter+1, 4*e-3:4*e-1] = posdiff/r

         weightMatrix[rowIter, rowIter] = 1e3      #Weight for phase m.
         weightMatrix[rowIter+1, rowIter+1] = 1    #Weight for code m.

         append!(measurements, [phaseData[e, s], rangeData[e, s]])
         m_pr = r + apri[4]               #model for pseudo range
         m_ph = r + apri[4] - aPriBias[s] #model for phase
         append!(model, [m_ph, m_pr])

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

function pointPositionIteration(aPriEst, rangeData, navconPos)
   nNavsat = length(navconPos)
   aPriPos = Tuple(aPriEst[1:3])
   aPriRange = map(x -> norm(aPriPos .- navconPos[x])+aPriEst[4], 1:length(navconPos))
   designMat = NaviSimu.findPPDesignMatrix(aPriPos, navconPos)

   #Least squares estimation of position correction
   deltaPos = pinv(designMat' * designMat) * designMat' * (rangeData - aPriRange)
   return aPriEst .+ deltaPos
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


codeSig = reshape([], 0, 32)
phaseSig = reshape([], 0, 32)
iter = 1:96
timevec = (iter.-1) * 15 * 60
trueIss = [position(propagateKeplerOrbit(iss, t)) for t in timevec]

for ti in iter
   gpsPos = gpsPositionData(ti)

   curSignals = instaSignal(trueIss[ti], gpsPos, timevec[ti])
   global codeSig = vcat(codeSig, curSignals.code')
   global phaseSig = vcat(phaseSig, curSignals.phase')

end

Random.seed!(1)

#Convert code and phase measurements to ranges and phase lengths (convert to meters)
freq = 1575.42e6
waveLen = lightConst / freq
pseudoRanges = codeSig * lightConst
pseudoRanges += randn(size(pseudoRanges)) * 1   #1m normal errors
phaseLens = phaseSig * waveLen
phaseLens += randn(size(phaseLens)) *1e-2       #1mm normal errors

#Perform point position estimation and compare position error (magnitude seems to be 1e-9)
ppesti = [pointPosition(pseudoRanges[i,:], gpsPositionData(i); niter=5) for i in iter]
ppPosErrors = [norm(trueIss[i] .- ppesti[i][1:3]) for i in iter]


#Perform kinematic position estimation
apriPosTime = reshape(hcat(map(x-> collect(ppesti[x]), iter)...), 4*96, 1)
apriBias = (pseudoRanges[1,:] - phaseLens[1,:])
apriEst = vcat(apriPosTime, apriBias)

kinEstim = apriEst
corrections = []

@progress for i in 1:3
      #weights, singulariteit
   corr = kinematicIter(kinEstim, pseudoRanges,  phaseLens)
   append!(corrections, norm(corr))
   global kinEstim += corr
end
kinPosErrors = [norm(trueIss[i] .- kinEstim[4*i-3:4*i-1]) for i in iter]


acc1 = sum(ppPosErrors) / length(ppPosErrors)      #pdop * normal error
acc2 = sum(kinPosErrors) / length(kinPosErrors)
