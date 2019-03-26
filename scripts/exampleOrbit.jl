using Plots, ProgressMeter, LinearAlgebra
include("../src/NaviSimu.jl")
using Main.NaviSimu

function distances(refPoint, pointVector)
   distVec = map(x -> norm(pointVector[x] .- refPoint), 1:length(pointVector))
   return distVec
end

function pointPositionIteration(aPriEst, rangeData, navconPos)
   nNavsat = length(navconPos)
   aPriPos = Tuple(aPriEst[1:3])
   aPriRange = map(x -> norm(aPriPos .- navconPos[x])+aPriEst[4], 1:length(navconPos))
   designMat = NaviSimu.findPPDesignMatrix(aPriPos, navconPos)

   #Least squares estimation of position correction
   deltaPos = inv(designMat' * designMat) * designMat' * (rangeData - aPriRange)
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

#Calculate 'measured' ranges from satellite position and gps data. Neglect shadowing
truePos = position(iss)
ranges = distances(truePos, gps1)
ranges += randn(length(ranges)) * 5

esti = [0, 0, 0, 0]
estis = esti'
for i in 1:100
   global esti = pointPositionIteration(esti, ranges, gps1)
end

estPos = esti[1:3]
posError = truePos .- estPos
display(norm(posError))
