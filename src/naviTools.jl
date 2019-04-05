# Find the point position design matrix for a navigation solution
function findPPDesignMatrix(
    refPoint::Tuple{Number, Number, Number},
    constelData::Array{Tuple{Float64, Float64, Float64}, 1} )

    n_navsats = size(constelData, 1)       #Number of navigation satellites
    mat_A = zeros(n_navsats, 4)     #Create empty A matrix
    mat_A[:, 4] .= 1.0            #Set -1 for time positions

    #Loop over the navigation satellites in the constellation
    for isat in 1:n_navsats
        posnavsat = constelData[isat]               #Nav satellite position
        posdiff = collect(refPoint .- posnavsat)    #Position difference
        r  = norm(posdiff)                          #Total distance
        mat_A[isat, 1:3] = posdiff ./ r             #Write unit vector distances
    end
    # Return matrix of covariances
    return mat_A
end
findPPDesignMatrix( refPoint::Tuple{Number, Number, Number},
     constel::KeplerConstellation ) = findPPDesignMatrix(refPoint, position(constel))


 # Evaluate Covariance matrix of the Navigation system
 # For a reference point and a constellation (in given position)\
 function findNavCovmat(
     refPoint::Tuple{Number, Number, Number},
     constelData::Array{Tuple{Float64, Float64, Float64}, 1} )
     mat_A = findPPDesignMatrix(refPoint, constelData)
     return inv(mat_A' * mat_A)
 end
 findNavCovmat( refPoint::Tuple{Number, Number, Number},
      constel::KeplerConstellation ) = findNavCovmat(refPoint, position(constel))

# # Evaluate Covariance matrix of the Navigation system
# # For a reference point and a constellation (in given position)\
# function findNavCovmat(
#     refPoint::Tuple{Number, Number, Number},
#     constelData::Array{Tuple{Float64, Float64, Float64}, 1} )
#
#     n_navsats = size(constelData, 1)       #Number of navigation satellites
#     mat_A = zeros(n_navsats, 4)     #Create empty A matrix
#     mat_A[:, 4] .= -1.0;            #Set -1 for time positions
#
#     #Loop over the navigation satellites in the constellation
#     for isat in 1:n_navsats
#         posnavsat = constelData[isat]        #Nav satellite position
#         posdiff = collect(posnavsat .- refPoint)    #Position difference
#         r  = norm(posdiff)             #Total distance
#         mat_A[isat, 1:3] = posdiff ./ r             #Write unit vector distances
#     end
#     # Return matrix of covariances
#     return inv(mat_A' * mat_A)
# end
# findNavCovmat( refPoint::Tuple{Number, Number, Number},
#      constel::KeplerConstellation ) = findNavCovmat(refPoint, position(constel))



#Find the GDOP for reference point and constellation
#Uses Earth position at t=0 and Earth radius for shadowing
function findNavGDOP(
    refPoint::Tuple{Number, Number, Number},
    constelData::Array{Tuple{Float64, Float64, Float64}, 1} )
    los = hasLineOfSight(refPoint, constelData, bodyPosition(earth, 0), earth.radius)
    GDOP = 1e9   #If gdop cannot be calculated (too little satellites or bad geometry)
    try
        mat_Q = findNavCovmat(refPoint, constelData[los])
        GDOP = sqrt(sum(diag(mat_Q)))   #Calculate gdop
    catch
    end

    return GDOP
end
findNavGDOP(refPoint::Tuple{Number, Number, Number},
    constel::KeplerConstellation) = findNavGDOP(refPoint, position(constel))

function findNavPDOP(
    refPoint::Tuple{Number, Number, Number},
    constelData::Array{Tuple{Float64, Float64, Float64}, 1} )
    los = hasLineOfSight(refPoint, constelData, bodyPosition(earth, 0), earth.radius)
    PDOP = 1e9   #If gdop cannot be calculated (too little satellites or bad geometry)
    try
        global mat_Q = findNavCovmat(refPoint, constelData[los])
        PDOP = sqrt(sum(diag(mat_Q)[1:3]))   #Calculate gdop
    catch
    end

    return PDOP
end
findNavPDOP(refPoint::Tuple{Number, Number, Number},
    constel::KeplerConstellation) = findNavGDOP(refPoint, position(constel))

function hasLineOfSight(receiverLocation::Tuple{Float64, Float64, Float64},
    transmitterLocation::Tuple{Float64, Float64, Float64},
    bodyLocation::Tuple{Float64, Float64, Float64},
    bodyRadius::Float64)

    #Set up the vector from receiver to transmitter
    transDir = transmitterLocation .- receiverLocation
    r = norm(transDir)      #Distance between trans. and rec.
    transDir = transDir ./ r           #Unit vector to transmitter

    #Body
    bodyToRec = receiverLocation .- bodyLocation
    bodyToTrans = transmitterLocation .- bodyLocation
    signalDistance = norm(bodyToRec .- (dot(bodyToRec, transDir)).*transDir )

    isBodyInline = signalDistance < bodyRadius
    isBodyBetween = max(sqrt( sum(bodyToTrans.^2) - bodyRadius^2),
        sqrt( sum(bodyToRec.^2) - bodyRadius^2)) < r


    return !(isBodyInline && isBodyBetween)
end

hasLineOfSight(receiverLoc::Tuple{Float64, Float64, Float64},
    transmitterLoc::Array{Tuple{Float64, Float64, Float64}, 1},
    bodyLoc::Tuple{Float64, Float64, Float64},
    bodyRadius::Float64) =
    [hasLineOfSight(receiverLoc, transmitterLoc[satid], bodyLoc, bodyRadius)
        for satid in 1:size(transmitterLoc, 1)]

hasLineOfSightEarth(receiverLoc, transmitterLoc) =
    hasLineOfSight(receiverLoc, transmitterLoc, bodyPosition(earth.name, 0), earth.radius)


function pointPosition(ranges, navSats; niter = 10)
   esti = [0, 0, 0, 0]
   for i in 1:niter
      esti = pointPositionIteration(esti, ranges, navSats)
   end
   return Tuple(esti)
end

function sequentialPointPosition(epochTimes, navigationConstellation::KeplerConstellation, pseudoRanges, availability; niter = 10)
    n_epochs = length(epochTimes)
    navCon = navigationConstellation
    return [pointPosition(pseudoRanges[epoch, availability[epoch, :]],
                position(propagateKeplerOrbit(gpsKepler, timevec[epoch]))[availability[epoch, :]]; niter=niter) for epoch in 1:n_epochs]
end



function pointPositionIteration(aPriEst, rangeData, navconPos)
   nNavsat = length(navconPos)
   aPriPos = Tuple(aPriEst[1:3])
   aPriRange = map(x -> norm(aPriPos .- navconPos[x])+aPriEst[4], 1:length(navconPos))
   designMat = findPPDesignMatrix(aPriPos, navconPos)

   #Least squares estimation of position correction
   deltaPos = pinv(designMat' * designMat) * designMat' * (rangeData - aPriRange)
   return aPriEst .+ deltaPos
end
