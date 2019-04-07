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
     constel::KeplerConstellation ) = findPPDesignMatrix(refPoint, globalPosition(constel))


 # Evaluate Covariance matrix of the Navigation system
 # For a reference point and a constellation (in given position)\
 function findNavCovmat(
     refPoint::Tuple{Number, Number, Number},
     constelData::Array{Tuple{Float64, Float64, Float64}, 1} )
     mat_A = findPPDesignMatrix(refPoint, constelData)
     return inv(mat_A' * mat_A)
 end
 findNavCovmat( refPoint::Tuple{Number, Number, Number},
      constel::KeplerConstellation ) = findNavCovmat(refPoint, globalPosition(constel))

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
    constel::KeplerConstellation) = findNavGDOP(refPoint, globalPosition(constel))

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
    constel::KeplerConstellation) = findNavGDOP(refPoint, globalPosition(constel))

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

# Perform a single point position
function pointPosition(ranges, navSats; maxIter = 10, correctionLimit = 1e-8)
   estimation = [0.0, 0.0, 0.0, 0.0]
   niter = 0
   correction = correctionLimit * 10

   while (niter < maxIter && correction > correctionLimit)
      estimation_new = pointPositionIteration(estimation, ranges, navSats)
      correction = maximum(abs.(estimation_new.-estimation))
      estimation = estimation_new
      niter += 1
   end
   return Tuple(estimation)
end

# Perform a single point position estimation iteration
function pointPositionIteration(aPriEst, rangeData, navconPos)
   nNavsat = length(navconPos)
   aPriPos = Tuple(aPriEst[1:3])
   aPriRange = map(x -> norm(aPriPos .- navconPos[x])+aPriEst[4], 1:length(navconPos))
   designMat = findPPDesignMatrix(aPriPos, navconPos)

   #Least squares estimation of position correction
   deltaPos = pinv(designMat' * designMat) * designMat' * (rangeData - aPriRange)
   return aPriEst .+ deltaPos
end

# Perform batch of point position estimations for a single satellite on various epochs
function sequentialPointPosition(epochTimes, navigationConstellation::KeplerConstellation,
    pseudoRanges, availability; maxIter = 10, correctionLimit = 1e-8)

    n_epochs = length(epochTimes)
    navCon = navigationConstellation
    return [pointPosition(pseudoRanges[epoch, availability[epoch, :]],
                globalPosition(propagateKeplerOrbit(navCon, epochTimes[epoch]))[availability[epoch, :]];
                    maxIter=maxIter, correctionLimit=correctionLimit)
                for epoch in 1:n_epochs]
end

# Perform the entire kinematic estmiation
function kinematicEstimation(navCon::KeplerConstellation, epochTimes,
   rangeData, phaseData, availability; maxIter::Number = 5, correctionLimit::Number = 1e-8, maxIter_pp = 20,
   codeWeight = 1, phaseWeight = 1e6)

   n_epochs = length(epochTimes)    #Number of epochs
   # Create apriori position and clock error estimation through point positioning
   global ppesti = sequentialPointPosition(epochTimes, navCon, rangeData, availability; maxIter=maxIter_pp, correctionLimit = 1e-3)
   global apriPosTime = vcat(map(x-> collect(ppesti[x]), 1:n_epochs)...)  #apriori position and time are those from point positioning

   kinEstim = apriPosTime  #Dont add bias estimation -> The kinematic iterator adds those dynamically when needed
   iter = 0
   correction = correctionLimit * 10   #Ensure the iterations are started

   # Perform iterate the kinematic estimation
   while (iter < maxIter && correction > correctionLimit)
      # Perform single iteration
      kinIter = kinematicIter(navCon, epochTimes, kinEstim,
         rangeData, phaseData, availability;
         codeWeight = codeWeight, phaseWeight = phaseWeight)
      kinEstim = kinIter.estimation                   #Apply estimation
      correction = maximum(abs.(kinIter.correction))  #Check the maximum correction
      iter += 1
   end

   # Restructure the ouput
   positionTimeEstim = [Tuple(kinEstim[epoch*4-3 : epoch*4])
      for epoch in 1:n_epochs]

   biasEstim = kinEstim[n_epochs*4+1 : end]


   return (positionTimeEstimation = positionTimeEstim, biasEstimation = biasEstim,
      iterations = iter, lastCorrection = correction)
end


# Perform a single iteration of kinematic estimation of satellite position, clock errors and bias
function kinematicIter(navCon::KeplerConstellation, epochTimes,
   aPriEst_c::Array{Float64, 1}, rangeData, phaseData, availability;
   codeWeight = 1.0, phaseWeight = 1e6)

   avail = availability
   n_epochs = size(avail, 1)
   n_sats = size(avail, 2)
   n_measurements = sum(avail) *2
   prns = collect(1:size(avail)[2])

   #Allocate bias numbers to the various satellite arcs
   biasNumber = zeros(Int16, n_epochs, maximum(prns))
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

   n_bias = maximum(biasNumber)  #Number of biases to be estimated

   # Pre-allocate vectors and matrices
   designMat = zeros(Float64,  n_measurements, n_epochs *4 + n_bias)
   measurements = zeros(Float64, n_measurements)
   model = zeros(Float64, n_measurements)
   weightMatrix = Matrix{Float64}(I, n_measurements, n_measurements)

   # Dynamically grow apriori estimation vector
   aPriEst = aPriEst_c
   while length(aPriEst) < 4*n_epochs
      #Risky zero-apriori estimations added to fill apriori vector.
      aPriEst.append(0.0)
   end

   aPrioriLength = 4*n_epochs + n_bias    #'Right' apriori vector length
   while length(aPriEst) < aPrioriLength
      #Add (estimated) bias to grow apriori vector size
      iAddedBias = length(aPriEst) - 4*n_epochs +1
      dataFilter = biasNumber.==iAddedBias
      aPriBias = sum((rangeData[dataFilter] - phaseData[dataFilter] )) / sum(dataFilter)
      append!(aPriEst, aPriBias)
   end

   rowIter = 1
   #Loop over all epochs
   for epoch in 1:n_epochs
      gpspos_e = globalPosition(propagateKeplerOrbit(navCon, epochTimes[epoch])) #navigation constellation positions
      apri_e = aPriEst_c[4*epoch-3:4*epoch]   #apriori pos en t estimation for this epoch

      #Loop over available satellites
      for prn in prns[avail[epoch, :]]
         navsatpos = gpspos_e[prn]
         posdiff = apri_e[1:3].- collect(navsatpos)
         r = norm(posdiff) #Geometric range
         sats_uvec = posdiff/r

         #Add phase measurement to design matrix and weight matrix
         designMat[rowIter, 4*epoch-3:4*epoch] = vcat(sats_uvec, 1.0)
         designMat[rowIter, n_epochs*4 + biasNumber[epoch, prn]] = -1   #bias
         weightMatrix[rowIter, rowIter] = phaseWeight      #Weight for phase measurement = sigma^2

         #Add code measurement to design matrix and weight matrix
         designMat[rowIter+1, 4*epoch-3:4*epoch] = vcat(sats_uvec, 1.0)
         weightMatrix[rowIter+1, rowIter+1] = codeWeight    #Weight for code measurement = sigma^2

         # Add code and phase measurements to measurements vector
         measurements[rowIter] = phaseData[epoch, prn]
         measurements[rowIter+1] = rangeData[epoch, prn]

         # Find mode code and phase data, and add to model vector
         model_phase = r + apri_e[4] - aPriEst[n_epochs*4 + biasNumber[epoch, prn]] #model for phase
         model[rowIter] = model_phase
         model_pseudorange = r + apri_e[4]               #model for pseudo range
         model[rowIter+1] = model_pseudorange

         #Add 2 (1x phase + 1x code) to the row iteration number
         rowIter += 2
      end
   end

   # Current error between measurements and model
   curr_error = measurements-model
   # Weighted least squares correction to estimation vector
   # correction = inv(designMat' *weightMatrix* designMat) * designMat' *weightMatrix* curr_error
   correction = (designMat' *weightMatrix* designMat) \ (designMat' *weightMatrix* curr_error)
   # ff = svd(designMat' *weightMatrix* designMat)
   # correction = (ff.V * inv(diagm(0=> ff.S)) * ff.U') * (designMat' *weightMatrix* curr_error)

   return (correction = correction, estimation = aPriEst+correction)
end
