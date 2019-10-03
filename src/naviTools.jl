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
    refPoint::Tup3d,
    constelData::Array{<:Tup3d, 1}; checkLineOfSightEarth = true, checkLOSMoon = false, time = 0 )
    # Check which satellites are within range
    if checkLineOfSightEarth && checkLOSMoon
        los = hasLineOfSightEarthMoon(refPoint,constelData, time)
    elseif checkLineOfSightEarth
        los = hasLineOfSight(refPoint, constelData, bodyPosition(earth, 0), earth.radius)
    elseif checkLOSMoon
        los = hasLineOfSightMoon(refPoint, constelData, time)
    else
        los = [true for i in 1:length(constelData)]
    end
    GDOP = 1e9   #If gdop cannot be calculated (too little satellites or bad geometry)
    try
        mat_Q = findNavCovmat(refPoint, constelData[los])
        GDOP = sqrt(sum(diag(mat_Q)))   #Calculate gdop
    catch
    end

    return GDOP
end
findNavGDOP(refPoint::Tup3d,
    constel::KeplerConstellation) = findNavGDOP(refPoint, globalPosition(constel))

function findNavTDOP(
    refPoint::Tup3d,
    constelData::Array{<:Tup3d, 1}; checkLineOfSightEarth = true, checkLOSMoon = false, time = 0 )
    # Check which satellites are within range
    if checkLineOfSightEarth && checkLOSMoon
        los = hasLineOfSightEarthMoon(refPoint,constelData, time)
    elseif checkLineOfSightEarth
        los = hasLineOfSight(refPoint, constelData, bodyPosition(earth, 0), earth.radius)
    elseif checkLOSMoon
        los = hasLineOfSightMoon(refPoint, constelData, time)
    else
        los = [true for i in 1:length(constelData)]
    end
    TDOP = 1e9   #If gdop cannot be calculated (too little satellites or bad geometry)
    try
        mat_Q = findNavCovmat(refPoint, constelData[los])
        TDOP = sqrt(mat_Q[4, 4])   #Calculate gdop
    catch
    end

    return TDOP
end
findNavTDOP(refPoint::Tup3d,
    constel::KeplerConstellation) = findNavTDOP(refPoint, globalPosition(constel))


function findNavPDOP(
    refPoint::Tup3d,
    constelData::Array{<:Tup3d, 1}; checkLineOfSightEarth = true, checkLOSMoon = false, time = 0 )
    # Check which satellites are within range
    if checkLineOfSightEarth && checkLOSMoon
        los = hasLineOfSightEarthMoon(refPoint,constelData, time)
    elseif checkLineOfSightEarth
        los = hasLineOfSight(refPoint, constelData, bodyPosition(earth, 0), earth.radius)
    elseif checkLOSMoon
        los = hasLineOfSightMoon(refPoint, constelData, time)
    else
        los = [true for i in 1:length(constelData)]
    end
    PDOP = 1e9   #If gdop cannot be calculated (too little satellites or bad geometry)
    try
        mat_Q = findNavCovmat(refPoint, constelData[los])
        PDOP = sqrt(sum(diag(mat_Q)[1:3]))   #Calculate gdop
    catch
    end

    return PDOP
end
findNavPDOP(refPoint::Tuple{Number, Number, Number},
    constel::KeplerConstellation) = findNavGDOP(refPoint, globalPosition(constel))

function hasLineOfSight(receiverLocation::Tup3d,
    transmitterLocation::Tup3d,
    bodyLocation::Tup3d,
    bodyRadius::Number)

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

hasLineOfSight(receiverLoc::Tup3d,
    transmitterLoc::Array{<:Tup3d, 1},
    bodyLoc::Tup3d,
    bodyRadius::Number) =
    [hasLineOfSight(receiverLoc, transmitterLoc[satid], bodyLoc, bodyRadius)
        for satid in 1:size(transmitterLoc, 1)]

hasLineOfSightEarth(receiverLoc, transmitterLoc) =
    hasLineOfSight(receiverLoc, transmitterLoc, bodyPosition(earth.name, 0), earth.radius)
hasLineOfSightMoon(receiverLoc, transmitterLoc, time) =
    hasLineOfSight(receiverLoc, transmitterLoc, bodyPosition(moon.name, time), moon.radius)
hasLineOfSightEarthMoon(receiverLoc, transmitterLoc, time) =
    hasLineOfSightEarth(receiverLoc, transmitterLoc) .* hasLineOfSightMoon(receiverLoc, transmitterLoc, time)

# Perform a single point position
function pointPosition(ranges, navSats, epochTime; maxIter = 10, correctionLimit = 1e-8,
    aprioriEstimation = [0.0, 0.0, 0.0, 0.0], lighttimeCorrection = true)

    estimation = aprioriEstimation
    nObservations = length(ranges)
    niter = 0

    if nObservations > 3
        sufficientObervables = true
            #TODO:: What happens with n < 4? neglect those?

        correction = correctionLimit * 10

        while (niter < maxIter && correction > correctionLimit)
          estimation_new = pointPositionIteration(estimation, ranges, navSats,
            epochTime; lighttimeCorrection = lighttimeCorrection)
          correction = maximum(abs.(estimation_new.-estimation))
          estimation = estimation_new
          niter += 1
        end
    else
        sufficientObervables = false
        correction = 0.0
    end

    return (estimation = Tuple(estimation), success = sufficientObervables, lastCorrection = correction, iterations = niter)
end

# Perform a single point position estimation iteration
# TODO:: ADD lightime correction warning
function pointPositionIteration(aPriEst, rangeData, navconPos::Array{<:Tup3d, 1}, epochTime; lighttimeCorrection = true)

    # println("Point position iteration based on static constellation coordinates")
    nNavsat = length(navconPos)
    aPriPos = Tuple(aPriEst[1:3])
    aPriRange = map(x -> norm(aPriPos .- navconPos[x])+aPriEst[4], 1:length(navconPos))
    designMat = findPPDesignMatrix(aPriPos, navconPos)

    #Least squares estimation of position correction
    deltaPos = pinv(designMat' * designMat) * designMat' * (rangeData - aPriRange)
    return aPriEst .+ deltaPos
end

# Perform a single point position estimation iteration based on constellation Ephemeris
function pointPositionIteration(aPriEst, rangeData, navigationEphemeris::Array{<:Ephemeris}, epochTime; lighttimeCorrection = true)

    nNavsat = length(navigationEphemeris)

    # println("Point position iteration based on constellation ephemeres")
    measurementTime = epochTime
    if (lighttimeCorrection)
        # Find positions of navigation satellites INCLUDING lighttime effect
        # Ephemeris are those selected from measurementTime
        # navConOrbits = [KeplerOrbit(navigationEphemeris[prn], measurementTime) for prn in 1:nNavsat]
        navconPos = [transmitterFinder(measurementTime-(aPriEst[4]/lightConst), Tuple(aPriEst[1:3]), KeplerOrbit(navigationEphemeris[prn], measurementTime)).pos for prn in 1:nNavsat]
    else
        # Find positions of navigation satellite at time of signal reception
        # e.g. not taking into account the light time effect
        navconPos = [globalPosition(navigationEphemeris[prn], measurementTime-(aPriEst[4]/lightConst)) for prn in 1:nNavsat]
    end

    if length(rangeData) == 0
        println("Zero-measurements case")
    end

    aPriPos = Tuple(aPriEst[1:3])
    aPriRange = map(x -> norm(aPriPos .- navconPos[x])+aPriEst[4], 1:length(navconPos))
    designMat = findPPDesignMatrix(aPriPos, navconPos)

    #Least squares estimation of position correction
    deltaPos = pinv(designMat' * designMat) * designMat' * (rangeData - aPriRange)
    return aPriEst .+ deltaPos
end

# Perform batch of point position estimations for a single satellite on various epochs
function sequentialPointPosition(epochTimes, navigationEphemeris::Array{<:Ephemeris},
    pseudoRanges, availability; aprioriEstimations = [0.0], maxIter = 10, correctionLimit = 1e-8,
    lighttimeCorrection = true)


    n_epochs = length(epochTimes)

    # Check validity of apriori estimations, and create zeros to reset when needed
    if (aprioriEstimations == [0.0])
        aprioriEstimations = zeros(n_epochs, 4)

    elseif (size(aprioriEstimations, 1) != n_epochs
        || size(aprioriEstimations, 2) != 4)

        print("sequentialPointPosition: Input apriori estimations of size ",
        size(aprioriEstimations), ", ", (n_epochs, 4), " expected.
        Zero vector used instead.")

        aprioriEstimations = zeros(n_epochs, 4)
    end


    # Ability to bipass lighttime correction with old point position estimation
    # It is recommended to keep ltc = true
    ltc = true


    if !ltc
        ppResult = [pointPosition(pseudoRanges[epoch, availability[epoch, :]],
                    [globalPosition(navigationEphemeris[i], epochTimes[epoch]) for i in 1:length(navigationEphemeris)][availability[epoch, :]],
                    epochTimes[epoch];
                    maxIter=maxIter, correctionLimit=correctionLimit,
                    aprioriEstimation = aprioriEstimations[epoch, :],
                    lighttimeCorrection = lighttimeCorrection)
                    for epoch in 1:n_epochs]
    else
        ppResult =  [pointPosition(pseudoRanges[epoch, availability[epoch, :]],
                    navigationEphemeris[availability[epoch, :]],
                    epochTimes[epoch];
                    maxIter=maxIter, correctionLimit=correctionLimit,
                    aprioriEstimation = aprioriEstimations[epoch, :],
                    lighttimeCorrection = lighttimeCorrection)
                    for epoch in 1:n_epochs]
    end
    ppEstimation = [i.estimation for i in ppResult]
    estimationSuccess = [i.success for i in ppResult]
    lastCorrections = [i.lastCorrection for i in ppResult]
    iterations = [i.iterations for i in ppResult]

    return (estimation = ppEstimation, resultValidity = estimationSuccess,
        lastCorrections = lastCorrections, iterations=iterations)
end


function findMeasurementArchs(availabilityMatrix)
    # arch=0
    inArch = false
    archRanges = []
    arch_startepoch = NaN
    nEpochs = size(availabilityMatrix)[1]
    for e in 1:nEpochs
        nObs_curr = sum(availabilityMatrix[e,:])
        if nObs_curr < 4 && inArch
            inArch = false
            push!( archRanges, arch_startepoch: e-1 )
        elseif nObs_curr < 4 && !inArch

        elseif inArch

        elseif !inArch
            inArch = true
            # arch += 1
            arch_startepoch = e
        end
    end
    if inArch
        push!( archRanges, arch_startepoch: nEpochs)
    end

    return archRanges
end

# Perform the entire kinematic estmiation
function kinematicEstimation(navigationEphemeris::Array{<:Ephemeris}, epochTimes,
    rangeData, phaseData, availability;
    ppApriori = [0], maxIter_kin::Int = 5, correctionLimit_kin::Number = 1e-8, maxIter_pp::Int = 20, correctionLimit_pp::Number = 1.0,
    codeWeight = 1, phaseWeight = 1e6, lighttimeCorrection = true, solver=2, verbose=false)

    n_epochs = length(epochTimes)    #Number of epochs
    prns = collect(1:size(availability)[2])

    # Create apriori position and clock error estimation through point positioning
    ppesti = sequentialPointPosition(epochTimes, navigationEphemeris, rangeData, availability;
    aprioriEstimations = ppApriori, maxIter=maxIter_pp, correctionLimit = correctionLimit_pp,
    lighttimeCorrection = lighttimeCorrection).estimation
    # apriPosTime = vcat(map(x-> collect(ppesti[x]), 1:n_epochs)...)  #apriori position and time are those from point positioning
    # println(typeof(ppesti))

    # Allocate bias numbers
    biasNumber = zeros(Int16, n_epochs, maximum(prns))
    for epoch in 1:n_epochs
        for prn in prns[availability[epoch, :]]
            if epoch > 1 && biasNumber[epoch-1, prn] > 0
               #Continue using old bias if a number has been assigned in previous epoch
               biasNumber[epoch, prn] = biasNumber[epoch-1, prn]
            else
               #Create new bias Number
               biasNumber[epoch, prn] = maximum(biasNumber) + 1
            end
        end
    end

    # Select function for solver iterator
    solverFuns = [kinematicIter, kinematicIterSmartSolve]
    kinIterFun = solverFuns[solver]

    # Output which function is used when verbose = true
    if (verbose)
        println("Kinematic iterator function: ", kinIterFun)
    end

    archs = findMeasurementArchs(availability)
    kinematic_xyztResults = copy(ppesti)
    kinematic_biasResults = zeros(maximum(biasNumber))
    kinematic_processed = BitArray(undef, length(epochTimes))
    kinematic_processed .= false
    iterations = Array{Int16, 1}(undef, 0)
    corrections = Array{Float64, 1}(undef, 0)

    variance_xyzt = zeros(n_epochs, 4)

    # Split into visibility archs and solev on a per-arch basis
    for archEpochs in archs
        if (verbose)
            println("Processing arch with epochs: ", archEpochs)
        end
        availability_arch = availability[archEpochs, :]
        nepochs_arch = length(archEpochs)
        range_arch = rangeData[archEpochs, :]
        phase_arch = phaseData[archEpochs, :]
        epochTimes_arch = epochTimes[archEpochs]
        apriPosTime_arch = vcat(map(x-> collect(ppesti[x]), archEpochs)...)
        biasNumber_arch = biasNumber[archEpochs, :]
        biasRange_arch = minimum(biasNumber_arch[biasNumber_arch.>0]):maximum(biasNumber_arch)
        for i in 1:length(biasNumber_arch)
            if biasNumber_arch[i] > 0
                biasNumber_arch[i] -= (first(biasRange_arch)-1)
            end
        end

        kinematic_processed[archEpochs] .= true

       kinEstim = apriPosTime_arch  #Dont add bias estimation -> The kinematic iterator adds those dynamically when needed

       iter = 0
       correction = correctionLimit_kin * 10   #Ensure the iterations are started

       # Perform iterate the kinematic estimation
       while (iter < maxIter_kin && correction > correctionLimit_kin)
          # Perform single iteration
          kinIter = kinIterFun(navigationEphemeris, epochTimes_arch, kinEstim,
             range_arch, phase_arch, availability_arch, biasNumber_arch;
             codeWeight = codeWeight, phaseWeight = phaseWeight,
             lighttimeCorrection = lighttimeCorrection)
          kinEstim = copy(kinIter.estimation)
          correction = maximum(abs.(kinIter.correction))  #Check the maximum correction
          iter += 1
       end

       # Restructure the ouput
       kinematic_xyztResults[archEpochs] = [Tuple(kinEstim[epoch*4-3 : epoch*4])
          for epoch in archEpochs.-(first(archEpochs)-1)]

       kinematic_biasResults[biasRange_arch] = kinEstim[nepochs_arch*4+1 : end]


       kincov = kinematic_covariance(navigationEphemeris, epochTimes_arch, kinEstim,
          range_arch, phase_arch, availability_arch, biasNumber_arch;
          codeWeight = codeWeight, phaseWeight = phaseWeight,
          lighttimeCorrection = lighttimeCorrection)

      var_xxs = diag(kincov.cov_xx)
      for epoch in archEpochs
          localEpochNumber = epoch - archEpochs[1] +1
          idx = 1+4*(localEpochNumber-1)
          variance_xyzt[epoch, :] = [var_xxs[idx+i] for i in 0:3 ]
      end

      append!(iterations, iter)
      append!(corrections, correction)
    end

   return (positionTimeEstimation = kinematic_xyztResults, biasEstimation = kinematic_biasResults,
      iterations = iterations, lastCorrection = corrections, kinTimes = epochTimes, archs = archs,
      boolarchs = kinematic_processed, variance_xyzt = variance_xyzt)
end

# Perform a single iteration of kinematic estimation of satellite position, clock errors and bias
function kinematicIter(navigationEphemeris::Array{<:Ephemeris}, epochTimes,
   aPriEst_c::Array{Float64, 1}, rangeData, phaseData, availability, biasNumbers;
   codeWeight = 1.0, phaseWeight = 1e6, lighttimeCorrection = true)

   avail = availability
   n_epochs = size(avail, 1)
   n_sats = size(avail, 2)
   n_measurements = sum(avail) *2
   prns = collect(1:size(avail)[2])
   n_bias = maximum(biasNumbers)  #Number of biases to be estimated

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
      dataFilter = biasNumbers.==iAddedBias
      aPriBias = sum((rangeData[dataFilter] - phaseData[dataFilter] )) / sum(dataFilter)
      append!(aPriEst, aPriBias)
   end

   rowIter = 1
   #Loop over all epochs
   for epoch in 1:n_epochs
      apri_e = aPriEst_c[4*epoch-3:4*epoch]   #apriori pos and time estimation for this epoch
      if (lighttimeCorrection)
          navconPos_e = [transmitterFinder(epochTimes[epoch], Tuple(apri_e[1:3]), KeplerOrbit(navigationEphemeris[prn], epochTimes[epoch])).pos for prn in 1:n_sats]
      else
          navconPos_e = [globalPosition(navigationEphemeris[i], epochTimes[epoch]) for i in 1:length(navigationEphemeris)] #navigation constellation positions
      end

      #Loop over available satellites
      for prn in prns[avail[epoch, :]]
         navsatpos = navconPos_e[prn]
         posdiff = apri_e[1:3].- collect(navsatpos)
         r = norm(posdiff) #Geometric range
         sats_uvec = posdiff/r

         #Add phase measurement to design matrix and weight matrix
         designMat[rowIter, 4*epoch-3:4*epoch] = vcat(sats_uvec, 1.0)
         designMat[rowIter, n_epochs*4 + biasNumbers[epoch, prn]] = -1   #bias
         weightMatrix[rowIter, rowIter] = phaseWeight      #Weight for phase measurement = sigma^2

         #Add code measurement to design matrix and weight matrix
         designMat[rowIter+1, 4*epoch-3:4*epoch] = vcat(sats_uvec, 1.0)
         weightMatrix[rowIter+1, rowIter+1] = codeWeight    #Weight for code measurement = sigma^2

         # Add code and phase measurements to measurements vector
         measurements[rowIter] = phaseData[epoch, prn]
         measurements[rowIter+1] = rangeData[epoch, prn]

         # Find mode code and phase data, and add to model vector
         model_phase = r + apri_e[4] - aPriEst[n_epochs*4 + biasNumbers[epoch, prn]] #model for phase
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

function kinematicIterSmartSolve(navigationEphemeris::Array{<:Ephemeris}, epochTimes,
   aPriEst_c::Array{Float64, 1}, rangeData, phaseData, availability, biasNumbers;
   codeWeight = 1.0, phaseWeight = 1e6, lighttimeCorrection = true)

   avail = availability
   n_epochs = size(avail, 1)
   n_sats = size(avail, 2)
   n_measurements = sum(avail) *2
   prns = collect(1:size(avail)[2])


    n_bias = maximum(biasNumbers)  #Number of biases to be estimated

    # Pre-allocate vectors and matrices
    Nxx_inv = spzeros(n_epochs*4, n_epochs*4)
    Nxb = spzeros(n_epochs*4, n_bias)
    Nbb = spzeros(n_bias, n_bias)
    nx = zeros(n_epochs*4, 1)
    nb = zeros(n_bias, 1)

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
        dataFilter = biasNumbers.==iAddedBias
        aPriBias = sum((rangeData[dataFilter] - phaseData[dataFilter] )) / sum(dataFilter)
        append!(aPriEst, aPriBias)
    end

    rowIter = 1
    #Loop over all epochs
    for epoch in 1:n_epochs
        apri_e = aPriEst_c[4*epoch-3:4*epoch]   #apriori pos and time estimation for this epoch
        if (lighttimeCorrection)
          navconPos_e = [transmitterFinder(epochTimes[epoch]-(apri_e[4]/lightConst), Tuple(apri_e[1:3]), KeplerOrbit(navigationEphemeris[prn], epochTimes[epoch])).pos for prn in 1:n_sats]
        else
          navconPos_e = [globalPosition(navigationEphemeris[i], epochTimes[epoch]-(apri_e[4]/lightConst)) for i in 1:length(navigationEphemeris)] #navigation constellation positions
        end

        eRange = 1+((epoch-1)*4):epoch*4
        Nxx_tmp = zeros(Float64, 4, 4)    #Temporary storage of local Nxx
        #Loop over available satellites = each 2 measurements
        for prn in prns[avail[epoch, :]]
            navsatpos = navconPos_e[prn]
            posdiff = apri_e[1:3].- collect(navsatpos)
            r = norm(posdiff) #Geometric range

            # Sats_uvec is [i; j; k; 1.0] = [r/|x|; r/|y|; r/|z|; dt/|dt|]
            sats_uvec = vcat(posdiff/r, 1.0)

            # Find mode code and phase data, and add to model vector
            biasNumber_current = biasNumbers[epoch, prn]
            model_phase = r + apri_e[4] - aPriEst[n_epochs*4 + biasNumber_current] #model for phase
            model_pseudorange = r + apri_e[4]               #model for pseudo range

            # Build onto error matrices
            nx[eRange] += sats_uvec .*codeWeight .* (rangeData[epoch, prn] - model_pseudorange)
            nx[eRange] += sats_uvec .*phaseWeight .* (phaseData[epoch, prn] - model_phase)
            nb[biasNumber_current] += -1 * phaseWeight .*(phaseData[epoch, prn] - model_phase)

            # Build onto normal equation matrices
            Nxx_tmp += (sats_uvec * sats_uvec') * (phaseWeight + codeWeight)
            Nbb[biasNumber_current, biasNumber_current] += phaseWeight
            Nxb[eRange, biasNumber_current] += -(sats_uvec) * phaseWeight
        end
        #Inverse the lcoal Nxx for faster inversion
        Nxx_inv[eRange, eRange] = inv(Nxx_tmp) #Inverse of local N_XX 4x4 matrix
    end

    Nbx = Nxb'

    # Bias and Position updates
    deltaB = (Nbb - Nxb'*Nxx_inv*Nxb) \ (nb - Nxb'*Nxx_inv*nx)
    deltaX = Nxx_inv*(nx - Nxb*deltaB)
    correction = vec(vcat(deltaX, deltaB))

    return (correction = correction, estimation = aPriEst+correction)
end


function kinematic_covariance(navigationEphemeris::Array{<:Ephemeris}, epochTimes,
   aPriEst_c::Array{Float64, 1}, rangeData, phaseData, availability, biasNumbers;
   codeWeight = 1.0, phaseWeight = 1e6, lighttimeCorrection = true)

   avail = availability
   n_epochs = size(avail, 1)
   n_sats = size(avail, 2)
   n_measurements = sum(avail) *2
   prns = collect(1:size(avail)[2])


    n_bias = maximum(biasNumbers)  #Number of biases to be estimated

    # Pre-allocate vectors and matrices
    Nxx_inv = zeros(n_epochs*4, n_epochs*4)
    Nxb = zeros(n_epochs*4, n_bias)
    Nbb = zeros(n_bias, n_bias)
    nx = zeros(n_epochs*4, 1)
    nb = zeros(n_bias, 1)

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
        dataFilter = biasNumbers.==iAddedBias
        aPriBias = sum((rangeData[dataFilter] - phaseData[dataFilter] )) / sum(dataFilter)
        append!(aPriEst, aPriBias)
    end

    rowIter = 1
    #Loop over all epochs
    for epoch in 1:n_epochs
        apri_e = aPriEst_c[4*epoch-3:4*epoch]   #apriori pos and time estimation for this epoch
        if (lighttimeCorrection)
          navconPos_e = [transmitterFinder(epochTimes[epoch]-(apri_e[4]/lightConst), Tuple(apri_e[1:3]), KeplerOrbit(navigationEphemeris[prn], epochTimes[epoch])).pos for prn in 1:n_sats]
        else
          navconPos_e = [globalPosition(navigationEphemeris[i], epochTimes[epoch]-(apri_e[4]/lightConst)) for i in 1:length(navigationEphemeris)] #navigation constellation positions
        end

        eRange = 1+((epoch-1)*4):epoch*4
        Nxx_tmp = zeros(Float64, 4, 4)    #Temporary storage of local Nxx
        #Loop over available satellites = each 2 measurements
        for prn in prns[avail[epoch, :]]
            navsatpos = navconPos_e[prn]
            posdiff = apri_e[1:3].- collect(navsatpos)
            r = norm(posdiff) #Geometric range

            # Sats_uvec is [i; j; k; 1.0] = [r/|x|; r/|y|; r/|z|; dt/|dt|]
            sats_uvec = vcat(posdiff/r, 1.0)

            # Find mode code and phase data, and add to model vector
            biasNumber_current = biasNumbers[epoch, prn]
            model_phase = r + apri_e[4] - aPriEst[n_epochs*4 + biasNumber_current] #model for phase
            model_pseudorange = r + apri_e[4]               #model for pseudo range

            # Build onto error matrices
            nx[eRange] += sats_uvec .*codeWeight .* (rangeData[epoch, prn] - model_pseudorange)
            nx[eRange] += sats_uvec .*phaseWeight .* (phaseData[epoch, prn] - model_phase)
            nb[biasNumber_current] += -1 * phaseWeight .*(phaseData[epoch, prn] - model_phase)

            # Build onto normal equation matrices
            Nxx_tmp += (sats_uvec * sats_uvec') * (phaseWeight + codeWeight)
            Nbb[biasNumber_current, biasNumber_current] += phaseWeight
            Nxb[eRange, biasNumber_current] += -(sats_uvec) * phaseWeight
        end
        #Inverse the lcoal Nxx for faster inversion
        Nxx_inv[eRange, eRange] = inv(Nxx_tmp) #Inverse of local N_XX 4x4 matrix
    end

    Nbx = Nxb'

    # Bias and Position updates
    cov_bb = inv(Nbb - Nxb'*Nxx_inv*Nxb)
    cov_xx = Nxx_inv + (Nxx_inv*Nxb) * cov_bb * ((Nxx_inv*Nxb)')

    return (cov_bb = cov_bb, cov_xx = cov_xx)
end
