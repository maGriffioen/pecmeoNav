# Obsolete light time effect function
function findSignalTravelTime(receiver::KeplerOrbit, transmitter::KeplerOrbit, time::Number)
    recState   = state(receiver, time)
    recPos = recState[1:3]
    return findSignalTravelTime(recPos, transmitter, time)
end
# Obsolete light time effect function
function findSignalTravelTime(recPos::Tup3d, transmitter::KeplerOrbit, time::Number)
    transState = state(transmitter, time)
    signalDistance = norm(recPos - transState[1:3])

    corr = 1.0
    i = 0
    while corr > 1e-6
        i += 1
        signalTravelTime = signalDistance / lightConst
        transState_new = state(transmitter, time - signalTravelTime)
        # Corrected receiver position
        corr = norm(transState_new[1:3] - transState[1:3])
        transState = transState_new
        signalDistance = norm(recPos .- transState[1:3])
    end

    return (signalTravelTime = signalDistance / lightConst, transmitterPosition = transState[1:3])
end

# Find transmitter position and true time during transmission (Light time effect)
function transmitterFinder(receptionTime::Number, receiverPos::Tup3d, transmitterOrbit::Orbit)
    travelTime = 0.0    #Assume instant signal as first estimate
    correction = 1      #Force loop to start
    transTime = receptionTime   #Initialization
    transPos = (0.0, 0.0, 0.0)  #Initialization

    nIter = 0   #Iteration counter
    # Iterate while no convergence of smaller than xx seconds is reached
    while (abs(correction) > 1e-14 && nIter < 10)   #TODO: Tweak value for balance between stability and numerical accuracy
        # Calculate new transmitter position & time
        transTime = receptionTime - travelTime      #Time of transmission of signal
        transPos = globalPosition(transmitterOrbit, transTime)  #Position of transmitting satellite
        distance = norm(transPos .- receiverPos)    #Travel distance of signal
        travelTime = distance/lightConst            #Travel time of signal

        # Find applied correction for the decision if convergence has been reached
        correction = receptionTime - transTime - travelTime
        nIter +=1
        # println(nIter," ", correction)
    end
    # println(nIter)

    return (transmissionTime = transTime, pos = transPos, travelTime = travelTime,
        lastCorrection = correction, iterations = nIter)
end

function transmitterFinder(receptionTime::Number, receiverPos::Tup3d, transmitterEpoch::Array{<:Ephemeris})

    return transmitterFinder(receptionTime, receiverPos, transmitterOrbit)
end


function instaSignal(recPos::Tup3d, transPos::Tup3d, time::Number)
    los = hasLineOfSightEarthMoon(recPos, transPos, time)
    if los
        #Neglect doppler shifting and relativistic effects
        freq = 1575.42e6
        waveLen = lightConst / freq
        distance = norm(recPos .- transPos)
        signalTravelTime = distance / lightConst

        #Neglect clock errors on receiver and transmitter satellite

        #phase at t0 is 0 for receiver and transmitter
        phaseRe0 = 0
        phaseTr0 = 1e7
        # phaseRe = phaseRe0 + freq*time
        # phaseTr = phaseTr0 + freq * (time - signalTravelTime)
        phaseRe = phaseRe0
        phaseTr = phaseTr0 -freq*signalTravelTime
        # Note that phase = phaseRe - phaseTR = phaseRe0 - phaseTr0 + freq*time - freq* (time-signalTravelTime)
        #   = phaseRe0 - phaseTr - freq * signalTravelTime
        #   This eliminates the time variable, the magnitude of which grows linearly with time
        #   Hence, its truncation error grows with time, and thus the numerical error in calculating
        #   Observations grows with time.

        phase = phaseRe - phaseTr #+ N + error
    else
        signalTravelTime = NaN64
        phase = NaN64
    end
    return (code = signalTravelTime, phase = phase, los = los)
end

function instaSignal(recPos::Tup3d, transPos::Array{<:Tup3d}, time::Number)
    nsats = length(transPos)
    codes = Array{Float64, 1}(undef, length(transPos))
    phases = Array{Float64, 1}(undef, length(transPos))
    avail = BitArray(undef, length(transPos))
    for isat in 1:nsats
        res = instaSignal(recPos, transPos[isat], time)
        codes[isat] = res.code
        phases[isat] = res.phase
        avail[isat] = res.los
    end
    return (code = codes, phase = phases, avail = avail)
end
