function findSignalTravelTime(receiver::KeplerOrbit, transmitter::KeplerOrbit, time::Number)
    recState   = state(receiver, time)
    recPos = recState[1:3]
    return findSignalTravelTime(recPos, transmitter, time)
end

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

    return signalDistance / lightConst
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
        phaseRe = phaseRe0 + freq * time                       #+ clock error stuff
        phaseTr = phaseTr0 + freq * (time - signalTravelTime) #+ clock error stuff

        phase = phaseRe - phaseTr #+ N + error
    else
        signalTravelTime = NaN64
        phase = NaN64
    end
    return (code = signalTravelTime, phase = phase, los = los)
end

function instaSignal(recPos::Tup3d, transPos::Array{Tup3d}, time::Number)
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
