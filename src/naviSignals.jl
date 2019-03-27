function findSignalTravelTime(receiver::KeplerOrbit, transmitter::KeplerOrbit, time::Number)
    recState   = state(receiver, time)
    transState = state(transmitter, time)
    signalDistance = norm(recState[1:3] - transState[1:3])

    corr = 1.0
    i = 0
    while corr > 1e-6
        i += 1
        signalTravelTime = signalDistance / lightConst
        transState_new = state(transmitter, time - signalTravelTime)
        corr = norm(transState_new[1:3] - transState[1:3])
        signalDistance = norm(recState[1:3] - transState[1:3])
    end

    return signalDistance / lightConst
end


function instaSignal(recPos::Tup3d, transPos::Tup3d, time::Number)
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

    return (code = signalTravelTime, phase = phase)
end

function instaSignal(recPos::Tup3d, transPos::Array{Tup3d}, time::Number)
    nsats = length(transPos)
    codes = zeros(nsats)
    phases = zeros(nsats)

    for isat in 1:nsats
        res = instaSignal(recPos, transPos[isat], time)
        codes[isat] = res.code
        phases[isat] = res.phase
    end
    return (code = codes, phase = phases)
end
