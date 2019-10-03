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

# Find transmitter position and time of transmission (due to light time effect)
# This is a higher accuracy version of transmitterFinder (using Double64 instead of Float64)
function transmitterFinder2(receptionTime::Number, receiverPos::Tup3d, transmitterOrbit::Orbit;
    maximumCorrection = 1e-20, maxIter = 10)
    receptionTime = Double64(receptionTime)
    receiverPos = Double64.(receiverPos)
    to = transmitterOrbit
    transmitterOrbit = KeplerOrbit(Double64(to.a), Double64(to.e), Double64(to.i),
        Double64(to.raan), Double64(to.aop), Double64(to.tanom), to.cbody)
    maximumCorrection = Double64(maximumCorrection)
    travelTime = df64"0.0"    #Assume instant signal as first estimate
    correction = df64"1"      #Force loop to start
    transTime = receptionTime   #Initialization
    transPos = (df64"0.0", df64"0.0", df64"0.0")  #Initialization
    distance = df64"0.0"

    nIter = 0   #Iteration counter
    # Iterate while no convergence of smaller than xx seconds is reached
    while (abs(correction) > maximumCorrection && nIter < maxIter)   #TODO: Tweak value for balance between stability and numerical accuracy
        # Calculate new transmitter position & time
        transTime = receptionTime - travelTime      #Time of transmission of signal
        transPos = globalPosition(transmitterOrbit, transTime)  #Position of transmitting satellite
        distance = norm(transPos .- receiverPos)    #Travel distance of signal
        travelTime = distance/lightConst            #Travel time of signal

        # Find applied correction for the decision if convergence has been reached
        correction = receptionTime - transTime - travelTime
        nIter +=1
    end

    return (transmissionTime = Float64(transTime), pos = Float64.(transPos), distance = Float64(distance), travelTime = Float64(travelTime),
        lastCorrection = Float64(correction), iterations = nIter)
end

function transmitterFinder(receptionTime::Number, receiverPos::Tup3d, transmitterOrbit::Orbit;
    maximumCorrection = 1e-15, maxIter = 10)
    travelTime = 0.0    #Assume instant signal as first estimate
    correction = 1      #Force loop to start
    transTime = receptionTime   #Initialization
    transPos = (0.0, 0.0, 0.0)  #Initialization
    distance = 0.0

    nIter = 0   #Iteration counter
    # Iterate while no convergence of smaller than xx seconds is reached
    while (abs(correction) > maximumCorrection && nIter < maxIter)   #TODO: Tweak value for balance between stability and numerical accuracy
        # Calculate new transmitter position & time
        transTime = receptionTime - travelTime      #Time of transmission of signal
        transPos = globalPosition(transmitterOrbit, transTime)  #Position of transmitting satellite
        distance = norm(transPos .- receiverPos)    #Travel distance of signal
        travelTime_new = distance/lightConst            #Travel time of signal

        # Find applied correction for the decision if convergence has been reached
        correction = travelTime - travelTime_new
        travelTime = travelTime_new
        nIter +=1
    end

    return (transmissionTime = transTime, pos = transPos, distance = distance, travelTime = travelTime,
        lastCorrection = correction, iterations = nIter)
end

# function transmitterFinder(receptionTime::Number, receiverPos::Tup3d, transmitterEpoch::Array{<:Ephemeris})
#
#     return transmitterFinder(receptionTime, receiverPos, transmitterOrbit)
# end

# Calculate 'fake' measurements without lighttime effect, noise, etc.
# Single epoch
function instantMeasurements(recPos::Tup3d, transPos::Tup3d, time::Number)
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

# Calculate 'fake' measurements without lighttime effect, noise, etc.
# Multiple epochs
function instantMeasurements(recPos::Tup3d, transPos::Array{<:Tup3d}, time::Number)
    nsats = length(transPos)
    codes = Array{Float64, 1}(undef, length(transPos))
    phases = Array{Float64, 1}(undef, length(transPos))
    avail = BitArray(undef, length(transPos))
    for isat in 1:nsats
        res = instantMeasurements(recPos, transPos[isat], time)
        codes[isat] = res.code
        phases[isat] = res.phase
        avail[isat] = res.los
    end
    return (code = codes, phase = phases, avail = avail)
end


function codeErrorVarianceFromNoise_TiberiusCoherent(CNR::Number, loopBandwidth::Number, earlyToLateSpacing::Number, prnChipLength::Number)
    sigma2 = (loopBandwidth * earlyToLateSpacing / (2*CNR))
    return sqrt(sigma2) * prnChipLength
end

function codeErrorVarianceFromNoise_TiberiusEarlyMinusLate(CNR::Number, loopBandwidth::Number, earlyToLateSpacing::Number, prnChipLength::Number, predetectionIntegrationInterval::Number)
    sigma2 = (loopBandwidth * earlyToLateSpacing / (2*CNR)) * (1+ (2/(predetectionIntegrationInterval*(2-earlyToLateSpacing)*CNR)))
    return sqrt(sigma2) * prnChipLength
end

function codeErrorVarianceFromNoise_TiberiusDot(CNR::Number, loopBandwidth::Number, earlyToLateSpacing::Number, prnChipLength::Number, predetectionIntegrationInterval::Number)
    sigma2 = (loopBandwidth * earlyToLateSpacing / (2*CNR)) * (1+ (2/(predetectionIntegrationInterval*CNR)))
    return sqrt(sigma2) * prnChipLength
end

function phaseErrorVarianceFromNoise_Tiberius(CNR::Number, loopBandwidth::Number, wavelength::Number, predetectionIntegrationInterval::Number)
    sigma2 = (loopBandwidth / CNR) * (1+ (1/(2*predetectionIntegrationInterval*CNR)))
    return (sqrt(sigma2)/ (2*pi)) * wavelength
end

function phaseErrorVarianceFromNoise_TiberiusSimplified(CNR::Number, loopBandwidth::Number, wavelength::Number)
    sigma2 = (loopBandwidth / CNR)
    return (sqrt(sigma2)/ (2*pi)) * wavelength
end

function phaseErrorVarianceFromNoise_GPSWorld(SNR::Number, sampleFreq::Number, wavelength::Number)
    sigma = wavelength/(pi * sqrt(2* SNR * 19200))
    return sigma
end


# Simulate the measurements at a single epoch
# Considers:
#   Lighttime effect
#   Signal noise due to CNR
#   Receiver clock error (from difference inertialTime and receiverTime inputs)
#
# INPUT:
#   inertialTime        Inertial (true) time when the measurements are taken
#   receiverTime        Epoch of the measurements, as stamped by receiver clock
#   receiverOrbit       Orbit of the recevier
#   navcon              Orbits of the navigation constallation
#   receiverSettings    RangeReceiverSettings object with info about the receiver
#   transmitterSettings Array of RangeTransmitterSettings object with info about the transmitters in the constellation
# INPUT (keyword arguments):
#   lighttimeEffect     True to apply lighttime effect, false for no lighttime effect
#   addNoise            True to apply random noise to measurements based on CNR
#   writeToFile         True to write to outputFile path
#   outputFile          Where to write the measurements to
#
# OUTPUT:
#   code                Code (pseudorange) meaasurements
#   phase               Phase measurements
#   avail               Availability of satellite for measurements
#   timeStamp           Timestamp (receiver time) of the measurement
#   codeSD              Standard Deviations of the code measurement noise
#   phaseSD             Standard Deviations of the phase measurement noise
function simulateMeasurements(inertialTime::Number, receiverTime::Number, receiverOrbit::Orbit, navcon::Orbit,
    receiverSettings::RangeReceiverSettings, transmitterSettings::Array{RangeTransmitterSettings,1};
    lighttimeEffect = true, addNoise = true, writeToFile = false, outputFile = "",
    addCarrierAmbiguity = true,
    lastAvailability = [], carrierAmbiguity = [])
    nsats = size(navcon)                        #Navigation constellation size
    codes = Array{Float64, 1}(undef, nsats)
    phases = Array{Float64, 1}(undef, nsats)
    codeSDs = Array{Float64, 1}(undef, nsats)
    phaseSDs = Array{Float64, 1}(undef, nsats)
    avail = BitArray(undef, nsats)

    freq = receiverSettings.frequency
    wavelength = lightConst / freq
    carrierAmbiguityRange = -10000:10000

    if lastAvailability == []
        lastAvailability = BitArray(undef, nsats)
        lastAvailability .= false
    end

    if addCarrierAmbiguity && carrierAmbiguity == []
        carrierAmbiguity = wavelength * rand(carrierAmbiguityRange, nsats)
    end

    #Loop over navigation satellites for this epoch
    for prn in 1:nsats
        receiverPosition = globalPosition(receiverOrbit, inertialTime)
        if lighttimeEffect
            # Transmitter position during transmission of received signal
            transmitter = transmitterFinder(inertialTime, receiverPosition, navcon[prn]; maxIter = 20)
        else
            # Transmitter position during signal reception (no light time effect)
            transmitter = (pos = globalPosition(navcon[prn], inertialTime),
            transmissionTime = inertialTime)
        end

        connection = true

        # Determine if in shadow of Earth or Moon
        # Check the line of sight both during reception and transmission. Assume no LOS when satellites are shadowed by body during either time.
        losEarly = hasLineOfSightEarthMoon(receiverPosition, transmitter.pos, transmitter.transmissionTime)
        losLate = hasLineOfSightEarthMoon(receiverPosition, transmitter.pos, inertialTime)
        los = (losEarly && losLate)

        # Check if satellites are connected and within line of sight
        available = los && connection


        avail[prn] = available
        if (available)
            if lighttimeEffect
                txTime = transmitter.transmissionTime
            else
                txTime = inertialTime
            end

            if (addNoise)
                #Check Link budget and use SNR & C/N0 to determine error SD -> generate errors from this
                txSettings = transmitterSettings[prn]
                snrCalc = snrCalculation(inertialTime, receiverTime, receiverPosition, transmitter,
                        receiverSettings, txSettings)
                cnr = snrCalc.cnr_straight
                codeSD = sqrt(((4 * 1) / (2*cnr)  )) * 293     # Simple model, allow replacement of model :)
                phaseSD =sqrt((receiverSettings.trackingLoopBandwidth / cnr) *
                    (1+ (1/ ( 2*cnr * receiverSettings.predetectionIntegrationInterval )))) * (wavelength/(2*pi))
                codeError = randn() * codeSD  # Error within receiver system, dependent on SNR
                phaseError = randn() * phaseSD  # Error within receiver system, dependent on C/N0
            else
                codeError = 0.0
                phaseError = 0.0
                codeSD = 0.0
                phaseSD = 0.0
            end
            codeSDs[prn] = codeSD
            phaseSDs[prn] = phaseSD

            if addCarrierAmbiguity
                if !lastAvailability[prn]
                    carrierAmbiguity[prn] = rand(carrierAmbiguityRange)
                end
                carrierAmbiguity_cur = carrierAmbiguity[prn]
            else
                carrierAmbiguity_cur = 0
            end

                # Neglecting transmitter clock error and relativistic effect -> will be largely recovered
            if lighttimeEffect# Please use this and apply the light time effect :)
                # Bad practice with an inaccurate method
                # codes[prn] = (receiverTime - transmitter.transmissionTime) * lightConst + codeError

                # Improved good method
                codes[prn] = (receiverTime - inertialTime) * lightConst + transmitter.distance + codeError

            else
                range = norm(globalPosition(receiverOrbit, inertialTime) .- transmitter.pos)
                # Unrealistically neglecting light time effect, which really shouldnt be neglected
                codes[prn] = (receiverTime - inertialTime) * lightConst + range + codeError
            end


            #Simulate Phase measurement
            if lighttimeEffect
                # Please use this and apply the light time effect :)
                # Bad practice method
                # phases[prn] = wavelength * (receiverSettings.zero_phase + freq*receiverTime -
                #     transmitterSettings[prn].zero_phase - freq * transmitter.transmissionTime) + phaseError

                # Slightly improved method
                # phases[prn] = transmitter.distance + wavelength * (receiverSettings.zero_phase + freq*(receiverTime- inertialTime) -
                #     transmitterSettings[prn].zero_phase) + phaseError

                # This is the final GOOD representation!
                phases[prn] = transmitter.distance + wavelength * (receiverSettings.zero_phase - transmitterSettings[prn].zero_phase) +
                    lightConst * (receiverTime- inertialTime) +phaseError + carrierAmbiguity_cur
            else
                # Unrealistically neglecting light time effect, which really shouldnt be neglected
                # phases[prn] = wavelength * (receiverSettings.zero_phase + freq*(receiverTime) -
                #     transmitterSettings[prn].zero_phase - freq * (inertialTime - range / lightConst)) + phaseError
                phases[prn] = range + wavelength * (receiverSettings.zero_phase + freq*(receiverTime - inertialTime) -
                    transmitterSettings[prn].zero_phase) + phaseError + carrierAmbiguity_cur


            end
                # Receiver Clock bias is taken into receiverTime.
                # Assuming no transmitter clock bias is present
                # (phi_u(t0) + f*(t-t0) + f*(clock bias u)) - (phi^s(t0) + f(t-tau-t0) + f*(clock bias s) )
        else
            # No measurement taken
            codes[prn] = NaN64
            phases[prn] = NaN64
            codeSDs[prn] = 0.0
            phaseSDs[prn] = 0.0
        end
    end
    return (code = codes, phase = phases, avail = avail, timeStamp = receiverTime,
    codeSD = codeSDs, phaseSD = phaseSDs, carrierAmbiguity = carrierAmbiguity)
end

# Simulate measurements at sequential times
# Considers:
#   Lighttime effect
#   Signal noise due to CNR
#   Receiver clock error (from difference inertialTime and receiverTime inputs)
#
# INPUT:
#   inertialTime        Inertial (true) times when the measurements are taken
#   receiverTime        Epochs of the measurements, as stamped by receiver clock
#   receiverOrbit       Orbit of the recevier
#   navcon              Orbits of the navigation constallation
#   receiverSettings    RangeReceiverSettings object with info about the receiver
#   transmitterSettings Array of RangeTransmitterSettings object with info about the transmitters in the constellation
# INPUT (keyword arguments):
#   lighttimeEffect     True to apply lighttime effect, false for no lighttime effect
#   addNoise            True to apply random noise to measurements based on CNR
#
# OUTPUT:
#   codeObs             Code (pseudorange) meaasurements
#   phaseObs            Phase measurements
#   availability        Availability of satellites for measurements, per epoch, per satellite
#   timeStamp           Timestamps (receiver time) of the measurements
#   codeSD              Standard Deviations of the code measurement noise, per epoch, per satellite
#   phaseSD             Standard Deviations of the phase measurement noise, per epoch, per satellite
function simulateMeasurements(inertialTime::Array{<:Number}, receiverTime::Array{<:Number}, receiverOrbit::Orbit, navcon::Orbit,
    receiverSettings::RangeReceiverSettings, transmitterSettings::Array{RangeTransmitterSettings,1}; lighttimeEffect = true, addNoise=true,
    addCarrierAmbiguity = true)
    #Raise error if number of measurements is not clear
    if (length(inertialTime) != length(receiverTime))
        error("simulateMeasurements: Number of true times not equal to number of receiver times")
    end

    nepochs = length(inertialTime)     # Number of epochs
    nsats = size(navcon) #Total number of satellites in the navigation constellation
    codeObs = Array{Float64}(undef, nepochs, nsats)    #Per-epoch, per-prn code observation
    phaseObs = Array{Float64}(undef, nepochs, nsats)   #Per-epoch, per-prn phase observation
    codeSD = Array{Float64}(undef, nepochs, nsats)
    phaseSD = Array{Float64}(undef, nepochs, nsats)
    availability = BitArray(undef, nepochs, nsats)     #Per-epoch, per-prn availability boolean

    availLast = BitArray(undef, nsats)
    availLast .= false
    carrierAmbiguity = zeros(Float64, nsats)

    for epoch = 1:nepochs
        msrmt = simulateMeasurements(inertialTime[epoch],
                    receiverTime[epoch],
                    receiverOrbit, navcon,
                    receiverSettings, transmitterSettings;
                    lighttimeEffect = lighttimeEffect,
                    addNoise = addNoise,
                    addCarrierAmbiguity = addCarrierAmbiguity,
                    lastAvailability = availLast,
                    carrierAmbiguity = carrierAmbiguity)
        codeObs[epoch, :] = msrmt.code
        phaseObs[epoch, :] = msrmt.phase
        codeSD[epoch, :] = msrmt.codeSD
        phaseSD[epoch, :] = msrmt.phaseSD
        availability[epoch, :] = msrmt.avail
        availLast = msrmt.avail
        carrierAmbiguity = msrmt.carrierAmbiguity
    end

    return (codeObs = codeObs, phaseObs = phaseObs, availability = availability, timeStamp = receiverTime, codeSD=codeSD, phaseSD=phaseSD)
end
