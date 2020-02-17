struct RangeTransmitterSettings
    zero_phase::Float64     #Radians theta(t0), where t_0 is the local Clock
    transmitPower::Number   #Watt. Unused for receivers
    frequency:: Number      #Hz. Ensure to match receiver and transmitter freq.
    antennaGainFunction::Function #Antenna Gain f(theta, phi)   in db
    pointingFunction::Function  #Pointing vector of antenna f(pos_s, t)
    pointingErrorSD::Number     #Standard deviation of pointing error in radians

end

struct RangeReceiverSettings
    zero_phase::Float64     #Radians theta(t0), where t_0 is the local Clock
    receptionNoiseTemperature::Number   #Temperature of receiver noise
    receptionBandwidth::Number  #Bandwidth of receiver
    trackingLoopBandwidth::Number
    predetectionIntegrationInterval::Number
    frequency:: Number      #Hz. Ensure to match receiver and transmitter freq.
    antennaGainFunction::Function #Antenna Gain f(theta, phi)  in db
    pointingFunction::Function  #Pointing vector of antenna f(pos_s, t)
    pointingErrorSD::Number     #Standard deviation of pointing error in radians
end

# Calculate the linkbudget and find SNR and C/N0
function snrCalculation(inertialTime, receiverTime, receiverPosition, transmitter,
        receiverSettings, transmitterSettings)

        # Calculate antenna gains from pointing and signal vectors
        relativePositionRec = receiverPosition .- transmitter.pos
        signalDistance = norm(relativePositionRec)
        pointingVectorTX = transmitterSettings.pointingFunction( transmitter.pos, transmitter.transmissionTime )
        pointingVectorRec = receiverSettings.pointingFunction( receiverPosition, receiverTime )

        pointingoffsetTX = acos(dot( relativePositionRec, pointingVectorTX ) / (norm(relativePositionRec) * norm(pointingVectorTX)))
        gainTX = transmitterSettings.antennaGainFunction(pointingoffsetTX, 0.0)

        pointingoffsetRec = acos(dot( .-relativePositionRec, pointingVectorRec ) / (norm(.-relativePositionRec) * norm(pointingVectorRec)))
        gainRec = receiverSettings.antennaGainFunction( pointingoffsetRec, 0.0 )

        # Calculate some fundamental properties
        freq = receiverSettings.frequency
        wavelength = lightConst / freq

        # EIRP and free space loss
        EIRP = lin2log( transmitterSettings.transmitPower ) + gainTX  # [dB] effective receiver power
        fsl = lin2log( ( wavelength/( 4*pi*signalDistance ))^2 )    # [db] free space loss
        attenuation = 0     # [dB] hardware noise, polarization offsets, etc.
        receivedPower = EIRP + fsl - attenuation        # [dB] signal Power

        # Noise level
        recTemp = receiverSettings.receptionNoiseTemperature
        recBw = receiverSettings.receptionBandwidth
        noise = lin2log( recTemp * recBw * bolzmanConst )

        #SNR
        snr = receivedPower - noise # [dB] signal to noise ratio
        snr_straight = log2lin( snr ) # [-] signal to noise (straight) ratio
        cnr_straight = snr_straight * recBw # [dB-Hz]
        cnr = lin2log( cnr_straight )


        # return (snr = snr, cnr = cnr, EIRP = EIRP, fsl = fsl, att = attenuation, receivedPower = receivedPower, noise = noise, wl = wavelength)
        return (snr = snr, cnr = cnr, cnr_straight = cnr_straight)
end
