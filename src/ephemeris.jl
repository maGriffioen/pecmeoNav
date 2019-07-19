abstract type Ephemeris end

# Kepler Ephemeris data structure
struct KeplerEphemeris <: Ephemeris
    timeReferences::Array{<:Number}
    orbits:: Array{<:Orbit}
end

# Standard deviations for Noisy Kepler Ephemeris
struct KeplerEphemerisSD
    a_sigma::Number
    e_sigma::Number
    i_sigma::Number
    raan_sigma::Number
    aop_sigma::Number
    tanom_sigma::Number
end

# Search for an orbit index in a Kepler Epehemeris corresponding to a given time
function ephemerisIndexSearch(ephemeris::KeplerEphemeris, time)
    i = 0   # Iterator to select correct epheremis

    # Find last reference time smaller than time
    for t_ref in ephemeris.timeReferences
        if (time > t_ref)
            i+=1
        end
    end
    # println("i= ",i)
    # If time is before first reference, extrapolate back in time
    if (i==0)
        i =1
    end

    return i
end

if (NaviSimuAsModule)
    # Interpolate in KeplerEphemeris to find position at given time
    function NaviSimu.globalPosition(ephemeris::KeplerEphemeris, time)
        i = ephemerisIndexSearch(ephemeris::KeplerEphemeris, time)

        timeOffset = time - ephemeris.timeReferences[i]
        # println("dt= ", timeOffset)
        return globalPosition(ephemeris.orbits[i], timeOffset)
    end
else
    # Interpolate in KeplerEphemeris to find position at given time
    function globalPosition(ephemeris::KeplerEphemeris, time)
        i = ephemerisIndexSearch(ephemeris::KeplerEphemeris, time)

        timeOffset = time - ephemeris.timeReferences[i]
        # println("dt= ", timeOffset)
        return globalPosition(ephemeris.orbits[i], timeOffset)
    end
end

if (NaviSimuAsModule)
    # Obtain keplerorbit from KeplerEphemeris with elements AT GIVEN TIME
    function NaviSimu.KeplerOrbit(ephemeris::KeplerEphemeris, time)
        i = ephemerisIndexSearch(ephemeris::KeplerEphemeris, time)
        epoch = ephemeris.timeReferences[i]
        # Move orbit to an epoch of 0
        orbit = propagateKeplerOrbit(ephemeris.orbits[i], -epoch)

        return orbit
    end
else
    # Obtain keplerorbit from KeplerEphemeris with elements AT GIVEN TIME
    function KeplerOrbit(ephemeris::KeplerEphemeris, time)
        i = ephemerisIndexSearch(ephemeris::KeplerEphemeris, time)
        epoch = ephemeris.timeReferences[i]
        # Move orbit to an epoch of 0
        orbit = propagateKeplerOrbit(ephemeris.orbits[i], -epoch)

        return orbit
    end
end

# function NaviSimu.globalPosition(Array{KeplerEphemeris, 1}) = [NaviSimu.globalPosition()]

# Generate a perfect Kepler Ephemeris
function trueKeplerEphemeris(times::Array{<:Number}, orbit::Orbit)
    return KeplerEphemeris(times, [propagateKeplerOrbit(orbit, t) for t in times])
end

# Generate a disturbed Kepler Ephemeris from Kepler element Noise
function noisyKeplerEphemeris(times::Array{<:Number}, orbit::KeplerOrbit, noise::KeplerEphemerisSD)
    orbit_array = Array{KeplerOrbit, 1}(undef, 0)
    normError = randn(6)
    for t in times
        curr_orbit = propagateKeplerOrbit(orbit, t)
        noisyOrbit = KeplerOrbit(
            curr_orbit.a + noise.a_sigma * normError[1],
            curr_orbit.e + noise.e_sigma * normError[2],
            curr_orbit.i + noise.i_sigma * normError[3],
            curr_orbit.raan + noise.raan_sigma * normError[4],
            curr_orbit.aop + noise.aop_sigma * normError[5],
            curr_orbit.tanom + noise.tanom_sigma * normError[6],
            curr_orbit.cbody
        )
        push!(orbit_array, noisyOrbit)
    end

    return KeplerEphemeris(times, orbit_array)
end


# Generate a disturbed Kepler Ephemeris from state Noise
#TODO: MAKE THIS FUNCTION
