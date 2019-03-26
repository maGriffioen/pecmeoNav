
# Object for Kepler Orbits
struct KeplerOrbit
    a::Float64      #Semi-major axis
    e::Float64      #Eccentricity
    i::Float64      #Inclination
    raan::Float64   #Right ascension of the ascending node
    aop::Float64    #Argument of perriapsis
    tanom::Float64  #True anomaly
    cbody::Body     #Central body of orbit
end
Base.copy(ko::KeplerOrbit)=
    KeplerOrbit(ko.a, ko.e, ko.i, ko.raan, ko.aop, ko.tanom, ko.cbody)

# Object for Cartesian states
struct CartesianState
    x::Float64      #X-position
    y::Float64      #Y-position
    z::Float64      #Z-position
    vx::Float64     #X-velocity
    vy::Float64     #Y-velocity
    vz::Float64     #Z-velocity
end
Base.copy(cs::CartesianState) =
    CartesianState(cs.x, cs.y, cs.z, cs.vx, cs.vy, cs.vz)

# Object for Cartesian orbits (=cartesian state + gravitational parameter)
struct CartesianOrbit
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    cbody::Body     #Central body of orbit
end
# Create Cartesian orbit directly from cartesian state
CartesianOrbit(cs::CartesianState, cbody) =
    CartesanOrbit(cs.x, cs.y, cs.z, cs.vx, cs.vy, cs.vz, cbody)
CartesianOrbit(kepOrbit::KeplerOrbit) = keplerToCartesian(kepOrbit::KeplerOrbit)
Base.copy(co::CartesianOrbit) =
    CartesianOrbit(co.x, co.y, co.z, co.vx, co.vy, co.vz, co.cbody)


# Object for a constellation of kepler orbits
struct KeplerConstellation
    orbits::Array{KeplerOrbit,1}
end
KeplerConstellation() = KeplerConstellation([])
KeplerConstellation(ko::KeplerOrbit) = KeplerConstellation([ko])
Base.copy(kc::KeplerConstellation) =
    KeplerConstellation([copy(orbit) for orbit in kc.orbits])
Base.getindex(constel::KeplerConstellation, i::Int) = constel.orbits[i]
Base.push!(constel::KeplerConstellation, orbit::KeplerOrbit) = push!(constel.orbits, orbit)
Base.size(constel::KeplerConstellation) = length(constel.orbits)

# Transform Kepler orbit to a Carthesian orbit type
function keplerToCartesian(keplerorbit::KeplerOrbit)
    ko = keplerorbit    #simplify notation
    mu = ko.cbody.gravitationalParameter

    #Find in-plane rectangular coordinates
    r = ko.a * (1- ko.e^2)/(1+ ko.e*cos(ko.tanom))
    pos2d = [r * cos(ko.tanom); r * sin(ko.tanom)]

    #Set up transformation matrix to go from in-plane 2d to 3d coordinates
    l1 = cos(ko.raan) * cos(ko.aop) - sin(ko.raan) * sin(ko.aop) * cos(ko.i)
    l2 = -cos(ko.raan) * sin(ko.aop) - sin(ko.raan) * cos(ko.aop) * cos(ko.i)
    m1 = sin(ko.raan) * cos(ko.aop) + cos(ko.raan) * sin(ko.aop) * cos(ko.i)
    m2 = -sin(ko.raan) * sin(ko.aop) + cos(ko.raan) * cos(ko.aop) * cos(ko.i)
    n1 = sin(ko.aop) * sin(ko.i)
    n2 = cos(ko.aop) * sin(ko.i)
    transforationMatrix = [l1 l2; m1 m2; n1 n2]

    #Find carthesian position
    pos3d = transforationMatrix * pos2d

    #Find angular momentum
    h = sqrt(mu * ko.a * (1- ko.e^2))

    #Find cathesian velocity vectors
    vel3d = (mu / h) * transforationMatrix * [-sin(ko.tanom); (ko.e+cos(ko.tanom))];

    CartesianOrbit(pos3d[1], pos3d[2], pos3d[3], vel3d[1], vel3d[2], vel3d[3], ko.cbody)
end

Base.position(cs::CartesianState) = (cs.x, cs.y, cs.z)
Base.position(co::CartesianOrbit) = (co.x, co.y, co.z)
Base.position(ko::KeplerOrbit) = position(keplerToCartesian(ko))
Base.position(kc::KeplerConstellation) = [position(orbit) for orbit in kc.orbits]

# Find orbital period of Kepler orbit
function orbitalPeriod(keplerorbit::KeplerOrbit)
    ko = keplerorbit
    return 2 * pi * sqrt((ko.a^3) / ko.cbody.gravitationalParameter)
end

# Find the eccentric anomaly for an orbit from true anomaly and eccentricity
function findEccentricAnomaly(trueAnomaly::Number, eccentricity::Number)
    #Equation 6.35 Astrodynamics reader
    eccentricAnomaly = 2 * atan(
        tan(trueAnomaly/2.0) *
        sqrt((1 - eccentricity)/(1 + eccentricity))
    )
    return eccentricAnomaly
end
findEccentricAnomaly(ko::KeplerOrbit) = findEccentricAnomaly(ko.tanom, ko.e)

# Find the mean anomaly for an orbit from true anomaly and eccentricity
function findMeanAnomaly(trueAnomaly::Number, eccentricity::Number)
    eccentricAnomaly = findEccentricAnomaly(trueAnomaly, eccentricity)
    #Equation 6.36-3 Astrodynamics reader
    meanAnomaly = eccentricAnomaly - eccentricity * sin( eccentricAnomaly )
    return meanAnomaly
end
findMeanAnomaly(ko::KeplerOrbit) = findMeanAnomaly(ko.tanom, ko.e)

function meanToTrueAnomaly(meanAnomaly::Number, eccentricity::Number,
    tolerance::Float64)

    eccAnomaly = meanAnomaly
    delta = 1

    if (eccentricity < 0.1)
        #Simple iterative method
        #6.37-5? Astrodynamics reader
        while(abs(delta) > tolerance)
            eccAnomalyNew = meanAnomaly + eccentricity * sin(eccAnomaly)
            delta =eccAnomalyNew - eccAnomaly
            eccAnomaly = eccAnomalyNew
        end

    elseif(eccentricity < 1.0)
        #Newton-Raphson method
        #6.37-4? Astrodynamics reader
        while(abs(delta) > tolerance)
            eccAnomalyNew = eccAnomaly -
            (( eccAnomaly - eccentricity * sin(eccentricAnomaly) - meanAnomaly ) /
            ( 1- eccentricity * cos(eccentricAnomaly) ))
            delta = eccAnomalyNew - eccAnomaly
            eccAnomaly = eccAnomalyNew
        end

    else
        #Not solvable (I think) with non-elliptical orbits
        error("Cannot find true anomaly with hyperbolic orbit")
    end

    #Eccentric anomaly to true anomaly conversion
    #Equation 6.35 Astrodynamics reader
    trueAnomaly = 2 * atan(
    sqrt((1 + eccentricity) / (1.0 - eccentricity)) *
    tan( eccAnomaly /2.0 )
    )

    return trueAnomaly
end
meanToTrueAnomaly(meanAnomaly::Number, eccentricity::Number) =
    meanToTrueAnomaly(meanAnomaly, eccentricity, 1e-7)

# Find the true anomaly of a Kepler orbit after a timeIncrement
function propagateKeplerOrbit(ko::KeplerOrbit, timeIncrement::Number)
    # Find initial mean anomaly of the orbit
    meanAnomaly = findMeanAnomaly(ko)

    # Find the mean motion of the orbit and progress it with timeIncrement
    meanMotion = sqrt(ko.cbody.gravitationalParameter / (ko.a^3))
    meanAnomaly += meanMotion * timeIncrement

    # Convert mean anomaly back to true anomaly
    newTrueAnomaly = meanToTrueAnomaly(meanAnomaly, ko.e)

    return KeplerOrbit(ko.a, ko.e, ko.i, ko.raan, ko.aop, newTrueAnomaly, ko.cbody)
end

# Propagate an entire constellation by looping over individual or bits
function propagateKeplerOrbit(constellation::KeplerConstellation, timeIncrement::Number)
    newConstellation = KeplerConstellation()
    for kOrbit in constellation.orbits
        push!( newConstellation, propagateKeplerOrbit(kOrbit, timeIncrement) )
    end
    return newConstellation
end



# Create a pecmeo constellation
#   Polar - Polar - Equatorial
function createCircPecmeo( radius::Float64, n_satellites::Tuple{Int, Int, Int}, cbody::Body,
    satSpacing::Tuple{Float64, Float64, Float64};
    initialOrbitShift::Tuple{Float64, Float64, Float64}= (0.0, 0.0, 0.0),
    equatorialRotation::Float64 = 0.0,
    inclination::Float64 = 0.0 )

    constellation = KeplerConstellation()
    inclinations = (pi/2 + inclination, pi/2, inclination)
    for i_orient in 1:3
        for i_sat in 1:n_satellites[i_orient]
            #Calculate orbital elements
            i = inclinations[i_orient]                              #inclination
            raan = equatorialRotation + (i_orient == 2 ? pi/2 : 0)  #right asc.
            aop = initialOrbitShift[i_orient]                       #arg. of periapsis
            tanom = satSpacing[i_orient] * (i_sat-1)     #true anomaly
            #Add Kepler orbit to the constellation
            push!( constellation,
                KeplerOrbit(radius, 0.0, i,
                raan, aop, tanom,
                cbody)
            )
        end
    end

    return constellation
end

createCircPecmeo( radius::Float64, n_satellites::Tuple{Int, Int, Int}, cbody::Body;
    initialOrbitShift::Tuple{Float64, Float64, Float64}= (0.0, 0.0, 0.0),
    equatorialRotation::Float64 = 0.0, inclination::Float64=0.0 )=
    createCircPecmeo(radius, n_satellites, cbody, 2*pi ./ n_satellites;
    initialOrbitShift = initialOrbitShift,
    equatorialRotation = equatorialRotation, inclination = inclination)
