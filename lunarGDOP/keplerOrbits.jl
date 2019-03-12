struct KeplerOrbit
    a::Float64
    e::Float64
    i::Float64
    raan::Float64   #Right ascension of the ascending node
    aop::Float64    #Argument of perriapsis
    tanom::Float64  #True anomaly
    mu::Float64
end

struct CartesianState
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
end

struct CartesianOrbit
    x::Float64
    y::Float64
    z::Float64
    vx::Float64
    vy::Float64
    vz::Float64
    mu::Float64
end

#Create Cartesian orbit directly from cartesian state
CartesianOrbit(cs::CartesianState, mu) =
    CartesanOrbit(cs.x, cs.y, cs.z, cs.vx, cs.vy, cs.vz, mu)

#Transform Kepler orbit to a Carthesian orbit type
function keplerToCartesian(keplerorbit::KeplerOrbit)
    ko = keplerorbit    #simplify notation

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
    h = sqrt(ko.mu * ko.a * (1- ko.e^2))

    #Find cathesian velocity vectors
    vel3d = (ko.mu / h) * transforationMatrix * [-sin(ko.tanom); (ko.e+cos(ko.tanom))];

    CartesianOrbit(pos3d[1], pos3d[2], pos3d[3], vel3d[1], vel3d[2], vel3d[3], ko.mu)
end

#Find the eccentric anomaly for an orbit from true anomaly and eccentricity
function findEccentricAnomaly(trueAnomaly::Number, eccentricity::Number)
    #Equation 6.35 Astrodynamics reader
    eccentricAnomaly = 2 * atan(
        tan(trueAnomaly/2.0) *
        sqrt((1 - eccentricity)/(1 + eccentricity))
    )
    return eccentricAnomaly
end
findEccentricAnomaly(ko::KeplerOrbit) = findEccentricAnomaly(ko.tanom, ko.e)

#Find the mean anomaly for an orbit from true anomaly and eccentricity
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

#Find the true anomaly of a Kepler orbit after a timeIncrement
function propagateKeplerOrbit(ko::KeplerOrbit, timeIncrement::Number)
    #Find initial mean anomaly of the orbit
    meanAnomaly = findMeanAnomaly(ko)

    #Find the mean motion of the orbit and progress it with timeIncrement
    meanMotion = sqrt(ko.mu / (ko.a^3))
    meanAnomaly += meanMotion * timeIncrement

    #Convert mean anomaly back to true anomaly
    newTrueAnomaly = meanToTrueAnomaly(meanAnomaly, ko.e)

    return KeplerOrbit(ko.a, ko.e, ko.i, ko.raan, ko.aop, newTrueAnomaly, ko.mu)
end
