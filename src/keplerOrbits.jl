# Define body type
# For the definition of several bodies with data, see bodies.jl
abstract type Orbit end
abstract type SingleOrbit <: Orbit end
abstract type Constellation <: Orbit end

struct Body
    name::String
    gravitationalParameter::Number
    radius::Number
    stateFunction::Function
end
Body(; name = "", gravitationalParameter, radius, stateFunction) =
    Body(name, gravitationalParameter, radius, stateFunction)

# Object for Kepler Orbits
struct KeplerOrbit <: SingleOrbit
    a::Number      #Semi-major axis
    e::Number      #Eccentricity
    i::Number      #Inclination
    raan::Number   #Right ascension of the ascending node
    aop::Number    #Argument of perriapsis
    tanom::Number  #True anomaly
    cbody::Body     #Central body of orbit
end
Base.copy(ko::KeplerOrbit)=
    KeplerOrbit(ko.a, ko.e, ko.i, ko.raan, ko.aop, ko.tanom, ko.cbody)

# Object for Cartesian orbits (=cartesian state + gravitational parameter)
struct CartesianOrbit <: SingleOrbit
    x::Number
    y::Number
    z::Number
    vx::Number
    vy::Number
    vz::Number
    cbody::Body     #Central body of orbit
end
# Create Cartesian orbit directly from cartesian state
CartesianOrbit(kepOrbit::KeplerOrbit) = keplerToCartesian(kepOrbit::KeplerOrbit)
KeplerOrbit(cartOrbit::CartesianOrbit) = cartesianToKepler(cartOrbit::CartesianOrbit)
Base.copy(co::CartesianOrbit) =
    CartesianOrbit(co.x, co.y, co.z, co.vx, co.vy, co.vz, co.cbody)


# Object for a constellation of kepler orbits
struct KeplerConstellation <: Constellation
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

# Transform Carthesian orbit into Kepler orbit type
function cartesianToKepler(cartesianorbit::CartesianOrbit)
    # Based on Noomen slides
    co = cartesianorbit
    r_vec = [co.x, co.y, co.z]
    v_vec = [co.vx, co.vy, co.vz]
    mu = co.cbody.gravitationalParameter
    r = norm(r_vec)     # orbital radius
    v = norm(v_vec)     # orbital velocity
    h = cross(r_vec, v_vec) # rotational energy
    N = cross([0, 0, 1], h)

    a = 1/((2/r) - (v*v/mu))    # semi-major axis
    e_vec = (cross(v_vec, h) / mu) - (r_vec / r)
    e = norm(e_vec)             # eccentricity

    i = acos(h[3] / norm(h))    # inclination

    Nx = N[1]
    Ny = N[2]
    Nxy = sqrt(Nx*Nx + Ny*Ny)
    N_unit = N / norm(N)
    e_unit = e_vec / e
    r_unit = r_vec / r
    # Final Kepler elements
    rightAscensionAscendingNode = atan(Ny / Nxy, Nx/Nxy)
    argumentOfPerigee = sign(dot(cross(N_unit, e_vec), h)) * acos(dot(e_unit, N_unit))
    trueAnomaly = sign(dot(cross(e_vec, r_vec), h)) * acos(dot(r_unit, e_unit))

    return KeplerOrbit(a, e, i, rightAscensionAscendingNode, argumentOfPerigee, trueAnomaly, co.cbody)
end

#Legacy positions
Base.position(co::CartesianOrbit) = (co.x, co.y, co.z)
Base.position(ko::KeplerOrbit) = position(keplerToCartesian(ko))
Base.position(kc::KeplerConstellation) = [position(orbit) for orbit in kc.orbits]

localPosition(co::CartesianOrbit; time::Number=0.0) = (co.x, co.y, co.z)
localPosition(ko::KeplerOrbit; time::Number=0.0) = localPosition(keplerToCartesian(ko); time=time)
localPosition(kc::KeplerConstellation; time::Number=0.0) = [localPosition(orbit; time=time) for orbit in kc.orbits]
localPosition(orbit::Orbit, time::Number) = localPosition(propagateKeplerOrbit(orbit, time); time=time)

globalPosition(co::CartesianOrbit; time::Number=0.0) = (co.x, co.y, co.z) .+ bodyPosition(co.cbody, time)
globalPosition(ko::KeplerOrbit; time::Number=0.0) = globalPosition(keplerToCartesian(ko); time=time)
globalPosition(kc::KeplerConstellation; time::Number=0.0) = [globalPosition(orbit; time=time) for orbit in kc.orbits]
globalPosition(orbit::Orbit, time::Number) = globalPosition(propagateKeplerOrbit(orbit, time); time=time)

localState(co::CartesianOrbit; time::Number=0.0) = [co.x, co.y, co.z, co.vx, co.vy, co.vz]
localState(ko::KeplerOrbit; time::Number=0.0) = localState(keplerToCartesian(ko); time=time)
localState(kc::KeplerConstellation; time::Number=0.0) = [localState(orbit; time=time) for orbit in kc.orbits]
localState(orbit::Orbit, time::Number) = localState(propagateKeplerOrbit(orbit, time); time=time)

globalState(co::CartesianOrbit; time::Number=0.0) = [co.x, co.y, co.z, co.vx, co.vy, co.vz] .+ globalState(co.cbody, time)
globalState(ko::KeplerOrbit; time::Number=0.0) = globalState(keplerToCartesian(ko); time=time)
globalState(kc::KeplerConstellation; time::Number=0.0) = [globalState(orbit; time=time) for orbit in kc.orbits]
globalState(orbit::Orbit, time::Number) = globalState(propagateKeplerOrbit(orbit, time); time=time)


#Legacy states
state(ko::KeplerOrbit; time::Number=0.0) = state(keplerToCartesian(ko))
state(co::CartesianOrbit; time::Number=0.0) = [co.x, co.y, co.z, co.vx, co.vy, co.vz]
state(kc::KeplerConstellation; time::Number=0.0) = [state(orbit) for orbit in kc.orbits]
state(orbit, time::Number) = state(propagateKeplerOrbit(orbit, time))

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
    tolerance::Number)
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
    elseif (eccentricity < 1.0)
        #Newton-Raphson method
        #6.37-4? Astrodynamics reader
        while(abs(delta) > tolerance)
            eccAnomalyNew = eccAnomaly -
            (( eccAnomaly - eccentricity * sin(eccAnomaly) - meanAnomaly ) /
            ( 1- eccentricity * cos(eccAnomaly) ))
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
    meanAnomaly += meanMotion * (timeIncrement % orbitalPeriod(ko))

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
function createCircPecmeo( radius::Number, n_satellites::Tuple{Int, Int, Int}, cbody::Body,
    satSpacing::Tuple{<:Number, <:Number, <:Number};
    initialOrbitShift::Tuple{<:Number, <:Number, <:Number}= (0.0, 0.0, 0.0),
    equatorialRotation::Number = 0.0,
    inclination::Number = 0.0 )

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

createCircPecmeo( radius::Number, n_satellites::Tuple{Int, Int, Int}, cbody::Body;
    initialOrbitShift::Tuple{<:Number, <:Number, <:Number}= (0.0, 0.0, 0.0),
    equatorialRotation::Number = 0.0, inclination::Number=0.0 )=
    createCircPecmeo(radius, n_satellites, cbody, 2*pi ./ n_satellites;
    initialOrbitShift = initialOrbitShift,
    equatorialRotation = equatorialRotation, inclination = inclination)
