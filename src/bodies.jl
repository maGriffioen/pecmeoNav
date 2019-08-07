#Quick numbers from google.
#Replace with well-sourced numbers!!
#https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
gravConst = 6.674e-11   #m^3 kg^-1 s^-2 Lissauer, de Pater
earthState(time::Number) = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
earth = Body(
    name = "Earth",
    gravitationalParameter =5.9724e24 * gravConst,
    radius = 6378.1e3, # Equatorial radius, pessimistic case for shadowing
    stateFunction = earthState
)
moonState(time::Number) = globalState(propagateKeplerOrbit(lunarOrbit, time))
moon = Body(
    name = "Moon",
    gravitationalParameter = 0.07346e24 * gravConst,
    radius = 1738.1e3,
    stateFunction = moonState
)
#https://ssd.jpl.nasa.gov/?sat_elem
lunarOrbit = KeplerOrbit(0.3844e9, 0.0554, deg2rad(5.16), deg2rad(318.15), deg2rad(125.08), deg2rad(135.27) , earth)
lightConst = 299792458          # Replace me with a proper source :(
bolzmanConst = 1.3806e-023      # Replace me with a proper source :(
println("Add proper speed of light & bolzman constant with sources!")



function bodyPosition(bodyName::String, time::Number)
    if (bodyName == "Earth")
        pos3d = (0.0, 0.0, 0.0)
    elseif (bodyName == "Moon")
        pos3d = position(propagateKeplerOrbit(lunarOrbit, time))
    else
        error("Body is unknown")
    end
    return pos3d

end
# bodyPosition(body::Body, time::Number) = bodyPosition(body.name, time)
bodyPosition(body::Body, time::Number) = Tup3d(body.stateFunction(time)[1:3])


function globalState(body::Body, time)
    # if (body.name == "Earth")
    #     state = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    # elseif (body.name == "Moon")
    #     state = globalState(propagateKeplerOrbit(lunarOrbit, time))
    # else
    #     error("Body is unknown")
    # end
    return body.stateFunction(time)
end
