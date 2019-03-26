#Quick numbers from google.
#Replace with well-sourced numbers

struct Body
    name::String
    gravitationalParameter::Float64
    radius::Float64
end

Body(; name = "", gravitationalParameter, radius) =
    Body(name, gravitationalParameter, radius)

earth = Body(
    name = "Earth",
    gravitationalParameter =3.986004418e14,
    radius = 6371e3
)
moon = Body(
    name = "Moon",
    gravitationalParameter = 4.9048695e12,
    radius = 1737e3
)
lunarOrbit = KeplerOrbit(384e6, 0.0, deg2rad(23.4 + 5.14), 0.0, 0.0, 0.0, earth)

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
bodyPosition(body::Body, time::Number) = bodyPosition(body.name, time)
