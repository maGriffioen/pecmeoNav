function sphereCoordinates(radius::Number, n::Int)
    #Define polar coordinates of points on sphere based on the n
    u = collect(range( 0, stop = 2*pi, length = n ))
    v = collect(range( 0, stop = pi, length = n))

    #Define coordinates of points on sphere
    x = radius * cos.(u) * sin.(v)'
    y = radius * sin.(u) * sin.(v)'
    z = radius * ones(n) * cos.(v)'

    #Return as a named tuple
    return ( x=x, y=y, z=z )
end
sphereCoordinates(n::Int) = sphereCoordinates(1::Int, n)

function plotConstellation(constellation::KeplerConstellation; orbits = false)
    earthSphere = sphereCoordinates(earth.radius/1e3, 25)
    Plots.plot(earthSphere.x, earthSphere.y, earthSphere.z, st=:surface, c=:Greens, zcolor=[-4,4], cb=false)
    nsats = size(constellation)
    for isat in 1:nsats
        cus_pos = position(constellation[isat])
        Plots.scatter!([cus_pos[1]/1e3], [cus_pos[2]/1e3], [cus_pos[3]/1e3], markersize=1, c="gray")
        if orbits
            x_pos = []
            y_pos = []
            z_pos = []
            for t in range(0, findOrbitalPeriod(iss), length = 40)
                cartOrbit = keplerToCartesian(propagateKeplerOrbit(constellation[isat], t))
                append!(x_pos, cartOrbit.x)
                append!(y_pos, cartOrbit.y)
                append!(z_pos, cartOrbit.z)
            end
            Plots.plot!(x_pos/1e3, y_pos/1e3, z_pos/1e3, w= 4.0)
        end
    end
end

function pointsInFullOrbit(kepOrbit::KeplerOrbit; n = 50)
    x_pos = zeros(n)
    y_pos = zeros(n)
    z_pos = zeros(n)

    for (i, t) in enumerate( range(0, findOrbitalPeriod(kepOrbit), length = n) )
        cartOrbit = keplerToCartesian(propagateKeplerOrbit(iss, t))
        x_pos[i] = cartOrbit.x
        y_pos[i] = cartOrbit.y
        z_pos[i] = cartOrbit.z
    end

    return (x = x_pos, y = y_pos, z = z_pos)
end
