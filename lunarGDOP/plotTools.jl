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

# Function to create a plot of the constellation in Plots
function plotConstellation(constellation::KeplerConstellation;
    body_radius= earth.radius/1e3, body_type=:lines,
    body_color = RGB(0.306, 0.51, 0.137),
    orbits = true,
    use_plot::Plots.Plot= Plots.plot())

    # Detail in the plotted sphere
    body_resolution = 36
    if (body_type==:lines)
        # Sphere constructed from grid lines
        plotLineSphere(body_radius, body_resolution; color = body_color)
    elseif (body_type==:solid)
        # Solid sphere
        earthSphere = sphereCoordinates(body_radius, 25)
        Plots.plot!(earthSphere.x, earthSphere.y, earthSphere.z,
            st=:surface, c=:Greens, cb=false)
    end

    # Plot the satellites and their orbits
    nsats = size(constellation)
    for isat in 1:nsats
        if orbits
            # Get the position series for the orbit
            x_pos = []
            y_pos = []
            z_pos = []
            for t in range(0, findOrbitalPeriod(constellation[isat]), length = 100)
                cartOrbit = keplerToCartesian(propagateKeplerOrbit(constellation[isat], t))
                append!(x_pos, cartOrbit.x/1e3)
                append!(y_pos, cartOrbit.y/1e3)
                append!(z_pos, cartOrbit.z/1e3)
            end
            Plots.plot!(x_pos, y_pos, z_pos, w= 2.0, c=RGB(0.569, 0.271, 0.153))
        end
    end
    # Satellite location marker size depends on backend
    sat_size = 10
    if Plots.backend_name() in [:plotly, :plotlyjs]
        sat_size = 2
    elseif Plots.backend_name() in [:gr]
        sat_size = 5
    end
    for isat in 1:nsats
        # Get the current positon only
        cartOrbit = keplerToCartesian(constellation[isat])
        x_pos = [cartOrbit.x/1e3]
        y_pos = [cartOrbit.y/1e3]
        z_pos = [cartOrbit.z/1e3]

        # Show the current satellite location
        Plots.scatter!([x_pos[1]], [y_pos[1]], [z_pos[1]], markersize=sat_size, c=RGB(0.51, 0.137, 0.251))
    end

    # Ensure the plot opens when running function by returning plot
    return use_plot
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

#Add the lines of a sphere to current plot
function plotLineSphere(radius::Number, n::Int; color=RGB(0.306, 0.51, 0.137))
    #Create sphere point coordinates
    sdata = sphereCoordinates(radius, n)

    for i in 1:n
        plot!(sdata.x[i,:], sdata.y[i,:], sdata.z[i,:], w=1.0, c=color, legend=false)
        plot!(sdata.x[:,i], sdata.y[:,i], sdata.z[:,i], w=1.0, c=color, legend=false)
    end
end


#Create animation of a constellation
function animateConstellation(constellation::KeplerConstellation, frames::Int;
    path::String = "tmp/", animation_time::Number=0,
    size::Tuple{Int, Int} = (600, 400)
    )
    #Check to see if we can use progress meter
    use_progressMeter = :ProgressMeter in names(Main, imported=true)
    if use_progressMeter
        prog = Progress(frames, 1)
    end

    if (animation_time == 0)
        timestep = findOrbitalPeriod(constellation[1]) / (frames+1)
    end

    time = 0
    for i_frame in 1:frames
        p = plot(size = size)
        con_tmp = propagateKeplerOrbit(constellation, time)
        plotConstellation(con_tmp; use_plot=p, body_type=:lines)

        time += timestep
        savefig("$(path)frame_$(i_frame).png")
        if use_progressMeter
            next!(prog)
        end
    end

    return 0
end


function plotConstellationConnections(receiver::Tuple{Float64,Float64,Float64},
    constellation::Array{Tuple{Float64, Float64, Float64}, 1},
    shadowBody::Body)
    con = constellation

    p = plot()
    los = hasLineOfSight(receiver, con, bodyPosition(shadowBody.name, 0), shadowBody.radius)

    for isat in 1:size(con, 1)
        scatter!(con[isat], c=:Gray)
        if los[isat]
            plot!([con[isat][1], receiver[1]],
                [con[isat][2], receiver[2]],
                [con[isat][3], receiver[3]], c=:Blue)
        end
    end
    scatter!(receiver, c=:Red)
    plotLineSphere(shadowBody.radius, 36)

    return p
end
