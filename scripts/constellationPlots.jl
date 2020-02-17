using Plots, LinearAlgebra, Random, Statistics, DoubleFloats
plotly()
include("../src/NaviSimu_lessModule.jl")
# using Main.NaviSimu

moonSat = KeplerOrbit(moon.radius + 100e3, 0.0,
    deg2rad(0), 0.0, 0.0, 0.0, moon)
pecmeo_nav = lunarnavPECMEODesigner( [0.431459575985078,
 0.999729596168082,
 0.962413062424563,
 0.999922037056928])

pecmeo_level = lunarnavPECMEODesigner( [0.431205002151572,   0.511002410029895,   0.285692220934846,   0.639705340836302])

#Define constellations
receiverOrbit = moonSat
navcon = pecmeo_nav

line_red = RGB(0.639, 0.192, 0.098)
line_blue = RGB(0.075, 0.31, 0.40)

# Simple PECMEO
begin
    orbit1 = KeplerOrbit(14e6, 0.0, 0.0, 0.0, 0.0, 0.0, earth)
    orbit2 = KeplerOrbit(14e6, 0.0, deg2rad(90), deg2rad(90), 0.0, 0.0, earth)
    orbit3 = KeplerOrbit(14e6, 0.0, deg2rad(90), 0.0, 0.0, 0.0, earth)
    PECMEO_simple = KeplerConstellation([orbit1, orbit2, orbit3])
    psimple = plotConstellation(PECMEO_simple; body_radius=6378.1, body_type=:solid, line_color = line_red, show_sats=false)
    plot!(xlabel="x [km]", ylabel="y [km]", zlabel= "z [km]", legend=false)
end


# Skewed constellation
begin
    psk = plotConstellation(pecmeo_nav; body_radius=6378.1, body_type=:solid, line_color = line_red)
    plot!(title = "Skewed PECMEO constellation", xlabel="x [km]", ylabel="y [km]", zlabel= "z [km]", legend=false)
end


# Level constellation
begin
    plotConstellation(pecmeo_level; body_radius=6378.1, body_type=:solid, line_color = line_red)
    plot!(title = "Level PECMEO constellation", xlabel="x [km]", ylabel="y [km]", zlabel= "z [km]", legend=false)
end


        # Optimization examples
# zero-vector
neutral_constel = lunarnavPECMEODesigner( [0.5,  0.5,   0.0,   0.0])
begin
    # neutral_constel = lunarnavPECMEODesigner( [0.0,  0.0,   0.0,   0.0])
    p_n = plotConstellation(neutral_constel; body_radius=6378.1, body_type=:solid, line_color = line_red)
    plot!(title = "Neutral PECMEO constellation", xlabel="x [km]", ylabel="y [km]", zlabel= "z [km]", legend=false)
end

# param 1
begin
    p_n = plotConstellation(neutral_constel; body_radius=6378.1, body_type=:solid, line_color = line_red)
    con1 = lunarnavPECMEODesigner( [0.5 + 0.5/3 , 0.5 ,   0.0,   0.0])
    # p1 = plotConstellation(neutral_constel; body_radius=6378.1, body_type=:solid, line_color = line_red)
    p1 = plotConstellation(con1; body_type = false, use_plot=p_n, line_color = line_blue)
    plot!(xlabel="x [km]", ylabel="y [km]", zlabel= "z [km]", legend=false)
end

# param 2
begin
    p_n = plotConstellation(neutral_constel; body_radius=6378.1, body_type=:solid, line_color = line_red)
    con2 = lunarnavPECMEODesigner( [0.5  , 0.5+ 0.5/3 ,   0.0,   0.0])
    # p1 = plotConstellation(neutral_constel; body_radius=6378.1, body_type=:solid, line_color = line_red)
    p2 = plotConstellation(con2; body_type = false, use_plot=p_n, line_color = line_blue)
    plot!(xlabel="x [km]", ylabel="y [km]", zlabel= "z [km]", legend=false)
end

# param 3
begin
    p_n = plotConstellation(neutral_constel; body_radius=6378.1, body_type=:solid, line_color = line_red)
    con3 = lunarnavPECMEODesigner( [0.5 , 0.5 ,   1/8,   0.0])
    # p1 = plotConstellation(neutral_constel; body_radius=6378.1, body_type=:solid, line_color = line_red)
    p3 = plotConstellation(con3; body_type = false, use_plot=p_n, line_color = line_blue)
    plot!(xlabel="x [km]", ylabel="y [km]", zlabel= "z [km]", legend=false)
end

# param 4
begin
    p_n = plotConstellation(neutral_constel; body_radius=6378.1, body_type=:solid, line_color = line_red)
    con4 = lunarnavPECMEODesigner( [0.5 , 0.5 ,   0.0,   1/8])
    # p1 = plotConstellation(neutral_constel; body_radius=6378.1, body_type=:solid, line_color = line_red)
    p4 = plotConstellation(con4; body_type = false, use_plot=p_n, line_color = line_blue)
    plot!(xlabel="x [km]", ylabel="y [km]", zlabel= "z [km]", legend=false)
end
