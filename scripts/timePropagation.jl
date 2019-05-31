using Plots, LinearAlgebra
gr()
include("../src/NaviSimu_adder.jl")
using Main.NaviSimu

moonSat = KeplerOrbit(moon.radius + 100e3, 0.0,
    deg2rad(0), 0.0, 0.0, 0.0, moon)
pecmeo32 = createCircPecmeo(26.4e6, (3, 3, 3), earth;
    initialOrbitShift=(0.0, (2/6)*pi, (2/6)*pi),
    equatorialRotation = 0*pi/8,
    inclination = 0.0)

#Define constellations
receiverOrbit = moonSat
navcon = pecmeo_333



# Find transmitter position and true time during transmission (Light time effect)
function transmitterFinder(receptionTime::Number, receiverPos::Tup3d, transmitterOrbit::Orbit)
    travelTime = 0.0    #Assume instant signal as first estimate
    correction = 1      #Force loop to start
    transTime = receptionTime   #Initialization
    transPos = (0.0, 0.0, 0.0)  #Initialization

    nIter = 0   #Iteration counter
    # Iterate while no convergence of smaller than xx seconds is reached
    while (abs(correction) > 1e-12)   #TODO: Tweak value for balance between stability and numerical accuracy
        # Calculate new transmitter position & time
        transTime = receptionTime - travelTime      #Time of transmission of signal
        transPos = globalPosition(transmitterOrbit, transTime)  #Position of transmitting satellite
        distance = norm(transPos .- receiverPos)    #Travel distance of signal
        travelTime = distance/lightConst            #Travel time of signal

        # Find applied correction for the decision if convergence has been reached
        correction = receptionTime - transTime - travelTime
        nIter +=1
        # println(nIter," ", correction)
    end
    # println(nIter)

    return (transmissionTime = transTime, pos = transPos, travelTime = travelTime,
        lastCorrection = correction, iterations = nIter)
end

# ODEs to describe time dilation of satellite with respect to Earth-fixed-rotating clock
function clockODE(x::Array{<:Number, 1}, t::Number)
    dxdy = zeros(length(x))
    state = globalState(receiverOrbit, t)
    velocity = norm(state[4:6])
    c = lightConst

    dxdy[1] = sqrt(1- (velocity^2 / c^2)) #Local time rate
    dxdy[2] = dxdy[1] * (1/sqrt(1- (1925.63^2 / c^2)))  #Clock time rate
    #x[1] = t',
    #x[2] = t'_c
    return dxdy
end

# Plain Runge-Kutta 4 integration scheme
# Integrates from t=0 to t=tspan with dt = stepSize
function rk4(odefun, tspan, y0, stepSize)
    h = stepSize
    t_curr = 0.0
    y_curr = y0
    dydt_curr = odefun(y_curr, t_curr)

    #Create Empty data vectors
    y = Array{Float64}(undef, 0)  #State vector, concentrated vertically to be 1D
    dydt = Array{Float64}(undef, 0)
    t = Array{Float64}(undef, 0)  #Time vector

    #Add initial conditions
    append!(y, y_curr)
    append!(dydt, dydt_curr)
    append!(t, t_curr)

    #Perform integration steps until time is beyond timespan
    while t_curr < tspan
        #See if it is a last reduced timestep
        if t_curr + stepSize > tspan
            h = tspan - t_curr
        end

        #Progress one step
        step=rk4Step(odefun, t_curr, y_curr, dydt_curr, h)

        #Store data
        append!(y, step.y)
        append!(dydt, step.dydt)
        append!(t, step.t)

        #Propagate step
        y_curr = step.y
        t_curr = step.t
        dydt_curr = step.dydt
    end
    y = reshape(y, size(y0)[1], size(t)[1]) #Return states as 2d array
    dydt = reshape(dydt, size(y0)[1], size(t)[1]) #Return states as 2d array

    return (y=y, dydt=dydt, t=t)
end

# RK4 integration to find data at target times.
# Integrates with integStepSize constantly.
# Propagation steps towards the targets are not used for further proapgation
# in order to keep the step size constant
function rk4Targets(odefun, targetTimes, y0, t0, integStepSize)
    h = integStepSize
    t_curr = t0
    y_curr = y0
    dydt_curr = odefun(y_curr, t_curr)

    target_iter = 1
    target_curr = targetTimes[target_iter]
    sort!(targetTimes)

    #Create Empty data vectors
    y = Array{Float64}(undef, 0)  #State vector, concentrated vertically to be 1D
    dydt = Array{Float64}(undef, 0)
    t = Array{Float64}(undef, 0)  #Time vector

    #Add initial conditions


    #Perform integration steps until time is beyond timespan
    while t_curr < targetTimes[end]

        #Propagate untill a target time is in reach of a timestep
        while (t_curr + h > target_curr)

            #Progress one step
            step=rk4Step(odefun, t_curr, y_curr, dydt_curr, h)

            #Propagate step
            y_curr = step.y
            t_curr = step.t
            dydt_curr = step.dydt
        end

        # Propagate towards the target Time
        target=rk4Step(odefun, t_curr, y_curr, dydt_curr, target_curr - t_curr)

        # Store data
        append!(y, target.y)
        append!(dydt, target.dydt)
        append!(t, target.t)

        # Set next target
        target_iter +=1
        target_curr = targetTimes[target_iter]


    end
    y = reshape(y, size(y0)[1], size(t)[1]) #Return states as 2d array
    dydt = reshape(dydt, size(y0)[1], size(t)[1]) #Return states as 2d array

    return (y=y, dydt=dydt, t=t)
end

#Interpolation within rk4 integration
function rk4measurementFinder(timeODE::Function, interpolValues::StepRangeLen,
    initTime::Number, propagationStepSize::Number; interParam::Int = 1)

    h = propagationStepSize
    odefun = timeODE

    t_curr = initTime
    y_curr = [initTime, initTime]
    dydt_curr = odefun(y_curr, t_curr)
    interpolInt_curr = 1
    interpolVal_tofind = interpolValues[interpolInt_curr]

    #Create Empty data vectors
    t = Array{Float64}(undef, 0)  #Time vector
    y = Array{Float64}(undef, 0) #State vector, concentrated vertically to be 1D
    dydt = Array{Float64}(undef, 0) #State derivative


    #Data vector for interpolated values
    t_inter = Array{Float64}(undef, 0)
    y_inter = Array{Float64}(undef, 0)
    dydt_inter = Array{Float64}(undef, 0)

    #Perform integration steps until time is beyond timespan
    while t_curr < interpolValues[end] && interpolInt_curr <= length(interpolValues)

        #Calculate propagation step
        new = rk4Step(odefun, t_curr, y_curr, dydt_curr, h)
            #new.y[1] -> Actual local time
            #new.y[2] -> Clock time
            #new.t -> Inertial time

        while (interpolVal_tofind >= y_curr[interParam] &&
            interpolVal_tofind < new.y[interParam] &&
            interpolInt_curr <= length(interpolValues))

            #Initial estimate for the root using linear interpolation
            root_dt = (new.t - t_curr)*((interpolVal_tofind - y_curr[interParam]) / (new.y[interParam]-y_curr[interParam]))
            root = rk4Step(odefun, t_curr, y_curr, dydt_curr, root_dt)
            error = abs(interpolVal_tofind - root.y[interParam])

            #Iterate With Newton-Raphson Method to find the root for measurement time
            iter = 0
            while error > 1e-13 && iter < 100
                root_dt -= (root.y[interParam] - interpolVal_tofind)/root.dydt[interParam]
                root = rk4Step(odefun, t_curr, y_curr, dydt_curr, root_dt)
                error = abs(interpolVal_tofind - root.y[interParam])
                iter += 1
            end
            if iter >1
                println("Many iterations: \t", iter, " @t_m= ", interpolVal_tofind)
            end
            # println(iter)

            #Append resulting root
            append!(t_inter, root.t)
            append!(y_inter, root.y)
            append!(dydt_inter, root.dydt)

            interpolInt_curr += 1
            # print(interpolInt_curr)
            # print(t_curr, "  ", root.t, "  ", new.t, "  ", root_dt, "  ", new.y[interParam], " ", interpolVal_tofind, "\n")
            if interpolInt_curr <= length(interpolValues)
                interpolVal_tofind = interpolValues[interpolInt_curr]
            end
        end

        #Propagate step
        y_curr = new.y
        t_curr = new.t
        dydt_curr = new.dydt
    end
    # y = reshape(y, size(y0)[1], size(t)[1]) #Return states as 2d array
    # dydt = reshape(dydt, size(y0)[1], size(t)[1]) #Return states as 2d array

    y_inter = reshape(y_inter, Int(length(y_inter)/length(t_inter)), length(t_inter))
    dydt_inter = reshape(dydt_inter, Int(length(dydt_inter)/length(t_inter)), length(t_inter))

    return (y = y_inter, t = t_inter, dydt = dydt_inter)
end

# Propagate single RK4 step, based on current time, y vector, dy/dt and step size.
# Returns new time, y vector and dy/dt
function rk4Step(odefun, t_curr, y_curr, dydt_curr, stepSize)
    h = stepSize
    k1 = h*dydt_curr
    k2 = h*odefun(y_curr + k1/2, t_curr + h/2)
    k3 = h*odefun(y_curr + k2/2, t_curr + h/2)
    k4 = h*odefun(y_curr + k3, t_curr + h)

    t = t_curr + stepSize
    y = y_curr + (k1 + 2*k2 + 2*k3 + k4) / 6
    dydt = odefun(y, t)

    return (y=y, t=t, dydt=dydt)
end

# Simulate the measurements sequentially
function simulateMeasurements(inertialTime::Number, receiverTime::Number, receiverOrbit::Orbit, navcon::Orbit;
    writeToFile = false, outputFile = "")
    nsats = size(navcon)
    codes = Array{Float64, 1}(undef, nsats)
    phases = Array{Float64, 1}(undef, nsats)
    avail = BitArray(undef, nsats)

    #Loop over navigation satellites for this epoch
    for prn in 1:nsats
        receiverPosition = globalPosition(receiverOrbit, inertialTime)
        transmitter = transmitterFinder(inertialTime, receiverPosition, navcon[prn])

        connection = true

        # Determine if in shadow of Earth or Moon
        # Check the line of sight both during reception and transmission. Assume no LOS when satellites are shadowed by body during either time.
        losEarly = NaviSimu.hasLineOfSightEarthMoon(receiverPosition, transmitter.pos, transmitter.transmissionTime)
        losLate = NaviSimu.hasLineOfSightEarthMoon(receiverPosition, transmitter.pos, inertialTime)
        los = (losEarly && losLate)

        # Check if satellites are connected and within line of sight
        available = los && connection

        avail[prn] = available
        if (available)
            #Simulate Code measurement
            measurementError = 0.0  # Error within receiver system, dependent on SNR
                # Neglecting transmitter clock error and relativistic effect -> will be largely recovered
            codes[prn] = (receiverTime - transmitter.transmissionTime) * lightConst + measurementError

            #Simulate Phase measurement
            phases[prn] = 0.0
        else
            # No measurement taken
            codes[prn] = NaN64
            phases[prn] = NaN64
        end
    end
    return (code = codes, phase = phases, avail = avail)
end

function simulateMeasurements(inertialTime::Array{<:Number}, receiverTime::Array{<:Number}, receiverOrbit::Orbit, navcon::Orbit)
    #Raise error if number of measurements is not clear
    if (length(inertialTime) != length(receiverTime))
        error("simulateMeasurements: Number of true times not equal to number of receiver times")
    end

    nepochs = length(inertialTime)     # Number of epochs
    nsats = size(navcon) #Total number of satellites in the navigation constellation
    codeObs = Array{Float64}(undef, nepochs, nsats)    #Per-epoch, per-prn code observation
    phaseObs = Array{Float64}(undef, nepochs, nsats)   #Per-epoch, per-prn phase observation
    availability = BitArray(undef, nepochs, nsats)     #Per-epoch, per-prn availability boolean

    for epoch = 1:nepochs
        msrmt = simulateMeasurements(inertialTime[epoch],
                    receiverTime[epoch],
                    receiverOrbit, navcon)
        codeObs[epoch, :] = msrmt.code
        phaseObs[epoch, :] = msrmt.phase
        availability[epoch, :] = msrmt.avail
    end

    return (codeObs = codeObs, phaseObs = phaseObs, availability = availability)
end

measurements = 5.0:100.2:10000.0
results = rk4(clockODE, 10000, [0.0, 0.0], 10)
inter = rk4measurementFinder(clockODE, measurements, 0.0, 1.0; interParam=2)
trueMeasurementTimes = inter.t  # Inertial times during measurements
localReceiverMeasurementTimes = inter.y[1, :]   # local time progression during measurements
receiverClockMeasurementTimes = inter.y[2, :]   # Time on receiver clock during measurements. Includes clock freq correction (and clock errors).
transmitter = transmitterFinder(inter.t[1], globalPosition(receiverOrbit, inter.t[1]), navcon[1])
singleMeasurement = simulateMeasurements(inter.t[20], inter.y[2, 20], moonSat, pecmeo_333)
allMeasurements = simulateMeasurements(inter.t, inter.y[2,:], moonSat, pecmeo_333)

trueOffset = results.y[1,:]- results.t
clockOffset= results.y[2,:]- results.t


p1 = plot(xlabel="Time [hours]", ylabel="Time offset [s]")
plot!(results.t / 3600, trueOffset, label="Local time offset")
plot!(results.t / 3600, clockOffset, label="Clock offset")

p2 = plot(xlabel="Time [hours]", ylabel="Time offset [s]")
plot!(inter.t / 3600, inter.y[2, :]-inter.t, label="Time offset from inertial time during measurement")
plot!(inter.t / 3600, inter.y[2,:]-measurements, label="Numerical measurement time error")
