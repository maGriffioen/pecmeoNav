using Plots, LinearAlgebra
gr()
include("../src/NaviSimu_adder.jl")
using Main.NaviSimu

moonSat = KeplerOrbit(moon.radius + 100e3, 0.0,
    deg2rad(0), 0.0, 0.0, 0.0, moon)
pecmeo_333 = createCircPecmeo(26.4e6, (3, 3, 3), earth,
    (5.372554642808982, 2.5348806135682063, 0.7211160718288704);
    initialOrbitShift=(5.372554642808982, 2.5348806135682063, 0.7211160718288704),
    equatorialRotation = 0*pi/8,
    inclination = 0.0651384236528878)

#Define constellations
receiverOrbit = moonSat
navcon = pecmeo_333

#Find transmitter position and true time during transmission (Light time effect)
function transmitterFinder(receptionTime::Number, receiverPos::Tup3d, transmitterOrbit::Orbit)
    travelTime = 0.0    #Assume instant signal as first estimate
    correction = 1      #Force loop to start
    transTime = receptionTime   #Initialization
    transPos = (0.0, 0.0, 0.0)  #Initialization

    nIter = 0   #Iteration counter
    # Iterate while no convergence of smaller than xx seconds is reached
    while (abs(correction) > 1e-15)   #TODO: Tweak value for balance between stability and numerical accuracy
        # Calculate new transmitter position & time
        transTime = receptionTime - travelTime      #Time of transmission of signal
        transPos = globalPosition(transmitterOrbit, transTime)  #Position of transmitting satellite
        distance = norm(transPos .- receiverPos)    #Travel distance of signal
        travelTime = distance/lightConst            #Travel time of signal

        # Find applied correction for the decision if convergence has been reached
        correction = receptionTime - transTime - travelTime
        nIter +=1
    end

    return (t_true = transTime, pos = transPos, travelTime = travelTime,
        lastCorrection = correction, iterations = nIter)
end

function clockODE(x::Array{<:Number, 1}, t::Number)
    dxdy = zeros(length(x))
    state = globalState(receiverOrbit, t)
    velocity = norm(state[4:6])
    c = lightConst

    dxdy[1] = 1/sqrt(1- (velocity^2 / c^2)) #Local time rate
    dxdy[2] = dxdy[1] + (1 - 1/sqrt(1- (1925.63^2 / c^2)))  #Clock time rate
    #x[1] = t',
    #x[2] = t'_c
    return dxdy
end

function rk4(odefun, tspan, y0, stepSize)
    h = stepSize
    t_curr = 0.0
    y_curr = y0
    dydt_curr = odefun(y_curr, t_curr)

    #Create Empty data vectors
    y = []  #State vector, concentrated vertically to be 1D
    dydt = []
    t = []  #Time vector

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
    y = []  #State vector, concentrated vertically to be 1D
    dydt = []
    t = []  #Time vector

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
    t = []  #Time vector
    y = [] #State vector, concentrated vertically to be 1D
    dydt = [] #State derivative


    #Data vector for interpolated values
    t_inter = []
    y_inter = []
    dydt_inter = []

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


measurements = 5.0:10.2:1000.0
results = rk4(clockODE, 1000, [0.0, 0.0], 10)
inter = rk4measurementFinder(clockODE, measurements, 0.0, 1.0; interParam=2)
transmitter = transmitterFinder(inter.t[1], globalPosition(receiverOrbit, inter.t[1]), navcon[1])

trueOffset = results.y[1,:]- results.t
clockOffset= results.y[2,:]- results.t


plot(xlabel="Time [hours]", ylabel="Time offset [s]")
plot!(results.t / 3600, trueOffset, label="Local time offset")
plot!(results.t / 3600, clockOffset, label="Clock offset")

plot(xlabel="Time [hours]", ylabel="Time offset [s]")
plot!(inter.t, inter.y[2, :]-inter.t, label="Time offset from inertial time during measurement")
plot!(inter.t, inter.y[2,:]-measurements, label="Numerical measurement time error")
