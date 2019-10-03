# Plain Runge-Kutta 4 integration function
# Integrates from t=0 to t=tspan with dt = stepSize
# odefun = function containing ODE(y, t)
# y0 = initial condition vector
function integratorRK4(odefun, tspan, y0, stepSize)
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

# # RK4 integration to find data at target times.
# # Integrates with integStepSize constantly.
# # Propagation steps towards the targets are not used for further proapgation
# # in order to keep the step size constant
# function rk4Targets(odefun::Function, targetTimes, y0, t0, integStepSize)
#     h = integStepSize
#     t_curr = t0
#     y_curr = y0
#     dydt_curr = odefun(y_curr, t_curr)
#
#     target_iter = 1
#     target_curr = targetTimes[target_iter]
#     sort!(targetTimes)
#
#     #Create Empty data vectors
#     y = Array{Float64}(undef, 0)  #State vector, concentrated vertically to be 1D
#     dydt = Array{Float64}(undef, 0)
#     t = Array{Float64}(undef, 0)  #Time vector
#
#     #Perform integration steps until time is beyond timespan
#     while t_curr < targetTimes[end]
#
#         #Propagate untill a target time is in reach of a timestep
#         while (t_curr + h > target_curr)
#
#             #Progress one step
#             step=rk4Step(odefun, t_curr, y_curr, dydt_curr, h)
#
#             #Propagate step
#             y_curr = step.y
#             t_curr = step.t
#             dydt_curr = step.dydt
#         end
#
#         # Propagate towards the target Time
#         target=rk4Step(odefun, t_curr, y_curr, dydt_curr, target_curr - t_curr)
#
#         # Store data
#         append!(y, target.y)
#         append!(dydt, target.dydt)
#         append!(t, target.t)
#
#         # Set next target
#         target_iter +=1
#         target_curr = targetTimes[target_iter]
#     end
#     y = reshape(y, size(y0)[1], size(t)[1]) #Return states as 2d array
#     dydt = reshape(dydt, size(y0)[1], size(t)[1]) #Return states as 2d array
#
#     return (y=y, dydt=dydt, t=t)
# end

# Interpolation within rk4 integration
# Finds inertial times at which the clock shows the times in interPolValues
# initTime is the start of the integration
# propagationStepSize: step size of the propagation of the timeODE
# interParam: location in vector of timeODE output with local clock time
function integration_interpolationRK4(timeODE::Function, interpolValues,
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
            while error > 1e-20 && iter < 100
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
# odefun: ODE(y, t)
# t_curr: current time (end time of last step)
# y_curr: y(t_curr)
# dydt_curr: dy/dt (t_curr) -> ode evaluation of last step to reduce number of evaluations
# stepSize: timestep -> dt
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
