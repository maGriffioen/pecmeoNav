using ProgressMeter, Random
include("../src/NaviSimu_lessModule.jl")

pecmeo_nav = lunarnavPECMEODesigner( [0.431459575985078,
 0.999729596168082,
 0.962413062424563,
 0.999922037056928])


meanerrors = []

orb_n = 1 #use 1, 4, 7
Random.seed!(orb_n)   #+1 for propagated case
reforb = pecmeo_nav[orb_n]
# reforb = KeplerOrbit(13e6, 0.0, deg2rad(45), 0.0, 0.0, 0.0, earth)

eph_error = 0.01
errors = []
x_errors = []
y_errors = []
z_errors = []
a_errors = []
theta_errors = []
i_errors = []
@showprogress for i in 1:10000
# sigma_a = eph_error * (1/sqrt(3))
# sigma_theta = eph_error * (1/sqrt(3)) * (1/reforb.a) * (1-cos(reforb.i))
# sigma_i = eph_error * (1/sqrt(3)) * (1/reforb.a) * (1+cos(reforb.i))
# sigma_raan = eph_error * (1/sqrt(3)) * (1/reforb.a) * sin(reforb.i)
# sigma_a = eph_error * (1/sqrt(3))*randn(Float64)
# sigma_theta = eph_error * (1/sqrt(3)) *(1/reforb.a) *randn(Float64)
# sigma_i = eph_error * (1/sqrt(3))*sqrt(2) * (1/reforb.a) * randn(Float64)
# sigma_raan = eph_error * (1/sqrt(3)) * (1/reforb.a) * sin(reforb.i)*0

sigma_a = eph_error  *randn(Float64)
sigma_theta = eph_error *(1/reforb.a) *randn(Float64)
sigma_i = eph_error *sqrt(2) * (1/reforb.a) * randn(Float64)
sigma_raan = eph_error * (1/sqrt(3)) * (1/reforb.a) * sin(reforb.i)*0
push!(a_errors, sigma_a)
push!(theta_errors, sigma_theta)
push!(i_errors, sigma_i)

rsserror = []
for t in range(0, orbitalPeriod(reforb); length=100)
    dt= 0   #0 or 1800
    orb = propagateKeplerOrbit(reforb, t)
    orb2 = KeplerOrbit(orb.a + sigma_a, orb.e, orb.i + sigma_i,
        orb.raan + sigma_raan, orb.aop, orb.tanom + sigma_theta, orb.cbody)
    orb = propagateKeplerOrbit(orb, dt)
    orb2 = propagateKeplerOrbit(orb2, dt)
    push!(errors, norm(globalPosition(orb) .- globalPosition(orb2)))
    push!(x_errors, globalPosition(orb)[1] - globalPosition(orb2)[1])
    push!(y_errors, globalPosition(orb)[2] - globalPosition(orb2)[2])
    push!(z_errors, globalPosition(orb)[3] - globalPosition(orb2)[3])
end
# push!(errors, mean(rsserror))

# println("Mean RSS error: ", mean(rsserror))
# plot(rsserror)
# push!(errors, mean(rsserror))
end
