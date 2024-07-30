using DifferentialEquations
using StaticArrays
using GLMakie
using Makie.GeometryBasics
using DataStructures: CircularBuffer

include("cartpole.jl")
# A. Kroll, H. Schulte / Applied Soft Computing 25 (2014) 496–513


linear_cartpole_system = linear_cartpole(cartpole(;Fc = 0.0))

pulse(t,dt) = .5 + .5tanh(1e6*t) - (.5 - .5tanh(-1e6*(t-dt)))
cl_cartpole_system = cl_cartpole(cartpole([0.3, 0.0, 4π/45, 0.0]), linear_cartpole_system, 1e1I, diagm([5e5,1e4,1e3,1e3]),t -> 100pulse(t-5,.1))
# cl_cartpole_system = cl_cartpole(cartpole([0.3, 0.0, π/9, 0.0]), linear_cartpole_system, 1e5I, diagm([5e4,1e-1,1e3,1e3]),t -> 100pulse(t-5,.1))
cp = cl_cartpole_system
u0 = cl_cartpole_system.u0
l = cartpole().p[4]

# Solve diffeq with constant step for smoother curves
# diffeq = (adaptdt = 1/30/5)
include("makefig_cartpole.jl")

# %% 7. Interactive application
# Makie.jl has tremendously strong capabilities for real-time
# interactivity. To learn all of this takes time of course,
# and you'll need to consult the online documentation.
# Here we will do two interactions: 1) a play/stop button
# 2) clicking on the screen and getting a new initial condition!

fig, integ, rod, ball, cart, traj = makefig(u0,cp)
# The run button is actually pretty simple, we'll add it below the plot
run = Button(fig[2,1]; label = "run", tellwidth = false)
# This button will start/stop an animation. It's actually surprisingly
# simple to do this. The magic code is:
isrunning = Observable(false)
on(run.clicks) do clicks; isrunning[] = !isrunning[]; end
on(run.clicks) do clicks
    @async while isrunning[]
        isopen(fig.scene) || break # ensures computations stop if closed window
        animstep!(integ, rod, ball, cart, traj)
        sleep(1/200) # or `yield()` instead
    end
end

# `on` an important Observables function when one starts
# doing advanced stuff. It triggers a piece of code once an observable
# triggers its update.

# We'll add one more interactive feature which will trigger once
# we click on the axis. Notice that by default makie performs a zoom
# once one clicks on the axis, so we'll disable this
ax = content(fig[1,1])
Makie.deactivate_interaction!(ax, :rectanglezoom)
# and we'll add a new trigger using the `select_point` function:
spoint = select_point(ax.scene, marker = :circle)

# Now we see that we can click on the screen and the `spoint` updates!

# Okay, let's do something useful when it triggers
# function θωcoords(x, y)
#     θ = atan(y,x) + π/2
#     return SVector(θ,0,0,0)
# end

on(spoint) do z
    x, y = z
    s = sign(integ[1] - x)
    # u = θωcoords(x, y)
    FΔt = s*200*integ.dt
    _cp = cartpole()
    M = _cp.p[2]
    v_cart = (FΔt + M*integ[2])/M
    reinit!(integ, [integ[1], v_cart, integ[3], integ[4]])
    # Reset tail and balls to new coordinates
    # x1,x2,y1,y2 = xycoords(integ)
    # traj[] .= fill(Point2f(x2, y2), length(traj[]))
    # traj[] = traj[]
    # rod[] = [Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)]
    # balls[] = [Point2f(x1, y1), Point2f(x2, y2)]
end
