using DifferentialEquations
using StaticArrays
using GLMakie
using Makie.GeometryBasics
using DataStructures: CircularBuffer

include("cartpole.jl")
# A. Kroll, H. Schulte / Applied Soft Computing 25 (2014) 496–513


linear_cartpole_system = linear_cartpole(cartpole(;Fc = 0.0))
cl_cartpole_system = cl_cartpole(cartpole([0.3, 0.0, 4π/45, 0.0]), linear_cartpole_system, 1e1I, diagm([5e5,1e4,1e3,1e3]),t -> 100pulse(t-5,.1))
# cl_cartpole_system = cl_cartpole(cartpole([0.3, 0.0, π/9, 0.0]), linear_cartpole_system, 1e5I, diagm([5e4,1e-1,1e3,1e3]),t -> 100pulse(t-5,.1))
cp = cl_cartpole_system
u0 = cl_cartpole_system.u0
l = cartpole().p[4]

# Solve diffeq with constant step for smoother curves
# diffeq = (adaptdt = 1/30/5)

function xycoords(state)
   θ1 = state[3]
   x1 = state[1]
   y1 = 0.0
   x2 = x1 - l * sin(θ1)
   y2 = y1 + l * cos(θ1)
   return x1,x2,y1,y2
end

function progress_for_one_step!(integ)
    step!(integ)
    return xycoords(integ)
end

function animstep!(integ, rod, ball, cart, traj)
    x1,x2,y1,y2 = progress_for_one_step!(integ)
    rod[] = [Point2f(x1, y1), Point2f(x2, y2)]
    ball[] = [Point2f(x1, y1), Point2f(x2, y2)]
    cart[] = [[Point2f(x1, y1)] .+ Point2f[(-.1, -.05), (-.1, 0.0), 
                                           ( .1,  .00), ( .1, -.05)]]
    rod  = [Point2f(x1, y1), Point2f(x2, y2)]
    push!(traj[], Point2f(x2, y2))
    traj[] = traj[] # <- important! Updating in-place the value of an
                    # `Observable` does not trigger an update!
end

# cool it works. Let's wrap up the creation of the observables
# and plots in a function (just to re-initialie everything)
function makefig(u0,cp)
    # cp = cartpole(;fₓ = one)
    integ = init(cp, Tsit5(); adaptive = false, dt = .005)
    
    x1,x2,y1,y2 = xycoords(u0)
    cart = Observable([[Point2f(x1, y1)] .+ Point2f[(-.1, -.05), (-.1, 0.0), 
                                        (.1, 0.0), (.1, -.05)]])
    rod  = Observable([Point2f(x1, y1), Point2f(x2, y2)])
    ball = Observable([Point2f(x1, y1), Point2f(x2, y2)])
    
    tail = 100   # length of plotted trajectory, in units of `dt`
    traj = CircularBuffer{Point2f}(tail)
    fill!(traj, Point2f(x2, y2)) # add correct values to the circular buffer
    traj = Observable(traj) # make it an observable
    fig = Figure(); display(fig)
    ax = Axis(fig[1,1])
    
    c = to_color(:purple)
    tailcol = [RGBAf(c.r, c.g, c.b, (i/tail)^2) for i in 1:tail]
    poly!(ax, cart, color = :blue)
    lines!(ax, rod; linewidth = 4, color = :purple)
    scatter!(ax, ball; marker = :circle, strokewidth = 2, 
        strokecolor = :purple,
        color = :black, markersize = [12,8]
    )
    lines!(ax, traj; linewidth = 3, color = tailcol)
    ax.title = "Cartpole"
    ax.aspect = DataAspect()
    xlims!(ax, -2l, 2l)
    ylims!(ax, -1.5l, 1.5l)
    # also return the figure object, we'll ues it!
    return fig, integ, rod, ball, cart, traj
end
    
# %% 6. Save animations to videos
fig, integ, rod, ball, cart, traj = makefig(u0,cp)
frames = 1:1200
record(fig, "video.mp4", frames; framerate = 40) do i # i = frame number
    for j in 1:5 # step 5 times per frame
        animstep!(integ, rod, ball, cart, traj)
    end
    # any other manipulation of the figure here...
end # for each step of this loop, a frame is recorded

# # %% 7. Interactive application
# # Makie.jl has tremendously strong capabilities for real-time
# # interactivity. To learn all of this takes time of course,
# # and you'll need to consult the online documentation.
# # Here we will do two interactions: 1) a play/stop button
# # 2) clicking on the screen and getting a new initial condition!

# fig, integ, rod, balls, traj = makefig(u0)
# # The run button is actually pretty simple, we'll add it below the plot
# run = Button(fig[2,1]; label = "run", tellwidth = false)
# # This button will start/stop an animation. It's actually surprisingly
# # simple to do this. The magic code is:
# isrunning = Observable(false)
# on(run.clicks) do clicks; isrunning[] = !isrunning[]; end
# on(run.clicks) do clicks
#     @async while isrunning[]
#         isopen(fig.scene) || break # ensures computations stop if closed window
#         animstep!(integ, rod, balls, traj)
#         sleep(1/60) # or `yield()` instead
#     end
# end

# `on` an important Observables function when one starts
# doing advanced stuff. It triggers a piece of code once an observable
# triggers its update.

# # We'll add one more interactive feature which will trigger once
# # we click on the axis. Notice that by default makie performs a zoom
# # once one clicks on the axis, so we'll disable this
# ax = content(fig[1,1])
# Makie.deactivate_interaction!(ax, :rectanglezoom)
# # and we'll add a new trigger using the `select_point` function:
# spoint = select_point(ax.scene, marker = :circle)

# # Now we see that we can click on the screen and the `spoint` updates!

# # Okay, let's do something useful when it triggers
# function θωcoords(x, y)
#     θ = atan(y,x) + π/2
#     return SVector(θ,0,0,0)
# end

# on(spoint) do z
#     x, y = z
#     u = θωcoords(x, y)
#     reinit!(integ, u)
#     # Reset tail and balls to new coordinates
#     x1,x2,y1,y2 = xycoords(u)
#     traj[] .= fill(Point2f(x2, y2), length(traj[]))
#     traj[] = traj[]
#     rod[] = [Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)]
#     balls[] = [Point2f(x1, y1), Point2f(x2, y2)]
# end
