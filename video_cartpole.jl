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

# %% 6. Save animations to videos
fig, integ, rod, ball, cart, traj, graph = makefig(u0,cp)
frames = 1:1200
record(fig, "video.mp4", frames; framerate = 40) do i # i = frame number
    for j in 1:5 # step 5 times per frame
        animstep!(integ, rod, ball, cart, traj, graph)
    end
    # any other manipulation of the figure here...
end # for each step of this loop, a frame is recorded
