using DifferentialEquations
using StaticArrays
using GLMakie
using Makie.GeometryBasics
using DataStructures: CircularBuffer

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

    cart_x
    cart_ẋ
    cart_θ
    cart_θ̇
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
    