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

function animstep!(integ, rod, ball, cart, traj, graph)
    x1,x2,y1,y2 = progress_for_one_step!(integ)
    rod[] = [Point2f(x1, y1), Point2f(x2, y2)]
    ball[] = [Point2f(x1, y1), Point2f(x2, y2)]
    cart[] = [[Point2f(x1, y1)] .+ Point2f[(-.1, -.05), (-.1, 0.0), 
                                           ( .1,  .00), ( .1, -.05)]]
    rod  = [Point2f(x1, y1), Point2f(x2, y2)]
    push!(traj[], Point2f(x2, y2))
    traj[] = traj[] # <- important! Updating in-place the value of an
                    # `Observable` does not trigger an update!

    K = integ.sol.prob.p[1]
    r = integ.sol.prob.p[2](integ.t)
    u = -K*(integ.u - r)    

    for (i,i_graph) in enumerate(graph[1:4])
        push!(i_graph[], Point2f(integ.t,integ.u[i]))
        i_graph[] = i_graph[]
    end

    push!(graph[5][], Point2f(integ.t,u[1]))
    graph[5][] = graph[5][]
    # graph[6][] = graph[6][]
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
    fig = Figure(size=(1000,500)); display(fig)
    ax = Axis(fig[1:3,1])
    
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

    nt = 600
    t = range(0.0, step=-integ.dt, length=nt)

    graph = [CircularBuffer{Point2f}(nt) for _ in 1:5]
    for i_graph in graph
        pushfirst!.([i_graph], Point2f.(t,0.0))
    end

    graph = Observable.(graph)

    ax_graph_1 = Axis(fig[1,2])
    ax_graph_1.title = "Cart position"
    ax_graph_1.ylabel = "Distance (m)"
    ax_graph_1.xlabel = "time (s)"
    ax_graph_2 = Axis(fig[2,2])
    ax_graph_2.title = "Pole position"
    ax_graph_2.ylabel = "Angle (rad)"
    ax_graph_2.xlabel = "time (s)"
    ax_graph_3 = Axis(fig[3,2])
    ax_graph_3.title = "Input"
    ax_graph_3.ylabel = "Force (N)"
    ax_graph_3.xlabel = "time (s)"

    lines!(ax_graph_1, graph[1], label=L"x")
    lines!(ax_graph_1, graph[2], label=L"\dot{x}")
    lines!(ax_graph_2, graph[3], label=L"\theta")
    lines!(ax_graph_2, graph[4], label=L"\dot{\theta}")
    lines!(ax_graph_3, graph[5], label=L"u")
    # lines!(ax_graph_3, graph[6], label=L"r")

    axislegend(ax_graph_1)
    axislegend(ax_graph_2)
    axislegend(ax_graph_3)

    lift(graph[1]) do g
        xlims!(ax_graph_1, g[1][1], g[end][1])
        xlims!(ax_graph_2, g[1][1], g[end][1])
        xlims!(ax_graph_3, g[1][1], g[end][1])
    end

    ylims!(ax_graph_3,-130,130)

    # also return the figure object, we'll ues it!
    return fig, integ, rod, ball, cart, traj, graph
end
    