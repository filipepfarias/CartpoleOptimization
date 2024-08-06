using JuMP, Optim

include("cartpole.jl")
cp = cartpole();


# Collocation methods
function swingup_cartpole(cp::ODEProblem=cartpole(),
            t=(0.0,2.0),
            x=([0,0,π,0],[0,0,0,0]),
            N=100)
    Tk = range(t...,N)
    Uk = zeros(N)
    d = size(cp.u0,1)
    model = Model(Optim.Optimizer);

    q = [@variable(model, [1:d]) for _ in 1:N]
    u = [@variable(model, [1:1]) for _ in 1:N]

    @inbounds function J(u,Tk,N)
        J = 0.0
        for k in 1:N-1
            J += 0.5(Tk[k+1]-Tk[k])*(u[k+1][1]+u[k][1])
        end
        return J
    end    
    @objective(model, Min, J(u,Tk,N))

    @inbounds function h(cp, uk, ukk, xk, xkk,tk,tkk)
        # d = size(x,1)
        # h = zeros(d,N-1);
        # p = cp.p[2:end]
        # for k in 1:N-1
        #     xk  = x[k]
        #     xkk = x[k+1]
            fk = cp.f.f(xk,(uk[1],p...),tk)
            fkk = cp.f.f(xkk,(ukk[1],p...),tk)

            return 0.5(tkk-tk)*(fkk+fk) - (xkk - xk)
        # end
    end
    @constraint(model, [k = 1:N-1], 0 .== h(cp,u[k],u[k+1],q[k],q[k+1],Tk[k],Tk[k+1]))
    @constraint(model, x[1] == q[1])
    @constraint(model, x[2] == q[N])
    # @constraint(model, [k = 1:N-1],  1 <= q[k][1] < 1)

    optimize!(model)
    return ([value.(uk)[1] for uk in u],[value.(qk) for qk in q])
end