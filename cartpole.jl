using DifferentialEquations
using StaticArrays
using ForwardDiff
using LinearAlgebra
using MatrixEquations
using JuMP, Ipopt
using DataInterpolations


function cartpole(
    q0::SVector=SVector(0.0, 0.0, 0.0, 0.0);
    u::Function = (x,t) -> zero(t), 
    M = 4.8, 
    m = 0.356, 
    l = 0.56, 
    g = 9.81,
    Fc = 4.9,
    dMf = 0.035)
    
    function f(q,p,t)
        p = (u=p.u([],t), M=M, m=m, l=l, g=g, Fc=Fc, dMf=dMf)
        return SVector{4}(cartpole_rule(q,p,t)...,)
    end
    
    return ODEProblem(f, 
    q0,(0.0,1.0),(u=u, M=M, m=m, l=l, g=g, Fc=Fc, dMf=dMf))
end


function cartpole(
    q0::Vector=[0.0, 0.0, 0.0, 0.0];
    u::Function = (x,t) -> zero(t), 
    M = 4.8, 
    m = 0.356, 
    l = 0.56, 
    g = 9.81,
    Fc = 4.9,
    dMf = 0.035)
    
    function f(dq,q,p,t)
        p = (u=p.u([],t), M=M, m=m, l=l, g=g, Fc=Fc, dMf=dMf)
        dq .= cartpole_rule(q,p,t)
    end
    
    return ODEProblem(f, 
    q0,(0.0,1.0),(u=u, M=M, m=m, l=l, g=g, Fc=Fc, dMf=dMf))
end


@inbounds function cartpole_rule(
    q::AbstractVector,
    p::@NamedTuple{u::T,
    M::Float64,
    m::Float64,
    l::Float64,
    g::Float64,
    Fc::Float64,
    dMf::Float64},
    t::Number) where T <: Union{Number,AbstractVariableRef}
    u, M, m, l, g, Fc, dMf = p
    x, ẋ, θ, θ̇ = q
    
    fc = Fc*tanh(1e6*ẋ)   
    
    ẍ = (fc - u + (-((-l*(fc - u + l*m*sin(θ)*(θ̇^2))*m*cos(θ)) / 
        (-M - m) + dMf*θ̇ - g*l*m*sin(θ))*l*m*cos(θ)) / 
        ((-(l^2)*(m^2)*(cos(θ)^2)) / (-M - m) - (l^2)*m) + 
        l*m*sin(θ)*(θ̇^2)) / (-M - m)
    θ̈ = ((-l*(fc - u + l*m*sin(θ)*(θ̇^2))*m*cos(θ)) / 
        (-M - m) + dMf*θ̇ - g*l*m*sin(θ)) / ((-(l^2)*(m^2)*(cos(θ)^2)) / 
        (-M - m) - (l^2)*m)
    
    return (ẋ, ẍ, θ̇, θ̈)
end

function linear_cartpole(
    cp::ODEProblem = cartpole(;Fc = 0.0), 
    q0::AbstractArray = [0.0, 0.0, 0.0, 0.0],
    t0 = 0.0)
    
    _p = (M=cp.p.M, m=cp.p.m, l=cp.p.l, g=cp.p.g, Fc=cp.p.Fc, dMf=cp.p.dMf)

    f(q,u) = cartpole_rule(q,(u=u, _p...,),t0)
    A = ForwardDiff.jacobian(x -> [f(x,1.0)...],q0) # t = 0.0 for Time-invariant systems
    B = ForwardDiff.derivative(x -> [f(q0,x)...],1.0)
    C = I
    D = 0I
    
    @inbounds function LTI(x, p, t)
        return p.A * x + p.B .* p.u([],t)
    end
    
    return ODEProblem(LTI,q0,(0.0,1.0),(A=A,B=B,C=C,D=D,u=cp.p.u))
end

function closed_loop_lqr(
    cp::ODEProblem = cartpole();
    lin_cp::ODEProblem = linear_cartpole(),
    x_ref::Function= t -> zero.([1.0,1.0,1.0,1.0]),
    Q = 1e1I,
    R = diagm([5e5,1e4,1e3,1e3]),
    W = zero
    )
    A,B,C,D,u = lin_cp.p
    S,_,K = arec(A,B,Q,R)
    u_lqr(x,t) = -(K*(x - x_ref(t)))[1] # + W(t); # [179.112  -31.4023  -240.162  -46.9267]
    
    return u_lqr
end


function trajectory_opt(cp::ODEProblem=cartpole();
    t=(0.0,2.0),
    x=([0.95,0,π,0],[0,0,0,0]),
    N=200)
    Tk = range(t...,N)
    # Uk = zeros(N)
    d = size(cp.u0,1)
    model = Model(Ipopt.Optimizer);
    
    q = [@variable(model, [1:d]) for _ in 1:N]
    u = [@variable(model, [1:1]) for _ in 1:N]
    
    @inbounds function J(u,Tk,N)
        J = 0.0
        for k in 1:N-1
            J += 0.5(Tk[k+1]-Tk[k])*(u[k+1][1]^2+u[k][1]^2)
        end
        return J
    end    
    @objective(model, Min, J(u,Tk,N))
    
    @inbounds function h(cp, uk, ukk, xk, xkk,tk,tkk)
        p = (M=cp.p.M, m=cp.p.m, l=cp.p.l, g=cp.p.g, Fc=cp.p.Fc, dMf=cp.p.dMf)
        fk = [cartpole_rule(xk,(u=uk[1],p...),tk)...]
        fkk = [cartpole_rule(xkk,(u=ukk[1],p...),tk)...]
        
        return 0.5(tkk-tk)*(fkk+fk) - (xkk - xk)
    end
    @constraint(model, [k = 1:N-1], 0 .== h(cp,u[k],u[k+1],q[k],q[k+1],Tk[k],Tk[k+1]))
    @constraint(model, [k = 1:N], -119 <= u[k][1] <= 119 )
    @constraint(model, [k = 1:N], -1 <= q[k][1] <= 1 )
    fix.(q[1],x[1]; force = true)
    fix.(q[N],x[2]; force = true)
    
    optimize!(model)
    
    H(t) = .5 - .5tanh(1e6t)
    u_opt(tt) = LinearInterpolation(
        [value.(uk)[1] for uk in u], Tk; 
        extrapolate = true)(tt)*H(tt-t[2])
    return u_opt
end

function swingup_and_stabilization()
    t_opt = (0.0,3.0)
    cp = cartpole(SVector(.95,0.0,π,0.0))
    u_lqr = closed_loop_lqr(cp)
    u_opt = trajectory_opt(cp; t=t_opt)

    _p = (M=cp.p.M, m=cp.p.m, l=cp.p.l, g=cp.p.g, Fc=cp.p.Fc, dMf=cp.p.dMf)

    H(t) = .5 - .5tanh(1e6t)
    T(t,th) = -H(t+th) + (1 - H(-t+th))
    u(x,t) = (1-T(x[3],.5))*u_opt(t) + (T(x[3],.5))*u_lqr(x,t)

    function f(x::SVector,_,t)
        return SVector(cartpole_rule(x,(u=u(x,t),_p...,),t)...)
    end

    function f(x::Vector,p,t)
        u = p.u
        return [cartpole_rule(x,(u=u(x,t),_p...,),t)...]
    end
    return ODEProblem(f,cp.u0,(0.0,10.0),(u=u,_p...))
end