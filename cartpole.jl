using DifferentialEquations
using StaticArrays
using ForwardDiff
using LinearAlgebra
using MatrixEquations


function cartpole(q0 = SVector([0.0, 0.0, 0.0, 0.0]...);
    u::Union{Number,Function} = zero, 
    M = 4.8, 
    m = 0.356, 
    l = 0.56, 
    g = 9.81,
    Fc = 4.9,
    dMf = 0.035)
    return ODEProblem(cartpole_rule, q0,(0.0,1.0),(u, M, m, l, g, Fc, dMf))
end
@inbounds function cartpole_rule(
                            q::SVector,
                            p::NTuple{7,Number},
                            t::Number)
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

    # du[1], du[2], du[3], du[4] = ẋ, ẍ, θ̇, θ̈
    return SVector{4}(ẋ, ẍ, θ̇, θ̈)
end

@inbounds function cartpole_rule(
    q::Vector,
    p::NTuple{7,Number},
    t::Number)

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

    # du[1], du[2], du[3], du[4] = ẋ, ẍ, θ̇, θ̈
    return [ẋ, ẍ, θ̇, θ̈]
end

@inbounds function cartpole_rule(
                            q::AbstractArray,
                            p::Tuple{Function,Number,Number,Number,Number,Number,Number},
                            t::Number)
    u, M, m, l, g, Fc, dMf = p
    return cartpole_rule(q,(u(t), M, m, l, g, Fc, dMf),t)
end

function linear_cartpole(
                        cp::ODEProblem = cartpole(), 
                        x0::AbstractArray = SVector{4}(0.0, 0.0, 0.0, 0.0))
    f = cp.f.f
    u = cp.p[1]
    A = ForwardDiff.jacobian(x -> f(x,cp.p,0.0),x0) # t = 0.0 for Time-invariant systems
    B = ForwardDiff.derivative(x -> f(x0,(x,cp.p[2:end]...,),1.0),1.0)
    C = I
    D = 0I

    @inbounds function LTI(x, p, t)
        A,B,C,D,u = p
        return A * x + B .* u(t)
    end

    return ODEProblem(LTI,x0,(0.0,1.0),(A,B,C,D,u))
end

function cl_cartpole(
                    cp::ODEProblem,
                    lin_cp::ODEProblem,
                    Q,
                    R,
                    W
)
    A,B,C,D,u = lin_cp.p
    S,_,K = arec(A,B,Q,R)
    
    @inbounds function cl_LTI(x, p, t)
        u, M, m, l, g, Fc, dMf = cp.p
        u =  -(K*x)[1] + W(t); # [179.112  -31.4023  -240.162  -46.9267]
        return cp.f.f(x, (u, M, m, l, g, Fc, dMf),t)
    end
    return ODEProblem(cl_LTI, cp.u0, cp.tspan,[])
end