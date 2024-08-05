using DifferentialEquations
using StaticArrays
using GLMakie
using Makie.GeometryBasics
using DataStructures: CircularBuffer

include("cartpole.jl")
cp = cartpole();

# Collocation methods