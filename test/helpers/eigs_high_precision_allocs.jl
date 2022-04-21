using LinearAlgebra
using Profile # to allow reseting allocation to get your target function 
using DoubleFloats, GenericArpack
using MultiFloats
using Quadmath
using StableRNGs

GenericArpack.@fix_doublefloats
T = Float128

rng = StableRNG(1)
A = Symmetric(randn(rng, 100,100))
eigs(T, A, 2; ritzvec=false)
Profile.clear_malloc_data()
eigs(T, A, 2; ritzvec=false)

##
include("allocations.jl")
lines = report_allocations(@__FILE__, system=true)
println.(lines);