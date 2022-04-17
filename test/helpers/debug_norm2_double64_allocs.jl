using LinearAlgebra
using Profile # to allow reseting allocation to get your target function 
using DoubleFloats, GenericArpack

a = randn(Double64, 10000)
GenericArpack.norm2(a)
Profile.clear_malloc_data()
GenericArpack.norm2(a)

##
include("allocations.jl")
lines = report_allocations(@__FILE__, system=true)
println.(lines);