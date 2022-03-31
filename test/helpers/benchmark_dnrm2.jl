using BenchmarkTools
@btime begin 
  ArpackInJulia._dnrm2_unroll_ext(x)
end setup=begin 
  x = randn(10000)
end 

##
using BenchmarkTools
using LinearAlgebra
@btime begin 
  norm(x)
end setup=begin 
  x = randn(10000)
end 

##
using LinearAlgebra
_dnrm2_blas(
  a::StridedVecOrMat{Float64},
) = ccall((LinearAlgebra.BLAS.@blasfunc("dnrm2_"), LinearAlgebra.BLAS.libblas), Float64,
  (Ref{LinearAlgebra.BlasInt}, Ptr{Float64}, Ref{LinearAlgebra.BlasInt}),
  length(a), a, stride(a,1))

using BenchmarkTools
@btime begin 
  _dnrm2_blas(x)
end setup=begin 
  x = randn(10000)
end   