import LinearAlgebra

#=
Our goal is to get rid of these functions, but they 
are helpful until we can get that one.
=# 

_dgemv_blas!(
  trans::Char,
  m::Int,
  n::Int,
  alpha::Float64,
  a::StridedVecOrMat{Float64},
  lda::Int,
  x::StridedVecOrMat{Float64},
  beta::Float64,
  y::StridedVecOrMat{Float64},
  ) =
  ccall((LinearAlgebra.BLAS.@blasfunc("dgemv_"), LinearAlgebra.BLAS.libblas), Cvoid,
   (Ref{UInt8}, Ref{LinearAlgebra.BlasInt}, Ref{LinearAlgebra.BlasInt},
    Ref{Float64}, Ptr{Float64}, Ref{LinearAlgebra.BlasInt},
      Ptr{Float64}, Ref{LinearAlgebra.BlasInt},
      Ref{Float64}, # beta
       Ptr{Float64}, Ref{LinearAlgebra.BlasInt}),
   trans, m, n, alpha, a, lda, x, stride(x,1), beta, y, stride(y,1))

