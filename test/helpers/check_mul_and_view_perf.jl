using LinearAlgebra
using BenchmarkTools
function mytransmul!(x::AbstractVector{T}, A::AbstractMatrix{T}, y::AbstractVector{T}) where T
  n = length(x)
  mul!(@view(y[1:n-1]), adjoint(@view(A[1:n-1,1:n-1])), @view(x[1:n-1]));
  #mul!(y, adjoint(A), x);
end 
@btime begin 
  mytransmul!(x, A, y);
end setup=begin
  n = 100
  A = ones(n,n)
  x = ones(n)
  y = zeros(n)
end

##

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


@btime begin
  _dgemv_blas!('T', n, n, 1.0, A, stride(A,2), x, 0.0, y)
end setup=begin
  n = 100
  A = ones(n,n)
  x = ones(n)
  y = zeros(n)
end