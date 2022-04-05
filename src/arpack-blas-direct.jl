

_pointer_to_array_offset(a::StridedVector{T}, offset::Int) where T =
  Base.unsafe_convert(Ptr{T}, a) + (offset-1)*sizeof(T)*strides(a)[1]
_pointer_to_array_offset(a::StridedMatrix{T}, offset1::Int, offset2::Int) where T =
  Base.unsafe_convert(Ptr{T}, a) + (offset1-1)*sizeof(T)*strides(a)[1] + (offset2-1)*sizeof(T)*strides(a)[2]



#= raw BLAS and LINPACK calls in dstqrb mapped to native Julia and also LAPACK/BLAS =#
import LinearAlgebra
_dlascl_blas!(cfrom::Float64, cto::Float64, a::StridedVecOrMat{Float64}) = begin
  ccall((LinearAlgebra.BLAS.@blasfunc("dlascl_"), LinearAlgebra.BLAS.libblas), Cvoid,
    (Ptr{UInt8}, # character
      Ref{LinearAlgebra.BlasInt},  # lower bandwidth
      Ref{LinearAlgebra.BlasInt}, # upper bandwidth
      Ref{Float64}, Ref{Float64}, # cfrom, cto
      Ref{LinearAlgebra.BlasInt}, Ref{LinearAlgebra.BlasInt}, # m, n
      Ptr{Float64}, # a
      Ref{LinearAlgebra.BlasInt}, # lda
      Ref{LinearAlgebra.BlasInt}) # info
      ,
      "G",
      1, 1,
      cfrom, cto,
      size(a,1), size(a,2),
      a,
      ndims(a) > 1 ? strides(a)[2] : 1, Ref{LinearAlgebra.BlasInt}(-1))
  a
end


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

_dnrm2_blas(
  a::StridedVecOrMat{Float64},
) = ccall((LinearAlgebra.BLAS.@blasfunc("dnrm2_"), LinearAlgebra.BLAS.libblas), Float64,
  (Ref{LinearAlgebra.BlasInt}, Ptr{Float64}, Ref{LinearAlgebra.BlasInt}),
  length(a), a, stride(a,1))


##
function _dlapy2_blas(x::Float64, y::Float64)
  return ccall((LinearAlgebra.BLAS.@blasfunc("dlapy2_"), LinearAlgebra.BLAS.libblas), Float64,
    (Ref{Float64}, Ref{Float64}),
    x, y)
end

##
_dlaev2_blas(a::Float64, b::Float64, c::Float64) = begin
  rt1 = Ref{Float64}(0.0)
  rt2 = Ref{Float64}(0.0)
  cs = Ref{Float64}(0.0)
  sn = Ref{Float64}(0.0)
  ccall((LinearAlgebra.BLAS.@blasfunc("dlaev2_"), LinearAlgebra.BLAS.libblas), Cvoid,
    (Ref{Float64}, Ref{Float64}, Ref{Float64},
     Ref{Float64}, Ref{Float64}, Ref{Float64}, Ref{Float64}),
    a, b, c, rt1, rt2, cs, sn)
  return rt1[], rt2[], cs[], sn[]
end

_dlaev2_blas(1.0, 1.0, 1.0)
_dlaev2_blas(5.0, 0.0, 2.0)

##
_dlasr_rvb_blas!(m::LinearAlgebra.BlasInt,
                 n::LinearAlgebra.BlasInt,
                 c::Union{Ptr{Float64},StridedVector{Float64}},
                 s::Union{Ptr{Float64},StridedVector{Float64}},
                 a::StridedVecOrMat{Float64},
                 lda::LinearAlgebra.BlasInt) = begin
  ccall((LinearAlgebra.BLAS.@blasfunc("dlasr_"), LinearAlgebra.BLAS.libblas), Cvoid,
    (Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8},
      Ref{LinearAlgebra.BlasInt}, Ref{LinearAlgebra.BlasInt},
      Ptr{Float64}, Ptr{Float64},
      Ptr{Float64},
      Ref{LinearAlgebra.BlasInt}),
    "r", "v", "b", m, n, c, s, a, lda)
  return a
end

##

_dlasr_rvf_blas!(m::LinearAlgebra.BlasInt,
                 n::LinearAlgebra.BlasInt,
                 c::Union{Ptr{Float64},StridedVector{Float64}},
                 s::Union{Ptr{Float64},StridedVector{Float64}},
                 a::StridedVecOrMat{Float64},
                 lda::LinearAlgebra.BlasInt) = begin
  ccall((LinearAlgebra.BLAS.@blasfunc("dlasr_"), LinearAlgebra.BLAS.libblas), Cvoid,
    (Ptr{UInt8}, Ptr{UInt8}, Ptr{UInt8},
      Ref{LinearAlgebra.BlasInt}, Ref{LinearAlgebra.BlasInt},
      Ptr{Float64}, Ptr{Float64},
      Ptr{Float64},
      Ref{LinearAlgebra.BlasInt}),
    "r", "v", "f", m, n, c, s, a, lda)
  return a
end

##
_dlartg(f::Float64, g::Float64) = begin
  c = Ref{Float64}(0.0)
  s = Ref{Float64}(0.0)
  r = Ref{Float64}(0.0)
  ccall((LinearAlgebra.BLAS.@blasfunc("dlartg_"), LinearAlgebra.BLAS.libblas), Cvoid,
    (Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64},Ref{Float64}),
    f, g, c, s, r)
  return c[], s[], r[]
end

##
_iladlc_blas(A::StridedMatrix) = begin
  ccall((LinearAlgebra.BLAS.@blasfunc("iladlc_"), LinearAlgebra.BLAS.libblas), LinearAlgebra.BlasInt,
    (Ref{LinearAlgebra.BlasInt},Ref{LinearAlgebra.BlasInt},Ref{Float64},Ref{LinearAlgebra.BlasInt}),
    size(A,1),size(A,2),A,stride(A,2))
end

##
function _dgeqr2_blas!(A::StridedMatrix{Float64}, 
  tau::StridedVector{Float64}=zeros(minimum(size(A))), 
  work::StridedVector{Float64}=zeros(size(A,2))) 
  info = Ref{LinearAlgebra.BlasInt}(0) 
  ccall((LinearAlgebra.BLAS.@blasfunc("dgeqr2_"), LinearAlgebra.BLAS.libblas), 
    Cvoid, 
    (Ref{LinearAlgebra.BlasInt},Ref{LinearAlgebra.BlasInt}, Ptr{Float64}, # m, n, A
    Ref{LinearAlgebra.BlasInt},  # lda 
    Ptr{Float64}, # tau 
    Ptr{Float64}, # work 
    Ref{LinearAlgebra.BlasInt}), # info
    size(A,1),size(A,2),A,stride(A,2),
    tau, work, info
    )
  return A, tau
end

##
_dorm2r_blas!(side::Char, trans::Char,
m::Int, n::Int, k::Int, A::StridedArray{Float64}, 
tau::StridedVector{Float64}, C::StridedMatrix{Float64},
work = zeros(maximum(size(A)))
) = begin
  info = Ref{LinearAlgebra.BlasInt}(0) 
  ccall((LinearAlgebra.BLAS.@blasfunc("dorm2r_"), LinearAlgebra.BLAS.libblas), 
    Cvoid, 
    (Ref{UInt8}, Ref{UInt8}, 
    Ref{LinearAlgebra.BlasInt},Ref{LinearAlgebra.BlasInt}, Ref{LinearAlgebra.BlasInt}, # m n k 
    Ptr{Float64}, # A, 
    Ref{LinearAlgebra.BlasInt},  # lda 
    Ptr{Float64}, # tau 
    Ptr{Float64}, # C
    Ref{LinearAlgebra.BlasInt},  # ldc 
    Ptr{Float64}, # work 
    Ref{LinearAlgebra.BlasInt}), # info, 1, 1 for character sizes
    side, trans, m, n, k, A, stride(A,2), tau, C, stride(C,2), work, info
    )
  return C
end

##
using LinearAlgebra: BLAS, BlasInt, LinearAlgebra
_dlarnv_blas!(idist::Int,
       iseed::Ref{NTuple{4,Int}},
      n::Int,
      x::Vector{Float64}) =
  ccall((LinearAlgebra.BLAS.@blasfunc("dlarnv_"), LinearAlgebra.BLAS.libblas), Cvoid,
    (Ref{LinearAlgebra.BlasInt}, Ptr{LinearAlgebra.BlasInt}, Ref{LinearAlgebra.BlasInt}, Ptr{Float64}),
    idist, iseed, n, x)#


##
function _dsteqr_blas!(d::StridedVecOrMat{Float64}, e::StridedVecOrMat{Float64}, 
  Z::StridedVecOrMat{Float64}, work::StridedVecOrMat{Float64})
  info = Ref{LinearAlgebra.BlasInt}(0) 
  ccall((LinearAlgebra.BLAS.@blasfunc("dsteqr_"), LinearAlgebra.BLAS.libblas), 
    Cvoid, 
    (Ref{UInt8}, # ICOMPZ
    Ref{LinearAlgebra.BlasInt}, # n 
    Ptr{Float64}, # d, 
    Ptr{Float64}, # e, 
    Ptr{Float64}, # z, 
    Ref{LinearAlgebra.BlasInt}, # ldz
    Ptr{Float64}, # work, 
    Ref{LinearAlgebra.BlasInt}, Int), # info, 1 for character sizes
    'I', length(d), d, e, Z, stride(Z,2), work, info, 1
    )
  if info[] != 0 
    error("dsteqr gave info $(info[])")
  end
  return Z 
end 