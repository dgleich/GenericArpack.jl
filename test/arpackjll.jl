import Arpack_jll, LinearAlgebra

# show the library if we are debugging...
@debug Arpack_jll.libarpack

function arpack_set_debug_high()
  # [Documentation and structure of debug block here.](https://github.com/opencollab/arpack-ng/blob/master/SRC/debug.h)
  arpack_debug = cglobal((:debug_, Arpack_jll.libarpack), Int64)
  unsafe_store!(arpack_debug, 6, 1) # logfile set logfil to 6, the default stdout
  unsafe_store!(arpack_debug, -14, 2) # ndigit - use 14 digits of precision (ndigit)
  unsafe_store!.(arpack_debug, 4, 3:24) # turn on most debugging
end

function arpack_set_debug_low()
  # [Documentation and structure of debug block here.](https://github.com/opencollab/arpack-ng/blob/master/SRC/debug.h)
  arpack_debug = cglobal((:debug_, Arpack_jll.libarpack), Int64)
  unsafe_store!(arpack_debug, 6, 1) # logfile set logfil to 6, the default stdout
  unsafe_store!(arpack_debug, -6, 2) # ndigit - use 14 digits of precision (ndigit)
  unsafe_store!.(arpack_debug, 1, 3:24) # turn on most debugging
end


function arpack_dsconv(n::Int, ritz::Vector{Float64}, bounds::Vector{Float64},
                      tol::Float64)
  nconv = Ref{LinearAlgebra.BlasInt}(-1)
  ccall((:dsconv_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ptr{Float64},
     Ref{Float64},
     Ref{LinearAlgebra.BlasInt}),
    n, ritz, bounds, tol, nconv)
  return nconv[]
end

function arpack_dsortr(
  which::Symbol, # Input
  apply::Bool, # Input
  n::Int, # Input
  x1::Vector{Float64}, # Input/Output
  x2::Vector{Float64}, # Input/Output
  )
  whichstr = string(which)
  ccall((:dsortr_, Arpack_jll.libarpack), Cvoid,
    (Ptr{UInt8},
     Ref{LinearAlgebra.BlasInt}, # bool
     Ref{LinearAlgebra.BlasInt}, # size
     Ptr{Float64},
     Ref{Float64}),
    whichstr, apply, n, x1, x2)
end

function arpack_dsgets(
  ishift::Int, # Input
  which::Symbol, # Input
  kev::Int, # Input
  np::Int, # Input
  ritz::Vector{Float64}, # Input/Output
  bounds::Vector{Float64}, # Input/Output
  shifts::Vector{Float64} # Input/Output
  )
  ccall((:dsgets_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt},
     Ptr{UInt8},
     Ref{LinearAlgebra.BlasInt},
     Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ptr{Float64},
     Ptr{Float64}),
   ishift, string(which), kev, np, ritz, bounds, shifts)
end

function arpack_dstqrb!(n::Int, d::Vector{Float64}, e::Vector{Float64},
                      z::Vector{Float64}, work::Vector{Float64})
  info = Ref{LinearAlgebra.BlasInt}(-1)
  ccall((:dstqrb_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ptr{Float64},
     Ptr{Float64},
     Ptr{Float64},
     Ref{LinearAlgebra.BlasInt}),
    n, d, e, z, work, info)
  return info[]
end

##
import Arpack_jll, LinearAlgebra
function arpack_dgetv0!(ido::Ref{LinearAlgebra.BlasInt}, bmat::Symbol, itry::Int, initv::Bool,
  n::Int, j::Int, v::StridedMatrix{Float64}, ldv::Int,
    resid::StridedVecOrMat{Float64}, rnorm::Ref{Float64},
    ipntr::StridedVecOrMat{LinearAlgebra.BlasInt}, workd::StridedVecOrMat{Float64})
  ierr = Ref{LinearAlgebra.BlasInt}(0)
  # NOTE, arpack doesn't touch ierr unless ido[] == 0 or there is
  # a restart failure.
  # @show ido[], itry, initv, n, j, string(bmat)
  ccall((:dgetv0_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt},
     Ptr{UInt8},
     Ref{LinearAlgebra.BlasInt},
     Ref{LinearAlgebra.BlasInt}, # really a logical...
     Ref{LinearAlgebra.BlasInt},
     Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ref{Float64},
     Ptr{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ref{LinearAlgebra.BlasInt}),
    ido, string(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd, ierr)
  return ierr[]
end
