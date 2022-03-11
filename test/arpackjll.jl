import Arpack_jll, LinearAlgebra, SHA

# show the library if we are debugging...
@debug Arpack_jll.libarpack, bytes2hex(open(SHA.sha256, Arpack_jll.libarpack))

function _reset_libarpack_dgetv0_iseed()
  # we need hardcoded offsets because we can't get these from cglobal :( )
  sha = open(SHA.sha256, Arpack_jll.libarpack)
  if bytes2hex(sha) == "5c1e951fad68bd7b180b83f2d821324efbdda7e00f5926698fad720249d6ac3f"
    iseedoffset = 0x0000000000059fc0
    dgetv0offset = 0x000000000001f940
  else
    errmsg = """Unknown libarpack sha $sha for $(Arpack_jll.libarpack)

    Please post a new issue on the github page for ArpackInJulia."""
    @error(errmsg)
  end
  libar = Base.Libc.dlopen(Arpack_jll.libarpack)
  dgetv0_real_offset = Base.Libc.dlsym(libar, "dgetv0_")
  base_offset = dgetv0_real_offset-dgetv0offset # this comes from the command above
  piseedoffset = Ptr{LinearAlgebra.BlasInt}(base_offset + iseedoffset)
  previseed = unsafe_load.(piseedoffset, (1,2,3,4))


  # store the values 1,3,5,6 at indices 1,2,3,4... which resets to the initial Arpack config.
  unsafe_store!.(piseedoffset, (1,3,5,7), (1,2,3,4))
  newiseed = unsafe_load.(piseedoffset, (1,2,3,4))

  return previseed, newiseed
end 

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

function arpack_dseigt!(rnorm::T,
  n::Int, H::AbstractMatrix{T}, ldh::Int,
  eig::AbstractVecOrMat{T}, bounds::AbstractVecOrMat{T}, 
  workl::AbstractVecOrMat{T}
) where T <: Float64 
  info = Ref{LinearAlgebra.BlasInt}(-1)
  ccall((:dseigt_, Arpack_jll.libarpack), Cvoid,
    (Ref{Float64},
      Ref{LinearAlgebra.BlasInt},
      Ptr{Float64},
      Ref{LinearAlgebra.BlasInt},
      Ptr{Float64},
      Ptr{Float64},
      Ptr{Float64},
      Ref{LinearAlgebra.BlasInt}),
    rnorm, n, H, ldh, eig, bounds, workl, info
  )
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
     Ref{LinearAlgebra.BlasInt},
     Int),
    ido, string(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd, ierr, 1)
  return ierr[]
end


##
import Arpack_jll, LinearAlgebra
function arpack_dsaitr!(
    ido::Ref{LinearAlgebra.BlasInt}, 
    bmat::Symbol, 
    n::Int, 
    k::Int,
    np::Int, 
    mode::Int, 
    resid::StridedVecOrMat{Float64},
    rnorm::Ref{Float64},
    v::StridedVecOrMat{Float64},
    ldv::Int,
    h::StridedVecOrMat{Float64},
    ldh::Int,
    ipntr::StridedVecOrMat{LinearAlgebra.BlasInt},
    workd::StridedVecOrMat{Float64}
)
  info = Ref{LinearAlgebra.BlasInt}(0)
  # NOTE, arpack doesn't touch ierr unless ido[] == 0 or there is
  # a restart failure.
  # @show ido[], itry, initv, n, j, string(bmat)
  ccall((:dsaitr_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt}, # ido
     Ptr{UInt8}, # bmat
     Ref{LinearAlgebra.BlasInt}, #n
     Ref{LinearAlgebra.BlasInt}, #k 
     Ref{LinearAlgebra.BlasInt}, #np
     Ref{LinearAlgebra.BlasInt}, #mode 
     Ptr{Float64}, #resid
     Ref{Float64}, #rnorm
     Ptr{Float64}, #v
     Ref{LinearAlgebra.BlasInt}, #ldv
     Ptr{Float64}, #h
     Ref{LinearAlgebra.BlasInt}, #lhd
     Ptr{LinearAlgebra.BlasInt}, #ipntr
     Ptr{Float64}, #workd
     Ref{LinearAlgebra.BlasInt}, Int), #info
    ido, string(bmat), n, k, np, mode, resid, rnorm, 
    v, ldv, h, ldh, ipntr, workd, info, 1)
  return info[]
end


##
import Arpack_jll, LinearAlgebra
function arpack_dsaup2!(
  ido::Ref{Int}, 
  bmat::Symbol,
  n::Int,
  which::Symbol,
  nev::Ref{Int},
  np::Ref{Int}, 
  tol::Float64,
  resid::StridedVecOrMat{Float64},
  mode::Int, 
  iupd::Int,
  ishift::Int,
  mxiter::Ref{Int},
  V::StridedMatrix{Float64},
  ldv::Int, 
  H::StridedMatrix{Float64}, 
  ldh::Int,
  ritz::StridedVecOrMat{Float64},
  bounds::StridedVecOrMat{Float64},
  Q::StridedMatrix{Float64},
  ldq::Int, 
  workl::StridedVecOrMat{Float64},
  ipntr::StridedVecOrMat{Int},
  workd::StridedVecOrMat{Float64},
  info_initv0::Int, # info in Arpack, but we return info... 
)
  info = Ref{LinearAlgebra.BlasInt}(info_initv0)
  ccall((:dsaup2_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt}, # ido
     Ptr{UInt8}, # bmat
     Ref{LinearAlgebra.BlasInt}, #n
     Ptr{UInt8}, # which
     Ref{LinearAlgebra.BlasInt}, # nev
     Ref{LinearAlgebra.BlasInt}, # np
     Ref{Float64}, # tol
     Ptr{Float64}, # resid
     Ref{LinearAlgebra.BlasInt}, # mode
     Ref{LinearAlgebra.BlasInt}, # iupd
     Ref{LinearAlgebra.BlasInt}, # ishift
     Ref{LinearAlgebra.BlasInt}, # mxiter
     Ptr{Float64}, # V
     Ref{LinearAlgebra.BlasInt}, # ldv
     Ptr{Float64}, # H
     Ref{LinearAlgebra.BlasInt}, # ldh
     Ptr{Float64}, # ritz
     Ptr{Float64}, # bounds 
     Ptr{Float64}, # Q
     Ref{LinearAlgebra.BlasInt}, # ldq
     Ptr{Float64}, # workl
     Ptr{LinearAlgebra.BlasInt}, # ipntr
     Ptr{Float64}, # workd
     Ref{LinearAlgebra.BlasInt},
     Int, Int), #info
    ido, string(bmat), n, string(which), nev, np, tol, 
    resid, 
    mode, iupd, ishift, mxiter, 
    V, ldv, H, ldh, ritz, bounds, Q, ldq, 
    workl, ipntr, workd, info, 1, 2)
  return info[]
end


##
import Arpack_jll, LinearAlgebra
function arpack_dsaupd!(
  ido::Ref{Int}, 
  bmat::Symbol,
  n::Int,
  which::Symbol,
  nev::Int,
  tol::Float64,
  resid::StridedVecOrMat{Float64},
  ncv::Int, 
  V::StridedMatrix{Float64},
  ldv::Int,
  iparam::StridedVecOrMat{Int},
  ipntr::StridedVecOrMat{Int},
  workd::StridedVecOrMat{Float64},
  workl::StridedVecOrMat{Float64},
  lworkl::Int,   
  info_initv0::Int, # info in Arpack, but we return info... 
)
  info = Ref{LinearAlgebra.BlasInt}(info_initv0)
  ccall((:dsaupd_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt}, # ido
     Ptr{UInt8}, # bmat
     Ref{LinearAlgebra.BlasInt}, #n
     Ptr{UInt8}, # which
     Ref{LinearAlgebra.BlasInt}, # nev
     Ref{Float64}, # tol
     Ptr{Float64}, # resid
     Ref{LinearAlgebra.BlasInt}, # ncv
     Ptr{Float64}, # V
     Ref{LinearAlgebra.BlasInt}, # ldv
     Ptr{LinearAlgebra.BlasInt}, # iparam
     Ptr{LinearAlgebra.BlasInt}, # ipntr
     Ptr{Float64}, # workd
     Ptr{Float64}, # workl
     Ref{LinearAlgebra.BlasInt}, # lworkl
     Ref{LinearAlgebra.BlasInt}, # info
     Int, Int), #info
    ido, string(bmat), n, string(which), nev, tol,
    resid, ncv, 
    V, ldv, iparam, ipntr, workd, workl, lworkl, info, 1, 2)
  return info[]
end