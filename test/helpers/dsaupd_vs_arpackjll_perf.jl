## test allocations
using BenchmarkTools
using ArpackInJulia
using LinearAlgebra
using Arpack_jll

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

function eigrun(op,ido, ::Val{BMAT}, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state; idonow=false)  where BMAT
  niter = 0 
  nbytes = 0 

  if idonow # just try everything in one go! 
    nbytes += @allocated ArpackInJulia.dsaupd!(ido, Val(BMAT), n, which, nev, tol, resid, ncv, V, ldv, iparam,
      ipntr, workd, workl, lworkl, info_initv;
      state, idonow = op
    )
    if ido[] != 99
      @warn("eigrun with idonow gave unfinished result")
    end
    return niter, nbytes
  end

  while ido[] != 99
    nbytes += @allocated ArpackInJulia.dsaupd!(ido, Val(BMAT), n, which, nev, tol, resid, ncv, V, ldv, iparam,
      ipntr, workd, workl, lworkl, info_initv;
      state 
    )
    if ido[] == 1 || ido[] == -1
      niter += 1
      ArpackInJulia._i_do_now_opx_1!(op, ipntr, workd, n)
    elseif ido[] == 99
      break
    else
      @error("this only supports standard eigenvalue problems")
    end 
  end
  return niter, nbytes
end 

function eigrun_arpackjll(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv)
  niter = 0 
  nbytes = 0 
  while ido[] != 99
    nbytes += @allocated arpack_dsaupd!(ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam,
      ipntr, workd, workl, lworkl, info_initv
    )
    
    info_initv = 0 # important for arpack... 
    if ido[] == 1 || ido[] == -1
      niter += 1
      ArpackInJulia._i_do_now_opx_1!(op, ipntr, workd, n)
    elseif ido[] == 99
      break
    else
      @error("this only supports standard eigenvalue problems")
    end 
  end
  return niter, nbytes
end 


##
println("Profiling arpack with reverse communcation")
begin 
  @btime begin
    eigrun(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state);
  end setup=begin
    op = ArpackInJulia.ArpackSimpleOp(Diagonal(1.0:10^3))
    nev = 6
    ido = Ref{Int}(0)
    bmat = Val(:I)
    n = size(op.A,1)
    which = :LM
    tol = eps(Float64)/2 # just use the default
    resid = zeros(n)
    ncv = min(2nev, n-1)
    V = zeros(n,ncv)
    ldv = n
    mode = 1 
    iparam = zeros(Int,11)
    iparam[1] = 1
    #iparam[3] = 300 # max iteration
    iparam[3] = 300 # 
    iparam[4] = 1
    iparam[7] = mode 
    ipntr = zeros(Int,11)
    workd = zeros(3n)
    lworkl = ncv*ncv + 8*ncv
    workl = zeros(lworkl)

    info_initv = 0

    # Note that we cannot run two sequences at once and check them where we start a whole
    # second arpack call because of the expected Arpack state. 
    state = ArpackInJulia.ArpackState{Float64}()

    niter = 0 
  end
end

##
println("Profiling arpack with idonow ops")
begin 
  @btime begin
    eigrun(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state; idonow=true);
  end setup=begin
    op = ArpackInJulia.ArpackSimpleOp(Diagonal(1.0:10^3))
    nev = 6
    ido = Ref{Int}(0)
    bmat = Val(:I)
    n = size(op.A,1)
    which = :LM
    tol = eps(Float64)/2 # just use the default
    resid = zeros(n)
    ncv = min(2nev, n-1)
    V = zeros(n,ncv)
    ldv = n
    mode = 1 
    iparam = zeros(Int,11)
    iparam[1] = 1
    #iparam[3] = 300 # max iteration
    iparam[3] = 300 # 
    iparam[4] = 1
    iparam[7] = mode 
    ipntr = zeros(Int,11)
    workd = zeros(3n)
    lworkl = ncv*ncv + 8*ncv
    workl = zeros(lworkl)

    info_initv = 0

    # Note that we cannot run two sequences at once and check them where we start a whole
    # second arpack call because of the expected Arpack state. 
    state = ArpackInJulia.ArpackState{Float64}()

    niter = 0 
  end
end

##
println("Profiling arpack with Arpack_jll")
begin 
  @btime begin
    eigrun_arpackjll(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv);
  end setup=begin
    op = ArpackInJulia.ArpackSimpleOp(Diagonal(1.0:10^3))
    nev = 6
    ido = Ref{Int}(0)
    bmat = :I
    n = size(op.A,1)
    which = :LM
    tol = eps(Float64)/2 # just use the default
    resid = zeros(n)
    ncv = min(2nev, n-1)
    V = zeros(n,ncv)
    ldv = n
    mode = 1 
    iparam = zeros(Int,11)
    iparam[1] = 1
    iparam[3] = 300 # 
    iparam[4] = 1
    iparam[7] = mode 
    ipntr = zeros(Int,11)
    workd = zeros(3n)
    lworkl = ncv*ncv + 8*ncv
    workl = zeros(lworkl)
    info_initv = 0
  end
end

## profiling 
using Revise, ArpackInJulia, LinearAlgebra, BenchmarkTools
begin
  op = ArpackInJulia.ArpackSimpleOp(Diagonal(1.0:10^3))
  nev = 6
  ido = Ref{Int}(0)
  bmat = Val(:I)
  n = size(op.A,1)
  which = :LM
  tol = eps(Float64)/2 # just use the default
  resid = zeros(n)
  ncv = min(2nev, n-1)
  V = zeros(n,ncv)
  ldv = n
  mode = 1 
  iparam = zeros(Int,11)
  iparam[1] = 1
  #iparam[3] = 300 # max iteration
  iparam[3] = 300 # 
  iparam[4] = 1
  iparam[7] = mode 
  ipntr = zeros(Int,11)
  workd = zeros(3n)
  lworkl = ncv*ncv + 8*ncv
  workl = zeros(lworkl)

  info_initv = 0

  # Note that we cannot run two sequences at once and check them where we start a whole
  # second arpack call because of the expected Arpack state. 
  state = ArpackInJulia.ArpackState{Float64}()

  niter = 0 
  @profview eigrun(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state)
end