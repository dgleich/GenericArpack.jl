## test allocations
using BenchmarkTools
using GenericArpack
using LinearAlgebra
function eigrun(op,ido, ::Val{BMAT}, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state) where BMAT 
  niter = 0 
  while ido[] != 99
    GenericArpack.dsaupd!(ido, Val(BMAT), n, which, nev, tol, resid, ncv, V, ldv, iparam,
      ipntr, workd, workl, lworkl, info_initv;
      state 
    )
    if ido[] == 1 || ido[] == -1
      niter += 1
      GenericArpack._i_do_now_opx_1!(op, ipntr, workd, n)
    elseif ido[] == 99
      break
    else
      @error("this only supports standard eigenvalue problems")
    end 
  end
  return niter
end 

begin 
  @btime begin
    eigrun(op, ido, Val(bmat), n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state);
  end setup=begin
    op = GenericArpack.ArpackSimpleOp(Diagonal(1.0:10^3))
    nev = 6
    ido = Ref{Int}(0)
    bmat = :I
    n = size(op.A,1)
    which = :LM
    tol = 0.0 # just use the default
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
    workd = zeros(n,3)
    lworkl = ncv*ncv + 8*ncv
    workl = zeros(lworkl)

    info_initv = 0

    # Note that we cannot run two sequences at once and check them where we start a whole
    # second arpack call because of the expected Arpack state. 
    state = GenericArpack.ArpackState{Float64}()

    niter = 0 
  end
end



## test allocations
using BenchmarkTools
using GenericArpack
using LinearAlgebra
begin 
  @btime begin
    for i=1:10
      ierr, state = GenericArpack.dsaupd!(ido, Val(bmat), n, which, nev, tol, resid, ncv, V, ldv, iparam,
        ipntr, workd, workl, lworkl, info_initv;
        state 
      )
    end
    
  end setup=begin
    op = GenericArpack.ArpackSimpleOp(Diagonal(1.0:10))
    nev = 6
    ido = Ref{Int}(0)
    bmat = :I
    n = size(op.A,1)
    which = :LM
    tol = 0.0 # just use the default
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
    workd = zeros(n,3)
    lworkl = ncv*ncv + 8*ncv
    workl = zeros(lworkl)

    info_initv = 0

    # Note that we cannot run two sequences at once and check them where we start a whole
    # second arpack call because of the expected Arpack state. 
    state = GenericArpack.ArpackState{Float64}()

    niter = 0 
  end
end


## profiling
using Revise, GenericArpack, LinearAlgebra, BenchmarkTools
begin
  op = GenericArpack.ArpackSimpleOp(Diagonal(1.0:10^3))
  nev = 6
  ido = Ref{Int}(0)
  bmat = Val(:I)
  n = size(op.A,1)
  which = :LM
  tol = 0.0 # just use the default
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
  state = GenericArpack.ArpackState{Float64}()

  niter = 0 
  @profview eigrun(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state)
end