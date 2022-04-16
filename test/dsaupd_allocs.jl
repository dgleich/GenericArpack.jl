# This code came from helpers/dsaupd_vs_arpackjll_perf 
# but we use it here to test allocs
@testset "dsaupd allocations" begin 
  function eigrunallocs(op,ido, ::Val{BMAT}, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state; idonow=false)  where BMAT
    niter = 0 
    nbytes = 0 

    if idonow # just try everything in one go! 
      nbytes += @allocated GenericArpack.dsaupd!(ido, Val(BMAT), n, which, nev, tol, resid, ncv, V, ldv, iparam,
        ipntr, workd, workl, lworkl, info_initv;
        state, idonow = op
      )
      if ido[] != 99
        @warn("eigrun with idonow gave unfinished result")
      end
      return niter, nbytes
    end

    while ido[] != 99
      nbytes += @allocated GenericArpack.dsaupd!(ido, Val(BMAT), n, which, nev, tol, resid, ncv, V, ldv, iparam,
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
    return niter, nbytes
  end 

  @testset "dsaupd with reverse and idownow" begin
    op = GenericArpack.ArpackSimpleOp(Diagonal(1.0:10^3))
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
    state = GenericArpack.ArpackState{Float64}()

    niter = 0 

    eigrunallocs(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state)

    niter, nbytes = eigrunallocs(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state)
    @test nbytes == 0

    eigrunallocs(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state; idonow=true)

    niter, nbytes = eigrunallocs(op, ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state; idonow=true)
    @test nbytes == 0
  end
end 