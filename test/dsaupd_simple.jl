using LinearAlgebra

@testset "diagonal test" begin 
  function mysimpleeigvals(op::ArpackOp, nev::Int = 6)
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

    iparam0 = copy(iparam)
    info_initv = 0

    T = Float64
    histdata = Vector{
      NamedTuple{(:ido, :resid, :V, :iparam, :ipntr, :workd, :workl, :ierr), 
        Tuple{Int,Vector{T}, Matrix{T}, Vector{Int}, Vector{Int}, Matrix{T}, Vector{T}, Int}}
    }()

    # Note that we cannot run two sequences at once and check them where we start a whole
    # second arpack call because of the expected Arpack state. 
    state = ArpackInJulia.ArpackState{Float64}()
    stats = ArpackStats()
    debug = ArpackInJulia.ArpackDebug(logfile=IOBuffer())
    ArpackInJulia.set_debug_high!(debug)
    while ido[] != 99
      ierr, state = ArpackInJulia.dsaupd!(ido, Val(bmat), n, which, nev, tol, resid, ncv, V, ldv, iparam,
        ipntr, workd, workl, lworkl, info_initv;
        state, stats, debug 
      )
      #ierr = arpack_dsaupd!(ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, 
      #  ipntr, workd, workl, lworkl, info_initv)
      # iparam9..11 are stats that aren't tracked the same... 
      push!(histdata, (;ido=ido[], resid=copy(resid), V=copy(V), iparam=copy(iparam), 
        ipntr=copy(ipntr), workd=copy(workd), workl=copy(workl), ierr))

      if ido[] == 1 || ido[] == -1
        ArpackInJulia._i_do_now_opx_1!(op, ipntr, workd, n)
      elseif ido[] == 99
        break
      else
        @error("this only supports standard eigenvalue problems")
      end 
    end
    # make sure we got the ritz values right...
    @test sort(eigvals(op.A), by=abs, rev=true)[1:nev] ≈ sort(workl[2*ncv+1:2*ncv+nev], by=abs, rev=true)
  end

  mysimpleeigvals(ArpackInJulia.ArpackSimpleOp(Diagonal(1.0:10)))
end 