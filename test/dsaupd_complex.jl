@testset "complex hermitian simple test" begin 
  using LinearAlgebra
  function mysimpleeigvals(op::ArpackOp, nev::Int = 6)
    T = Complex{Float64}
    ido = Ref{Int}(0)
    bmat = :I
    n = size(op.A,1)
    which = :LM
    tol = 0.0 # just use the default
    resid = zeros(T, n)
    ncv = min(2nev, n-1)
    V = zeros(T, n,ncv)
    ldv = n
    mode = 1 
    iparam = zeros(Int,11)
    iparam[1] = 1
    #iparam[3] = 300 # max iteration
    iparam[3] = 300 # 
    iparam[4] = 1
    iparam[7] = mode 
    ipntr = zeros(Int,11)
    workd = zeros(T,3n)
    lworkl = ncv*ncv + 8*ncv
    workl = zeros(Float64,lworkl)

    iparam0 = copy(iparam)
    info_initv = 0

    histdata = Vector{
      NamedTuple{(:ido, :resid, :V, :iparam, :ipntr, :workd, :workl, :ierr), 
        Tuple{Int,Vector{T}, Matrix{T}, Vector{Int}, Vector{Int}, Vector{T}, Vector{T}, Int}}
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
    @test sort(eigvals(Matrix(op.A)), by=abs, rev=true)[1:nev] ≈ sort(workl[2*ncv+1:2*ncv+nev], by=abs, rev=true)

    # now test to make sure we get dseupd okay too
    Z = zeros(ComplexF64, n, nev)
    d = zeros(Float64, nev)
    select = zeros(Int, ncv)

    ierr = ArpackInJulia.simple_dseupd!(true, select, d, Z, 0.0, Val(bmat), n, which, nev, tol, resid, ncv, V, iparam, ipntr, workd, workl)

    @test sort(eigvals(Matrix(op.A)), by=abs, rev=true)[1:nev] ≈ sort(d, by=abs, rev=true)
    for i in 1:size(Z,2)
      lam = Z[:,i]'*op.A*Z[:,i]
      @test lam ≈ d[i]
      @test norm(Z[:,i]) ≈ 1
    end
    @test Z'*Z ≈ I
  end

  n = 20 
  A = Tridiagonal(1.0im*ones(n-1), collect(range(0.1, 2.0, length=n)) .+ 0.0*im, -1.0im*ones(n-1))
  mysimpleeigvals(ArpackInJulia.ArpackSimpleOp(A))
end
