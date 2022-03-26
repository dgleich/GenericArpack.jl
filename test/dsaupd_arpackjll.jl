@testset "dsaupd compare against arpackjll" begin 
  function _check_dsaupd_sequence!(M;
    B=1.0LinearAlgebra.I,
    nev=min(size(M,1)-2,6),
    ncv=min(2nev,size(M,1)-1),
    mode::Int=1,
    initv=nothing,
    bmat=:I,
    which=:LM,
    maxiter=300,
    tol=0.0, # use the default 
  )
    n = size(M,1)
    
    resid = zeros(n)
    info_initv = 0 
    if initv !== nothing
      copyto!(resid, initv)
      info_initv = 1 
    end 

    # allocations... 
    ido = Ref{Int}(0)
    V = zeros(n,ncv)
    ldv = n
    iparam = zeros(Int,11)
    iparam[1] = 1
    iparam[3] = maxiter # maxiter
    iparam[4] = 1 # not used 
    iparam[7] = mode 
    ipntr = zeros(Int,11)
    workd = zeros(n,3)
    lworkl = ncv*ncv + 8*ncv
    workl = zeros(lworkl)

    state = ArpackInJulia.ArpackState{Float64}()
    stats = ArpackInJulia.ArpackStats()
    iter = 0
    debug = ArpackInJulia.ArpackDebug()
    #debug.maitr = 1
    while ido[] != 99
      # make a copy of state for arpack
      arido = Ref(ido[])
      arV = copy(V)
      ariparam = copy(iparam)
      aripntr = copy(ipntr)
      arresid = copy(resid)
      arworkd = copy(workd)
      arworkl = copy(workl)

      ierr, state = ArpackInJulia.dsaupd!(ido, Val(bmat), n, which, nev, tol, resid, ncv, V, ldv, iparam,
        ipntr, workd, workl, lworkl, info_initv;
        state, stats, debug
      )
      #ierr = arpack_dsaupd!(ido, bmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, 
      #  ipntr, workd, workl, lworkl, info_initv)
      # iparam9..11 are stats that aren't tracked the same... 
      arierr = arpack_dsaupd!(arido, bmat, n, which, nev, tol, arresid, ncv, arV, ldv, ariparam, 
        aripntr, arworkd, arworkl, lworkl, info_initv)

      # set info_initv to zero after the first iteration to
      # make Arpack work...
      info_initv = 0 

      @test ido[] == arido[]
      @test workl == arworkl 
      @test resid == arresid
      @test ierr == arierr
      @test V == arV 
      @test ipntr == aripntr
      @test iparam == ariparam
      @test workd == arworkd 
 
      if ido[] == -1 || ido[] == 1
        if mode == 2
          # crazy interface, see remark 5 
          mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
          copyto!(@view(workd[ipntr[1]:ipntr[1]+n-1]),@view(workd[ipntr[2]:ipntr[2]+n-1]))
          ldiv!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
        else
          #mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
          ArpackInJulia._i_do_now_opx_1!(ArpackInJulia.ArpackSimpleOp(M), ipntr, workd, n)
        end 
      elseif ido[] == 2
        mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
      elseif ido[] == 99 
        # we are done... will exit...
      else
        @error("Wrong ido, $(ido[])")
      end
    end
    return nothing
  end


  @testset "dsaitr compare arpackjll simple call" begin
    using LinearAlgebra
    bmat = :I
    n = 10 
    nev = 3
    mode = 1
    resid = ones(n)/sqrt(n)

    M = Diagonal(1.0:n)
    
    @testset "run 1" begin 
      seqdata = _check_dsaupd_sequence!(M;bmat, nev, initv=resid, mode)
    end 

    @testset "run 2" begin 
      # make sure there aren't any weird duplicate scenarios. 
      seqdata = _check_dsaupd_sequence!(M;bmat, nev, initv=resid, mode)
    end 
  end 

  @testset "dsaitr compare arpackjll generalized call" begin
    using LinearAlgebra
    bmat = :G
    n = 10
    nev = 3
    ncv = 6
    mode = 2
    resid = ones(n)/sqrt(n)
    M = Diagonal(1.0:n)
    B = Diagonal(range(0.1, 1.0, length=n))
    bmat = :G

    @testset "run 1" begin 
      seqdata = _check_dsaupd_sequence!(M;B, bmat,nev,initv=resid, mode)
    end 

    # this one does need a random getv0
    _reset_libarpack_dgetv0_iseed()

    @testset "run 2" begin 
      # make sure there aren't any weird duplicate scenarios. 
      seqdata = _check_dsaupd_sequence!(M;B, bmat, nev,initv=resid, mode)
    end 
  end

  @testset "Generalized Eigenproblem with Diagonal that is too long" begin 
    using LinearAlgebra
    ido = Ref{Int}(0)
    bmat = :G
    n = 10
    ncv = 10 
    nev = 9 
    mode = 3

    A = Diagonal(1.0:n)
    B = Diagonal(collect(range(0.1, 1.0, length=n)))

    _reset_libarpack_dgetv0_iseed()
    _check_dsaupd_sequence!(A; B, mode, ncv, nev, bmat)
  end

  @testset "Generalized Eigenproblem with Diagonal that is too long" begin 
    using LinearAlgebra
    ido = Ref{Int}(0)
    bmat = :G
    n = 10
    ncv = 10 
    nev = 9 
    mode = 3

    A = Diagonal(1.0:n)
    B = Diagonal(collect(range(0.1, 1.0, length=n)))

    _reset_libarpack_dgetv0_iseed()
    _check_dsaupd_sequence!(A; B, mode, ncv, nev, bmat)
  end
end

#=
  


  

  @testset "Generalized Eigenproblem with Diagonal" begin 
    using LinearAlgebra
    using Random 
    Random.seed!(0)
    ido = Ref{Int}(0)
    bmat = :G
    n = 10
    k = 0 # number of current columns in V
    np = 9
    mode = 2
    resid = randn(n)

    A = Diagonal(1.0:n)
    B = Diagonal(collect(range(0.1, 1.0, length=n)))
    #M = inv(B)*A # okay since they are diagonal....

    V = zeros(n,k+np)
    ldv = n 
    H = zeros(n,2) # full h
    ldh = n 

    rnorm = Ref{Float64}(sqrt(abs(resid'*B*resid)))
    _reset_libarpack_dgetv0_iseed()
    seqdata = _check_dsaitr_sequence!(A; B, idostart=0, 
      bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
    )
    
  end 

  @testset "Generalized Eigenproblem with Diagonal that is too long" begin 
    using LinearAlgebra
    using Random 
    Random.seed!(0)
    ido = Ref{Int}(0)
    bmat = :G
    n = 10
    k = 0 # number of current columns in V
    np = 12
    mode = 3
    resid = randn(n)

    A = Diagonal(1.0:n)
    B = Diagonal(collect(range(0.1, 1.0, length=n)))
    #M = inv(B)*A # okay since they are diagonal....

    V = zeros(n,k+np)
    ldv = n 
    H = zeros(k+np,2) # full h
    ldh = k+np 

    rnorm = Ref{Float64}(sqrt(abs(resid'*B*resid)))
    _reset_libarpack_dgetv0_iseed()
    seqdata = _check_dsaitr_sequence!(A; B, idostart=0, 
      bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
    )
    
  end 
end 

=#

#=
Codes that were made obselete by new tests, but may still be useful if things go south.



@testset "dsaitr arpackjll call" begin
  using LinearAlgebra
  using Random 
  Random.seed!(0)
  ido = Ref{Int}(0)
  bmat = :I
  n = 10
  k = 0 # number of current columns in V
  np = 3
  #V = randn(n,k+np) # total memory for v
  #V[:,1:k] = qr(V).Q[:,1:k] # random orthogonal init
  mode = 1
  resid = randn(n)
  rnorm = Ref{Float64}(norm(resid))

  V = zeros(n,k+np)
  ldv = n 
  H = zeros(n,2) # full h
  ldh = n 
  ipntr = zeros(Int,3)
  workd = zeros(3n)

  arpack_set_debug_high()

  M = Diagonal(1.0:n)
  B = Diagonal(range(0.1, n, length=n))

  handle_ido = (ido, workd, ipntr) -> begin 
    if ido[] == 1 || ido[] == -1
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    elseif ido[] == 2
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    end
  end 

  resid0 = copy(resid)

  info = arpack_dsaitr!(
    ido, bmat, n, k, np, mode, resid, rnorm, 
    V, ldv, H, ldh, ipntr, workd
  )
  @show ido[]
  @test workd[ipntr[1]:ipntr[1]+n-1] ≈ resid0/norm(resid0)

  niter = 1
  while ido[] != 99 
    handle_ido(ido, workd, ipntr)
    niter += 1
    info = arpack_dsaitr!(
      ido, bmat, n, k, np, mode, resid, rnorm, 
      V, ldv, H, ldh, ipntr, workd
    )
    println("Call $niter, ido=$(ido[])")
  end 
end
=#