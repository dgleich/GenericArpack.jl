@testset "dsaupd compare against arpackjll" begin 
  function _check_dsaupd_sequence!(M;
    B=1.0LinearAlgebra.I,
    invBop=B, 
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

    if tol <= 0 
      # need to force Tol for arpackjll
      tol = eps(Float64)/2
    end 

    state = GenericArpack.ArpackState{Float64}()
    stats = GenericArpack.ArpackStats()
    iter = 0
    debug = GenericArpack.ArpackDebug()
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

      ierr, state = GenericArpack.dsaupd!(ido, Val(bmat), n, which, nev, tol, resid, ncv, V, ldv, iparam,
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

      if true
        # add some debugging code
        # that can help track down failures...
        test_diffs("V on iter $(iter) ", V, arV)
        test_diffs("workl on iter $(iter) ", workl, arworkl)
        test_diffs("workd on iter $(iter) ", workd, arworkd)
      end 

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
          #ldiv!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
          ldiv!(invBop, @view(workd[ipntr[2]:ipntr[2]+n-1]))
        else
          #mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
          GenericArpack._i_do_now_opx_1!(GenericArpack.ArpackSimpleOp(M), ipntr, workd, n)
        end 
      elseif ido[] == 2
        mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
      elseif ido[] == 99 
        # we are done... will exit...
      else
        @error("Wrong ido, $(ido[])")
      end
      iter += 1
    end
    #@show iparam[5], workl[2*ncv+1:2*ncv+nev], stats.nopx
    return iter
  end


  @testset "dsaupd compare arpackjll simple call" begin
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


  @testset "dsaupd compare arpackjll simple call long diagonal" begin
    using LinearAlgebra
    bmat = :I
    n = 1000
    nev = 6
    mode = 1
    resid = ones(n)/sqrt(n)

    M = Diagonal(1.0:n)
    
    @testset "run 1" begin 
      seqdata = _check_dsaupd_sequence!(M;bmat, nev, initv=resid, mode)
      #println("Long diagonal run takes $(seqdata) iterations")
    end 
  end 

  @testset "dsaupd compare arpackjll generalized call" begin
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

    # this one does need a random getv0
    _reset_libarpack_dgetv0_iseed()

    @testset "run 1" begin 
      seqdata = _check_dsaupd_sequence!(M;B,invBop=B, bmat,nev,initv=resid, mode)
    end 

    # this one does need a random getv0
    _reset_libarpack_dgetv0_iseed()

    @testset "run 2" begin 
      # make sure there aren't any weird duplicate scenarios. 
      seqdata = _check_dsaupd_sequence!(M;B,invBop=B, bmat, nev,initv=resid, mode)
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
    _check_dsaupd_sequence!(A; B,invBop=B, mode, ncv, nev, bmat)
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
    _check_dsaupd_sequence!(A; B,invBop=B, mode, ncv, nev, bmat)
  end

  @testset "Generalized Eigenproblem from dsdrv3" begin 

    n = 100
    ncv = 10
    nev = 4
    A = Tridiagonal(-ones(n-1),2*ones(n),-ones(n-1)).*((n+1))
    B = Tridiagonal(ones(n-1),4*ones(n),ones(n-1)).*(1/(6*(n+1)))
    F = lu(B)
    mode = 2
    bmat = :G
    ncv = 10
    nev = 4
    _reset_libarpack_dgetv0_iseed()
    _check_dsaupd_sequence!(A; B,invBop=F, mode, ncv, nev, bmat)

    _reset_libarpack_dgetv0_iseed()
    _check_dsaupd_sequence!(A; B,invBop=F, mode, ncv, nev, bmat, which=:SA)


    
  end

  @testset "Problem that needs early shift and deflation" begin
    using SparseArrays
    function loadmat(colptr, rowval, nzval)
      n = length(colptr)-1
      SparseMatrixCSC(n,n,colptr, rowval, nzval)
    end
    A = loadmat([1, 9, 16, 22, 27, 33, 38, 47, 51, 58, 62, 71, 78, 85, 89, 95], [2, 3, 4, 7, 8, 9, 11, 15, 1, 3, 5, 7, 11, 12, 14, 1, 2, 4, 5, 6, 9, 1, 3, 4, 6, 13, 2, 3, 6, 7, 12, 15, 3, 4, 5, 10, 13, 1, 2, 5, 8, 9, 10, 11, 12, 13, 1, 7, 9, 11, 1, 3, 7, 8, 13, 14, 15, 6, 7, 11, 15, 1, 2, 7, 8, 10, 12, 13, 14, 15, 2, 5, 7, 11, 12, 13, 14, 4, 6, 7, 9, 11, 12, 15, 2, 9, 11, 12, 1, 5, 9, 10, 11, 13], [1.0619322548149914, 0.8766879002200858, -0.9064189303728223, 0.21600317366983124, -1.283210593569618, 0.94558609230275, 1.9155428985937315, 0.9974159980590636, 1.0619322548149914, 0.07774792930492884, 0.2649798493656264, 0.6491820935668625, -0.6188565210918642, 0.33167352440104847, 0.06581708949429992, 0.8766879002200858, 0.07774792930492884, -0.831030685702723, -0.9593728263464558, -1.0783982816827704, 1.5537090724942313, -0.9064189303728223, -0.831030685702723, -2.73568483736268, -0.2601662488800133, -2.0138536457007397, 0.2649798493656264, -0.9593728263464558, 1.2063734934686803, -0.23345707234182503, -0.598980265258237, -3.235531538166721, -1.0783982816827704, -0.2601662488800133, 1.2063734934686803, -1.2968161204379096, -0.5037853135483386, 0.21600317366983124, 0.6491820935668625, -0.23345707234182503, -1.635181006451923, 1.2505504541996117, -1.8979584955867908, 0.8078340897940672, 0.3488529638883331, -0.8853907377562747, -1.283210593569618, -1.635181006451923, -1.0402874079575923, -0.023802716431538706, 0.94558609230275, 1.5537090724942313, 1.2505504541996117, -1.0402874079575923, -0.4150736664751938, -0.7608617919970803, -0.1985227463733881, -1.2968161204379096, -1.8979584955867908, 0.4120271156295672, -0.033303822949929375, 1.9155428985937315, -0.6188565210918642, 0.8078340897940672, -0.023802716431538706, 0.4120271156295672, -1.7060217210513384, 0.6923442199350044, 0.5061616454994609, 0.5200080839628756, 0.33167352440104847, -0.598980265258237, 0.3488529638883331, -1.7060217210513384, -0.4259291978517505, -0.7347260397652922, -0.45812838793146876, -2.0138536457007397, -0.5037853135483386, -0.8853907377562747, -0.4150736664751938, 0.6923442199350044, -0.7347260397652922, -0.28756151077851555, 0.06581708949429992, -0.7608617919970803, 0.5061616454994609, -0.45812838793146876, 0.9974159980590636, -3.235531538166721, -0.1985227463733881, -0.033303822949929375, 0.5200080839628756, -0.28756151077851555])
    _reset_libarpack_dgetv0_iseed()
    _check_dsaupd_sequence!(A; mode=1, ncv=6, nev=3, bmat=:I)
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