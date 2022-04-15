function setup_saitr_problem(M; 
  np, 
  V = nothing, 
  H = nothing, 
  idonow = false,
  B = nothing,
  v0 = nothing, 
)
  n = size(M)
  bmat = :I
  if B !== nothing
    bmat = :G
  end
  k = 0 
  if V !== nothing
    @assert H !== nothing # need an H to go with V....
    k = size(V,2)
    V = [V zeros(n,np)]
  else
    H = zeros(n,2)
    V = zeros(n,k+np)
  end 
  return (; n, V, H, bmat, idonow)
end 

function _run_saitr_sequence!(M; 
  B=1.0LinearAlgebra.I,
  idostart::Int,
  bmat::Symbol,
  n::Int,
  k::Int, 
  np::Int,
  mode::Int,
  resid::AbstractVecOrMat{T},
  rnorm::Ref{T},
  V::AbstractMatrix{T},
  ldv::Int,
  H::AbstractMatrix{T},
  ldh::Int,
  stats = nothing,
  debug = nothing, 
  state = nothing, 
  idonow = nothing, 
) where T
  resid0 = copy(resid)

  @assert(size(M,1) == n)

  ido = Ref{Int}(idostart)
  ipntr = zeros(Int, 3)
  workd = zeros(3n)

  histdata = Vector{
      NamedTuple{(:info,:ido,:rnorm), Tuple{Int64,Int64,T}}
  }()

  if state === nothing
    state = GenericArpack.ArpackState{Float64}()
  end 
  while ido[] != 99
    info = GenericArpack.dsaitr!(
      ido, Val(bmat), n, k, np, mode, resid, rnorm, V, ldv, H, ldh, ipntr, workd, state; 
      stats, debug, idonow)

    if ido[] == -1 || ido[] == 1
      if mode == 2
        # crazy interface, see remark 5 
        mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
        copyto!(@view(workd[ipntr[1]:ipntr[1]+n-1]),@view(workd[ipntr[2]:ipntr[2]+n-1]))
        ldiv!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
      else
        mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
      end 
    elseif ido[] == 2
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    elseif ido[] == 99 
      # we are done... will exit...
    else
      @error("Wrong ido, $(ido[])")
    end

    push!(histdata, (;info,ido=ido[],rnorm=rnorm[]))
  end
  return histdata
end

@testset "dsaitr Simple Checks" begin
  @testset "Simple Standard Eigenproblem" begin 
    using LinearAlgebra
    using Random 
    Random.seed!(0)
    ido = Ref{Int}(0)
    bmat = :I
    n = 10
    k = 0 # number of current columns in V
    np = 3
    mode = 1
    resid = randn(n)
    rnorm = Ref{Float64}(norm(resid))

    V = zeros(n,k+np)
    ldv = n 
    H = zeros(n,2) # full h
    ldh = n 

    M = Diagonal(1.0:n)

    debug = GenericArpack.ArpackDebug(logfile=IOBuffer())
    GenericArpack.set_debug_high!(debug)
    stats = ArpackStats()
    _run_saitr_sequence!(M; idostart=0,
      n, k, np, mode, resid, rnorm, V, H, ldv, ldh, stats, debug, bmat
    )
    @test stats.nrstrt == 0 
    @test stats.tgetv0 == 0 
    @test stats.nopx > 0 
    @test stats.taitr > 0 
    @test stats.tmvopx > 0 
    @test V'*V ≈ Matrix(1.0I,np,np)
    @test norm(V'*resid) ≈ 0 atol=eps(1.0)

    k = 3 
    np = 6
    V = [V zeros(n, np)]
    # now test extending the factorization...
    stats = ArpackStats()
    _run_saitr_sequence!(M; idostart=0,
      n, k, np, mode, resid, rnorm, V, H, ldv, ldh, stats, debug, bmat
    )

    @test stats.nrstrt == 0 
    @test stats.tgetv0 == 0 
    @test stats.nopx > 0 
    @test stats.taitr > 0 
    @test stats.tmvopx > 0 
    @test V'*V ≈ Matrix(1.0I,k+np,k+np)
    @test norm(V'*resid) ≈ 0 atol=eps(1.0)
  end 

  @testset "Check Identity with Restart" begin
    using LinearAlgebra
    using Random 
    Random.seed!(0)
    ido = Ref{Int}(0)
    bmat = :I
    n = 10
    k = 0 # number of current columns in V
    np = 1
    mode = 1
    resid = collect(1.0:n)
    rnorm = Ref{Float64}(norm(resid))
    V = zeros(n,k+np+1)
    ldv = n 
    H = zeros(n,2) # full h
    ldh = n 
    M = Diagonal(1.0I, 10) # 10x10 identity 


    stats = ArpackStats()
    @test_nowarn rval = _run_saitr_sequence!(M; idostart=0,
      n, k, np, mode, resid, rnorm, V, H, ldv, ldh, bmat, stats
    )
    @test stats.nrstrt == 0 
    # resid should be zero... 
    @test norm(resid) ≈ 0 atol=n*eps(1.0)

    k = 0 # number of current columns in V
    np = 2
    resid = collect(1.0:n)
    rnorm = Ref{Float64}(norm(resid))
    @test_nowarn rval = _run_saitr_sequence!(M; idostart=0,
      n, k, np, mode, resid, rnorm, V, H, ldv, ldh, bmat, stats
    )
    @test stats.nrstrt > 0 
    @test stats.tgetv0 > 0 
    # because it's the identity, resid should still be zero...
    @test norm(resid) ≈ 0 atol=n*eps(1.0)
  end 

  @testset "Check zero initial resid" begin 
    using LinearAlgebra
    using Random 
    Random.seed!(0)
    ido = Ref{Int}(0)
    bmat = :I
    n = 10
    k = 0 # number of current columns in V
    np = 3
    mode = 1
    resid = randn(n)
    rnorm = Ref{Float64}(0.0)

    V = zeros(n,k+np)
    ldv = n 
    H = zeros(n,2) # full h
    ldh = n 

    M = Diagonal(1.0:n)

    debug = GenericArpack.ArpackDebug(logfile=IOBuffer())
    stats = ArpackStats()
    _run_saitr_sequence!(M; idostart=0,
      n, k, np, mode, resid, rnorm, V, H, ldv, ldh, stats, bmat
    )
    @test stats.nrstrt > 0 
    @test stats.tgetv0 > 0 
    @test stats.nopx > 0 
    @test stats.taitr > 0 
    @test stats.tmvopx > 0 
    @test V'*V ≈ Matrix(1.0I,np,np)
    @test norm(V'*resid) ≈ 0 atol=2*eps(1.0)
  end

  @testset "Simple Generalized Eigenproblem with Identity" begin 
    using LinearAlgebra
    using Random 
    Random.seed!(0)
    ido = Ref{Int}(0)
    bmat = :G
    n = 10
    k = 0 # number of current columns in V
    np = 3
    mode = 1
    resid = randn(n)
    rnorm = Ref{Float64}(norm(resid))

    V = zeros(n,k+np)
    ldv = n 
    H = zeros(n,2) # full h
    ldh = n 

    M = Diagonal(1.0:n)
    B = 1.0*I

    debug = GenericArpack.ArpackDebug(logfile=IOBuffer())
    stats = ArpackStats()
    _run_saitr_sequence!(M; B, idostart=0,
      n, k, np, mode, resid, rnorm, V, H, ldv, ldh, stats, debug, bmat
    )
    @test stats.nrstrt == 0 
    @test stats.tgetv0 == 0 
    @test stats.nopx > 0 
    @test stats.nbx > 0
    @test stats.taitr > 0 
    @test stats.tmvopx > 0 
    @test stats.tmvbx > 0 
    @test V'*V ≈ Matrix(1.0I,np,np)
    @test norm(V'*resid) ≈ 0 atol=eps(1.0)
  end 

  @testset "Simple Generalized Eigenproblem with Diagonal" begin 
    using LinearAlgebra
    using Random 
    Random.seed!(0)
    ido = Ref{Int}(0)
    bmat = :G
    n = 10
    k = 0 # number of current columns in V
    np = 3
    mode = 2
    resid = randn(n)

    A = Diagonal(1.0:n)
    B = Diagonal(collect(range(0.1, n, length=n)))
    M = inv(B)*A # okay since they are diagonal....

    V = zeros(n,k+np)
    ldv = n 
    H = zeros(n,2) # full h
    ldh = n 

    rnorm = Ref{Float64}(sqrt(abs(resid'*B*resid)))

    #debug = GenericArpack.ArpackDebug(logfile=IOBuffer())
    debug = GenericArpack.ArpackDebug()
    #GenericArpack.set_debug_high!(debug)
    stats = ArpackStats()
    _run_saitr_sequence!(A; B, idostart=0,
      n, k, np, mode, resid, rnorm, V, H, ldv, ldh, stats, debug, bmat
    )
    @test rnorm[] ≈ sqrt(abs(resid'*B*resid))
    @test stats.nrstrt == 0 
    @test stats.tgetv0 == 0 
    @test stats.nopx > 0 
    @test stats.nbx > 0
    @test stats.taitr > 0 
    @test stats.tmvopx > 0 
    @test stats.tmvbx > 0 
    @test V'*B*V ≈ Matrix(1.0I,np,np)
    @test norm(V'*B*resid) ≈ 0 atol=eps(1.0)
  end 

  @testset "Check initial workd" begin 
    @testset "dsaitr simple call" begin
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

      state = GenericArpack.ArpackState{Float64}()

      resid0 = copy(resid)

      info = GenericArpack.dsaitr!(
        ido, Val(bmat), n, k, np, mode, resid, rnorm, 
        V, ldv, H, ldh, ipntr, workd,
        state 
      )
      @test workd[ipntr[1]:ipntr[1]+n-1] ≈ resid0/norm(resid0)

    end
  end 

  @testset "Simple Standard Eigenproblem with idonow" begin 
    using LinearAlgebra
    using Random 
    Random.seed!(0)
    ido = Ref{Int}(0)
    bmat = :I
    n = 10
    k = 0 # number of current columns in V
    np = 3
    mode = 1
    resid = randn(n)
    rnorm = Ref{Float64}(norm(resid))

    V = zeros(n,k+np)
    ldv = n 
    H = zeros(n,2) # full h
    ldh = n 

    M = Diagonal(1.0:n)
    idonow = GenericArpack.ArpackSimpleOp(M)
    stats = ArpackStats()
    hist = _run_saitr_sequence!(M; idostart=0,
      n, k, np, mode, resid, rnorm, V, H, ldv, ldh, stats, bmat, idonow
    )
    @test length(hist) == 1 
    @test stats.nrstrt == 0 
    @test stats.tgetv0 == 0 
    @test stats.nopx > 0 
    @test stats.taitr > 0 
    @test stats.tmvopx > 0 
    @test V'*V ≈ Matrix(1.0I,np,np)
    @test norm(V'*resid) ≈ 0 atol=eps(1.0)

    
  end 
end 