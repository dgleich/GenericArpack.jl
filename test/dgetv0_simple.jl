@testset "dgetv0" begin
  using LinearAlgebra
  @testset "Initial vector" begin
    # Arpack v0 
    # Get this by running helpers/get_arpack_jll_initial_random.jl
    v0 = [0.3957424639187579
    0.0008649603975001696
    -0.9227205789982591
    -0.9165671495278005
    0.1175963848841306
    -0.2996262520371218
    0.9038269570258635
    -0.25045104802183715
    0.33224741301423677
    -0.29023922021963955]
    ido = Ref{Int}(0)
    bmat = :I
    itry = 1
    initv = false
    n = 10
    j = 1
    v = zeros(n, j+1)
    ldv = n
    resid = 2*rand(n).-1
    rnorm = Ref{Float64}(0.0)
    ipntr = zeros(Int, 3)
    workd = zeros(2n)

    resid0 = copy(resid)

    state=ArpackInJulia.ArpackState{Float64}()
    ierr = ArpackInJulia.dgetv0!(
      ido, Val(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state)
    @test resid == v0
    @test state.getv0.iseed[] != tuple(1,3,5,7)
    @test ido[] == 99
    @test ierr==0
  end
  @testset "Standard Simple Eigenproblem" begin
    ido = Ref{Int}(0)
    bmat = :I
    itry = 1
    initv = true
    n = 10
    j = 1
    v = zeros(n, j+1)
    ldv = n
    resid = 2*rand(n).-1
    rnorm = Ref{Float64}(0.0)
    ipntr = zeros(Int, 3)
    workd = zeros(2n)
    M = Diagonal(1.0:10.0)

    resid0 = copy(resid)

    state=ArpackInJulia.ArpackState{Float64}()
    ierr = ArpackInJulia.dgetv0!(
      ido, Val(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state)
    @test ierr == 0
    @test state.getv0.first == false
    @test state.getv0.orth == false
    @test state.getv0.iseed[] == tuple(1,3,5,7)
    @test ido[] == 99
    @test resid == resid0
    @test rnorm[] == norm(resid0)

    #= Old code for testing the updated dgetv0.jl
    @test ido[] == -1
    @test state.getv0.first == false
    @test state.getv0.orth == false
    @test state.getv0.iseed[] == tuple(1,3,5,7)
    if ido[]==-1
      @test @view(workd[ipntr[1]:ipntr[1]+n-1]) == resid0
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    end
    ierr = ArpackInJulia.dgetv0!(
      ido, Val(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state)
    @test state.getv0.first == false
    @test state.getv0.orth == false
    @test state.getv0.iseed[] == tuple(1,3,5,7)
    @test ido[] == 99
    @test ierr==0
    @test resid == M*resid0
    =#
  end

  @testset "Standard Generalized Eigenproblem" begin
    ido = Ref{Int}(0)
    bmat = :G
    itry = 1
    initv = true
    n = 10
    j = 1
    v = zeros(n, j+1)
    ldv = n
    resid = 2*rand(n).-1
    rnorm = Ref{Float64}(0.0)
    ipntr = zeros(Int, 3)
    workd = zeros(2n)
    M = Diagonal(1.0:10.0)
    B = Diagonal(0.1:0.1:1)

    resid0 = copy(resid)

    state=ArpackInJulia.ArpackState{Float64}()
    ierr = ArpackInJulia.dgetv0!(
      ido, Val(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state)
    @test ierr == 0
    @test ido[] == -1
    @test state.getv0.first == false
    @test state.getv0.orth == false
    @test state.getv0.iseed[] == tuple(1,3,5,7)
    if ido[]==-1
      @test @view(workd[ipntr[1]:ipntr[1]+n-1]) == resid0
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    end
    ierr = ArpackInJulia.dgetv0!(
      ido, Val(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state)
    @test state.getv0.first == true
    @test state.getv0.orth == false
    @test state.getv0.iseed[] == tuple(1,3,5,7)
    @test ido[] == 2
    @test ierr==0
    @test resid == M*resid0

    if ido[]==2
      @test @view(workd[ipntr[1]:ipntr[1]+n-1]) == M*resid0
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    end


    ierr = ArpackInJulia.dgetv0!(
      ido, Val(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state)
    
    @test state.getv0.first == false
    @test state.getv0.orth == false
    @test ido[] == 99
    @test ierr==0
    @test rnorm[] ≈ sqrt(abs((M*resid0)'*(B*M*resid0)))
  end

  @testset "Restarted case where bmat is I" begin
    using LinearAlgebra
    using Random
    Random.seed!(0)
    ido = Ref{Int}(0)
    bmat = :I
    itry = 1
    initv = true
    n = 5
    j = 6
    v = Matrix(1.0I, n, j+1)
    ldv = n
    resid = 2*rand(n).-1
    rnorm = Ref{Float64}(0.0)
    ipntr = zeros(Int, 3)
    workd = zeros(2n)
    M = Diagonal(1.0:n)

    resid0 = copy(resid)

    state=ArpackInJulia.ArpackState{Float64}()
    ierr = ArpackInJulia.dgetv0!(
      ido, Val(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state)
    @test ierr == -1

    # try another case that should work
    initv = true
    v = zeros(n,3)
    v[:,1] .= 1/sqrt(n)
    j = 2
    fill!(workd, 0)
    fill!(ipntr, 0)
    ido[] = 0
    itry = 1
    resid = 2*rand(n).-1

    state=ArpackInJulia.ArpackState{Float64}()
    debug=ArpackInJulia.ArpackDebug(;logfile=IOBuffer())
    stats=ArpackInJulia.ArpackStats()
    debug.mgetv0 = 4
    ierr = ArpackInJulia.dgetv0!(
      ido, Val(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state, debug, stats)
    @test ierr == 0
    @test isapprox(sum(resid), 0.0; atol=n*eps(1.0))
    @test rnorm[] > 0
  end

  @testset "Restarted case where bmat is G" begin
    using LinearAlgebra
    using Random
    ido = Ref{Int}(0)
    bmat = :G
    itry = 1
    initv = true
    n = 8
    j = 3
    v = Matrix(1.0I,n,j+1)

    Random.seed!(0)
    ldv = n
    resid = 2*rand(n).-1
    rnorm = Ref{Float64}(0.0)
    ipntr = zeros(Int, 3)
    workd = zeros(2n)
    M = Diagonal(1.0:n)
    B = Diagonal(zeros(n))
    B[1,1] = 1
    B[2,2] = 1

    resid0 = copy(resid)
    state=ArpackInJulia.ArpackState{Float64}()
    ierr = -1
    niter = 0
    niterB = 0 
    stats=ArpackInJulia.ArpackStats()
    while ido[] != 99
      ierr = ArpackInJulia.dgetv0!(
        ido, Val(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
        state,stats)
      niter += 1
      if ido[] == -1
        mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
      elseif ido[] == 2
        niterB += 1
        mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
      end
    end
    @test ierr==-1
    @test niterB == 7 # initial + 6 (in dgetv0 iter = 0...5) more... 
    @test niter == niterB + 2 # initial + final when ido == 99
    @test stats.nbx == niterB # initial 
    @test stats.nopx == 1 # initial 
    @test stats.tgetv0 > 0 
    @test stats.tmvbx > 0 
    @test stats.tmvopx > 0 
  end
end