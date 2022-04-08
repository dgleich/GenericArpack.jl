using LinearAlgebra
@testset "dsaupd with idonow for simple" begin 
  op = ArpackInJulia.ArpackSimpleOp(Diagonal(1.0:30))
  prob1 = _allocate_symproblem(op, 12)
  prob2 = deepcopy(prob1)

  niter1, ierr1 = _eigrun!(prob1, 6)
  niter2, ierr2 = _eigrun!(prob2, 6; idonow=true)

  @test niter2 == 1
  @test niter1 > niter2 # needs at least one matvec...
  @test ierr1 == ierr2
  _compare_probs(prob1, prob2) # run tests on each field
  @test prob2.op == prob1.op
end

@testset "dsaupd with idonow for generalized" begin 
  n = 100
  ncv = 10
  nev = 4
  A = Tridiagonal(-ones(n-1),2*ones(n),-ones(n-1)).*((n+1))
  B = Tridiagonal(ones(n-1),4*ones(n),ones(n-1)).*(1/(6*(n+1)))
  op = ArpackSymmetricGeneralizedOp(A, lu!(copy(B)), B)
  prob1 = _allocate_symproblem(op, 10)
  prob2 = _allocate_symproblem(op, 10)

  stats1 = ArpackInJulia.ArpackStats()
  stats2 = ArpackInJulia.ArpackStats()
  debug = ArpackInJulia.ArpackDebug()
  #debug.maupd = 1
  #debug.maup2 = 1
  #ArpackInJulia.set_debug_high!(debug)
  
  niter1, ierr1 = _eigrun!(prob1, nev; mode=2, bmat=Val(:G), stats=stats1)
  niter2, ierr2 = _eigrun!(prob2, nev; mode=2, bmat=Val(:G), idonow=true, stats=stats2)

  #@show ierr1, ierr2
  #@show prob1.workl[2*ncv+1:2*ncv+nev]
  @test ierr1 == ierr2
  @test niter2 == 1
  @test niter1 > niter2 # needs at least one matvec...
  _compare_probs(prob1, prob2) # run tests on each field

  @test sort(eigvals(B\A), by=abs, rev=true)[1:nev] ≈ sort(prob1.workl[2*ncv+1:2*ncv+nev], by=abs, rev=true)
  @test sort(eigvals(B\A), by=abs, rev=true)[1:nev] ≈ sort(prob2.workl[2*ncv+1:2*ncv+nev], by=abs, rev=true)
end  