@testset "dsaupd errors" begin 

  #=
  = -1: N must be positive.
         = -3: NCV must be greater than NEV and less than or equal to N.
         = -4: The maximum number of Arnoldi update iterations allowed
               must be greater than zero.
         = -8: Error return from trid. eigenvalue calculation;
               Informatinal error from LAPACK routine dsteqr .
         = -12: IPARAM(1) must be equal to 0 or 1.
         =#
  
  @testset "initial vector is zero" begin 
    n = 10 
    op = ArpackSimpleOp(Diagonal(1.0:n))
    prob = _allocate_symproblem(op, n-1)
    niter, ierr = _eigrun!(prob, 6; initv=zeros(n))
    @test ierr==-9 # initial vector is 0
  end 

  @testset "nev must be positive" begin 
    n = 5
    op = ArpackSimpleOp(Diagonal(1.0:n))
    prob = _allocate_symproblem(op, 4)
    niter, ierr = _eigrun!(prob, 0)
    @test ierr==-2 # nev needs to be positive 
  end

  @testset "n must be positive" begin 
    faken = 5
    ncv = 4
    op = ArpackSimpleOp(Diagonal(1.0:faken))
    prob = _allocate_symproblem(op, ncv)
    #prob.op = ArpackSimpleOp(Diagonal(1.0:0))
    # allocate much too large.
    prob = ArpackSymProblem(ArpackSimpleOp(Diagonal(1.0:0)), 
      zeros(0,ncv), prob.resid, prob.workd, prob.workl, prob.iparam, prob.ipntr, prob.ido)
    niter, ierr = _eigrun!(prob, 1)
    @test ierr==-1 # n must be positive 
  end 

  @testset "workl isn't long enough" begin 
    n = 5 
    ncv = 3
    op = ArpackSimpleOp(Diagonal(1.0:n))
    prob = _allocate_symproblem(op, ncv)
    prob = ArpackSymProblem(prob.op, 
      prob.V, prob.resid, prob.workd, zeros(ncv*ncv), prob.iparam, prob.ipntr, prob.ido)
    
    niter, ierr = _eigrun!(prob, 2)
    @test ierr==-7 # initial vector is 0
  end 

  @testset "simple parameter problems" begin 
    n = 5 
    ncv = 3
    op = ArpackSimpleOp(Diagonal(1.0:n))
    prob = _allocate_symproblem(op, ncv)
    niter, ierr = _eigrun!(prob, 1; which=:BE)
    @test ierr==-13 # which=:BE and also odd
    niter, ierr = _eigrun!(prob, 2; which=:SS)
    @test ierr==-5 # WHICH must be one of 'LM", 'SM", 'LA", 'SA' or 'BE'.
    niter, ierr = _eigrun!(prob, 2; which=:SS, bmat=Val(:B))
    @test ierr== -6 #: BMAT must be one of 'I' or 'G'.
    niter, ierr = _eigrun!(prob, 2; which=:SS, bmat=Val(:G), mode=1)
    @test ierr== -11 # -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
    niter, ierr = _eigrun!(prob, 2; which=:SS, bmat=Val(:G), mode=6)
    @test ierr== -10 # -10: IPARAM(7) must be 1,2,3,4,5.  


    niter, ierr = _eigrun!(prob, 2; maxiter=-1)
    @test ierr == -4 # maxiter needs to be positive


    
  end 

  @testset "ncv must be greater than nev and less than or equal to n" begin 

    # annoying to test because symproblem/eigrun enforces these properties :) 
    n = 5 
    ncv = 4
    op = ArpackSimpleOp(Diagonal(1.0:n))
    prob = _allocate_symproblem(op, ncv)
    niter, ierr = _eigrun!(prob, 3; ncv=3)
    @test ierr == -3 # NCV must be greater than NEV and less than or equal to N.

    # another case for ncv 
    prob = ArpackSymProblem(op,
      zeros(n,n+1), prob.resid, prob.workd, zeros((10+(n+1))*(n+1)), prob.iparam, prob.ipntr, prob.ido)
    niter, ierr = _eigrun!(prob, 2)
    @test ierr == -3 # NCV must be greater than NEV and less than or equal to N.
  end 
end 

