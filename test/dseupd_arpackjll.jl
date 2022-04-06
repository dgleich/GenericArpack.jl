@testset "dseupd compare to arpackjll" begin
  using LinearAlgebra
  @testset "dseupd simple" begin
    T= Float64 
    bmat = :I
    BMAT = Val(bmat)
    which = :LM
    n = 30 
    tol = eps(T)/2
    op = ArpackInJulia.ArpackSimpleOp(Diagonal(1.0:n))
    ncv = 12
    prob = _allocate_symproblem(op, ncv)
    nev = 6
    niter, ierr = _eigrun!(prob, nev; which, bmat=BMAT)

    rvec = true
    sigma = 0.0

    Z = zeros(n, nev)
    d = zeros(nev) 
    select = zeros(Int, ncv)

    arprob = deepcopy(prob)
    ard = copy(d)
    arZ = copy(Z)
    arselect = copy(select)

    debug = ArpackInJulia.ArpackDebug()
    #ArpackInJulia.set_debug_high!(debug)
    ierr = ArpackInJulia.simple_dseupd!(rvec, select, d, Z, sigma,  BMAT, n, which, nev, 
      tol, 
      prob.resid, ncv, prob.V, prob.iparam, prob.ipntr, 
      prob.workd, prob.workl; debug
    )

    #arpack_set_debug_high()
    arierr = arpack_dseupd!(true, arselect, ard, arZ, 
      stride(arZ,2), sigma, bmat, n, which, nev, tol, 
      arprob.resid, ncv, arprob.V, stride(arprob.V,2), arprob.iparam, arprob.ipntr, 
      arprob.workd, arprob.workl, length(arprob.workl))

    @test ierr == arierr 
    @test ierr == 0 
    @test d == ard 
    @test Z == arZ
    @test prob.workd == arprob.workd 
    @test prob.workl == arprob.workl 
    @test prob.ipntr == arprob.ipntr # just a few workspace updates
    @test prob.iparam == arprob.iparam # shouldn't change ... 
    @test prob.V == arprob.V # will get hit by orthogonalization... 
  end
  
end
