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

  function parse_workl_dseupd(workl, ncv, nconv, ipntr)
    ih = ipntr[5]
    ritz = ipntr[6]
    bounds = ipntr[7]
    ldh = ncv
    ldq = ncv
    ihd = bounds + ldh # diagonal of H
    ihb = ihd + ldh    # subdiagonal of H 
    iq = ihb + ldh
    iw = iq + ldh*ncv 
    next = iw + 2*ncv

    # julia helpers in new motif...
    ncv0 = 0:ncv-1
    nconv0 = 0:nconv-1

    H = reshape(@view(workl[1:2*ncv]), ncv, 2)
    ritzvec = @view(workl[ritz .+ ncv0])
    boundsvec = @view(workl[bounds .+ ncv0])
    hdvec = @view(workl[ihd .+ ncv0])
    hbvec = @view(workl[ihb .+ ncv0])
    #workvec = @view(workl[iw .+ 0:(2*ncv-1)])
    workvec = @view(workl[iw:iw+2*ncv-1])
    Q = reshape(@view(workl[iq:iq+ldq*ncv-1]), ldq, ncv)
    Q1 = reshape(@view(workl[iq:iq+ldq*nconv-1]), ldq, nconv)

    @show length(workvec)
    @show iw:iw+2*ncv-1
    @show iw .+ 0:(2*ncv-1)
    @show length(workl)
    
    return H, ritzvec, boundsvec, hdvec, hbvec, workvec, Q, Q1
  end 

  parse_workl_dseupd(prob) = parse_workl_dseupd(prob.workl, size(prob.V,2), prob.iparam[5], prob.ipntr)

  @testset "dseupd with alt spectrum" begin
    T= Float64 
    bmat = :I
    BMAT = Val(bmat)
    which = :LM
    n = 25 
    tol = eps(T)/2
    op = ArpackInJulia.ArpackSimpleOp(Diagonal(range(-10, 9, length=n)))
    ncv = 12
    prob = _allocate_symproblem(op, ncv)
    nev = 6
    niter, ierr = _eigrun!(prob, nev; which, bmat=BMAT)

    rvec = false
    sigma = 0.0

    Z = zeros(n, nev)
    d = zeros(nev) 
    select = zeros(Int, ncv)

    arprob = deepcopy(prob)
    ard = copy(d)
    arZ = copy(Z)
    arselect = copy(select)

    #test_diffs("before workl: ", prob.workl, arprob.workl)

    #H, ritzvec, boundsvec, hdvec, hbvec, workvec, Q, Q1 = parse_workl_dseupd(prob)
    #@show boundsvec

    debug = ArpackInJulia.ArpackDebug()
    #ArpackInJulia.set_debug_high!(debug)
    ierr = ArpackInJulia.simple_dseupd!(rvec, select, d, Z, sigma,  BMAT, n, which, nev, 
      tol, 
      prob.resid, ncv, prob.V, prob.iparam, prob.ipntr, 
      prob.workd, prob.workl; debug
    )

    #arpack_set_debug_high()
    arierr = arpack_dseupd!(rvec, arselect, ard, arZ, 
      stride(arZ,2), sigma, bmat, n, which, nev, tol, 
      arprob.resid, ncv, arprob.V, stride(arprob.V,2), arprob.iparam, arprob.ipntr, 
      arprob.workd, arprob.workl, length(arprob.workl))

    
    # H, ritzvec, boundsvec, hdvec, hbvec, workvec, Q, Q1 = parse_workl_dseupd(prob)
    # @show H
    # @show ritzvec 
    # @show boundsvec
    # @show hdvec
    # @show hbvec
    # @show workvec


    # test_diffs("workl: ", prob.workl, arprob.workl)

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



  @testset "Generalized Eigenproblem from dsdrv3" begin 
    n = 100
    A = Tridiagonal(-ones(n-1),2*ones(n),-ones(n-1)).*((n+1))
    B = Tridiagonal(ones(n-1),4*ones(n),ones(n-1)).*(1/(6*(n+1)))
    F = lu(B)
    op = ArpackInJulia.ArpackSymmetricGeneralizedOp(A,F,B)
    bmat = :G
    BMAT = Val(bmat)
    mode = ArpackInJulia.arpack_mode(op)
    which = :SA
    sigma = 0.0
    tol = eps(Float64)/2

    ncv = 10
    nev = 4
    prob = _allocate_symproblem(op, ncv)
    niter, ierr = _eigrun!(prob, nev; which, bmat=BMAT, maxiter=400)

    Z = zeros(n, nev)
    d = zeros(nev) 
    select = zeros(Int, ncv)

    arprob = deepcopy(prob)
    backupprob = deepcopy(prob)

    ard = copy(d)
    arZ = copy(Z)
    arselect = copy(select)

    rvec = false

    debug = ArpackInJulia.ArpackDebug()
    #ArpackInJulia.set_debug_high!(debug)
    ierr = ArpackInJulia.simple_dseupd!(rvec, select, d, Z, sigma,  BMAT, n, which, nev, 
      tol, 
      prob.resid, ncv, prob.V, prob.iparam, prob.ipntr, 
      prob.workd, prob.workl; debug
    )

    #arpack_set_debug_high()
    arierr = arpack_dseupd!(rvec, arselect, ard, arZ, 
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

    rvec = true
    prob = deepcopy(backupprob)
    arprob = deepcopy(backupprob)
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
    arierr = arpack_dseupd!(rvec, arselect, ard, arZ, 
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



  @testset "Generalized Eigenproblem from dsdrv3" begin 
    n = 100
    A = Tridiagonal(-ones(n-1),2*ones(n),-ones(n-1)).*((n+1))
    B = Tridiagonal(ones(n-1),4*ones(n),ones(n-1)).*(1/(6*(n+1)))
    F = lu(B)
    op = ArpackInJulia.ArpackSymmetricGeneralizedOp(A,F,B)
    bmat = :G
    BMAT = Val(bmat)
    mode = ArpackInJulia.arpack_mode(op)
    which = :LM
    sigma = 0.0
    tol = eps(Float64)/2

    ncv = 10
    nev = 4
    prob = _allocate_symproblem(op, ncv)
    niter, ierr = _eigrun!(prob, nev; which, bmat=BMAT, maxiter=400)

    Z = zeros(n, nev)
    d = zeros(nev) 
    select = zeros(Int, ncv)

    arprob = deepcopy(prob)
    backupprob = deepcopy(prob)

    ard = copy(d)
    arZ = copy(Z)
    arselect = copy(select)

    rvec = false

    debug = ArpackInJulia.ArpackDebug()
    #ArpackInJulia.set_debug_high!(debug)
    ierr = ArpackInJulia.simple_dseupd!(rvec, select, d, Z, sigma,  BMAT, n, which, nev, 
      tol, 
      prob.resid, ncv, prob.V, prob.iparam, prob.ipntr, 
      prob.workd, prob.workl; debug
    )

    #arpack_set_debug_high()
    arierr = arpack_dseupd!(rvec, arselect, ard, arZ, 
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

    rvec = true
    prob = deepcopy(backupprob)
    arprob = deepcopy(backupprob)
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
    arierr = arpack_dseupd!(rvec, arselect, ard, arZ, 
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
