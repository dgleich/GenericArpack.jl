@testset "dsaitr compare against arpackjll" begin 
  function _check_dsaitr_sequence!(M; 
    B=1.0LinearAlgebra.I,
    idostart::Int,
    bmat::Symbol,
    n::Int,
    k::Int, 
    np::Int,
    mode::Int,
    resid::StridedVecOrMat{T},
    rnorm::Ref{T},
    V::StridedMatrix{T},
    ldv::Int,
    H::StridedMatrix{T},
    ldh::Int,
  ) where T
    resid0 = copy(resid)

    ido = Ref{Int}(idostart)
    ipntr = zeros(Int, 3)
    workd = zeros(3n)

    arido = Ref{Int}(idostart)
    arv = copy(V)
    arh = copy(H)
    arresid = copy(resid0)
    arrnorm = Ref{T}(rnorm[])
    aripntr = copy(ipntr)
    arworkd = copy(workd)

    histdata = Vector{
        NamedTuple{(:info,:ido,:rnorm), Tuple{Int64,Int64,T}}
    }()

    state = ArpackInJulia.ArpackState{Float64}()
    while ido[] != 99
      info = ArpackInJulia.dsaitr!(
        ido, Val(bmat), n, k, np, mode, resid, rnorm, V, ldv, H, ldh, ipntr, workd,
        state)

        
      arinfo = arpack_dsaitr!(
        arido, bmat, n, k, np, mode, arresid, arrnorm, arv, ldv, arh, ldh, 
        aripntr, arworkd)

      @debug arido[], arinfo, ido[], info

      @test arinfo == info
      @test arido[] == ido[]
      @test arv == V
      @test arh == H
      @test arresid == resid
      @test floatsbetween(arrnorm[], rnorm[]) <= 1
      @test aripntr == ipntr
      @test arworkd == workd  

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

      # everything should be same...
      copyto!(arworkd, workd)
      push!(histdata, (;info,ido=ido[],rnorm=rnorm[]))
    end
    return histdata
  end


  @testset "dsaitr compare arpackjll simple call" begin
    using LinearAlgebra
    using Random 
    Random.seed!(0)
    bmat = :I
    n = 10
    k = 0 # number of current columns in V
    np = 3
    #V = randn(n,k+np) # total memory for v
    #V[:,1:k] = qr(V).Q[:,1:k] # random orthogonal init
    mode = 1
    resid = ones(n)/sqrt(n)
    rnorm = Ref{Float64}(norm(resid))

    V = zeros(n,k+np)
    ldv = n 
    H = zeros(n,2) # full h
    ldh = n 

    M = Diagonal(1.0:n)
    
    @testset "run 1" begin 
      seqdata = _check_dsaitr_sequence!(M;  idostart=0, 
        bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
      )
    end 

    V = zeros(n,k+np)
    resid = ones(n)/sqrt(n)
    rnorm = Ref{Float64}(norm(resid))

    
    @testset "run 2" begin 
      # make sure there aren't any weird duplicate scenarios. 
      seqdata = _check_dsaitr_sequence!(M;  idostart=0, 
        bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
      )
    end 

    # now extend this...
    k += np 
    V = [V zeros(n,np)]

    @testset "extend Arnoldi 3->6" begin
      seqdata = _check_dsaitr_sequence!(M;  idostart=0, 
       bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
      )
    end

    # now extend this...
    V = [V zeros(n,np)]
    k += np 

    @testset "extend Arnoldi more 6->9" begin
      seqdata = _check_dsaitr_sequence!(M;  idostart=0, 
       bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
      )
    end

    # now extend this... beyond n... 
    V = [V zeros(n,np)]
    k += np 

    # need to adjust H at this point
    H = [H; zeros(5,2)]
    ldh = size(H,1)

    _reset_libarpack_dgetv0_iseed()
    @testset "extend Arnoldi beyond 9" begin

      seqdata = _check_dsaitr_sequence!(M;  idostart=0, 
       bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
      )
      @test seqdata[end].info != 0 
      @test seqdata[end].info == n 
      @test seqdata[end].rnorm == 0
    end
  end 

  @testset "dsaitr compare arpackjll generalized call" begin
    using LinearAlgebra
    using Random 
    Random.seed!(0)
    bmat = :G
    n = 10
    k = 0 # number of current columns in V
    np = 3
    #V = randn(n,k+np) # total memory for v
    #V[:,1:k] = qr(V).Q[:,1:k] # random orthogonal init
    mode = 1
    resid = ones(n)/sqrt(n)
    rnorm = Ref{Float64}(norm(resid))

    V = zeros(n,k+np)
    ldv = n 
    H = zeros(n,2) # full h
    ldh = n 

    M = Diagonal(1.0:n)
    B = Diagonal(range(0.1, 1.0, length=n))
    
    @testset "run 1" begin 
      seqdata = _check_dsaitr_sequence!(M; B, idostart=0, 
        bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
      )
    end 

    V = zeros(n,k+np)
    H = zeros(n,2) # full h
    resid = ones(n)/sqrt(n)
    rnorm = Ref{Float64}(norm(resid))

    
    @testset "run 2" begin 
      # make sure there aren't any weird duplicate scenarios. 
      seqdata = _check_dsaitr_sequence!(M; B, idostart=0, 
        bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
      )
    end 
  end 

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