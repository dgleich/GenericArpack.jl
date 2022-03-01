
# a helper function to check a sequence of calls and 
# make sure all intermediate data is the same 
function _check_dgetv0_sequence!(M; 
  B=1.0LinearAlgebra.I,
  idostart::Int,
  bmat::Symbol,
  itry::Int, 
  initv::Bool,
  n::Int,
  j::Int,
  v::AbstractVecOrMat{T},
  ldv::Int,
  resid::AbstractVecOrMat{T},
  rnorm::Ref{T}
) where T

  resid0 = copy(resid)

  ido = Ref{Int}(idostart)
  ipntr = zeros(Int, 3)
  workd = zeros(2n)

  # allocate a separate copy of Arpack
  arido = Ref{Int}(0)
  arv = copy(v)
  arresid = copy(resid0)
  arrnorm = Ref{Float64}(rnorm[])
  aripntr = copy(ipntr)
  arworkd = copy(workd)

  histdata = Vector{
      NamedTuple{(:ierr,:ido,:rnorm), Tuple{Int64,Int64,T}}
  }()

  # run the sequence
  state = ArpackInJulia.ArpackState{Float64}()
  while ido[] != 99
    ierr = ArpackInJulia.dgetv0!(
      ido, Val{bmat}, itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state)
    
    arierr = arpack_dgetv0!(
      arido, bmat, itry, initv, n, j, arv, ldv, arresid, arrnorm, aripntr, arworkd)

    @debug arido[], arierr, ido[], ierr

    @test arierr == ierr
    @test arido[] == ido[]
    @test arv == v
    @test arresid == resid
    @test arrnorm[] ≈ rnorm[]
    @test aripntr == ipntr
    @test arworkd == workd  
    if ido[] == -1
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    elseif ido[] == 2
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    end
    # everything should be same...
    copyto!(arworkd, workd)
    push!(histdata, (;ierr,ido=ido[],rnorm=rnorm[]))
  end
  return histdata
end

@testset "Simple case to check against ArpackJLL" begin
  using LinearAlgebra
  using Random
  Random.seed!(0)
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

  resid0 = copy(resid)

  state=ArpackInJulia.ArpackState{Float64}()
  ierr = ArpackInJulia.dgetv0!(
    ido, Val{bmat}, itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
    state)

  arido = Ref{Int}(0)
  #bmat = :I # these are static... 
  #itry = 1
  #initv = true
  #n = 10
  #j = 1
  arv = zeros(n, j+1)
  #ldv = n
  arresid = copy(resid0)
  arrnorm = Ref{Float64}(0.0)
  aripntr = zeros(Int, 3)
  arworkd = zeros(2n)
  
  arierr = arpack_dgetv0!(
    arido, bmat, itry, initv, n, j, arv, ldv, arresid, arrnorm, aripntr, arworkd)
  
  @test arierr == ierr
  @test arido[] == ido[]
  @test arv == v
  @test arresid == resid
  @test arrnorm[] == rnorm[]
  @test aripntr == ipntr
  @test arworkd == workd
end

@testset "Simple Generalized Eigenproblem" begin 
  using LinearAlgebra
  using Random
  Random.seed!(0)
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

  seqdata = _check_dgetv0_sequence!(1.0*I;  idostart=0, 
    rnorm, resid, n, j, v, bmat, ldv, itry, initv)

  @test seqdata[end].ierr == 0
  oneresid = sqrt(abs((M*resid0)'*(B*M*resid0)))
  rval = begin 
    for v in seqdata # one of these should be true...
      if v.rnorm ≈ oneresid 
        return true 
      end
    end
    return false
  end
  @test rval == true 
end 

@testset "Restarted case where bmat is G and singular" begin
  using LinearAlgebra
  using Random
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
  M = Diagonal(1.0:n)
  B = Diagonal(zeros(n))
  B[1,1] = 1
  B[2,2] = 1

  seqdata = _check_dgetv0_sequence!(M; B, idostart=0, 
    rnorm, resid, n, j, v, bmat, ldv, itry, initv )

  @test seqdata[end].ierr==-1
end

##
@testset "Simple Restarted Eigenproblem that succeeds" begin
  using LinearAlgebra
  Random.seed!(0)
  ido = Ref{Int}(0)
  bmat = :I
  itry = 1
  initv = true
  n = 5
  j = 2
  v = zeros(n,j+1)
  v[:,1] .= 1/sqrt(n)
  ldv = n
  resid = 2*rand(n).-1
  rnorm = Ref{Float64}(0.0)
  ipntr = zeros(Int, 3)
  workd = zeros(2n)
  
  seqdata = _check_dgetv0_sequence!(1.0*I;  idostart=0, 
    rnorm, resid, n, j, v, bmat, ldv, itry, initv)

  @test seqdata[end].ierr == 0    

  @test isapprox(sum(resid), 0.0; atol=n*eps(1.0))
  @test rnorm[] > 0
end

@testset "Simple Restarted Eigenproblem that fails" begin
  using LinearAlgebra
  using Random
  Random.seed!(0)
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

  seqdata = _check_dgetv0_sequence!(1.0*I;  idostart=0, 
    rnorm, resid, n, j, v, bmat, ldv, itry, initv)

  @test seqdata[end].ierr == -1    
  @test length(seqdata) == 1
end


## OLD CODE with unrolled 

#=
@testset "Restarted case where bmat is G and singular" begin
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

  # allocate a separate copy of Arpack
  arido = Ref{Int}(0)
  arv = copy(v)
  arresid = copy(resid0)
  arrnorm = Ref{Float64}(0.0)
  aripntr = zeros(Int, 3)
  arworkd = zeros(2n)

  ierr = -1
  state = ArpackInJulia.ArpackState{Float64}()
  while ido[] != 99
    #ierr = ArpackInJulia.dgetv0!(

    ierr = ArpackInJulia.dgetv0!(
      ido, Val{bmat}, itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state)
    
    arierr = arpack_dgetv0!(
      arido, bmat, itry, initv, n, j, arv, ldv, arresid, arrnorm, aripntr, arworkd)

    @debug arido[], arierr, ido[], ierr

    @test arierr == ierr
    @test arido[] == ido[]
    @test arv == v
    @test arresid == resid
    @test arrnorm[] ≈ rnorm[]
    @test aripntr == ipntr
    @test arworkd == workd  
    if ido[] == -1
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    elseif ido[] == 2
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    end
    # everything should be same...
    copyto!(arworkd, workd)
  end
  @test ierr==-1
end
=#