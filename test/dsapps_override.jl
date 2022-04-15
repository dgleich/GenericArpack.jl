module CheckDsappsWithArpackjll
  using Test
  using GenericArpack
  import GenericArpack: AitrState, Getv0State, Saup2State, AbstractArpackState
  using LinearAlgebra

  # this 
  Base.@kwdef mutable struct CheckDsappsState{T} <: AbstractArpackState{T}
    aitr::AitrState{T} = AitrState{T}()
    getv0::Getv0State{T} = Getv0State{T}()
    saup2::Saup2State{T} = Saup2State{T}()
    aupd_nev0 = Ref{Int}(0)
    aupd_np = Ref{Int}(0)
    aupd_mxiter = Ref{Int}(0)
    aup2_rnorm = Ref{T}(zero(T))
  end 

  function GenericArpack.dsapps!(
    n::Int,
    kev::Int,
    np::Int,
    shift::AbstractVecOrMat{T},
    V::AbstractMatrix{T},
    ldv::Int,
    H::AbstractMatrix{T},
    ldh::Int,
    resid::AbstractVecOrMat{T},
    Q::AbstractMatrix{T},
    ldq::Int,
    workd::AbstractVecOrMat{T},
    state::CheckDsappsState{T},
    ;
    stats::Union{ArpackStats,Nothing}=nothing,
    debug::Union{ArpackDebug,Nothing}=nothing,
  ) where {T, BMAT}
    @debug "In override dsapps"
    # these codes won't work with idonow. 
    normalstate = GenericArpack.ArpackState{T}()

    arshift = copy(shift)
    arV = copy(V)
    arH = copy(H)
    arQ = copy(Q)
    arresid = copy(resid)
    arworkd = copy(workd) 

    # run arpack... 
    arrval = Main.arpack_dsapps!(n, kev, np, arshift, arV, ldv, arH, ldh, arresid, arQ, ldq, arworkd)

    # run ours
    rval = GenericArpack.dsapps!(n, kev, np, shift, V, ldv, H, ldh, resid, Q, ldq, workd, normalstate; stats, debug)

    @test arrval == rval
    @test arshift == shift 
    @test arresid == resid 
    @test arQ == Q
    @test arworkd == workd 
    @test arH == H
    @test arV == V

    return rval 
  end
end

using LinearAlgebra
function mysimpleeigvals_check_dsapps(op::ArpackOp, nev::Int = 6)
  ido = Ref{Int}(0)
  bmat = :I
  n = size(op.A,1)
  which = :LM
  tol = eps(Float64)/2 # just use the default
  resid = zeros(n)
  ncv = min(2nev, n-1)
  V = zeros(n,ncv)
  ldv = n
  mode = 1 
  iparam = zeros(Int,11)
  iparam[1] = 1
  iparam[3] = 300 # max iteration
  iparam[4] = 1
  iparam[7] = mode 
  ipntr = zeros(Int,11)
  workd = zeros(n,3)
  lworkl = ncv*ncv + 8*ncv
  workl = zeros(lworkl)

  iparam0 = copy(iparam)
  info_initv = 0

  T = Float64
  histdata = Vector{
    NamedTuple{(:ido, :resid, :V, :iparam, :ipntr, :workd, :workl, :ierr), 
      Tuple{Int,Vector{T}, Matrix{T}, Vector{Int}, Vector{Int}, Matrix{T}, Vector{T}, Int}}
  }()

  
  # make sure we reset this seed...
  _reset_libarpack_dgetv0_iseed()

  # Note that we cannot run two sequences at once and check them where we start a whole
  # second arpack call because of the expected Arpack state. 
  state = CheckDsappsWithArpackjll.CheckDsappsState{Float64}()
  stats = ArpackStats()
  debug = ArpackDebug()
  #debug = GenericArpack.set_debug_high!(ArpackDebug())
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
    push!(histdata, (;ido=ido[], resid=copy(resid), V=copy(V), iparam=copy(iparam), 
      ipntr=copy(ipntr), workd=copy(workd), workl=copy(workl), ierr))



    arierr = arpack_dsaupd!(arido, bmat, n, which, nev, tol, arresid, ncv, arV, ldv, ariparam, 
      aripntr, arworkd, arworkl, lworkl, info_initv)

    @test ido[] == arido[]
    @test workl == arworkl 
    @test resid == arresid
    @test ierr == arierr
    @test V == arV 
    @test ipntr == aripntr
    @test iparam == ariparam
    @test workd == arworkd 

    if ido[] == 1 || ido[] == -1
      GenericArpack._i_do_now_opx_1!(op, ipntr, workd, n)
    elseif ido[] == 99
      break
    else
      @error("this only supports standard eigenvalue problems")
    end 
  end
end

mysimpleeigvals_check_dsapps(GenericArpack.ArpackSimpleOp(Diagonal(1.0:10)))