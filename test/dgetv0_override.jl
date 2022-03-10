## Here, we are going to try overriding dgetv0 with a state variable.
using Test
@testset "Checked state" begin 
  import ArpackInJulia

  Base.@kwdef mutable struct Override_dgetv0_ArpackState{T} <: ArpackInJulia.AbstractArpackState{T}
    aitr::ArpackInJulia.AitrState{T} = ArpackInJulia.AitrState{T}()
    getv0::ArpackInJulia.Getv0State{T} = ArpackInJulia.Getv0State{T}()
    saup2::ArpackInJulia.Saup2State{T} = ArpackInJulia.Saup2State{T}()
    aupd_nev0 = Ref{Int}(0)
    aupd_np = Ref{Int}(0)
    aupd_mxiter = Ref{Int}(0)
    extravar = Int(52)
  end


  ntimes = 0 


  function ArpackInJulia.dgetv0!(
    ido::Ref{Int}, # input/output
    ::Val{BMAT},
    itry::Int, # input
    initv::Bool,
    n::Int,
    j::Int,
    V::AbstractMatrix{T},
    ldv::Int, # TODO, try and remove
    resid::AbstractVecOrMat{T},
    rnorm::Ref{T}, # output
    ipntr::AbstractVector{Int}, # output
    workd::AbstractVector{T}, # output
    state::Override_dgetv0_ArpackState{T};
    stats::Union{ArpackInJulia.ArpackStats,Nothing}=nothing,
    debug::Union{ArpackInJulia.ArpackDebug,Nothing}=nothing,
    idonow::Union{ArpackInJulia.ArpackOp,Nothing}=nothing
    ) where {T, BMAT}

    ntimes::Int += 1 # use global variable 
    normalstate = ArpackInJulia.ArpackState{T}()
    normalstate.getv0 = state.getv0

    return ArpackInJulia.dgetv0!(ido, Val(BMAT), itry, initv, n, j, V, ldv, resid, rnorm, ipntr, workd,
      normalstate; debug, idonow, stats
    )
    state.getv0 = normalstate.getv0
  end 

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

  @test ntimes == 0 
  state=Override_dgetv0_ArpackState{Float64}()
  #state=ArpackInJulia.ArpackState{Float64}()
  ierr = ArpackInJulia.dgetv0!(
    ido, Val(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd, state
  )

  @test ntimes == 1
  @test state.getv0.iseed[] != tuple(1,3,5,7)
  @test ido[] == 99
  @test ierr==0
end 