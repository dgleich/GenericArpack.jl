using LinearAlgebra: Symmetric

struct ArpackException <: Exception
  msg::String
end 

eigs(A::Symmetric, k::Integer; kwargs...) = eigs(eltype(A), A, k; kwargs...)
eigs(T::Type, A::Symmetric, k::Integer; kwargs...) = symeigs(T, ArpackSimpleOp(A), k; kwargs...)
eigs(T::Type, A::Symmetric, B::Symmetric, k::Integer; kwargs...) = symeigs(T, ArpackSymmetricGeneralizedOp(A), k; kwargs...)
symeigs(op::ArpackOp, k::Integer; kwargs...) = symeigs(Float64, op, k; kwargs...)

struct ArpackEigen{TL,TV,OpT,StateT,BMAT}
  ipntr::Vector{Int}
  iparam::Vector{Int}
  V::Matrix{TV}
  workd::Vector{TV}
  workl::Vector{TV} 
  resid::Vector{TV}
  values::Vector{TL}
  vectors::Matrix{TV}
  bmat::Val{BMAT}
  op::OpT
  state::StateT
end
Base.iterate(S::ArpackEigen) = (S.values, Val(:vectors))
Base.iterate(S::ArpackEigen, ::Val{:vectors}) = (S.vectors, Val(:done))
Base.iterate(S::ArpackEigen, ::Val{:done}) = nothing

function symeigs(::Type{T}, op::ArpackOp, nev::Integer; 
  which::Symbol=:LM, 
  maxiter::Int=max(300, ceil(Int,sqrt(size(op)))),
  bmat=bmat(op),
  ncv::Int=min(size(op),2nev),
  mode=arpack_mode(op),
  v0=nothing,
  stats=nothing,
  debug=nothing,
  state=ArpackState{T}(),
  tol::T=zero(T),
  ritzvec::Bool = true,
  failonmaxiter::Bool = true) where T 

  n = size(op) 

  ncv <= n || throw(ArgumentError("ncv=$ncv must be at most n=$n"))
  nev < ncv || throw(ArgumentError("nev=$nev must be strictly less than ncv=$ncv"))
  
  isreal = T <: Real 
  iscomplex = T <: Complex
  isother = !isreal && !iscomplex 
  @assert(iscomplex == false)
  @assert(isother == false) # not yet implemented... 

  # goal is to run arpack with type T for a symmetric problem
  if isreal
    TV = T
    TL = T # lambda are 
  end
  @assert(sizeof(TV) >= sizeof(TL)) # this means we can always store a TL in an TV spot... 

  nvectors = ritzvec ? nev : 0 

  # allocate workspace... 
  vectors = Matrix{TV}(undef, n, nvectors)
  values = Vector{TL}(undef, nev)
  V = Matrix{TV}(undef, n, ncv)
  resid = Vector{TV}(undef, n)
  iparam = zeros(Int, 11)
  ipntr = zeros(Int, 11)
  workd = Vector{TV}(undef, 3n)
  lworkl = ncv*ncv + 8*ncv
  workl = Vector{TV}(undef, lworkl)

  iparam[1] = 1
  iparam[4] = 1
  iparam[3] = maxiter 
  iparam[7] = mode 

  info_initv = 0 
  if v0 !== nothing
    copyto!(resid, v0)
    info_initv = 1 # v is already initialized!
  end

  ido = Ref{Int}(0)
  info = dsaupd!(ido, bmat, n, which, nev, tol, resid, ncv, V, size(V,1), 
    iparam, ipntr, workd, workl, lworkl, info_initv;
    state, stats, debug, idonow=op)[1]

  if (info == 1 && failonmaxiter) || (info != 0 || ido[] != 99)
    throw(ArpackException("symmetric aupd gave error code info=$info with ido=$(ido[])"))
  end 

  select = Vector{Int}(undef, ncv)
  ierr = simple_dseupd!(ritzvec, select, values, vectors, shift(T, op), 
    bmat, n, which, nev, tol, resid, ncv, V, iparam, ipntr, 
    workd, workl; debug, stats)

  if ierr != 0 
    throw(ArpackException("symmetric eupd gave error code ierr=$ierr"))
  end 

  return ArpackEigen(ipntr, iparam, V, workd, workl, resid, values, vectors, bmat, op, state)
end

#=
eigs(A::Symmetric, k; kwargs...)
eigs(T, A::Symmetric, k; kwargs...) # force Arpack to use type T. 
eigs(T, A::Symmetric, B::Symmetric, k; kwargs...) # force Arpack to use type T. 
eigs(A::ArpackSymmetricOp, k; kwargs...) -> ArpackEigen()

svds(A::AbstractMatrix{T}, k; kwargs...) where {T <: AbstractFloat}
svds(A::ArpackSVDOp, k; kwargs...)
svds(T, A::ArpackSVDOp, k; kwargs...)

##
struct ArpackEigen{T,TL,TV}
  ipntr::Vector{Int}
  iparam::Vector{Int}
  V::Matrix{TV}
  state::ArpackState{T}
  workd::Vector{TV}
  workl::Vector{TV} 
  resid::Vector{TV}
  values::Vector{TL}
  vectors::Matrix{TV}
end
isconverged()
Base.iterate(S::ArpackEigen) = (S.values, Val(:vectors))
Base.iterate(S::ArpackEigen, ::Val{:vectors}) = (S.vectors, Val(:done))
Base.iterate(S::ArpackEigen, ::Val{:done}) = nothing

# in the future... 
eigs(A::Hermitian, k; kwargs...)
eigs(T, A::Hermitian, k; kwargs...) # force Arpack to use type T. 
svds()
=#