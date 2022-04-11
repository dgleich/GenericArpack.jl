using LinearAlgebra: Symmetric, Hermitian, HermOrSym

struct ArpackException <: Exception
  msg::String
end 

_float_type(::Type{T}) where T = typeof(one(T)/one(T))
_real_type(::Type{Complex{T}}) where T = _float_type(T)
_real_type(::Type{T}) where T = _float_type(T)

eigs(A::HermOrSym, k::Integer; kwargs...) = eigs(_float_type(eltype(A)), A, k; kwargs...)
# the next call must be special cased to make sure TV <: Complex for hermitian...
# this gives a super helpful error message!
eigs(::Type{TV}, A::Hermitian, k::Integer; kwargs...) where {TV <: Complex} = eigs(
  TV, _real_type(TV), A, k; kwargs...)
eigs(TV::Type, A::Symmetric, k::Integer; kwargs...) = eigs(TV, _real_type(TV), A, k; kwargs...)
eigs(TV::Type, TF::Type, A::HermOrSym, k::Integer; kwargs...) = symeigs(TV, TF, ArpackSimpleOp(A), k; kwargs...)

_gen_eigs_vtype(A, B) = promote_type(_float_type(eltype(A)), _float_type(eltype(B)))
_gen_eigs_ftype(A, B) = promote_type(_real_type(eltype(A)), _real_type(eltype(B)))
eigs(A::HermOrSym, B::HermOrSym, k::Integer; kwargs...) = symeigs(
  _gen_eigs_vtype(A,B), _gen_eigs_ftype(A,B), ArpackSymmetricGeneralizedOp(A,B), k; kwargs...)
  
eigs(A::Symmetric, B::Symmetric, k::Integer; kwargs...) = symeigs(T, ArpackSymmetricGeneralizedOp(A), k; kwargs...)
eigs(T::Type, A::Symmetric, B::Symmetric, k::Integer; kwargs...) = symeigs(T, ArpackSymmetricGeneralizedOp(A), k; kwargs...)

symeigs(A::AbstractMatrix,k::Integer; kwargs...) = symeigs(ArpackSimpleOp(A),k; kwargs...)
symeigs(TV::Type, A::AbstractMatrix,k::Integer; kwargs...) = symeigs(TV, ArpackSimpleOp(A),k; kwargs...)
symeigs(A::AbstractMatrix,B::AbstractMatrix,k::Integer; kwargs...) = symeigs(ArpackSymmetricGeneralizedOp(A,B),k; kwargs...)
symeigs(TV::Type, A::AbstractMatrix,B::AbstractMatrix,k::Integer; kwargs...) = symeigs(
    TV, ArpackSymmetricGeneralizedOp(A,B),k; kwargs...)

hermeigs(A::AbstractMatrix,k::Integer; kwargs...) = hermeigs(ArpackSimpleOp(A),k; kwargs...)
hermeigs(TV::Type, A::AbstractMatrix,k::Integer; kwargs...) = hermeigs(TV, ArpackSimpleOp(A),k; kwargs...)
hermeigs(A::AbstractMatrix,B::AbstractMatrix,k::Integer; kwargs...) = hermeigs(ArpackSymmetricGeneralizedOp(A,B),k; kwargs...)
hermeigs(TV::Type, A::AbstractMatrix,B::AbstractMatrix,k::Integer; kwargs...) = hermeigs(
    TV, ArpackSymmetricGeneralizedOp(A,B),k; kwargs...)

symeigs(op::ArpackOp, k::Integer; kwargs...) = symeigs(Float64, op, k; kwargs...)
hermeigs(op::ArpackOp, k::Integer; kwargs...) = symeigs(ComplexF64, op, k; kwargs...)
symeigs(TV::Type, op::ArpackOp, k::Integer; kwargs...) = symeigs(TV, _real_type(TV), op, k; kwargs...)
hermeigs(TV::Type, op::ArpackOp, k::Integer; kwargs...) = symeigs(TV, _real_type(TV), op, k; kwargs...)

struct ArpackEigen{TL,TV,OpT,StateT,BMAT}
  which::Symbol
  ipntr::Vector{Int}
  iparam::Vector{Int}
  V::Matrix{TV}
  workd::Vector{TV}
  workl::Vector{TL} 
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

function symeigs(::Type{TV}, ::Type{TF}, op::ArpackOp, nev::Integer; 
  which::Symbol=:LM, 
  maxiter::Int=max(300, ceil(Int,sqrt(size(op)))),
  bmat=bmat(op),
  ncv::Int=min(size(op),2nev),
  mode=arpack_mode(op),
  v0=nothing,
  stats=nothing,
  debug=nothing,
  state=ArpackState{TF}(),
  tol::TF=zero(TF),
  ritzvec::Bool = true,
  failonmaxiter::Bool = true) where {TV,TF} 

  n = size(op) 

  ncv <= n || throw(ArgumentError("ncv=$ncv must be at most n=$n"))
  nev < ncv || throw(ArgumentError("nev=$nev must be strictly less than ncv=$ncv"))
  
  isreal = TV <: Real 
  iscomplex = TV <: Complex
  isother = !isreal && !iscomplex 
  @assert(isother == false) # not yet implemented... 

  if !is_arpack_mode_valid_for_op(mode,op)
    @warn("mode $mode is not valid for op with type $(typeof(op))")
  end 

  TL = TF 

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
  workl = Vector{TF}(undef, lworkl)

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
  ierr = simple_dseupd!(ritzvec, select, values, vectors, shift(TF, op), 
    bmat, n, which, nev, tol, resid, ncv, V, iparam, ipntr, 
    workd, workl; debug, stats)

  if ierr != 0 
    throw(ArpackException("symmetric eupd gave error code ierr=$ierr"))
  end 

  return ArpackEigen(which, ipntr, iparam, V, workd, workl, resid, values, vectors, bmat, op, state)
end

svds(A::AbstractMatrix{T}, k::Integer; kwargs...) where T = svds(T, A, k; kwargs...)
svds(T::Type, A::AbstractMatrix, k::Integer; kwargs...) = svds(T, ArpackNormalOp(A), k; kwargs...)
svds(op::ArpackOp, k::Integer; kwargs...) = svds(Float64, op, k; kwargs...)
function svds(T::Type, op::ArpackOp, k::Integer; kwargs...)
  einfo = symeigs(T, op, k; kwargs...) # eigeninfo
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