using LinearAlgebra: Symmetric, Hermitian, HermOrSym, factorize, SVD
import Base: show

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
#_gen_eigs_ftype(A, B) = promote_type(_real_type(eltype(A)), _real_type(eltype(B)))
eigs(A::HermOrSym, B::HermOrSym, k::Integer; kwargs...) = symeigs(
  _gen_eigs_vtype(A,B), ArpackSymmetricGeneralizedOp(A,factorize(B),B), k; kwargs...)
  
eigs(T::Type, A::Symmetric, B::Symmetric, k::Integer; kwargs...) = symeigs(
    T, ArpackSymmetricGeneralizedOp(A,factorize(B),B), k; kwargs...)

symeigs(A::AbstractMatrix{T},k::Integer; kwargs...) where {T <: Complex} = symeigs(ComplexF64, ArpackSimpleOp(A),k; kwargs...)
symeigs(A::AbstractMatrix{T},k::Integer; kwargs...) where {T <: Real} = symeigs(Float64, ArpackSimpleOp(A),k; kwargs...)

symeigs(TV::Type, A::AbstractMatrix,k::Integer; kwargs...) = symeigs(TV, ArpackSimpleOp(A),k; kwargs...)
symeigs(A::AbstractMatrix,B::AbstractMatrix,k::Integer; kwargs...) = symeigs(
    ArpackSymmetricGeneralizedOp(A,factorize(B),B),k; kwargs...)
symeigs(TV::Type, A::AbstractMatrix, B::AbstractMatrix,k::Integer; kwargs...) = symeigs(
    TV, ArpackSymmetricGeneralizedOp(A,factorize(B),B),k; kwargs...)

hermeigs(A::AbstractMatrix,k::Integer; kwargs...) = hermeigs(ArpackSimpleOp(A),k; kwargs...)
hermeigs(TV::Type, A::AbstractMatrix,k::Integer; kwargs...) = hermeigs(TV, ArpackSimpleOp(A),k; kwargs...)
hermeigs(A::AbstractMatrix,B::AbstractMatrix,k::Integer; kwargs...) = hermeigs(
    ArpackSymmetricGeneralizedOp(A,factorize(B),B),k; kwargs...)
hermeigs(TV::Type, A::AbstractMatrix,B::AbstractMatrix,k::Integer; kwargs...) = hermeigs(
    TV, ArpackSymmetricGeneralizedOp(A,factorize(B),B),k; kwargs...)

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

function Base.getproperty(S::ArpackEigen, d::Symbol)
  if d === :bounds
    ncv = size(S.V,2)
    nconv = length(S.values) 
    return @view(S.workl[5ncv .+ (1:nconv)])
  else
    return getfield(S, d)
  end
end

Base.propertynames(F::ArpackEigen, private::Bool=false) =
  private ? (:bounds, fieldnames(typeof(F))...) : (:values, :vectors, :bounds)

function show(io::IO, mime::MIME{Symbol("text/plain")}, F::ArpackEigen)
  summary(io, F); println(io)
  println("eigenspace: ", F.which, )
  println(io, "values:")
  show(io, mime, F.values)
  if size(F.vectors,2) > 0
    println(io, "\nvectors:")
    show(io, mime, F.vectors)
  end 
end

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

  if (info == 1 && failonmaxiter)
    throw(ArpackException("symmetric aupd hit maxiter=$maxiter with $(iparam[5]) < $nev eigenvalues converged"))
  elseif ((info != 0 && info != 1)|| ido[] != 99) 
    throw(ArpackException("symmetric aupd gave error code info=$info with ido=$(ido[])"))
  end 

  # allocate output... 
  nconv = iparam[5]
  nvectors = ritzvec ? min(nev,nconv) : 0
  vectors = Matrix{TV}(undef, n, nvectors)
  values = Vector{TL}(undef, nconv)
  select = Vector{Int}(undef, ncv)
  ierr = simple_dseupd!(ritzvec, select, values, vectors, shift(TF, op), 
    bmat, n, which, nev, tol, resid, ncv, V, iparam, ipntr, 
    workd, workl; debug, stats)

  if ierr != 0 
    throw(ArpackException("symmetric eupd gave error code ierr=$ierr"))
  end 

  return ArpackEigen(which, ipntr, iparam, V, workd, workl, resid, values, vectors, bmat, op, state)
end

svds(A::AbstractMatrix{T}, k::Integer; kwargs...) where T = svds(_float_type(T), A, k; kwargs...)
svds(T::Type, A::AbstractMatrix, k::Integer; kwargs...) = svds(T, ArpackNormalOp(T,A), k; kwargs...)
svds(TV::Type, TF::Type, A::AbstractMatrix, k::Integer; kwargs...) = svds(TV, TF, ArpackNormalOp(TV, A), k; kwargs...)

svds(op::ArpackSVDOp, k::Integer; kwargs...) = svds(Float64, op, k; kwargs...)
complexsvds(op::ArpackSVDOp, k::Integer; kwargs...) = svds(ComplexF64, op, k; kwargs...)
svds(T::Type, op::ArpackSVDOp, k::Integer; kwargs...) = svds(T, _real_type(T), op, k; kwargs...)

function svds(TV::Type, TF::Type, op::ArpackSVDOp, k::Integer; which::Symbol=:LM, kwargs...)
  # this translates into the appropriate arpack parameters for the svd computation
  arnev, arwhich = _adjust_nev_which_for_svd(op, k, which)
  m, n = _uv_size(op)
  einfo = symeigs(TV, TF, op, arnev; kwargs..., which=arwhich) # eigeninfo
  _adjust_eigenvalues_to_singular_values!(einfo.values, op, arwhich, einfo.values)
  # re-sort the info... 
  nconv = length(einfo.values)
  nvecs = size(einfo.vectors,2)
  applysort = nvecs > 0 # I guess  
  if which==:LM || which==:LA # they want largest first, so we sort arpack info by smallest to reverse
    dsortr(:SM, applysort, nconv, einfo.values, einfo.vectors)
  elseif which == :SA || which == :LM # they want smallest first, so sort arpack info by largest
    dsortr(:LM, applysort, nconv, einfo.values, einfo.vectors)
  elseif which == :BE # ugh... this is annoying
    # don't sort here, arpack sorted them 
    # in smallest-to-largest for us. 
  end 
  
  U = Matrix{TV}(undef, m, nvecs)
  V = Matrix{TV}(undef, n, nvecs)
  if size(einfo.vectors,2) > 0 # then we have eigenvector info...
    eigenvecs_to_singvecs!(U, V, op, einfo.vectors)
  end
  return SVD(U, einfo.values, V')
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

"""
    svd_residuals(A, U, s, V)
    svd_residuals(A, SVD)
    svd_residuals(A, U, s, V, k) # compute for only the top k 
    svd_residuals!(r, A, U, s, V) # write result in place

Compute the residuals of an SVD computation \$||A*V[:,i] - U[:,i]*sigma[i]||\$,
and return the result in a vector. 

Using a matvec function 
-----------------------
Note that A can also be a function A(y,x), where y = A*x is updated 
in place. e.g. `svd(A,U,s,V) == svd((y,x)->mul!(y,A,x), U,s,V)`
are equivalent.
"""
:svd_residuals, :svd_residuals!

function svd_residuals!(r, A, U, s, V)
  k = length(r) 
  @assert(k <= size(U,2))
  @assert(k <= size(V,2))
  @assert(k <= length(s))
  m = size(U,1)
  av = Vector{eltype(U)}(undef, m) 
  for i in eachindex(r) 
    if A isa Function 
      A(av, @view(V[:,i]))
    else
      mul!(av, A, @view(V[:,i]))
    end 
    LinearAlgebra.BLAS.axpy!(-s[i], @view(U[:,i]), av)
    r[i] = norm(av)
  end 
  return r
end 
svd_residuals(A, U, s, V, k=length(s)) = svd_residuals!(Vector{eltype(U)}(undef, k), A, U, s, V)
svd_residuals(A, info::SVD, k=length(info.S)) = svd_residuals(A, info..., k)