#= 
This file handles functions related to 
Arpack Operators that can handle matvec ops
without the reverse communication interface
in Julia. 

Cases:  
-- symmetric A: mode 1  --> Ax = lambda x, OP = A, B = I 
-- symmetric A: mode 2 (special!) --> Ax = lambda Mx, OP = inv(M)*A, B = M, (M spd) 
   but need to write Ax into memory during OP*x too. y = OP*x -> x -> A*x, y -> inv(M)*x 
   so no need for extra vector :) thanks Arpack team! 
-- symmetric A: mode 3 --> K*x = lambda*M*x, OP = (inv[K - sigma*M])*M, B = M (M spsd)
-- symmetric A: mode 4 --> K*x = lambda KGx, OP = (inv(K-sigma*KG))*K, B = K, (K spsd, KG indef)
-- symmetric A: mode 5 --> A*x = lambda*M*x, OP = inv(A-sigma*M)*(A + sigma M), B = M (M spsd)

-- complex, non-sym: mode 1 --> Ax = lambda x, OP = A, B = I 
-- complex, non-sym: mode 2 --> Ax = lambda Mx, OP = inv(M)*A, B = M (M hpd)
-- complex, non-sym: mode 3 --> Ax = lambda Mx, OP = inv(A - sigma M)*M, B = M (M hpsd)

Cases: 
ido[] == -1 --> idonow_opx_neg1 --> opx!(y, OP, x) # only to overwrite x here... Bx not avail.
ido[] ==  1 --> idonow_opx_1 -> opx!(y, OP, x)
ido[] ==  1 --> idonow_opx_mode2_1! -> genopx!(y, OP, x)     # sym, mode 2 only! 
ido[] ==  1 --> idonow_shiftinvert_opx! -> opx!(y,OP,x,Bx)   # mode 3, 4, 5 only! 
ido[] ==  2 --> idonow_bx -> bx!(y, OP, x)
ido[] ==  3 --> idonow_shifts -> shifts!(lams, OP))
=#

import Base: Matrix

abstract type ArpackOp end
abstract type ArpackSVDOp <: ArpackOp end

# general abstract cases... 
shift(T::Type, ::ArpackOp) = zero(T) # default shift is none!
opx!(y,OP::ArpackOp,x,Bx) = opx!(y,OP,x) # default shift-invert operator. 
shifts!(lams, OP::ArpackOp) = nothing # this is just a no-op unless you implement it yourself! 
function Matrix(T::Type, op::ArpackOp) 
  n = size(op)
  v = Vector{T}(undef, n) 
  A = Matrix{T}(undef, n, n)
  for i=1:n 
    fill!(v, 0) # in some cases, opx! can overwrite v... 
    v[i] = one(T) 
    opx!(@view(A[:,i]), op, v)
  end 
  return A
end 

## These map from the raw codes to more user-friendly functions. 
function _i_do_now_opx_neg1!(idonow::OpT, ipntr, workd, n) where {OpT <: ArpackOp} # handle the operation
  opx!(@view(workd[ipntr[2]:ipntr[2]+n-1]),idonow,@view(workd[ipntr[1]:ipntr[1]+n-1]))
end

function _i_do_now_opx_1!(idonow::OpT, ipntr, workd, n) where {OpT <: ArpackOp} # handle the operation
  opx!(@view(workd[ipntr[2]:ipntr[2]+n-1]),idonow,@view(workd[ipntr[1]:ipntr[1]+n-1]))
end

function _i_do_now_opx_mode2_1!(idonow::OpT, ipntr, workd, n) where {OpT <: ArpackOp} # handle the operation
  genopx!(@view(workd[ipntr[2]:ipntr[2]+n-1]),idonow,@view(workd[ipntr[1]:ipntr[1]+n-1]))
end

function _i_do_now_opx_shiftinvert_1!(idonow::OpT, ipntr, workd, n) where {OpT <: ArpackOp} # handle the operation
  opx!(@view(workd[ipntr[2]:ipntr[2]+n-1]),idonow,@view(workd[ipntr[1]:ipntr[1]+n-1]),@view(workd[ipntr[3]:ipntr[3]+n-1]))
end

function _i_do_now_bx!(idonow::OpT, ipntr, workd, n) where {OpT <: ArpackOp} # handle the operation
  bx!(@view(workd[ipntr[2]:ipntr[2]+n-1]),idonow,@view(workd[ipntr[1]:ipntr[1]+n-1]))
end 

using LinearAlgebra: mul!, ldiv!

#=
c  Mode 1:  A*x = lambda*x, A symmetric
c           ===> OP = A  and  B = I.
=#

"""
    ArpackSimpleOp(A)

This corresponds to a simple eigenvalue problem Ax = lambda x, and builds
a Julia object that represents the minimal information Arpack needs about
the matrix to run an eigenvalue problem. 

# Arguments
- `A`: Anything that implements `Base.size`, `LinearAlgebra.mul!`

# Examples
```julia-repl
julia> op = ArpackSimpleOp(A)
juila> size(op) == size(A,1)
```

See also [`ArpackSymmetricGeneralizedOp`](@ref), [`ArpackSimpleFunctionOp`](@ref)
"""
struct ArpackSimpleOp{MatT} <: ArpackOp
  A::MatT
end
arpack_mode(::ArpackSimpleOp) = 1
Base.size(op::ArpackSimpleOp) = Base.size(op.A, 1)
bmat(::ArpackSimpleOp) = Val(:I)

opx!(y,OP::ArpackSimpleOp,x) = mul!(y,OP.A,x)      
is_arpack_mode_valid_for_op(mode::Int, ::ArpackSimpleOp) = mode == 1 


"""
    ArpackAugmentedOp

This corresponds to a single eigenvalue problem on the augmented matrix [0 A; A' 0]
to mirror ARPACK's SVD drivers.
"""
struct ArpackAugmentedOp{MatT} <: ArpackSVDOp
  A::MatT    
end
arpack_mode(::ArpackAugmentedOp) = 1
Base.size(op::ArpackAugmentedOp) = +(Base.size(op.A)...)
bmat(::ArpackAugmentedOp) = Val(:I)

function _adjust_nev_which_for_svd(op::ArpackAugmentedOp, nev::Integer, svdwhich::Symbol)
  if svdwhich==:LM || svdwhich == :LA
    return nev, :LA 
  elseif which == :SM || svdwhich == :SA 
    # need to offset for all the zero eigenvalues... 
    # m+n eigenvalues, with 2*min(m,n) sing. pairs
    # so m+n-2*min(m,n)
    return size(op)-2*minimum(Base.size(op.A)) + 2*nev, :SM
  else
    throw(ArgumentError("which = $svdwhich is not valid for SVD of ArpackAugmentedOp"))
  end 
end 

function _adjust_eigenvalues_to_singular_values!(svdvals, op::ArpackAugmentedOp, which::Symbol, eigenvals)
  if which == :LA
    length(svdvals) >= length(eigenvals) || throw(
        ArgumentError("incompatible lengths of svdvals and eigenvals for which=$which"))

    for i in eachindex(eigenvals)
      svdvals[i] = max(eigenvals[i],0) # truncate negative values to zero 
    end 

  elseif which == :SM
    for i in eachindex(eigenvals)
      svdvals[i] = max(eigenvals[i],0) # truncate negative values to zero 
    end 

  else
    throw(ArgumentError("which = $svdwhich is not valid for SVD of ArpackAugmentedOp"))
  end 
end

function opx!(y,OP::ArpackAugmentedOp,x)
  m,n = Base.size(OP.A)
  x1 = view(x,1:m)
  x2 = view(x,m+1:m+n)
  y1 = view(y,1:m)
  y2 = view(y,m+1:m+n)
  mul!(y1,OP.A,x2)
  mul!(y2,adjoint(OP.A),x1)
end
is_arpack_mode_valid_for_op(mode::Int, ::ArpackAugmentedOp) = mode == 1 

function eigenvecs_to_singvecs!(U::AbstractMatrix, V::AbstractMatrix, 
  OP::ArpackAugmentedOp, Z::AbstractMatrix
)  
  m,n = Base.size(OP.A)
  copyto!(U, Z[1:m,:])
  copyto!(V, Z[m+1:m+n,:])
  # TODO, check if these are allocating or not... 
  normalize!.(eachcol(V))
  normalize!.(eachcol(U))
end 

"""
    ArpackNormalOp(A)

This provides an [`ArpackSVDOp`](@ref) that is used as a bridge in `svds` 
to call the symeigs routine. It maniuplates the so-called normal matrix 
\$AA^H\$ or \$A^H A\$ depending on which 
is smaller. Note that this function allocates memory to compute the operation,
so it will not be safe to use with multiple threads. 

# Arguments
- `A`: Anything that implements `Base.size`, `LinearAlgebra.mul!(y,A,x)`,
  `LinearAlgebra.adjoint` and `LinearAlgebra.mul!(y,adjoint(A),x)``

# Examples
```julia-repl
julia> op = ArpackNormalOp(A)
```

See also [`ArpackNormalFunctionOp`](@ref)
"""
abstract type ArpackNormalSVDOp <: ArpackSVDOp end 

struct ArpackNormalOp{MatT,ET} <: ArpackNormalSVDOp
  A::MatT    
  n::Int
  At_side::Symbol 
  extra::Vector{ET}
end
function ArpackNormalOp(T,A)
  m,n = Base.size(A)
  if m <= n
    At_side = :right
    extradim = n 
    sz = m 
  else
    At_side = :left 
    extradim = m 
    sz = n 
  end 
  extra = Vector{T}(undef, extradim)
  return ArpackNormalOp(A, sz, At_side, extra)
end 
arpack_mode(::ArpackNormalOp) = 1
Base.size(op::ArpackNormalOp) = op.n
bmat(::ArpackNormalOp) = Val(:I)

function _adjust_nev_which_for_svd(op::ArpackNormalSVDOp, nev::Integer, svdwhich::Symbol)
  if svdwhich==:LM || svdwhich==:LA
    return nev, :LM
  elseif svdwhich == :SM || svdwhich == :SA
    # need to offset for all the zero eigenvalues... 
    # m+n eigenvalues, with 2*min(m,n) sing. pairs
    # so m+n-2*min(m,n)
    return nev, :SA 
  elseif svdwhich == :BE 
    return nev, :BE 
  else
    throw(ArgumentError("which = $svdwhich is not valid for SVD of ArpackNormalOp "))
  end 
end 

function _adjust_eigenvalues_to_singular_values!(svdvals, op::ArpackNormalSVDOp, which::Symbol, eigenvals)
  length(svdvals) >= length(eigenvals) || throw(ArgumentError("incompatible lengths of svdvals and eigenvals"))
  for i in eachindex(eigenvals)
    svdvals[i] = sqrt(max(eigenvals[i], 0))
  end 
end

function _uv_size(op::ArpackNormalOp)
  return size(op.A)
end 


function opx!(y,OP::ArpackNormalOp,x)
  if OP.At_side == :left
    # AtA
    mul!(OP.extra, OP.A, x)
    mul!(y, adjoint(OP.A), OP.extra)
  else
    mul!(OP.extra, adjoint(OP.A), x)
    mul!(y, OP.A, OP.extra)
  end 
end
is_arpack_mode_valid_for_op(mode::Int, ::ArpackNormalSVDOp) = mode == 1 

function orth!(X::AbstractMatrix)
  # should really run gram-schmidt?
  # run QR, then dorgqr
  #http://www.netlib.org/lapack/explore-html/da/dba/group__double_o_t_h_e_rcomputational_ga14b45f7374dc8654073aa06879c1c459.html#ga14b45f7374dc8654073aa06879c1c459
  # or dorg2r?

  # Qinfo = LinearAlgebra.qr!(X)
  # copyto!(X, Matrix(Qinfo))

  #signinfo = Vector{eltype(X)}(undef, size(X,2))
  tau = Vector{eltype(X)}(undef, size(X,2))
  work = Vector{eltype(X)}(undef, size(X,2))
  _dgeqr2!(X, tau, work)
  _sign(x) = iszero(x) ? one(x) : sign(x)
  rsign = [_sign(X[i,i]) for i in 1:size(X,2)]
  dorg2r!(X, tau, work)
  # scale the columns... can be done ourselves too
  LinearAlgebra.rmul!(X, LinearAlgebra.Diagonal(rsign))
end 

function eigenvecs_to_singvecs!(U::AbstractMatrix, V::AbstractMatrix, 
  OP::ArpackNormalOp, Z::AbstractMatrix
)
  # TODO, better handling for singular subspace... 
  warnthresh = _eps23(real(promote_type(eltype(U),eltype(V))))
  errorthresh = eps(real(promote_type(eltype(U),eltype(V))))/2
  check_norm(nrm, var) = begin
    if nrm <= errorthresh
      throw(ArpackException("encountered σ ≈ $(nrm) when computing $var subspace"))
    elseif nrm <= warnthresh
      @warn "$var subspace may be inaccurate due to small singular value (σ ≈ $(nrm)) ≤ ε^(2/3) = $warnthresh"
    end 
  end
  if OP.At_side == :left
    copyto!(V, Z)
    for i=1:size(V,2)
      mul!(@view(U[:,i]), OP.A, @view(V[:,i]))
      #LinearAlgebra.normalize!(@view(U[:,i]))
      nrm = norm(@view(U[:,i]))
      check_norm(nrm, "U")
      _scale_from_to(nrm, one(eltype(nrm)),@view(U[:,i]))
    end 
    orth!(U)
  else
    copyto!(U, Z)
    for i=1:size(V,2)
      mul!(@view(V[:,i]), adjoint(OP.A), @view(U[:,i]))
      #LinearAlgebra.normalize!(@view(V[:,i]))
      nrm = norm(@view(V[:,i]))
      check_norm(nrm, "V")
      _scale_from_to(nrm, one(eltype(nrm)),@view(V[:,i]))
    end 
    orth!(V)
  end 
end 

"""
    ArpackSimpleFunctionOp(F::Function, n::Integer)

This corresponds to a simple eigenvalue problem Ax = lambda x, but 
takes a functional operator that we apply.

# Arguments
- `F::Function` this is a function (y,x) -> mul!(y,A,x), i.e. a function 
  that writes \$A*x\$ into \$y\$
- `n::Integer` the dimension of the problem 

# Examples
```julia-repl
julia> using GenericArpack, SparseArrays
julia> function myf(y,x) 
       fill!(y, 0)
       y[1] += x[2] + x[100]
       for i in 2:99
         y[i] += x[i+1] + x[i-1]
       end
       y[100] += x[99] + x[1]
       end
julia> op = ArpackSimpleFunctionOp(myf, 100)
"""
struct ArpackSimpleFunctionOp <: ArpackOp
  F::Function 
  n::Int
end 
arpack_mode(::ArpackSimpleFunctionOp) = 1
Base.size(op::ArpackSimpleFunctionOp) = op.n
bmat(::ArpackSimpleFunctionOp) = Val(:I)
opx!(y,op::ArpackSimpleFunctionOp,x) = op.F(y,x)
is_arpack_mode_valid_for_op(mode::Int, ::ArpackSimpleFunctionOp) = mode == 1 


# This will try and determine if your matvec function is real or complex
# valued... 
function _autotype_function(T::Type, F::Function, m, n)
  y = Vector{T}(undef, m)
  x = ones(T, n) 
  try
    F(y,x) 
    return T
  catch e
    if e isa InexactError
      # see if we can get this to work with a ComplexT type...
      try
        y = Vector{Complex{T}}(undef, m)
        x = ones(Complex{T}, n)
        F(y,x) 
        return Complex{T}
      catch
        throw(ArgumentError("unable to auto-type the function as real or complex, please specify an explicit type"))
      end
    else
      rethrow(e) 
    end 
  end
end 
"""
    ArpackNormalFunctionOp(ax::Function, atx::Function, m::Integer, n::Integer)

# Arguments
- `ax` A function to compute (y,x) and write y = A*x into the memory for y
- `atx` A function to compute (y,x) and write y = A^H*x into the memory for y
- `m` the number of rows of A
- `n` the number of rows of A

See also [`ArpackNormalOp`](@ref), [`ArpackSimpleFunctionOp`](@ref)
"""
struct ArpackNormalFunctionOp{T} <: ArpackNormalSVDOp
  av::Function
  atv::Function
  m::Int
  n::Int 
  extra::Vector{T}
end 

ArpackNormalFunctionOp(av::Function, atv::Function, m::Int, n::Int) = ArpackNormalFunctionOp(
  _autotype_function(Float64, av, m, n), av, atv, m, n)

ArpackNormalFunctionOp(T::Type, av::Function, atv::Function, m::Int, n::Int) = ArpackNormalFunctionOp{T}(
  av, atv, m, n, Vector{T}(undef, m))

arpack_mode(::ArpackNormalFunctionOp) = 1
Base.size(op::ArpackNormalFunctionOp) = op.n
bmat(::ArpackNormalFunctionOp) = Val(:I)

function _uv_size(op::ArpackNormalFunctionOp)
  return op.m, op.n 
end 
function opx!(y,OP::ArpackNormalFunctionOp,x)
  OP.av(OP.extra, x)
  OP.atv(y, OP.extra) 
end
is_arpack_mode_valid_for_op(mode::Int, ::ArpackNormalFunctionOp) = mode == 1 

function eigenvecs_to_singvecs!(U::AbstractMatrix, V::AbstractMatrix, 
  OP::ArpackNormalFunctionOp, Z::AbstractMatrix
)
  copyto!(V, Z)
  for i=1:size(V,2)
    OP.av(@view(U[:,i]), @view(V[:,i]))
    LinearAlgebra.normalize!(@view(U[:,i]))
  end 
  orth!(U)
end 


"""
    ArpackSymmetricGeneralizedOp(A,invB,B)    

We need three operations: A*x, invB*x, B*x    
B must also be symmetric, pos. def. 
Note that if B can be factorized efficiently via Cholesky, then there is a better way to proceed. 

# Examples
```julia-repl
julia> using GenericArpack, LinearAlgebra
julia> n = 100 
julia> A = Tridiagonal(-ones(n-1),2*ones(n),-ones(n-1)).*(n+1)
julia> B = Tridiagonal(ones(n-1),4*ones(n),ones(n-1)).*(1/(6*(n+1)))
julia> op = ArpackSymmetricGeneralizedOp(A, lu!(copy(B)), B)

See also [`ArpackSimpleOp`](@ref)
"""    
struct ArpackSymmetricGeneralizedOp{MatT, SolveType, BType} <: ArpackOp
  A::MatT
  S::SolveType
  B::BType
end 
arpack_mode(::ArpackSymmetricGeneralizedOp) = 2
is_arpack_mode_valid_for_op(mode::Int, ::ArpackSymmetricGeneralizedOp) = mode == 2
bmat(::ArpackSymmetricGeneralizedOp) = Val(:G)

function genopx!(y, OP::ArpackSymmetricGeneralizedOp, x)
  mul!(y, OP.A, x) # compute A*x
  x[:] = y         # save A*x in x per remark 5 in
  try 
    ldiv!(OP.S, y)   # overwrite y with inv(B)*y 
  catch e
    #throw(e)
    if e isa MethodError
      # try one more time...
      # this fixes an issue with Cholmod factorizations. 
      y[:] = OP.S \ x
    else
      rethrow(e)
    end
  end 
  #ldiv!(y, OP.S, x) # overwrite y with inv(B)*y = inv(B)*x = inv(B)*A*x
end 
opx!(y, OP::ArpackSymmetricGeneralizedOp, x) = genopx!(y, OP, x) # we can use that same one 
Base.size(op::ArpackSymmetricGeneralizedOp) = Base.size(op.A, 1)
function bx!(y, OP::ArpackSymmetricGeneralizedOp, x)
  mul!(y, OP.B, x)
end 

"""
    ArpackShiftInvertOp

"""    
struct ArpackShiftInvertOp{MatT, SolveType, BType, FType} <: ArpackOp
  A::MatT
  S::SolveType
  B::BType
  sigma::FType
end 
Base.size(op::ArpackShiftInvertOp) = Base.size(op.A, 1)
is_arpack_mode_valid_for_op(mode::Int, ::ArpackShiftInvertOp) = mode == 3
shift(T, op::ArpackShiftInvertOp) = T(op.sigma) 

arpack_mode(::ArpackShiftInvertOp) = 3
function bmat(::ArpackShiftInvertOp{MatT,SolveT,BType,FType}
) where {MatT, SolveT, BType, FType}
  if BType <: UniformScaling
    return Val(:I)
  else
    return Val(:G)
  end 
end

function shiftinvert_op(A,B,sigma)
  S = factorize(A-sigma*B)
  return ArpackShiftInvertOp(A, S, B, sigma)
end 

function opx!(y, OP::ArpackShiftInvertOp, x, Bx)
  copyto!(y, Bx)
  ldiv!(OP.S, y)   # overwrite y with inv(A-sigma*B)*y 
end 

function opx!(y, OP::ArpackShiftInvertOp, x)
  mul!(y, OP.B, x) # compute B*x
  ldiv!(OP.S, y)   # overwrite y with inv(A-sigma*B)*y 
end 

function bx!(y, OP::ArpackShiftInvertOp, x)
  mul!(y, OP.B, x)
end 

#=
c          IDO =  3: compute the IPARAM(8) shifts where
c                    IPNTR(11) is the pointer into WORKL for
c                    placing the shifts. See remark 6 below.
c
c  6. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the
c     NP = IPARAM(8) shifts in locations:
c     1   WORKL(IPNTR(11))
c     2   WORKL(IPNTR(11)+1)
c                        .
c                        .
c                        .
c     NP  WORKL(IPNTR(11)+NP-1).
c
c     The eigenvalues of the current tridiagonal matrix are located in
c     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are in the
c     order defined by WHICH. The associated Ritz estimates are located in
c     WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... , WORKL(IPNTR(8)+NCV-1).
=#
function _i_do_now_shifts!(idonow::OpT, np::Int, ipntr, workl) where {OpT <: ArpackOp}
  offsets = ipntr[11]:ipntr[11]+np-1
  # TODO, handle eigenvalues and ritz estimates
  shifts!(@view(workl[offsets]), idonow)
end 