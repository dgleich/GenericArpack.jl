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

abstract type ArpackOp end
# general abstract cases... 
shift(T::Type, ::ArpackOp) = zero(T) # default shift is none!
opx!(y,OP::ArpackOp,x,Bx) = opx!(y,OP,x) # default shift-invert operator. 
shifts!(lams, OP::ArpackOp) = nothing # this is just a no-op unless you implement it yourself! 

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
    ArpackSimpleOp

This corresponds to a simple eigenvalue problem Ax = lambda x.

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
    ArpackSymmetricGeneralizedOp(A,invB,B)    

We need three operations: A*x, invB*x, B*x    
B must also be symmetric, pos. def. 
Note that if B can be factorized efficiently via Cholesky, then there is a better way to proceed. 

Example: 
    using LinearAlgebra
    n = 100 
    A = Tridiagonal(-ones(n-1),2*ones(n),-ones(n-1)).*(n+1)
    B = Tridiagonal(ones(n-1),4*ones(n),ones(n-1)).*(1/(6*(n+1)))
    ArpackSymmetricGeneralizedOp(A, lu!(copy(B)), B)
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
  ldiv!(OP.S, y)   # overwrite y with inv(B)*y 
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
shift(T, op::ArpackShiftInvertOp) = T(sigma) 

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
  opM = A-sigma*B
  S = lu!(opM)
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