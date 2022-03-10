#= 
This file handles functions related to 
Arpack Operators that can handle matvec ops
without the reverse communication interface
in Julia. 
=#

using LinearAlgebra: mul!

#=
c  Mode 1:  A*x = lambda*x, A symmetric
c           ===> OP = A  and  B = I.
=#
abstract type ArpackOp end

struct ArpackSimpleOp{MatT} <: ArpackOp
  A::MatT
end
opx!(y,OP::ArpackSimpleOp,x) = mul!(y,OP.A,x)      
is_arpack_mode_valid_for_op(mode::Int, ::ArpackSimpleOp) = mode == 1 

function _i_do_now_opx_1!(idonow::OpT, ipntr, workd, n) where {OpT <: ArpackOp} # handle the operation
  opx!(@view(workd[ipntr[2]:ipntr[2]+n-1]),idonow,@view(workd[ipntr[1]:ipntr[1]+n-1]))
end

function _i_do_now_bx!(idonow::OpT, ipntr, workd, n) where {OpT <: ArpackOp} # handle the operation
  bx!(@view(workd[ipntr[2]:ipntr[2]+n-1]),idonow,@view(workd[ipntr[1]:ipntr[1]+n-1]))
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
_i_do_now_shift!(idonow, np[], ipntr, workl)
function _i_do_now_shifts!(idonow::OpT, np::Int, ipntr, workl) where {OpT <: ArpackOp}
  offsets = ipntr[11]:ipntr[11]+np-1
  # TODO, handle eigenvalues and ritz estimates
  shifts!(@view(workl[offsets]), idonow)
end 