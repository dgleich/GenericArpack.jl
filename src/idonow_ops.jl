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
