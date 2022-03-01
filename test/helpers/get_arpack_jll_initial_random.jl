## This code is helpful to get the initial residual vector out of Arpack_jll on the very first call.
# Note that because this has state, this will return different things, so it should only be called once.
# This results in the exit(1) at the end of the script!
# Don't remove this or you'll regret it.
# you can use the little script at the very end to run it!
import Arpack_jll, LinearAlgebra
function arpack_dgetv0!(ido::Ref{LinearAlgebra.BlasInt}, bmat::Symbol, itry::Int, initv::Bool,
  n::Int, j::Int, v::StridedMatrix{Float64}, ldv::Int,
    resid::StridedVecOrMat{Float64}, rnorm::Ref{Float64},
    ipntr::StridedVecOrMat{LinearAlgebra.BlasInt}, workd::StridedVecOrMat{Float64})
  ierr = Ref{LinearAlgebra.BlasInt}(0)
  # NOTE, arpack doesn't touch ierr unless ido[] == 0 or there is
  # a restart failure.
  @show ido[], itry, initv, n, j, string(bmat)
  ccall((:dgetv0_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt},
     Ptr{UInt8},
     Ref{LinearAlgebra.BlasInt},
     Ref{LinearAlgebra.BlasInt}, # really a logical...
     Ref{LinearAlgebra.BlasInt},
     Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ref{Float64},
     Ptr{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ref{LinearAlgebra.BlasInt}),
    ido, string(bmat), itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd, ierr)
  return ierr[]
end
function get_initial_vector(n::Int=10)
  ido = Ref{LinearAlgebra.BlasInt}(0)
  bmat = :I
  itry = 1
  initv = false
  j = 1
  v = zeros(n, j+1)
  ldv = n
  resid = zeros(n)
  rnorm = Ref{Float64}(0.0)
  ipntr = zeros(Int, 3)
  workd = zeros(2n)

  resid0 = copy(resid)

  ierr = arpack_dgetv0!(ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd)
  println.(resid)
  return nothing
end
get_initial_vector()
exit()
##
file = @__FILE__
run(`/Applications/Julia-1.7.app/Contents/Resources/julia/bin/julia $file`)