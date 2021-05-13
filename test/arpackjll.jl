import Arpack_jll, LinearAlgebra

function arpack_dsconv(n::Int, ritz::Vector{Float64}, bounds::Vector{Float64},
                      tol::Float64)
  nconv::LinearAlgebra.BlasInt = -1
  ccall((:dsconv_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ptr{Float64},
     Ref{Float64},
     Ref{LinearAlgebra.BlasInt}),
    n, ritz, bounds, tol, nconv)
  return nconv
end
