import Arpack_jll, LinearAlgebra

function arpack_dsconv(n::Int, ritz::Vector{Float64}, bounds::Vector{Float64},
                      tol::Float64)
  nconv = Ref{LinearAlgebra.BlasInt}(-1)
  ccall((:dsconv_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ptr{Float64},
     Ref{Float64},
     Ref{LinearAlgebra.BlasInt}),
    n, ritz, bounds, tol, nconv)
  return nconv[]
end

function arpack_dsortr(
  which::Symbol, # Input
  apply::Bool, # Input
  n::Int, # Input
  x1::Vector{Float64}, # Input/Output
  x2::Vector{Float64}, # Input/Output
  )
  if which==:LM; whichstr="LM"
  elseif which==:SM; whichstr="SM"
  elseif which==:LA; whichstr="LA"
  elseif which==:SA; whichstr="SA"
  end
  ccall((:dsortr_, Arpack_jll.libarpack), Cvoid,
    (Ptr{UInt8},
     Ref{LinearAlgebra.BlasInt}, # bool
     Ref{LinearAlgebra.BlasInt}, # size
     Ptr{Float64},
     Ref{Float64}),
    whichstr, apply, n, x1, x2)
end
