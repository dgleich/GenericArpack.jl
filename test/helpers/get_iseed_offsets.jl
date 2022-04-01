## This code is helpful to find the offsets to iseed
# so that we can reset the Arpack random number generator.
## get the offsets

## Try and automatically get them... 
import Arpack_jll, SHA
run(pipeline(`nm $(Arpack_jll.libarpack_path)`,`grep dgetv0`))
run(pipeline(`nm $(Arpack_jll.libarpack_path)`,`grep iseed`))
@show bytes2hex(open(SHA.sha256, Arpack_jll.libarpack_path))
iseed_offsets_str = readchomp(pipeline(`nm $(Arpack_jll.libarpack_path)`,`grep iseed`))
iseed_offsets = iseed_offsets_str |> 
          x->split(x,"\n") |>
          x->map(y->y[1:16],x) |>
          x->parse.(UInt64, x;base=16)
dgetv0_offset = readchomp(pipeline(`nm $(Arpack_jll.libarpack_path)`,`grep dgetv0`)) |> x->parse(UInt64, x[1:16];base=16)


## This code will update iseed, and allow you to check if you get it right/wrong...
import Arpack_jll, LinearAlgebra
function arpack_dgetv0!(ido::Ref{LinearAlgebra.BlasInt}, bmat::Symbol, itry::Int, initv::Bool,
  n::Int, j::Int, v::StridedMatrix{Float64}, ldv::Int,
    resid::StridedVecOrMat{Float64}, rnorm::Ref{Float64},
    ipntr::StridedVecOrMat{LinearAlgebra.BlasInt}, workd::StridedVecOrMat{Float64})
  ierr = Ref{LinearAlgebra.BlasInt}(0)
  # NOTE, arpack doesn't touch ierr unless ido[] == 0 or there is
  # a restart failure.
  # 
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
function check_iseed_offset(dgetv0offset, iseed_offsets)

  libar = Base.Libc.dlopen(Arpack_jll.libarpack)
  dgetv0_real_offset = Base.Libc.dlsym(libar, "dgetv0_")
  dgetv0_nm_offset = dgetv0offset # this comes from the command above
  base_offset = dgetv0_real_offset - dgetv0_nm_offset



  iseeds = Ptr{LinearAlgebra.BlasInt}.(base_offset .+ iseed_offsets)
  iseedvals = [unsafe_load(iseeds[i], j) for i=1:length(iseeds), j=1:4]
  @show iseedvals

  n = 10

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
  #println.(resid)

  iseedvals = [unsafe_load(iseeds[i], j) for i=1:length(iseeds), j=1:4]
  @show iseedvals
  @show iseed_offsets
  @show dgetv0_nm_offset
  return resid
  
end
resid1 = check_iseed_offset(dgetv0_offset, iseed_offsets)
resid2 = check_iseed_offset(dgetv0_offset, iseed_offsets)

include("../arpackjll.jl")

_reset_libarpack_dgetv0_iseed()

resid3 = check_iseed_offset(dgetv0_offset, iseed_offsets)
@show resid1 == resid3
@show resid1 == resid2


exit()
##
file = @__FILE__
run(`/Applications/Julia-1.7.app/Contents/Resources/julia/bin/julia $file`)