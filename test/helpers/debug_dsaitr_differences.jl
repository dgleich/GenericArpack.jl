## Issue: debugging
# while developing dsaitr! I found a case where we
# get answers very close to Julia, but are off by
# a few floating point operations. I wanted to 
# understand this case in depth, so I could figure 
# out why the final answer seemed to be off by
# even more. 

## The code here was just to run the darn thing
# twice with debug mode on to see what was happening.

# This was causing it to get 
using Revise
using GenericArpack
using Test
using LinearAlgebra
include("../arpackjll.jl")
include("../utility.jl")
arpack_set_debug_high()

function store_arpackjll_dsaitr_sequence(M;
  B=1.0LinearAlgebra.I,
  idostart::Int,
  bmat::Symbol,
  n::Int,
  k::Int, 
  np::Int,
  mode::Int,
  resid::StridedVecOrMat{T},
  rnorm::Ref{T},
  V::StridedMatrix{T},
  ldv::Int,
  H::StridedMatrix{T},
  ldh::Int, 
) where T

  ido = Ref{Int}(idostart)
  ipntr = zeros(Int, 3)
  workd = zeros(10, 3)
  rnorm = Ref{T}(rnorm[])
  ipntr = copy(ipntr)
  
  histdata = Vector{
      NamedTuple{(:V,:H,:resid,:workd), Tuple{Matrix{T},Matrix{T},Vector{T},Matrix{T}}}
  }()


  if bmat==:G
    mul!(@view(workd[n+1:2n]),M,@view(resid[1:n]))
  end 

  state = GenericArpack.ArpackState{Float64}()
  while ido[] != 99
    info = arpack_dsaitr!(
      ido, bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh, 
      ipntr, workd)
    @info "Arpack_jll: Step $(length(histdata)): ", ido[], info, mode, bmat  
    if ido[] == -1 || ido[] == 1
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    elseif ido[] == 2
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    end
    # everything should be same...
    push!(histdata, (;V=copy(V),H=copy(H),resid=copy(resid),workd=copy(workd)))
  end
  return histdata
end

function store_GenericArpack_dsaitr_sequence(M;
  B=1.0LinearAlgebra.I,
  idostart::Int,
  bmat::Symbol,
  n::Int,
  k::Int, 
  np::Int,
  mode::Int,
  resid::StridedVecOrMat{T},
  rnorm::Ref{T},
  V::StridedMatrix{T},
  ldv::Int,
  H::StridedMatrix{T},
  ldh::Int, 
) where T

  ido = Ref{Int}(idostart)
  ipntr = zeros(Int, 3)
  workd = zeros(10, 3)
  rnorm = Ref{T}(rnorm[])
  ipntr = copy(ipntr)

  histdata = Vector{
      NamedTuple{(:V,:H,:resid,:workd), Tuple{Matrix{T},Matrix{T},Vector{T},Matrix{T}}}
  }()

  state = GenericArpack.ArpackState{Float64}()
  debug = GenericArpack.set_debug_high!(GenericArpack.ArpackDebug())
  #debug = GenericArpack.ArpackDebug()
  while ido[] != 99
    info = GenericArpack.dsaitr!(
      ido, Val(bmat), n, k, np, mode, resid, rnorm, V, ldv, H, ldh, 
      ipntr, workd;
      state, debug)
    @info "GenericArpack: Step $(length(histdata)): ", ido[], info, mode, bmat
    if ido[] == -1 || ido[] == 1
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    elseif ido[] == 2
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    end
    # everything should be same...
    push!(histdata, (;V=copy(V),H=copy(H),resid=copy(resid),workd=copy(workd)))
  end
  return histdata
end


function run_example_jll()
  bmat = :I
  n = 10
  k = 0 # number of current columns in V
  np = 3
  #V = randn(n,k+np) # total memory for v
  #V[:,1:k] = qr(V).Q[:,1:k] # random orthogonal init
  mode = 1
  resid = collect(1.0:n)
  rnorm = Ref{Float64}(norm(resid))

  V = zeros(n,k+np)
  ldv = n 
  H = zeros(n,2) # full h
  ldh = n 

  M = Diagonal(1.0:n)
  
  histinfo = store_arpackjll_dsaitr_sequence(M;  idostart=0, 
    bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
  )
  return histinfo
end

function run_example_AiJ()
  bmat = :I
  n = 10
  k = 0 # number of current columns in V
  np = 3
  #V = randn(n,k+np) # total memory for v
  #V[:,1:k] = qr(V).Q[:,1:k] # random orthogonal init
  mode = 1
  resid = collect(1.0:n)
  rnorm = Ref{Float64}(norm(resid))

  V = zeros(n,k+np)
  ldv = n 
  H = zeros(n,2) # full h
  ldh = n 

  M = Diagonal(1.0:n)
  
  histinfo = store_GenericArpack_dsaitr_sequence(M;  idostart=0, 
    bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
  )
  return histinfo
end

map(x->println(), 1:5)
println("Run Example - Try 1")
hist1 = run_example_jll()

map(x->println(), 1:5)
println("Run Example - Try 2")
hist2 = run_example_AiJ()

##
@show hist1[1] == hist2[1]
## That was true! so one iteration is the same ;) 

##
@show hist1[2] == hist2[2]

## That was false! :( so two iterations are not the same.
# Where are the differences?

##
A = hist1[2]
B = hist2[2]
@show A.resid == B.resid
@show A.H == B.H
##
floatsbetween.(A.V, B.V) 
##
floatsbetween.(A.workd, B.workd)
##
@test A.workd ≈ B.workd atol=0

##
@show hist1[3] == hist2[3]
##
@show hist1[4] == hist2[4]


## Now try a generalized example
function run_example_jll()
  bmat = :I
  n = 10
  k = 0 # number of current columns in V
  np = 3
  #V = randn(n,k+np) # total memory for v
  #V[:,1:k] = qr(V).Q[:,1:k] # random orthogonal init
  mode = 1
  resid = collect(1.0:n)
  rnorm = Ref{Float64}(norm(resid))

  V = zeros(n,k+np)
  ldv = n 
  H = zeros(n,2) # full h
  ldh = n 

  M = Diagonal(1.0:n)
  
  histinfo = store_arpackjll_dsaitr_sequence(M;  idostart=0, 
    bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
  )
  return histinfo
end

##
function run_example_generalized()
  bmat = :G # this isn't valid, but seems to indicate something else wrong...
  # as it still _should_ be the same :) 
  n = 10
  k = 0 # number of current columns in V
  np = 3
  #V = randn(n,k+np) # total memory for v
  #V[:,1:k] = qr(V).Q[:,1:k] # random orthogonal init
  mode = 1
  resid = ones(n)/sqrt(n)
  rnorm = Ref{Float64}(norm(resid))

  V = zeros(n,k+np)
  ldv = n 
  H = zeros(n,2) # full h
  ldh = n 

  M = Diagonal(1.0:n)
  B = Diagonal(range(0.1, 1.0, length=n))
  
  histinfo2 = store_GenericArpack_dsaitr_sequence(M; B, idostart=0, 
    bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
  )

  resid = ones(n)/sqrt(n)
  rnorm = Ref{Float64}(norm(resid))
  V = zeros(n,k+np)
  ldv = n 
  H = zeros(n,2) # full h
  ldh = n 

  histinfo1 = store_arpackjll_dsaitr_sequence(M; B, idostart=0, 
    bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
  )
  return histinfo1, histinfo2
end
hist1,hist2 = run_example_generalized()

##
hist1[1] == hist2[1]
##
hist1[2] == hist2[2]
##
hist1[3] == hist2[3]
##
hist1[4] == hist2[4]

##
jll = hist1[2]
aij = hist2[2]
@show jll.resid == aij.resid
@show jll.H == aij.H
@show jll.V == aij.V

## Study the matrix V from a generalized example
# Ahh, that's right, mode 2 is the worst!
function study_generalized_example()
  Random.seed!(0)
  ido = Ref{Int}(0)
  bmat = :G
  n = 10
  k = 0 # number of current columns in V
  np = 3
  mode = 2
  #resid = randn(n)
  resid = zeros(n)

  A = Diagonal(1.0:n)
  B = Diagonal(collect(range(0.1, 1.0, length=n)))
  M = inv(sqrt.(B))*A # okay since they are diagonal....
  @show

  V = zeros(n,k+np)
  ldv = n 
  H = zeros(n,2) # full h
  ldh = n 

  rnorm = Ref{Float64}(sqrt(abs(resid'*B*resid)))

  #return store_GenericArpack_dsaitr_sequence(M; B, idostart=0, 
  return store_arpackjll_dsaitr_sequence(M; B, idostart=0,
    bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
  )
end
hist = study_generalized_example() 