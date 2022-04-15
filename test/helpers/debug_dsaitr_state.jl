## Issue: debugging
# while developing dsaitr! I found a case where it seemed like
# on run 1, I would get the same floating point values as Julia
# but on run 2, I would not. This is annoying because... well, 
# we should be more determinstic and this case (with a diagonal
# matrix with distinct entries initialized with the vector of
# all ones, should produce a complete orthogonal factorization)

## Okay, we figured this one out
# This happened because I had a typo in
# the idonow handler that had a == 2 instead of == 1
# so we weren't doing the matvec. So the alg
# thought it was getting the identity matrix, and then
# callding dgetv0... 
# which used the random state... 
# hence, non-determistic. 

## The code here was just to run the darn thing
# twice with debug mode on to see what was happening.

# This was causing it to get 
using Revise
using GenericArpack
using Test
using LinearAlgebra
include("../arpackjll.jl")
arpack_set_debug_high()

function store_arpackjll_dsaitr_sequence(M;
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
  workd = zeros(3n)
  rnorm = Ref{T}(rnorm[])
  ipntr = copy(ipntr)
  workd = copy(workd)

  histdata = Vector{
      NamedTuple{(:V,:H,:resid,:workd), Tuple{Matrix{T},Matrix{T},Vector{T},Vector{T}}}
  }()

  state = GenericArpack.ArpackState{Float64}()
  while ido[] != 99
    arinfo = arpack_dsaitr!(
      ido, bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh, 
      ipntr, workd)
      
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

function run_example()
  bmat = :I
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
  
  histinfo = store_arpackjll_dsaitr_sequence(M;  idostart=0, 
    bmat, n, k, np, mode, resid, rnorm, V, ldv, H, ldh
  )
  return histinfo
end

map(x->println(), 1:5)
println("Run Example - Try 1")
hist1 = run_example()

map(x->println(), 1:5)
println("Run Example - Try 2")
hist2 = run_example()

##

