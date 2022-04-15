# https://discourse.julialang.org/t/memory-allocation-and-profile/47573/5
using GenericArpack
using LinearAlgebra

#=
We used this to track down where copyto! was making copies with 
what were actually  non-aliased arrays Julia couldn't detect.

function Base.mightalias(A::AbstractArray, B::AbstractArray) 
  if !Base._isdisjoint(Base.dataids(A), Base.dataids(B))

    println(Base.dataids(A))
    println(Base.dataids(B))
    error("Non disjoint copyto")
  end 

  !isbits(A) && !isbits(B) && !Base._isdisjoint(Base.dataids(A), Base.dataids(B))
end 
=#

function eigrun(op,ido, ::Val{BMAT}, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state) where BMAT 
  niter = 0 
  nbytes = 0 
  while ido[] != 99
    nbytes += @allocated GenericArpack.dsaupd!(ido, Val(BMAT), n, which, nev, tol, resid, ncv, V, ldv, iparam,
      ipntr, workd, workl, lworkl, info_initv;
      state 
    )
    if ido[] == 1 || ido[] == -1
      niter += 1
      GenericArpack._i_do_now_opx_1!(op, ipntr, workd, n)
    elseif ido[] == 99
      break
    else
      @error("this only supports standard eigenvalue problems")
    end 
  end
  return niter, nbytes
end 

op = GenericArpack.ArpackSimpleOp(Diagonal(1.0:10^3))
nev = 6
ido = Ref{Int}(0)
bmat = :I
n = size(op.A,1)
which = :LM
tol = 0.0 # just use the default
resid = zeros(n)
ncv = min(2nev, n-1)
V = zeros(n,ncv)
ldv = n
mode = 1 
iparam = zeros(Int,11)
iparam[1] = 1
#iparam[3] = 300 # max iteration
iparam[3] = 300 # 
iparam[4] = 1
iparam[7] = mode 
ipntr = zeros(Int,11)
workd = zeros(3n)
lworkl = ncv*ncv + 8*ncv
workl = zeros(lworkl)

info_initv = 0

# Note that we cannot run two sequences at once and check them where we start a whole
# second arpack call because of the expected Arpack state. 
state = GenericArpack.ArpackState{Float64}()


eigrun(op, ido, Val(bmat), n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state);

# reset state
state = GenericArpack.ArpackState{Float64}()
ido[] = 0 

valbmat = Val(bmat)

# rest profiling after compile
using Profile
Profile.clear_malloc_data()

niter = eigrun(op, ido, valbmat, n, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, info_initv, state);
println(niter)

exit()
##
include("allocations.jl")
#lines = report_allocations(@__FILE__; system=false)
lines = report_allocations(@__FILE__)
println.(lines);
