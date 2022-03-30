using LinearAlgebra
using DoubleFloats

# specialize for Double64
#floatmin2(::Type{T}) where {T} = (twopar = 2one(T); twopar^trunc(Integer,log(floatmin(T)/eps(T))/log(twopar)/twopar))
import LinearAlgebra.floatmin2
LinearAlgebra.floatmin2(::Type{Double64}) = Double64(reinterpret(Float64, 0x2350000000000000), 0.0)

using Random

function ArpackInJulia._dlarnv_idist_2!(iseed::Base.RefValue{NTuple{4,Int}}, n::Int, x::Vector{Double64})
  # just a hacky Implementation

  rng = Random.MersenneTwister(UInt32[iseed[][1],iseed[][2],iseed[][3],iseed[][4]])
  rand!(rng, @view(x[1:n]))
  iseed[] = tuple(rand(rng, UInt32),rand(rng, UInt32),rand(rng, UInt32),rand(rng, UInt32))
end 

function ArpackInJulia._dnrm2_unroll_ext(a::AbstractVector{T}) where {T <: Double64}
  norm(a,2)
end 

n = 25 
T = Double64
op = ArpackInJulia.ArpackSimpleOp(Diagonal(one(T):n))

ido = Ref{Int}(0)
bmat = :I

which = :LM
tol = zero(T) # just use the default
resid = zeros(T, n)
ncv = min(2nev, n-1)
V = zeros(T, n,ncv)
ldv = n
mode = 1 
iparam = zeros(Int,11)
iparam[1] = 1
iparam[3] = 300 # 
iparam[4] = 1
iparam[7] = mode 
ipntr = zeros(Int,11)
workd = zeros(T,3n)
lworkl = ncv*ncv + 8*ncv
workl = zeros(T,lworkl)

iparam0 = copy(iparam)
info_initv = 0

state = ArpackInJulia.ArpackState{T}()
stats = ArpackStats()
debug = ArpackInJulia.ArpackDebug()
ArpackInJulia.set_debug_high!(debug)

iter = 1
while ido[] != 99
  if ido[] == 1 || ido[] == -1
    ArpackInJulia._i_do_now_opx_1!(op, ipntr, workd, n)
  end 
  local ierr = ArpackInJulia.dsaupd!(ido, Val(bmat), n, which, nev, tol, resid, ncv, V, ldv, iparam,
    ipntr, workd, workl, lworkl, info_initv;
    state,  debug 
  ).ierr 

  println("Call number $iter, ierr=$ierr, ido=$(ido[])")

  global iter = iter + 1
end 
