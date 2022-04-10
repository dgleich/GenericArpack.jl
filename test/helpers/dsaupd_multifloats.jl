using LinearAlgebra
using MultiFloats
using Random
using ArpackInJulia

ArpackInJulia._eps23(Float64x4) = Float64x4((1.8078980427655233e-42, 5.038360808408627e-59, -4.0108760346726413e-75, 1.2660248040679727e-91))
ArpackInJulia._eps23(::Type{Float64x4}) = Float64x4((1.8078980427655233e-42, 5.038360808408627e-59, -4.0108760346726413e-75, 1.2660248040679727e-91))
ArpackInJulia._dstqrb_maxit(::Type{Float64x4}) = 150

n = 25 
T = Float64x4
#op = ArpackInJulia.ArpackSimpleOp(Diagonal(1.0:n))
op = ArpackInJulia.ArpackSimpleOp(Float64x4.(Diagonal(1.0:n)))

ido = Ref{Int}(0)
bmat = :I

which = :LM
nev = 6
tol = eps(T)/10 # just use the default
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
