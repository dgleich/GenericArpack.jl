## Test accuracy...
using DoubleFloats
using MultiFloats
using GenericArpack
using LinearAlgebra

function arpack_dsdrv3(::Type{T}, n::Int) where T 
  A = Tridiagonal(-ones(T, n-1),2*ones(T,n),-ones(T,n-1)).*(n+1)
  B = SymTridiagonal(4*ones(T, n),ones(T,n-1)).*(one(T)/(6*(n+1)))
  return A, B
end 

# compute 
T = Float64x4
n = 100 
setprecision(BigFloat, 2048) # very high precision... 
lams = 101*(2 .-2*cos.(pi.*BigFloat.(collect(1:n))./(n+1)))
A = SymTridiagonal(arpack_dsdrv3(T, n)[1])
##
l1 = eigvals(SymTridiagonal(A))
##
l2 = eigs(Symmetric(A), 4)

##
vals = Any[l2.values[end], lams[end], 0,  l1[end],lams[end], 0, l2.values[end], l1[end]]
##
l2.values[end] - lams[end]
##
l1[end] - lams[end]
##

##
ind = n-1
l2ind = ind - n + length(l2.values)
vals = Any[l2.values[l2ind], lams[ind], 0,  l1[ind],lams[ind], 0, l2.values[l2ind], l1[ind]]
##
l2.values[l2ind] - lams[ind]
##
l1[ind] - lams[ind]
##

##
T = BigFloat
n = 100 
setprecision(BigFloat, 2048) # very high precision... 
lams = 101*(2 .-2*cos.(pi.*BigFloat.(collect(1:n))./(n+1)))
A = SymTridiagonal(arpack_dsdrv3(T, n)[1])
##
@time l2 = eigs(Symmetric(A), 4; ncv=40, maxiter=1000)