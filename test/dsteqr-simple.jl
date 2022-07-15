using GenericArpack
using LinearAlgebra
using Test

##
@testset "simple tridiag" begin 
  n = 12
  dvals = collect(range(1.0, 0.1; length=n)) # collect(range(1.0, 0.1, n)
  A = SymTridiagonal(dvals, ones(n-1))
  Z = zeros(n,n)
  lams = GenericArpack.flexible_dsteqr!(copy(A), Z)[1]
  E = eigen(A)
  @test abs.(E.vectors) ≈ abs.(Z)
  @test E.values ≈ diag(lams)
end 

##
@testset "reverse simple tridiag" begin 
  n = 12
  dvals = collect(reverse(range(1.0, 0.1; length=n)))
  A = SymTridiagonal(dvals, ones(n-1))
  Z = zeros(n,n)
  lams = GenericArpack.flexible_dsteqr!(copy(A), Z)[1]
  E = eigen(A)
  @test abs.(E.vectors) ≈ abs.(Z)
  @test E.values ≈ diag(lams)
end 
