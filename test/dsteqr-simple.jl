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

@testset "issue #5 from GenericArpack.jl" begin
  Adiag = [0.08599443783810437, 0.9140055621618955, 0.9257211782861121, 0.07427882171388793, 0.9927790910828268, 0.007220908917173123, 0.9565054060449746, 0.043494593955025376]
  Asubdiag = [0.28035583550019566, 9.277803709232427e-17, 0.26222409950018777, 8.063486413292979e-17, 0.08466857381332825, 1.7712507888072247e-18, 0.20396767942915087]
  A = SymTridiagonal(Adiag, Asubdiag)
  n = length(Adiag)
  Z = zeros(n,n)
  lams = GenericArpack.flexible_dsteqr!(copy(A), Z)[1]
  @test(A*Z ≈ Z*Diagonal(lams))
  @test(Z'*Z ≈ Diagonal(ones(n)))
  E = eigen(A)
  @test E.values ≈ diag(lams)
end
