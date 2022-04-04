using ArpackInJulia
using LinearAlgebra

##
@testset "simple tridiag" begin 
  n = 12
  A = SymTridiagonal(collect(range(1.0, 0.1, n)), ones(n-1))
  Z = zeros(n,n)
  lams = ArpackInJulia.flexible_dsteqr!(copy(A), Z)[1]
  @test 
  display(diag(lams)')
  display(Z)
  display(eigen(A).values')
  display(eigen(A).vectors)
end 

##
@testset "reverse simple tridiag" begin 
  n = 12
  A = SymTridiagonal(collect(reverse(range(1.0, 0.1, n))), ones(n-1))
  Z = zeros(n,n)
  lams = ArpackInJulia.flexible_dsteqr!(copy(A), Z)[1]
  display(diag(lams)')
  display(Z)
  display(eigen(A).values')
  display(eigen(A).vectors)
end 
