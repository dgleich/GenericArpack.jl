using GenericArpack
using LinearAlgebra
using Test
include("../../src/arpack-blas-direct.jl")

##
function test_diffs(tagstr, a, b)
  @assert(size(a) == size(b)) # otherwise, use view to make them the same...
  first = false 
  for i in eachindex(a)
    if a[i] != b[i]
      if first == false
        first = true
        println("Difference between $tagstr")
      end 
      println(" ", rpad(a[i], 30), " , ", rpad(b[i], 30), " # ", i)
    end 
  end
end 
function compare_dsteqr(A::SymTridiagonal{Float64})
  d1 = A.dv
  e1 = A.ev
  Z1 = zeros(length(d1),length(d1))
  work1 = zeros(2*length(d1))

  d2 = copy(d1)
  e2 = copy(e1)
  Z2 = copy(Z1)
  work2 = copy(work1)

  GenericArpack.simple_dsteqr!(d1, e1, Z1; work=work1)
  _dsteqr_blas!(d2, e2, Z2, work2)
  @show d2

  test_diffs("dsteqr: d1,d2", d1, d2)
  test_diffs("dsteqr: e1,e2", e1, e2)
  test_diffs("dsteqr: Z1,Z2", Z1, Z2)
  test_diffs("dsteqr: work1,work2", work1, work2)
  @test d1 == d2
  @test e1 == e2
  @test Z1 == Z2
  @test work1 == work2 
end 

compare_dsteqr(SymTridiagonal(collect(range(1.0, 0.1, n)), ones(n-1)))

##
@testset "simple tridiag" begin 
  n = 12
  A = SymTridiagonal(collect(range(1.0, 0.1, n)), ones(n-1))
  Z = zeros(n,n)
  lams = GenericArpack.flexible_dsteqr!(copy(A), Z)[1]
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
  lams = GenericArpack.flexible_dsteqr!(copy(A), Z)[1]
  display(diag(lams)')
  display(Z)
  display(eigen(A).values')
  display(eigen(A).vectors)
end 
