using LinearAlgebra
include("../src/arpack-blas-direct.jl")

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

  test_diffs("dsteqr: d1,d2", d1, d2)
  test_diffs("dsteqr: e1,e2", e1, e2)
  test_diffs("dsteqr: Z1,Z2", Z1, Z2)
  test_diffs("dsteqr: work1,work2", work1, work2)
  @test d1 == d2
  @test e1 == e2
  @test Z1 == Z2
  @test work1 == work2 
end 

@testset "dsteqr force qr" begin
  n = 12
  compare_dsteqr(SymTridiagonal(collect(range(1.0, 0.1, n)), ones(n-1)))
end

@testset "dsteqr forced ql" begin
  n = 15
  compare_dsteqr(SymTridiagonal(collect(range(1.0, 2.0, n)), ones(n-1)))
end

@testset "dsteqr split" begin
  Hdata = [7.078500017390861 29.999999976645363; 0.00015282224133827933 28.999999984019453; 0.00019833099664768876 28.00000003918011; 1.2453099363261543e-5 26.999999999305622; 2.9144607362751737e-5 26.000000000849127; 5.476034823168761e-7 25.00000000000023; 4.9362983911508065e-15 12.856546991821302; 8.536389676058763 12.201877285145336; 5.696747114072514 12.741129075106617; 6.219630101820137 11.83319317846694; 5.4770323243957435 13.437916643126739; 5.271530633768442 13.161663328320165]
  H = SymTridiagonal(Hdata[:,2], Hdata[2:end, 1])
  compare_dsteqr(H)
end


