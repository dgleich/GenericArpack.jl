@testset "blas qr compare" begin 
  include("../src/arpack-blas-direct.jl")
  using LinearAlgebra

  function _simple_dgeqr2(A)
    tau = zeros(minimum(size(A)))
    work = zeros(size(A,2))
    cA = copy(A)
    ArpackInJulia._dgeqr2!(cA, tau, work)
    return cA, tau
  end 

  function compare_qr(A)

    B, tauB = _dgeqr2_blas!(copy(A))
    tau = zeros(minimum(size(A)))
    work = zeros(size(A,2))
    ArpackInJulia._dgeqr2!(A, tau, work)

    @test A == B
    @test tau == tauB 
  end 

  compare_qr(zeros(10,8))
  compare_qr(zeros(8,10))
  compare_qr(Matrix(1.0I, 5,8))
  compare_qr(Matrix(1.0I, 8,8))
  compare_qr(Matrix(1.0I, 5,5))
  
  n = 10
  compare_qr(Matrix(SymTridiagonal(zeros(n), ones(n-1))))
  compare_qr(Matrix(Tridiagonal(-ones(n-1), zeros(n), ones(n-1))))

  n = 21
  compare_qr(Matrix(Tridiagonal(-ones(n-1), 2*ones(n), ones(n-1))))

  compare_qr(Matrix(reshape(1:10.0, 10, 1)))
  compare_qr(Matrix(reshape(1:10.0, 1, 10)))
end 
