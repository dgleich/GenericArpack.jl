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

  @testset "qr apply with dorm2r" begin 
    function compare_qr_apply(A::StridedMatrix)
      mn = maximum(size(A))
      QH1, tau1 = _simple_dgeqr2(copy(A))
      QH2, tau2 = _dgeqr2_blas!(copy(A))

      Q1b = _dorm2r_blas!('L', 'N', mn, mn, length(tau1), QH1, tau1, Matrix(1.0I, mn, mn))
      Q2b = _dorm2r_blas!('R', 'N', mn, mn, length(tau1), QH1, tau1, Matrix(1.0I, mn, mn))
      Q3b = _dorm2r_blas!('L', 'T', mn, mn, length(tau1), QH1, tau1, Matrix(1.0I, mn, mn))
      Q4b = _dorm2r_blas!('R', 'T', mn, mn, length(tau1), QH1, tau1, Matrix(1.0I, mn, mn))

      work = zeros(mn)
      Q1j = ArpackInJulia.dorm2r(Val(:L), Val(:N), mn, mn, length(tau1), QH1, tau1, Matrix(1.0I, mn, mn), work)
      Q2j = ArpackInJulia.dorm2r(Val(:R), Val(:N), mn, mn, length(tau1), QH1, tau1, Matrix(1.0I, mn, mn), work)
      Q3j = ArpackInJulia.dorm2r(Val(:L), Val(:T), mn, mn, length(tau1), QH1, tau1, Matrix(1.0I, mn, mn), work)
      Q4j = ArpackInJulia.dorm2r(Val(:R), Val(:T), mn, mn, length(tau1), QH1, tau1, Matrix(1.0I, mn, mn), work)

      @test Q1b == Q1j
      @test Q2b == Q2j
      @test Q3b == Q3j
      @test Q4b == Q4j
    end

    n = 7
    compare_qr_apply(Matrix(SymTridiagonal(zeros(n), ones(n-1))))
  
  end 
end 
