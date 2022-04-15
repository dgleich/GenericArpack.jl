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
      mn  = maximum(size(A))
      m,n = size(A)
      QH1, tau1 = _simple_dgeqr2(copy(A))
      QH2, tau2 = _dgeqr2_blas!(copy(A))

      @test QH1 == QH2
      @test tau1 == tau2
      Q1b = _dorm2r_blas!('L', 'N', m, n, length(tau1), QH1, tau1, Matrix(1.0I, m, n))
      Q2b = _dorm2r_blas!('R', 'N', m, m, length(tau1), QH1, tau1, Matrix(1.0I, m, m))
      Q3b = _dorm2r_blas!('L', 'T', n, n, length(tau1), QH1, tau1, Matrix(1.0I, n, n))
      Q4b = _dorm2r_blas!('R', 'T', m, m, length(tau1), QH1, tau1, Matrix(1.0I, m, m))



      work = zeros(mn)
      Q1j = ArpackInJulia.dorm2r!(Val(:L), Val(:N), m, n, length(tau1), QH1, tau1, Matrix(1.0I, m, n), work)
      Q2j = ArpackInJulia.dorm2r!(Val(:R), Val(:N), m, m, length(tau1), QH1, tau1, Matrix(1.0I, m, m), work)
      Q3j = ArpackInJulia.dorm2r!(Val(:L), Val(:T), n, n, length(tau1), QH1, tau1, Matrix(1.0I, n, n), work)
      Q4j = ArpackInJulia.dorm2r!(Val(:R), Val(:T), m, m, length(tau1), QH1, tau1, Matrix(1.0I, m, m), work)

      @test Q1b == Q1j
      @test Q2b == Q2j
      @test Q3b == Q3j
      @test Q4b == Q4j
    end

    n = 7
    compare_qr_apply(Matrix(SymTridiagonal(zeros(n), ones(n-1))))

    Q1 = [-5.019931333600887e-26 4.583545848260025e-19 -6.290763553850451e-14 1.5154694014158623e-8 0.00015282224011718752 -0.9999999883226816; 1.6424086134953873e-21 -1.1997064776330149e-14 1.234917793991677e-9 -0.00019833099821702998 -0.9999999686550888 -0.00015282224011718755; -3.312459742101567e-17 1.8147034339346287e-10 -1.2453099114015112e-5 0.9999999802548672 -0.0001983309982170299 -1.5154693893460058e-8; 7.979844240888431e-12 -2.9144607372867256e-5 0.999999999497756 1.245309911401511e-5 -1.2349178135518348e-9 -6.290763628158253e-14; -5.476034820843436e-7 0.999999999575146 2.9144607372867262e-5 1.8147034214075507e-10 -1.1997064937141322e-14 -4.5835458998368585e-19; 0.99999999999985 5.476034820843435e-7 7.979844244071235e-12 3.3124597097839987e-17 -1.6424086342901649e-21 -5.019931392220031e-26; 3.760421589892889e-16 1.919443790093301e-22 2.61862134362123e-27 1.0215964062303585e-32 -4.777060177258836e-37 -1.3812568763549187e-41; 1.5061066509517348e-16 7.343137757532626e-23 9.588246944325921e-28 3.586789342676401e-33 -1.6109517000400644e-37 -4.480947425847092e-42; 0.0 0.0 0.0 0.0 0.0 -0.0; 0.0 0.0 0.0 0.0 0.0 -0.0; 0.0 0.0 0.0 0.0 0.0 -0.0; 0.0 0.0 0.0 0.0 0.0 -0.0]
    compare_qr_apply(Q1)
    Q2 = [-5.0199313336008365e-26 0.0 -6.290763553850389e-14 1.5154694014158603e-8 -0.00015282224011718812 0.9999999883226813; 1.6424086134953723e-21 -1.1997064776330143e-14 1.2349177939916661e-9 -0.00019833099821702952 0.9999999686550886 0.00015282224011718804; -3.3124597421015475e-17 1.8147034339346326e-10 -1.2453099114015055e-5 0.9999999802548674 0.0001983309982170295 1.5154693893460117e-8; 7.979844240888405e-12 -2.9144607372867357e-5 0.9999999994977561 1.245309911401506e-5 1.2349178135518293e-9 6.290763628158279e-14; -5.476034820843422e-7 0.9999999995751461 2.9144607372867377e-5 1.8147034214075442e-10 1.1997064937141275e-14 4.583545899836887e-19; 0.9999999999998501 5.476034820843429e-7 7.979844244071283e-12 0.0 0.0 5.019931392220078e-26; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
    compare_qr_apply(Q2)
    Q3 = reinterpret(Float64, UInt64[0xbaaf126611decddc 0x3c20e909625f8160 0xbd31b4f90d9567aa 0x3e5045b0c942a3fb 0x3f2407dd0bfbb616 0xbfeffffff9bb14e2; 0x3b9f0633622dea91 0xbd0b03d644bd8ecb 0x3e153739d4475678 0xbf29fee24d604783 0xbfefffffef2bfcfb 0xbf2407dd0bfbb617; 0xbc831853e8b15e7d 0x3de8f0eb8ab1844a 0xbeea1db4e8e2f8fe 0x3feffffff5663fe7 0xbf29fee24d604780 0xbe5045b0c7160444; 0x3da18c40c0eb6240 0xbefe8f722b184dec 0x3fefffffffbaf8dc 0x3eea1db4e8e2f8fd 0xbe153739d9ea9f67 0xbd31b4f91117bb2c; 0xbea25fe0cabb9429 0x3fefffffffc59bc7 0x3efe8f722b184dee 0x3de8f0eb87ce0c70 0xbd0b03d64ad0d3f9 0xbc20e9096590c5de; 0x3feffffffffffab9 0x3ea25fe0cabb9428 0x3da18c40c10971cc 0x3c831853e5913990 0xbb9f063368c4fd33 0xbaaf126617f52b15; 0x3cbb18c0e3440a36 0x3b6d017d0e85f204 0x3a69eefe6c09a3d9 0x394a85ac373b6398 0xb86451c205cba689 0xb773407d3dd19a98; 0x3ca5b48d837c5bb0 0x3b56317a476f9016 0x3a52fdd19cb8f127 0x39329fa92acbd4e1 0xb84b68af2b637a50 0xb758fb6bfffa949b; 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x8000000000000000; 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x8000000000000000; 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x8000000000000000; 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x8000000000000000])
    compare_qr_apply(Q3)


  end 

  @testset "qr with views" begin 
    Qdata = reinterpret(Float64, UInt64[0xbaaf126611decddc 0x3c20e909625f8160 0xbd31b4f90d9567aa 0x3e5045b0c942a3fb 0x3f2407dd0bfbb616 0xbfeffffff9bb14e2; 0x3b9f0633622dea91 0xbd0b03d644bd8ecb 0x3e153739d4475678 0xbf29fee24d604783 0xbfefffffef2bfcfb 0xbf2407dd0bfbb617; 0xbc831853e8b15e7d 0x3de8f0eb8ab1844a 0xbeea1db4e8e2f8fe 0x3feffffff5663fe7 0xbf29fee24d604780 0xbe5045b0c7160444; 0x3da18c40c0eb6240 0xbefe8f722b184dec 0x3fefffffffbaf8dc 0x3eea1db4e8e2f8fd 0xbe153739d9ea9f67 0xbd31b4f91117bb2c; 0xbea25fe0cabb9429 0x3fefffffffc59bc7 0x3efe8f722b184dee 0x3de8f0eb87ce0c70 0xbd0b03d64ad0d3f9 0xbc20e9096590c5de; 0x3feffffffffffab9 0x3ea25fe0cabb9428 0x3da18c40c10971cc 0x3c831853e5913990 0xbb9f063368c4fd33 0xbaaf126617f52b15; 0x3cbb18c0e3440a36 0x3b6d017d0e85f204 0x3a69eefe6c09a3d9 0x394a85ac373b6398 0xb86451c205cba689 0xb773407d3dd19a98; 0x3ca5b48d837c5bb0 0x3b56317a476f9016 0x3a52fdd19cb8f127 0x39329fa92acbd4e1 0xb84b68af2b637a50 0xb758fb6bfffa949b; 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x8000000000000000; 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x8000000000000000; 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x8000000000000000; 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x0000000000000000 0x8000000000000000])
    m, n = size(Qdata)

    qdatavec = copy(vec(Qdata))

    worklen = n # number of columns
    taulen = n # number of cols

    qdatalen = length(qdatavec)
    append!(qdatavec, zeros(2n))
    tau = @view(qdatavec[qdatalen+1:qdatalen+n])
    work = @view(qdatavec[qdatalen+n+1:qdatalen+2n])
    Q = reshape(@view(qdatavec[1:qdatalen]), m, n)

    Q1 = copy(Q)
    tau1 = copy(tau)
    work1 = copy(work)

    Q2 = copy(Q)
    tau2 = copy(tau)
    work2 = copy(work)

    ArpackInJulia._dgeqr2!(Q, tau, work)
    _dgeqr2_blas!(Q1, tau1, work1)

    test_diffs("Q, Q1", Q, Q1)
    test_diffs("tau, tau1", tau, tau1)
    test_diffs("work, work1", work, work1)
    @test Q == Q1
    @test tau == tau1
    @test work == work1 

    ArpackInJulia._dgeqr2!(Q2, tau2, work2)

    test_diffs("Q1, Q2", Q1, Q2)
    test_diffs("tau1, tau2", tau1, tau2)
    test_diffs("work1, work2", work1, work2)
    @test Q1 == Q2
    @test tau1 == tau2
    @test work1 == work2

    # now test applying this...
    longn = 100
    ncv = m
    nconv = n 
    Vfull = randn(longn, ncv)
    #Vfull = Matrix(1.0I, longn, ncv)
    V = copy(Vfull)
    V1 = copy(Vfull) 
    V2 = copy(Vfull)
    longwork = zeros(longn) 
    longwork1 = copy(longwork)
    longwork2 = copy(longwork)
  
    ArpackInJulia.dorm2r!(Val(:R), Val(:N), longn, ncv, nconv, 
      Q, tau, @view(V[1:longn, 1:ncv]), @view(longwork[1:longn]))

    _dorm2r_blas!('R', 'N', longn, ncv, nconv, Q, tau, @view(V1[1:longn, 1:ncv]), 
      @view(longwork1[1:longn]))

    ArpackInJulia.dorm2r!(Val(:R), Val(:N), longn, ncv, nconv, 
      copy(Q), copy(tau), V2, longwork2)

    test_diffs("V, V1", V, V1)
    test_diffs("V1, V2", V1, V2)
    @test V == V1 
    @test V1 == V2 
  end 

  @testset "dorg2r inplace orthogonalize" begin 
    function compare_qr_orth(A::StridedMatrix)
      m, n = size(A)
      @assert(m >= n)
      work = zeros(n) 
      work2 = zeros(n) 
      Qinfo, tau = _simple_dgeqr2(A)
      Q2 = copy(Qinfo) # copy for blas
      Q = ArpackInJulia.dorg2r!(Qinfo, tau, work)
      Q2 = _dorg2r_blas!(m, n, length(tau), Q2, stride(Q2,2), tau, work2)
      @test Q == Q2
      @test work == work2 
    end 
    n = 7 
    compare_qr_orth(Matrix(SymTridiagonal(zeros(n), ones(n-1))))    
  end 

end 
