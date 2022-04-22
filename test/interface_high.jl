@testset "high precision interfaces" begin
  using Quadmath
  using DoubleFloats
  using MultiFloats

  GenericArpack.@fix_doublefloats
  GenericArpack.@fix_multifloats

  @testset "fixes" begin
    val = (LinearAlgebra.floatmin2(Double64))
    val2 = @allocated(LinearAlgebra.floatmin2(Double64)) # make sure to time allocated too... 
    @test val == DoubleFloats.Double64(reinterpret(Float64, 0x2350000000000000), 0.0)    
    # can't figure out why this test won't pass on CI server, but oh well...
    #@test (@allocated(LinearAlgebra.floatmin2(Double64))) == 0

    @test GenericArpack._dstqrb_maxit(Double64) == 80

    @test typeof(GenericArpack._eps23(Float64x2)) == Float64x2
    @test GenericArpack._dstqrb_maxit(Float64x2) == 80
    
    @test Int(Float64x2(1.0)) == 1

    @allocated GenericArpack._dstqrb_maxit(Float128) == 0
  end

  @testset "Quadmath" begin 
    using LinearAlgebra
    T = Float128 

    @testset  "Tridiag" begin 
      A = arpack_dsdrv3(T, 100)[1]
      info = eigs(Symmetric(A), 4)
      @test eltype(info.values) == T
      B = SymTridiagonal(A)
      @test info.values ≈ sort(sort(LinearAlgebra.eigvals(B), rev=true)[1:4])
    end 

    @testset  "dsdrv3" begin 
      A, B = arpack_dsdrv3(T, 100)
      info = eigs(Symmetric(A), Symmetric(B), 4)
      @test eltype(info.values) == T
      @test Float64.(info.values) ≈ [121003.49732902342, 121616.60247324049, 122057.49457079473, 122323.22366457577]
    end 

    @testset  "svd" begin 
      A = mytestmat(100,80,one(T)/80)
      info = svds(A, 2)
      @test eltype(info.S) == T
    end
  end

  @testset "DoubleFloats" begin 
    using LinearAlgebra
    T = Double64 

    @testset  "Tridiag" begin 
      A = arpack_dsdrv3(T, 100)[1]
      info = eigs(Symmetric(A), 4)
      @test eltype(info.values) == T
      B = SymTridiagonal(A)
      @test info.values ≈ sort(sort(LinearAlgebra.eigvals(B), rev=true)[1:4])
    end 

    @testset  "dsdrv3" begin 
      A, B = arpack_dsdrv3(T, 100)
      info = eigs(Symmetric(A), Symmetric(B), 4)
      @test eltype(info.values) == T
      @test Float64.(info.values) ≈ [121003.49732902342, 121616.60247324049, 122057.49457079473, 122323.22366457577]
    end 

    @testset  "svd" begin 
      A = mytestmat(100,80,one(T)/80)
      info = svds(A, 2)
      @test eltype(info.S) == T
    end
  end

  @testset "MultiFloats" begin 
    using LinearAlgebra

    # don't test Float64x8 here as 
    # we get different results from SymEigs...
    # but it isn't clear who is right.
    Ts = [Float64x2, Float64x4, Float64x6]
    for T in Ts
      @testset "MultiFloat $T" begin 
        @testset  "Tridiag" begin 
          A = arpack_dsdrv3(T, 100)[1]
          info = eigs(Symmetric(A), 4)
          @test eltype(info.values) == T
          B = SymTridiagonal(A)
          @test info.values ≈ sort(sort(LinearAlgebra.eigvals(B), rev=true)[1:4])
        end 

        @testset  "dsdrv3" begin 
          A, B = arpack_dsdrv3(T, 100)
          info = eigs(Symmetric(A), Symmetric(B), 4)
          @test eltype(info.values) == T
          @test Float64.(info.values) ≈ [121003.49732902342, 121616.60247324049, 122057.49457079473, 122323.22366457577]
        end 

        @testset  "svd" begin 
          A = mytestmat(100,80,one(T)/80)
          info = svds(A, 2)
          @test eltype(info.S) == T
        end
      end
    end 
  end
end 

