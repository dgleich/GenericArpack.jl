@testset "ArpackOp" begin 
  using LinearAlgebra
  function exercise_op(T, op, A, B=I)
    n = size(op)  
    if B isa UniformScaling
      # standard op
      x = ones(T, n) 
      y = zeros(T, n)
      opx!(y, op, x)
      @test y == A*x
    else
      # generalized op... or shift-invert op
      #shift = GenericArpack.shift(T, op) 
    end 
    
  end 
  @testset "ArpackSimpleOp" begin 
    using LinearAlgebra
    using Random
    Random.seed!(0)
    A = Diagonal(1.0:10.0)
    idonow = GenericArpack.ArpackSimpleOp(A)
    y = zeros(size(A,1))
    x = randn(size(A,1))
    GenericArpack.opx!(y,idonow,x)
    @test y == A*x 
  end 

  @testset "ArpackSymmetricGeneralizedOp" begin 
    n = 100
    A = Tridiagonal(-ones(n-1),2*ones(n),-ones(n-1)).*(n+1)
    B = Tridiagonal(ones(n-1),4*ones(n),ones(n-1))*(1/(6*(n+1)))
    #B = SymTridiagonal(4*ones(n),ones(n-1))*(1/(6*(n+1)))

    op1 = GenericArpack.ArpackSymmetricGeneralizedOp(A,factorize(B),B)
    # TODO, this would ideally work with factorize
    # factorize(Symmetric(B)) 
    # but that hits issue https://github.com/JuliaLang/julia/issues/44973
    op2 = GenericArpack.ArpackSymmetricGeneralizedOp(Symmetric(A),lu!(copy(B)),Symmetric(B))

    using Random
    Random.seed!(0)

    y = zeros(Float32, size(A,1))
    x = randn(Float32, size(A,1))
    x1 = copy(x)

    x = copy(x1) 
    GenericArpack.genopx!(y,op1,x)
    @test y ≈ B\(A*x1)

    x = copy(x1)
    GenericArpack.genopx!(y,op2,x)
    @test y ≈ B\(A*x1)

    x = copy(x1)
    GenericArpack.bx!(y,op1,x)
    @test y ≈ B*x1

    x = copy(x1)
    GenericArpack.bx!(y,op2,x)
    @test y ≈ B*x1
  end 

  @testset "ArpackNormalOp" begin 
    A = randn(20,8)
    op = ArpackNormalOp(Float64, A)
    B = Matrix(Float64, op) 
    @test size(op) == 8 
    @test GenericArpack._uv_size(op) == (20,8)
    @test A'*A ≈ B 
  end 

  @testset "ArpackNormalFunctionOp" begin 
    A = randn(20,8)
    myav = (y,x) -> mul!(y, A, x)
    myatv = (y,x) -> mul!(y, adjoint(A), x)
    op = ArpackNormalFunctionOp(myav, myatv, size(A)...)
    B = Matrix(Float64, op) 
    @test size(op) == 8 
    @test GenericArpack._uv_size(op) == (20,8)
    @test A'*A ≈ B 
  end 
end 