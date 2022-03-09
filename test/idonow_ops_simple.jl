@testset "ArpackOp" begin 
  @testset "ArpackSimpleOp" begin 
    using LinearAlgebra
    using Random
    Random.seed!(0)
    A = Diagonal(1.0:10.0)
    idonow = ArpackInJulia.ArpackSimpleOp(A)
    y = zeros(size(A,1))
    x = randn(size(A,1))
    ArpackInJulia.opx!(y,idonow,x)
    @test y == A*x 
  end 
end 