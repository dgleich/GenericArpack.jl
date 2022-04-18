@testset "interface complex" begin 
  using SparseArrays, LinearAlgebra, StableRNGs
  n = 15 
  rng = StableRNG(1234)
  A = sprandn(rng, ComplexF64, n, n, 0.25) |> A -> A + A' 
  Avals,Avecs = eigen(Matrix(A))
  #@show A.colptr, A.rowval, A.nzval

  k = 3 
  #vals, vecs = eigs(Symmetric(A), k; debug=set_debug_high!(ArpackDebug()))
  vals, vecs = eigs(Hermitian(A), k; debug=(ArpackDebug()))
  @test vals ≈ sort(sort(Avals, by=abs,rev=true)[1:k])
  vals, vecs = eigs(Hermitian(A), k; ritzvec=false, which=:LA)
  @test vals ≈ sort(sort(Avals, rev=true)[1:k])

  vals, vecs = eigs(ComplexF32, Hermitian(A), k; ritzvec=false, which=:LA)
  @test eltype(vals) == Float32
  @test vals ≈ Float32.(sort(sort(Avals, rev=true)[1:k]))

  fop = ArpackSimpleFunctionOp((y,x) -> mul!(y,A,x), n)
  vals, vecs = hermeigs(fop, k)
  @test vals ≈ sort(sort(Avals, by=abs,rev=true)[1:k])
  #vals, vecs = eigs(Float64, Symmetric(A), k; which=:SM)
  #@test vals ≈ sort(sort(Avals, by=abs)[1:k])
  
  #vals, vecs = eigs(Symmetric(A), k; which=:LA)
end 

@testset "complex svd" begin 
  @testset "using matrix" begin 
    A = mytestmat(10,8,-3.0im)
    U,s,V = svds(A, 2; which=:BE)
    @test s ≈ [0.9640220546797924, 6.0316491543925785]
    check_svd(A,U,s,V;tol=30) 
  end 

  @testset "using ArpackNormalFunctionOp" begin 
    A = mytestmat(10,8,-3.0im)
    fop = ArpackNormalFunctionOp((y,x) -> mul!(y, A, x), (y,x) -> mul!(y, adjoint(A), x), size(A)...)
    U,s,V = complexsvds(fop, 2; which=:LM)
    @test s ≈ [6.0316491543925785, 1.9374903374733126]
    check_svd(A,U,s,V;tol=30) 
  end 
end 