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