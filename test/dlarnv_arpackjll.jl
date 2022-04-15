@testset "dlarnv" begin # not strictly arpackjll, but this is a good place
  _dlarnv_check_blas!(idist::Int,
     iseed::Ref{NTuple{4,Int}},
    n::Int,
    x::Vector{Float64}) =
  ccall((LinearAlgebra.BLAS.@blasfunc("dlarnv_"), LinearAlgebra.BLAS.libblas), Cvoid,
    (Ref{LinearAlgebra.BlasInt}, Ptr{LinearAlgebra.BlasInt}, Ref{LinearAlgebra.BlasInt}, Ptr{Float64}),
    idist, iseed, n, x)

  _dlarnv_check_blas!(idist::Int,
    iseed::Ref{NTuple{4,Int}},
   n::Int,
   x::Vector{Float32}) =
  ccall((LinearAlgebra.BLAS.@blasfunc("slarnv_"), LinearAlgebra.BLAS.libblas), Cvoid,
   (Ref{LinearAlgebra.BlasInt}, Ptr{LinearAlgebra.BlasInt}, Ref{LinearAlgebra.BlasInt}, Ptr{Float32}),
   idist, iseed, n, x)    

  _dlarnv_check_blas!(idist::Int,
   iseed::Ref{NTuple{4,Int}},
    n::Int,
    x::Vector{ComplexF64}) =
  ccall((LinearAlgebra.BLAS.@blasfunc("zlarnv_"), LinearAlgebra.BLAS.libblas), Cvoid,
    (Ref{LinearAlgebra.BlasInt}, Ptr{LinearAlgebra.BlasInt}, Ref{LinearAlgebra.BlasInt}, Ptr{ComplexF64}),
    idist, iseed, n, x)       

  seed1 = Base.RefValue(tuple(1,3,5,7))
  seed2 = Base.RefValue(tuple(1,3,5,7))
  for n=0:4097
    z1 = zeros(n)
    z2 = zeros(n)
    _dlarnv_check_blas!(2, seed1, length(z1), z1)
    GenericArpack._dlarnv_idist_2!(seed2, length(z2), z2)
    @test seed1[] == seed2[]
    @test z1 == z2
  end

  seed1 = Base.RefValue(tuple(1,3,5,7))
  seed2 = Base.RefValue(tuple(1,3,5,7))
  for n=0:4097
    z1 = zeros(Float32, n)
    z2 = zeros(Float32, n)
    _dlarnv_check_blas!(2, seed1, length(z1), z1)
    GenericArpack._dlarnv_idist_2!(seed2, length(z2), z2)
    @test seed1[] == seed2[]
    @test z1 == z2
  end  

  seed1 = Base.RefValue(tuple(1,3,5,7))
  seed2 = Base.RefValue(tuple(1,3,5,7))
  for n=0:4097
    z1 = zeros(ComplexF64, n)
    z2 = zeros(ComplexF64, n)
    _dlarnv_check_blas!(2, seed1, length(z1), z1)
    GenericArpack._dlarnv_idist_2!(seed2, length(z2), z2)
    @test seed1[] == seed2[]
    @test z1 == z2
  end  

  using Random
  Random.seed!(0)
  for t=1:5000 
    s1 = Base.RefValue(tuple(rand(0:4095),rand(0:4095),rand(0:4095),rand(0:4095)))
    s2 = Base.RefValue(s1[])
    n = rand(1:5)
    z1 = zeros(n)
    z2 = zeros(n)
    _dlarnv_check_blas!(2, s1, length(z1), z1)
    GenericArpack._dlarnv_idist_2!(s2, length(z2), z2)
    @test seed1[] == seed2[]
    @test z1 == z2
  end
end