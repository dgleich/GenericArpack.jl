#=
Goal: Compare our _dnrm2_unroll_ext to _drnm2 in OpenBLAS
=# 

#=
# Codes for standalone testing... 
using Test
using GenericArpack

""" Return the number of floating point values between two values, including
one of the endpoints. """
function floatsbetween(a::T,b::T) where {T <: Float64}
  # based on https://github.com/JuliaLang/julia/commit/79920db3ac9f1e1915d836c82bc401ba0cdc322f
  # https://stackoverflow.com/questions/3587880/number-of-floats-between-two-floats
  IntT = Int64 # inttype(T)
  ia = reinterpret(IntT, a)
  ib = reinterpret(IntT, b)
  ia = ifelse(ia < zero(IntT), ia ⊻ typemax(IntT), ia)
  ib = ifelse(ib < zero(IntT), ib ⊻ typemax(IntT), ib)
  return ib-ia
end

_dnrm2_blas(
  a::StridedVecOrMat{Float64},
) = ccall((LinearAlgebra.BLAS.@blasfunc("dnrm2_"), LinearAlgebra.BLAS.libblas), Float64,
  (Ref{LinearAlgebra.BlasInt}, Ptr{Float64}, Ref{LinearAlgebra.BlasInt}),
  length(a), a, stride(a,1))

# study a case where we get different nroms.
begin 
  using Random
  Random.seed!(1)
  vals = exp.(exp.(randn(2)))
  n1 = _dnrm2_blas(vals)
  n2 = GenericArpack._dnrm2_unroll_ext(vals)
  @show vals
  @show n1, n2, floatsbetween(n1, n2)
end
=#

##
if Sys.ARCH==:x86_64 || Sys.ARCH==:x86
  @testset "_dnrm2 Compare to OpenBLAS on x86 / x86_64" begin
    import LinearAlgebra
    _dnrm2_blas(
      a::StridedVecOrMat{Float64},
    ) = ccall((LinearAlgebra.BLAS.@blasfunc("dnrm2_"), LinearAlgebra.BLAS.libblas), Float64,
      (Ref{LinearAlgebra.BlasInt}, Ptr{Float64}, Ref{LinearAlgebra.BlasInt}),
      length(a), a, stride(a,1))

    @testset "length 1000 vectors" begin
      x = zeros(1000); fill!(x, 0.03162277660168379)
      n1 = _dnrm2_blas(vals)
      n2 = GenericArpack.norm2(Float64, vals)
      @test floatsbetween(n1,n2) == 0 
      if floatsbetween(n1,n2) != 0
        @show n1, n2, floatsbetween(n1,n2)
      end
    end      

    @testset "Random vectors" begin
      import StableRNGs
      #rng = Random.MersenneTwister(2)
      rng = StableRNG(2)
      sizes = [1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,32,33,63,100,151,252,1023,2042,8123,32191,121892,1000001]
      ntrials = 10000
      nfloatsoff = 0 
      maxtotaldiff = 10
      maxdiff = 1 
      nrand = 0 
      for t =1:ntrials
        sz = rand(rng, sizes)
        vals = exp.(exp.(randn(rng, sz))) # very wide dynamic range... 
        nrand += sz + 1
        n1 = _dnrm2_blas(vals)
        n2 = GenericArpack._dnrm2_unroll_ext(vals)
        ndiff = abs(floatsbetween(n1,n2))
        nfloatsoff += ndiff 
        if !isfinite(n2)
          @show nrand, length(vals), n1, n2, norm(vals)
        end
        @test abs(floatsbetween(n1,n2)) <= maxdiff 
        if ndiff != 0
          @show n1, n2, sz, norm(vals)
        end
        

        vals2 = vals ./ maximum(vals)
        n3 = _dnrm2_blas(vals2)
        n4 = GenericArpack._dnrm2_unroll_ext(vals2)
        ndiff = abs(floatsbetween(n3,n4))
        nfloatsoff += ndiff 
        @test abs(floatsbetween(n3,n4)) <= maxdiff 

        if ndiff != 0
          @show n3, n4, sz, norm(vals2)
        end
      end
      @test nfloatsoff <= maxtotaldiff
    end
  end 
end 