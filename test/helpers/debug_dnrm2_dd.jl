##
# There is a case where we get NaN/Inf from our implementation 
# and not from Blas
using Test
using ArpackInJulia
using LinearAlgebra

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

_eval_dnrm2_blas(
  a::StridedVecOrMat{Float64},
) = ccall((LinearAlgebra.BLAS.@blasfunc("dnrm2_"), LinearAlgebra.BLAS.libblas), Float64,
  (Ref{LinearAlgebra.BlasInt}, Ptr{Float64}, Ref{LinearAlgebra.BlasInt}),
  length(a), a, stride(a,1))


function eval_nrm2_dd(a::AbstractVector{T}) where T
  min,max = extrema(abs, a)
  # here's the question, do we have enough precision in strict double-double,
  # or do we need to scale?
  if max >= sqrt(prevfloat(Inf))
    #scale = (2*one(T))^(-511)
    scale = one(T)
  else
    scale = one(T)
  end 
  len = length(a) 
  offset = 1
  ss1 = (zero(T), zero(T))
  ss2 = (zero(T), zero(T))
  ss3 = (zero(T), zero(T))
  ss4 = (zero(T), zero(T))
  @inbounds while offset+8 <= len
    ss1 = ArpackInJulia.add_sum_sq(ss1, a[offset]*scale)
    ss2 = ArpackInJulia.add_sum_sq(ss2, a[offset+1]*scale)
    ss3 = ArpackInJulia.add_sum_sq(ss3, a[offset+2]*scale)
    ss4 = ArpackInJulia.add_sum_sq(ss4, a[offset+3]*scale)
    ss1 = ArpackInJulia.add_sum_sq(ss1, a[offset+4]*scale)
    ss2 = ArpackInJulia.add_sum_sq(ss2, a[offset+5]*scale)
    ss3 = ArpackInJulia.add_sum_sq(ss3, a[offset+6]*scale)
    ss4 = ArpackInJulia.add_sum_sq(ss4, a[offset+7]*scale)
    offset += 8 
  end
  #@show ss1, ss2, ss3, ss4
  @inbounds for i=offset:len
    ss1 = ArpackInJulia.add_sum_sq(ss1, a[i]*scale)
  end
  ss1 = ArpackInJulia.add_dddd_dd(ss1, ss3)
  ss1 = ArpackInJulia.add_dddd_dd(ss1, ss2) 
  ss1 = ArpackInJulia.add_dddd_dd(ss1, ss4) 

  
  
  # lots of work to get the last bit right...
  # inline sqrt_dd_dd from DoubleFloats.jl
  r = inv(sqrt(ss1[1]))
  h = (ss1[1]*0.5, ss1[2]*0.5)
  r2 = ArpackInJulia.two_prod(r, r) 
  hr2 = ArpackInJulia.mul_dddd_dd(h, r2)
  radj = ArpackInJulia.sub_fpdd_dd(0.5, hr2)
  radj = ArpackInJulia.mul_dddd_dd(radj, (r, 0.0))
  r = ArpackInJulia.add_fpdd_dd(r, radj)

  r2 = ArpackInJulia.mul_dddd_dd(r,r)
  hr2 = ArpackInJulia.mul_dddd_dd(h, r2)
  radj = ArpackInJulia.sub_fpdd_dd(0.5, hr2)
  radj = ArpackInJulia.mul_dddd_dd(radj, r)
  r = ArpackInJulia.add_dddd_dd(r, radj)

  r = ArpackInJulia.mul_dddd_dd(r, ss1)
  return r[1]/scale 
end  

@testset "Long vector failure" begin
  import Random
  rng = Random.MersenneTwister(1)
  sz = 1000001
  vals = exp.(exp.(1.3*randn(rng, sz))) # very wide dynamic range... 
  n1 = _eval_dnrm2_blas(vals)
  n2 = eval_nrm2_dd(vals)
  @show norm(vals), n1, n2
  @test isfinite(n2)
end  

##
@testset "_dnrm2 Compare to OpenBLAS on x86 / x86_64" begin
  @testset "Random vectors" begin
    import Random
    rng = Random.MersenneTwister(2)
    sizes = [1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,32,33,63,100,151,252,1023,2042,8123,32191,121892,1000001]
    ntrials = 10000
    nfloatsoff = 0 
    nrand = 0 
    for t =1:ntrials
      sz = rand(rng, sizes)
      vals = exp.(exp.(randn(rng, sz))) # very wide dynamic range... 
      nrand += sz + 1
      n1 = _dnrm2_blas(vals)
      n2 = eval_nrm2_dd(vals)
      ndiff = abs(floatsbetween(n1,n2))
      nfloatsoff += ndiff 
      if !isfinite(n2)
        @show nrand, length(vals), n1, n2, norm(vals)
      end
      @test abs(floatsbetween(n1,n2)) <= 0
      

      vals2 = vals ./ maximum(vals)
      n3 = _dnrm2_blas(vals2)
      n4 = eval_nrm2_dd(vals2)
      ndiff = abs(floatsbetween(n3,n4))
      nfloatsoff += ndiff 
      @test abs(floatsbetween(n3,n4)) <= 0
    end
    @test nfloatsoff <= 1
  end
end 