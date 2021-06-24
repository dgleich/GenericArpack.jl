using LinearAlgebra
conceptual_dstqrb! = ArpackInJulia.conceptual_dstqrb!
## Helper functions
# all subidagonals are zero.
function reverse_T(T::SymTridiagonal)
  return SymTridiagonal(reverse(T.dv), reverse(T.ev))
end
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

""" Count the number of floats between but treat the range [-ref,ref] as two
fractional floating point values,

relfloatsbetween(-ref,ref;ref) = 2.0
relfloatsbetween(-prevfloat(ref),nextfloat(ref);ref) = 4.0
relfloatsbetween(ref/4,ref/2;ref) = 0.25
"""
function relfloatsbetween(a::T,b::T;ref=eps(T)) where {T <: Float64}
  # this is awful code at the moment :(
  if b < a
    flipsign = -1
    a,b = b,a # swap so we are in order to minimize cases below
  else
    flipsign = 1
  end

  rval = 0.0
  if a < -ref && b < -ref # then since a < b, we can just return this directly
    rval += floatsbetween(a,b)
  elseif a <= -ref && b <= ref
    rval += floatsbetween(a, -ref)
    rval += (b-(-ref))/ref
  elseif a <= -ref && b > ref
    rval += floatsbetween(a, -ref)
    rval += floatsbetween(ref, b)
    rval += 2.0
  elseif a <= ref && b <= ref
    # a in [-ref,ref], b in [-ref,ref]
    rval += (b-a)/ref
  elseif a <= ref && b > ref
    rval += (ref-a)/ref
    rval += floatsbetween(ref, b)
  else
    # a >= ref, b >= ref
    rval += floatsbetween(a,b)
  end

  return flipsign*rval
end
##
@testset "small" begin
  T = SymTridiagonal(ones(1,1))
  conceptual_dstqrb!(T)
  @test T.dv[1] == 1

  T = SymTridiagonal(ones(2,2))
  z = zeros(2)
  conceptual_dstqrb!(T;z)
  @test T.dv[1] == 0
  @test T.dv[2] == 2
  @test abs(z[1]) == 1/sqrt(2)
  @test abs(z[2]) == 1/sqrt(2)
end
##
@testset "subdiagonals" begin
  T = SymTridiagonal(ones(15,15))
  fill!(T.ev, 0) #
  T.ev[8] = 10
  conceptual_dstqrb!(T)
  @test all(T.ev .== 0)
  fill!(T.ev, 0)
  fill!(T.dv, 0)
  fill!(T.ev, 2*eps(1.0))
  conceptual_dstqrb!(T)
  @test all(T.ev .== 0)

  fill!(T.ev, 2*sqrt(eps(1.0)))
  fill!(T.dv, 0)
  T.dv[1] = 1.0
  conceptual_dstqrb!(T)
  @test all(T.ev .== 0)

  fill!(T.ev, 2*sqrt(eps(1.0)))
  fill!(T.dv, 0)
  T.dv[end] = 1.0
  conceptual_dstqrb!(T)
  @test all(T.ev .== 0)
end

@testset "breakpoints" begin
  # need even to get two even splits to avoid zero...
  T = SymTridiagonal(ones(14,14))
  T.ev[4] = eps(1.0)/10
  T.dv[5:end] .+= 6
  L,z = conceptual_dstqrb!(copy(T))
  @test all(z[1:4] .== 0)
  @test maximum(map(floatsbetween, eigvals(T), L.dv) .<= 10)
end

##
@testset "random-breakpoints" begin
  using Random
  Random.seed!(0)
  for trial=1:100
    n = 15
    T = SymTridiagonal(randn(129,129) |> A->A*A') # sym pos def...
    for subtrial = 1:10
      S = copy(T)
      bps = rand(1:n-1)
      for bpi = 1:bps
        S.ev[bpi] = eps(1.0)*S.dv[bpi]*S.dv[bpi+1]/rand(1:3)
      end
      L,z = conceptual_dstqrb!(copy(S))
      @test all(L.dv .>= 0)
      @test maximum(map(floatsbetween, eigvals(S), L.dv) .<= 10)
    end
  end
end
##
@testset "zeros" begin
  for n=1:10
    T,z = conceptual_dstqrb!(SymTridiagonal(zeros(n), zeros(n-1)))
    @test all(T.dv .== 0)
    @test z[end] == 1.0
    @test all(z[1:end-1] .== 0.0)
  end
end
@testset "ones" begin
  # skip odd values because they given an eigenvalue near 0 that is problematic for precision
  for n=[2,8,16,128,500,512,1000,1024]
    T,z = conceptual_dstqrb!(SymTridiagonal(zeros(n), ones(n-1)))
    lams = 2*cos.(pi.*(1:n)./(n+1))
    sort!(lams) # sort lams
    @test maximum(abs.(map(relfloatsbetween, lams, T.dv))) .<= 2*n
    @test norm(abs.(z) - abs.(sqrt(2/(n+1)).*sin.(n.*(1:n).*pi/(n+1)))) .<= 2*n*eps(1.0)
  end
  # check for the accuracy in the eigenvalue near zero
  for n=[3,9,15,17,127,129,501,511,513,1001,1025,2001]
    T,z = conceptual_dstqrb!(SymTridiagonal(zeros(n), ones(n-1)))
    #@show abs(T.dv[div(n,2)+1])
    @test abs(T.dv[div(n,2)+1]) .<= 10*eps(1.0)
    @test norm(abs.(z) - abs.(sqrt(2/(n+1)).*sin.(n.*(1:n).*pi/(n+1)))) .<= 2*n*eps(1.0)
  end
end
##
@testset "laplacian" begin
  for n=[1,2,3,5,8,11,13,16,17,99,127,128,129,500,511,512,513,1000,1023,1024,1025]
    T,z = conceptual_dstqrb!(SymTridiagonal(2*ones(n), -ones(n-1)))
    lams = -2*cos.(pi.*(1:n)./(n+1)).+2.0
    sort!(lams) # sort lams
    #@show norm(T.dv .- lams,Inf)
    @test norm(lams-T.dv,Inf) .<= 2*sqrt(n)*eps(1.0)
    @test norm(abs.(z) - abs.(sqrt(2/(n+1)).*sin.(n.*(1:n).*pi/(n+1)))) .<= 2*n*eps(1.0)
  end
end

## TODO, look into why 7 is so inaccurate in the above case.
