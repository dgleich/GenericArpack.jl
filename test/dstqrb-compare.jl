
function checked_dstqrb(T::SymTridiagonal)
  T1 = copy(T)
  T2 = copy(T)
  work = zeros(2*size(T,1))
  z1 = zeros(size(T,1))
  z2 = zeros(size(T,2))
  conceptual_dstqrb!(T1; z=z1, work)
  arpack_dstqrb!(size(T,1), T2.dv, T2.ev, z2, work)
  return T1, T2, z1, z2
end
function simple_eval(T1,T2,z1,z2)
  @show maximum(abs.(map(relfloatsbetween, T1.dv, T2.dv)))
  @show maximum(abs.(map(relfloatsbetween, z1, z2)))
  @show norm(T1-T2,Inf)
  @show norm((z1)-(z2),Inf)
  @show norm(abs.(z1)-abs.(z2),Inf)
end

## This is one of the problems
@testset "sign-probs" begin
  simple_z_sign_test(T1,T2,z1,z2) = @test norm(z1-z2,Inf) < 1e-8
  using Random
  Random.seed!(76102)
  T = SymTridiagonal(randn(9,9) |>  A-> A+A') # random symmetric matrix
  simple_z_sign_test(checked_dstqrb(T)...)

  Random.seed!(10973)
  T = SymTridiagonal(randn(11,11) |>  A-> A+A') # random symmetric matrix
  simple_z_sign_test(checked_dstqrb(T)...)

  Random.seed!(796)
  T = SymTridiagonal(randn(15,15) |>  A-> A+A') # random symmetric matrix
  simple_z_sign_test(checked_dstqrb(T)...)

  Random.seed!(8)
  T = SymTridiagonal(randn(1025,1025) |>  A-> A+A')
  simple_z_sign_test(checked_dstqrb(T)...)
end

## Let's repeat the above tests but check against arpack
@testset "simple_arpack_checks" begin
  function strict_eval(T1,T2,z1,z2)
    @test maximum(map(floatsbetween,T1.dv,T2.dv)) .<= 1
    @test norm(z1-z2,Inf) < 1e-8
  end
  @testset "subdiagonals" begin
    T = SymTridiagonal(ones(15,15))
    fill!(T.ev, 0) #
    T.ev[8] = 10
    conceptual_dstqrb!(T)
    @test all(T.ev .== 0)
    fill!(T.ev, 0)
    fill!(T.dv, 0)
    fill!(T.ev, 2*eps(1.0))
    strict_eval(checked_dstqrb(T)...)

    fill!(T.ev, 2*sqrt(eps(1.0)))
    fill!(T.dv, 0)
    T.dv[1] = 1.0
    strict_eval(checked_dstqrb(T)...)

    fill!(T.ev, 2*sqrt(eps(1.0)))
    fill!(T.dv, 0)
    T.dv[end] = 1.0
    strict_eval(checked_dstqrb(T)...)
  end

  @testset "breakpoints" begin
    # need even to get two even splits to avoid zero...
    T = SymTridiagonal(ones(14,14))
    T.ev[4] = eps(1.0)/10
    T.dv[5:end] .+= 6
    strict_eval(checked_dstqrb(T)...)
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
        strict_eval(checked_dstqrb(S)...)
      end
    end
  end
  ##
  @testset "zeros" begin
    for n=1:10
      strict_eval(checked_dstqrb(
        SymTridiagonal(zeros(n), zeros(n-1)))...)
    end
  end
  @testset "ones" begin
    # skip odd values because they given an eigenvalue near 0 that is problematic for precision
    for n=[2,8,16,128,500,512,1000,1024]
      strict_eval(checked_dstqrb(SymTridiagonal(zeros(n), ones(n-1)))...)
    end
    # check for the accuracy in the eigenvalue near zero
    for n=[3,9,15,17,127,129,501,511,513,1001,1025,2001]
      strict_eval(checked_dstqrb(SymTridiagonal(zeros(n), ones(n-1)))...)
    end
  end
  ##
  @testset "laplacian" begin
    for n=[1,2,3,5,8,11,13,16,17,99,127,128,129,500,511,512,513,1000,1023,1024,1025]
      strict_eval(checked_dstqrb(SymTridiagonal(2*ones(n), -ones(n-1)))...)
    end
  end
end
