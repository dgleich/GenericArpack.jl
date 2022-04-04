#=
running tests

Options:
debug - Enable debugging macros (@debug statements)
arpackjll - Check against the Arpack library itself
full - Check a wider range of options

Example:

In ArparkInJulia environment...

    Pkg.test(;test_args=["debug"])

Elsewhere

    using Pkg; Pkg.test("ArpackInJulia"; test_args=["arpackjll"])

=#    

# This needs to be before ArpackInJulia to get debug macros enabled
if "debug" in ARGS
  ENV["JULIA_DEBUG"] = "ArpackInJulia,Main"
end

using Test
using ArpackInJulia

# utility can depend on test... 
include("utility.jl")

# want to run these first... 
# Generally, the strategy is to write "method_simple" and "method_arpackjll" here...
if false # switch to true while developing
  @testset "development..." begin
    # include("arpackjll.jl") # uncomment to develop arpackjll tests
    # arpack_set_debug_high()

  end
  #exit(0)
end

#@testset "external interfaces" begin 
#end 

include("macros.jl")

@testset "simple features" begin
  include("idonow_ops_simple.jl")

  @testset "dsaupd" begin 
    include("dsaupd_simple.jl")
    include("dsaupd_idonow.jl")
  end 

  include("dsgets_simple.jl")
  
  @testset "dsconv" begin
    @test_throws ArgumentError ArpackInJulia.dsconv(6, zeros(5), ones(5), 1e-8)
    @test_nowarn ArpackInJulia.dsconv(5, zeros(5), ones(5), 1e-8)
    @test_nowarn ArpackInJulia.dsconv(5, zeros(5), ones(5), 1e-8;
                    stats=ArpackInJulia.ArpackStats())
    stats = ArpackInJulia.ArpackStats()
    @test ArpackInJulia.dsconv(15, ones(15), zeros(15), 1e-8; stats)==15
    #@test stats.tconv > 0  # we did record some time
    t1 = stats.tconv
    @test ArpackInJulia.dsconv(25, ones(25), zeros(25), 1e-8; stats)==25
    #@test stats.tconv > t1  # we did record some more time

    # TODO, add more relevant tests here.
  end

  @testset "dsortr" begin
    x = randn(10)
    # check obvious error
    @test_throws ArgumentError ArpackInJulia.dsortr(:LA, false, length(x)+1, x, zeros(0))
    @test_throws ArgumentError ArpackInJulia.dsortr(:LA, true, length(x), x, zeros(0))
    @test_nowarn ArpackInJulia.dsortr(:LA, false, length(x), x, zeros(0))
    @test issorted(x)

    x = randn(8)
    @test_nowarn ArpackInJulia.dsortr(:LA, false, length(x), x, zeros(0))
    @test issorted(x)
  end

  @testset "dseigt" begin
    using LinearAlgebra
    n = 34 # don't use odd values as they have a zero eigenvalue that is annoying.
    lams = 2*cos.(pi.*(1:n)./(n+1))
    sort!(lams) # sort lams
    T = SymTridiagonal(zeros(n), ones(n-1))
    H = zeros(n, 2)
    H[2:end,1] .= T.ev
    H[:,2] .= T.dv
    rnorm = 128.0
    eigval = randn(n)
    bounds = randn(n) 
    workl = randn(n,3)
    debug = ArpackDebug(logfile=IOBuffer())
    ArpackInJulia.set_debug_high!(debug)
    stats = ArpackStats()
    ArpackInJulia.dseigt!(rnorm, n, H, n, eigval, bounds, workl, nothing; debug, stats)

    # check against the analytical results...
    @test maximum(abs.(map(relfloatsbetween, lams, eigval))) .<= 2*n
    @test norm(bounds - abs(rnorm)*abs.(sqrt(2/(n+1)).*sin.(n.*(1:n).*pi/(n+1)))) <= 2*n*eps(1.0)*rnorm  
  end

  @testset "dnrm2" begin 
    @test ArpackInJulia._dnrm2_unroll_ext(zeros(5)) == 0 # make sure we get zero
  end 

  include("dgetv0_simple.jl")
  include("dsaitr_simple.jl")
end

@testset "fancy features" begin
  @testset "override methods" begin 
    include("dgetv0_override.jl")
  end 
end

@testset "dstqrb" begin
include("dstqrb.jl")
end

##
# To run, use
# using Pkg; Pkg.test("ArpackInJulia"; test_args=["arpackjll"])
#
if "arpackjll" in ARGS
  @testset "blas-lapack" begin
    # comparisons with blas/lapack
    include("dlarnv_arpackjll.jl")
    include("blas-qr-compare.jl")
  end
  

  include("arpackjll.jl") # get the helpful routines to call stuff in "arpackjll"
  @testset "arpackjll" begin

    @testset "dsconv" begin
      soln = arpack_dsconv(10, ones(10), zeros(10), 1e-8)
      @test ArpackInJulia.dsconv(10, ones(10), zeros(10), 1e-8)==soln

      soln = arpack_dsconv(10, zeros(10), ones(10), 1e-8)
      @test ArpackInJulia.dsconv(10, zeros(10), ones(10), 1e-8)==soln
    end

    @testset "dsortr" begin
      #soln = arpack_dsconv(10, ones(10), zeros(10), 1e-8)
      #@test ArpackInJulia.dsconv(10, ones(10), zeros(10), 1e-8)==soln
      x = collect(9.0:-1:1)
      x1 = copy(x); x2 = copy(x)
      ArpackInJulia.dsortr(:LM, false, length(x), x1, zeros(0))
      arpack_dsortr(:LM, false, length(x), x2, zeros(0))
      @test x1==x2

      x = collect(9.0:-1:1)
      y = collect(1.0:9)
      x1 = copy(x); x2 = copy(x)
      y1 = copy(y); y2 = copy(y)
      ArpackInJulia.dsortr(:LA, true, length(x), x1, y1)
      arpack_dsortr(:LA, true, length(x), x2, y2)
      @test x1==x2
      @test y1==y2
    end

    @testset "dseigt" begin
      using LinearAlgebra
      n = 30
      lams = 2*cos.(pi.*(1:n)./(n+1))
      sort!(lams) # sort lams
      T1 = SymTridiagonal(zeros(n), ones(n-1))
      H = zeros(n, 2)
      H[2:end,1] .= T1.ev
      H[:,2] .= T1.dv
      H2 = copy(H)
      rnorm = 64.0
      # when we compare, we want zeros...
      eigval1 = randn(n)
      bounds1 = randn(n) 
      workl1 = randn(n,3)

      eigval2 = copy(eigval1)
      bounds2 = copy(bounds1)
      workl2 = copy(workl1)

      info1 = ArpackInJulia.dseigt!(rnorm, n, H, n, eigval1, bounds1, workl1, nothing)
      info2 = arpack_dseigt!(rnorm, n, H2, n, eigval2, bounds2, workl2)
      @test H == H2
      @test eigval1 == eigval2 
      @test bounds1 == bounds2 
      @test workl1 == workl2 
      @test info1 == info2
    end

    
    include("dsgets_arpackjll.jl")

    @testset "dsaitr" begin 
      include("dsaitr_arpackjll.jl")
    end 

    include("dsaupd_arpackjll.jl")

    include("dseupd_arpackjll.jl")
    
    @testset "dstqrb" begin
      include("dstqrb-compare.jl")
    end
    @testset "dgetv0" begin
      include("dgetv0_arpackjll.jl")
    end

    

    # The next set of tests delves deeply into
    # comparisons between Arpackjll and our code
    include("arpackjll_check_state.jl")
    include("dsapps_override.jl")

    
  end
end
##
function all_permutations(v::Vector{Float64};
    all::Vector{Vector{Float64}}=Vector{Vector{Float64}}(),
    prefix::Vector{Float64}=zeros(0))

  if length(v) == 0
    push!(all, prefix)
  end

  for i in eachindex(v)
    all_permutations(deleteat!(copy(v), i); prefix=push!(copy(prefix), v[i]), all)
  end
  return all
end
##
if "full" in ARGS
  include("dnrm2_openblas.jl")

  @testset "dsortr" begin
    for range in [0:0.0,0:1.0,-1:1.0,-2:1.0,-2:2.0,-2:3.0,-3:3.0,-4:3.0,-4:4.0,-4:5.0]
      X = all_permutations(collect(range))
      for x in X
        y = copy(x)
        @test_nowarn ArpackInJulia.dsortr(:LA, false, length(y), y, zeros(0))
        @test issorted(y)

        y = copy(x)
        @test_nowarn ArpackInJulia.dsortr(:LM, false, length(y), y, zeros(0))
        @test issorted(y, by=abs)

        y = copy(x)
        @test_nowarn ArpackInJulia.dsortr(:SM, false, length(y), y, zeros(0))
        @test issorted(y, by=abs, rev=true)

        y = copy(x)
        @test_nowarn ArpackInJulia.dsortr(:SA, false, length(y), y, zeros(0))
        @test issorted(y, rev=true)
      end
    end
  end
end
