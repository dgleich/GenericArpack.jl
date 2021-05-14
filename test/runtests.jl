using Test
using ArpackInJulia

@testset "macros" begin

  function test_arpack_macros(;stats=nothing)
    t0 = ArpackInJulia.@jl_arpack_time
    sleep(5.0)
    ArpackInJulia.@jl_update_time(taupd, t0)
    return 5
  end
  @test test_arpack_macros()==5 # this checks with no timing information

  arpackstat = ArpackStats()
  test_arpack_macros(stats=arpackstat)
  @test arpackstat.taupd > 0

end


@testset "simple" begin

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
end

if "arpackjll" in ARGS
  include("arpackjll.jl") # get the helpful routines to call stuff in "arpackjll"
  @testset "arpackjll" begin
    @testset "dsconv" begin
      soln = arpack_dsconv(10, ones(10), zeros(10), 1e-8)
      @test ArpackInJulia.dsconv(10, ones(10), zeros(10), 1e-8)==soln

      soln = arpack_dsconv(10, zeros(10), ones(10), 1e-8)
      @test ArpackInJulia.dsconv(10, zeros(10), ones(10), 1e-8)==soln

    end
  end
end
