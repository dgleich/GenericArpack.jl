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
    @test ArpackInJulia.dsconv(5, ones(5), zeros(5), 1e-8; stats)==5
    @test stats.tconv > 0  # we did record some time
    t1 = stats.tconv
    @test ArpackInJulia.dsconv(10, ones(10), zeros(10), 1e-8; stats)==10
    @test stats.tconv > t1  # we did record some more time

    # TODO, add more relevant tests here.
  end
end
