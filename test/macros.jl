@testset "macros" begin
  @testset "docs" begin
    mutable struct MyTestType
      x::Int64
      y::Float64
    end
    @test """# Fields
    ```
    x :: Int64
    y :: Float64
    ```
    """ == ArpackInJulia._output_fields_in_markdown(MyTestType)
  end

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

  function test_arpack_debug(x;debug=nothing)
    msglvl = ArpackInJulia.@jl_arpack_debug(mgets,0)
    if msglvl > 0
      x[1] += 1
    else
      x[1] = 0
    end
  end

  debug=ArpackDebug()

  x = ones(1)
  test_arpack_debug(x)
  @test x[1] == 0
  x = ones(1)
  test_arpack_debug(x;debug)
  @test x[1] == 0
  debug.mgets = 1
  test_arpack_debug(x;debug)
  @test x[1] == 1

  function test_arpack_debug_output(x;debug=nothing)
    msglvl = ArpackInJulia.@jl_arpack_debug(mgets,0)
    if msglvl > 0
      println(debug.logfile, "xval=", x[1])
    end
  end

  io = IOBuffer()
  debug=ArpackDebug(logfile=io)
  x[1] = 0.5
  @test_nowarn test_arpack_debug_output(x)
  @test_nowarn test_arpack_debug_output(x;debug=nothing)
  @test_nowarn test_arpack_debug_output(x;debug)
  @test String(take!(io)) == ""
  debug.mgets = 1
  @test_nowarn test_arpack_debug_output(x;debug)
  @test String(take!(io)) == "xval=0.5\n"
end
