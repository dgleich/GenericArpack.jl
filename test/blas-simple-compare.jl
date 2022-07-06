@testset "_dlaev2" begin 
  include("../src/arpack-blas-direct.jl")
  function _checked_dlaev2(a::Float64, b::Float64, c::Float64)
    @testset "a=$a, b=$v, c=$c" begin 
      @test GenericArpack._dlaev2(a,b,c)==_dlaev2_blas(a,b,c)
    end 
  end 

  _checked_dlaev2(1.0, 0.0, 2.0)
  _checked_dlaev2(2.0, 0.0, 1.0)
  _checked_dlaev2(1.0, 0.0, -2.0)
  _checked_dlaev2(-1.0, 0.0, -2.0)
  _checked_dlaev2(1.0, 0.0, -2.0)
end 