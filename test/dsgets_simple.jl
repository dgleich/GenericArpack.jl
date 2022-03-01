@testset "dsgets" begin
  ishift = 1
  for which = [:BE,:LA,:LM,:SM,:SA]
    kev = 3
    np = 2
    ritz = [-1.0, 0.05, -0.001, 1.0, 2.0]
    bounds = [0.0, 0.01, 0.00001, 2.0, 0.5]
    pairs = Set(zip(ritz, bounds))
    shifts = [0.0, 0.0, 0.0]
    @test_nowarn ArpackInJulia.dsgets(ishift, which, kev, np, ritz, bounds, shifts)
    @test pairs == Set(zip(ritz, bounds)) # should get the same pairs in both cases

    kev = 2
    np = 3
    ritz = [-1.5, 0.05, -0.001, 1.0, 2.0]
    bounds = [0.0, 0.01, 0.00001, 2.0, 0.5]
    shifts = [0.0, 0.0, 0.0]
    pairs = Set(zip(ritz, bounds))
    @test_nowarn ArpackInJulia.dsgets(ishift, which, kev, np, ritz, bounds, shifts)
    @test pairs == Set(zip(ritz, bounds)) # should get the same pairs in both cases

    @test union!(Set(shifts[1:np]), ritz[np+1:end]) == Set(ritz)
    if which == :LM
      @test ritz[end] == 2.0
      @test ritz[end-1] == -1.5
    elseif which == :SA
      @test ritz[end] == -1.5 # smallest algebraic entry
      @test ritz[end-1] == -0.001
    elseif which == :SM
      @test ritz[end] == -0.001 # smallest magnitude is last
      @test ritz[end-1] == 0.05
    elseif which == :BE
      @test ritz[end] == 2.0 # largest magnitude
      @test ritz[end-1] == -1.5 # largest at other end
    elseif which == :LA
      @test ritz[end] == 2.0 # largest magnitude
      @test ritz[end-1] == 1.0 # largest at other end
    end
  end

  ishift = 0
  # in this case shifts should not really be changed.
end