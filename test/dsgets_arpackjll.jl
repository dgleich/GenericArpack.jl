@testset "dsgets" begin
  ishift = 1
  for which = [:BE,:LA,:LM,:SM,:SA]
    kev = 3
    np = 2
    ritz = [-1.0, 0.05, -0.001, 1.0, 2.0]
    bounds = [0.0, 0.01, 0.00001, 2.0, 0.5]
    shifts = [0.0, 0.0, 0.0]

    aritz = copy(ritz)
    abounds = copy(bounds)
    ashifts = copy(shifts)
    arpack_dsgets(ishift, which, kev, np, aritz, abounds, ashifts)
    GenericArpack.dsgets(ishift, which, kev, np, ritz, bounds, shifts)
    @test ritz == aritz
    @test bounds == abounds
    @test shifts == ashifts

    kev = 2
    np = 3
    ritz = [-1.0, 0.05, -0.001, 1.0, 2.0]
    bounds = [0.0, 0.01, 0.00001, 2.0, 0.5]
    shifts = [0.0, 0.0, 0.0]

    aritz = copy(ritz)
    abounds = copy(bounds)
    ashifts = copy(shifts)
    arpack_dsgets(ishift, which, kev, np, aritz, abounds, ashifts)
    GenericArpack.dsgets(ishift, which, kev, np, ritz, bounds, shifts)
    @test ritz == aritz
    @test bounds == abounds
    @test shifts == ashifts

    kev = 5
    np = 0
    ritz = [-1.0, 0.05, -0.001, 1.0, 2.0]
    bounds = [0.0, 0.01, 0.00001, 2.0, 0.5]
    shifts = [0.0, 0.0, 0.0]

    aritz = copy(ritz)
    abounds = copy(bounds)
    ashifts = copy(shifts)
    arpack_dsgets(ishift, which, kev, np, aritz, abounds, ashifts)
    GenericArpack.dsgets(ishift, which, kev, np, ritz, bounds, shifts)
    @test ritz == aritz
    @test bounds == abounds
    @test shifts == ashifts
  end

  ishift=0
  for which = [:BE,:LA,:LM,:SM,:SA]
    kev = 2
    np = 3
    ritz = [-1.0, 0.05, -0.001, 1.0, 2.0]
    bounds = [0.0, 0.01, 0.00001, 2.0, 0.5]
    shifts = randn(3)
    orig_shifts = copy(shifts)

    aritz = copy(ritz)
    abounds = copy(bounds)
    ashifts = copy(shifts)
    arpack_dsgets(ishift, which, kev, np, aritz, abounds, ashifts)
    GenericArpack.dsgets(ishift, which, kev, np, ritz, bounds, shifts)
    @test ritz == aritz
    @test bounds == abounds
    @test shifts == ashifts == orig_shifts # should be unchanged for ishift=0
  end
end