using LinearAlgebra
using SparseArrays

@testset "BMAT parameter check (bmat vs BMAT bug)" begin
  # Bug: line 1070 used `bmat` (a module-level function) instead of `BMAT`
  # (the Val type parameter). Since `bmat` is a function, `bmat == :G` is
  # always false, so the error check for mode=1 + BMAT=:G was silently skipped.
  n = 10
  A = spdiagm(0 => 2*ones(n), -1 => -ones(n-1), 1 => -ones(n-1))
  info = eigs(Symmetric(A), 3)

  ncv = size(info.V, 2)
  nconv = info.iparam[5]
  select = zeros(Int, ncv)

  # mode=1 with BMAT=:G is invalid and should return ierr=-11
  ierr = GenericArpack.simple_dseupd!(
    true, select, copy(info.values), Matrix{Float64}(undef, n, nconv), 0.0,
    Val(:G), n, info.which, 3, eps(Float64)/2, copy(info.resid),
    ncv, copy(info.V), copy(info.iparam), copy(info.ipntr),
    copy(info.workd), copy(info.workl))
  @test ierr == -11
end

@testset "bnorm2 typo (bnrom2 bug)" begin
  # Bug: line 1402 used `bnrom2` (typo) instead of `bnorm2`, causing
  # UndefVarError in the shift-invert mode (mode=3) with rvec=false path.
  n = 10
  A = spdiagm(0 => 2*ones(n), -1 => -ones(n-1), 1 => -ones(n-1))
  info = eigs(Symmetric(A), 3)

  ncv = size(info.V, 2)
  select = zeros(Int, ncv)
  iparam = copy(info.iparam)
  iparam[7] = 3  # force mode=3 (shift-invert) to trigger the buggy path

  # rvec=false triggers the bnorm2 path; should not crash
  ierr = GenericArpack.simple_dseupd!(
    false, select, copy(info.values), Matrix{Float64}(undef, n, 0), 0.0,
    Val(:I), n, info.which, 3, eps(Float64)/2, copy(info.resid),
    ncv, copy(info.V), iparam, copy(info.ipntr),
    copy(info.workd), copy(info.workl))
  @test ierr == 0
end
