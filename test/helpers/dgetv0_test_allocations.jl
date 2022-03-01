## test allocations
using BenchmarkTools
using ArpackInJulia
using LinearAlgebra
begin
  @btime begin
    ierr = ArpackInJulia.dgetv0!(
      ido, Val{bmat}, itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state, stats)
    if ido[] == -1
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    end
    ierr = ArpackInJulia.dgetv0!(
      ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state, stats)
    if ido[]==2
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),B,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    end
    ierr = ArpackInJulia.dgetv0!(
      ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state, stats)
  end setup=begin
    ido = Ref{Int}(0)
    bmat = Val{:G}
    itry = 1
    initv = true
    n = 10
    j = 1
    v = zeros(n, j+1)
    ldv = n
    resid = 2*rand(n).-1
    rnorm = Ref{Float64}(0.0)
    ipntr = zeros(Int, 3)
    workd = zeros(2n)
    M = Diagonal(1.0:10.0)
    B = Diagonal(0.1:0.1:1.0)
    stats = ArpackStats()
    state=ArpackInJulia.ArpackState{Float64}()
  end
end

##
using ArpackInJulia
using LinearAlgebra
begin
  @btime begin
    ierr = ArpackInJulia.dgetv0!(
      ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state, stats)
    if ido[] == -1
      mul!(@view(workd[ipntr[2]:ipntr[2]+n-1]),M,@view(workd[ipntr[1]:ipntr[1]+n-1]))
    end
    ierr = ArpackInJulia.dgetv0!(
      ido, bmat, itry, initv, n, j, v, ldv, resid, rnorm, ipntr, workd;
      state, stats)
  end setup=begin
    ido = Ref{Int}(0)
    bmat = Val{:I}
    itry = 1
    initv = true
    n = 10
    j = 1
    v = zeros(n, j+1)
    ldv = n
    resid = 2*rand(n).-1
    rnorm = Ref{Float64}(0.0)
    ipntr = zeros(Int, 3)
    workd = zeros(2n)
    M = Diagonal(1.0:10.0)
    B = Diagonal(0.1:0.1:1.0)
    stats = ArpackStats()
    state=ArpackInJulia.ArpackState{Float64}()
  end
end
