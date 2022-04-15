using ArpackInJulia

function test_diffs(tagstr, a, b)
  @assert(size(a) == size(b)) # otherwise, use view to make them the same...
  first = false 
  for i in eachindex(a)
    if a[i] != b[i]
      if first == false
        first = true
        println("Difference between $tagstr")
      end 
      println(" ", rpad(a[i], 30), " , ", rpad(b[i], 30), " # diff=", relfloatsbetween(a[i], b[i]), " index=", i)
    end 
  end
end 

""" Return the number of floating point values between two values, including
one of the endpoints. """
function floatsbetween(a::T,b::T) where {T <: Float64}
  # based on https://github.com/JuliaLang/julia/commit/79920db3ac9f1e1915d836c82bc401ba0cdc322f
  # https://stackoverflow.com/questions/3587880/number-of-floats-between-two-floats
  IntT = Int64 # inttype(T)
  ia = reinterpret(IntT, a)
  ib = reinterpret(IntT, b)
  ia = ifelse(ia < zero(IntT), ia ⊻ typemax(IntT), ia)
  ib = ifelse(ib < zero(IntT), ib ⊻ typemax(IntT), ib)
  return ib-ia
end

""" Count the number of floats between but treat the range [-ref,ref] as two
fractional floating point values,

relfloatsbetween(-ref,ref;ref) = 2.0
relfloatsbetween(-prevfloat(ref),nextfloat(ref);ref) = 4.0
relfloatsbetween(ref/4,ref/2;ref) = 0.25
"""
function relfloatsbetween(a::T,b::T;ref=eps(T)/2) where {T <: Float64}
  # this is awful code at the moment :(
  if b < a
    flipsign = -1
    a,b = b,a # swap so we are in order to minimize cases below
  else
    flipsign = 1
  end

  rval = 0.0
  if a < -ref && b < -ref # then since a < b, we can just return this directly
    rval += floatsbetween(a,b)
  elseif a <= -ref && b <= ref
    rval += floatsbetween(a, -ref)
    rval += (b-(-ref))/ref
  elseif a <= -ref && b > ref
    rval += floatsbetween(a, -ref)
    rval += floatsbetween(ref, b)
    rval += 2.0
  elseif a <= ref && b <= ref
    # a in [-ref,ref], b in [-ref,ref]
    rval += (b-a)/ref
  elseif a <= ref && b > ref
    rval += (ref-a)/ref
    rval += floatsbetween(ref, b)
  else
    # a >= ref, b >= ref
    rval += floatsbetween(a,b)
  end

  return flipsign*rval
end

struct ArpackSymProblem{Op <: ArpackInJulia.ArpackOp,T}
  op::Op
  V::Matrix{T}
  resid::Vector{T}
  workd::Vector{T}
  workl::Vector{T}
  iparam::Vector{Int}
  ipntr::Vector{Int}
  ido::Base.RefValue{Int} 
end 

function _compare_probs(prob1, prob2; compare_ido=true)
  @test prob1.V == prob2.V
  @test prob1.resid == prob2.resid
  @test prob1.workd == prob2.workd
  @test prob1.workl == prob2.workl
  @test prob1.iparam == prob2.iparam
  @test prob1.ipntr == prob2.ipntr
  if compare_ido
    @test prob2.ido[] == prob1.ido[] 
  end 
end 

function _reset!(prob) 
  fill!(prob.iparam, 0)
  prob.iparam[1] = 1
  prob.iparam[4] = 1
  fill!(prob.V, 0)
  fill!(prob.resid, 0)
  fill!(prob.workd, 0)
  fill!(prob.workl, 0)
  fill!(prob.ipntr, 0)
  prob.ido[] = 0 
  return prob 
end 

function _allocate_symproblem(TV, TF, op, ncv::Int)
  n = size(op)
  @assert(ncv <= n)
  resid = zeros(TV, n)
  V = zeros(TV, n,ncv)
  iparam = zeros(Int,11)

  ipntr = zeros(Int,11)
  workd = zeros(TV, 3n)
  lworkl = ncv*ncv + 8*ncv
  workl = zeros(TF, lworkl)

  ido = Ref{Int}(0)
  prob = ArpackSymProblem(op, V, resid, workd, workl, iparam, ipntr, ido)
  return _reset!(prob)
end 

_allocate_symproblem(op, ncv::Int) = _allocate_symproblem(Float64, Float64, op, ncv)
  

function __val2sym(::Val{BMAT}) where BMAT
  return BMAT
end 

function _eigrun!(prob,nev; which=:LM, bmat=ArpackInJulia.bmat(prob.op), 
    mode=ArpackInJulia.arpack_mode(prob.op), 
    iterfunc=nothing, state=nothing, ncv=size(prob.V,2), tol=0, initv=nothing,
    stats=nothing, debug=nothing, maxiter=300, idonow::Bool=false, 
    arpackjllfunc = nothing)

  V = prob.V
  _reset!(prob)
  prob.iparam[3] = maxiter 
  prob.iparam[7] = mode 

  @assert(nev < size(V,2))
  @assert(ncv <= size(V,2))
  TF = eltype(prob.workl)
  n = size(V, 1)
  tol = tol*one(TF)
  
  
  if state === nothing 
    state = ArpackInJulia.ArpackState{TF}()
  end 
  ido = prob.ido

  lworkl = length(prob.workl)

  info_initv = 0
  if initv !== nothing
    info_initv = 1
    copyto!(prob.resid, initv)
  end 

  ierr = 0 
  niter = 0 
  if arpackjllfunc !== nothing
    @assert(idonow == false)
  end

  if idonow == false 
    op = prob.op
    aupdfunc = () -> ArpackInJulia.dsaupd!(ido, bmat, n, which, nev, tol, prob.resid, ncv, V, size(V,1), 
        prob.iparam,
        prob.ipntr, prob.workd, prob.workl, lworkl, info_initv;
        state, stats, debug 
      ).ierr 
    if arpackjllfunc !== nothing
      assert(TF == Float64)
      aupdfunc = () -> arpackjllfunc(ido, __val2sym(bmat), n, which, nev, tol, prob.resid, ncv, V, size(V,1),
        prob.iparam, prob.ipntr, prob.workd, prob.workl, lworkl, info_initv )
    end 

    while ido[] != 99
      ierr = aupdfunc()
      niter += 1
      if iterfunc !== nothing
        iterfunc(prob, ierr, state)
      end 
      info_initv = 0 

      if ido[] == 1 && mode == 1
        ArpackInJulia._i_do_now_opx_1!(op, prob.ipntr, prob.workd, n)
      elseif ido[] == 1 && mode == 2
        ArpackInJulia._i_do_now_opx_mode2_1!(op, prob.ipntr, prob.workd, n)
      elseif ido[] == 1 # mode 3, 4, 5
        ArpackInJulia._i_do_now_opx_shiftinvert_1!(op, prob.ipntr, prob.workd, n)
      elseif ido[] == -1 
        ArpackInJulia._i_do_now_opx_neg1!(op, prob.ipntr, prob.workd, n)
      elseif ido[] == 2
        ArpackInJulia._i_do_now_bx!(op, prob.ipntr, prob.workd, n)
      elseif ido[] == 3
        ArpackInJulia._i_do_now_shifts!(op, prob.iparam[8], prob.ipntr, prob.workl, n)
      end 
    end
    return niter, ierr  
  else
    ierr = ArpackInJulia.dsaupd!(ido, bmat, n, which, nev, tol, prob.resid, ncv, V, size(V,1), 
        prob.iparam,
        prob.ipntr, prob.workd, prob.workl, lworkl, info_initv;
        state, stats, debug, idonow= prob.op 
      ).ierr

    # with idonow, everything should be handled in-place! 
    @assert ido[] == 99 
    return niter+1, ierr
  end 
end 

function _eigenvecs(prob,nev;which=:LM, bmat=ArpackInJulia.bmat(prob.op), 
  iterfunc=nothing, state=nothing, ncv=size(prob.V,2), tol=0, 
  #stats=nothing, 
  debug=nothing)

  # make sure we are done! 
  @assert(prob.ido[] == 99)

  n = size(prob.V, 1)
  @assert(ncv > nev)
  Z = zeros(n, nev)
  d = zeros(nev)
  select = zeros(Int, ncv)

  rvec = true
  sigma = shift(eltype(prob.V), prob.op)

  ierr = ArpackInJulia.dseupd!(rvec, select, d, Z, sigma, 
    bmat, n, which, nev, tol, prob.resid, ncv, prob.V, prob.iparam, prob.ipntr, 
    prob.workd, prob.workl; debug)

  if ierr != 0
    throw(ErrorException("dseupd returns $ierr"))
  end 

  return d, V
end 

function check_svd(A, U, s, V; tol=2)
  @test U'*U ≈ I 
  @test V'*V ≈ I 
  r = ArpackInJulia.svd_residuals(A, U, s, V)
  if tol <= 0
    @test_broken maximum(r) <= eps(real(eltype(U)))*abs(tol)
  else 
    @test maximum(r) <= eps(real(eltype(U)))*tol
  end 
end 
##
# use something inspired by the lauchli matrix from the test matrix toolkit from Higham
# http://www.ma.man.ac.uk/~higham/mctoolbox/toolbox.pdf
# this has eigenvalues
# I'm sure we can work out a closed form solution for the singular values here...
function mytestmat(m,n,minval=1/m, maxdiagval=2)
  # the matrix is
  # [v, B]
  # where B is n-1 cols, m rows... 
  B = sparse(1:n-1, 1:n-1, range(1, maxdiagval, length=n-1), m, n-1)
  return [collect(range(1, minval, length=m)) B]
end 