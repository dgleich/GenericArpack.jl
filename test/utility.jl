using ArpackInJulia
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

function _allocate_symproblem(op, ncv::Int)
  n = size(op)
  @assert(ncv <= n)
  resid = zeros(n)
  V = zeros(n,ncv)
  iparam = zeros(Int,11)

  ipntr = zeros(Int,11)
  workd = zeros(3n)
  lworkl = ncv*ncv + 8*ncv
  workl = zeros(lworkl)

  ido = Ref{Int}(0)
  prob = ArpackSymProblem(op, V, resid, workd, workl, iparam, ipntr, ido)
  return _reset!(prob)
end 

function _eigrun!(prob,nev; which=:LM, bmat=Val(:I), mode=1, 
    iterfunc=nothing, state=nothing, ncv=size(prob.V,2), tol=0, initv=nothing,
    stats=nothing, debug=nothing, maxiter=300, idonow::Bool=false)

  V = prob.V
  _reset!(prob)
  prob.iparam[4] = maxiter 
  prob.iparam[7] = mode 

  @assert(nev < size(V,2))
  @assert(ncv <= size(V,2))
  T = eltype(V)
  n = size(V, 1)
  tol = tol*one(T)
  
  
  if state === nothing 
    state = ArpackInJulia.ArpackState{T}()
  end 
  ido = prob.ido

  lworkl = length(prob.workl)

  info_initv = 0
  if initv !==nothing
    info_initv = 1
    copyto!(prob.resid, initv)
  end 
  
  niter = 0 
  if idonow == false 
    op = prob.op
    while ido[] != 99
      ierr = ArpackInJulia.dsaupd!(ido, bmat, n, which, nev, tol, prob.resid, ncv, V, size(V,1), 
        prob.iparam,
        prob.ipntr, prob.workd, prob.workl, lworkl, info_initv;
        state, stats, debug 
      ).ierr 
      niter += 1
      if iterfunc !== nothing
        iterfunc(prob, ierr, state)
      end 

      if ido[] == 1 
        ArpackInJulia._i_do_now_opx_1!(op, prob.ipntr, prob.workd, n)
      elseif ido[] == -1 
        ArpackInJulia._i_do_now_opx_neg1!(op, prob.ipntr, prob.workd, n)
      elseif ido[] == 2
        ArpackInJulia._i_do_now_bx!(op, prob.ipntr, prob.workd, n)
      end 
    end
    return niter 
  else
    ierr = ArpackInJulia.dsaupd!(ido, bmat, n, which, nev, tol, prob.resid, ncv, V, size(V,1), 
        prob.iparam,
        prob.ipntr, prob.workd, prob.workl, lworkl, info_initv;
        state, stats, debug, idonow= prob.op 
      ).ierr

    # with idonow, everything should be handled in-place! 
    @assert ido[] == 99 
    return niter+1
  end 
end 

