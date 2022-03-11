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