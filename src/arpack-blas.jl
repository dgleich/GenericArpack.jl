#= raw BLAS and LINPACK calls in ARPACK mapped to Julia code =#
using LinearAlgebra: givens
##
function plane_rotation(f::T, g::T) where T
  g, r = LinearAlgebra.givens(f, g, 1, 2)
  return g.c, g.s, r
end


##
""" Compute a scaling factor that is conceptually to/from,
but handles overflow/underflow. use in the following fashion

mul, cfromc, ctoc, done = _dlascal_factor(cfrom, cto, done)
mul!(v, mul)
while !done
  mul, cfromc, ctoc, done = _dlascal_factor(cfromc, ctoc, done)
  mul!(v, mul)
end
"""
function _dlascl_factor(cfrom::T, cto::T, done::Bool=false) where T
  # code from dlascl.f in LAPACK
  smlnum = floatmin(T)
  bignum = 1 / smlnum
  cfromc = cfrom
  ctoc = cto

  cfrom1 = cfromc*smlnum
  if (cfrom1 == cfromc)
    # cfromc is an inf... this will give a signed zero
    mul = ctoc / cfromc
    cto1 = ctoc
    done = true
  else
    cto1 = ctoc/bignum
    if cto1 == ctoc
      # ctoc is either 0 or an inf, so we get the correct factor
      mul = ctoc
      done = true
      cfromc = one(T)
    elseif (abs(cfrom1) > abs(ctoc) && ctoc != 0)
      mul = smlnum # multiply evrything by mul and restart
      done = false
      cfromc = cfrom1
    elseif (abs(cto1) > abs(cfromc))
      mul = bignum
      done = false
      ctoc = cto1
    else
      mul = ctoc / cfromc
      done = true
    end
  end
  return mul, cfromc, ctoc, done
end

function _scale_from_to(from::T, to::T, v) where T
  if from==to # test easy out
    return v
  end
  mul, cfromc, ctoc, done = _dlascl_factor(from, to)
  broadcast!(*, v, v, mul) # compute v .*= mul
  while !done
    mul, cfromc, ctoc, done = _dlascl_factor(cfromc, ctoc, done)
    broadcast!(*, v, v, mul)
  end
  return v
end

function _dscal(a::T, v) where T
  broadcast!(*, v, v, mul)
end
##
function _dlassq(v::AbstractVector{T}) where T
  scale::T = zero(T)
  ssq::T = one(T)
  for xi in v
    if !iszero(xi)
      absxi = abs(xi)
      if scale < absxi
        ssq = 1+ssq*(scale/absxi)^2
        scale = absxi
      else
        ssq += (absxi/scale)^2
      end
    end
  end
  return scale, ssq
end

function _dnrm2(v::AbstractVector{T}) where T
  n = length(v)
  norm::T = zero(T)
  if n == 1
    norm = abs(v[1])
  else
    scale, ssq = _dlassq(v)
    norm = scale*sqrt(ssq)
  end
  return norm::T
end

##
""" 2x2 eigenvalue problem
given
[a b
 b c]
compute the eigenvalue decomp
[ c s]*[a b]*[c -s] = [lam1 0    ]
[-s c] [b c] [s  c]   [0    lam2 ]
and return
(lam1,lam2,c,s)

Based on the LAPACK code
  """
function _dlaev2(a::T, b::T, c::T) where T
  sm = a + c
  df = a - c
  adf = abs( df )
  tb = b + b
  ab = abs( tb )
  half = T(1/2) # 1/2 is fully accurate in Float64
  if abs(a) > abs(c)
    acmx = a
    acmn = c
  else
    acmx = c
    acmn = a
  end
  if adf > ab
    rt = adf*sqrt(1+(ab/adf)^2)
  elseif adf < ab
    rt = ab*sqrt(1 + (adf/ab)^2)
  else
    rt = ab*sqrt(T(2))
  end

  if sm < 0
    rt1 = half*(sm-rt)
    sgn1 = -1
    rt2 = (acmx/rt1)*acmn - (b/rt1)*b
  elseif sm > 0
    rt1 = half*(sm+rt)
    sgn1 = 1
    rt2 = (acmx/rt1)*acmn - (b/rt1)*b
  else
    rt1 = half*rt
    rt2 = -half*rt
    sgn1 = 1
  end

  # compute the eigenvector
  if df >= 0
    cs = df + rt
    sgn2 = 1
  else
    cs = df - rt
    sgn2 = -1
  end
  acs = abs(cs)
  if (acs>ab)
    ct = -tb/cs
    sn1 = 1/sqrt(1+ct*ct)
    cs1 = ct*sn1
  else
    if ab == 0
      cs1::T = 1
      sn1::T = 0
    else
      tn = -cs/tb
      cs1 = 1/sqrt(1+tn*tn)
      sn1 = tn*cs1
    end
  end
  if sgn1 == sgn2
    tn = cs1
    cs1 = -sn1
    sn1 = tn
  end
  return rt1, rt2, cs1, sn1
end

##
function _dlapy2_julia(x::T, y::T) where T
  !isnan(x) || return x
  !isnan(y) || return y
  w = max(abs(x),abs(y))
  z = min(abs(x),abs(y))
  if z == 0
    return w
  else
    return w*sqrt(1+(z/w)^2)
  end
end

##
function _apply_plane_rotations_right!(a::AbstractVecOrMat{T},
    c::AbstractVector{T}, s::AbstractVector{T};
    n=size(a,2), m=size(a,1), order=1:n-1, rev::Bool=false) where T
  if rev
    _apply_plane_rotations_right!(a, c, s; n, m, rev=false, order=reverse(order))
  else
    @inbounds for j in order # forward or reverse order
      ctemp = c[j]
      stemp = s[j]
      if ctemp != one(T) || stemp != zero(T) # there is work to do!
        for i=1:m
          temp = a[i,j+1]
          a[i,j+1] = ctemp*temp - stemp*a[i,j]
          a[i,j] = stemp*temp + ctemp*a[i,j]
        end
      end
    end
  end
  return a
end

##
function _dlasr_right_side_variable_pivot_backward!(m::Integer, n::Integer,
    c::AbstractVector{T}, s::AbstractVector{T}, a::AbstractVecOrMat{T}, ) where T
#=
    ELSE IF( lsame( direct, 'B' ) ) THEN
362                DO 160 j = n - 1, 1, -1
363                   ctemp = c( j )
364                   stemp = s( j )
365                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
366                      DO 150 i = 1, m
367                         temp = a( i, j+1 )
368                         a( i, j+1 ) = ctemp*temp - stemp*a( i, j )
369                         a( i, j ) = stemp*temp + ctemp*a( i, j )
370   150                CONTINUE
371                   END IF
372   160          CONTINUE
=#

  @inbounds for j=n-1:-1:1
    ctemp = c[j]
    stemp = s[j]
    if ctemp != one(T) || stemp != one(T) # there is work to do!
      for i=1:m
        temp = a[i,j+1]
        a[i,j+1] = ctemp*temp - stemp*a[i,j]
        a[i,j] = stemp*temp + ctemp*a[i,j]
      end
    end
  end
  return a
end



function _dlasr_right_side_variable_pivot_forward!(m::Integer, n::Integer,
    c::AbstractVector{T}, s::AbstractVector{T}, a::AbstractVecOrMat{T}) where T

#= from dlasr.f
    IF( lsame( direct, 'F' ) ) THEN
      350                DO 140 j = 1, n - 1
      351                   ctemp = c( j )
      352                   stemp = s( j )
      353                   IF( ( ctemp.NE.one ) .OR. ( stemp.NE.zero ) ) THEN
      354                      DO 130 i = 1, m
      355                         temp = a( i, j+1 )
      356                         a( i, j+1 ) = ctemp*temp - stemp*a( i, j )
      357                         a( i, j ) = stemp*temp + ctemp*a( i, j )
      358   130                CONTINUE
      359                   END IF
=#

  @inbounds for j=1:n-1
    ctemp = c[j]
    stemp = s[j]
    if ctemp != one(T) || stemp != one(T) # there is work to do!
      @simd for i=1:m
        temp = a[i,j+1]
        a[i,j+1] = ctemp*temp - stemp*a[i,j]
        a[i,j] = stemp*temp - ctemp*a[i,j]
      end
    end
  end
  return a
end
