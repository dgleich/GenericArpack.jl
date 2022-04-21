#= raw BLAS and LINPACK calls in ARPACK mapped to Julia code =#
using LinearAlgebra: givens

## include the better _drnm2 
include("arpack-blas-nrm2.jl")

##
function plane_rotation(f::T, g::T) where T
  g, r = givens(f, g, 1, 2)
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

function _dscal!(a::T, v) where T
  broadcast!(*, v, v, a)
end
##
@inline function _dlassq(TA::Type, v::AbstractVector{T}) where T
  scale::TA = zero(TA)
  ssq::TA = one(TA)
  for xi in v
    if !iszero(xi)
      absxi = TA(abs(xi))
      if scale < absxi
        #ssq = one(TA)+ssq*(scale/absxi)^2
        val = scale/absxi
        ssq = one(TA)+ssq*(val)*(val)
        scale = TA(absxi)
      else
        #ssq += (absxi/scale)^2
        val = (absxi/scale)
        ssq += val*val
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
    scale, ssq = _dlassq(T, v)
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

# this is really the function _dlapy3 in Fortran... 
function _dlapy2_julia(x::Complex{T}, y::T) where T
  _dlapy3_julia(real(x), imag(x), y)
end 

function _dlapy3_julia(x::T, y::T, z::T) where T
  xabs = abs(x)
  yabs = abs(y)
  zabs = abs(z)
  w = max(xabs, yabs, zabs)
  if w == 0 || w > floatmax(T)
    #  *     W can be zero for max(0,nan,0)
    # *     adding all three entries together will make sure
    # *     NaN will not disappear.
    return xabs + yabs + zabs
  else
    return w*sqrt( ( xabs / w )^2+( yabs / w )^2+ ( zabs / w )^2 )
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
    if ctemp != one(T) || stemp != zero(T) # there is work to do!
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
    if ctemp != one(T) || stemp != zero(T) # there is work to do!
      for i=1:m
        temp = a[i,j+1]
        a[i,j+1] = ctemp*temp - stemp*a[i,j]
        a[i,j] = stemp*temp + ctemp*a[i,j]
      end
    end
  end
  return a
end


##
const _dlaruv_mm=tuple(tuple(494,322,2508,2549),
    tuple(2637,789,3754,1145),
    tuple(255,1440,1766,2253),
    tuple(2008,752,3572,305),
    tuple(1253,2859,2893,3301),
    tuple(3344,123,307,1065),
    tuple(4084,1848,1297,3133),
    tuple(1739,643,3966,2913),
    tuple(3143,2405,758,3285),
    tuple(3468,2638,2598,1241),
    tuple(688,2344,3406,1197),
    tuple(1657,46,2922,3729),
    tuple(1238,3814,1038,2501),
    tuple(3166,913,2934,1673),
    tuple(1292,3649,2091,541),
    tuple(3422,339,2451,2753),
    tuple(1270,3808,1580,949),
    tuple(2016,822,1958,2361),
    tuple(154,2832,2055,1165),
    tuple(2862,3078,1507,4081),
    tuple(697,3633,1078,2725),
    tuple(1706,2970,3273,3305),
    tuple(491,637,17,3069),
    tuple(931,2249,854,3617),
    tuple(1444,2081,2916,3733),
    tuple(444,4019,3971,409),
    tuple(3577,1478,2889,2157),
    tuple(3944,242,3831,1361),
    tuple(2184,481,2621,3973),
    tuple(1661,2075,1541,1865),
    tuple(3482,4058,893,2525),
    tuple(657,622,736,1409),
    tuple(3023,3376,3992,3445),
    tuple(3618,812,787,3577),
    tuple(1267,234,2125,77),
    tuple(1828,641,2364,3761),
    tuple(164,4005,2460,2149),
    tuple(3798,1122,257,1449),
    tuple(3087,3135,1574,3005),
    tuple(2400,2640,3912,225),
    tuple(2870,2302,1216,85),
    tuple(3876,40,3248,3673),
    tuple(1905,1832,3401,3117),
    tuple(1593,2247,2124,3089),
    tuple(1797,2034,2762,1349),
    tuple(1234,2637,149,2057),
    tuple(3460,1287,2245,413),
    tuple(328,1691,166,65),
    tuple(2861,496,466,1845),
    tuple(1950,1597,4018,697),
    tuple(617,2394,1399,3085),
    tuple(2070,2584,190,3441),
    tuple(3331,1843,2879,1573),
    tuple(769,336,153,3689),
    tuple(1558,1472,2320,2941),
    tuple(2412,2407,18,929),
    tuple(2800,433,712,533),
    tuple(189,2096,2159,2841),
    tuple(287,1761,2318,4077),
    tuple(2045,2810,2091,721),
    tuple(1227,566,3443,2821),
    tuple(2838,442,1510,2249),
    tuple(209,41,449,2397),
    tuple(2770,1238,1956,2817),
    tuple(3654,1086,2201,245),
    tuple(3993,603,3137,1913),
    tuple(192,840,3399,1997),
    tuple(2253,3168,1321,3121),
    tuple(3491,1499,2271,997),
    tuple(2889,1084,3667,1833),
    tuple(2857,3438,2703,2877),
    tuple(2094,2408,629,1633),
    tuple(1818,1589,2365,981),
    tuple(688,2391,2431,2009),
    tuple(1407,288,1113,941),
    tuple(634,26,3922,2449),
    tuple(3231,512,2554,197),
    tuple(815,1456,184,2441),
    tuple(3524,171,2099,285),
    tuple(1914,1677,3228,1473),
    tuple(516,2657,4012,2741),
    tuple(164,2270,1921,3129),
    tuple(303,2587,3452,909),
    tuple(2144,2961,3901,2801),
    tuple(3480,1970,572,421),
    tuple(119,1817,3309,4073),
    tuple(3357,676,3171,2813),
    tuple(837,1410,817,2337),
    tuple(2826,3723,3039,1429),
    tuple(2332,2803,1696,1177),
    tuple(2089,3185,1256,1901),
    tuple(3780,184,3715,81),
    tuple(1700,663,2077,1669),
    tuple(3712,499,3019,2633),
    tuple(150,3784,1497,2269),
    tuple(2000,1631,1101,129),
    tuple(3375,1925,717,1141),
    tuple(1621,3912,51,249),
    tuple(3090,1398,981,3917),
    tuple(3765,1349,1978,2481),
    tuple(1149,1441,1813,3941),
    tuple(3146,2224,3881,2217),
    tuple(33,2411,76,2749),
    tuple(3082,1907,3846,3041),
    tuple(2741,3192,3694,1877),
    tuple(359,2786,1682,345),
    tuple(3316,382,124,2861),
    tuple(1749,37,1660,1809),
    tuple(185,759,3997,3141),
    tuple(2784,2948,479,2825),
    tuple(2202,1862,1141,157),
    tuple(2199,3802,886,2881),
    tuple(1364,2423,3514,3637),
    tuple(1244,2051,1301,1465),
    tuple(2020,2295,3604,2829),
    tuple(3160,1332,1888,2161),
    tuple(2785,1832,1836,3365),
    tuple(2772,2405,1990,361),
    tuple(1217,3638,2058,2685),
    tuple(1822,3661,692,3745),
    tuple(1245,327,1194,2325),
    tuple(2252,3660,20,3609),
    tuple(3904,716,3285,3821),
    tuple(2774,1842,2046,3537),
    tuple(997,3987,2107,517),
    tuple(2573,1368,3508,3017),
    tuple(1148,1848,3525,2141),
    tuple(545,2366,3801,1537))
function _dlaruv!(iseed::Base.RefValue{NTuple{4,Int}}, n::Int, x::Vector{Float64})
  IPW2=4096
  R = 1/IPW2
  i1 = iseed[][1]
  i2 = iseed[][2]
  i3 = iseed[][3]
  i4 = iseed[][4]
  # setup scope for final call in julia
  IT1 = i1
  IT2 = i2
  IT3 = i3
  IT4 = i4
  @inbounds for i=1:min(n,length(_dlaruv_mm))
    # directly copied/pasted from fortran, hence the caps
    while true
      #        Multiply the seed by i-th power of the multiplier modulo 2**48
      IT4 = i4*_dlaruv_mm[i][4]
      IT3 = IT4 ÷ IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + i3*_dlaruv_mm[i][4] + i4*_dlaruv_mm[i][3]
      IT2 = IT3 ÷ IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + i2*_dlaruv_mm[i][4] + i3*_dlaruv_mm[i][3] + i4*_dlaruv_mm[i][2]
      IT1 = IT2 ÷ IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + i1*_dlaruv_mm[i][4] + i2*_dlaruv_mm[i][3] + i3*_dlaruv_mm[i][2] + i4*_dlaruv_mm[i][1]
      IT1 = mod( IT1, IPW2 )

      #        Convert 48-bit integer to a real number in the interval (0,1)
      x[i] = R*( Float64( IT1 )+R*( Float64( IT2 )+R*( Float64( IT3 )+R*
                 Float64( IT4 ) ) ) )
      if x[i] == 1
        # handle the case where the value is different
        i1 += 2
        i2 += 2
        i3 += 2
        i4 += 2
      else
        break
      end
    end
  end
  iseed[] = tuple(IT1,IT2,IT3,IT4) # update the seed
end

function _dlarnv_idist_2!(iseed::Base.RefValue{NTuple{4,Int}}, n::Int, x::Vector{T}) where {
  T <: Union{AbstractFloat,Complex{<: AbstractFloat}}}

  iscomplex = eltype(x) <: Complex 
  if iscomplex
    RT = typeof(real(one(T)))
    maxloop = 2n 
  else
    RT = T 
    maxloop = n 
  end 

  # manually inline dlaruv to avoid excess cruft...
  IPW2=4096
  R = one(RT)/IPW2
  i1 = iseed[][1]
  i2 = iseed[][2]
  i3 = iseed[][3]
  i4 = iseed[][4]
  # setup scope for final call in julia
  IT1 = i1
  IT2 = i2
  IT3 = i3
  IT4 = i4

  nrv = 0

  

  # there is some subpar code here in terms of how 
  # the loop is handled for mixed complex and floats.
  # basically, we store the last rv to make sure 
  # that the complex code path can get two rvs
  # and we double the loop. 
  lastrv = zero(T)
  @inbounds for i=1:maxloop
    rv = 1.0
    nrv += 1
    if nrv > length(_dlaruv_mm)
      nrv = 1  # update seed information
      i1 = IT1
      i2 = IT2
      i3 = IT3
      i4 = IT4
    end

    while true
      #        Multiply the seed by i-th power of the multiplier modulo 2**48
      IT4 = i4*_dlaruv_mm[nrv][4]
      IT3 = IT4 ÷ IPW2
      IT4 = IT4 - IPW2*IT3
      IT3 = IT3 + i3*_dlaruv_mm[nrv][4] + i4*_dlaruv_mm[nrv][3]
      IT2 = IT3 ÷ IPW2
      IT3 = IT3 - IPW2*IT2
      IT2 = IT2 + i2*_dlaruv_mm[nrv][4] + i3*_dlaruv_mm[nrv][3] + i4*_dlaruv_mm[nrv][2]
      IT1 = IT2 ÷ IPW2
      IT2 = IT2 - IPW2*IT1
      IT1 = IT1 + i1*_dlaruv_mm[nrv][4] + i2*_dlaruv_mm[nrv][3] + i3*_dlaruv_mm[nrv][2] + i4*_dlaruv_mm[nrv][1]
      IT1 = mod( IT1, IPW2 )

      #        Convert 48-bit integer to a real number in the interval (0,1)
      rv = R*( RT( IT1 )+R*( RT( IT2 )+R*( RT( IT3 )+R*
               RT( IT4 ) ) ) )
      if rv == 1
        # handle the case where the value is different
        i1 += 2
        i2 += 2
        i3 += 2
        i4 += 2
      else
        break
      end
    end
    # now, rv is a random value
    if iscomplex
      # we index on 2, 4, 6... -> exact with i div 2
      ind,off = divrem(i,2)
      if off == 0 
        x[ind] = T(2*lastrv-1, 2*rv-1)
      end
    else
      x[i] = 2*rv - 1
    end 

    lastrv = rv 
  end
  iseed[] = tuple(IT1,IT2,IT3,IT4) # update the seed
end

function ddot(T::Type, a::AbstractVector, b::AbstractVector)
  rval = zero(T)
  length(a) == length(b) || throw(ArgumentError("ddot needs both a and b to have the same length"))
  @simd for i in eachindex(a)
    rval += T(conj(a[i])*b[i])
  end 
  return rval 
end 