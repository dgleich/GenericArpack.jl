
#=
The goal of this file is to implement
the blas dnrm2 function as closely as
possible to that which is in OpenBLAS.
This is very tricky becasue the OpenBLAS
version uses 80-bit extended precision.
So when ARPACK calls dnrm2, it uses
this 80-bit version. And the ARPACK/
Lanczos/Arnoldi codes are very sensitive
to the precise norm computation. So we 
want to get this as correct as possible. 

Our idea is to use doubledouble algs.
But this isn't quite good enough because
80-bit Intel stuff uses 15-bits for
exponent, so they don't have to deal
with scaling. We do :( 

So basically, we check if we need to scale.
Then scale by something good for Float64.
=# 

# Need a few ops from doublefloats.jl
# to get high-precision 

#=
MIT License

Copyright (c) 2018 Julia Math

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
=#

@inline function two_sum(a::T, b::T) where {T<:Float64}
  hi = a + b
  a1 = hi - b
  b1 = hi - a1
  lo = (a - a1) + (b - b1)
  return hi, lo
end

"""
    two_sum(a, b, c)
Computes `hi = fl(a+b+c)` and `lo = err(a+b+c)`.
"""
function two_sum(a::T,b::T,c::T) where {T<:Float64}
    t0, t1 = two_sum(a,  b)
    hi, t2 = two_sum(t0, c)
    lo = t2 + t1
    hi, lo = two_hilo_sum(hi, lo)
    return hi, lo
end

@inline function two_hilo_sum(a::T, b::T) where {T<:Float64}
  s = a + b
  e = b - (s - a)
  return s, e
end

# Algorithm 6 from Tight and rigourous error bounds. relative error < 3u²
@inline function add_dddd_dd(x::Tuple{T,T}, y::Tuple{T,T}) where T<:Float64
  xhi, xlo = x
  yhi, ylo = y
  hi, lo = two_sum(xhi, yhi)
  thi, tlo = two_sum(xlo, ylo)
  c = lo + thi
  hi, lo = two_hilo_sum(hi, c)
  c = tlo + lo
  hi, lo = two_hilo_sum(hi, c)
  return hi, lo
end

# Algorithm 12 from Tight and rigourous error bounds.  relative error <= 5u²
@inline function mul_dddd_dd(x::Tuple{T,T}, y::Tuple{T,T}) where T<:Float64
  xhi, xlo = x
  yhi, ylo = y
  hi, lo = two_prod(xhi, yhi)
  t = xlo * ylo
  t = fma(xhi, ylo, t)
  t = fma(xlo, yhi, t)
  t = lo + t
  hi, lo = two_hilo_sum(hi, t)
  return hi, lo
end

@inline function add_fpdd_dd(x::T, y::Tuple{T,T}) where {T<:Float64}
  yhi, ylo = y
  yhi, ylo = two_sum(x, yhi, ylo)
  return yhi, ylo
end

"""
    two_diff(a, b)
Computes `hi = fl(a-b)` and `lo = err(a-b)`.
"""
function two_diff(a::T, b::T) where {T<:Float64}
    hi = a - b
    a1 = hi + b
    b1 = hi - a1
    lo = (a - a1) - (b + b1)
    return hi, lo
end


"""
    two_diff(a, b, c)
Computes `hi = fl(a-b-c)` and `lo = err(a-b-c)`.
"""

function two_diff(a::T, b::T, c::T) where {T<:Float64}
    s, t = two_diff(-b, c)
    x, u = two_sum(a, s)
       y = u + t
    x, y = two_sum(x, y)
    return x, y
end

@inline function sub_fpdd_dd(x::T, y::Tuple{T,T}) where {T<:Float64}
  yhi, ylo = y
  yhi, ylo = two_diff(x, yhi, ylo)
  return yhi, ylo
end

@inline function two_prod(a::T, b::T) where {T<:Float64}
  s = a * b
  t = fma(a, b, -s)
  return s, t
end

@inline function add_sum_sq(a::Tuple{T,T}, b::T) where {T <: Float64}
  return add_dddd_dd(a, two_prod(b,b))
end

# This function is super annoying because 
# the OpenBLAS routines uses the 80-bit FP
# registers on x86_64, so here, we are
# going to use Kahan summation to simulation that...
# (using the 80-bit registers was a neat trick
# to )
function _dnrm2_unroll_ext(a::AbstractVector{T}) where {T <: Float64}
  max = maximum(abs, a)
  #min,max = extrema(abs, a)
  # here's the question, do we have enough precision in strict double-double,
  # or do we need to scale?
  if max >= sqrt(prevfloat(Inf))
    scale = (2*one(T))^(-511)
    #scale = one(T)
  else
    scale = one(T)
  end 
  len = length(a) 
  offset = 1
  ss1 = (zero(T), zero(T))
  ss2 = (zero(T), zero(T))
  ss3 = (zero(T), zero(T))
  ss4 = (zero(T), zero(T))
  @inbounds while offset+8 <= len
    ss1 = add_sum_sq(ss1, a[offset]*scale)
    ss2 = add_sum_sq(ss2, a[offset+1]*scale)
    ss3 = add_sum_sq(ss3, a[offset+2]*scale)
    ss4 = add_sum_sq(ss4, a[offset+3]*scale)
    ss1 = add_sum_sq(ss1, a[offset+4]*scale)
    ss2 = add_sum_sq(ss2, a[offset+5]*scale)
    ss3 = add_sum_sq(ss3, a[offset+6]*scale)
    ss4 = add_sum_sq(ss4, a[offset+7]*scale)
    offset += 8 
  end
  @inbounds for i=offset:len
    ss1 = add_sum_sq(ss1, a[i]*scale)
  end
  ss1 = add_dddd_dd(ss1, ss3)
  ss1 = add_dddd_dd(ss1, ss2) 
  ss1 = add_dddd_dd(ss1, ss4) 
  
  # lots of work to get the last bit right...
  # inline sqrt_dd_dd from DoubleFloats.jl
  r = inv(sqrt(ss1[1]))
  h = (ss1[1]*0.5, ss1[2]*0.5)
  r2 = two_prod(r, r) 
  hr2 = mul_dddd_dd(h, r2)
  radj = sub_fpdd_dd(0.5, hr2)
  radj = mul_dddd_dd(radj, (r, 0.0))
  r = add_fpdd_dd(r, radj)

  r2 = mul_dddd_dd(r,r)
  hr2 = mul_dddd_dd(h, r2)
  radj = sub_fpdd_dd(0.5, hr2)
  radj = mul_dddd_dd(radj, r)
  r = add_dddd_dd(r, radj)

  r = mul_dddd_dd(r, ss1)
  r = mul_dddd_dd(r, (1/scale, zero(T)))
  return r[1]
end  