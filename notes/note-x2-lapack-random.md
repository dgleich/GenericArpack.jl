LAPACK's random number generators in Julia (part ... blah ... in work to port ARPACK to Julia.)

ARPACK depends on LAPACK's internal random number generators

  - `dlarnv`
  - `slarnv`
  - `zlarnv`
  - `clarnv`

which use the raw routines   

  - `dlaruv`
  - `slaruv`

This is true even if you use Julia to initialize the starting
vector using Julia's random number generator. The reason is because
ARPACK needs to be able to get random data to restart the
process is it detects an invariant subspace. This suggests two
approaches: call LAPACK's RNG from Julia, or port LAPACK's RNG
to Julia. Why do the simple thing when the more complicated thing
is so much more fun?

> Aside, best guess on LAPACK names: _larnv = lapack random normal vector and _laruv = lapack random uniform vector 

# Calling LAPACK's random number generator from julia

The tricky bit about calling LAPACK's RNG from Julia is
the seed value. This is just 4 ints, which we'd like to be
stack allocated. Unfortunately, we won't be able to get that
for reasons that are internal to Julia. (Blah blah,
stack vs. heap and gc and guarantees...)

Here's some code that works though.

    import LinearAlgebra

    _dlarnv_blas!(idist::Int,
           iseed::Ref{NTuple{4,Int}},
          n::Int,
          x::Vector{Float64}) =
      ccall((LinearAlgebra.BLAS.@blasfunc("dlarnv_"), LinearAlgebra.BLAS.libblas), Cvoid,
        (Ref{LinearAlgebra.BlasInt}, Ptr{LinearAlgebra.BlasInt}, Ref{LinearAlgebra.BlasInt}, Ptr{Float64}),
        idist, iseed, n, x)

    x = zeros(5); _dlarnv!(1,Base.Ref(tuple(1,2,3,4)),length(x),x)

    using BenchmarkTools
    @btime _dlarnv_blas!(1, tseed, length(z), z) setup=begin
        z=zeros(128);     
        tseed=Base.Ref(tuple(1,3,5,7));
    end

This uses `NTuple` as the holder for the 4 ints we use as the seed. At this
point, I thought the following should work...

    mutable struct MySeedState
      seed::NTuple{4,Int}
    end

    function myrand!(x::Vector{Float64}, state::MySeedState)
      iseed = Base.Ref(state.seed) # this will allocate :(

      _dlarnv_blas!(2, iseed, length(x), x)
      state.seed = iseed[]
    end

    using BenchmarkTools
    @btime myrand!(z, tstate) setup=begin z=zeros(128); tstate=MySeedState(tuple(1,3,5,7)); end

But -- see above -- I learned when we want to use this is that allocating
a Ref can actually allocate in Julia, even when it's an itsbits type.

    julia> isbitstype(NTuple{4,Int})
    true

To avoid this, we can just allocate the ref as part of the mutable struct.

    mutable struct MyRefSeedState
      seed::Ref{NTuple{4,Int}}
    end

    function myrand!(x::Vector{Float64}, state::MyRefSeedState)
      _dlarnv_blas!(2, state.seed, length(x), x)
    end

    using BenchmarkTools
    @btime myrand!(z, tstate) setup=begin
        z=zeros(128);
        tstate=MyRefSeedState(Base.Ref(tuple(1,3,5,7)));
    end

Either way, it probably isn't that big of a deal, but I'd like
to be clear on where the allocations are happening. e.g. one-time
thread local setup is okay. Random internal allocations are not
okay.

It's puzzling to me why this causes an allocation as this
seems like it can be safely done on the stack and the
same allocation doesn't occur for types like Int.

# Porting LAPACK's random number generator to Julia.

Since the routine is sufficiently simple, we can actually just port
the random number generator directly.

> Aside -- that this is actually
a weakness as Julia's RNGs are much better. On the other hand,
these are just used to generate data that is
not-orthogonal to a desired invariant subspace (with very very high
probability). So the quality
of the RNG does not need to be statistical.

The real double-precision random number generator is `dlaruv`
not `dlarnv`. The second function `dlarnv` is a wrapper that
does seemingly odd things like `work in 128-length` segments,
until you realize that `dlaruv` only works in 128-length segments.
(Weird!) Anyhoo.

[Here's a link to the double precision generator dlaruv.f in OpenBLAS](https://github.com/xianyi/OpenBLAS/blob/974acb39ff86121a5a94be4853f58bd728b56b81/lapack-netlib/SRC/dlaruv.f)

The routine is fairly simple. Have a set of fixed values and
do some integer operations/etc. to get pseudorandom values out.
Specifically
> ```  
  This routine uses a multiplicative congruential method with modulus
  2**48 and multiplier 33952834046453 (see G.S.Fishman,
  'Multiplicative congruential random number generators with modulus
  2**b: an exhaustive analysis for b = 32 and a partial analysis for
  b = 48', Math. Comp. 189, pp 331-344, 1990). ```

See the appendix for code to extract this fixed array of values.

We set these values into a constant global array

    const _dlaruv_mm=tuple(tuple(494,322,2508,2549), ...

    julia> typeof(_dlaruv_mm)
    NTuple{128, NTuple{4, Int64}}

Next up, we port the actual RNG. In this case, I was actually able to
copy/paste much of the FORTRAN code.

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
          # Multiply the seed by i-th power of the multiplier modulo 2**48
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

          # Convert 48-bit integer to a real number in the interval (0,1)
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

This uses the global constant array to avoid dealing
with that big array of values each function call; otherwise, the
only difference is that we use a while loop instead of a
GOTO loop in fortran.

This gives exactly the same results as the FORTRAN call


    _dlaruv_blas!(
          iseed::Ref{NTuple{4,Int}},
          n::Int,
          x::Vector{Float64}) =
      ccall((LinearAlgebra.BLAS.@blasfunc("dlaruv_"), LinearAlgebra.BLAS.libblas), Cvoid,
        (Ptr{LinearAlgebra.BlasInt}, Ref{LinearAlgebra.BlasInt}, Ptr{Float64}),
        iseed, n, x)

    seed1 = Base.RefValue(tuple(1,3,5,7))
    seed2 = Base.RefValue(tuple(1,3,5,7))

    x1 = zeros(97); _dlaruv_blas!(seed1, length(x1), x1)
    x2 = similar(x1); _dlaruv!(seed2, length(x2), x2)

    @assert x1 == x2
    @assert seed1[] == seed2[]

    x3 = zeros(143); _dlaruv_blas!(seed1, length(x3), x3)
    x4 = zeros(axes(x3)...); _dlaruv!(seed2, length(x4), x4)

    @assert x3 == x4
    @assert seed1[] == seed2[]

Awesome. But if you look at `dlarnv`, it also has a statically
allocated length 128 array that is not going to be fun in Julia.
Doable, but not fun... nor at all necessary. We can actually
streamline things by just manually inlining the RNG directly
into `dlarnv`. The only thing we have to remember is that
we need to "refresh" the seed every 128 RV generations.

Since ARPACK only ever calls dlarnv with idist=2, we are just
going to port that code.

    function _dlarnv_idist_2!(iseed::Base.RefValue{NTuple{4,Int}}, n::Int, x::Vector{Float64})
      # manually inline dlaruv to avoid excess cruft...
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

      nrv = 0

      @inbounds for i=1:n
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
          # Multiply the seed by i-th power of the multiplier modulo 2**48
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

          # Convert 48-bit integer to a real number in the interval (0,1)
          rv = R*( Float64( IT1 )+R*( Float64( IT2 )+R*( Float64( IT3 )+R*
                   Float64( IT4 ) ) ) )
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
        # port the line from `dlarnv` for idist=2
        x[i] = 2*rv - 1
      end
      iseed[] = tuple(IT1,IT2,IT3,IT4) # update the seed
    end

This works great!

    seed1 = Base.RefValue(tuple(1,3,5,7))
    seed2 = Base.RefValue(tuple(1,3,5,7))

    x1 = zeros(501); _dlarnv_blas!(2, seed1, length(x1), x1)
    x2 = similar(x1); _dlarnv_idist_2!(seed2, length(x2), x2)

    @assert x1 == x2
    @assert seed1[] == seed2[]

    ##
    x3 = zeros(143); _dlarnv_blas!(2, seed1, length(x3), x3)
    x4 = zeros(axes(x3)...); _dlarnv_idist_2!(seed2, length(x4), x4)

    @assert x3 == x4
    @assert seed1[] == seed2[]    

Exactly the same values and exactly the same seeds.


Things that might happen in the future... I may pull out the key piece
of the random number generator into a separate function...

    rv,IT1,IT2,IT3,IT4,i1,i2,i3,i4 = _dlaruv_rng_step(
        nrv,IT1,IT2,IT3,IT4,i1,i2,i3,i4)

Well... I tried that, and it turns out slower, which brings us to ...

# Julia is faster than Fortran.

What's neat about inlining that code is we get a small performance
boost.

    using BenchmarkTools
    @btime _dlarnv_idist_2!(tseed, length(z), z) setup=begin
      z=zeros(501);
      tseed=Base.Ref(tuple(1,3,5,7)); end
    @btime _dlarnv_blas!(2, tseed, length(z), z) setup=begin
      z=zeros(501);
      tseed=Base.Ref(tuple(1,3,5,7)); end

Gives ...

*Julia*  6.686 μs (0 allocations: 0 bytes)

*Fortran* 7.333 μs (0 allocations: 0 bytes)

So Julia for the win!      

(For reference, with the `_dlaruv_rng_step` idea, it was
  7.0 microseconds, so just a nudge slower. Not sure why as
  it would seem to inline to exactly the same code. )

# Appendix
## Getting the values from the fortran code.

Here is how we extracted that big array of constant values
from the FORTRAN code for `dlaruv`. I believe it's the same array
for `slaruv` too, but I'd need to double-check that. So this code
is helpful there.  

    code = raw"""
    DATA               ( MM( 1, J ), J = 1, 4 ) / 494, 322, 2508,
    $                   2549 /
    DATA               ( MM( 2, J ), J = 1, 4 ) / 2637, 789, 3754,
    $                   1145 /
    DATA               ( MM( 3, J ), J = 1, 4 ) / 255, 1440, 1766,
    $                   2253 /
    DATA               ( MM( 4, J ), J = 1, 4 ) / 2008, 752, 3572,
    $                   305 /
    DATA               ( MM( 5, J ), J = 1, 4 ) / 1253, 2859, 2893,
    $                   3301 /
    DATA               ( MM( 6, J ), J = 1, 4 ) / 3344, 123, 307,
    $                   1065 /
    DATA               ( MM( 7, J ), J = 1, 4 ) / 4084, 1848, 1297,
    $                   3133 /
    DATA               ( MM( 8, J ), J = 1, 4 ) / 1739, 643, 3966,
    $                   2913 /
    DATA               ( MM( 9, J ), J = 1, 4 ) / 3143, 2405, 758,
    $                   3285 /
    DATA               ( MM( 10, J ), J = 1, 4 ) / 3468, 2638, 2598,
    $                   1241 /
    DATA               ( MM( 11, J ), J = 1, 4 ) / 688, 2344, 3406,
    $                   1197 /
    DATA               ( MM( 12, J ), J = 1, 4 ) / 1657, 46, 2922,
    $                   3729 /
    DATA               ( MM( 13, J ), J = 1, 4 ) / 1238, 3814, 1038,
    $                   2501 /
    DATA               ( MM( 14, J ), J = 1, 4 ) / 3166, 913, 2934,
    $                   1673 /
    DATA               ( MM( 15, J ), J = 1, 4 ) / 1292, 3649, 2091,
    $                   541 /
    DATA               ( MM( 16, J ), J = 1, 4 ) / 3422, 339, 2451,
    $                   2753 /
    DATA               ( MM( 17, J ), J = 1, 4 ) / 1270, 3808, 1580,
    $                   949 /
    DATA               ( MM( 18, J ), J = 1, 4 ) / 2016, 822, 1958,
    $                   2361 /
    DATA               ( MM( 19, J ), J = 1, 4 ) / 154, 2832, 2055,
    $                   1165 /
    DATA               ( MM( 20, J ), J = 1, 4 ) / 2862, 3078, 1507,
    $                   4081 /
    DATA               ( MM( 21, J ), J = 1, 4 ) / 697, 3633, 1078,
    $                   2725 /
    DATA               ( MM( 22, J ), J = 1, 4 ) / 1706, 2970, 3273,
    $                   3305 /
    DATA               ( MM( 23, J ), J = 1, 4 ) / 491, 637, 17,
    $                   3069 /
    DATA               ( MM( 24, J ), J = 1, 4 ) / 931, 2249, 854,
    $                   3617 /
    DATA               ( MM( 25, J ), J = 1, 4 ) / 1444, 2081, 2916,
    $                   3733 /
    DATA               ( MM( 26, J ), J = 1, 4 ) / 444, 4019, 3971,
    $                   409 /
    DATA               ( MM( 27, J ), J = 1, 4 ) / 3577, 1478, 2889,
    $                   2157 /
    DATA               ( MM( 28, J ), J = 1, 4 ) / 3944, 242, 3831,
    $                   1361 /
    DATA               ( MM( 29, J ), J = 1, 4 ) / 2184, 481, 2621,
    $                   3973 /
    DATA               ( MM( 30, J ), J = 1, 4 ) / 1661, 2075, 1541,
    $                   1865 /
    DATA               ( MM( 31, J ), J = 1, 4 ) / 3482, 4058, 893,
    $                   2525 /
    DATA               ( MM( 32, J ), J = 1, 4 ) / 657, 622, 736,
    $                   1409 /
    DATA               ( MM( 33, J ), J = 1, 4 ) / 3023, 3376, 3992,
    $                   3445 /
    DATA               ( MM( 34, J ), J = 1, 4 ) / 3618, 812, 787,
    $                   3577 /
    DATA               ( MM( 35, J ), J = 1, 4 ) / 1267, 234, 2125,
    $                   77 /
    DATA               ( MM( 36, J ), J = 1, 4 ) / 1828, 641, 2364,
    $                   3761 /
    DATA               ( MM( 37, J ), J = 1, 4 ) / 164, 4005, 2460,
    $                   2149 /
    DATA               ( MM( 38, J ), J = 1, 4 ) / 3798, 1122, 257,
    $                   1449 /
    DATA               ( MM( 39, J ), J = 1, 4 ) / 3087, 3135, 1574,
    $                   3005 /
    DATA               ( MM( 40, J ), J = 1, 4 ) / 2400, 2640, 3912,
    $                   225 /
    DATA               ( MM( 41, J ), J = 1, 4 ) / 2870, 2302, 1216,
    $                   85 /
    DATA               ( MM( 42, J ), J = 1, 4 ) / 3876, 40, 3248,
    $                   3673 /
    DATA               ( MM( 43, J ), J = 1, 4 ) / 1905, 1832, 3401,
    $                   3117 /
    DATA               ( MM( 44, J ), J = 1, 4 ) / 1593, 2247, 2124,
    $                   3089 /
    DATA               ( MM( 45, J ), J = 1, 4 ) / 1797, 2034, 2762,
    $                   1349 /
    DATA               ( MM( 46, J ), J = 1, 4 ) / 1234, 2637, 149,
    $                   2057 /
    DATA               ( MM( 47, J ), J = 1, 4 ) / 3460, 1287, 2245,
    $                   413 /
    DATA               ( MM( 48, J ), J = 1, 4 ) / 328, 1691, 166,
    $                   65 /
    DATA               ( MM( 49, J ), J = 1, 4 ) / 2861, 496, 466,
    $                   1845 /
    DATA               ( MM( 50, J ), J = 1, 4 ) / 1950, 1597, 4018,
    $                   697 /
    DATA               ( MM( 51, J ), J = 1, 4 ) / 617, 2394, 1399,
    $                   3085 /
    DATA               ( MM( 52, J ), J = 1, 4 ) / 2070, 2584, 190,
    $                   3441 /
    DATA               ( MM( 53, J ), J = 1, 4 ) / 3331, 1843, 2879,
    $                   1573 /
    DATA               ( MM( 54, J ), J = 1, 4 ) / 769, 336, 153,
    $                   3689 /
    DATA               ( MM( 55, J ), J = 1, 4 ) / 1558, 1472, 2320,
    $                   2941 /
    DATA               ( MM( 56, J ), J = 1, 4 ) / 2412, 2407, 18,
    $                   929 /
    DATA               ( MM( 57, J ), J = 1, 4 ) / 2800, 433, 712,
    $                   533 /
    DATA               ( MM( 58, J ), J = 1, 4 ) / 189, 2096, 2159,
    $                   2841 /
    DATA               ( MM( 59, J ), J = 1, 4 ) / 287, 1761, 2318,
    $                   4077 /
    DATA               ( MM( 60, J ), J = 1, 4 ) / 2045, 2810, 2091,
    $                   721 /
    DATA               ( MM( 61, J ), J = 1, 4 ) / 1227, 566, 3443,
    $                   2821 /
    DATA               ( MM( 62, J ), J = 1, 4 ) / 2838, 442, 1510,
    $                   2249 /
    DATA               ( MM( 63, J ), J = 1, 4 ) / 209, 41, 449,
    $                   2397 /
    DATA               ( MM( 64, J ), J = 1, 4 ) / 2770, 1238, 1956,
    $                   2817 /
    DATA               ( MM( 65, J ), J = 1, 4 ) / 3654, 1086, 2201,
    $                   245 /
    DATA               ( MM( 66, J ), J = 1, 4 ) / 3993, 603, 3137,
    $                   1913 /
    DATA               ( MM( 67, J ), J = 1, 4 ) / 192, 840, 3399,
    $                   1997 /
    DATA               ( MM( 68, J ), J = 1, 4 ) / 2253, 3168, 1321,
    $                   3121 /
    DATA               ( MM( 69, J ), J = 1, 4 ) / 3491, 1499, 2271,
    $                   997 /
    DATA               ( MM( 70, J ), J = 1, 4 ) / 2889, 1084, 3667,
    $                   1833 /
    DATA               ( MM( 71, J ), J = 1, 4 ) / 2857, 3438, 2703,
    $                   2877 /
    DATA               ( MM( 72, J ), J = 1, 4 ) / 2094, 2408, 629,
    $                   1633 /
    DATA               ( MM( 73, J ), J = 1, 4 ) / 1818, 1589, 2365,
    $                   981 /
    DATA               ( MM( 74, J ), J = 1, 4 ) / 688, 2391, 2431,
    $                   2009 /
    DATA               ( MM( 75, J ), J = 1, 4 ) / 1407, 288, 1113,
    $                   941 /
    DATA               ( MM( 76, J ), J = 1, 4 ) / 634, 26, 3922,
    $                   2449 /
    DATA               ( MM( 77, J ), J = 1, 4 ) / 3231, 512, 2554,
    $                   197 /
    DATA               ( MM( 78, J ), J = 1, 4 ) / 815, 1456, 184,
    $                   2441 /
    DATA               ( MM( 79, J ), J = 1, 4 ) / 3524, 171, 2099,
    $                   285 /
    DATA               ( MM( 80, J ), J = 1, 4 ) / 1914, 1677, 3228,
    $                   1473 /
    DATA               ( MM( 81, J ), J = 1, 4 ) / 516, 2657, 4012,
    $                   2741 /
    DATA               ( MM( 82, J ), J = 1, 4 ) / 164, 2270, 1921,
    $                   3129 /
    DATA               ( MM( 83, J ), J = 1, 4 ) / 303, 2587, 3452,
    $                   909 /
    DATA               ( MM( 84, J ), J = 1, 4 ) / 2144, 2961, 3901,
    $                   2801 /
    DATA               ( MM( 85, J ), J = 1, 4 ) / 3480, 1970, 572,
    $                   421 /
    DATA               ( MM( 86, J ), J = 1, 4 ) / 119, 1817, 3309,
    $                   4073 /
    DATA               ( MM( 87, J ), J = 1, 4 ) / 3357, 676, 3171,
    $                   2813 /
    DATA               ( MM( 88, J ), J = 1, 4 ) / 837, 1410, 817,
    $                   2337 /
    DATA               ( MM( 89, J ), J = 1, 4 ) / 2826, 3723, 3039,
    $                   1429 /
    DATA               ( MM( 90, J ), J = 1, 4 ) / 2332, 2803, 1696,
    $                   1177 /
    DATA               ( MM( 91, J ), J = 1, 4 ) / 2089, 3185, 1256,
    $                   1901 /
    DATA               ( MM( 92, J ), J = 1, 4 ) / 3780, 184, 3715,
    $                   81 /
    DATA               ( MM( 93, J ), J = 1, 4 ) / 1700, 663, 2077,
    $                   1669 /
    DATA               ( MM( 94, J ), J = 1, 4 ) / 3712, 499, 3019,
    $                   2633 /
    DATA               ( MM( 95, J ), J = 1, 4 ) / 150, 3784, 1497,
    $                   2269 /
    DATA               ( MM( 96, J ), J = 1, 4 ) / 2000, 1631, 1101,
    $                   129 /
    DATA               ( MM( 97, J ), J = 1, 4 ) / 3375, 1925, 717,
    $                   1141 /
    DATA               ( MM( 98, J ), J = 1, 4 ) / 1621, 3912, 51,
    $                   249 /
    DATA               ( MM( 99, J ), J = 1, 4 ) / 3090, 1398, 981,
    $                   3917 /
    DATA               ( MM( 100, J ), J = 1, 4 ) / 3765, 1349, 1978,
    $                   2481 /
    DATA               ( MM( 101, J ), J = 1, 4 ) / 1149, 1441, 1813,
    $                   3941 /
    DATA               ( MM( 102, J ), J = 1, 4 ) / 3146, 2224, 3881,
    $                   2217 /
    DATA               ( MM( 103, J ), J = 1, 4 ) / 33, 2411, 76,
    $                   2749 /
    DATA               ( MM( 104, J ), J = 1, 4 ) / 3082, 1907, 3846,
    $                   3041 /
    DATA               ( MM( 105, J ), J = 1, 4 ) / 2741, 3192, 3694,
    $                   1877 /
    DATA               ( MM( 106, J ), J = 1, 4 ) / 359, 2786, 1682,
    $                   345 /
    DATA               ( MM( 107, J ), J = 1, 4 ) / 3316, 382, 124,
    $                   2861 /
    DATA               ( MM( 108, J ), J = 1, 4 ) / 1749, 37, 1660,
    $                   1809 /
    DATA               ( MM( 109, J ), J = 1, 4 ) / 185, 759, 3997,
    $                   3141 /
    DATA               ( MM( 110, J ), J = 1, 4 ) / 2784, 2948, 479,
    $                   2825 /
    DATA               ( MM( 111, J ), J = 1, 4 ) / 2202, 1862, 1141,
    $                   157 /
    DATA               ( MM( 112, J ), J = 1, 4 ) / 2199, 3802, 886,
    $                   2881 /
    DATA               ( MM( 113, J ), J = 1, 4 ) / 1364, 2423, 3514,
    $                   3637 /
    DATA               ( MM( 114, J ), J = 1, 4 ) / 1244, 2051, 1301,
    $                   1465 /
    DATA               ( MM( 115, J ), J = 1, 4 ) / 2020, 2295, 3604,
    $                   2829 /
    DATA               ( MM( 116, J ), J = 1, 4 ) / 3160, 1332, 1888,
    $                   2161 /
    DATA               ( MM( 117, J ), J = 1, 4 ) / 2785, 1832, 1836,
    $                   3365 /
    DATA               ( MM( 118, J ), J = 1, 4 ) / 2772, 2405, 1990,
    $                   361 /
    DATA               ( MM( 119, J ), J = 1, 4 ) / 1217, 3638, 2058,
    $                   2685 /
    DATA               ( MM( 120, J ), J = 1, 4 ) / 1822, 3661, 692,
    $                   3745 /
    DATA               ( MM( 121, J ), J = 1, 4 ) / 1245, 327, 1194,
    $                   2325 /
    DATA               ( MM( 122, J ), J = 1, 4 ) / 2252, 3660, 20,
    $                   3609 /
    DATA               ( MM( 123, J ), J = 1, 4 ) / 3904, 716, 3285,
    $                   3821 /
    DATA               ( MM( 124, J ), J = 1, 4 ) / 2774, 1842, 2046,
    $                   3537 /
    DATA               ( MM( 125, J ), J = 1, 4 ) / 997, 3987, 2107,
    $                   517 /
    DATA               ( MM( 126, J ), J = 1, 4 ) / 2573, 1368, 3508,
    $                   3017 /
    DATA               ( MM( 127, J ), J = 1, 4 ) / 1148, 1848, 3525,
    $                   2141 /
    DATA               ( MM( 128, J ), J = 1, 4 ) / 545, 2366, 3801,
    $                   1537 /
    """

    lines = split(code, "\n"; keepempty=false)
    for i = 1:2:length(lines)
      l1 = lines[i]
      l2 = lines[i+1] # next line  with last piece
      vals13str = split(l1, "/")[2] # values 1 to 3
      vals13 = parse.(Int, split(vals13str, ","; keepempty=false))
      vals14 = push!(vals13, parse(Int, split(l2[2:end])[1]))
      println("tuple(", join(vals14, ","), "),")
    end
