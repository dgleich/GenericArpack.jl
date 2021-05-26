Next up on our list of routines to convert is the sorting routine `dsortr`.
Again, this _could_ be replaced with internal Julia stuff. However, this sorting
routine does something that shows up quite often in practice. You want to sort
array `x`, and take array `y` along for the ride!

The pseudocode is simple:

      function sort_permute(x,y)
         p = sortperm(x)
         return x[p], y[p]
      end
      a = [5, 8, 9, 1, 4, 3, 2]
      b = [3., 4., 5., 2., 1., 9., 8.]
      sort_permute(a,b)
      ([1, 2, 3, 4, 5, 8, 9], [2.0, 8.0, 9.0, 1.0, 3.0, 4.0, 5.0])

This is a recurring pattern.
[Here's a link to a question on stack overflow.](https://stackoverflow.com/questions/5872415/c-sort-vector-of-t-according-to-vector-of-double).
That post references an [old blog post](https://web.archive.org/web/20090618074917/http://www.stanford.edu/~dgleich/notebook/2006/03/sorting_two_arrays_simultaneou.html) (updated to Wayback machine).
And hey, I think I know that fellow. Good to know that things haven't changed
too much in 15 years! (And that I'm still working on the same old stuff...
hmm... maybe that isn't so good.)

So let's solve this problem in Julia again. Or wait, let's not.
Off the top of my head, I don't know how to do this in Julia without allocating
the extra permutation vector. After glancing around sort.jl in base, it seems
like you could probably do this by ... something complicated by creating a new
array type that wraps the pair of arrays.

Uh oh, the problem with having an idea...

    struct SortPermutePair{U,V} <: AbstractArray{Tuple{U,V},1}
      x::Vector{U}
      y::Vector{V}
    end

... trying to resist urge to test it ...

    Base.size(z::SortPermutePair{U,V}) where {U,V} = (min(length(z.x),length(z.y)),)
    Base.getindex(z::SortPermutePair{U,V}, i::Integer) where {U,V} = (z.x[i], z.y[i])
    Base.setindex!(z::SortPermutePair{U,V}, v::Tuple, i::Integer) where {U,V} = begin z.x[i] = v[1];  z.y[i] = v[2] end

... yeah, not gonna work. Let's see how to do this in Julia before
actually porting the ARPACK routine.    

    a = [5, 8, 9, 1, 4, 3, 2]
    b = [3., 4., 5., 2., 1., 9., 8.]
    sort!(SortPermutePair(a,b), by=x->x[1])

So that wasn't so bad! No extra memory allocated... (right?!?) This shouldn't
allocate extra memory...

    a = [5, 8, 9, 1, 4, 3, 2]
    b = [3., 4., 5., 2., 1., 9., 8.]
    @time sort!(SortPermutePair(a,b), by=x->x[1]);
    @time sort!(SortPermutePair(a,b), by=x->x[1]);

    julia> @time sort!(SortPermutePair(a,b), by=x->x[1]);
      0.050316 seconds (58.15 k allocations: 3.460 MiB, 99.04% compilation time)

    julia> @time sort!(SortPermutePair(a,b), by=x->x[1]);
      0.050027 seconds (58.15 k allocations: 3.460 MiB, 99.08% compilation time)

Hmm. Well, that'll have to wait for a future blog post to debug. My guess
is that I need type annotations in all the functions, which is something
we'll return to in a future post. (Okay -- I had to at least try `@btime`...)
I'm very bad at tangents here! That's the fun of this post series. I'm going
to go into all of those and not edit them out.

    using BenchmarkTools
    @btime sort!(SortPermutePair($a,$b), by=x->x[1]); setup=begin
      a=[5, 8, 9, 1, 4, 3, 2]
      b=[3., 4., 5., 2., 1., 9., 8.]
    end

This gives

    52.661 ns (1 allocation: 80 bytes)
    7-element Vector{Float64}:
    ...

Well, that is slightly better. I'm not sure why we see one allocation. This
should be doable with zero allocations. (I think it might somewhere
be allocating that last output vector???)

    @btime sort!(a); setup=(a=randn(10))
      44.800 ns (0 allocations: 0 bytes)

Maybe it's a function boundary thing?

    function sort_permute!(x,y;kwargs...)
      sort!(SortPermutePair(x,y);kwargs...,by=first)
      x,y
    end   

    julia> @time sort_permute!(a,b)
      0.000006 seconds (2 allocations: 112 bytes)
    ([1, 2, 3, 4, 5, 8, 9], [2.0, 8.0, 9.0, 1.0, 3.0, 4.0, 5.0])   

And this clearly has a problem with allocating extra memory somewhere

    julia> a = randn(10^5); b=randn(10^5); @time sort_permute!(a,b)
      0.014885 seconds (3 allocations: 781.422 KiB)   

And now you have one of my least favorite things in Julia. No mental model
for where/how these decisions/impacts get made. This all is designed
to be allocated on the stack, so there should be no need to go into the weeds.
Anyway, a few iterations later, we end up with one fix. Using that
new function x->x[1] is what causes the huge hit to the compiler time.

    # this occurs no matter how many times you run it as I'm running from
    # the repl. as it seems to need to precompile this function x->x[1] everytime.
    julia> @time sort!(SortPermutePair(a,b), by=x->x[1]);
      0.050027 seconds (58.15 k allocations: 3.460 MiB, 99.08% compilation time)

    julia> @time sort!(SortPermutePair(a,b), by=first);
      0.004317 seconds (152 allocations: 7.844 KiB, 99.72% compilation time)

    julia> @time sort!(SortPermutePair(a,b), by=first);
      0.000011 seconds (3 allocations: 144 bytes)

Big long diversion before we get to the actual Arpack Fortran sorting routine
I wanted to get to.  That'll have to wait for another day! So where are we.
We demonstrated it's fairly easy to do this in Julia although there is one
allocation I'd like to get rid of that still seems to be there. Maybe I
should use the code linter or something like that. (We'll see `track allocations`
in a future post.)
