Here, we pick up with debugging the call to `dsconv` within the ARPACK library
itself from the last post. This is when the repo was at
[this commit](https://github.com/dgleich/ArpackInJulia.jl/tree/9597f00dddd507d0fc57a8145e954761e049f07d)

When we run the tests against the Arpack wrapped function directly, it errors.

    using Pkg; Pkg.test("ArpackInJulia"; test_args=["arpackjll"])

In particular, the following code always returns "-1"

    import Arpack_jll, LinearAlgebra

    function arpack_dsconv(n::Int, ritz::Vector{Float64}, bounds::Vector{Float64},
                          tol::Float64)
      nconv::LinearAlgebra.BlasInt = -1
      ccall((:dsconv_, Arpack_jll.libarpack), Cvoid,
        (Ref{LinearAlgebra.BlasInt},
         Ptr{Float64},
         Ptr{Float64},
         Ref{Float64},
         Ref{LinearAlgebra.BlasInt}),
        n, ritz, bounds, tol, nconv)
      return nconv
    end

Ugh. This is tough to debug. It _seems_ like nconv always isn't getting updated.
One of my common refrains is that problems occur at _boundaries_ between code.
Interfaces are very hard to get right because slight differences in models
can manifest as bugs or just outright incorrect answers. This is exactly what
we are seeing here.

Now, I could just jump to the solution, which is blindingly obvious in retrospect.
On the other hand, there's value in walking through the process I used. (Or
so I tell myself.)

Given that the code above doesn't work, what were my working hypotheses about
the causes?

- So maybe I don't understand how the fortran calling stuff works and I shouldn't
have expected this to get updated.
- Maybe Julia is doing something I don't understand?
- There is some odd interaction for this particular piece of code?
- The Fortran / Arpack code doesn't do what I think it does.

Well, the last one is the the most important as everything follows from my
assumption the fortran code works. So that seems like a good place to start.
(In truth, I did futz around with some of the others, but then quickly
settled on this one.)

Step 1. Do something you know or strongly suspect will work.  
----------

Arpack is designed to work with C and C++. There are just fewer places for
bugs to hide there. So I wanted to make sure I could call this from C and
get the response I expected.

This gave me the following C program.

    #include <inttypes.h>
    typedef int64_t blasint;
    void dsconv_(blasint* n, double* ritz, double* bounds, double* tol, blasint *nconv);

    #include <printf.h>

    int main(int argc, char** argv) {
      double ritz[] = {1.0, 1.0, 1.0};
      double bounds[] = {0.0, 0.0, 0.0};
      double tol = 1e-8;
      blasint n = 3;
      blasint nconv = 4;

      dsconv_(&n, ritz, bounds, &tol, &nconv);

      printf("nconv = %i\n", (int)nconv);
      return 0;
    }

Which I saved into a file: `test_dsconv.c`. This should print 3 if I've understood
this correctly. (And the other Julia code is correct.)

I wanted to compile this against the Julia libraries. Here's how I got it
working. Get the path to libarpack from Julia

    @show Arpack_jll.libarpack_path
    /Users/dgleich/.julia/artifacts/096eefeb84e5b1a463aefac5a561d6109f6e9467/lib/libarpack.2.0.0.dylib    

(Ugh, security breach... yes, my local username is dgleich. Entirely surprising!)
After much fussing with library paths, etc. here's how to get this all compiled.

    gcc test_dsconv.c \
      /Users/dgleich/.julia/artifacts/096eefeb84e5b1a463aefac5a561d6109f6e9467/lib/libarpack.2.0.0.dylib \
      -Wl,-rpath /Users/dgleich/.julia/artifacts/096eefeb84e5b1a463aefac5a561d6109f6e9467/lib/ \
      -Wl,-rpath /Applications/Julia-1.6.app/Contents/Resources/julia/lib/julia

Let's pick this apart.

    gcc test_dsconv.c
    # above is the call to compile this piece of code

    /Users/dgleich/.julia/artifacts/096eefeb84e5b1a463aefac5a561d6109f6e9467/lib/libarpack.2.0.0.dylib
    # above is telling it we want to use this library.

    -Wl,-rpath /Users/dgleich/.julia/artifacts/096eefeb84e5b1a463aefac5a561d6109f6e9467/lib/
    # above is telling gcc to past an option to the linker command. Specifically
    # rpath means where to get runtime libraries ("runtime path")
    # this can also be set with LD_LIBRARY_PATH, but that's... controversial

Let's try without the last link to Julia's libraries to see what happens!

    dgleich@circulant test_dsconv %     gcc test_dsconv.c \
      /Users/dgleich/.julia/artifacts/096eefeb84e5b1a463aefac5a561d6109f6e9467/lib/libarpack.2.0.0.dylib \
      -Wl,-rpath /Users/dgleich/.julia/artifacts/096eefeb84e5b1a463aefac5a561d6109f6e9467/lib/ \

    dgleich@circulant test_dsconv % ./a.out
    dyld: Library not loaded: @rpath/libopenblas64_.dylib
    Referenced from: /Users/dgleich/.julia/artifacts/096eefeb84e5b1a463aefac5a561d6109f6e9467/lib/libarpack.2.0.0.dylib
    Reason: image not found
    zsh: abort      ./a.out

Which is because it can't find libopenblas, referenced by libarpack. So we need
it tell it about all the other common Julia libraries. That's what the final

    -Wl,-rpath /Applications/Julia-1.6.app/Contents/Resources/julia/lib/julia
    # tell it about Julia's libraries.

And now it works great.

    dgleich@circulant test_dsconv % gcc test_dsconv.c \
          /Users/dgleich/.julia/artifacts/096eefeb84e5b1a463aefac5a561d6109f6e9467/lib/libarpack.2.0.0.dylib \
          -Wl,-rpath /Users/dgleich/.julia/artifacts/096eefeb84e5b1a463aefac5a561d6109f6e9467/lib/ \
          -Wl,-rpath /Applications/Julia-1.6.app/Contents/Resources/julia/lib/julia
    dgleich@circulant test_dsconv % ./a.out
    nconv = 3

Okay, so my understanding of the Arpack function was correct. Ugh. That means
I need to learn something about Julia.

Step 2: Better understanding Julia.
---------------

Okay, let's go back to Julia because what we are doing ought to work. Here's
where @code_native and some hints about assembly code come in handy.

    julia> @code_native arpack_dsconv(3, ones(3), zeros(3), 1e-8)
        .section    __TEXT,__text,regular,pure_instructions
    ; ┌ @ none:1 within `arpack_dsconv'
        subq    $40, %rsp
    ; │ @ none:4 within `arpack_dsconv'
    ; │┌ @ essentials.jl:396 within `cconvert'
    ; ││┌ @ refpointer.jl:104 within `convert'
    ; │││┌ @ refvalue.jl:8 within `RefValue'
        movq    %rdi, 32(%rsp)
        vmovsd  %xmm0, 16(%rsp)
        movq    $-1, (%rsp)
    ; │└└└
    ; │┌ @ pointer.jl:65 within `unsafe_convert'
        movq    (%rsi), %rsi
        movq    (%rdx), %rdx
        leaq    32(%rsp), %rdi
        leaq    16(%rsp), %rcx
        movq    %rsp, %r8
        movabsq $dsconv_, %rax
    ; │└
        callq   *%rax
    ; │ @ none:11 within `arpack_dsconv'
        movq    $-1, %rax
        addq    $40, %rsp
        retq
        nopw    %cs:(%rax,%rax)
    ; └

If you look at this code long enough and google around about how `amd64` assembly
works (or `x86_64` if you prefer, although it's really amd64...) then the return
value of the function is in register `rax`. This code

    movq    $-1, %rax

sets the return register to -1 right before returning (`retq`)    

If you want to learn more about calling conventions and register conventions
on `amd64`, see [this powerpoint from cs217 at Princeton](https://www.cs.princeton.edu/courses/archive/spr18/cos217/lectures/15_AssemblyFunctions.pdf)

Okay, so what's happening is that we are just always returning `-1` independently
of the function call. This means that the Julia compiler is treating the value
of `nconv` as a constant and optimizing the code.

Another debugging tenant is try things you think will likely work. Here's something
I strongly suspect will work. And just for fun, I tried the following code.

    import Arpack_jll, LinearAlgebra

    function arpack_dsconv(n::Int, ritz::Vector{Float64}, bounds::Vector{Float64},
                          tol::Float64)
      nconv::Vector{LinearAlgebra.BlasInt} = [-1]
      ccall((:dsconv_, Arpack_jll.libarpack), Cvoid,
        (Ref{LinearAlgebra.BlasInt},
         Ptr{Float64},
         Ptr{Float64},
         Ref{Float64},
         Ref{LinearAlgebra.BlasInt}),
        n, ritz, bounds, tol, nconv)
      return nconv[1]
    end

Why do I think this will work? Because I'm passing an array and then explicitly
looking at an entry. This is what I want to do, just without the array allocation
that creating the vector causes. (Julia's memory allocater is good, but not
that good.)

    julia> arpack_dsconv(3, ones(3), zeros(3), 1e-8)
    1-element Vector{Int64}:
     3

Hmm... so I need [to learn more about Ref](https://docs.julialang.org/en/v1/base/c/#Core.Ref)
to understand why the other call didn't work.

This line stuck out to me.

> In Julia, Ref objects are dereferenced (loaded or stored) with [].

Anyway, at this point, I realized my mental model was wrong and got to ...

The fix
=======

Now, here's the fix. The issue was the `ccall` will auto wrap
an intrinsic type (i.e. int, float, etc.) with a Ref if you tell
it the ref type. But this doesn't refer to the original value, it
makes a new copy of it and refers to that.  So here's working code.

import Arpack_jll, LinearAlgebra

function arpack_dsconv(n::Int, ritz::Vector{Float64}, bounds::Vector{Float64},
                      tol::Float64)
  nconv = Ref{LinearAlgebra.BlasInt}(-1)
  ccall((:dsconv_, Arpack_jll.libarpack), Cvoid,
    (Ref{LinearAlgebra.BlasInt},
     Ptr{Float64},
     Ptr{Float64},
     Ref{Float64},
     Ref{LinearAlgebra.BlasInt}),
    n, ritz, bounds, tol, nconv)
  return nconv[]
end

The error in my mental model
----------------------------

As I said above, bugs exist on boundaries. Here's the origin of this
bug and my incorrect mental model of the code.

My mental model based on C/C++

    nconv::Int = -1
    ccall((:func, libfunc), Cvoid, (Ref{Int}), nconv) # <-> func(&nconv) in C

But that's not right! Here's what the code actually does.

    nconv::Int = -1
    ccall((:func, libfunc), Cvoid, (Ref{Int}), nconv) # <-> int newint = nconv; func(&int(newint)) in C

Which of course explains why nconv isn't updated. It wasn't even passed in!

Suggestion
----------

- Add a note to `Ref` documentation on this
- Add a note to `ccall` documentation on this

Is it worth it for one person with an incorrect mental model? I think yes
because the `ccall` area is likely to be a source of similar bugs
and it's needed to document the assumptions! This is already partially done
in this paragraph

> When passed as a `ccall` argument (either as a `Ptr` or `Ref` type), a `Ref`
> object will be converted to a native pointer to the data it references.
> For most `T`, or when converted to a `Ptr{Cvoid}`, this is a pointer to the
> object data. When `T` is an `isbits` type, this value may be safely mutated,
> otherwise mutation is strictly undefined behavior.
