See [Overview](https://dgleich.micro.blog/2021/04/21/a-small-project.html) for a
high-level picture, also see [ArpackInJulia.jl on Github](https://github.com/dgleich/ArpackInJulia.jl)

One of the things I wanted to do was check against the true ARPACK functions.
The idea is that we can swap out true ARPACK routines for our routines easily.
So I wanted to be able to call the true ARPACK routines individually myself.

This is remarkably easy and how Arpack.jl works with the functions `aupd` and
`eupd` in ARPACK to implement eigs. So this isn't surprisingly. It's still fun
to see that ease of use translate into new codes.

The idea is that `Arpack.jl` uses `Arpack_jll` to access the compiled FORTRAN
code as a dynamic library. All the binary compiling, etc. is handled by
`Arpack_jll` and `Arpack.jl` can just use that. So we are going to

    import Arpack_jll

so we can get access to

    @show Arpack_jll.libarpack

which gives the path to the ARPACK library in a system independent fashion.
(At some point, this will work on the Apple M1 too, but not yet.) [Probably
will be about the same time as I finish posting these notes though and getting
all of the routines ported.]

Once we have that, we can simply setup a ccall to get what we need! See
the blas.jl functions for many examples of how to do this. Also see Arpack.jl
for even more examples that are specific to ARPACK.

    import Arpack_jll, LinearAlgebra

    function arpack_dsconv(n::Int, ritz::Vector{Float64}, bounds::Vector{Float64},
                          tol::Float64)
      nconv::Int = 0
      ccall((:dsconv_, Arpack_jll.libarpack), Cvoid,
        (Ref{LinearAlgebra.BlasInt},
         Ptr{Float64},
         Ptr{Float64},
         Ref{Float64},
         Ref{LinearAlgebra.BlasInt}),
        n, ritz, bounds, tol, nconv)
      return nconv
    end

Of course, I only want to depend on `Arpack_jll` in the test routines.
(Actually, I'd really only like to depend on it in a subset of the
test routines. But I don't see how to do that right now and it isn't
relevant.)

Put another way, I want to enable the tests to run without `Arpack_jll`
working on a particular platform. And where we don't run the tests that
need it.

To control which tests are run, we can check `ARGS` as described
in this pull request <https://github.com/JuliaLang/Pkg.jl/pull/1226>
where this functionality was implemented in Julia.


    if "arpackjll" in ARGS
      include("arpackjll.jl") # get the helpful routines to call stuff in "arpackjll"
      @testset "arpackjll" begin
        # tests that check against libarpack directly
      end
    end    

So I stuck the call to dsconv into `arpackjll.jl`    

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

This does require me to add `Arpack_jll` and `LinearAlgebra` to the
Test package itself. We can do this by activating the test package
and then adding them.

    # change to the test subdirectory
    ] activate .
    add LinearAlgebra
    add Arpack_jll

And our full piece of new test code

    if "arpackjll" in ARGS
      include("arpackjll.jl") # get the helpful routines to call stuff in "arpackjll"
      @testset "arpackjll" begin
        @testset "dsconv" begin
          soln = arpack_dsconv(10, ones(10), zeros(10), 1e-8)
          @test ArpackInJulia.dsconv(10, ones(10), zeros(10), 1e-8)==soln

          soln = arpack_dsconv(10, zeros(10), ones(10), 1e-8)
          @test ArpackInJulia.dsconv(10, zeros(10), ones(10), 1e-8)==soln
        end
      end
    end

This flips ones and zeros, which should show converged in one of the cases.

And what's weird, it fails. I always get back arpack_dsconv==-1. Hmm... will debug
in the next post!
