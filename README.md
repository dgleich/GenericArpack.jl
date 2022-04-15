GenericArpack.jl
===============

This is a pure-Julia translation of the Arpack library. As such, it works
with generic real-valued types in Julia. It's only dependency is `LinearAlgebra.jl`
in the standard library and includes self-contained BLAS routines for the generic
tools it requires. 

```
using GenericArpack
using DoubleFloats
A = Symmetric(sprand(Double64, 100000, 100000, 5/100000) |> A -> A + A')
... 
julia> eigs(A, 2; ncv=10)
ArpackInJulia.ArpackEigen{...}
eigenspace: LM
values:
2-element Vector{Double64}:
 -4.156547290415474
  5.757828927650802
...
```

The library also supports mixed precision, although, this needs to be used carefully. 
```
julia> eigs(Float16, Float64, A, 2; ncv=10) # use Float16 for vectors, Float64 for Arnoldi info, Double64 for A
values:
2-element Vector{Float64}:
 -5.043265045062479
  6.522275837391058
```
The library is better used with use high precision types, like Double64.

Because it was fairly trivial to do in Julia, `GenericArpack.jl` also has a specialized Hermitian eigensolver.
This has had more limited testing, but should work for many cases. This is useful because it gives a
complex-valued SVD without using the non-symmetric solver as in `Arpack.jl`.

See also the packages
- [Arpack.jl](https://github.com/JuliaLinearAlgebra/Arpack.jl)
- [KrylovKit.jl](https://github.com/Jutho/KrylovKit.jl)
- [ArnoldiMethod.jl](https://github.com/haampie/ArnoldiMethod.jl)

Key Functionality
=================

Arpack has a tremendous amount of functionality. Right now, only the symmetric (and by Julia magic, the new Hermitian)
solvers are implemented. See below for more on the list of functionality. Not everything has been testing, but
- Real Symmetric, Generalized Symmetric, Shift and Invert, with mixed/high precision
- Complex Hermitian, Generalized Symmetric, Shift and Invert, with mixed/high precision
- Real and Complex SVD via the Normal equations with mixed/high precision. 
- Smallest singular subspace estimation. 

More details on functionality
-----------------------------

The Arpack symmetric eigensolver functionality is all here and this is ported 
in its entirety. No work has yet been done on the non-symmetric eigensolver (yet).
Right now, 

|	Status	|	Information	|
|	--------	|	--------	|
|	post-beta	|	Most cases should work, there may be edge cases	|
|	beta	|	A few edges cases are likely to appear	|
|	alpha	|	Limited testing, simple things will probably work	|
|	pre-alpha	|	Limited or virtually no testing, likely to have issues	|
|	coded	|	No idea, but the code is there. 	|


|	Functionality	|	Types	|	Status	|	Notes	|
|	--------	|	--------	|	--------	|	--------	|
|	simple real symmetric eigenvalues	|	Float64	|	post-beta	|	Bitwise matches arpack_jll	|
|	generalized real symmetric eigenvalues 	|	Float64	|	beta	|	Bitwise matches arpack_jll	|
|	singular value decomposition 	|	Float64	|	beta	|	Uses normal equations, sorry Gene	|
|	simple shift-invert symmetric eigenvalues	|	Float64	|	pre-alpha	|	not yet tested	|
|	generalized shift-invert symmetric eigenvalues	|	Float64	|	pre-alpha	|	not yet tested	|
|	generalized buckle symmetric eigenvalues	|	Float64	|	pre-alpha	|	not yet tested	|
|	generalized cayley symmetric eigenvalues 	|	Float64	|	pre-alpha	|	not yet tested	|
|	**Complex**	|		|		|		|
|	simple complex hermitian eigenvalues	|	ComplexF64	|	alpha 	|	no specialized solver in Arpack.jl	|
|	generalized complex hermitian eigenvalues 	|	ComplexF64	|	alpha 	|		|
|	complex singular value decomposition	|	ComplexF64	|	alpha 	|		|
|	**High precision**	|		|		|		|
|	simple real symmetric eigenvalues	|	Double64, Float64x2, etc	|	alpha 	|	no specialized solver in Arpack.jl	|
|	generalized real symmetric eigenvalues 	|	Double64, Float64x2, etc	|	alpha 	|		|
|	singular value decomposition 	|	Double64, Float64x2, etc	|	alpha 	|		|
|	**Mixed type**	|		|		|		|
| 	all previous cases	|	Allows Mixed Types	|	alpha 	|	no specialized solver in Arpack.jl	|
|	**Exotic features**	|		|	--------	|		|
|	user-computed shifts	|	 	|	coded	|		|
|	shift invert, buckling, cayley	|	ComplexF64, Float64	|	coded	|		|
|	shift invert, buckling, cayley	|	high-precision	|	coded	|		|

History
-------
This started as a "me-project" to work on to I can learn something and see how various ideas work.

The goal of this exercise is to port the double-precision ARPACK
for symmetric matrices in Julia. Including all ARPACK stuff. So this should
give "exactly" what ARPACK does but be a pure Julia implementation.
(Where exactly is ... it should be executing roughly the same sequence of
floating point operations and can differ on levels that would be expected
for different compilers compiling the same code.)

- not a goal to "Julia-ize" the package; I want to keep as close to the FORTRAN
  as possible so that I might be able to replace calls to Julia's Arpack.saupd /
  Arpack.seupd (which call the Fortran library) with this code; 
  while this is possible, it was easier to use new features from `GenericArpack.jl`
  to implement superior interfaces.
- small internal function changes are okay, e.g. ARPACK has various debugging
  and timing stuff that would need to be done differently in Julia.
- small simplifications, e.g. if a function computes a single Int, we can
  rewrite that to return the Int rather than writing it into an array like in
  FORTRAN.
- Why? Why not implement my own ImplicitRestart/Eigensolver/Etc.? Simple: I trust
  ARPACK. Also, I want to understand exactly what the symmetric ARPACK solver is doing.
- Why not use a Fortran->Julia compiler? Well, I could. But I could also do
  this and learn intimate details of how it works to help out in teaching :)
- I want to document some of the cool stuff in ARPACK!

Along the way, it seemed like in many cases it was possible to get
_bitwise equivalent floating point_ results for Float64 types. So 
we seek to do that. 
