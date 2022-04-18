```@setup using-pkgs
using GenericArpack, LinearAlgebra
```

# GenericArpack.jl Documentation 

This library is a julia translation of the Arpack library. 
I strongly recommend browsing the documentation on the 
functions in that library: 

- [dsaupd](https://github.com/opencollab/arpack-ng/blob/9233f7f86f063ca6ca3793cb54dec590eb146e10/SRC/dsaupd.f)
- [dseupd](https://github.com/opencollab/arpack-ng/blob/9233f7f86f063ca6ca3793cb54dec590eb146e10/SRC/dseupd.f)

The major portion of the Julia code is a translation of these functions. 

However, we also provide our own user-friendly interface to these codes. This is inspired by the `eigs` and `svds`
interface in the Julia wrappers [Arpack.jl](https://github.com/JuliaLinearAlgebra/Arpack.jl), however it has 
a number of simplifications.

# Eigenvalue, eigenvectors, standard, generalized. 

The documentation needs to talk about matrices and eigenvectors, of course. 

A standard eigenvalue and eigenvector is a pair: 
- ``Ax = \lambda x`` where ``x`` is not-zero and \$||x|| = 1\$ by convention

If ``A``is symmetic, then all of the ``\lambda`` are real. Note that 
the vector ``x`` is only determined up to sign, as well as rotation 
if there are multliple linearly independent 
eigenvectors with the same eigenvalue.

Likewise if ``A`` is Hermitian, then all of the ``\\lambda`` are still real
whereas the vectors are complex-valued. 
Also, the vectors ``x`` are only unique up to scaling by complex-value 
of magnitude 1. (i.e. ``e^{i\theta}`` for an angle ``\theta``). 

A generalized eigenvalue and eigenvector is a pair: 
- ``Ax = \lambda Bx`` where ``x`` is not-zero and ``||Bx|| = 1`` by convention.

The vectors are just as unique as in the other cases. 

There are many complexities with eigenvalues and subtleties in the defintion.
Please see the ARPACK manual for exact detail on what we mean here as
this is a translation of ARPACK. 

# ARPACK Notes and Key Parameters 

The three key parameters of the Arpack methods are:
- `k`
- `which`
- `ncv` 

The goal of Arpack is to compute a partial _eigenspace_ of a matrix. 

The size of this eigenspace is the number of dimension `k` requested. 

In order to decide _which_ eigenspace we are targeting, we need to tell it! 
Arpack uses the natural variable: `which` in order to provide this information.

- `which`: which part of the spectrum to compute. There are a few choices: 
    + `:LM` _(the default)_ find the largest magnitude eigenspace 
    + `:SM` find the smallest magnitude eigenspace
    + `:LA` find the largest algebraic value eigenspace 
    + `:SA` find the smallest algebraic value eigenspace     
    + `:BE` find both smallest and largest algebraic eigenspace

It turns out that ARPACK uses a strategy based on looking for a larger eigenspace than only `k` vectors.
The parameter `ncv` (number of compute vectors) determines this size. These parameters
must satisfy ``k < ncv \le =n`` (where ``n`` is the dimension of the matrix).

!!! tip "Convergence issues"
    If your eigenproblem isn't converging, it is often useful to increase `ncv` instead
    as well as increasing `maxiter`. This is because if the eigenspace you are trying to find is
    poorly separated, it's easier to find if you can find _any_ separation of a larger eigenspace. 

# Examples of real symmetric, complex Hermitian, and singular value decomposition 

The best way to understand the library is just to see a variety of examples.

!!! info "There is more than one way ... "
    There is more than one way to do something in `GenericArpack.jl`. This
    is by design and there are various limits to the interfaces. Also, there
    exist multiple synonymous calls. This is make it easier for someone reading
    your code to understand what you meant, and often because there is a small
    hint about what a natural default type would be. For instance: 
    `hermeigs` causes us to use `ComplexF64` as a default type, whereas
    `symeigs` will use `Float64` instead. Of course, if you all `symeigs`
    and specify `ComplexF64`, it'll work. Likewise, if you call `hermeigs` 
    and specify `Float64`, but on the other hand, your readers will probably
    be annoyed by that confusion, so use sparingly. 

## Find the largest eigenvalues of a real symmetric matrix.
Let `A` be any type in Julia that implements `Base.size` and `LineareAlgebra.mul!(y,A,x)`.
For instance, `A` can be a `SparseMatrixCSC`, a `Matrix`, a `Tridiagonal`, a `Diagonal`
Moreover, we assume that `A` is symmetric (and we won't check). This can be a
common error.

```@repl using-pkgs
n = 100 
A = Tridiagonal(-ones(n-1), 2*ones(n), -ones(n-1))
symeigs(A, 6)
```

```@setup tridiag
using GenericArpack, LinearAlgebra
n = 100 
A = Tridiagonal(-ones(n-1), 2*ones(n), -ones(n-1))
```

We can also use any of the following equivalent calls
```@repl tridiag
eigs(Symmetric(A), 6)
hermeigs(Symmetric(A), 6)
```

## Finding the largest eigenvalues with an operator
All `GenericArpack` calls are mapped to a single computational interface with an `ArpackOp` type.
```@setup tridiag
op = ArpackSimpleOp(A);

symeigs(op, 6)
```

!!! note "On limits of eigs call"
    In the previous block, we can't use `eigs` and must use `symeigs`. That's
    because right now, we only have the symmetric eigensolvers ported. However,
    since the op type has no way to tell eigs that it should pick the symmetric
    or non-symmetric version, we need to tell it. 

    So if you see 
    ```
    MethodError: no method matching eigs(::ArpackSimpleOp{ ... })
    ```
    That means you just need to call `symeigs` instead! 
    This will probaby get fixed at some point. 


## Complex Hermitian
We also extend Arpack's symmetric solvers with the ability to solve Hermitian problems.

```@repl using-pkgs
n = 100 
A = Tridiagonal(-ones(n-1)*1im, 2*ones(n)+0im, ones(n-1)*1im)
hermeigs(A, 6)
```

## SVD (Singular Value Decomposition) via the Normal Equations
Gene Golub often railed against the use of the normal equations for SVD computations. 
However, he would also use them when appropriate. For computing singular 
values via ARPACK, the normal equations offer a few advantages. 
- smallest singular values 
- they reduce computation as the dimension of the vectors is smaller
The downside is that they also reduce maximum obtainable accuracy.
So please consider using the high-precision types we enable in `GenericArpack.jl`,
if accuracy is desired.

(I need a better SVD example here, but alas.)

```@repl 
using GenericArpack, SparseArrays, StableRNGs 
A = sprand(StableRNG(1), 200, 100, 5/200)
svds(A, 3)
```

```@setup realsvd
using GenericArpack, SparseArrays, StableRNGs 
A = sprand(StableRNG(1), 200, 100, 5/200)
```

We can also find the smallest subspace 
```@repl realsvd 
svds(A, 3; which=:SM) # smallest 
```

Or only the singular values (although this doesn't
really save as much work as it might sound like)
and use them to compute the matrix 2-norm. 
```@repl realsvd 
vals = svds(A, 4; which=:BE, ritzvec=false, ncv=12).S # both smallest and largest
[vals[end]/vals[1] cond(Matrix(A))]
```
So we can estimate the matrix 2-norm with some reasonable 
accuracy (in a well-conditioned case). 


## Complex SVD (Singular Value Decomposition)
This operation is not in the pure ARPACK library.

```@repl
using GenericArpack, SparseArrays, StableRNGs 
A = sprand(StableRNG(1), ComplexF64, 200, 100, 5/200)
svds(A, 3)
```

## Real and Complex SVD via an ArpackNormalOp
All of the complex SVD calls for any Arpack matrix
type are mapped to a single call with an `ArpackNormalOp`

```@repl realsvd 
svds(ArpackNormalOp(A), 3; which=:SM) # smallest 
```

# Examples of Generalized Eigenvalue Problems

## Generalized Eigenvalue Example (from Arpack `dsdrv3` sample).

```@repl using-pkgs 
n = 100
A = Tridiagonal(-ones(n-1),2*ones(n),-ones(n-1)).*(n+1)
B = SymTridiagonal(4*ones(n),ones(n-1)).*(1/(6*(n+1)))
vals, vecs = eigs(Symmetric(A), Symmetric(B), 4)
```

## Complex Hermitian Generalized Eigenvalue Example 

Not here yet! Send me one if you have one! 

## Shift-Invert Example

Not here yet! Send me one if you have one! 

## Buckling Mode Example

Not here yet! Send me one if you have one! 

## 

# Examples with high-precision (or low-precision) types

Both `svds` and `symeigs`/`hermeigs` take in a set of types.

## With just one type

Consider our initial tridiagonal matrix ``A``. 
By default, we use a float-type based on the element-type of
the input matrix ``A``. But when the input includes a specific
type, we use that: 

```@repl tridiag
eigs(Float32, Symmetric(A), 6) # use Float32 for all computations
```

If `A` has a different float type, then we may use that. 

```@repl tridiag
Af32 = Float32.(A) # make a Float32 copy
symeigs(Float32.(A), 6) # will switch default to Float32. 
```

But if you give an `op`, that has no natural type, so we default to `Float64`

```@repl tridiag
op = ArpackSimpleOp(Float32.(A))
symeigs(op, 6)
```

(This is because an `ArpackOp` may wrap a variety of information and we don't
need the type information.)

```@repl tridiag
op = ArpackSimpleOp(Float32.(A))
symeigs(Float32, op, 6)
```

## With two types

In fact, the `symeigs` functions all take in _two_ types: the type
of the Eigeninformation (which must be be non-complex) and
the type of the vector information (which can be complex).

```@repl tridiag
op = ArpackSimpleOp(Float16.(A))
symeigs(ComplexF64, Float32, op, 6)
```

!!! tip "These are correct" 
    These eigenvectors are correct, even though they aren't real-valued.
    It's just that they are less unique in the complex plane compared to the
    real-plane. 

This allows us to use higher or lower-precision in the vector compared
with the eigeninformation (which also is used to represent the 
Arnoldi factorization). 

This uses `BigFloat` for the Eigenvalue information and `Float64` for
the vectors. This will realistically limit you to Float64 accuracy
but might handle some edge cases better. (And `BigFloat` should
really be a different type, but this was handy to write as it's 
built into Julia.)

```@repl tridiag
symeigs(Float64, BigFloat, A, 3; ncv=36)
```

# Advanced usage

## Writing your own `ArpackOp`
It is "easy" (says the developer) to write your own `ArpackOp` type. 

Here is the code to wrap a function `F` in an `ArpackOp`. This is
the complete implementation of the `ArpackFunctionOp`

```
struct ArpackSimpleFunctionOp <: ArpackOp
  F::Function 
  n::Int
end 
arpack_mode(::ArpackSimpleFunctionOp) = 1
Base.size(op::ArpackSimpleFunctionOp) = op.n
bmat(::ArpackSimpleFunctionOp) = Val(:I)
opx!(y,op::ArpackSimpleFunctionOp,x) = op.F(y,x)
is_arpack_mode_valid_for_op(mode::Int, ::ArpackSimpleFunctionOp) = mode == 1 
```

## The Arpack Drivers

All of the Arpack interface is through the functions `dsaupd` and `dseupd`.
_(Yes, Lapack-heads, I know that it should be `_saupd` / `_seupd` in proper nomenclature, but I'll get around to fixing that at some point.)_


The key differences from the standard interface are four new
types of arguments for Julia. 

- `ArpackDebug` controls the Arpack debug messages (from `debug.h`)
- `ArpackStats` controls the Arpack statistics collected (from `stats.h`)
- `idonow` is the Julia information for _avoiding_ the reverse communication 
  interface. It's _slightly_ faster. But likely to make compile times 
  longer. I'm leaving it in even though it'll be a source of bugs, ugh. 
- `ArpackState` tracks the computation specific state that was in the
  Fortran `save` variables. This includes things like the random number
  seeds. 

For `stats` if the type is `nothing`, then we don't track stats. 
Likewise for debug. 

The other key difference is that all of the functions _return_ the info 
value instead of setting the _info_ parameter. (I didn't want a needed
Ref type hanging around in the Julia code.)

Here is how a dsaupd call maps between the two libraries. This is
one we use in the test code. 

```
ido = Ref{Int}(0); bmat=:I; which=:LM; 
state = ArpackState{Float64}() # (or whatever the Eigenvector info is...)
stats = nothing # okay to not track...
ierr, state = GenericArpack.dsaupd!(ido, Val(bmat), n, which, nev, tol, resid, ncv, V, ldv, iparam,
  ipntr, workd, workl, lworkl, info_initv;
  state, stats, debug # these are the new parameters, which are Julia keywoard params, so any order is okay! 
)
arierr = arpack_dsaupd!(arido, bmat, n, which, nev, tol, arresid, ncv, arV, ldv, ariparam, 
  aripntr, arworkd, arworkl, lworkl, info_initv)
```

If `tol=0`, the `tol` parameter is typically initalized by Arpack to `dlamch("E")` which is `eps(Float64)/2`.
This parameter is then propagated to future calls since everything in Fortran is pass-by-reference. 

Instead of worrying about this, we just set `tol` on each call to `GenericArpack.dsaupd!`; however
this can cause some small issues if you try and compare our calls directly. So in those cases, 
I recommend setting tol to be anything non-zero so Arpack won't touch it. 

## Getting debug information
You can also pass `ArpackStats` and `ArpackDebug` to any of the higher-level drivers.

```@repl tridiag
eigs(Symmetric(A), 3; debug=ArpackDebug(maupd=2, maup2=1), stats=ArpackStats()) # lots of convergence info
```

