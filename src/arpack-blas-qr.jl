import ArpackInJulia: _dnrm2_unroll_ext, _dlapy2_julia, _dscal!, @jl_arpack_check_length
"""
*> DGEQR2 computes a QR factorization of a real m-by-n matrix A:
*>
*>    A = Q * ( R ),
*>            ( 0 )
*>
*> where:
*>
*>    Q is a m-by-m orthogonal matrix;
*>    R is an upper-triangular n-by-n matrix;
*>    0 is a (m-n)-by-n zero matrix, if m > n.
*> param[in] M
*>          M is INTEGER
*>          The number of rows of the matrix A.  M >= 0.
*>
*> param[in] N
*>          N is INTEGER
*>          The number of columns of the matrix A.  N >= 0.
*>
*> param[in,out] A
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the m by n matrix A.
*>          On exit, the elements on and above the diagonal of the array
*>          contain the min(m,n) by n upper trapezoidal matrix R (R is
*>          upper triangular if m >= n); the elements below the diagonal,
*>          with the array TAU, represent the orthogonal matrix Q as a
*>          product of elementary reflectors (see Further Details).
*>
*> param[in] LDA
*>          LDA is INTEGER
*>          The leading dimension of the array A.  LDA >= max(1,M).
*>
*> param[out] TAU
*>          TAU is DOUBLE PRECISION array, dimension (min(M,N))
*>          The scalar factors of the elementary reflectors (see Further
*>          Details).
*>
*> param[out] WORK
*>          WORK is DOUBLE PRECISION array, dimension (N)
*>
Return 
*> param[out] INFO
*>          INFO is INTEGER
*>          = 0: successful exit
*>          < 0: if INFO = -i, the i-th argument had an illegal value
"""
function _dgeqr2!(
    A::AbstractMatrix{T},
    tau::AbstractVecOrMat{T},
    work::AbstractVecOrMat{T},
) where {T <: AbstractFloat}
  info = 0 
  m, n = size(A)
  @assert m >= 0 
  @assert n >= 0 
  k = min(m,n)
  @jl_arpack_check_length(tau, k)
  @jl_arpack_check_length(work, n)

  for i=1:k
    tau[i] = _dlarfg!(@view(A[i:m,i]))
    if i < n 
      # apply H[i] to A[i:m, i+1:n] from the left
      aii = A[i,i]
      A[i,i] = one(T)
      #_dlarf_left()
      _dlarf_left!(@view(A[i:m,i]), tau[i], @view(A[i:m,i+1:end]), work)
      A[i,i] = aii
    end 
  end

  return info
end

"""
37 *> DLARFG generates a real elementary reflector H of order n, such
38 *> that
39 *>
40 *>       H * ( alpha ) = ( beta ),   H**T * H = I.
41 *>           (   x   )   (   0  )
42 *>
43 *> where alpha and beta are scalars, and x is an (n-1)-element real
44 *> vector. H is represented in the form
45 *>
46 *>       H = I - tau * ( 1 ) * ( 1 v**T ) ,
47 *>                     ( v )
48 *>
49 *> where tau is a real scalar and v is a real (n-1)-element
50 *> vector.
51 *>
52 *> If the elements of x are all zero, then tau = 0 and H is taken to be
53 *> the unit matrix.
54 *>
55 *> Otherwise  1 <= tau <= 2.
56 
57 *
58 *  Arguments:
59 *  ==========
60 *
61 *> param[in] N
62 
63 *>          N is INTEGER
64 *>          The order of the elementary reflector.
65 
66 *>
67 *> param[in,out] ALPHA
68 
69 *>          ALPHA is DOUBLE PRECISION
70 *>          On entry, the value alpha.
71 *>          On exit, it is overwritten with the value beta.
73 *>
74 *> param[in,out] X
76 *>          X is DOUBLE PRECISION array, dimension
77 *>                         (1+(N-2)*abs(INCX))
78 *>          On entry, the vector x.
79 *>          On exit, it is overwritten with the vector v.
81 *>
82 *> param[in] INCX
84 *>          INCX is INTEGER
85 *>          The increment between elements of X. INCX > 0.
86 
87 *>
88 *> param[out] TAU
90 *>          TAU is DOUBLE PRECISION
91 *>          The value tau.

Julia Arguments
---------------
- `x` the trailing vector. (Can be a view of something else); this will be modified in place
- `alpha` the first element of the full x. 

Julia Return value
------------------
- (beta, tau): 
"""
function _dlarfg!(alpha::T, x::AbstractVector{T}) where T
  if length(x) <= 0
    return (alpha, zero(T))
  end
  xnorm = _dnrm2_unroll_ext(x)
  tau = zero(T) # 
  beta = alpha
  if xnorm == 0 
    # tau = zero(T) initialized above
    # just return 
  else
    # beta = -sign( dlapy2( alpha, xnorm ), alpha )
    beta = -copysign(_dlapy2_julia(alpha, xnorm), alpha)
    # safmin = dlamch( 'S' ) / dlamch( 'E' )
    safmin = floatmin(T)/(eps(T)/2) 
    knt = 0 
    if abs(beta) < safmin
      # xnorm, beta may be inaccurate, scale X and recompute them
      rsafmn = one(T)/safmin
      while knt < 20
        knt += 1
        _dscal!(rsafmn, x)
        beta = beta*rsafmn
        alpha = alpha*rsafmn
        if abs(beta) >= safmin 
          break
        end
      end
      # new beta is at most 1, at least safmin
      xnorm = ArpackInJulia._dnrm2_unroll_ext(x)
      beta = -copysign(_dlapy2_julia(alpha, xnorm), alpha)
    end
    tau = (beta-alpha)/beta
    _dscal!(one(T)/(alpha-beta), x)
    # if alpha is subnormal, it may lose relative accuracy
    for j=1:knt
      beta = beta*safmin
    end
  end
  return (beta, tau)
end

function _dlarfg!(x::AbstractVector{T}) where T
  beta, tau = _dlarfg!(x[1], @view(x[2:end]))
  x[1] = beta 
  return tau
end 

"""*>
38 *> DLARF applies a real elementary reflector H to a real m by n matrix
39 *> C, from either the left or the right. H is represented in the form
40 *>
41 *>       H = I - tau * v * v**T
42 *>
43 *> where tau is a real scalar and v is a real vector.
44 *>
45 *> If tau = 0, then H is taken to be the unit matrix.
47 *
48 *  Arguments:
49 *  ==========
50 *
51 *> param[in] SIDE
52 
53 *>          SIDE is CHARACTER*1
54 *>          = 'L': form  H * C
55 *>          = 'R': form  C * H
56 
57 *>
58 *> param[in] M
59 
60 *>          M is INTEGER
61 *>          The number of rows of the matrix C.
62 
63 *>
64 *> param[in] N
65 
66 *>          N is INTEGER
67 *>          The number of columns of the matrix C.
68 
69 *>
70 *> param[in] V
71 
72 *>          V is DOUBLE PRECISION array, dimension
73 *>                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
74 *>                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
75 *>          The vector v in the representation of H. V is not used if
76 *>          TAU = 0.
77 
78 *>
79 *> param[in] INCV
80 
81 *>          INCV is INTEGER
82 *>          The increment between elements of v. INCV <> 0.
83 
84 *>
85 *> param[in] TAU
86 
87 *>          TAU is DOUBLE PRECISION
88 *>          The value tau in the representation of H.
89 
90 *>
91 *> param[in,out] C
92 
93 *>          C is DOUBLE PRECISION array, dimension (LDC,N)
94 *>          On entry, the m by n matrix C.
95 *>          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
96 *>          or C * H if SIDE = 'R'.
97 
98 *>
99 *> param[in] LDC
100 
101 *>          LDC is INTEGER
102 *>          The leading dimension of the array C. LDC >= max(1,M).
103 
104 *>
105 *> param[out] WORK
106 
107 *>          WORK is DOUBLE PRECISION array, dimension
108 *>                         (N) if SIDE = 'L'
109 *>                      or (M) if SIDE = 'R'

Julia Arguments
---------------
- tau, 
"""
function _dlarf_left!(
  v::AbstractVecOrMat{T}, 
  tau::T, 
  C::AbstractMatrix{T}, 
  work::AbstractVecOrMat{T}
) where T 
  lastv = 0 
  lastc = 0 
  m = length(v)
  if tau != 0
    # set up variables for scanning V, lastv begins 
    lastv = m
    i = 1 + (lastv-1)
    while lastv > 0 && v[i] == 0 
      lastv -= 1
      i -= 1
    end 
    # scan for the last non-zero column in C[1:lastv,:]
    lastc = _last_nonzero_column(@view(C[1:lastv,:]))
    # form H*c
    if lastv > 0
      # w(1:lastc,1) := C(1:lastv,1:lastc)**T * v(1:lastv,1)
      mul!(@view(work[1:lastc]), adjoint(@view(C[1:lastv,1:lastc])), @view(v[1:lastv]))
      # C(1:lastv,1:lastc) := C(...) - v(1:lastv,1) * w(1:lastc,1)**T
      _dger!(-tau, @view(v[1:lastv]), @view(work[1:lastc]), @view(C[1:lastv,1:lastc]))
    end
  end 
end

function _dlarf_right!(
  v::AbstractVecOrMat{T}, 
  tau::T, 
  C::AbstractMatrix{T}, 
  work::AbstractVecOrMat{T}
) where T 
  lastv = 0 
  lastc = 0 
  m = length(v)
  if tau != 0
    # set up variables for scanning V, lastv begins 
    lastv = m
    i = 1 + (lastv-1)
    while lastv > 0 && v[i] == 0 
      lastv -= 1
      i -= 1
    end 
    # Scan for the last non-zero row in C(:,1:lastv).
    lastc = _last_nonzero_column(transpose(@view(C[:, 1:lastv])))
    # Form  C * H
    if lastv > 0
      # w(1:lastc,1) := C(1:lastc,1:lastv) * v(1:lastv,1)
      mul!(@view(work[1:lastc]), @view(C[1:lastc,1:lastv]), @view(v[1:lastv]))
      # C(1:lastc,1:lastv) := C(...) - w(1:lastc,1) * v(1:lastv,1)**T
      _dger!(-tau, @view(work[1:lastc]), @view(v[1:lastv]), @view(C[1:lastc,1:lastv]))
    end
  end 
end

function _last_nonzero_column(A::AbstractMatrix{T}) where T
  m,n = size(A)
  if n == 0 
    return n 
  elseif A[1,n] != 0 || A[m,n] != 0 # quick return 
    return n
  else
    # scan 
    for lastc = n:-1:1
      for i=1:m
        if A[i,lastc] != 0 
          return lastc 
        end
      end
    end
    return zero(typeof(m))
  end
end 
#=
direct call
_iladlc(A::StridedMatrix) = begin
  ccall((LinearAlgebra.BLAS.@blasfunc("iladlc_"), LinearAlgebra.BLAS.libblas), LinearAlgebra.BlasInt,
    (Ref{LinearAlgebra.BlasInt},Ref{LinearAlgebra.BlasInt},Ref{Float64},Ref{LinearAlgebra.BlasInt}),
    size(A,1),size(A,2),A,stride(A,2))
end
_iladlc(zeros(10,10)) # does give zero
=#


"""
DGER   performs the rank 1 operation

A := alpha*x*y**T + A,

where alpha is a scalar, x is an m element vector, y is an n element
vector and A is an m by n matrix.
"""
function _dger!(alpha::T, x::AbstractVecOrMat{T}, y::AbstractVecOrMat{T}, A::AbstractMatrix{T}) where T
  m,n = size(A)
  for j=1:n
    if y[j] != 0 
      temp = alpha*y[j]
      for i=1:m
        A[i,j] += x[i]*temp
      end
    end
  end
end


"""
DORM2R overwrites the general real m by n matrix C with

Q * C  if SIDE = 'L' and TRANS = 'N', or

Q**T* C  if SIDE = 'L' and TRANS = 'T', or

C * Q  if SIDE = 'R' and TRANS = 'N', or

C * Q**T if SIDE = 'R' and TRANS = 'T',

where Q is a real orthogonal matrix defined as the product of k
elementary reflectors

Q = H(1) H(2) . . . H(k)

as returned by DGEQRF. Q is of order m if SIDE = 'L' and of order n
if SIDE = 'R'.
"""
function dorm2r(::Val{SIDE},::Val{TRANS},
    m::Integer, n::Integer, k::Integer, 
    A::AbstractMatrix{T}, 
    tau::AbstractVector{T}, 
    C::AbstractMatrix{T}, 
    work::AbstractVector{T}
) where {SIDE, TRANS, T}

  info = 0 
  left = SIDE==:L
  notrans = TRANS==:N
  if left 
    nq = m
  else
    nq = n
  end
  if !left && !(SIDE==:R)
    throw(ArgumentError("Unknown parameter for side, $(SIDE), use Val(:L) or Val(:R)"))
  elseif !notrans && !(TRANS==:T)
    throw(ArgumentError("Unknown parameter for side, $(TRANS), use Val(:N) or Val(:T)"))
  end
  @assert(m >= 0)
  @assert(n >= 0)
  @assert(0 <= k <= nq)

  if left && !notrans || (!left && notrans)
    i1 = 1
    i2 = k
    i3 = 1
  else
    i1 = k
    i2 = 1
    i3 = -1
  end

  if left 
    ni = n
    jc = 1
  else
    mi = m
    ic = 1
  end

  for i in range(start=i1, stop=i2, step=i3)
    if left 
      # H(i) is applied to C(i:m,1:n)
      mi = m - i + 1
      ic = i
    else
      # H(i) is applied to C(1:m,i:n)
      ni = n - i + 1
      jc = i 
    end 
    aii = A[i,i]
    A[i,i] = one(T)
    if left
      _dlarf_left!(@view(A[i:(i+mi-1),i]), tau[i], 
        @view(C[ic:ic+mi-1, jc:jc+ni-1]), work)
    else
      _dlarf_right!(@view(A[i:(i+ni-1),i]), tau[i], 
        @view(C[ic:ic+mi-1, jc:jc+ni-1]), work)
    end
    A[i,i] = aii
  end 
  return C
end 
