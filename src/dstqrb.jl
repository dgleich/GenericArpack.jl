#=
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dstqrb
c
c\Description:
c  Computes all eigenvalues and the last component of the eigenvectors
c  of a symmetric tridiagonal matrix using the implicit QL or QR method.
c
c  This is mostly a modification of the LAPACK routine dsteqr.
c  See Remarks.
c
c\Usage:
c  call dstqrb
c     ( N, D, E, Z, WORK, INFO )
c
c\Arguments
c  N       Integer.  (INPUT)
c          The number of rows and columns in the matrix.  N >= 0.
c
c  D       Double precision array, dimension (N).  (INPUT/OUTPUT)
c          On entry, D contains the diagonal elements of the
c          tridiagonal matrix.
c          On exit, D contains the eigenvalues, in ascending order.
c          If an error exit is made, the eigenvalues are correct
c          for indices 1,2,...,INFO-1, but they are unordered and
c          may not be the smallest eigenvalues of the matrix.
c
c  E       Double precision array, dimension (N-1).  (INPUT/OUTPUT)
c          On entry, E contains the subdiagonal elements of the
c          tridiagonal matrix in positions 1 through N-1.
c          On exit, E has been destroyed.
c
c  Z       Double precision array, dimension (N).  (OUTPUT)
c          On exit, Z contains the last row of the orthonormal
c          eigenvector matrix of the symmetric tridiagonal matrix.
c          If an error exit is made, Z contains the last row of the
c          eigenvector matrix associated with the stored eigenvalues.
c
c  WORK    Double precision array, dimension (max(1,2*N-2)).  (WORKSPACE)
c          Workspace used in accumulating the transformation for
c          computing the last components of the eigenvectors.
c
c  INFO    Integer.  (OUTPUT)
c          = 0:  normal return.
c          < 0:  if INFO = -i, the i-th argument had an illegal value.
c          > 0:  if INFO = +i, the i-th eigenvalue has not converged
c                              after a total of  30*N  iterations.
c
c\Remarks
c  1. None.
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dswap   Level 1 BLAS that swaps the contents of two vectors.
c     lsame   LAPACK character comparison routine.
c     dlae2   LAPACK routine that computes the eigenvalues of a 2-by-2
c             symmetric matrix.
c     dlaev2  LAPACK routine that eigendecomposition of a 2-by-2 symmetric
c             matrix.
c     dlamch  LAPACK routine that determines machine constants.
c     dlanst  LAPACK routine that computes the norm of a matrix.
c     dlapy2  LAPACK routine to compute sqrt(x**2+y**2) carefully.
c     dlartg  LAPACK Givens rotation construction routine.
c     dlascl  LAPACK routine for careful scaling of a matrix.
c     dlaset  LAPACK matrix initialization routine.
c     dlasr   LAPACK routine that applies an orthogonal transformation to
c             a matrix.
c     dlasrt  LAPACK sorting routine.
c     dsteqr  LAPACK routine that computes eigenvalues and eigenvectors
c             of a symmetric tridiagonal matrix.
c     xerbla  LAPACK error handler routine.
c
c\Authors
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: stqrb.F   SID: 2.5   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c     1. Starting with version 2.5, this routine is a modified version
c        of LAPACK version 2.0 subroutine SSTEQR. No lines are deleted,
c        only commeted out and new lines inserted.
c        All lines commented out have "c$$$" at the beginning.
c        Note that the LAPACK version 1.0 subroutine SSTEQR contained
c        bugs.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dstqrb ( n, d, e, z, work, info )
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    info, n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           d( n ), e( n-1 ), z( n ), work( 2*n-2 )
c
c     .. parameters ..
      Double precision
     &                   zero, one, two, three
      parameter          ( zero = 0.0D+0, one = 1.0D+0,
     &                     two = 2.0D+0, three = 3.0D+0 )
      integer            maxit
      parameter          ( maxit = 30 )
c     ..
c     .. local scalars ..
      integer            i, icompz, ii, iscale, j, jtot, k, l, l1, lend,
     &                   lendm1, lendp1, lendsv, lm1, lsv, m, mm, mm1,
     &                   nm1, nmaxit
      Double precision
     &                   anorm, b, c, eps, eps2, f, g, p, r, rt1, rt2,
     &                   s, safmax, safmin, ssfmax, ssfmin, tst
c     ..
c     .. external functions ..
      logical            lsame
      Double precision
     &                   dlamch, dlanst, dlapy2
      external           lsame, dlamch, dlanst, dlapy2
c     ..
c     .. external subroutines ..
      external           dlae2, dlaev2, dlartg, dlascl, dlaset, dlasr,
     &                   dlasrt, dswap, xerbla
c     ..
c     .. intrinsic functions ..
      intrinsic          abs, max, sign, sqrt
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
c
c$$$      IF( LSAME( COMPZ, 'N' ) ) THEN
c$$$         ICOMPZ = 0
c$$$      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
c$$$         ICOMPZ = 1
c$$$      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
c$$$         ICOMPZ = 2
c$$$      ELSE
c$$$         ICOMPZ = -1
c$$$      END IF
c$$$      IF( ICOMPZ.LT.0 ) THEN
c$$$         INFO = -1
c$$$      ELSE IF( N.LT.0 ) THEN
c$$$         INFO = -2
c$$$      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
c$$$     $         N ) ) ) THEN
c$$$         INFO = -6
c$$$      END IF
c$$$      IF( INFO.NE.0 ) THEN
c$$$         CALL XERBLA( 'SSTEQR', -INFO )
c$$$         RETURN
c$$$      END IF
c
c    *** New starting with version 2.5 ***
c
      icompz = 2
c    *************************************
c
c     quick return if possible
c
      if( n.eq.0 )
     $   return
c
      if( n.eq.1 ) then
         if( icompz.eq.2 )  z( 1 ) = one
         return
      end if
c
c     determine the unit roundoff and over/underflow thresholds.
c
      eps = dlamch( 'e' )
      eps2 = eps**2
      safmin = dlamch( 's' )
      safmax = one / safmin
      ssfmax = sqrt( safmax ) / three
      ssfmin = sqrt( safmin ) / eps2
c
c     compute the eigenvalues and eigenvectors of the tridiagonal
c     matrix.
c
c$$      if( icompz.eq.2 )
c$$$     $   call dlaset( 'full', n, n, zero, one, z, ldz )
c
c     *** New starting with version 2.5 ***
c
      if ( icompz .eq. 2 ) then
         do 5 j = 1, n-1
            z(j) = zero
  5      continue
         z( n ) = one
      end if
c     *************************************
c
      nmaxit = n*maxit
      jtot = 0
c
c     determine where the matrix splits and choose ql or qr iteration
c     for each block, according to whether top or bottom diagonal
c     element is smaller.
c
      l1 = 1
      nm1 = n - 1
c
   10 continue
      if( l1.gt.n )
     $   go to 160
      if( l1.gt.1 )
     $   e( l1-1 ) = zero
      if( l1.le.nm1 ) then
         do 20 m = l1, nm1
            tst = abs( e( m ) )
            if( tst.eq.zero )
     $         go to 30
            if( tst.le.( sqrt( abs( d( m ) ) )*sqrt( abs( d( m+
     $          1 ) ) ) )*eps ) then
               e( m ) = zero
               go to 30
            end if
   20    continue
      end if
      m = n
c
   30 continue
      l = l1
      lsv = l
      lend = m
      lendsv = lend
      l1 = m + 1
      if( lend.eq.l )
     $   go to 10
c
c     scale submatrix in rows and columns l to lend
c
      anorm = dlanst( 'i', lend-l+1, d( l ), e( l ) )
      iscale = 0
      if( anorm.eq.zero )
     $   go to 10
      if( anorm.gt.ssfmax ) then
         iscale = 1
         call dlascl( 'g', 0, 0, anorm, ssfmax, lend-l+1, 1, d( l ), n,
     $                info )
         call dlascl( 'g', 0, 0, anorm, ssfmax, lend-l, 1, e( l ), n,
     $                info )
      else if( anorm.lt.ssfmin ) then
         iscale = 2
         call dlascl( 'g', 0, 0, anorm, ssfmin, lend-l+1, 1, d( l ), n,
     $                info )
         call dlascl( 'g', 0, 0, anorm, ssfmin, lend-l, 1, e( l ), n,
     $                info )
      end if
c
c     choose between ql and qr iteration
c
      if( abs( d( lend ) ).lt.abs( d( l ) ) ) then
         lend = lsv
         l = lendsv
      end if
c
      if( lend.gt.l ) then
c
c        ql iteration
c
c        look for small subdiagonal element.
c
   40    continue
         if( l.ne.lend ) then
            lendm1 = lend - 1
            do 50 m = l, lendm1
               tst = abs( e( m ) )**2
               if( tst.le.( eps2*abs( d( m ) ) )*abs( d( m+1 ) )+
     $             safmin )go to 60
   50       continue
         end if
c
         m = lend
c
   60    continue
         if( m.lt.lend )
     $      e( m ) = zero
         p = d( l )
         if( m.eq.l )
     $      go to 80
c
c        if remaining matrix is 2-by-2, use dlae2 or dlaev2
c        to compute its eigensystem.
c
         if( m.eq.l+1 ) then
            if( icompz.gt.0 ) then
               call dlaev2( d( l ), e( l ), d( l+1 ), rt1, rt2, c, s )
               work( l ) = c
               work( n-1+l ) = s
c$$$               call dlasr( 'r', 'v', 'b', n, 2, work( l ),
c$$$     $                     work( n-1+l ), z( 1, l ), ldz )
c
c              *** New starting with version 2.5 ***
c
               tst      = z(l+1)
               z(l+1) = c*tst - s*z(l)
               z(l)   = s*tst + c*z(l)
c              *************************************
            else
               call dlae2( d( l ), e( l ), d( l+1 ), rt1, rt2 )
            end if
            d( l ) = rt1
            d( l+1 ) = rt2
            e( l ) = zero
            l = l + 2
            if( l.le.lend )
     $         go to 40
            go to 140
         end if
c
         if( jtot.eq.nmaxit )
     $      go to 140
         jtot = jtot + 1
c
c        form shift.
c
         g = ( d( l+1 )-p ) / ( two*e( l ) )
         r = dlapy2( g, one )
         g = d( m ) - p + ( e( l ) / ( g+sign( r, g ) ) )
c
         s = one
         c = one
         p = zero
c
c        inner loop
c
         mm1 = m - 1
         do 70 i = mm1, l, -1
            f = s*e( i )
            b = c*e( i )
            call dlartg( g, f, c, s, r )
            if( i.ne.m-1 )
     $         e( i+1 ) = r
            g = d( i+1 ) - p
            r = ( d( i )-g )*s + two*c*b
            p = s*r
            d( i+1 ) = g + p
            g = c*r - b
c
c           if eigenvectors are desired, then save rotations.
c
            if( icompz.gt.0 ) then
               work( i ) = c
               work( n-1+i ) = -s
            end if
c
   70    continue
c
c        if eigenvectors are desired, then apply saved rotations.
c
         if( icompz.gt.0 ) then
            mm = m - l + 1
c$$$            call dlasr( 'r', 'v', 'b', n, mm, work( l ), work( n-1+l ),
c$$$     $                  z( 1, l ), ldz )
c
c             *** New starting with version 2.5 ***
c
              call dlasr( 'r', 'v', 'b', 1, mm, work( l ),
     &                    work( n-1+l ), z( l ), 1 )
c             *************************************
         end if
c
         d( l ) = d( l ) - p
         e( l ) = g
         go to 40
c
c        eigenvalue found.
c
   80    continue
         d( l ) = p
c
         l = l + 1
         if( l.le.lend )
     $      go to 40
         go to 140
c
      else
c
c        qr iteration
c
c        look for small superdiagonal element.
c
   90    continue
         if( l.ne.lend ) then
            lendp1 = lend + 1
            do 100 m = l, lendp1, -1
               tst = abs( e( m-1 ) )**2
               if( tst.le.( eps2*abs( d( m ) ) )*abs( d( m-1 ) )+
     $             safmin )go to 110
  100       continue
         end if
c
         m = lend
c
  110    continue
         if( m.gt.lend )
     $      e( m-1 ) = zero
         p = d( l )
         if( m.eq.l )
     $      go to 130
c
c        if remaining matrix is 2-by-2, use dlae2 or dlaev2
c        to compute its eigensystem.
c
         if( m.eq.l-1 ) then
            if( icompz.gt.0 ) then
               call dlaev2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2, c, s )
c$$$               work( m ) = c
c$$$               work( n-1+m ) = s
c$$$               call dlasr( 'r', 'v', 'f', n, 2, work( m ),
c$$$     $                     work( n-1+m ), z( 1, l-1 ), ldz )
c
c               *** New starting with version 2.5 ***
c
                tst      = z(l)
                z(l)   = c*tst - s*z(l-1)
                z(l-1) = s*tst + c*z(l-1)
c               *************************************
            else
               call dlae2( d( l-1 ), e( l-1 ), d( l ), rt1, rt2 )
            end if
            d( l-1 ) = rt1
            d( l ) = rt2
            e( l-1 ) = zero
            l = l - 2
            if( l.ge.lend )
     $         go to 90
            go to 140
         end if
c
         if( jtot.eq.nmaxit )
     $      go to 140
         jtot = jtot + 1
c
c        form shift.
c
         g = ( d( l-1 )-p ) / ( two*e( l-1 ) )
         r = dlapy2( g, one )
         g = d( m ) - p + ( e( l-1 ) / ( g+sign( r, g ) ) )
c
         s = one
         c = one
         p = zero
c
c        inner loop
c
         lm1 = l - 1
         do 120 i = m, lm1
            f = s*e( i )
            b = c*e( i )
            call dlartg( g, f, c, s, r )
            if( i.ne.m )
     $         e( i-1 ) = r
            g = d( i ) - p
            r = ( d( i+1 )-g )*s + two*c*b
            p = s*r
            d( i ) = g + p
            g = c*r - b
c
c           if eigenvectors are desired, then save rotations.
c
            if( icompz.gt.0 ) then
               work( i ) = c
               work( n-1+i ) = s
            end if
c
  120    continue
c
c        if eigenvectors are desired, then apply saved rotations.
c
         if( icompz.gt.0 ) then
            mm = l - m + 1
c$$$            call dlasr( 'r', 'v', 'f', n, mm, work( m ), work( n-1+m ),
c$$$     $                  z( 1, m ), ldz )
c
c           *** New starting with version 2.5 ***
c
            call dlasr( 'r', 'v', 'f', 1, mm, work( m ), work( n-1+m ),
     &                  z( m ), 1 )
c           *************************************
         end if
c
         d( l ) = d( l ) - p
         e( lm1 ) = g
         go to 90
c
c        eigenvalue found.
c
  130    continue
         d( l ) = p
c
         l = l - 1
         if( l.ge.lend )
     $      go to 90
         go to 140
c
      end if
c
c     undo scaling if necessary
c
  140 continue
      if( iscale.eq.1 ) then
         call dlascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv+1, 1,
     $                d( lsv ), n, info )
         call dlascl( 'g', 0, 0, ssfmax, anorm, lendsv-lsv, 1, e( lsv ),
     $                n, info )
      else if( iscale.eq.2 ) then
         call dlascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv+1, 1,
     $                d( lsv ), n, info )
         call dlascl( 'g', 0, 0, ssfmin, anorm, lendsv-lsv, 1, e( lsv ),
     $                n, info )
      end if
c
c     check for no convergence to an eigenvalue after a total
c     of n*maxit iterations.
c
      if( jtot.lt.nmaxit )
     $   go to 10
      do 150 i = 1, n - 1
         if( e( i ).ne.zero )
     $      info = info + 1
  150 continue
      go to 190
c
c     order eigenvalues and eigenvectors.
c
  160 continue
      if( icompz.eq.0 ) then
c
c        use quick sort
c
         call dlasrt( 'i', n, d, info )
c
      else
c
c        use selection sort to minimize swaps of eigenvectors
c
         do 180 ii = 2, n
            i = ii - 1
            k = i
            p = d( i )
            do 170 j = ii, n
               if( d( j ).lt.p ) then
                  k = j
                  p = d( j )
               end if
  170       continue
            if( k.ne.i ) then
               d( k ) = d( i )
               d( i ) = p
c$$$               call dswap( n, z( 1, i ), 1, z( 1, k ), 1 )
c           *** New starting with version 2.5 ***
c
               p    = z(k)
               z(k) = z(i)
               z(i) = p
c           *************************************
            end if
  180    continue
      end if
c
  190 continue
      return
c
c     %---------------%
c     | End of dstqrb |
c     %---------------%
c
      end
=#

using LinearAlgebra: SymTridiagonal, norm, UniformScaling

"""
  Computes all eigenvalues and the last component of the eigenvectors
  of a symmetric tridiagonal matrix using the implicit QL or QR method.

  This is mostly a modification of the LAPACK routine dsteqr.
  See Remarks.
"""
function dstqrb!(
  n::Int,
  d::AbstractVector{T},
  e::AbstractVector{T},
  z::AbstractVector{T},
  work::AbstractVector{T},
  ::Union{AbstractArpackState{T},Nothing} # we have state for dispatch
) where {T}
  @jl_arpack_check_length(d,n)
  @jl_arpack_check_length(e,n-1)
  @jl_arpack_check_length(z,n)
  @jl_arpack_check_length(work, max(2*n-2,1))

  conceptual_dstqrb!(SymTridiagonal(@view(d[1:n]), @view(e[1:n-1]));
    z, work)
end

function simple_dsteqr!(
  d::AbstractVector{T},
  e::AbstractVector{T},
  Z::AbstractMatrix{T};
  work=Vector{T}(undef,max(1,2*size(A,1)-2))
) where T
  n = length(d) 
  @jl_arpack_check_length(d,n)
  @jl_arpack_check_length(e,n-1)
  @jl_arpack_check_size(Z, n, n)
  @jl_arpack_check_length(work, max(2*n-2,1))

  flexible_dsteqr!(SymTridiagonal(@view(d[1:n]), @view(e[1:n-1])), Z; 
    work)
end 

function opnorm1(A::SymTridiagonal{T}) where T
  m, n = size(A)
  d, e = A.dv, A.ev
  Tnorm = typeof(float(real(zero(T))))
  Tsum = promote_type(Float64, Tnorm)
  nrm::Tsum = n >= 1 ?
                max(norm(d[1]) + norm(e[1]),
                  norm(d[end]) + norm(e[end])) :
              0

  @inbounds begin
    for j = 2:n-1
      nrmj::Tsum = norm(d[j]) + norm(e[j]) + norm(e[j-1])
      nrm = max(nrm,nrmj)
    end
  end
  return convert(Tnorm, nrm)
end

_dstqrb_maxit(::Type{Float64}) = 30
_dstqrb_maxit(::Type{Float32}) = 30
_dstqrb_maxit(::Type{T}) where T = max(30, 30*ceil(Int, sqrt(log(eps(T))/log(2^-52))))

function conceptual_dstqrb!(A::SymTridiagonal{T};
  z=Vector{T}(undef,size(A,1)),
  work=Vector{T}(undef,max(1,2*size(A,1)-2))) where T

  flexible_dsteqr!(A, z; work)
end 

function flexible_dsteqr!(A::SymTridiagonal{T},
    z::AbstractVecOrMat{T};
    work=Vector{T}(undef,max(1,2*size(A,1)-2)),initz=true) where T
  d = A.dv # same variables at dstqrb
  e = A.ev

  isvector = ndims(z) == 1 
  n = size(A,1)
  #println()
  #println("Starting dstqrb")
  #display(A)
  if initz 
    if isvector 
      fill!(z, zero(T))
      z[n] = one(T)
    else
      copyto!(z, UniformScaling(one(T)))
    end 
  end 

  jtot = 0
  nmaxit = _dstqrb_maxit(T)*n
  l1 = 1
  while l1 <= n
    if l1 > 1
      e[l1-1]=zero(T) # zero out last component
    end
    l, lend = _find_next_block!(l1, d, e)
    if lend > l # so there is actually work to do!
      #println("Working on block $(l:lend)")
      #display(A)
      if isvector 
        info, niter = _process_block(@view(d[l:lend]), @view(e[l:lend-1]),
            @view(z[l:lend]), work, nmaxit - jtot)
      else
        info, niter = _process_block(@view(d[l:lend]), @view(e[l:lend-1]),
            @view(z[1:n,l:lend]), work, nmaxit - jtot)
      end 

      jtot += niter
      # note, we are going to throw an error in the inner-most function
    else
    end
    l1 = lend + 1
  end

  # sort eigenvectors with selection sort
  for ii=2:n
    i = ii-1
    k = i
    p = d[i]
    for j=ii:n # find the smallest element and assocated index.
      if d[j] < p
        k = j
        p = d[j]
      end
    end
    if k != i # swap if need be.
      d[k] = d[i]
      d[i] = p
      if isvector 
        p = z[k]
        z[k] = z[i]
        z[i] = p
      else
        # swap the columns... 
        for r=1:n
          p = z[r,k]
          z[r,k] = z[r,i]
          z[r,i] = p 
        end 
      end 
    end
  end
  return A,z
end

function _find_next_block!(l::Integer, d, e)
  eps = Base.eps(eltype(d))/2
  for m=l:length(e)
    tst = abs(e[m])
    if iszero(tst)
      return l, m
    elseif tst <= sqrt(abs(d[m]))*sqrt(abs(d[m+1]))*eps
      e[m] = 0
      return l, m
    end
  end
  return l, length(d) # total length of block
end

function _process_block(d, e, z, work, maxit)
  T = Float64
  eps = Base.eps(T)/2 # dlamch("E") is eps/2
  eps2 = eps^2
  safmin = floatmin(T) # dlamch("S")
  safmax = 1/safmin
  ssfmax = sqrt(safmax)/3
  ssfmin = sqrt(safmin)/eps2

  A = SymTridiagonal(d,e)
  anorm = opnorm1(A)

  scaleto = anorm > ssfmax ? ssfmax : anorm < ssfmin ? ssfmin : anorm
  _scale_from_to(anorm, scaleto, d)
  _scale_from_to(anorm, scaleto, e)

  ## TODO, handle iterations...
  if abs(d[end]) < abs(d[1])
    #println("qr")
    info, niter = _qr_iteration(d, e, z, work; maxit)
  else
    #println("ql")
    info, niter = _ql_iteration(d, e, z, work; maxit)
  end

  # reverse scaling
  _scale_from_to(scaleto, anorm, d)
  _scale_from_to(scaleto, anorm, e)

  return info, niter
end

function _ql_iteration(d::AbstractVecOrMat{T},
  e::AbstractVecOrMat{T}, z::AbstractVecOrMat{T}, work::AbstractVecOrMat{T};
  maxit::Int = length(d)*30) where T
  eps = Base.eps(T)/2 # dlamch("E") is eps/2
  eps2 = eps^2
  safmin = floatmin(T) # dlamch("S")
  safmax = 1/safmin
  ssfmax = sqrt(safmax)/3
  ssfmin = sqrt(safmin)/eps2
  two = 2*one(T)

  isvector = ndims(z) == 1 

  n = length(d)
  iend = length(d)
  is = 1
  jtot = 0

  while is < iend # l != lend
    #println("At itereration $(jtot) is=$is iend=$iend")
    #display(SymTridiagonal(d, e))
    # look for a small subdiagonal element to see if we can split it
    lastm = iend
    for m=is:iend-1
      tst = abs(e[m])^2
      if tst <= (eps2*abs(d[m]))*abs(d[m+1])+safmin
        lastm = m
        break
      end
    end
    if lastm < iend # found a splitting element
      e[lastm] = zero(T)
    end

    if lastm == is
      # we have an isolated eigenvalue
      is += 1
    elseif lastm == is+1
      # we have a little 2x2 pair
      rt1,rt2,c,s = _dlaev2(d[is], e[is], d[is+1])
      #rt1,rt2,c,s = _dlaev2_blas(d[is], e[is], d[is+1])
      #rt1,rt2,c,s = _checked_dlaev2(d[is], e[is], d[is+1])
      # note that these are here in dstqrb 
      work[is] = c
      work[is+n-1] = s
      if isvector 
        tst = z[is+1]
        # apply transformation to eigenvector
        z[is+1] = c*tst - s*z[is]
        z[is] = s*tst + c*z[is]
      else
        #println("Split piece")
        _dlasr_right_side_variable_pivot_backward!(size(z,1), 2, 
          @view(work[is:is]), 
          @view(work[is+n-1:is+n-1]), 
          @view(z[:,is:is+1]))
      end 
      # save eigenvalue info
      d[is] = rt1
      d[is+1] = rt2
      e[is] = zero(T) # note that e[is+1=lastm] was already set to zero
      is += 2
    else
      if jtot == maxit
        # failure case handled in jtot == maxit below...
        break
      end
      jtot += 1

      # form shift.
      p = d[is] # save the pivot element
      g = (d[is+1]-p) / (2*e[is])
      r = _dlapy2_julia(g, one(T))
      g = d[lastm] - p + (e[is] / (g+copysign(r, g)))
      s = one(T)
      c = one(T)
      p = zero(T)

      # inner-loop
      for i=reverse(is:lastm-1)
        f = s*e[i]
        b = c*e[i]
        c, s, r = plane_rotation(g, f)
        #c, s, r = _dlartg(g, f)
        if i != lastm-1
          e[i+1] = r
        end
        g = d[i+1] - p
        r = (d[i]-g)*s + two*c*b
        p = s*r
        d[i+1] = g+p
        g = c*r-b
        # save stuff for eigenvecs
        work[i] = c
        work[n-1+i] = -s
      end
      #=_dlasr_right_side_variable_pivot_backward!(
        1, length(is:lastm),
        @view(work[is:lastm-1]), @view(work[(n-1).+(is:lastm-1)]),
        @view(z[is:lastm])') # note don't need "1" here...=#
      #_dlasr_rvb_blas!(1, length(is:lastm),
      #  @view(work[is:lastm-1]), @view(work[(n-1).+(is:lastm-1)]),
      #  @view(z[is:lastm]), 1)
      #_checked_dlasr_right_side_variable_pivot_backward!(
      #  1, length(is:lastm),
      #  @view(work[is:lastm-1]), @view(work[(n-1).+(is:lastm-1)]),
      #  @view(z[is:lastm])') # note don't need "1" here...
       #  @view(z[is:lastm]), 1)
      if isvector
        _apply_plane_rotations_right!(adjoint(@view(z[is:lastm])),
           @view(work[is:lastm-1]), @view(work[(n-1).+(is:lastm-1)]); rev=true)
      else
        _apply_plane_rotations_right!(@view(z[:,is:lastm]),
          @view(work[is:lastm-1]), @view(work[(n-1).+(is:lastm-1)]); rev=true)
      end
      d[is] -= p
      e[is] = g
    end
  end
  info = 0
  if jtot == maxit # early termination, count non-converged vectors
    for v in e
      if !iszero(v)
        info += 1
      end
    end
    error(string("failed to converge eigenvalues\n",
      "$maxit iterations for $n-length tridiagonal problem;",
      "$info not converged"))
  end
  return info, jtot
end

function _qr_iteration(d::AbstractVecOrMat{T},
    e::AbstractVecOrMat{T}, z::AbstractVecOrMat{T}, work::AbstractVecOrMat{T};
    maxit::Int = length(d)*30) where T
  eps = Base.eps(T)/2 # dlamch("E") is eps/2
  eps2 = eps^2
  safmin = floatmin(T) # dlamch("S")
  safmax = 1/safmin
  ssfmax = sqrt(safmax)/3
  ssfmin = sqrt(safmin)/eps2
  two = 2*one(T)

  isvector = ndims(z) == 1 

  n = length(d)
  is = length(d) # we work backwards here... so many loops are reversed...
  iend = 1
  jtot = 0

  while iend < is # l != lend
    #println("At itereration $(jtot) is=$is iend=$iend")
    #display(SymTridiagonal(d, e))
    # look for a small superdiagonal element to see if we can split it
    lastm = iend
    for m=reverse(iend+1:is)
      tst = abs(e[m-1])^2
      if tst <= (eps2*abs(d[m]))*abs(d[m-1])+safmin
        lastm = m
        break
      end
    end
    if lastm > iend # found a splitting element
      e[lastm-1] = zero(T)
    end

    if lastm == is
      # we have an isolated eigenvalue
      is -= 1
    elseif lastm == is-1
      # we have a little 2x2 pair
      rt1,rt2,c,s = _dlaev2(d[is-1], e[is-1], d[is])
      #rt1,rt2,c,s = _dlaev2_blas(d[is-1], e[is-1], d[is])
      #rt1,rt2,c,s = _checked_dlaev2(d[is-1], e[is-1], d[is])
      # apply transformation to eigenvector
      if isvector 
        # no need to save work info as this is now irrelevant...
        tst = z[is]
        # apply transformation to eigenvector
        z[is] = c*tst - s*z[is-1]
        z[is-1] = s*tst + c*z[is-1]
      else
        work[lastm] = c
        work[lastm+n-1] = s
        _dlasr_right_side_variable_pivot_forward!(size(z,1), 2, 
          @view(work[lastm:lastm]), 
          @view(work[lastm+n-1:lastm+n-1]), 
          @view(z[:,is-1:is]))
      end 
      # save eigenvalue info
      d[is-1] = rt1
      d[is] = rt2
      e[is-1] = zero(T) # note that e[is+1=lastm] was already set to zero
      is -= 2
    else
      if jtot == maxit
        # failure case handled in jtot == maxit below...
        break
      end
      jtot += 1

      # form shift.
      p = d[is] # save the pivot element
      g = (d[is-1]-p) / (2*e[is-1])
      r = _dlapy2_julia(g, one(T))
      g = d[lastm] - p + (e[is-1] / (g+copysign(r, g)))
      s = one(T)
      c = one(T)
      p = zero(T)

      # inner-loop
      for i=lastm:is-1
        f = s*e[i]
        b = c*e[i]
        c, s, r = plane_rotation(g, f)
        #c, s, r = _dlartg(g, f)
        if i != lastm
          e[i-1] = r
        end
        g = d[i] - p
        r = (d[i+1]-g)*s + two*c*b
        p = s*r
        d[i] = g+p
        g = c*r-b
        # save stuff for eigenvecs
        work[i] = c
        work[n-1+i] = s
      end
      #@show z
      #_dlasr_rvf_blas!(1, length(lastm:is),
      #  @view(work[lastm:is-1]), @view(work[(n-1).+(lastm:is-1)]),
      #  @view(z[lastm:is]), 1)
      #@show z
      if isvector
        _apply_plane_rotations_right!(adjoint(@view(z[lastm:is])),
          @view(work[lastm:is-1]),
          @view(work[(n-1).+(lastm:is-1)]);rev=false) # forwards in qr
      else
        _apply_plane_rotations_right!(@view(z[:,lastm:is]),
          @view(work[lastm:is-1]),
          @view(work[(n-1).+(lastm:is-1)]);rev=false) # forwards in qr
      end 
      d[is] -= p
      e[is-1] = g
    end
  end
  info = 0
  if jtot == maxit # early termination, count non-converged vectors
    for v in e
      if !iszero(v)
        info += 1
      end
    end
    error(string("failed to converge eigenvalues\n",
      "$maxit iterations for $n-length tridiagonal problem;",
      "$info not converged"))
  end
  return info, jtot
end
