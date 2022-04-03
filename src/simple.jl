"""
Convergence testing for the symmetric Arnoldi eigenvalue routine.

In Julia, this can be replaced with
  sum( )

Julia port of ARPACK dsconv routine.
Usage:
 call dsconv
    ( N, RITZ, BOUNDS, TOL; debug, stats )

Arguments
 N       Integer.  (INPUT)
         Number of Ritz values to check for convergence.

 RITZ    Double precision array of length N.  (INPUT)
         The Ritz values to be checked for convergence.

 BOUNDS  Double precision array of length N.  (INPUT)
         Ritz estimates associated with the Ritz values in RITZ.

 TOL     Double precision scalar.  (INPUT)
         Desired relative accuracy for a Ritz value to be considered
         "converged".

Optional values (they default to nothing)
 debug (INPUT/OUTPUT) nothing or an ArpackStats struct
 stats (INPUT/OUTPUT) nothing or an ArpackStats struct
These are used to collect debugging and performance stats and are
mostly internal.


Return value
 NCONV   Integer scalar.  (RETURN)
         Number of "converged" Ritz values.
"""
function dsconv(
            n::Int,
            ritz::AbstractVector{T},
            bounds::AbstractVector{T},
            tol::T;
            stats=nothing,
            debug=nothing
  ) where T
#=
  c\Author
  c     Danny Sorensen               Phuong Vu
  c     Richard Lehoucq              CRPC / Rice University
  c     Dept. of Computational &     Houston, Texas
  c     Applied Mathematics
  c     Rice University
  c     Houston, Texas

  Julia port by David F. Gleich, Purdue University
=#

  # c     | Executable Statements |
  t0 = @jl_arpack_time

  @jl_arpack_check_length(ritz, n)
  @jl_arpack_check_length(bounds, n)
  
  eps23::T = _eps23(T) # dlamch("E")**(2.0D+0 / 3.0D+0)
  nconv::Int = 0
  @inbounds for i=1:n
    #=
    c        | The i-th Ritz value is considered "converged"       |
    c        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   |
    =#
    tmp = max( eps23, abs(ritz[i]) )
    if ( bounds[i] <= tol*tmp )
      nconv += 1
    end
  end

  #  call arscnd (t1)
  #  tsconv = tsconv + (t1 - t0)
  @jl_update_time(tconv, t0)

  return nconv
end

#= ORIGINAL CODE
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsconv
c
c\Description:
c  Convergence testing for the symmetric Arnoldi eigenvalue routine.
c
c\Usage:
c  call dsconv
c     ( N, RITZ, BOUNDS, TOL, NCONV )
c
c\Arguments
c  N       Integer.  (INPUT)
c          Number of Ritz values to check for convergence.
c
c  RITZ    Double precision array of length N.  (INPUT)
c          The Ritz values to be checked for convergence.
c
c  BOUNDS  Double precision array of length N.  (INPUT)
c          Ritz estimates associated with the Ritz values in RITZ.
c
c  TOL     Double precision scalar.  (INPUT)
c          Desired relative accuracy for a Ritz value to be considered
c          "converged".
c
c  NCONV   Integer scalar.  (OUTPUT)
c          Number of "converged" Ritz values.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Routines called:
c     arscnd  ARPACK utility routine for timing.
c     dlamch  LAPACK routine that determines machine constants.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\SCCS Information: @(#)
c FILE: sconv.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
c
c\Remarks
c     1. Starting with version 2.4, this routine no longer uses the
c        Parlett strategy using the gap conditions.
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsconv (n, ritz, bounds, tol, nconv)
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    n, nconv
      Double precision
     &           tol
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           ritz(n), bounds(n)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i
      Double precision
     &           temp, eps23
c
c     %-------------------%
c     | External routines |
c     %-------------------%
c
      Double precision
     &           dlamch
      external   dlamch

c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    abs
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      call arscnd (t0)
c
      eps23 = dlamch('Epsilon-Machine')
      eps23 = eps23**(2.0D+0 / 3.0D+0)
c
      nconv  = 0
      do 10 i = 1, n
c
c        %-----------------------------------------------------%
c        | The i-th Ritz value is considered "converged"       |
c        | when: bounds(i) .le. TOL*max(eps23, abs(ritz(i)))   |
c        %-----------------------------------------------------%
c
         temp = max( eps23, abs(ritz(i)) )
         if ( bounds(i) .le. tol*temp ) then
            nconv = nconv + 1
         end if
c
   10 continue
c
      call arscnd (t1)
      tsconv = tsconv + (t1 - t0)
c
      return
c
c     %---------------%
c     | End of dsconv |
c     %---------------%
c
      end
=#

"""
Sort the array X1 in the order specified by WHICH and optionally
applies the permutation to the array X2.

Usage:
call dsortr
   ( WHICH, APPLY, N, X1, X2 )

Arguments
WHICH   Character*2.  (Input)
        'LM' -> X1 is sorted into increasing order of magnitude.
        'SM' -> X1 is sorted into decreasing order of magnitude.
        'LA' -> X1 is sorted into increasing order of algebraic.
        'SA' -> X1 is sorted into decreasing order of algebraic.

APPLY   Logical.  (Input)
        APPLY = .TRUE.  -> apply the sorted order to X2.
        APPLY = .FALSE. -> do not apply the sorted order to X2.

N       Integer.  (INPUT)
        Size of the arrays.

X1      Double precision array of length N.  (INPUT/OUTPUT)
        The array to be sorted.

X2      Double precision array of length N.  (INPUT/OUTPUT)
        Only referenced if APPLY = .TRUE.

Returns: nothing

"""
function dsortr(
  which::Symbol, # Input
  apply::Bool, # Input
  n::Int, # Input
  x1::AbstractVector{T}, # Input/Output
  x2::AbstractVecOrMat{T}, # Input/Output
) where {T <: AbstractFloat}

  @jl_arpack_check_length(x1, n)
  if apply
    @jl_arpack_check_length(x2, n)
  end

  # the ARPACK sorting routine is all the same, just changes the
  # comparison, Julia gives us a way of stating that very easily.
  if which == :SA    
    _dsortr_loop((x::T, y::T) -> (x < y), apply, n, x1, x2)
  elseif which == :SM
    _dsortr_loop((x::T, y::T) -> (abs(x) < abs(y)), apply, n, x1, x2)
  elseif which == :LA
    _dsortr_loop((x::T, y::T) -> (x > y), apply, n, x1, x2)
  elseif which == :LM
    _dsortr_loop((x::T, y::T) -> (abs(x) > abs(y)), apply, n, x1, x2)
  end
end

function _dsortr_loop(cmp, apply, n::Int, x1, x2)
  igap = div(n,2) # integer division
  @inbounds while igap != 0
    for i=igap:n-1
      j = i-igap
      while j >= 0
        # in the fortran code, it uses indices 0, n-1
        if cmp(x1[j+1],x1[j+igap+1])
          x1[j+1],x1[j+igap+1] = x1[j+igap+1],x1[j+1]
          if apply
            _sort_swap(x2, j+1, j+igap+1)
          end
        else
          break # go to 30, means go to continue branch
        end
        j = j - igap
      end
    end
    igap = div(igap, 2)
  end
end

function _sort_swap(x2::AbstractVector, i1::Integer, i2::Integer) 
  x2[i1],x2[i2] = x2[i2],x2[i1]
end 
function _sort_swap(x2::AbstractMatrix, i1::Integer, i2::Integer)  # swap columns
  @simd for i in 1:size(x2,1)
    x2[i,i1],x2[i,i2] = x2[i,i2],x2[i,i1]
  end 
end 
    
##
#=
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsortr
c
c\Description:
c  Sort the array X1 in the order specified by WHICH and optionally
c  applies the permutation to the array X2.
c
c\Usage:
c  call dsortr
c     ( WHICH, APPLY, N, X1, X2 )
c
c\Arguments
c  WHICH   Character*2.  (Input)
c          'LM' -> X1 is sorted into increasing order of magnitude.
c          'SM' -> X1 is sorted into decreasing order of magnitude.
c          'LA' -> X1 is sorted into increasing order of algebraic.
c          'SA' -> X1 is sorted into decreasing order of algebraic.
c
c  APPLY   Logical.  (Input)
c          APPLY = .TRUE.  -> apply the sorted order to X2.
c          APPLY = .FALSE. -> do not apply the sorted order to X2.
c
c  N       Integer.  (INPUT)
c          Size of the arrays.
c
c  X1      Double precision array of length N.  (INPUT/OUTPUT)
c          The array to be sorted.
c
c  X2      Double precision array of length N.  (INPUT/OUTPUT)
c          Only referenced if APPLY = .TRUE.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     12/16/93: Version ' 2.1'.
c               Adapted from the sort routine in LANSO.
c
c\SCCS Information: @(#)
c FILE: sortr.F   SID: 2.3   DATE OF SID: 4/19/96   RELEASE: 2
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsortr (which, apply, n, x1, x2)
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      logical    apply
      integer    n
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           x1(0:n-1), x2(0:n-1)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, igap, j
      Double precision
     &           temp
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      igap = n / 2
c
      if (which .eq. 'SA') then
c
c        X1 is sorted into decreasing order of algebraic.
c
   10    continue
         if (igap .eq. 0) go to 9000
         do 30 i = igap, n-1
            j = i-igap
   20       continue
c
            if (j.lt.0) go to 30
c
            if (x1(j).lt.x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 30
            endif
            j = j-igap
            go to 20
   30    continue
         igap = igap / 2
         go to 10
c
      else if (which .eq. 'SM') then
c
c        X1 is sorted into decreasing order of magnitude.
c
   40    continue
         if (igap .eq. 0) go to 9000
         do 60 i = igap, n-1
            j = i-igap
   50       continue
c
            if (j.lt.0) go to 60
c
            if (abs(x1(j)).lt.abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 60
            endif
            j = j-igap
            go to 50
   60    continue
         igap = igap / 2
         go to 40
c
      else if (which .eq. 'LA') then
c
c        X1 is sorted into increasing order of algebraic.
c
   70    continue
         if (igap .eq. 0) go to 9000
         do 90 i = igap, n-1
            j = i-igap
   80       continue
c
            if (j.lt.0) go to 90
c
            if (x1(j).gt.x1(j+igap)) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 90
            endif
            j = j-igap
            go to 80
   90    continue
         igap = igap / 2
         go to 70
c
      else if (which .eq. 'LM') then
c
c        X1 is sorted into increasing order of magnitude.
c
  100    continue
         if (igap .eq. 0) go to 9000
         do 120 i = igap, n-1
            j = i-igap
  110       continue
c
            if (j.lt.0) go to 120
c
            if (abs(x1(j)).gt.abs(x1(j+igap))) then
               temp = x1(j)
               x1(j) = x1(j+igap)
               x1(j+igap) = temp
               if (apply) then
                  temp = x2(j)
                  x2(j) = x2(j+igap)
                  x2(j+igap) = temp
               end if
            else
               go to 120
            endif
            j = j-igap
            go to 110
  120    continue
         igap = igap / 2
         go to 100
      end if
c
 9000 continue
      return
c
c     %---------------%
c     | End of dsortr |
c     %---------------%
c
      end
=#

function _swap_within_array(n::Int, v::AbstractVector{T}, offset::Int, ioffset::Int=1) where T
  @jl_arpack_check_length(v, ioffset+n-1)
  @jl_arpack_check_length(v, n+offset-1)
  for i=1:n #
    tmp=v[offset-1+i]
    v[offset-1+i]=v[i+ioffset-1]
    v[i+ioffset-1] = tmp
  end
end

function _copyn!(n::Int, dst::AbstractVecOrMat{T}, src::AbstractVecOrMat{T}) where T
  #@assert n <= length(dst) && n <= length(src)
  @jl_arpack_check_length(dst, n)
  @jl_arpack_check_length(src, n)
  @inbounds @simd for i=1:n
    dst[i] = src[i]
  end
end

"""
 Given the eigenvalues of the symmetritridiagonal matrix H,
 computes the NP shifts AMU that are zeros of the polynomial of
 degree NP which filters out components of the unwanted eigenvectors
 corresponding to the AMU's based on some given criteria.

 NOTE: This is called even in the case of user specified shifts in
 order to sort the eigenvalues, and error bounds of H for later use.

Usage:
 call dsgets
    ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS, SHIFTS )

Arguments
 ISHIFT  Integer.  (INPUT)
         Method for selecting the implicit shifts at each iteration.
         ISHIFT = 0: user specified shifts
         ISHIFT = 1: exact shift with respect to the matrix H.

 WHICH   Character*2.  (INPUT)
         Shift selection criteria.
         'LM' -> KEV eigenvalues of largest magnitude are retained.
         'SM' -> KEV eigenvalues of smallest magnitude are retained.
         'LA' -> KEV eigenvalues of largest value are retained.
         'SA' -> KEV eigenvalues of smallest value are retained.
         'BE' -> KEV eigenvalues, half from each end of the spectrum.
                 If KEV is odd, compute one more from the high end.

 KEV      Integer.  (INPUT)
         KEV+NP is the size of the matrix H.

 NP      Integer.  (INPUT)
         Number of implicit shifts to be computed.

 RITZ    Double precision array of length KEV+NP.  (INPUT/OUTPUT)
         On INPUT, RITZ contains the eigenvalues of H.
         On OUTPUT, RITZ are sorted so that the unwanted eigenvalues
         are in the first NP locations and the wanted part is in
         the last KEV locations.  When exact shifts are selected, the
         unwanted part corresponds to the shifts to be applied.

 BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
         Error bounds corresponding to the ordering in RITZ.

 SHIFTS  Double precision array of length NP.  (INPUT/OUTPUT)
         On INPUT:  contains the user specified shifts if ISHIFT = 0.
         On OUTPUT: contains the shifts sorted into decreasing order
         of magnitude with respect to the Ritz estimates contained in
         BOUNDS. If ISHIFT = 0, SHIFTS is not modified on exit.
"""
function dsgets(
    ishift::Int,
    which::Symbol,
    kev::Int,
    np::Int,
    ritz::AbstractVecOrMat{T},
    bounds::AbstractVecOrMat{T},
    shifts::AbstractVecOrMat{T};
    stats::Union{ArpackStats,Nothing}=nothing,
    debug::Union{ArpackDebug,Nothing}=nothing
    ) where T 

  @jl_arpack_check_length(ritz, kev+np)
  @jl_arpack_check_length(bounds, kev+np)
  @jl_arpack_check_length(shifts, np)

  #=
    c     | Initialize timing statistics  |
    c     | & message level for debugging |
  =#
  t0 = @jl_arpack_time
  msglvl = @jl_arpack_debug(mgets,0)

  #=
    c     | Local Scalars |
          integer    kevd2, msglvl
  =#
  if which == :BE
    #=
    c        | Both ends of the spectrum are requested.            |
    c        | Sort the eigenvalues into algebraically increasing  |
    c        | order first then swap high end of the spectrum next |
    c        | to low end in appropriate locations.                |
    c        | NOTE: when np < floor(kev/2) be careful not to swap |
    c        | overlapping locations.                              |
    =#

    #=
    call dsortr ('LA', .true., kev+np, ritz, bounds)
    kevd2 = kev / 2
    if ( kev .gt. 1 ) then
       call dswap ( min(kevd2,np), ritz, 1,
&                   ritz( max(kevd2,np)+1 ), 1)
       call dswap ( min(kevd2,np), bounds, 1,
&                   bounds( max(kevd2,np)+1 ), 1)
    end if
    =#
    dsortr(:LA, true, kev+np, ritz, bounds)
    kevd2::Int = div(kev,2)
    if kev > 1
      # maps to the dswap command above
      _swap_within_array( min(kevd2,np), ritz, max(kevd2,np)+1)
      _swap_within_array( min(kevd2,np), bounds, max(kevd2,np)+1)
    end
  else # :LM, :SM, :LA, :SA
    #=
    c        | LM, SM, LA, SA case.                               |
    c        | Sort the eigenvalues of H into the desired order   |
    c        | and apply the resulting order to BOUNDS.           |
    c        | The eigenvalues are sorted so that the wanted part |
    c        | are always in the last KEV locations.               |
    =#
    dsortr(which, true, kev+np, ritz, bounds)
  end

  if ishift == 1 && np > 0
    #=
    c        | Sort the unwanted Ritz values used as shifts so that  |
    c        | the ones with largest Ritz estimates are first.       |
    c        | This will tend to minimize the effects of the         |
    c        | forward instability of the iteration when the shifts  |
    c        | are applied in subroutine dsapps.                     |
    =#
    dsortr(:SM, true, np, bounds, ritz )
    _copyn!(np, shifts, ritz)
  end

  @jl_update_time(tgets, t0)

  if msglvl > 0
    #@jl_arpack_ivout()
    @assert debug !== nothing
    println(debug.logfile, "_sgets: KEV is ", kev)
    println(debug.logfile, "_sgets: NP is ", np)
    _arpack_vout(debug,"_sgets: Eigenvalues of current H matrix", @view ritz[1:kev+np])
    _arpack_vout(debug,"_sgets: Associated Ritz estimates", @view bounds[1:kev+np])
  end
    #=
          if (msglvl .gt. 0) then
             call ivout (logfil, 1, [kev], ndigit, '_sgets: KEV is')
             call ivout (logfil, 1, [np], ndigit, '_sgets: NP is')
             call dvout (logfil, kev+np, ritz, ndigit,
         &        '_sgets: Eigenvalues of current H matrix')
             call dvout (logfil, kev+np, bounds, ndigit,
         &        '_sgets: Associated Ritz estimates')
          end if
    c
          return
    =#
end

#=
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsgets
c
c\Description:
c  Given the eigenvalues of the symmetric tridiagonal matrix H,
c  computes the NP shifts AMU that are zeros of the polynomial of
c  degree NP which filters out components of the unwanted eigenvectors
c  corresponding to the AMU's based on some given criteria.
c
c  NOTE: This is called even in the case of user specified shifts in
c  order to sort the eigenvalues, and error bounds of H for later use.
c
c\Usage:
c  call dsgets
c     ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS, SHIFTS )
c
c\Arguments
c  ISHIFT  Integer.  (INPUT)
c          Method for selecting the implicit shifts at each iteration.
c          ISHIFT = 0: user specified shifts
c          ISHIFT = 1: exact shift with respect to the matrix H.
c
c  WHICH   Character*2.  (INPUT)
c          Shift selection criteria.
c          'LM' -> KEV eigenvalues of largest magnitude are retained.
c          'SM' -> KEV eigenvalues of smallest magnitude are retained.
c          'LA' -> KEV eigenvalues of largest value are retained.
c          'SA' -> KEV eigenvalues of smallest value are retained.
c          'BE' -> KEV eigenvalues, half from each end of the spectrum.
c                  If KEV is odd, compute one more from the high end.
c
c  KEV      Integer.  (INPUT)
c          KEV+NP is the size of the matrix H.
c
c  NP      Integer.  (INPUT)
c          Number of implicit shifts to be computed.
c
c  RITZ    Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c          On INPUT, RITZ contains the eigenvalues of H.
c          On OUTPUT, RITZ are sorted so that the unwanted eigenvalues
c          are in the first NP locations and the wanted part is in
c          the last KEV locations.  When exact shifts are selected, the
c          unwanted part corresponds to the shifts to be applied.
c
c  BOUNDS  Double precision array of length KEV+NP.  (INPUT/OUTPUT)
c          Error bounds corresponding to the ordering in RITZ.
c
c  SHIFTS  Double precision array of length NP.  (INPUT/OUTPUT)
c          On INPUT:  contains the user specified shifts if ISHIFT = 0.
c          On OUTPUT: contains the shifts sorted into decreasing order
c          of magnitude with respect to the Ritz estimates contained in
c          BOUNDS. If ISHIFT = 0, SHIFTS is not modified on exit.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dsortr  ARPACK utility sorting routine.
c     ivout   ARPACK utility routine that prints integers.
c     arscnd  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dcopy   Level 1 BLAS that copies one vector to another.
c     dswap   Level 1 BLAS that swaps the contents of two vectors.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University
c     Dept. of Computational &     Houston, Texas
c     Applied Mathematics
c     Rice University
c     Houston, Texas
c
c\Revision history:
c     xx/xx/93: Version ' 2.1'
c
c\SCCS Information: @(#)
c FILE: sgets.F   SID: 2.4   DATE OF SID: 4/19/96   RELEASE: 2
c
c\Remarks
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsgets ( ishift, which, kev, np, ritz, bounds, shifts )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character*2 which
      integer    ishift, kev, np
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           bounds(kev+np), ritz(kev+np), shifts(np)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           one, zero
      parameter (one = 1.0D+0, zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    kevd2, msglvl
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dswap, dcopy, dsortr, arscnd
c
c     %---------------------%
c     | Intrinsic Functions |
c     %---------------------%
c
      intrinsic    max, min
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------%
c
      call arscnd (t0)
      msglvl = msgets
c
      if (which .eq. 'BE') then
c
c        %-----------------------------------------------------%
c        | Both ends of the spectrum are requested.            |
c        | Sort the eigenvalues into algebraically increasing  |
c        | order first then swap high end of the spectrum next |
c        | to low end in appropriate locations.                |
c        | NOTE: when np < floor(kev/2) be careful not to swap |
c        | overlapping locations.                              |
c        %-----------------------------------------------------%
c
         call dsortr ('LA', .true., kev+np, ritz, bounds)
         kevd2 = kev / 2
         if ( kev .gt. 1 ) then
            call dswap ( min(kevd2,np), ritz, 1,
     &                   ritz( max(kevd2,np)+1 ), 1)
            call dswap ( min(kevd2,np), bounds, 1,
     &                   bounds( max(kevd2,np)+1 ), 1)
         end if
c
      else
c
c        %----------------------------------------------------%
c        | LM, SM, LA, SA case.                               |
c        | Sort the eigenvalues of H into the desired order   |
c        | and apply the resulting order to BOUNDS.           |
c        | The eigenvalues are sorted so that the wanted part |
c        | are always in the last KEV locations.               |
c        %----------------------------------------------------%
c
         call dsortr (which, .true., kev+np, ritz, bounds)
      end if
c
      if (ishift .eq. 1 .and. np .gt. 0) then
c
c        %-------------------------------------------------------%
c        | Sort the unwanted Ritz values used as shifts so that  |
c        | the ones with largest Ritz estimates are first.       |
c        | This will tend to minimize the effects of the         |
c        | forward instability of the iteration when the shifts  |
c        | are applied in subroutine dsapps.                     |
c        %-------------------------------------------------------%
c
         call dsortr ('SM', .true., np, bounds, ritz)
         call dcopy (np, ritz, 1, shifts, 1)
      end if
c
      call arscnd (t1)
      tsgets = tsgets + (t1 - t0)
c
      if (msglvl .gt. 0) then
         call ivout (logfil, 1, [kev], ndigit, '_sgets: KEV is')
         call ivout (logfil, 1, [np], ndigit, '_sgets: NP is')
         call dvout (logfil, kev+np, ritz, ndigit,
     &        '_sgets: Eigenvalues of current H matrix')
         call dvout (logfil, kev+np, bounds, ndigit,
     &        '_sgets: Associated Ritz estimates')
      end if
c
      return
c
c     %---------------%
c     | End of dsgets |
c     %---------------%
c
      end
=#


"""
call dseigt
( RNORM, N, H, LDH, EIG, BOUNDS, WORKL, IERR )

Arguments
RNORM   Double precision scalar.  (INPUT)
     RNORM contains the residual norm corresponding to the current
     symmetri tridiagonal matrix H.

N       Integer.  (INPUT)
     Size of the symmetri tridiagonal matrix H.

H       Double precision N by 2 array.  (INPUT)
     H contains the symmetri tridiagonal matrix with the 
     subdiagonal in the first column starting at H(2,1) and the 
     main diagonal in second column.

LDH     Integer.  (INPUT)
     Leading dimension of H exactly as declared in the calling 
     program.

EIG     Double precision array of length N.  (OUTPUT)
     On output, EIG contains the N eigenvalues of H possibly 
     unsorted.  The BOUNDS arrays are returned in the
     same sorted order as EIG.

BOUNDS  Double precision array of length N.  (OUTPUT)
     On output, BOUNDS contains the error estimates corresponding
     to the eigenvalues EIG.  This is equal to RNORM times the
     last components of the eigenvectors corresponding to the
     eigenvalues in EIG.

WORKL   Double precision work array of length 3*N.  (WORKSPACE)
     Private (replicated) array on each PE or array allocated on
     the front end.

Return value
------------
IERR    Integer.  (OUTPUT)
     Error exit flag from dstqrb.
"""
function dseigt!(
  rnorm::T,
  n::Int,
  H::AbstractMatrix{T},
  ldh::Int,
  eig::AbstractVecOrMat{T},
  bounds::AbstractVecOrMat{T},
  workl::AbstractVecOrMat{T},
  state::Union{AbstractArpackState{T},Nothing}
  ; 
  stats::Union{ArpackStats,Nothing}=nothing,
  debug::Union{ArpackDebug,Nothing}=nothing,
) where T 
  
  @jl_arpack_check_length(eig, n)
  @jl_arpack_check_length(bounds, n)
  @jl_arpack_check_length(workl, 3n)
  @jl_arpack_check_size(H,n,2)
  @jl_arpack_check_size(H,ldh,2) # correct input size...

  # c     | Initialize timing statistics  |
  # c     | & message level for debugging |
  t0 = @jl_arpack_time()
  msglvl = @jl_arpack_debug(meigt,0)

  if msglvl > 0
    _arpack_vout(debug, "_seigt: main diagonal of matrix H", @view(H[1:n,2]))
    if n > 1
      _arpack_vout(debug, "_seigt: sub diagonal of matrix H", @view(H[2:n,1]))
    end 
  end

  # copyto!(@view(eig[1:n]), @view(H[1:n,2]))
  # copyto!(@view(workl[1:n-1]), @view(H[2:n,1]))
  # Julia's conservative nature detects that these arrays may 
  # overlap because H is a reshaped version of a view into workl
  _copyn!(n, @view(eig[1:n]), @view(H[1:n,2]))
  _copyn!(n-1, @view(workl[1:n-1]), @view(H[2:n,1]))
  # 
  
  dstqrb!(
    n, @view(eig[1:n]), @view(workl[1:n-1]), @view(bounds[1:n]), @view(workl[n+1:3n]), state
  )
  # Julia note:
  # We have dstqrb! throw a hard error if the eigenvalue 
  # routine fails as this is a very serious problem that
  # really should cause a bigger error than just a bad
  # return value. So this is a slight deviation from Arpack
  ierr = 0 
  if ierr != 0 
    # failure, don't do anything and just return 
  else
    if msglvl > 1
      _arpack_vout(debug, "_seigt: last row of the eigenvector matrix for H", 
        @view(bounds[1:n]))
    end
  end
  # c     | Finally determine the error bounds associated |
  # c     | with the n Ritz values of H.                  |
  for k=1:n
    bounds[k] = rnorm*abs(bounds[k])
  end

  # note, arpack doesn't update this on a failure, but I think we should...
  @jl_update_time(teigt, t0)
  return ierr
end 

#=
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dseigt
c
c\Description: 
c  Compute the eigenvalues of the current symmetric tridiagonal matrix
c  and the corresponding error bounds given the current residual norm.
c
c\Usage:
c  call dseigt
c     ( RNORM, N, H, LDH, EIG, BOUNDS, WORKL, IERR )
c
c\Arguments
c  RNORM   Double precision scalar.  (INPUT)
c          RNORM contains the residual norm corresponding to the current
c          symmetric tridiagonal matrix H.
c
c  N       Integer.  (INPUT)
c          Size of the symmetric tridiagonal matrix H.
c
c  H       Double precision N by 2 array.  (INPUT)
c          H contains the symmetric tridiagonal matrix with the 
c          subdiagonal in the first column starting at H(2,1) and the 
c          main diagonal in second column.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling 
c          program.
c
c  EIG     Double precision array of length N.  (OUTPUT)
c          On output, EIG contains the N eigenvalues of H possibly 
c          unsorted.  The BOUNDS arrays are returned in the
c          same sorted order as EIG.
c
c  BOUNDS  Double precision array of length N.  (OUTPUT)
c          On output, BOUNDS contains the error estimates corresponding
c          to the eigenvalues EIG.  This is equal to RNORM times the
c          last components of the eigenvectors corresponding to the
c          eigenvalues in EIG.
c
c  WORKL   Double precision work array of length 3*N.  (WORKSPACE)
c          Private (replicated) array on each PE or array allocated on
c          the front end.
c
c  IERR    Integer.  (OUTPUT)
c          Error exit flag from dstqrb.
c
c\EndDoc
c
c-----------------------------------------------------------------------
c
c\BeginLib
c
c\Local variables:
c     xxxxxx  real
c
c\Routines called:
c     dstqrb  ARPACK routine that computes the eigenvalues and the
c             last components of the eigenvectors of a symmetric
c             and tridiagonal matrix.
c     arscnd  ARPACK utility routine for timing.
c     dvout   ARPACK utility routine that prints vectors.
c     dcopy   Level 1 BLAS that copies one vector to another.
c
c\Author
c     Danny Sorensen               Phuong Vu
c     Richard Lehoucq              CRPC / Rice University 
c     Dept. of Computational &     Houston, Texas 
c     Applied Mathematics
c     Rice University           
c     Houston, Texas            
c
c\Revision history:
c     xx/xx/92: Version ' 2.4'
c
c\SCCS Information: @(#) 
c FILE: seigt.F   SID: 2.4   DATE OF SID: 8/27/96   RELEASE: 2
c
c\Remarks
c     None
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dseigt 
     &   ( rnorm, n, h, ldh, eig, bounds, workl, ierr )
c
c     %----------------------------------------------------%
c     | Include files for debugging and timing information |
c     %----------------------------------------------------%
c
      include   'debug.h'
      include   'stat.h'
c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      integer    ierr, ldh, n
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      Double precision
     &           eig(n), bounds(n), h(ldh,2), workl(3*n)
c
c     %------------%
c     | Parameters |
c     %------------%
c
      Double precision
     &           zero
      parameter (zero = 0.0D+0)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer    i, k, msglvl
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   dcopy, dstqrb, dvout, arscnd
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
c     %-------------------------------%
c     | Initialize timing statistics  |
c     | & message level for debugging |
c     %-------------------------------% 
c
      call arscnd (t0)
      msglvl = mseigt
c
      if (msglvl .gt. 0) then
         call dvout (logfil, n, h(1,2), ndigit,
     &              '_seigt: main diagonal of matrix H')
         if (n .gt. 1) then
         call dvout (logfil, n-1, h(2,1), ndigit,
     &              '_seigt: sub diagonal of matrix H')
         end if
      end if
c
      call dcopy  (n, h(1,2), 1, eig, 1)
      call dcopy  (n-1, h(2,1), 1, workl, 1)
      call dstqrb (n, eig, workl, bounds, workl(n+1), ierr)
      if (ierr .ne. 0) go to 9000
      if (msglvl .gt. 1) then
         call dvout (logfil, n, bounds, ndigit,
     &              '_seigt: last row of the eigenvector matrix for H')
      end if
c
c     %-----------------------------------------------%
c     | Finally determine the error bounds associated |
c     | with the n Ritz values of H.                  |
c     %-----------------------------------------------%
c
      do 30 k = 1, n
         bounds(k) = rnorm*abs(bounds(k))
   30 continue
c 
      call arscnd (t1)
      tseigt = tseigt + (t1 - t0)
c
 9000 continue
      return
c
c     %---------------%
c     | End of dseigt |
c     %---------------%
c
      end
=#