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
            ritz::Vector{Float64},
            bounds::Vector{Float64},
            tol::Float64;
            stats=nothing,
            debug=nothing
  )
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

#=
c     | Executable Statements |
=#
  t0 = @jl_arpack_time

  if n > length(ritz)
    throw(ArgumentError("range 1:n=$n out of bounds for ritz, of length $(length(ritz))"))
  end
  if n > length(bounds)
    throw(ArgumentError("range 1:n=$n out of bounds for bounds, of length $(length(bounds))"))
  end

  eps23::Float64 = eps(1.0)^(2.0/3.0)
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
  x1::Vector{Float64}, # Input/Output
  x2::Vector{Float64}, # Input/Output
  )

  @jl_arpack_check_length(x1, n)
  if apply
    @jl_arpack_check_length(x2, n)
  end

  igap = div(n,2) # integer division

  # the ARPACK sorting routine is all the same, just changes the
  # comparison, Julia gives us a way of stating that very easily.
  if which == :SA
    cmp = (x::Float64, y::Float64) -> (x < y)
  elseif which == :SM
    cmp = (x::Float64, y::Float64) -> (abs(x) < abs(y))
  elseif which == :LA
    cmp = (x::Float64, y::Float64) -> (x > y)
  elseif which == :LM
    cmp = (x::Float64, y::Float64) -> (abs(x) > abs(y))
  end

  @inbounds while igap != 0
    for i=igap:n-1
      j = i-igap
      while j >= 0
        # in the fortran code, it uses indices 0, n-1
        if cmp(x1[j+1],x1[j+igap+1])
          x1[j+1],x1[j+igap+1] = x1[j+igap+1],x1[j+1]
          if apply
            x2[j+1],x2[j+igap+1] = x2[j+igap+1],x2[j+1]
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
