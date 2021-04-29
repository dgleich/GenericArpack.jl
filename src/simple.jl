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
