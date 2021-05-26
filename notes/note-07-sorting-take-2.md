A small project to do a pure port of ARPACK Symmetric solver to Julia - Part 7 - Sorting Take 2

Last post was a disaster. I tried to do something quick and simple in Julia
and ran into all sorts of issues with things I didn't quite understand and
breaks in my mental model about how Julia works. This was just to try and
write the `real` julia code for the following task: sort an array and take
another along for the ride as is done in the following code.

This post is straightforward. We are simply going to port this code.

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

This code has the same routine four times with different compare functions.
Let's pull out one routine.

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

Things to note. First, this uses Fortran's flexible indexing to switch to
C-based indices. While this is possible in Julia, we are just going to
do it manually. Again, we are trying to keep this simple. (Despite the
last blog post.)

First step is just to unwrap the loops and GOTOs. In this case, I find
an outer to inner strategy is helpful.

The outer loop (line number 10) seems to be the following

    while igap != 0
      # stuff
      igap = div(igap,2)
    end

The next loops (line number 20 and 30) seems to be encoding the pattern.

    for i=igap:(n-1)
      j = i-igap
      while j > 0
        if cmp(j,j+igap)
          # swap
        else
          break # out of while loop
        end
        j = j-igap
      end
    end

Now, we assemble these into a simple function. Other changes made are:
* Predefined strings become symbols.
* Rather than repeating the code 4 times, we use functions to
  handle the comparison
* Add a new macro to handle checking the length of input vectors.




      """
      Sort the array X1 in the order specified by WHICH and optionally
      applies the permutation to the array X2.

      Usage:
      call dsortr
         ( WHICH, APPLY, N, X1, X2 )

      Arguments
      WHICH   Symbol  (Input)
              :LM -> X1 is sorted into increasing order of magnitude.
              :SM -> X1 is sorted into decreasing order of magnitude.
              :LA -> X1 is sorted into increasing order of algebraic.
              :SA -> X1 is sorted into decreasing order of algebraic.

      APPLY   Boolean.  (Input)
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

        while igap != 0
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
                break # go to 30, which means to exit while loop
              end
              j = j - igap
            end
          end
          igap = div(igap, 2)
        end
      end


The macro for parameter checking is pretty simple, although it took me
about an hour to get all the various options correct.

    macro jl_arpack_check_length(var,len)
      varstr = string(var)
      lenstr = string(len)
      errstr1 = "range 1:$lenstr="
      errstr2 = " out of bounds for $varstr of length "

      return esc( quote
        if $len > length($var)
          throw(ArgumentError(string($errstr1, $len, $errstr2, length($var))))
        end
      end)
    end

To get the "string" command, I had to look at

    x = "Hi"
    mytest(x) = "My string $x"
    @code_lowered mytest(x)

To find how the string interpolation actually gets map to the string command.    

Some history
============

There is a comment in the code that references LANSO. Here is the original
routine from LANSO (found in `svdpack` by Mike Berry).

C
C @(#)DSORT2.F	3.2 (BNP) 12/9/88
C
      SUBROUTINE DSORT2(N,ARRAY1,ARRAY2)
      INTEGER N
      REAL*8 ARRAY1(0:N-1),ARRAY2(0:N-1)
C
C.... SORT ARRAY1 AND ARRAY2 INTO INCREASING ORDER FOR ARRAY1
C
      INTEGER IGAP,I,J
      REAL*8 TEMP
C
      IGAP = N/2
 10   IF (IGAP.GT.0) THEN
        DO 200 I = IGAP,N-1
          J = I-IGAP
 50       IF (J.GE.0.AND.ARRAY1(J).GT.ARRAY1(J+IGAP)) THEN
              TEMP = ARRAY1(J)
              ARRAY1(J) = ARRAY1(J+IGAP)
              ARRAY1(J+IGAP) = TEMP
              TEMP = ARRAY2(J)
              ARRAY2(J) = ARRAY2(J+IGAP)
              ARRAY2(J+IGAP) = TEMP
          ELSE
            GO TO 200
          ENDIF
          J = J-IGAP
          GO TO 50
 200    CONTINUE
      ELSE
        RETURN
      ENDIF
      IGAP = IGAP/2
      GO TO 10
      END
