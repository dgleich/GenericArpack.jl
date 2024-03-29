#=
c-----------------------------------------------------------------------
c\BeginDoc
c
c\Name: dsaitr
c
c\Description:
c  Reverse communication interface for applying NP additional steps to
c  a K step symmetric Arnoldi factorization.
c
c  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T
c
c          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.
c
c  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T
c
c          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.
c
c  where OP and B are as in dsaupd.  The B-norm of r_{k+p} is also
c  computed and returned.
c
c\Usage:
c  call dsaitr
c     ( IDO, BMAT, N, K, NP, MODE, RESID, RNORM, V, LDV, H, LDH,
c       IPNTR, WORKD, INFO )
c
c\Arguments
c  IDO     Integer.  (INPUT/OUTPUT)
c          Reverse communication flag.
c          -------------------------------------------------------------
c          IDO =  0: first call to the reverse communication interface
c          IDO = -1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c                    This is for the restart phase to force the new
c                    starting vector into the range of OP.
c          IDO =  1: compute  Y = OP * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y,
c                    IPNTR(3) is the pointer into WORK for B * X.
c          IDO =  2: compute  Y = B * X  where
c                    IPNTR(1) is the pointer into WORK for X,
c                    IPNTR(2) is the pointer into WORK for Y.
c          IDO = 99: done
c          -------------------------------------------------------------
c          When the routine is used in the "shift-and-invert" mode, the
c          vector B * Q is already available and does not need to be
c          recomputed in forming OP * Q.
c
c  BMAT    Character*1.  (INPUT)
c          BMAT specifies the type of matrix B that defines the
c          semi-inner product for the operator OP.  See dsaupd.
c          B = 'I' -> standard eigenvalue problem A*x = lambda*x
c          B = 'G' -> generalized eigenvalue problem A*x = lambda*M*x
c
c  N       Integer.  (INPUT)
c          Dimension of the eigenproblem.
c
c  K       Integer.  (INPUT)
c          Current order of H and the number of columns of V.
c
c  NP      Integer.  (INPUT)
c          Number of additional Arnoldi steps to take.
c
c  MODE    Integer.  (INPUT)
c          Signifies which form for "OP". If MODE=2 then
c          a reduction in the number of B matrix vector multiplies
c          is possible since the B-norm of OP*x is equivalent to
c          the inv(B)-norm of A*x.
c
c  RESID   Double precision array of length N.  (INPUT/OUTPUT)
c          On INPUT:  RESID contains the residual vector r_{k}.
c          On OUTPUT: RESID contains the residual vector r_{k+p}.
c
c  RNORM   Double precision scalar.  (INPUT/OUTPUT)
c          On INPUT the B-norm of r_{k}.
c          On OUTPUT the B-norm of the updated residual r_{k+p}.
c
c  V       Double precision N by K+NP array.  (INPUT/OUTPUT)
c          On INPUT:  V contains the Arnoldi vectors in the first K
c          columns.
c          On OUTPUT: V contains the new NP Arnoldi vectors in the next
c          NP columns.  The first K columns are unchanged.
c
c  LDV     Integer.  (INPUT)
c          Leading dimension of V exactly as declared in the calling
c          program.
c
c  H       Double precision (K+NP) by 2 array.  (INPUT/OUTPUT)
c          H is used to store the generated symmetric tridiagonal matrix
c          with the subdiagonal in the first column starting at H(2,1)
c          and the main diagonal in the second column.
c
c  LDH     Integer.  (INPUT)
c          Leading dimension of H exactly as declared in the calling
c          program.
c
c  IPNTR   Integer array of length 3.  (OUTPUT)
c          Pointer to mark the starting locations in the WORK for
c          vectors used by the Arnoldi iteration.
c          -------------------------------------------------------------
c          IPNTR(1): pointer to the current operand vector X.
c          IPNTR(2): pointer to the current result vector Y.
c          IPNTR(3): pointer to the vector B * X when used in the
c                    shift-and-invert mode.  X is the current operand.
c          -------------------------------------------------------------
c
c  WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
c          Distributed array to be used in the basic Arnoldi iteration
c          for reverse communication.  The calling program should not
c          use WORKD as temporary workspace during the iteration !!!!!!
c          On INPUT, WORKD(1:N) = B*RESID where RESID is associated
c          with the K step Arnoldi factorization. Used to save some
c          computation at the first step.
c          On OUTPUT, WORKD(1:N) = B*RESID where RESID is associated
c          with the K+NP step Arnoldi factorization.
c
c  INFO    Integer.  (OUTPUT)
c          = 0: Normal exit.
c          > 0: Size of an invariant subspace of OP is found that is
c               less than K + NP.
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
c     dgetv0  ARPACK routine to generate the initial vector.
c     ivout   ARPACK utility routine that prints integers.
c     dmout   ARPACK utility routine that prints matrices.
c     dvout   ARPACK utility routine that prints vectors.
c     dlamch  LAPACK routine that determines machine constants.
c     dlascl  LAPACK routine for careful scaling of a matrix.
c     dgemv   Level 2 BLAS routine for matrix vector multiplication.
c     daxpy   Level 1 BLAS that computes a vector triad.
c     dscal   Level 1 BLAS that scales a vector.
c     dcopy   Level 1 BLAS that copies one vector to another .
c     ddot    Level 1 BLAS that computes the scalar product of two vectors.
c     dnrm2   Level 1 BLAS that computes the norm of a vector.
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
c     xx/xx/93: Version ' 2.4'
c
c\SCCS Information: @(#)
c FILE: saitr.F   SID: 2.6   DATE OF SID: 8/28/96   RELEASE: 2
c
c\Remarks
c  The algorithm implemented is:
c
c  restart = .false.
c  Given V_{k} = [v_{1}, ..., v_{k}], r_{k};
c  r_{k} contains the initial residual vector even for k = 0;
c  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already
c  computed by the calling program.
c
c  betaj = rnorm ; p_{k+1} = B*r_{k} ;
c  For  j = k+1, ..., k+np  Do
c     1) if ( betaj < tol ) stop or restart depending on j.
c        if ( restart ) generate a new starting vector.
c     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];
c        p_{j} = p_{j}/betaj
c     3) r_{j} = OP*v_{j} where OP is defined as in dsaupd
c        For shift-invert mode p_{j} = B*v_{j} is already available.
c        wnorm = || OP*v_{j} ||
c     4) Compute the j-th step residual vector.
c        w_{j} =  V_{j}^T * B * OP * v_{j}
c        r_{j} =  OP*v_{j} - V_{j} * w_{j}
c        alphaj <- j-th component of w_{j}
c        rnorm = || r_{j} ||
c        betaj+1 = rnorm
c        If (rnorm > 0.717*wnorm) accept step and go back to 1)
c     5) Re-orthogonalization step:
c        s = V_{j}'*B*r_{j}
c        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
c        alphaj = alphaj + s_{j};
c     6) Iterative refinement step:
c        If (rnorm1 > 0.717*rnorm) then
c           rnorm = rnorm1
c           accept step and go back to 1)
c        Else
c           rnorm = rnorm1
c           If this is the first time in step 6), go to 5)
c           Else r_{j} lies in the span of V_{j} numerically.
c              Set r_{j} = 0 and rnorm = 0; go to 1)
c        EndIf
c  End Do
c
c\EndLib
c
c-----------------------------------------------------------------------
c
      subroutine dsaitr
     &   (ido, bmat, n, k, np, mode, resid, rnorm, v, ldv, h, ldh,
     &    ipntr, workd, info)
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
      character  bmat*1
      integer    ido, info, k, ldh, ldv, n, mode, np
      Double precision
     &           rnorm
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer    ipntr(3)
      Double precision
     &           h(ldh,2), resid(n), v(ldv,k+np), workd(3*n)
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
      logical    first, orth1, orth2, rstart, step3, step4
      integer    i, ierr, ipj, irj, ivj, iter, itry, j, msglvl,
     &           infol, jj
      Double precision
     &           rnorm1, wnorm, safmin, temp1
      save       orth1, orth2, rstart, step3, step4,
     &           ierr, ipj, irj, ivj, iter, itry, j, msglvl,
     &           rnorm1, safmin, wnorm
c
c     %-----------------------%
c     | Local Array Arguments |
c     %-----------------------%
c
      Double precision
     &           xtemp(2)
c
c     %----------------------%
c     | External Subroutines |
c     %----------------------%
c
      external   daxpy, dcopy, dscal, dgemv, dgetv0, dvout, dmout,
     &           dlascl, ivout, arscnd
c
c     %--------------------%
c     | External Functions |
c     %--------------------%
c
      Double precision
     &           ddot, dnrm2, dlamch
      external   ddot, dnrm2, dlamch
c
c     %-----------------%
c     | Data statements |
c     %-----------------%
c
      data      first / .true. /
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c
      if (first) then
         first = .false.
c
c        %--------------------------------%
c        | safmin = safe minimum is such  |
c        | that 1/sfmin does not overflow |
c        %--------------------------------%
c
         safmin = dlamch('safmin')
      end if
c
      if (ido .eq. 0) then
c
c        %-------------------------------%
c        | Initialize timing statistics  |
c        | & message level for debugging |
c        %-------------------------------%
c
         call arscnd (t0)
         msglvl = msaitr
c
c        %------------------------------%
c        | Initial call to this routine |
c        %------------------------------%
c
         info   = 0
         step3  = .false.
         step4  = .false.
         rstart = .false.
         orth1  = .false.
         orth2  = .false.
c
c        %--------------------------------%
c        | Pointer to the current step of |
c        | the factorization to build     |
c        %--------------------------------%
c
         j      = k + 1
c
c        %------------------------------------------%
c        | Pointers used for reverse communication  |
c        | when using WORKD.                        |
c        %------------------------------------------%
c
         ipj    = 1
         irj    = ipj   + n
         ivj    = irj   + n
      end if
c
c     %-------------------------------------------------%
c     | When in reverse communication mode one of:      |
c     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
c     | will be .true.                                  |
c     | STEP3: return from computing OP*v_{j}.          |
c     | STEP4: return from computing B-norm of OP*v_{j} |
c     | ORTH1: return from computing B-norm of r_{j+1}  |
c     | ORTH2: return from computing B-norm of          |
c     |        correction to the residual vector.       |
c     | RSTART: return from OP computations needed by   |
c     |         dgetv0.                                 |
c     %-------------------------------------------------%
c
      if (step3)  go to 50
      if (step4)  go to 60
      if (orth1)  go to 70
      if (orth2)  go to 90
      if (rstart) go to 30
c
c     %------------------------------%
c     | Else this is the first step. |
c     %------------------------------%
c
c     %--------------------------------------------------------------%
c     |                                                              |
c     |        A R N O L D I     I T E R A T I O N     L O O P       |
c     |                                                              |
c     | Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1) |
c     %--------------------------------------------------------------%
c
 1000 continue
c
         if (msglvl .gt. 2) then
            call ivout (logfil, 1, [j], ndigit,
     &                  '_saitr: generating Arnoldi vector no.')
            call dvout (logfil, 1, [rnorm], ndigit,
     &                  '_saitr: B-norm of the current residual =')
         end if
c
c        %---------------------------------------------------------%
c        | Check for exact zero. Equivalent to determining whether |
c        | a j-step Arnoldi factorization is present.              |
c        %---------------------------------------------------------%
c
         if (rnorm .gt. zero) go to 40
c
c           %---------------------------------------------------%
c           | Invariant subspace found, generate a new starting |
c           | vector which is orthogonal to the current Arnoldi |
c           | basis and continue the iteration.                 |
c           %---------------------------------------------------%
c
            if (msglvl .gt. 0) then
               call ivout (logfil, 1, [j], ndigit,
     &                     '_saitr: ****** restart at step ******')
            end if
c
c           %---------------------------------------------%
c           | ITRY is the loop variable that controls the |
c           | maximum amount of times that a restart is   |
c           | attempted. NRSTRT is used by stat.h         |
c           %---------------------------------------------%
c
            nrstrt = nrstrt + 1
            itry   = 1
   20       continue
            rstart = .true.
            ido    = 0
   30       continue
c
c           %--------------------------------------%
c           | If in reverse communication mode and |
c           | RSTART = .true. flow returns here.   |
c           %--------------------------------------%
c
            call dgetv0 (ido, bmat, itry, .false., n, j, v, ldv,
     &                   resid, rnorm, ipntr, workd, ierr)
            if (ido .ne. 99) go to 9000
            if (ierr .lt. 0) then
               itry = itry + 1
               if (itry .le. 3) go to 20
c
c              %------------------------------------------------%
c              | Give up after several restart attempts.        |
c              | Set INFO to the size of the invariant subspace |
c              | which spans OP and exit.                       |
c              %------------------------------------------------%
c
               info = j - 1
               call arscnd (t1)
               tsaitr = tsaitr + (t1 - t0)
               ido = 99
               go to 9000
            end if
c
   40    continue
c
c        %---------------------------------------------------------%
c        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
c        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
c        | when reciprocating a small RNORM, test against lower    |
c        | machine bound.                                          |
c        %---------------------------------------------------------%
c
         call dcopy (n, resid, 1, v(1,j), 1)
         if (rnorm .ge. safmin) then
             temp1 = one / rnorm
             call dscal (n, temp1, v(1,j), 1)
             call dscal (n, temp1, workd(ipj), 1)
         else
c
c            %-----------------------------------------%
c            | To scale both v_{j} and p_{j} carefully |
c            | use LAPACK routine SLASCL               |
c            %-----------------------------------------%
c
             call dlascl ('General', i, i, rnorm, one, n, 1,
     &                    v(1,j), n, infol)
             call dlascl ('General', i, i, rnorm, one, n, 1,
     &                    workd(ipj), n, infol)
         end if
c
c        %------------------------------------------------------%
c        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
c        | Note that this is not quite yet r_{j}. See STEP 4    |
c        %------------------------------------------------------%
c
         step3 = .true.
         nopx  = nopx + 1
         call arscnd (t2)
         call dcopy (n, v(1,j), 1, workd(ivj), 1)
         ipntr(1) = ivj
         ipntr(2) = irj
         ipntr(3) = ipj
         ido = 1
c
c        %-----------------------------------%
c        | Exit in order to compute OP*v_{j} |
c        %-----------------------------------%
c
         go to 9000
   50    continue
c
c        %-----------------------------------%
c        | Back from reverse communication;  |
c        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}.   |
c        %-----------------------------------%
c
         call arscnd (t3)
         tmvopx = tmvopx + (t3 - t2)
c
         step3 = .false.
c
c        %------------------------------------------%
c        | Put another copy of OP*v_{j} into RESID. |
c        %------------------------------------------%
c
         call dcopy (n, workd(irj), 1, resid, 1)
c
c        %-------------------------------------------%
c        | STEP 4:  Finish extending the symmetric   |
c        |          Arnoldi to length j. If MODE = 2 |
c        |          then B*OP = B*inv(B)*A = A and   |
c        |          we don't need to compute B*OP.   |
c        | NOTE: If MODE = 2 WORKD(IVJ:IVJ+N-1) is   |
c        | assumed to have A*v_{j}.                  |
c        %-------------------------------------------%
c
         if (mode .eq. 2) go to 65
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            step4 = .true.
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c
c           %-------------------------------------%
c           | Exit in order to compute B*OP*v_{j} |
c           %-------------------------------------%
c
            go to 9000
         else if (bmat .eq. 'I') then
              call dcopy(n, resid, 1 , workd(ipj), 1)
         end if
   60    continue
c
c        %-----------------------------------%
c        | Back from reverse communication;  |
c        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j}. |
c        %-----------------------------------%
c
         if (bmat .eq. 'G') then
            call arscnd (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c
         step4 = .false.
c
c        %-------------------------------------%
c        | The following is needed for STEP 5. |
c        | Compute the B-norm of OP*v_{j}.     |
c        %-------------------------------------%
c
   65    continue
         if (mode .eq. 2) then
c
c           %----------------------------------%
c           | Note that the B-norm of OP*v_{j} |
c           | is the inv(B)-norm of A*v_{j}.   |
c           %----------------------------------%
c
            wnorm = ddot (n, resid, 1, workd(ivj), 1)
            wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'G') then
            wnorm = ddot (n, resid, 1, workd(ipj), 1)
            wnorm = sqrt(abs(wnorm))
         else if (bmat .eq. 'I') then
            wnorm = dnrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------%
c        | Compute the j-th residual corresponding |
c        | to the j step factorization.            |
c        | Use Classical Gram Schmidt and compute: |
c        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
c        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
c        %-----------------------------------------%
c
c
c        %------------------------------------------%
c        | Compute the j Fourier coefficients w_{j} |
c        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
c        %------------------------------------------%
c
         if (mode .ne. 2 ) then
            call dgemv('T', n, j, one, v, ldv, workd(ipj), 1, zero,
     &                  workd(irj), 1)
         else if (mode .eq. 2) then
            call dgemv('T', n, j, one, v, ldv, workd(ivj), 1, zero,
     &                  workd(irj), 1)
         end if
c
c        %--------------------------------------%
c        | Orthgonalize r_{j} against V_{j}.    |
c        | RESID contains OP*v_{j}. See STEP 3. |
c        %--------------------------------------%
c
         call dgemv('N', n, j, -one, v, ldv, workd(irj), 1, one,
     &               resid, 1)
c
c        %--------------------------------------%
c        | Extend H to have j rows and columns. |
c        %--------------------------------------%
c
         h(j,2) = workd(irj + j - 1)
         if (j .eq. 1  .or.  rstart) then
            h(j,1) = zero
         else
            h(j,1) = rnorm
         end if
         call arscnd (t4)
c
         orth1 = .true.
         iter  = 0
c
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c
c           %----------------------------------%
c           | Exit in order to compute B*r_{j} |
c           %----------------------------------%
c
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if
   70    continue
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH1 = .true. |
c        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
c        %---------------------------------------------------%
c
         if (bmat .eq. 'G') then
            call arscnd (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c
         orth1 = .false.
c
c        %------------------------------%
c        | Compute the B-norm of r_{j}. |
c        %------------------------------%
c
         if (bmat .eq. 'G') then
            rnorm = ddot (n, resid, 1, workd(ipj), 1)
            rnorm = sqrt(abs(rnorm))
         else if (bmat .eq. 'I') then
            rnorm = dnrm2(n, resid, 1)
         end if
c
c        %-----------------------------------------------------------%
c        | STEP 5: Re-orthogonalization / Iterative refinement phase |
c        | Maximum NITER_ITREF tries.                                |
c        |                                                           |
c        |          s      = V_{j}^T * B * r_{j}                     |
c        |          r_{j}  = r_{j} - V_{j}*s                         |
c        |          alphaj = alphaj + s_{j}                          |
c        |                                                           |
c        | The stopping criteria used for iterative refinement is    |
c        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
c        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
c        | Determine if we need to correct the residual. The goal is |
c        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
c        %-----------------------------------------------------------%
c
         if (rnorm .gt. 0.717*wnorm) go to 100
         nrorth = nrorth + 1
c
c        %---------------------------------------------------%
c        | Enter the Iterative refinement phase. If further  |
c        | refinement is necessary, loop back here. The loop |
c        | variable is ITER. Perform a step of Classical     |
c        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
c        %---------------------------------------------------%
c
   80    continue
c
         if (msglvl .gt. 2) then
            xtemp(1) = wnorm
            xtemp(2) = rnorm
            call dvout (logfil, 2, xtemp, ndigit,
     &           '_saitr: re-orthonalization ; wnorm and rnorm are')
         end if
c
c        %----------------------------------------------------%
c        | Compute V_{j}^T * B * r_{j}.                       |
c        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
c        %----------------------------------------------------%
c
         call dgemv ('T', n, j, one, v, ldv, workd(ipj), 1,
     &               zero, workd(irj), 1)
c
c        %----------------------------------------------%
c        | Compute the correction to the residual:      |
c        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1).  |
c        | The correction to H is v(:,1:J)*H(1:J,1:J) + |
c        | v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j, but only   |
c        | H(j,j) is updated.                           |
c        %----------------------------------------------%
c
         call dgemv ('N', n, j, -one, v, ldv, workd(irj), 1,
     &               one, resid, 1)
c
         if (j .eq. 1  .or.  rstart) h(j,1) = zero
         h(j,2) = h(j,2) + workd(irj + j - 1)
c
         orth2 = .true.
         call arscnd (t2)
         if (bmat .eq. 'G') then
            nbx = nbx + 1
            call dcopy (n, resid, 1, workd(irj), 1)
            ipntr(1) = irj
            ipntr(2) = ipj
            ido = 2
c
c           %-----------------------------------%
c           | Exit in order to compute B*r_{j}. |
c           | r_{j} is the corrected residual.  |
c           %-----------------------------------%
c
            go to 9000
         else if (bmat .eq. 'I') then
            call dcopy (n, resid, 1, workd(ipj), 1)
         end if
   90    continue
c
c        %---------------------------------------------------%
c        | Back from reverse communication if ORTH2 = .true. |
c        %---------------------------------------------------%
c
         if (bmat .eq. 'G') then
            call arscnd (t3)
            tmvbx = tmvbx + (t3 - t2)
         end if
c
c        %-----------------------------------------------------%
c        | Compute the B-norm of the corrected residual r_{j}. |
c        %-----------------------------------------------------%
c
         if (bmat .eq. 'G') then
             rnorm1 = ddot (n, resid, 1, workd(ipj), 1)
             rnorm1 = sqrt(abs(rnorm1))
         else if (bmat .eq. 'I') then
             rnorm1 = dnrm2(n, resid, 1)
         end if
c
         if (msglvl .gt. 0 .and. iter .gt. 0) then
            call ivout (logfil, 1, [j], ndigit,
     &           '_saitr: Iterative refinement for Arnoldi residual')
            if (msglvl .gt. 2) then
                xtemp(1) = rnorm
                xtemp(2) = rnorm1
                call dvout (logfil, 2, xtemp, ndigit,
     &           '_saitr: iterative refinement ; rnorm and rnorm1 are')
            end if
         end if
c
c        %-----------------------------------------%
c        | Determine if we need to perform another |
c        | step of re-orthogonalization.           |
c        %-----------------------------------------%
c
         if (rnorm1 .gt. 0.717*rnorm) then
c
c           %--------------------------------%
c           | No need for further refinement |
c           %--------------------------------%
c
            rnorm = rnorm1
c
         else
c
c           %-------------------------------------------%
c           | Another step of iterative refinement step |
c           | is required. NITREF is used by stat.h     |
c           %-------------------------------------------%
c
            nitref = nitref + 1
            rnorm  = rnorm1
            iter   = iter + 1
            if (iter .le. 1) go to 80
c
c           %-------------------------------------------------%
c           | Otherwise RESID is numerically in the span of V |
c           %-------------------------------------------------%
c
            do 95 jj = 1, n
               resid(jj) = zero
  95        continue
            rnorm = zero
         end if
c
c        %----------------------------------------------%
c        | Branch here directly if iterative refinement |
c        | wasn't necessary or after at most NITER_REF  |
c        | steps of iterative refinement.               |
c        %----------------------------------------------%
c
  100    continue
c
         rstart = .false.
         orth2  = .false.
c
         call arscnd (t5)
         titref = titref + (t5 - t4)
c
c        %----------------------------------------------------------%
c        | Make sure the last off-diagonal element is non negative  |
c        | If not perform a similarity transformation on H(1:j,1:j) |
c        | and scale v(:,j) by -1.                                  |
c        %----------------------------------------------------------%
c
         if (h(j,1) .lt. zero) then
            h(j,1) = -h(j,1)
            if ( j .lt. k+np) then
               call dscal(n, -one, v(1,j+1), 1)
            else
               call dscal(n, -one, resid, 1)
            end if
         end if
c
c        %------------------------------------%
c        | STEP 6: Update  j = j+1;  Continue |
c        %------------------------------------%
c
         j = j + 1
         if (j .gt. k+np) then
            call arscnd (t1)
            tsaitr = tsaitr + (t1 - t0)
            ido = 99
c
            if (msglvl .gt. 1) then
               call dvout (logfil, k+np, h(1,2), ndigit,
     &         '_saitr: main diagonal of matrix H of step K+NP.')
               if (k+np .gt. 1) then
               call dvout (logfil, k+np-1, h(2,1), ndigit,
     &         '_saitr: sub diagonal of matrix H of step K+NP.')
               end if
            end if
c
            go to 9000
         end if
c
c        %--------------------------------------------------------%
c        | Loop back to extend the factorization by another step. |
c        %--------------------------------------------------------%
c
      go to 1000
c
c     %---------------------------------------------------------------%
c     |                                                               |
c     |  E N D     O F     M A I N     I T E R A T I O N     L O O P  |
c     |                                                               |
c     %---------------------------------------------------------------%
c
 9000 continue
      return
c
c     %---------------%
c     | End of dsaitr |
c     %---------------%
c
      end
=#


"""
Usage:
 call dsaitr
    ( IDO, BMAT, N, K, NP, MODE, RESID, RNORM, V, LDV, H, LDH,
      IPNTR, WORKD, INFO )

Arguments
 IDO     Integer.  (INPUT/OUTPUT)
         Reverse communication flag.
         -------------------------------------------------------------
         IDO =  0: first call to the reverse communication interface
         IDO = -1: compute  Y = OP * X  where
                   IPNTR(1) is the pointer into WORK for X,
                   IPNTR(2) is the pointer into WORK for Y.
                   This is for the restart phase to force the new
                   starting vector into the range of OP.
         IDO =  1: compute  Y = OP * X  where
                   IPNTR(1) is the pointer into WORK for X,
                   IPNTR(2) is the pointer into WORK for Y,
                   IPNTR(3) is the pointer into WORK for B * X.
         IDO =  2: compute  Y = B * X  where
                   IPNTR(1) is the pointer into WORK for X,
                   IPNTR(2) is the pointer into WORK for Y.
         IDO = 99: done
         -------------------------------------------------------------
         When the routine is used in the "shift-and-invert" mode, the
         vector B * Q is already available and does not need to be
         recomputed in forming OP * Q.

 BMAT    Character*1.  (INPUT)
         BMAT specifies the type of matrix B that defines the
         semi-inner product for the operator OP.  See dsaupd.
         B = 'I' -> standard eigenvalue problem A*x = lambda*x
         B = 'G' -> generalized eigenvalue problem A*x = lambda*M*x

 N       Integer.  (INPUT)
         Dimension of the eigenproblem.

 K       Integer.  (INPUT)
         Current order of H and the number of columns of V.

 NP      Integer.  (INPUT)
         Number of additional Arnoldi steps to take.

 MODE    Integer.  (INPUT)
         Signifies which form for "OP". If MODE=2 then
         a reduction in the number of B matrix vector multiplies
         is possible since the B-norm of OP*x is equivalent to
         the inv(B)-norm of A*x.

 RESID   Double precision array of length N.  (INPUT/OUTPUT)
         On INPUT:  RESID contains the residual vector r_{k}.
         On OUTPUT: RESID contains the residual vector r_{k+p}.

 RNORM   Double precision scalar.  (INPUT/OUTPUT)
         On INPUT the B-norm of r_{k}.
         On OUTPUT the B-norm of the updated residual r_{k+p}.

 V       Double precision N by K+NP array.  (INPUT/OUTPUT)
         On INPUT:  V contains the Arnoldi vectors in the first K
         columns.
         On OUTPUT: V contains the new NP Arnoldi vectors in the next
         NP columns.  The first K columns are unchanged.

 LDV     Integer.  (INPUT)
         Leading dimension of V exactly as declared in the calling
         program.

 H       Double precision (K+NP) by 2 array.  (INPUT/OUTPUT)
         H is used to store the generated symmetritridiagonal matrix
         with the subdiagonal in the first column starting at H(2,1)
         and the main diagonal in the second column.

 LDH     Integer.  (INPUT)
         Leading dimension of H exactly as declared in the calling
         program.

 IPNTR   Integer array of length 3.  (OUTPUT)
         Pointer to mark the starting locations in the WORK for
         vectors used by the Arnoldi iteration.
         -------------------------------------------------------------
         IPNTR(1): pointer to the current operand vector X.
         IPNTR(2): pointer to the current result vector Y.
         IPNTR(3): pointer to the vector B * X when used in the
                   shift-and-invert mode.  X is the current operand.
         -------------------------------------------------------------

 WORKD   Double precision work array of length 3*N.  (REVERSE COMMUNICATION)
         Distributed array to be used in the basiArnoldi iteration
         for reverse communication.  The calling program should not
         use WORKD as temporary workspace during the iteration !!!!!!
         On INPUT, WORKD(1:N) = B*RESID where RESID is associated
         with the K step Arnoldi factorization. Used to save some
         computation at the first step.
         On OUTPUT, WORKD(1:N) = B*RESID where RESID is associated
         with the K+NP step Arnoldi factorization.


Return value:
 INFO    Integer.  (OUTPUT)
         = 0: Normal exit.
         > 0: Size of an invariant subspace of OP is found that is
              less than K + NP.
"""
function dsaitr!(
  ido::Ref{Int},
  ::Val{BMAT},
  n::Int,
  k::Int,
  np::Int,
  mode::Int,
  resid::AbstractVecOrMat{T},
  rnorm::Ref{TR},
  V::AbstractMatrix{T},
  ldv::Int, # TODO, try and remove
  H::AbstractMatrix{TR},
  ldh::Int, 
  ipntr::AbstractVecOrMat{Int},
  workd::AbstractVecOrMat{T},
  state::AbstractArpackState{TR};
  stats::Union{ArpackStats,Nothing}=nothing,
  debug::Union{ArpackDebug,Nothing}=nothing,
  idonow::Union{ArpackOp,Nothing}=nothing
) where {T, TR, BMAT}

  @attach_aitr_state(state)

  # check lengths
  # RESID   Double precision array of length N.  (INPUT/OUTPUT)
  # V       Double precision N by K+NP array.  (INPUT/OUTPUT)
  # H       Double precision (K+NP) by 2 array.  (INPUT/OUTPUT)
  @jl_arpack_check_length(resid, n)
  @jl_arpack_check_size(V, n, k+np)
  @jl_arpack_check_size(H, k+np, 2)
  @jl_arpack_check_bmat(BMAT)

  #=
  c        | Initialize timing statistics  |
  c        | & message level for debugging |
  =#
  # in Arpack, saved in 'save' variables, but just recomputed here ...
  safmin = floatmin(TR)
  msglvl = @jl_arpack_debug(maitr,0)

  info = 0
  firststep = false # can't get these from reverse comm, so no need to save 
  step2 = false 
  step5 = false 
  laststep = false 

  if ido[] == 0 # first run
    #=
    c        | Initial call to this routine |
    =#
    t0 = @jl_arpack_time()
    step3 = false
    step4 = false
    rstart = false
    orth1 = false
    orth2 = false
    # c        | Pointer to the current step of |
    # c        | the factorization to build     |
    j = k+1
    # c        | Pointers used for reverse communication  |
    # c        | when using WORKD.                        |
    ipj = 1
    irj = ipj + n
    ivj = irj + n
    #
    iter = 0 
    itry = 0 
    rnorm1 = zero(T) 
    wnorm = zero(T) 
    # add a variable for julia for the firststep
    firststep = true 
  else
    @assert(step3 || step4 || rstart || orth1 || orth2, "need a valid configuration")
  end

  #=
  c     | When in reverse communication mode one of:      |
  c     | STEP3, STEP4, ORTH1, ORTH2, RSTART              |
  c     | will be .true.                                  |
  c     | STEP3: return from computing OP*v_{j}.          |
  c     | STEP4: return from computing B-norm of OP*v_{j} |
  c     | ORTH1: return from computing B-norm of r_{j+1}  |
  c     | ORTH2: return from computing B-norm of          |
  c     |        correction to the residual vector.       |
  c     | RSTART: return from OP computations needed by   |
  c     |         dgetv0.                                 |
  =#
  #= in Julia, we implement this with a while loop to allow us
  to easily restart by running "continue" =#
  while true
    if firststep
      firststep = false 

      # this is label 1000
      if msglvl > 2
        @assert debug !== nothing
        println(debug.logfile, "_saitr: generating Arnoldi vector no. ", j)
        println(debug.logfile, "_saitr: B-norm of the current residual ", rnorm[])
      end
      if rnorm[] <= 0
        #=
        c           | Invariant subspace found, generate a new starting |
        c           | vector which is orthogonal to the current Arnoldi |
        c           | basis and continue the iteration.                 |
        =#
        if msglvl > 0
          println(debug.logfile, "_saitr: ****** restart at step ****** ", j)
        end
        #= 
        c           | ITRY is the loop variable that controls the |
        c           | maximum amount of times that a restart is   |
        c           | attempted. NRSTRT is used by stat.h         |
        =#
        @jl_arpack_increment_stat(nrstrt)
        itry = 1
        # label 20 in dsaitr.f 
        rstart = true
        ido[] = 0 # NOTE, dgetv0 will set ido[] != 0
        # need to get to label 30 in dsaitr.f,
        # which will happen after the first step... 
        continue
      else
        # need to get to label 40
        step2 = true 
        continue 
      end
    else # firststep is false ... either we are back from reverse communication or moving around
      if step2 
        @debug "start of step 2, label 40 in dsaitr.f"
        #=
        c        | STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm  |
        c        | Note that p_{j} = B*r_{j-1}. In order to avoid overflow |
        c        | when reciprocating a small RNORM, test against lower    |
        c        | machine bound.                                          |
        =#
        copyto!(@view(V[1:n,j]), @view(resid[1:n]))
        if rnorm[] >= safmin
          temp1 = one(T)/rnorm[]
          _dscal!(temp1, @view(V[1:n,j]))
          _dscal!(temp1, @view(workd[ipj:ipj+n-1]))
        else
          _scale_from_to(rnorm[], one(T), @view(V[1:n,j]))
          _scale_from_to(rnorm[], one(T), @view(workd[ipj:ipj+n-1]))
        end
        step2 = false 
        # setup for step 3
        #=
        c        | STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j} |
        c        | Note that this is not quite yet r_{j}. See STEP 4    |
        =#
        step3 = true
        @jl_arpack_set_stat(nopx, stats.nopx+1)
        t2 = @jl_arpack_time()
        copyto!(@view(workd[ivj:ivj+n-1]), @view(V[1:n,j]))
        ipntr[1] = ivj
        ipntr[2] = irj
        ipntr[3] = ipj
        ido[] = 1
        if idonow === nothing
          # c        | Exit in order to compute OP*v_{j} |
          break # goto 9000, this will save the current state before return
        else
          # we do it right away in Julia...
          if mode == 1
            _i_do_now_opx_1!(idonow, ipntr, workd, n) # handle the operation
          elseif mode == 2
            _i_do_now_opx_mode2_1!(idonow, ipntr, workd, n) # handle the operation
          else # this means we are mode 3, 4, 5
            # so bx is available in ipntr[3]
            _i_do_now_opx_shiftinvert_1!(idonow, ipntr, workd, n) # handle the operation
          end 
        end
        # if we get here, step3 = true, step2 = false, so we will restart in step3
      elseif step3
        @debug "start of step 3, label 50 in dsaitr.f" # label 50 in dsaitr.f
        # c        | Back from reverse communication;  |
        # c        | WORKD(IRJ:IRJ+N-1) := OP*v_{j}.   |
        @jl_update_time(tmvopx, t2)
        step3 = false
        # c        | Put another copy of OP*v_{j} into RESID. |
        copyto!(@view(resid[1:n]), @view(workd[irj:irj+n-1]))

        # setup for step 4
        #=
        c        | STEP 4:  Finish extending the symmetric   |
        c        |          Arnoldi to length j. If MODE = 2 |
        c        |          then B*OP = B*inv(B)*A = A and   |
        c        |          we don't need to compute B*OP.   |
        c        | NOTE: If MODE = 2 WORKD(IVJ:IVJ+N-1) is   |
        c        | assumed to have A*v_{j}.                  |        
        =#
        step4 = true
        if mode == 2
          # do nothing...
          # FORTRAN if (mode .eq. 2) go to 65
          # this will just swing around and head on down 
          # to step 4... because step4 = true
        else
          t2 = @jl_arpack_time()
          if BMAT == :G
            # compute something with B*Op*v ...
            @jl_arpack_set_stat(nbx, stats.nbx+1)
            ipntr[1] = irj
            ipntr[2] = ipj
            ido[] = 2 
            if idonow === nothing
              break # break out of while loop to return ...
            else
              _i_do_now_bx!(idonow, ipntr, workd, n)
            end 
          elseif BMAT == :I
            copyto!(@view(workd[ipj:ipj+n-1]), @view(resid[1:n]))
          end
        end
      elseif step4
        # we aren't quite at label 60 here, because we should only
        # hit label 60 if mode isn't 2... 
        step4 = false 
        if mode == 2
          # this was in the case above, so we don't want to do more work here either.
        else
          @debug "label 60 in dsaitr.f" # label 60 in dsaitr.f
          # c        | Back from reverse communication;  |
          # c        | WORKD(IPJ:IPJ+N-1) := B*OP*v_{j}. |
          if BMAT == :G
            @jl_update_time(tmvbx, t2)
          end
        end

        # c        | The following is needed for STEP 5. |
        # c        | Compute the B-norm of OP*v_{j}.     |

        @debug "label 65 in dsaitr.f" # label 65 in dsaitr.f
        if mode == 2
          #=
          c           | Note that the B-norm of OP*v_{j} |
          c           | is the inv(B)-norm of A*v_{j}.   |
          =#
          wnorm = dot(@view(resid[1:n]), @view(workd[ivj:ivj+n-1]))
          wnorm = sqrt(abs(wnorm))
        elseif BMAT == :G
          wnorm = dot(@view(resid[1:n]), @view(workd[ipj:ipj+n-1]))
          wnorm = sqrt(abs(wnorm))
        elseif BMAT == :I
          #wnorm = _dnrm2(@view(resid[1:n])) # TODO, implement _dnrm2
          #wnorm = _dnrm2_unroll_ext(@view(resid[1:n]))
          wnorm = norm2(TR, @view(resid[1:n]))
        end

        #=
        c        | Compute the j-th residual corresponding |
        c        | to the j step factorization.            |
        c        | Use Classical Gram Schmidt and compute: |
        c        | w_{j} <-  V_{j}^T * B * OP * v_{j}      |
        c        | r_{j} <-  OP*v_{j} - V_{j} * w_{j}      |
        c
        c        | Compute the j Fourier coefficients w_{j} |
        c        | WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.  |
        =#

        # TODO check to make sure the workd[irj:irj+j-1] is correct here...
        if mode != 2
          #_dgemv_blas!('T', n, j, 1.0, V, ldv, @view(workd[ipj:ipj+n-1]), 
          #  0.0, @view(workd[irj:irj+j-1]))
          # A B α + C β ==> mul!(C, A, B, α, β)
          mul!(@view(workd[irj:irj+j-1]), adjoint(@view(V[1:n, 1:j])), @view(workd[ipj:ipj+n-1]))
        else 
          # the difference is ipj (mode != 2) and ivj (mode == 2)
          # _dgemv_blas!('T', n, j, 1.0, V, ldv, @view(workd[ivj:ivj+n-1]), 
          #  0.0, @view(workd[irj:irj+j-1]))
          mul!(@view(workd[irj:irj+j-1]), adjoint(@view(V[1:n, 1:j])), @view(workd[ivj:ivj+n-1]))
        end
      
        # c        | Orthgonalize r_{j} against V_{j}.    |
        # c        | RESID contains OP*v_{j}. See STEP 3. | 
      
        #_dgemv_blas!('N', n, j, -1.0, V, ldv, @view(workd[irj:irj+j-1]), 
        #  1.0, @view(resid[1:n]))
        mul!(@view(resid[1:n]), @view(V[1:n, 1:j]), @view(workd[irj:irj+j-1]), -one(T), one(T))

        # c        | Extend H to have j rows and columns. |
        H[j,2] = real(workd[irj+j-1])
        # Draft code to try and get better precision... 
        #=
        if sizeof(T) < sizeof(TR) 
          # compute this with higher-precision
          if mode != 2
            H[j,2] = ddot(TR, @view(V[:,j]), @view(workd[ipj:ipj+n-1]))
          else
            H[j,2] = ddot(TR, @view(V[:,j]), @view(workd[ivj:ivj+n-1]))
          end
          println("H[j,2] = $(H[j,2])  vs.  $(real(workd[irj+j-1]))  diff = $(real(workd[irj+j-1]) - H[j,2]) ")
        else 
          H[j,2] = real(workd[irj+j-1])
        end 
        =#

        if j == 1 || rstart
          H[j,1] = 0
        else
          H[j,1] = rnorm[]
        end 

        # start the orthgonalization step
        t4 = @jl_arpack_time()
        orth1 = true
        iter = 0 
        # TODO, figure out how to save this sequence...
        t2 = @jl_arpack_time()
        if BMAT == :G
          # compute something with B*Op*v ...
          @jl_arpack_set_stat(nbx, stats.nbx+1)
          copyto!(@view(workd[irj:irj+n-1]), @view(resid[1:n]))
          ipntr[1] = irj
          ipntr[2] = ipj
          ido[] = 2 
          if idonow === nothing
            break # break out of while loop to return ...
          else
            _i_do_now_bx!(idonow, ipntr, workd, n)
          end 
        elseif BMAT == :I
          copyto!(@view(workd[ipj:ipj+n-1]), @view(resid[1:n]))
        end
      elseif orth1
        @debug "orth1==true case, label 70 in dsaitr.f" # label 70 in dsaitr.f
        orth1 = false 
        # c        | Back from reverse communication if ORTH1 = .true. |
        # c        | WORKD(IPJ:IPJ+N-1) := B*r_{j}.                    |
        if BMAT == :G
          @jl_update_time(tmvbx, t2)
        end

        # c        | Compute the B-norm of r_{j}. |
        # TODO, this is a shared function we could out-source...
        # it occurs one other place...
        if BMAT == :G
          rnorm[] = dot(@view(resid[1:n]), @view(workd[ipj:ipj+n-1]))
          rnorm[] = sqrt(abs(rnorm[]))
        elseif BMAT == :I
          #wnorm = _dnrm2(@view(resid[1:n])) # TODO, implement _dnrm2
          #rnorm[] = _dnrm2_unroll_ext(@view(resid[1:n]))
          rnorm[] = norm2(TR, @view(resid[1:n]))
        end

        # end of step4 move to step5 if necessary...

        #=
        c        | STEP 5: Re-orthogonalization / Iterative refinement phase |
        c        | Maximum NITER_ITREF tries.                                |
        c        |                                                           |
        c        |          s      = V_{j}^T * B * r_{j}                     |
        c        |          r_{j}  = r_{j} - V_{j}*s                         |
        c        |          alphaj = alphaj + s_{j}                          |
        c        |                                                           |
        c        | The stopping criteria used for iterative refinement is    |
        c        | discussed in Parlett's book SEP, page 107 and in Gragg &  |
        c        | Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.         |
        c        | Determine if we need to correct the residual. The goal is |
        c        | to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||  |
        =#

        if rnorm[] > 0.717*wnorm
          # GOTO label 100
          laststep = true 
        else 
          step5 = true
        end
      elseif step5 == true
        @debug "label 80 in dsaitr.f"
        step5 = false 
        
        # see LOOP ISSUE below with respect to label 80
        @jl_arpack_increment_stat(nrorth)
        #=
        c        | Enter the Iterative refinement phase. If further  |
        c        | refinement is necessary, loop back here. The loop |
        c        | variable is ITER. Perform a step of Classical     |
        c        | Gram-Schmidt using all the Arnoldi vectors V_{j}  |
        =#
          
        if msglvl > 2
          # no xtemp here... just use two print statements...
          println(debug.logfile, "_saitr: re-orthonalization ; wnorm is ", wnorm)
          println(debug.logfile, "_saitr: re-orthonalization ; rnorm is ", rnorm[])
        end
        # c        | Compute V_{j}^T * B * r_{j}.                       |
        # c        | WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1). |
        #_dgemv_blas!('T', n, j, one(T), V, ldv, @view(workd[ipj:ipj+n-1]), 
        #  zero(T), @view(workd[irj:irj+j-1]))
        mul!(@view(workd[irj:irj+j-1]), adjoint(@view(V[1:n,1:j])), @view(workd[ipj:ipj+n-1]))
        #=
        c        | Compute the correction to the residual:      |
        c        | r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1).  |
        c        | The correction to H is v(:,1:J)*H(1:J,1:J) + |
        c        | v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j, but only   |
        c        | H(j,j) is updated.                           |
        =#
        #_dgemv_blas!('N', n, j, -one(T), V, ldv, 
        #    @view(workd[irj:irj+j-1]), 
        #    one(T), @view(resid[1:n]))
        mul!(@view(resid[1:n]), @view(V[1:n, 1:j]), @view(workd[irj:irj+j-1]), -one(T), one(T))
        if j == 1 || rstart 
          H[j,1] = 0
        end
        H[j,2] = H[j,2] + real(workd[irj+j-1])

        orth2 = true

        # TODO, figure out how to save this sequence... t is reused...
        t2 = @jl_arpack_time()
        if BMAT == :G
          # compute something with B*Op*v ...
          @jl_arpack_increment_stat(nbx)
          copyto!(@view(workd[irj:irj+n-1]), @view(resid[1:n]))
          ipntr[1] = irj
          ipntr[2] = ipj
          ido[] = 2 
          if idonow === nothing
            break # break out of while loop to return ...
          else
            _i_do_now_bx!(idonow, ipntr, workd, n)
          end 
        elseif BMAT == :I
          copyto!(@view(workd[ipj:ipj+n-1]), @view(resid[1:n]))
        end
        # c           | Exit in order to compute B*r_{j}. |
        # c           | r_{j} is the corrected residual.  |
      elseif orth2
        @debug "orth2 == true, label 90 in dsaitr.f"

        # This line isn't in dsaitr.f, but given our looping logic, 
        # we need it here because we need to get down to laststep
        # below if rnorm1 > rnorm or if we aren't doing more
        # iterative refinement
        orth2 = false 

        # c        | Back from reverse communication if ORTH2 = .true. |
        if BMAT == :G
          @jl_update_time(tmvbx, t2)
        end

        # c        | Compute the B-norm of the corrected residual r_{j}. |
        if BMAT == :G
          rnorm1 = dot(@view(resid[1:n]), @view(workd[ipj:ipj+n-1]))
          rnorm1 = sqrt(abs(rnorm1))
        elseif BMAT == :I
          # rnorm1 = _dnrm2_unroll_ext(@view(resid[1:n]))
          rnorm1 = norm2(TR, @view(resid[1:n]))
        end

        if msglvl > 0 && iter > 0
          println(debug.logfile, "_saitr: iterative refinement for Arnoldi residual ", j)
          if msglvl > 2
            # no xtemp here... just use two print statements...
            println(debug.logfile, "_saitr: iterative refinement ; rnorm is ", rnorm[])
            println(debug.logfile, "_saitr: iterative refinement ; rnorm1 is ", rnorm1)
          end
        end

        # c        | Determine if we need to perform another |
        # c        | step of re-orthogonalization.           |
        if rnorm1 > 0.717*rnorm[]
          rnorm[] = rnorm1 
          laststep = true 
        else
          # c           | Another step of iterative refinement step |
          # c           | is required. NITREF is used by stat.h     |
          @jl_arpack_increment_stat(nitref)
          rnorm[] = rnorm1 
          iter += 1
          if iter <= 1
            step5 = true # go back up to label 80...
            # LOOP ISSUE
            # but _decrement_ nrorth, this is because we always 
            # increase north in "our" label 80, but Arpack
            # didn't if we "go to" label 80...
            @jl_arpack_set_stat(nrorth, stats.nrorth - 1)
          else
            # c           | Otherwise RESID is numerically in the span of V |
            fill!(@view(resid[1:n]), 0)
            rnorm[] = 0 
            # need to get to label 100 now, so we are going to set laststep = true
            laststep = true
          end 
        end
      elseif laststep # this is the end of the iteration
        #=
        c        | Branch here directly if iterative refinement |
        c        | wasn't necessary or after at most NITER_REF  |
        c        | steps of iterative refinement.               |
        =# 
        @debug "laststep == true, label 100 in dsaitr.f"
        laststep = false 

        # label 100 in the Fortran code! 
        rstart = false
        orth2 = false 
        @jl_update_time(titref, t4)
        #=
        c        | Make sure the last off-diagonal element is non negative  |
        c        | If not perform a similarity transformation on H(1:j,1:j) |
        c        | and scale v(:,j) by -1.                                  |
        =#
        if H[j,1] < 0
          H[j,1] = - H[j,1]
          if j < k+np 
            _dscal!(-one(T), @view(V[1:n,j+1]))
          else
            _dscal!(-one(T), @view(resid[1:n]))
          end
        end 

        #=
        c        | STEP 6: Update  j = j+1;  Continue |
        =#
        j = j+1

        if j > k+np
          # we are done... prep for exit from the routine...
          @jl_update_time(taitr, t0)
          ido[] = 99
          if msglvl > 1
            _arpack_vout(debug, "_saitr: main diagonal of matrix H of step K+NP", 
              @view(H[1:k+np, 2]))
            if k+np > 1
              _arpack_vout(debug, "_saitr: sub diagonal of matrix H of step K+NP", 
                @view(H[2:k+np, 1]))
            end
          end
          break 
        else
          # c        | Loop back to extend the factorization by another step. |
          firststep = true
        end
      elseif rstart # rstart often stays on for a while, so it should come last...
        # label 30 in dsaitr.f
        @debug "calling dgetv0"
        ierr = dgetv0!(ido, Val(BMAT), itry, false, n, j, V, 
                      ldv, resid, rnorm, ipntr, workd, state; 
                      debug, stats, idonow)
        if ido[] != 99
          break # we want to exit the while loop and return
        end
        if ierr < 0 
          itry += 1

          if itry <= 3
            # we just inline the code for label 20 here 20 again here...
            rstart = true
            ido[] = 0 
            continue # this will move us back to label 30
          end
          #=
          c              | Give up after several restart attempts.        |
          c              | Set INFO to the size of the invariant subspace |
          c              | which spans OP and exit.                       |
          =#
          info = j-1 
          @jl_update_time(taitr, t2)
          ido[] = 99 
          break #
        end
        step2 = true
        # no real need for continue here... 
      else
        error("impossible situtation")
      end
    end
  end

  state.aitr = AitrState{TR}(@aitr_state_vars)
  return info
end
