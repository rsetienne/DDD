! Helper function: 
! fill vec with N elements from parms, starting at position ii
!==========================================================================

      SUBROUTINE dd_fill1d (vec, DIMP, parms, II)
      IMPLICIT NONE
      INTEGER DIMP, II, I
      DOUBLE PRECISION vec(DIMP), parms(*)
      II = II
        DO I = 1, DIMP
          II = II + 1
          vec(I) = parms(II)
        ENDDO
        
      END SUBROUTINE dd_fill1d

!==========================================================================
! module with declarations
!==========================================================================

      MODULE dd_dimmod

      ! length of the vector -  decided in R-code
      INTEGER  :: N
      INTEGER  :: kk
      
      ! 1 parameter vectors with unknown length
      DOUBLE PRECISION, ALLOCATABLE  :: P(:)  
      
      ! Boolean: will become TRUE if the parameters have a value
      LOGICAL :: initialised = .FALSE.

      END MODULE dd_dimmod

!==========================================================================
!==========================================================================
! Initialisation: name of this function as passed by "initfunc" argument
! Sets the fixed parameter vector, and allocates memory
!==========================================================================
!==========================================================================

      SUBROUTINE dd_initmod (steadyparms)
      USE dd_dimmod 

      IMPLICIT NONE
      EXTERNAL steadyparms

      INTEGER, PARAMETER :: nparsmall = 2  ! constant-length parameters
      
      DOUBLE PRECISION parms(nparsmall)
      COMMON /XCBPar/parms                 ! common block 

! Set the fixed parameters obtained from R
      CALL steadyparms(nparsmall, parms)

! first parameter has the length of the vector       
      N = INT(parms(1) + 1e-6)
      kk = INT(parms(2) + 1e-6)

! Allocate variable size arrays (state variables, derivatives and parameters)

      IF (ALLOCATED(P)) DEALLOCATE(P)  
      ALLOCATE(P(3 * (N + 2 + 2 * kk)))

      initialised = .FALSE.
       
      END SUBROUTINE dd_initmod
      
!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================
       
      SUBROUTINE dd_runmod (neq, t, Conc, dConc, yout, ip)
      USE dd_dimmod
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii
      DOUBLE PRECISION  :: t, Conc(N), dConc(N), yout(*)
      DOUBLE PRECISION  :: V(N + 2)
      DOUBLE PRECISION  :: lavec(N + 2 + 2 * kk),muvec(N + 2 + 2 * kk)
      DOUBLE PRECISION  :: nn(N + 2 + 2 * kk)
      DOUBLE PRECISION  :: FF1, FF2, FF3

! parameters - named here
      DOUBLE PRECISION rn(2)
      COMMON /XCBPar/rn

! local variables
      CHARACTER(len=100) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough") 

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL dd_fill1d(P, 3 * (N + 2 + 2 * kk), yout, ii)   ! ii is updated in fill1d
        Initialised = .TRUE.          ! to prevent from initialising more than once
      ENDIF

! dynamics

 !  dx = lavec[(2:(lx+1))+kk-1] * nn[(2:(lx+1))+2*kk-1] * xx[(2:(lx+1))-1] + muvec[(2:(lx+1))+kk+1] 
    !  * nn[(2:(lx+1))+1] * xx[(2:(lx+1))+1] - (lavec[(2:(lx+1))+kk] + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]

      V(1) = 0
      DO I = 2, N + 1 
        V(I) = Conc(I - 1)
      ENDDO
      V(N + 2) = 0
      DO I = 1, N + 2 + 2 * kk
       lavec(I) = P(I)
       muvec(I) = P(I + N + 2 + 2 * kk)
       nn(I)    = P(I + 2 * (N + 2 + 2 * kk))
      ENDDO

      DO I = 2, N + 1 
        FF1 = lavec(I + kk - 1) * nn(I + 2 * kk - 1) * V(I - 1)
        FF2 = muvec(I + kk + 1) * nn(I + 1) * V(I + 1)
        FF3 = (lavec(I + kk) + muvec(I + kk)) * nn(I + kk) * V(I)
        dConc(I - 1) = FF1 + FF2 - FF3
      ENDDO
  
      END SUBROUTINE dd_runmod
      
!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
!==========================================================================
!==========================================================================
       
      SUBROUTINE dd_runmodbw (neq, t, Conc, dConc, yout, ip)
      USE dd_dimmod
      IMPLICIT NONE

!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii
      DOUBLE PRECISION  :: t, Conc(N), dConc(N), yout(*)
      DOUBLE PRECISION  :: V(N + 2)
      DOUBLE PRECISION  :: lavec(N + 2 + 2 * kk),muvec(N + 2 + 2 * kk)
      DOUBLE PRECISION  :: nn(N + 2 + 2 * kk)
      DOUBLE PRECISION  :: FF1, FF2, FF3

! parameters - named here
      DOUBLE PRECISION rn(2)
      COMMON /XCBPar/rn

! local variables
      CHARACTER(len=100) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough") 

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL dd_fill1d(P, 3 * (N + 2 + 2 * kk), yout, ii)   ! ii is updated in fill1d
        Initialised = .TRUE.          ! to prevent from initialising more than once
      ENDIF

! dynamics

! dx = lavec[(2:(lx+1))+kk] * nn[(2:(lx+1))+2*kk] * xx[(2:(lx+1))+1] + muvec[(2:(lx+1))+kk] * nn[(2:(lx+1))] * xx[(2:(lx+1))-1] - (c(lavec[(2:(lx))+kk],0) + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]
! dG = x[1 + (kk == 0)]


      V(1) = 0
      DO I = 2, N 
        V(I) = Conc(I - 1)
      ENDDO
      V(N + 1) = 0
      V(N + 2) = 0
      DO I = 1, N + 2 + 2 * kk
       lavec(I) = P(I)
       muvec(I) = P(I + N + 2 + 2 * kk)
       nn(I)    = P(I + 2 * (N + 2 + 2 * kk))
      ENDDO

      DO I = 2, N - 1
        FF1 = lavec(I + kk) * nn(I + 2 * kk) * V(I + 1)
        FF2 = muvec(I + kk) * nn(I) * V(I - 1)
        FF3 = (lavec(I + kk) + muvec(I + kk)) * nn(I + kk) * V(I)
        dConc(I - 1) = FF1 + FF2 - FF3
      ENDDO
      I = N
      FF1 = lavec(I + kk) * nn(I + 2 * kk) * V(I + 1)
      FF2 = muvec(I + kk) * nn(I) * V(I - 1)
      FF3 = (0 + muvec(I + kk)) * nn(I + kk) * V(I)
      dConc(I - 1) = FF1 + FF2 - FF3
      I = N + 1
      IF (kk .eq. 0) THEN
         dConc(I - 1) = Conc(2)
      ELSE
         dConc(I - 1) = Conc(1)
      ENDIF      

      END SUBROUTINE dd_runmodbw
      
!==========================================================================
!==========================================================================
! Dynamic routine: name of this function as passed by "func" argument
! variable parameter values are passed via yout
! Reinterpretation: N is the length of probs and kk is the number of sigmas
!==========================================================================
!==========================================================================
       
      SUBROUTINE dd_runmodtd (neq, t, Conc, dConc, yout, ip)
      USE dd_dimmod
      IMPLICIT NONE
!......................... declaration section.............................
      INTEGER           :: neq, ip(*), i, ii, M
      DOUBLE PRECISION  :: t, Conc(N), dConc(N), yout(*)
      DOUBLE PRECISION  :: V(N + 2)
      DOUBLE PRECISION  :: lavec((N - kk) + 2),muvec((N - kk) + 2)
      DOUBLE PRECISION  :: nn((N - kk) + 2)
      DOUBLE PRECISION  :: FF1, FF2, FF3, En
      DOUBLE PRECISION  :: c, t1, y
      REAL(16)          :: Envec(N - kk)

! parameters - named here
      DOUBLE PRECISION rn(2)
      COMMON /XCBPar/rn

! local variables
      CHARACTER(len=100) msg

!............................ statements ..................................

      IF (.NOT. Initialised) THEN
        ! check memory allocated to output variables
        IF (ip(1) < 1) CALL rexit("nout not large enough") 

        ! save parameter values in yout
        ii = ip(1)   ! Start of parameter values
        CALL dd_fill1d(P, 3 * ((N - kk) + 2), yout, ii)   ! ii is updated in fill1d
        Initialised = .TRUE.          ! to prevent from initialising more than once
      ENDIF

! dynamics

 !  dp = lavec[(2:(lp+1))-1] * nn[(2:(lp + 1))-1] * p[(2:(lp + 1))-1] + muvec[(2:(lp+1))+1] * nn[(2:(lp+1))+1] * p[(2:(lp+1))+1] - (lavec[(2:(lp+1))] + muvec[(2:(lp+1))]) * nn[(2:(lp+1))] * p[(2:(lp+1))]
 !  This is the same as setting kk = 0 in
 !  dx = lavec[(2:(lx+1))+kk-1] * nn[(2:(lx+1))+2*kk-1] * xx[(2:(lx+1))-1] + muvec[(2:(lx+1))+kk+1] 
    !  * nn[(2:(lx+1))+1] * xx[(2:(lx+1))+1] - (lavec[(2:(lx+1))+kk] + muvec[(2:(lx+1))+kk]) * nn[(2:(lx+1))+kk] * xx[2:(lx+1)]
 !  mutd = rep(mu,lrs)
 !  En = sum((0:(lx - 1)) * x[1:lx] )
 !  dsigdiv = mutd / En

      M = N - kk
      V(1) = 0
      DO I = 2, M + 1 
        V(I) = Conc(I - 1)
      ENDDO
      V(M + 2) = 0
      DO I = 1, M + 2
       lavec(I) = P(I)
       muvec(I) = P(I + M + 2)
       nn(I)    = P(I + 2 * (M + 2))
      ENDDO
      
      Envec(1) = 0
      DO I = 2, M
       Envec(I)       = (I - 1) * Conc(I)
      ENDDO
      
      !Naive Sum
      !En = Envec(1)
      !DO I = 2,M
      !  En = En + Envec(I)
      !ENDDO
      !
      !Kahan Sum
      !En = Envec(1)
      !c = 0.0
      !DO I = 2, M
      !  y = Envec(I) - c
      !  t1 = En + y
      !  c = (t1 - En) - y
      !  En = t1
      !ENDDO
      !
      !Improved Kahan Sum
      !En = Envec(1)
      !c = 0.0
      !DO I = 2, M
      !   t1 = En + Envec(I)
      !   IF (ABS(En) >= ABS(Envec(I))) THEN
      !      c = c + (En - t1) + Envec(I)
      !   ELSE
      !      c = c + (Envec(I) - t1) + En
      !   ENDIF
      !   En = t1
      !ENDDO
      !En = En + c
      
      En = SUM(Envec)

      DO I = 2, M + 1
        FF1 = lavec(I + 0 - 1) * nn(I + 2 * 0 - 1) * V(I - 1)
        FF2 = muvec(I + 0 + 1) * nn(I + 1) * V(I + 1)
        FF3 = (lavec(I + 0) + muvec(I + 0)) * nn(I + 0) * V(I)
        dConc(I - 1) = FF1 + FF2 - FF3
      ENDDO
      
      DO I = 1,kk
        dConc(M + I) = muvec(I)/En
      ENDDO
  
      END SUBROUTINE dd_runmodtd
