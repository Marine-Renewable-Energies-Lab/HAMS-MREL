!  ------------------------------------------------------------------------------------------------------
!                                                               
!    Program HAMS-MREL for the diffraction and radiation of waves 
!    by 3D structures.
! 
!  License:
! 
!    This routine is part of HAMS and HAMS-MREL.
!
!    HAMS-MREL is a free software framework: you can redistribute it and/or modify it 
!    under the terms of the Apache License, Version 2.0 (the "License"); you may 
!    not use these subroutines except in compliance with the License. The software
!    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND.
!
!    You should have received a copy of the Apache License, Version 2.0, along with 
!    HAMS-MREL. If not, see <http://www.apache.org/licenses/LICENSE-2.0>.
!
!  Code Original Author:
!
!    Yingyi Liu
!
!  ------------------------------------------------------------------------------------------------------
MODULE ImplementSubs

   USE HAMS_mod
   USE Const_mod 
   USE WaveDyn_mod

   IMPLICIT NONE

   ! ..... Public Subroutines ...................................................................................................
   
   PUBLIC :: CalWaveProperts
   PUBLIC :: SetWaveFromOmega
   PUBLIC :: InputOmega
   PUBLIC :: WarnLargeKH
   PUBLIC :: CheckFiniteResults

   ! Number of frequencies with a NaN or infinite result, and the first such frequency in rad/s
   INTEGER, PUBLIC :: NONFINITE_FREQS = 0
   REAL*8, PUBLIC :: FIRST_NONFINITE_FREQ = 0.D0

CONTAINS

! ------------------------------------------------------------------- 
!    Calculate the incident wave properties
! ------------------------------------------------------------------- 
      SUBROUTINE CalWaveProperts(KK)
      IMPLICIT  NONE
!        
      INTEGER,INTENT(IN):: KK

        IF (SYBO.EQ.1) THEN
          IF (KK.EQ.1) THEN
            V=0.D0
            WK=0.D0
            W1=0.D0
            TP=-1.D0
            WL=-1.D0
            IF (INFT.LE.3) THEN
             INFR=0.D0
            ELSE
             INFR=-1.D0
            ENDIF
            GOTO 100
            ELSEIF (KK.EQ.2) THEN
            V=-1.D0
            WK=-1.D0
            W1=-1.D0
            TP=0.D0
            WL=0.D0
            IF (INFT.LE.3) THEN
             INFR=-1.D0
            ELSE
             INFR=0.D0
            ENDIF
            GOTO 100
          ENDIF
        ENDIF
        
        INFR=WVNB(KK)
        IF (INFT .EQ. 1)  THEN
            V=WVNB(KK)
          IF (H .LE. 0.0D0) THEN
            WK=V
            W1=SQRT(G*V)
          ELSE
            W1=SQRT(G*V)
            CALL DISPERSION(WVN,NK,W1,H)
            WK=WVN(1)
          END IF
            WL=2.D0*PI/WK
            TP=2.D0*PI/W1
        ELSEIF (INFT .EQ. 2)  THEN
            WK=WVNB(KK)
          IF (H .LE. 0.0D0) THEN
            V=WK
            W1=SQRT(G*V)
          ELSE
            V=WK*TANH(WK*H)
            W1=SQRT(G*V)
            CALL DISPERSION(WVN,NK,W1,H)
          END IF
            WL=2.D0*PI/WK
            TP=2.D0*PI/W1
        ELSEIF (INFT .EQ. 3)  THEN
            CALL SetWaveFromOmega(WVNB(KK))
        ELSEIF (INFT .EQ. 4)  THEN
            CALL SetWaveFromOmega(2.D0*PI/WVNB(KK))
            TP=WVNB(KK)
        ELSEIF (INFT .EQ. 5)  THEN
            WL=WVNB(KK)
            WK=2.D0*PI/WL
          IF(H .LE. 0.0D0) THEN
            V=WK
            W1=SQRT(G*V)
          ELSE
            V=WK*TANH(WK*H)
            W1=SQRT(G*V)
            CALL DISPERSION(WVN,NK,W1,H)
          END IF
        END IF

100     CONTINUE
        WVNB(KK)=WK
        WVFQ(KK)=W1
        
        IF (OUFT .EQ. 1)  THEN
         OUFR=V
        ELSEIF (OUFT .EQ. 2)  THEN
         OUFR=WK
        ELSEIF (OUFT .EQ. 3)  THEN
         OUFR=W1
        ELSEIF (OUFT .EQ. 4)  THEN
         OUFR=TP
        ELSEIF (OUFT .EQ. 5)  THEN
         OUFR=WL
        END IF

      END SUBROUTINE CalWaveProperts

! -------------------------------------------------------------------
!    Wave properties of a finite frequency W (rad/s)
! -------------------------------------------------------------------
!    Sets W1, V, WK (and the roots WVN of the finite-depth Green function), WL and TP, as CalWaveProperts does for
!    INFT = 3, and nothing else: not INFR, OUFR or the per-frequency arrays WVNB and WVFQ. The radiation-only solves of the
!    QTF difference frequencies use it after the frequency loop.
      SUBROUTINE SetWaveFromOmega(W)
      IMPLICIT  NONE
      REAL*8,INTENT(IN):: W

            W1=W
          IF(H .LE. 0.0D0) THEN
            V=W1*W1/G
            WK=V
          ELSE
            V=W1*W1/G
            CALL DISPERSION(WVN,NK,W1,H)
            WK=WVN(1)
          END IF
            WL=2.D0*PI/WK
            TP=2.D0*PI/W1

      END SUBROUTINE SetWaveFromOmega

! -------------------------------------------------------------------
!    Frequency in rad/s of an input frequency value X of type INFT
! -------------------------------------------------------------------
!    The W1 that CalWaveProperts computes from WVNB(KK) = X before the frequency loop overwrites WVNB, without side
!    effects. Not for the limit frequencies of SYBO = 1.
      REAL*8 FUNCTION InputOmega(X)
      IMPLICIT  NONE
      REAL*8,INTENT(IN):: X
      REAL*8:: VX,KX

        IF (INFT .EQ. 1) THEN
          InputOmega=SQRT(G*X)
        ELSEIF (INFT .EQ. 2 .OR. INFT .EQ. 5) THEN
          KX=X
          IF (INFT .EQ. 5) KX=2.D0*PI/X
          IF (H .LE. 0.0D0) THEN
            VX=KX
          ELSE
            VX=KX*TANH(KX*H)
          END IF
          InputOmega=SQRT(G*VX)
        ELSEIF (INFT .EQ. 4) THEN
          InputOmega=2.D0*PI/X
        ELSE
          InputOmega=X
        END IF

      END FUNCTION InputOmega

! -------------------------------------------------------------------
!    Warn when kh is large in finite depth
! -------------------------------------------------------------------
!    The series of the finite-depth Green function lose accuracy at large kh on the small panels a mesh
!    needs there: measured on a cylinder at depth 300, against the deep-water solver, A11 is 2.5% off
!    at kh = 300 and 20% off at kh = 600, and B33 more. At such kh the waves do not reach the seabed,
!    so deep water gives the same result unless the body is close to the seabed.
      SUBROUTINE WarnLargeKH()
      IMPLICIT  NONE
      REAL*8,PARAMETER:: KH_WARN=100.D0

        IF (H.GT.0.D0 .AND. WK.GT.0.D0 .AND. WK*H.GT.KH_WARN) THEN
          WRITE(6,1000) WK*H
        ENDIF
1000    FORMAT(10x,'Warning: kh =',F9.1,' is above 100. The finite depth Green function is', &
               ' inaccurate at large kh;',/,10x,'unless the body is near the seabed,', &
               ' set Waterdepth < 0 for deep water, which gives the same result.')

      END SUBROUTINE WarnLargeKH

! -------------------------------------------------------------------
!    Check the results of one frequency for NaN or Inf
! -------------------------------------------------------------------
      SUBROUTINE CheckFiniteResults(AMS,BDP,EXF)
      USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
      IMPLICIT  NONE
      REAL*8,INTENT(IN):: AMS(:,:),BDP(:,:)
      COMPLEX*16,INTENT(IN):: EXF(:,:)

        IF (ALL(IEEE_IS_FINITE(AMS)) .AND. ALL(IEEE_IS_FINITE(BDP)) .AND. &
            ALL(IEEE_IS_FINITE(DBLE(EXF))) .AND. ALL(IEEE_IS_FINITE(DIMAG(EXF)))) RETURN

        IF (NONFINITE_FREQS.EQ.0) FIRST_NONFINITE_FREQ=W1
        NONFINITE_FREQS=NONFINITE_FREQS+1
        WRITE(6,1000)
1000    FORMAT(10x,'Error: the results at this frequency contain NaN or Inf.')

      END SUBROUTINE CheckFiniteResults
!-------------------------------------------------------------------------------
END MODULE ImplementSubs
!*******************************************************************************