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
MODULE AssbMatx_irr

   USE HAMS_mod
   USE Body_mod
   USE Const_mod
   USE WaveDyn_mod
   USE PanelMesh_mod
   USE Inerfs_mod
   
   USE BodyIntgr_irr
   USE PatcVelct
   USE IterativeSolvers, ONLY: ZITER_SOLVE_MULTI
   IMPLICIT NONE
   
      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: ASSB_LEFT_IRR
   PUBLIC :: ASSB_RBC_IRR
   
   PUBLIC :: ASSB_DBC_IRR
   PUBLIC :: RADIATION_SOLVER_IRR
   PUBLIC :: DIFFRACTION_SOLVER_IRR
   PUBLIC :: NORMAL_MATRIX
   PUBLIC :: NORMAL_RHS
   
CONTAINS
! ------------------------------------------------------------------- 
!    Calculate the element contribution, assembly the left-hand
!    side coefficient matrix [A]
!    This subroutine is for suppression of irregular frequencies
! ------------------------------------------------------------------- 

      SUBROUTINE ASSB_LEFT_IRR(AMAT,CMAT,IPIV,NELEM,TNELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM,TNELEM,NSYS
      INTEGER,INTENT(OUT):: IPIV(NELEM,NSYS)
      COMPLEX*16,INTENT(OUT):: AMAT(TNELEM,TNELEM,NSYS),CMAT(NELEM,NELEM,NSYS)
      
      INTEGER IEL,JEL,KEL,IS,IP,IRR,FLAG,INFO
      REAL*8 DIST
      COMPLEX*16 TINDP(4),SUM

      AMAT=CMPLX(0.0D0,0.0D0)

      IRR=1

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,JEL,IP,IS,FLAG,DIST,TINDP) SCHEDULE(DYNAMIC,16)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.

      DO  1000 IEL=1,  NELEM
        
       DO 100  IP=1,  NSYS

        AMAT(IEL,IEL,IP)=2.D0*PI
        
100    CONTINUE
        
       DO 200 JEL=1,  NELEM

        DIST=SQRT((XYZ_P(IEL,1)-XYZ_P(JEL,1))**2+(XYZ_P(IEL,2)-XYZ_P(JEL,2))**2+(XYZ_P(IEL,3)-XYZ_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF

        TINDP=CMPLX(0.0D0, 0.0D0)

        DO  200   IS=1,  NSYS
        
         CALL BODINT_LEFT_IRR(IS,IEL,JEL,TINDP,IRR,FLAG)

         DO IP=1, NSYS
          ! Both ISX/ISY branches were identical — collapsed.
          AMAT(IEL,JEL,IP)=AMAT(IEL,JEL,IP)+RXY(IS,IP)*TINDP(IS)
         ENDDO
        
200    CONTINUE
     
1000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
       
      IRR=3

!$OMP PARALLEL NUM_THREADS(NTHREAD)  
!$OMP DO PRIVATE(IEL,JEL,IP,IS,FLAG,DIST,TINDP) SCHEDULE(DYNAMIC,16)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
      
      DO  3000 IEL=NELEM+1,  NELEM+iNELEM
       
       DO 500  IP=1,  NSYS

        AMAT(IEL,IEL,IP)=0.D0*PI
        
500    CONTINUE

       DO 600 JEL=1,  NELEM
            
        DIST=SQRT((iXYZ_P(IEL-NELEM,1)-XYZ_P(JEL,1))**2+(iXYZ_P(IEL-NELEM,2)-XYZ_P(JEL,2))**2+(iXYZ_P(IEL-NELEM,3)-XYZ_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
        
        TINDP=CMPLX(0.0D0, 0.0D0)
        
        DO  600   IS=1,  NSYS            
        
         CALL BODINT_LEFT_IRR(IS,IEL-NELEM,JEL,TINDP,IRR,FLAG)

         DO IP=1, NSYS
          ! Both ISX/ISY branches were identical — collapsed.
          AMAT(IEL,JEL,IP)=AMAT(IEL,JEL,IP)+RXY(IS,IP)*TINDP(IS)
         ENDDO
        
600    CONTINUE

3000  CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
!
!--------------------------------------------------------------
!    Matrix set-up for the least-square problem

      DO IP=1, NSYS
       CALL NORMAL_MATRIX(AMAT(:,:,IP),TNELEM,NELEM,CMAT(:,:,IP))
      ENDDO

! NSYS is 1 or 2 — let MKL's internal threading do the work in ZGETRF.
! Only factorize when direct-LU will be used; ISOLV>1 uses the unfactored normal-equations matrix.
       IF (ISOLV .EQ. 1) THEN
          DO IP=1, NSYS
             CALL ZGETRF( NELEM, NELEM, CMAT(:,:,IP), NELEM, IPIV(:,IP), INFO )
          ENDDO
       ENDIF

      RETURN
      END SUBROUTINE ASSB_LEFT_IRR
      
! ------------------------------------------------------------------- 
!    Calculate the element contribution, assembly the right-hand
!    side vector[B]
!    This subroutine is for suppression of irregular frequencies
! ------------------------------------------------------------------- 

      SUBROUTINE ASSB_RBC_IRR(BRMAT,DRMAT,AMAT,NELEM,TNELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM,TNELEM,NSYS
      COMPLEX*16,INTENT(IN):: AMAT(TNELEM,TNELEM,NSYS)
      COMPLEX*16,INTENT(OUT):: BRMAT(TNELEM,6,NSYS),DRMAT(NELEM,6,NSYS)
      
      INTEGER IEL,JEL,KEL,IS,IP,IRR,MD,FLAG
      REAL*8 DIST
      COMPLEX*16 TINRD(4,6,4),BTMP(6,4)

      BRMAT=CMPLX(0.0D0,0.0D0)
!
      IRR=1
  
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,JEL,MD,IP,IS,FLAG,DIST,TINRD,BTMP) SCHEDULE(DYNAMIC,16)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
      
      DO  1000 IEL=1,  NELEM
            
       BTMP=CMPLX(0.0D0,0.0D0)

        DO 200 JEL=1,  NELEM

         DIST=SQRT((XYZ_P(IEL,1)-XYZ_P(JEL,1))**2+(XYZ_P(IEL,2)-XYZ_P(JEL,2))**2+(XYZ_P(IEL,3)-XYZ_P(JEL,3))**2)
         IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
          FLAG=1
         ELSE
          FLAG=0
         ENDIF

         TINRD=CMPLX(0.0D0, 0.0D0)
        
         DO  200   IS=1,  NSYS
        
          CALL RBC_IRR(IS,IEL,JEL,TINRD,IRR,FLAG)
        
         DO MD=1,  6
         DO IP=1, NSYS
          BTMP(MD,IP)=BTMP(MD,IP)+TINRD(IS,MD,IP)
         ENDDO
         ENDDO
        
200     CONTINUE
       
        DO  300  MD=1,  6
        DO  300  IP=1, NSYS

        DO  300  IS=1, NSYS
          ! Both ISX/ISY branches were identical — collapsed.
          BRMAT(IEL,MD,IP)=BRMAT(IEL,MD,IP)+RXY(IP,IS)*BTMP(MD,IS)
300     CONTINUE
     
1000    CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
       
      IRR=3

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,JEL,MD,IP,IS,FLAG,DIST,TINRD,BTMP) SCHEDULE(DYNAMIC,16)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
      
      DO  3000 IEL=NELEM+1,  NELEM+iNELEM
            
       BTMP=CMPLX(0.0D0,0.0D0)

        DO 600 JEL=1,  NELEM

           DIST=SQRT((iXYZ_P(IEL-NELEM,1)-XYZ_P(JEL,1))**2 +&
                (iXYZ_P(IEL-NELEM,2)-XYZ_P(JEL,2))**2 +&
                (iXYZ_P(IEL-NELEM,3)-XYZ_P(JEL,3))**2)
         IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
          FLAG=1
         ELSE
          FLAG=0
         ENDIF

         TINRD=CMPLX(0.0D0, 0.0D0)
        
         DO  600   IS=1,  NSYS
        
          CALL RBC_IRR(IS,IEL-NELEM,JEL,TINRD,IRR,FLAG)
        
         DO MD=1,  6
         DO IP=1, NSYS
          BTMP(MD,IP)=BTMP(MD,IP)+TINRD(IS,MD,IP)
         ENDDO
         ENDDO
        
600     CONTINUE

       DO  700  MD=1,  6
       DO  700  IP=1, NSYS

       DO  700  IS=1, NSYS
          ! Both ISX/ISY branches were identical — collapsed.
          BRMAT(IEL,MD,IP)=BRMAT(IEL,MD,IP)+RXY(IP,IS)*BTMP(MD,IS)
700    CONTINUE
     
3000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
!
!--------------------------------------------------------------
!    Matrix set-up for the least-square problem
!
      !WRITE(6, *) ' Right-hand side matrix setting-up for the least-square problem...'

      DO IP=1, NSYS
       CALL NORMAL_RHS(AMAT(:,:,IP),TNELEM,NELEM,BRMAT(:,:,IP),6,DRMAT(:,:,IP))
      ENDDO

      RETURN
      END SUBROUTINE ASSB_RBC_IRR

! ------------------------------------------------------------------- 
!    Calculate the element contribution, assembly the right-hand
!    side vector[B]
!    This subroutine is for suppression of irregular frequencies
! ------------------------------------------------------------------- 

      SUBROUTINE ASSB_DBC_IRR(BDMAT,DDMAT,AMAT,NELEM,TNELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM,TNELEM,NSYS
      COMPLEX*16,INTENT(IN):: AMAT(TNELEM,TNELEM,NSYS)
      COMPLEX*16,INTENT(OUT):: BDMAT(TNELEM,NSYS),DDMAT(NELEM,NSYS)
      
      INTEGER IEL,JEL,KEL,IS,IP,IRR,MD,FLAG
      REAL*8  XP,YP,ZP,DIST
      COMPLEX*16 TINRD(4,4),BTMP(4)

      MD=7
      BDMAT=CMPLX(0.0D0,0.0D0)
!
      IRR=1

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(XP,YP,ZP,IEL,JEL,IP,IS,FLAG,DIST,TINRD,BTMP) SCHEDULE(DYNAMIC,16)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
      DO  1000 IEL=1,  NELEM

       BTMP=CMPLX(0.0D0,0.0D0)
       
       IF (ISOL.NE.1.AND.ISOL.NE.2) THEN
        PRINT*,"  Error: The input for ISOL should be either 1 or 2."
        STOP
       ENDIF

       ! Both ISOL=1 (total) and ISOL=2 (scattered) solve the same validated
       ! direct BIE for the scattered potential; the total/scattered choice is
       ! applied later by adding back the incident potential F0 (ISOL=1 only).
       DO 200 JEL=1,  NELEM

          DIST=SQRT((XYZ_P(IEL,1)-XYZ_P(JEL,1))**2+(XYZ_P(IEL,2)-XYZ_P(JEL,2))**2+(XYZ_P(IEL,3)-XYZ_P(JEL,3))**2)
          IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
           FLAG=1
          ELSE
           FLAG=0
          ENDIF

          TINRD=CMPLX(0.0D0, 0.0D0)

          DO  200   IS=1,  NSYS

            CALL DBC_IRR(IS,IEL,JEL,TINRD,IRR,FLAG)

          DO IP=1, NSYS
            BTMP(IP)=BTMP(IP)+TINRD(IS,IP)
          ENDDO

200     CONTINUE
 
        DO  300  IP=1, NSYS

        DO  300  IS=1, NSYS
          ! Both ISX/ISY branches were identical — collapsed.
          BDMAT(IEL,IP)=BDMAT(IEL,IP)+RXY(IP,IS)*BTMP(IS)
300     CONTINUE
     
1000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
       
      IRR=3

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(XP,YP,ZP,IEL,JEL,IP,IS,FLAG,DIST,TINRD,BTMP) SCHEDULE(DYNAMIC,16)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
      
      DO  3000 IEL=NELEM+1,  NELEM+iNELEM

       BTMP=CMPLX(0.0D0,0.0D0)
       
       IF (ISOL.NE.1.AND.ISOL.NE.2) THEN
        PRINT*,"  Error: The input for ISOL should be either 1 or 2."
        STOP
       ENDIF

       ! Both ISOL=1 (total) and ISOL=2 (scattered) solve the same validated
       ! direct BIE for the scattered potential; the total/scattered choice is
       ! applied later by adding back the incident potential F0 (ISOL=1 only).
       DO 600 JEL=1,  NELEM

        DIST=SQRT((iXYZ_P(IEL-NELEM,1)-XYZ_P(JEL,1))**2+(iXYZ_P(IEL-NELEM,2)-XYZ_P(JEL,2))**2+(iXYZ_P(IEL-NELEM,3)-XYZ_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF

          TINRD=CMPLX(0.0D0, 0.0D0)

          DO  600   IS=1,  NSYS

            CALL DBC_IRR(IS,IEL-NELEM,JEL,TINRD,IRR,FLAG)

          DO IP=1, NSYS
            BTMP(IP)=BTMP(IP)+TINRD(IS,IP)
          ENDDO

600     CONTINUE
       
       DO  700  IP=1, NSYS

       DO  700  IS=1, NSYS
          ! Both ISX/ISY branches were identical — collapsed.
          BDMAT(IEL,IP)=BDMAT(IEL,IP)+RXY(IP,IS)*BTMP(IS)
700    CONTINUE

3000   CONTINUE

!$OMP END DO NOWAIT
!$OMP END PARALLEL
!
!--------------------------------------------------------------
!    Matrix set-up for the least-square problem

      DO IP=1, NSYS
       CALL NORMAL_RHS(AMAT(:,:,IP),TNELEM,NELEM,BDMAT(:,IP),1,DDMAT(:,IP))
      ENDDO

      RETURN
      END SUBROUTINE ASSB_DBC_IRR
      
! -----------------------------------------------------------------------------
!     Solve the radiation potentials after removing irregular frequencies
! -----------------------------------------------------------------------------

      SUBROUTINE RADIATION_SOLVER_IRR(CMAT,DRMAT,IPIV,MXPOT,NELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM,NSYS
      INTEGER,INTENT(IN):: IPIV(NELEM,NSYS)
      COMPLEX*16,INTENT(IN):: CMAT(NELEM,NELEM,NSYS),DRMAT(NELEM,6,NSYS)
      COMPLEX*16,INTENT(OUT):: MXPOT(NELEM,7,NSYS)
      
      INTEGER IEL,IS,IP,MD,INFO
      COMPLEX*16,ALLOCATABLE:: DRTMAT(:,:,:)
      ! ZGETRS does not modify its A argument — pass CMAT (the LU factorization) directly,
      ! eliminating the previous full-size CTMAT copy (Phase 2.4).

      ALLOCATE(DRTMAT(NELEM,6,NSYS))
      DRTMAT=DRMAT

      IF (ISOLV .EQ. 1) THEN
! Direct LU (ISOLV=1). NSYS is 1 or 2 — let MKL's internal threading do the work in ZGETRS.
         DO IP=1, NSYS
              CALL ZGETRS( 'No transpose', NELEM, 6, CMAT(:,:,IP), NELEM, IPIV(:,IP), DRTMAT(:,:,IP), NELEM, INFO )
         ENDDO
      ELSE
! GMRES (ISOLV=2) — iterate on the unfactored normal-equations matrix directly.
! H-matrix compression is skipped here because normal-equations blocks are not admissible.
         DO IP=1, NSYS
              CALL ZITER_SOLVE_MULTI(ISOLV, NELEM, CMAT(:,:,IP), DRMAT(:,:,IP), DRTMAT(:,:,IP), 6, INFO)
         ENDDO
      ENDIF
      
      DO 1400 MD=1, 6
      DO 1400 IP=1, NSYS
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,IS)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
      DO 1450 IEL=1, NELEM
         MXPOT(IEL,MD,IP)=CMPLX(0.0D0, 0.0D0)
      DO 1460 IS=1, NSYS
         ! Both ISX/ISY branches were identical — collapsed.
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+DRTMAT(IEL,MD,IS)*RXY(IP,IS)
1460  CONTINUE
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)/NSYS
1450  CONTINUE
!$OMP END DO NOWAIT
!$OMP END PARALLEL
1400  CONTINUE
      
      DEALLOCATE(DRTMAT)

      RETURN
      END SUBROUTINE RADIATION_SOLVER_IRR
      
! -----------------------------------------------------------------------------
!     Solve the diffraction potentials after removing irregular frequencies
! -----------------------------------------------------------------------------

      SUBROUTINE DIFFRACTION_SOLVER_IRR(CMAT,DDMAT,IPIV,MXPOT,NELEM,NSYS)
      IMPLICIT  NONE
!
      INTEGER,INTENT(IN):: NELEM,NSYS
      INTEGER,INTENT(IN):: IPIV(NELEM,NSYS)
      COMPLEX*16,INTENT(IN):: CMAT(NELEM,NELEM,NSYS),DDMAT(NELEM,NSYS)
      COMPLEX*16,INTENT(OUT):: MXPOT(NELEM,7,NSYS)
      
      INTEGER IEL,IS,IP,MD,INFO
      COMPLEX*16,ALLOCATABLE:: DDTMAT(:,:)
      COMPLEX*16,ALLOCATABLE:: BSING(:,:), XSING(:,:)
      ! ZGETRS does not modify its A argument — pass CMAT (the LU factorization) directly,
      ! eliminating the previous full-size CTMAT copy (Phase 2.4).

      ALLOCATE(DDTMAT(NELEM,NSYS))
      DDTMAT=DDMAT

      MD=7

      IF (ISOLV .EQ. 1) THEN
! Direct LU (ISOLV=1). NSYS is 1 or 2 — let MKL's internal threading do the work in ZGETRS.
         DO IP=1, NSYS
              CALL ZGETRS( 'No transpose', NELEM, 1, CMAT(:,:,IP), NELEM, IPIV(:,IP), DDTMAT(:,IP), NELEM, INFO )
         ENDDO
      ELSE
! GMRES (ISOLV=2) — iterate on the unfactored normal-equations matrix directly.
         ALLOCATE(BSING(NELEM,1), XSING(NELEM,1))
         DO IP=1, NSYS
              BSING(:,1) = DDMAT(:,IP)
              CALL ZITER_SOLVE_MULTI(ISOLV, NELEM, CMAT(:,:,IP), BSING, XSING, 1, INFO)
              DDTMAT(:,IP) = XSING(:,1)
         ENDDO
         DEALLOCATE(BSING, XSING)
      ENDIF
      
       DO 400 IP=1, NSYS
!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(IEL,IS)
! Row-disjoint writes on the IEL axis — no REDUCTION needed.
       DO 450 IEL=1, NELEM
         MXPOT(IEL,MD,IP)=CMPLX(0.0D0, 0.0D0)
       DO 460 IS=1, NSYS
         ! Both ISX/ISY branches were identical — collapsed.
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)+DDTMAT(IEL,IS)*RXY(IP,IS)
460    CONTINUE
         MXPOT(IEL,MD,IP)=MXPOT(IEL,MD,IP)/NSYS
450    CONTINUE
!$OMP END DO NOWAIT
!$OMP END PARALLEL
400    CONTINUE

      DEALLOCATE(DDTMAT)

      RETURN
      END SUBROUTINE DIFFRACTION_SOLVER_IRR

!-------------------------------------------------------------------------------
!     Normal equations of the overdetermined hull + lid system: C = A^H A over the N hull columns of A(1:M,1:N).
!     A is complex at a real frequency, and the least-squares solution, the minimum of ||A x - b||, solves
!     A^H A x = A^H b (Golub and Van Loan, Matrix Computations, section 5.3). The plain transpose of Liu (2019),
!     J. Mar. Sci. Eng. 7(3) 81, Eq. (27), is the same only for a real A, at the zero- and infinite-frequency limits.
!     Used by ASSB_LEFT_IRR and ASSB_LEFT_IRR_MULTI. Rows are summed in order, each (I,J) by one thread.
!-------------------------------------------------------------------------------
      SUBROUTINE NORMAL_MATRIX(A,M,N,C)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: M,N
      COMPLEX*16,INTENT(IN):: A(M,*)
      COMPLEX*16,INTENT(OUT):: C(N,N)
      INTEGER:: I,J,K
      COMPLEX*16:: S

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(I,J,K,S)
      DO I=1,N
       DO J=1,N
        S=CMPLX(0.0D0,0.0D0)
        DO K=1,M
         S=S+CONJG(A(K,I))*A(K,J)
        ENDDO
        C(I,J)=S
       ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      RETURN
      END SUBROUTINE NORMAL_MATRIX

!-------------------------------------------------------------------------------
!     Right sides of the normal equations: D = A^H B for NRHS columns of B(1:M,:), with the same A as NORMAL_MATRIX.
!     Used by the radiation and diffraction right sides of both paths. D is overwritten.
!-------------------------------------------------------------------------------
      SUBROUTINE NORMAL_RHS(A,M,N,B,NRHS,D)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: M,N,NRHS
      COMPLEX*16,INTENT(IN):: A(M,*),B(M,NRHS)
      COMPLEX*16,INTENT(OUT):: D(N,NRHS)
      INTEGER:: I,K,MD

!$OMP PARALLEL NUM_THREADS(NTHREAD)
!$OMP DO PRIVATE(I,K,MD)
      DO I=1,N
       DO MD=1,NRHS
        D(I,MD)=CMPLX(0.0D0,0.0D0)
        DO K=1,M
         D(I,MD)=D(I,MD)+B(K,MD)*CONJG(A(K,I))
        ENDDO
       ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

      RETURN
      END SUBROUTINE NORMAL_RHS

!-------------------------------------------------------------------------------
END MODULE AssbMatx_irr
!*******************************************************************************
