!  ------------------------------------------------------------------------------------------------------
!
!    Program HAMS-MREL for the diffraction and radiation of waves
!    by 3D structures.
!
!  License:
!
!    This routine is part of HAMS-MREL.
!
!    HAMS-MREL is a free software framework: you can redistribute it and/or modify it
!    under the terms of the Apache License, Version 2.0 (the "License"); you may
!    not use these subroutines except in compliance with the License. The software
!    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND.
!
!    You should have received a copy of the Apache License, Version 2.0, along with
!    HAMS-MREL. If not, see <http://www.apache.org/licenses/LICENSE-2.0>.
!
!  ------------------------------------------------------------------------------------------------------
!
!  Assembly and solve with dipole (zero-thickness) panels, first order, one body, no symmetry, direct LU
!  (notes/specs/2026-11-24-dipole-panels-design.md, Sec. 3 and 4). Used only when DipoleMesh.pnl gives DNELEM > 0; the
!  path without dipole panels (AssbMatx, AssbMatx_irr) is not touched.
!
!  Unknowns: [phi on the N = NELEM hull panels; mu = phi(+n) - phi(-n) on the M = DNELEM dipole panels], NTOT = N + M.
!    Hull row i, hull columns: as ASSB_LEFT (2 pi + BODINT_LEFT from RKBN and CGRN of CALGREEN).
!    Hull row i, dipole column j: -W_ij, W_ij = Omega_j(x_i) + s_V Omega_j(xbar_i) + n_Dj . grad_Q G_w ΔS_j    (5)
!    Dipole row i (times DPNSZ_i): hull columns H_ij, dipole columns -E_ij; right side sum_j v_j F_ij - 4 pi v_D,i (6)
!  with s_V = -1 at the infinite-frequency limit (V = -1) and +1 otherwise, Omega the solid angle and the Rankine parts
!  of H and E the vortex-ring gradients of DipoleIntgr, and the wave parts from INFGREEN3D / FINGREEN3D through
!  GREEN_WAVE_GRAD, GREEN_WAVE_FIELD_GRAD and GREEN_WAVE_NN (one point per source panel, as the hull rows).
!  F_ij = int_Sj dG/dn_F dS (hull panel j, derivative in the dipole collocation point) is kept for the right sides.
!-------------------------------------------------------------------------------
MODULE DipoleAssb

   USE HAMS_mod
   USE Body_mod, ONLY: XW
   USE Const_mod
   USE WaveDyn_mod
   USE PanelMesh_mod
   USE Potentials_mod
   USE DipoleMesh_mod
   USE DipoleIntgr
   USE BodyIntgr, ONLY: BODINT_LEFT, RBC_RIGHT, DBC_RIGHT
   USE SingularIntgr, ONLY: SGLINTBD_TRI2, SGLINTBD_QUAD2
   USE PatcVelct, ONLY: DINP
   USE omp_lib
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: DIPOLE_ASSB_LEFT
   PUBLIC :: DIPOLE_ASSB_RBC
   PUBLIC :: DIPOLE_ASSB_DBC
   PUBLIC :: DIPOLE_SOLVE
   PUBLIC :: DIPOLE_FORCE_SUM
   PUBLIC :: DIPOLE_RELEASE_ASSB

   COMPLEX*16,ALLOCATABLE:: FDIP(:,:)          ! F_ij (DNELEM, NELEM) of the last DIPOLE_ASSB_LEFT
   REAL*8,PARAMETER:: FOURPI=12.566370614359172D0

CONTAINS

!-------------------------------------------------------------------------------
!     s_V of the image term 1/r1: -1 at the infinite-frequency limit, +1 otherwise (design note item 1.2)
!-------------------------------------------------------------------------------
      REAL*8 FUNCTION IMAGE_SIGN()
      IMAGE_SIGN=1.D0
      IF (ABS(V+1.D0).LT.1.E-8) IMAGE_SIGN=-1.D0
      END FUNCTION IMAGE_SIGN

      SUBROUTINE HULL_NODES(J,A,NN)
      INTEGER,INTENT(IN):: J
      REAL*8,INTENT(OUT):: A(4,3)
      INTEGER,INTENT(OUT):: NN
      INTEGER:: K
      NN=NCN(J)
      A=0.D0
      DO K=1,NN
       A(K,:)=XYZ(NCON(J,K),:)
      ENDDO
      END SUBROUTINE HULL_NODES

      SUBROUTINE SHEET_NODES(J,A,NN)
      INTEGER,INTENT(IN):: J
      REAL*8,INTENT(OUT):: A(4,3)
      INTEGER,INTENT(OUT):: NN
      INTEGER:: K
      NN=DNCN(J)
      A=0.D0
      DO K=1,NN
       A(K,:)=DPXYZ(DNCON(J,K),:)
      ENDDO
      END SUBROUTINE SHEET_NODES

!-------------------------------------------------------------------------------
!     W of (5): the dipole column of a row at the point XF (a hull or lid point), dipole panel J.
!-------------------------------------------------------------------------------
      COMPLEX*16 FUNCTION SHEET_COLUMN(XF,J)
      REAL*8,INTENT(IN):: XF(3)
      INTEGER,INTENT(IN):: J
      REAL*8:: A(4,3),XI(3),NQ(3)
      INTEGER:: NN
      COMPLEX*16:: GW(4)

      CALL SHEET_NODES(J,A,NN)
      NQ(:)=DNRM(J,1:3)
      XI(1)=XF(1)
      XI(2)=XF(2)
      XI(3)=-XF(3)
      CALL GREEN_WAVE_GRAD(XF,DPXYZ_C(J,:),V,WVN,NK,H,GW)
      SHEET_COLUMN=DIPOLE_SOLID_ANGLE(XF,A,NN,NQ)+IMAGE_SIGN()*DIPOLE_SOLID_ANGLE(XI,A,NN,NQ)                    &
                  +(GW(2)*NQ(1)+GW(3)*NQ(2)+GW(4)*NQ(3))*DDS(J)
      END FUNCTION SHEET_COLUMN

!-------------------------------------------------------------------------------
!     The matrix on NTOT = NELEM + DNELEM unknowns, its LU factors (ZGETRF), and F_ij for the right sides. Needs
!     CALGREEN at this frequency (RKBN, CGRN of the hull rows).
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ASSB_LEFT(AMAT,IPIV,NTOT)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NTOT
      COMPLEX*16,INTENT(OUT):: AMAT(NTOT,NTOT)
      INTEGER,INTENT(OUT):: IPIV(NTOT)
      INTEGER:: IEL,JEL,FLAG,INFO,NN,ID
      REAL*8:: DIST,SV,XF(3),XI(3),NF(3),A(4,3),ENN,RK(4),RKI(4)
      COMPLEX*16:: TINDP(4),GNN,GW(4)

      IF (ALLOCATED(FDIP)) DEALLOCATE(FDIP)
      ALLOCATE(FDIP(DNELEM,NELEM))
      AMAT=CMPLX(0.D0,0.D0)
      SV=IMAGE_SIGN()

!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(IEL,JEL,FLAG,DIST,TINDP,XF) SCHEDULE(DYNAMIC,16)
      DO IEL=1,NELEM
       ! hull columns, as ASSB_LEFT with NSYS = 1
       AMAT(IEL,IEL)=2.D0*PI
       DO JEL=1,NELEM
        DIST=SQRT((XYZ_P(IEL,1)-XYZ_P(JEL,1))**2+(XYZ_P(IEL,2)-XYZ_P(JEL,2))**2+(XYZ_P(IEL,3)-XYZ_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
        TINDP=CMPLX(0.D0,0.D0)
        CALL BODINT_LEFT(1,IEL,JEL,TINDP,FLAG)
        AMAT(IEL,JEL)=AMAT(IEL,JEL)+RXY(1,1)*TINDP(1)
       ENDDO
       ! dipole columns
       XF(:)=XYZ_P(IEL,:)
       DO JEL=1,DNELEM
        AMAT(IEL,NELEM+JEL)=-SHEET_COLUMN(XF,JEL)
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(ID,JEL,XF,XI,NF,A,NN,ENN,GNN,RK,RKI,GW) SCHEDULE(DYNAMIC,4)
      DO ID=1,DNELEM
       XF(:)=DPXYZ_C(ID,:)
       XI(1)=XF(1)
       XI(2)=XF(2)
       XI(3)=-XF(3)
       NF(:)=DNRM(ID,1:3)
       ! hull columns: H_ij, and F_ij for the right sides
       DO JEL=1,NELEM
        CALL HULL_NODES(JEL,A,NN)
        CALL DIPOLE_RANKINE_NN(XF,NF,A,NN,GEONRM_P(JEL,:),SV,ENN)
        CALL GREEN_WAVE_NN(XF,NF,XYZ_P(JEL,:),GEONRM_P(JEL,:),V,WVN,NK,H,GNN)
        AMAT(NELEM+ID,JEL)=DPNSZ(ID)*(ENN+GNN*DS(JEL))
        IF (NN.EQ.3) THEN
         CALL SGLINTBD_TRI2(XF,JEL,RK,1)
         CALL SGLINTBD_TRI2(XI,JEL,RKI,1)
        ELSE
         CALL SGLINTBD_QUAD2(XF,JEL,RK,1)
         CALL SGLINTBD_QUAD2(XI,JEL,RKI,1)
        ENDIF
        CALL GREEN_WAVE_FIELD_GRAD(XF,XYZ_P(JEL,:),V,WVN,NK,H,GW)
        ! grad_x int 1/r = -RK(2:4); grad_x int 1/r1 = -R RKI(2:4), R = diag(1, 1, -1) (IMAGE_TERM_EXACT convention)
        FDIP(ID,JEL)=-(NF(1)*RK(2)+NF(2)*RK(3)+NF(3)*RK(4))-SV*(NF(1)*RKI(2)+NF(2)*RKI(3)-NF(3)*RKI(4))         &
                     +(GW(2)*NF(1)+GW(3)*NF(2)+GW(4)*NF(3))*DS(JEL)
       ENDDO
       ! dipole columns: -E_ij
       DO JEL=1,DNELEM
        CALL SHEET_NODES(JEL,A,NN)
        CALL DIPOLE_RANKINE_NN(XF,NF,A,NN,DNRM(JEL,1:3),SV,ENN)
        CALL GREEN_WAVE_NN(XF,NF,DPXYZ_C(JEL,:),DNRM(JEL,1:3),V,WVN,NK,H,GNN)
        AMAT(NELEM+ID,NELEM+JEL)=-DPNSZ(ID)*(ENN+GNN*DDS(JEL))
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL ZGETRF(NTOT,NTOT,AMAT,NTOT,IPIV,INFO)
      IF (INFO.NE.0) THEN
       PRINT *,' Error: the matrix with dipole panels is singular, ZGETRF INFO =',INFO
       STOP 1
      ENDIF

      RETURN
      END SUBROUTINE DIPOLE_ASSB_LEFT

!-------------------------------------------------------------------------------
!     Radiation right sides, modes 1-6: hull rows as ASSB_RBC (NSYS = 1); dipole rows (6), times DPNSZ.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ASSB_RBC(BRMAT,NTOT)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NTOT
      COMPLEX*16,INTENT(OUT):: BRMAT(NTOT,6)
      INTEGER:: IEL,JEL,MD,FLAG,ID
      REAL*8:: DIST
      COMPLEX*16:: TINRD(4,6,4),BTMP(6)

      BRMAT=CMPLX(0.D0,0.D0)
!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(IEL,JEL,MD,FLAG,DIST,TINRD,BTMP) SCHEDULE(DYNAMIC,16)
      DO IEL=1,NELEM
       BTMP=CMPLX(0.D0,0.D0)
       DO JEL=1,NELEM
        DIST=SQRT((XYZ_P(IEL,1)-XYZ_P(JEL,1))**2+(XYZ_P(IEL,2)-XYZ_P(JEL,2))**2+(XYZ_P(IEL,3)-XYZ_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
        TINRD=CMPLX(0.D0,0.D0)
        CALL RBC_RIGHT(1,IEL,JEL,TINRD,FLAG)
        BTMP(:)=BTMP(:)+TINRD(1,:,1)
       ENDDO
       BRMAT(IEL,:)=RXY(1,1)*BTMP(:)
      ENDDO
!$OMP END PARALLEL DO
      DO ID=1,DNELEM
       DO MD=1,6
        BTMP(MD)=CMPLX(0.D0,0.D0)
        DO JEL=1,NELEM
         BTMP(MD)=BTMP(MD)+FDIP(ID,JEL)*DXYZ_P(JEL,MD)
        ENDDO
        BRMAT(NELEM+ID,MD)=DPNSZ(ID)*(BTMP(MD)-FOURPI*DNRM(ID,MD))
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DIPOLE_ASSB_RBC

!-------------------------------------------------------------------------------
!     Diffraction right side (scattered potential): hull rows as ASSB_DBC (NSYS = 1); dipole rows (6) with
!     v = -dphi_I/dn, times DPNSZ.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ASSB_DBC(BDMAT,NTOT)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NTOT
      COMPLEX*16,INTENT(OUT):: BDMAT(NTOT)
      INTEGER:: IEL,JEL,FLAG,ID
      REAL*8:: DIST
      COMPLEX*16:: TINRD(4,4),BTMP,DPOX,DPOY,DPOZ,VH,VD

      BDMAT=CMPLX(0.D0,0.D0)
!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(IEL,JEL,FLAG,DIST,TINRD,BTMP) SCHEDULE(DYNAMIC,16)
      DO IEL=1,NELEM
       BTMP=CMPLX(0.D0,0.D0)
       DO JEL=1,NELEM
        DIST=SQRT((XYZ_P(IEL,1)-XYZ_P(JEL,1))**2+(XYZ_P(IEL,2)-XYZ_P(JEL,2))**2+(XYZ_P(IEL,3)-XYZ_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
        TINRD=CMPLX(0.D0,0.D0)
        CALL DBC_RIGHT(1,IEL,JEL,TINRD,FLAG)
        BTMP=BTMP+TINRD(1,1)
       ENDDO
       BDMAT(IEL)=RXY(1,1)*BTMP
      ENDDO
!$OMP END PARALLEL DO
      DO ID=1,DNELEM
       BTMP=CMPLX(0.D0,0.D0)
       DO JEL=1,NELEM
        CALL DINP(XYZ_P(JEL,1),XYZ_P(JEL,2),XYZ_P(JEL,3),XW(1),XW(2),BETA,DPOX,DPOY,DPOZ)
        VH=-(DPOX*GEONRM_P(JEL,1)+DPOY*GEONRM_P(JEL,2)+DPOZ*GEONRM_P(JEL,3))
        BTMP=BTMP+FDIP(ID,JEL)*VH
       ENDDO
       CALL DINP(DPXYZ_C(ID,1),DPXYZ_C(ID,2),DPXYZ_C(ID,3),XW(1),XW(2),BETA,DPOX,DPOY,DPOZ)
       VD=-(DPOX*DNRM(ID,1)+DPOY*DNRM(ID,2)+DPOZ*DNRM(ID,3))
       BDMAT(NELEM+ID)=DPNSZ(ID)*(BTMP-FOURPI*VD)
      ENDDO

      RETURN
      END SUBROUTINE DIPOLE_ASSB_DBC

!-------------------------------------------------------------------------------
!     Solve with the LU factors of DIPOLE_ASSB_LEFT: radiation (6 right sides) into MXPOT(1:NTOT, 1:6, 1), diffraction
!     into MXPOT(1:NTOT, 7, 1). NRHS = 6 or 1 picks which.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_SOLVE(AMAT,IPIV,B,NTOT,NRHS,MDFIRST)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NTOT,NRHS,MDFIRST
      COMPLEX*16,INTENT(IN):: AMAT(NTOT,NTOT)
      INTEGER,INTENT(IN):: IPIV(NTOT)
      COMPLEX*16,INTENT(IN):: B(NTOT,NRHS)
      COMPLEX*16,ALLOCATABLE:: X(:,:)
      INTEGER:: INFO

      ALLOCATE(X(NTOT,NRHS))
      X=B
      CALL ZGETRS('No transpose',NTOT,NRHS,AMAT,NTOT,IPIV,X,NTOT,INFO)
      MXPOT(1:NTOT,MDFIRST:MDFIRST+NRHS-1,1)=X(:,:)
      DEALLOCATE(X)

      RETURN
      END SUBROUTINE DIPOLE_SOLVE

!-------------------------------------------------------------------------------
!     The sheet's share of the force sums of RFORCE and EFORCE: F(k) = sum_m (-mu_m) nu_Dk,m ΔS_m for the potential in
!     MXPOT(:, MD, 1) (design note Sec. 4). Zero without dipole panels.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_FORCE_SUM(MD,F)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: MD
      COMPLEX*16,INTENT(OUT):: F(6)
      INTEGER:: M

      F=CMPLX(0.D0,0.D0)
      DO M=1,DNELEM
       F(:)=F(:)-MXPOT(NELEM+M,MD,1)*DDS(M)*DNRM(M,1:6)
      ENDDO

      RETURN
      END SUBROUTINE DIPOLE_FORCE_SUM

      SUBROUTINE DIPOLE_RELEASE_ASSB
      IF (ALLOCATED(FDIP)) DEALLOCATE(FDIP)
      END SUBROUTINE DIPOLE_RELEASE_ASSB
!-------------------------------------------------------------------------------
END MODULE DipoleAssb
!*******************************************************************************
