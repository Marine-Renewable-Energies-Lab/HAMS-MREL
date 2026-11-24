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
   USE BodyIntgr_irr, ONLY: BODINT_LEFT_IRR, RBC_IRR, DBC_IRR
   USE Inerfs_mod, ONLY: iXYZ_P
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
   PUBLIC :: DIPOLE_ASSB_LEFT_IRR
   PUBLIC :: DIPOLE_ASSB_RHS_IRR
   PUBLIC :: DIPOLE_LID_RESIDUAL
   PUBLIC :: DIPOLE_ROW_WEIGHTS

   COMPLEX*16,ALLOCATABLE:: FDIP(:,:)          ! F_ij (DNELEM, NELEM) of the last DIPOLE_ASSB_LEFT
   REAL*8,PARAMETER:: FOURPI=12.566370614359172D0
   REAL*8,ALLOCATABLE:: DIPOLE_ROW_WEIGHTS(:)  ! W of the lid least squares, rows 1..NTOT+NLID (design note Sec. 5)

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
      INTEGER:: INFO

      CALL FILL_SQUARE(AMAT,NTOT,NTOT,.FALSE.)
      CALL ZGETRF(NTOT,NTOT,AMAT,NTOT,IPIV,INFO)
      IF (INFO.NE.0) THEN
       PRINT *,' Error: the matrix with dipole panels is singular, ZGETRF INFO =',INFO
       STOP 1
      ENDIF

      RETURN
      END SUBROUTINE DIPOLE_ASSB_LEFT

!-------------------------------------------------------------------------------
!     Rows 1 to NTOT of the matrix (hull rows, then dipole rows) into AMAT(1:NTOT, 1:NTOT), AMAT dimensioned (LDA, NTOT),
!     and F_ij. IRRP: the hull block from CALGREEN_IRR (BODINT_LEFT_IRR, IRR = 1) instead of CALGREEN (BODINT_LEFT).
!-------------------------------------------------------------------------------
      SUBROUTINE FILL_SQUARE(AMAT,LDA,NTOT,IRRP)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: LDA,NTOT
      LOGICAL,INTENT(IN):: IRRP
      COMPLEX*16,INTENT(INOUT):: AMAT(LDA,NTOT)
      INTEGER:: IEL,JEL,FLAG,NN,ID
      REAL*8:: DIST,SV,XF(3),XI(3),NF(3),A(4,3),ENN,RK(4),RKI(4)
      COMPLEX*16:: TINDP(4),GNN,GW(4)

      IF (ALLOCATED(FDIP)) DEALLOCATE(FDIP)
      ALLOCATE(FDIP(DNELEM,NELEM))
      AMAT(1:NTOT,:)=CMPLX(0.D0,0.D0)
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
        IF (IRRP) THEN
         CALL BODINT_LEFT_IRR(1,IEL,JEL,TINDP,1,FLAG)
        ELSE
         CALL BODINT_LEFT(1,IEL,JEL,TINDP,FLAG)
        ENDIF
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

      RETURN
      END SUBROUTINE FILL_SQUARE

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
      CALL DIPOLE_ROWS_RBC(BRMAT)

      RETURN
      END SUBROUTINE DIPOLE_ASSB_RBC

!-------------------------------------------------------------------------------
!     The dipole rows NELEM+1..NELEM+DNELEM of the radiation right sides, (6) times DPNSZ, from FDIP.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ROWS_RBC(BRMAT)
      IMPLICIT NONE
      COMPLEX*16,INTENT(INOUT):: BRMAT(NELEM+DNELEM,6)
      INTEGER:: ID,JEL,MD
      COMPLEX*16:: BTMP

      DO ID=1,DNELEM
       DO MD=1,6
        BTMP=CMPLX(0.D0,0.D0)
        DO JEL=1,NELEM
         BTMP=BTMP+FDIP(ID,JEL)*DXYZ_P(JEL,MD)
        ENDDO
        BRMAT(NELEM+ID,MD)=DPNSZ(ID)*(BTMP-FOURPI*DNRM(ID,MD))
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DIPOLE_ROWS_RBC

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
      CALL DIPOLE_ROWS_DBC(BDMAT)

      RETURN
      END SUBROUTINE DIPOLE_ASSB_DBC

!-------------------------------------------------------------------------------
!     The dipole rows of the diffraction right side, (6) with v = -dphi_I/dn, times DPNSZ, from FDIP.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ROWS_DBC(BDMAT)
      IMPLICIT NONE
      COMPLEX*16,INTENT(INOUT):: BDMAT(NELEM+DNELEM)
      INTEGER:: ID,JEL
      COMPLEX*16:: BTMP,DPOX,DPOY,DPOZ,VH,VD

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
      END SUBROUTINE DIPOLE_ROWS_DBC

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

!-------------------------------------------------------------------------------
!     Weights of the rows of the lid least squares A^H W A x = A^H W b (design note Sec. 5): hull, dipole and lid rows,
!     length NTOT + NLID. On development_2 the IRSP = 1 solve without dipoles has no row weights, so W = 1 everywhere:
!     a dipole row, scaled by DPNSZ, then weighs like a hull row. (With the row weights of fix/field-limes-lid-weight
!     the dipole rows would take DDS/S on the scaled row.)
!-------------------------------------------------------------------------------
      SUBROUTINE SET_DIPOLE_ROW_WEIGHTS(NROW)
      INTEGER,INTENT(IN):: NROW
      IF (ALLOCATED(DIPOLE_ROW_WEIGHTS)) DEALLOCATE(DIPOLE_ROW_WEIGHTS)
      ALLOCATE(DIPOLE_ROW_WEIGHTS(NROW))
      DIPOLE_ROW_WEIGHTS=1.D0
      END SUBROUTINE SET_DIPOLE_ROW_WEIGHTS

      SUBROUTINE CHECK_DIPOLE_ROW_WEIGHTS(NROW)
      INTEGER,INTENT(IN):: NROW
      IF (.NOT.ALLOCATED(DIPOLE_ROW_WEIGHTS)) THEN
       PRINT *,' Error: the dipole row weights are not set.'
       STOP 1
      ENDIF
      IF (SIZE(DIPOLE_ROW_WEIGHTS).NE.NROW) THEN
       PRINT *,' Error: the dipole row weights have length',SIZE(DIPOLE_ROW_WEIGHTS),', not NTOT + NLID =',NROW
       STOP 1
      ENDIF
      END SUBROUTINE CHECK_DIPOLE_ROW_WEIGHTS

!-------------------------------------------------------------------------------
!     IRSP = 1 with dipole panels: the NTOT + NL rows (hull, dipole, then lid) into AFULL, and the normal matrix
!     C = A^H W A with its LU factors. The lid rows are those of ASSB_LEFT_IRR (BODINT_LEFT_IRR, IRR = 3) plus -W_lj (5)
!     in the dipole columns (design note Sec. 5). Needs CALGREEN_IRR at this frequency.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ASSB_LEFT_IRR(AFULL,CMAT,IPIV,NTOT,NL)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NTOT,NL
      COMPLEX*16,INTENT(OUT):: AFULL(NTOT+NL,NTOT),CMAT(NTOT,NTOT)
      INTEGER,INTENT(OUT):: IPIV(NTOT)
      INTEGER:: L,JEL,FLAG,INFO,I,J,K
      REAL*8:: DIST,XF(3)
      COMPLEX*16:: TINDP(4),S

      CALL FILL_SQUARE(AFULL,NTOT+NL,NTOT,.TRUE.)
!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(L,JEL,FLAG,DIST,TINDP,XF) SCHEDULE(DYNAMIC,16)
      DO L=1,NL
       XF(:)=iXYZ_P(L,:)
       DO JEL=1,NELEM
        DIST=SQRT((XF(1)-XYZ_P(JEL,1))**2+(XF(2)-XYZ_P(JEL,2))**2+(XF(3)-XYZ_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
        TINDP=CMPLX(0.D0,0.D0)
        CALL BODINT_LEFT_IRR(1,L,JEL,TINDP,3,FLAG)
        AFULL(NTOT+L,JEL)=RXY(1,1)*TINDP(1)
       ENDDO
       DO JEL=1,DNELEM
        AFULL(NTOT+L,NELEM+JEL)=-SHEET_COLUMN(XF,JEL)
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL SET_DIPOLE_ROW_WEIGHTS(NTOT+NL)
      CALL CHECK_DIPOLE_ROW_WEIGHTS(NTOT+NL)
!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(I,J,K,S)
      DO I=1,NTOT
       DO J=1,NTOT
        S=CMPLX(0.D0,0.D0)
        DO K=1,NTOT+NL
         S=S+CONJG(AFULL(K,I))*DIPOLE_ROW_WEIGHTS(K)*AFULL(K,J)
        ENDDO
        CMAT(I,J)=S
       ENDDO
      ENDDO
!$OMP END PARALLEL DO
      CALL ZGETRF(NTOT,NTOT,CMAT,NTOT,IPIV,INFO)
      IF (INFO.NE.0) THEN
       PRINT *,' Error: the normal matrix with dipole panels is singular, ZGETRF INFO =',INFO
       STOP 1
      ENDIF

      RETURN
      END SUBROUTINE DIPOLE_ASSB_LEFT_IRR

!-------------------------------------------------------------------------------
!     IRSP = 1 right sides B(1:NTOT+NL, NRHS): radiation (KIND = 1, NRHS = 6) or diffraction (KIND = 2, NRHS = 1).
!     Hull rows as RBC_IRR / DBC_IRR (IRR = 1), dipole rows as DIPOLE_ROWS_RBC / DBC, lid rows as RBC_IRR / DBC_IRR
!     (IRR = 3); then D = A^H W B. B stays with the caller for DIPOLE_LID_RESIDUAL.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ASSB_RHS_IRR(AFULL,B,D,NTOT,NL,NRHS,KIND)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NTOT,NL,NRHS,KIND
      COMPLEX*16,INTENT(IN):: AFULL(NTOT+NL,NTOT)
      COMPLEX*16,INTENT(OUT):: B(NTOT+NL,NRHS),D(NTOT,NRHS)
      INTEGER:: IEL,JEL,FLAG,IRR,IROW,I,K,MD
      REAL*8:: DIST,XF(3)
      COMPLEX*16:: TINRD(4,6,4),TINRDD(4,4),BTMP(6)
      COMPLEX*16,ALLOCATABLE:: BSQ(:,:)

      B=CMPLX(0.D0,0.D0)
!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(IROW,IEL,JEL,FLAG,IRR,DIST,XF,TINRD,TINRDD,BTMP) SCHEDULE(DYNAMIC,16)
      DO IROW=1,NELEM+NL
       IF (IROW.LE.NELEM) THEN
        IEL=IROW
        IRR=1
        XF(:)=XYZ_P(IEL,:)
       ELSE
        IEL=IROW-NELEM
        IRR=3
        XF(:)=iXYZ_P(IEL,:)
       ENDIF
       BTMP=CMPLX(0.D0,0.D0)
       DO JEL=1,NELEM
        DIST=SQRT((XF(1)-XYZ_P(JEL,1))**2+(XF(2)-XYZ_P(JEL,2))**2+(XF(3)-XYZ_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
        IF (KIND.EQ.1) THEN
         TINRD=CMPLX(0.D0,0.D0)
         CALL RBC_IRR(1,IEL,JEL,TINRD,IRR,FLAG)
         BTMP(:)=BTMP(:)+TINRD(1,:,1)
        ELSE
         TINRDD=CMPLX(0.D0,0.D0)
         CALL DBC_IRR(1,IEL,JEL,TINRDD,IRR,FLAG)
         BTMP(1)=BTMP(1)+TINRDD(1,1)
        ENDIF
       ENDDO
       IF (IROW.LE.NELEM) THEN
        B(IEL,1:NRHS)=RXY(1,1)*BTMP(1:NRHS)
       ELSE
        B(NTOT+IEL,1:NRHS)=RXY(1,1)*BTMP(1:NRHS)
       ENDIF
      ENDDO
!$OMP END PARALLEL DO
      ALLOCATE(BSQ(NTOT,NRHS))
      BSQ=CMPLX(0.D0,0.D0)
      IF (KIND.EQ.1) THEN
       CALL DIPOLE_ROWS_RBC(BSQ)
      ELSE
       CALL DIPOLE_ROWS_DBC(BSQ(:,1))
      ENDIF
      B(NELEM+1:NTOT,:)=BSQ(NELEM+1:NTOT,:)
      DEALLOCATE(BSQ)

      CALL CHECK_DIPOLE_ROW_WEIGHTS(NTOT+NL)
      DO MD=1,NRHS
       DO I=1,NTOT
        D(I,MD)=CMPLX(0.D0,0.D0)
        DO K=1,NTOT+NL
         D(I,MD)=D(I,MD)+CONJG(AFULL(K,I))*DIPOLE_ROW_WEIGHTS(K)*B(K,MD)
        ENDDO
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DIPOLE_ASSB_RHS_IRR

!-------------------------------------------------------------------------------
!     How well the least-squares solution X satisfies the lid rows (as the hub's J13 diagnostic): the rms over the lid
!     points of |sum_j A_lj X_j - b_l| / (sum_j |A_lj X_j| + |b_l|).
!-------------------------------------------------------------------------------
      REAL*8 FUNCTION DIPOLE_LID_RESIDUAL(AFULL,B,X,NTOT,NL)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NTOT,NL
      COMPLEX*16,INTENT(IN):: AFULL(NTOT+NL,NTOT),B(NTOT+NL),X(NTOT)
      INTEGER:: L,J
      COMPLEX*16:: R
      REAL*8:: T,S

      S=0.D0
      DO L=1,NL
       R=-B(NTOT+L)
       T=ABS(B(NTOT+L))
       DO J=1,NTOT
        R=R+AFULL(NTOT+L,J)*X(J)
        T=T+ABS(AFULL(NTOT+L,J)*X(J))
       ENDDO
       S=S+(ABS(R)/MAX(T,1.D-300))**2
      ENDDO
      DIPOLE_LID_RESIDUAL=SQRT(S/MAX(NL,1))

      RETURN
      END FUNCTION DIPOLE_LID_RESIDUAL

      SUBROUTINE DIPOLE_RELEASE_ASSB
      IF (ALLOCATED(FDIP)) DEALLOCATE(FDIP)
      IF (ALLOCATED(DIPOLE_ROW_WEIGHTS)) DEALLOCATE(DIPOLE_ROW_WEIGHTS)
      END SUBROUTINE DIPOLE_RELEASE_ASSB
!-------------------------------------------------------------------------------
END MODULE DipoleAssb
!*******************************************************************************
