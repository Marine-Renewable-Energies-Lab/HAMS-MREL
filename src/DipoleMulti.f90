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
!  Dipole (zero-thickness) panels on several bodies (NBODY > 1), first order, no symmetry, direct LU: the multi-body
!  twin of DipoleAssb (docs/developer_documentation.md, "Dipole panels on several bodies"). Used only when a
!  HullMesh_N.pnl has dipole panels (DNELEM_TOTAL > 0); the multi-body path without them (AssbMatxMulti,
!  AssbMatx_irrMulti) is not touched.
!
!  Unknowns: [phi on the NH hull panels of all bodies; mu on the ND dipole panels of all bodies], NTOT = NH + ND, in the
!  order of the combined arrays (body by body). The rows are those of DipoleAssb on the combined arrays:
!    hull row i, hull columns: 2 pi + BODINT_MULTI_LEFT (or BODINT_LEFT_IRR_MULTI, IRR = 1, with the lid);
!    hull row i, dipole column j: -W_ij (solid angle, image and wave gradient of sheet panel j at x_i);
!    dipole row i (times DPNSZ_i): hull columns H_ij, dipole columns -E_ij.
!  Radiation of body b, mode k: the normal velocity is n_k on b's hull panels and nu_Dk on b's sheet panels, zero on
!  the other bodies, so dipole row i gets DPNSZ_i (sum over b's hull panels j of F_ij n_jk - 4 pi [owner(i) = b] nu_Dk,i).
!  F_ij = int_Sj dG/dn_F dS (hull panel j, derivative at the collocation point of dipole panel i).
!-------------------------------------------------------------------------------
MODULE DipoleMulti

   USE HAMS_mod
   USE Body_mod, ONLY: XW,XR_MULTI
   USE Const_mod
   USE WaveDyn_mod
   USE PanelMesh_mod
   USE Inerfs_mod, ONLY: tNELEM
   USE Potentials_mod, ONLY: MXPOT_MULTI_COMB
   USE DipolePanelMulti_mod
   USE DipoleIntgr, ONLY: DIPOLE_SOLID_ANGLE,DIPOLE_RANKINE_NN,GREEN_WAVE_GRAD,GREEN_WAVE_FIELD_GRAD,GREEN_WAVE_NN
   USE BodyIntgrMulti, ONLY: BODINT_MULTI_LEFT,RBC_MULTI_RIGHT,DBC_MULTI_RIGHT
   USE BodyIntgr_irrMulti, ONLY: BODINT_LEFT_IRR_MULTI,RBC_IRR_MULTI,DBC_IRR_MULTI
   USE SingularIntgrMulti, ONLY: SGLINTBD_MULTI_TRI2,SGLINTBD_MULTI_QUAD2
   USE PatcVelct, ONLY: DINP
   USE NormalProcess, ONLY: CalPanelCentre,CalPanelArea,CalPanelChartLength,CalTransNormals,CalRotNormals
   USE ReadPanelMesh, ONLY: DIPOLE_CHECK_GEOMETRY,DIPOLE_CHECK_CLEARANCE,DIPOLE_CHECK_COPLANAR
   USE omp_lib
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: DIPOLE_SETUP_MULTI
   PUBLIC :: DIPOLE_CHECK_OPTIONS_MULTI
   PUBLIC :: DIPOLE_ASSB_LEFT_MULTI
   PUBLIC :: DIPOLE_ASSB_RBC_MULTI
   PUBLIC :: DIPOLE_ASSB_DBC_MULTI
   PUBLIC :: DIPOLE_SOLVE_MULTI
   PUBLIC :: DIPOLE_ASSB_LEFT_IRR_MULTI
   PUBLIC :: DIPOLE_ASSB_RHS_IRR_MULTI
   PUBLIC :: DIPOLE_FORCE_SUM_MULTI
   PUBLIC :: DIPOLE_RELEASE_MULTI
   PUBLIC :: DIPOLE_NH_MULTI

   COMPLEX*16,ALLOCATABLE:: FDIPM(:,:)          ! F_ij (ND, NH) of the last DIPOLE_ASSB_LEFT_MULTI or _LEFT_IRR_MULTI
   REAL*8,PARAMETER:: FOURPI=12.566370614359172D0

CONTAINS

!-------------------------------------------------------------------------------
!     NH: the hull panels of all bodies (the hull unknowns come first, then the dipole jumps).
!-------------------------------------------------------------------------------
      INTEGER FUNCTION DIPOLE_NH_MULTI()
      DIPOLE_NH_MULTI=SUM(NELEM_MULTI(1:NBODY))
      END FUNCTION DIPOLE_NH_MULTI

!-------------------------------------------------------------------------------
!     First and last hull panel of body B in the combined arrays.
!-------------------------------------------------------------------------------
      SUBROUTINE HULL_RANGE(B,J0,J1)
      INTEGER,INTENT(IN):: B
      INTEGER,INTENT(OUT):: J0,J1
      J0=SUM(NELEM_MULTI(1:B-1))+1
      J1=J0+NELEM_MULTI(B)-1
      END SUBROUTINE HULL_RANGE

!-------------------------------------------------------------------------------
!     s_V of the image term 1/r1: -1 at the infinite-frequency limit, +1 otherwise (as DipoleAssb).
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
      NN=NCN_MULTI_COMB(J)
      A=0.D0
      DO K=1,NN
       A(K,:)=XYZ_GLOBAL_MULTI_COMB(NCON_MULTI_COMB(J,K),:)
      ENDDO
      END SUBROUTINE HULL_NODES

      SUBROUTINE SHEET_NODES(J,A,NN)
      INTEGER,INTENT(IN):: J
      REAL*8,INTENT(OUT):: A(4,3)
      INTEGER,INTENT(OUT):: NN
      INTEGER:: K
      NN=DNCN_COMB(J)
      A=0.D0
      DO K=1,NN
       A(K,:)=DPXYZ_COMB(DNCON_COMB(J,K),:)
      ENDDO
      END SUBROUTINE SHEET_NODES

!-------------------------------------------------------------------------------
!     W_ij: the dipole column of a row at the point XF (a hull or lid point), dipole panel J.
!-------------------------------------------------------------------------------
      COMPLEX*16 FUNCTION SHEET_COLUMN(XF,J)
      REAL*8,INTENT(IN):: XF(3)
      INTEGER,INTENT(IN):: J
      REAL*8:: A(4,3),XI(3),NQ(3)
      INTEGER:: NN
      COMPLEX*16:: GW(4)

      CALL SHEET_NODES(J,A,NN)
      NQ(:)=DNRM_COMB(J,1:3)
      XI(1)=XF(1)
      XI(2)=XF(2)
      XI(3)=-XF(3)
      CALL GREEN_WAVE_GRAD(XF,DPXYZ_C_COMB(J,:),V,WVN,NK,H,GW)
      SHEET_COLUMN=DIPOLE_SOLID_ANGLE(XF,A,NN,NQ)+IMAGE_SIGN()*DIPOLE_SOLID_ANGLE(XI,A,NN,NQ)                    &
                  +(GW(2)*NQ(1)+GW(3)*NQ(2)+GW(4)*NQ(3))*DDS_COMB(J)
      END FUNCTION SHEET_COLUMN

!=======================================================================
!     The limits of the multi-body dipole path: no symmetry, direct LU, first order only, no field points. OK false and
!     a one-line MSG naming the option and the first body with dipole panels otherwise.
!=======================================================================
      SUBROUTINE DIPOLE_CHECK_OPTIONS_MULTI(NSYS_IN,ISOLV_IN,SO_ACTIVE,NFP_IN,OK,MSG)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSYS_IN,ISOLV_IN,NFP_IN
      LOGICAL,INTENT(IN):: SO_ACTIVE
      LOGICAL,INTENT(OUT):: OK
      CHARACTER(LEN=*),INTENT(OUT):: MSG
      INTEGER:: B,B1
      CHARACTER(LEN=24):: WHO

      B1=0
      IF (ALLOCATED(DNELEM_MULTI)) THEN
       DO B=SIZE(DNELEM_MULTI),1,-1
        IF (DNELEM_MULTI(B).GT.0) B1=B
       ENDDO
      ENDIF
      WRITE(WHO,'(A,I0,A)') ' (body ',B1,')'
      OK=.FALSE.
      IF (NSYS_IN.NE.1) THEN
       MSG=' Error: dipole panels'//TRIM(WHO)//' do not support symmetry: set the hull symmetry flags to 0.'
      ELSEIF (ISOLV_IN.NE.1) THEN
       MSG=' Error: dipole panels'//TRIM(WHO)//' need the direct solver, ISOLV = 1.'
      ELSEIF (SO_ACTIVE) THEN
       MSG=' Error: dipole panels'//TRIM(WHO)//' are first order only: remove SecondOrder.in or switch it off.'
      ELSEIF (NFP_IN.GT.0) THEN
       MSG=' Error: field-point pressures and elevations are not supported with dipole panels'//TRIM(WHO)//'.'
      ELSE
       OK=.TRUE.
       MSG=''
      ENDIF

      RETURN
      END SUBROUTINE DIPOLE_CHECK_OPTIONS_MULTI

!=======================================================================
!     The geometry of the sheets of all bodies from the dipole blocks (DIPOLE_READ_PANELS_MULTI left their node numbers
!     in DNCON_COMB, numbers of the combined node list XYZ_GLOBAL_MULTI_COMB). Call after BuildCombinedArrays. Without
!     dipole panels (DNELEM_TOTAL = 0) nothing is done and OK is true. The nodes the sheets use are copied out, in the
!     order of the combined node list, into DPXYZ_COMB; normals follow the node order as for hull panels, rotational
!     normals are about the owning body's XR_MULTI. The checks of the single body (submerged, clear of z = 0, not on a
!     hull panel) run on all sheets against all hulls.
!=======================================================================
      SUBROUTINE DIPOLE_SETUP_MULTI(NSYS_IN,ISOLV_IN,SO_ACTIVE,NFP_IN,OK,MSG)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NSYS_IN,ISOLV_IN,NFP_IN
      LOGICAL,INTENT(IN):: SO_ACTIVE
      LOGICAL,INTENT(OUT):: OK
      CHARACTER(LEN=*),INTENT(OUT):: MSG
      INTEGER,ALLOCATABLE:: NEWNUM(:)
      INTEGER:: J,IND,IEL,NND,B,NH,I0,N

      OK=.TRUE.
      MSG=''
      IF (DNELEM_TOTAL.LE.0) RETURN

      CALL DIPOLE_CHECK_OPTIONS_MULTI(NSYS_IN,ISOLV_IN,SO_ACTIVE,NFP_IN,OK,MSG)
      IF (.NOT.OK) RETURN

      NND=SIZE(XYZ_GLOBAL_MULTI_COMB,1)
      NH=DIPOLE_NH_MULTI()
      ALLOCATE(NEWNUM(NND))
      NEWNUM=0
      DO IEL=1,DNELEM_TOTAL
       IF (DNCN_COMB(IEL).NE.3.AND.DNCN_COMB(IEL).NE.4) THEN
        OK=.FALSE.
        WRITE(MSG,'(A,I0,A,I0,A)') ' Error: dipole panel ',IEL-DSTART(DOWNER(IEL))+1,' of HullMesh_',DOWNER(IEL),         &
                                   '.pnl does not have 3 or 4 nodes.'
        RETURN
       ENDIF
       DO J=1,DNCN_COMB(IEL)
        IF (DNCON_COMB(IEL,J).LT.1.OR.DNCON_COMB(IEL,J).GT.NND) THEN
         OK=.FALSE.
         WRITE(MSG,'(A,I0,A,I0,A)') ' Error: dipole panel ',IEL-DSTART(DOWNER(IEL))+1,' of HullMesh_',DOWNER(IEL),        &
                                    '.pnl refers to a node outside the node list.'
         RETURN
        ENDIF
        NEWNUM(DNCON_COMB(IEL,J))=1
       ENDDO
      ENDDO
      DNTND_TOTAL=0
      DO IND=1,NND
       IF (NEWNUM(IND).GT.0) THEN
        DNTND_TOTAL=DNTND_TOTAL+1
        NEWNUM(IND)=DNTND_TOTAL
       ENDIF
      ENDDO
      IF (ALLOCATED(DPXYZ_COMB)) DEALLOCATE(DPXYZ_COMB)
      ALLOCATE(DPXYZ_COMB(DNTND_TOTAL,3))
      DO IND=1,NND
       IF (NEWNUM(IND).GT.0) DPXYZ_COMB(NEWNUM(IND),1:3)=XYZ_GLOBAL_MULTI_COMB(IND,1:3)
      ENDDO
      DO IEL=1,DNELEM_TOTAL
       DO J=1,DNCN_COMB(IEL)
        DNCON_COMB(IEL,J)=NEWNUM(DNCON_COMB(IEL,J))
       ENDDO
      ENDDO
      DEALLOCATE(NEWNUM)
      IF (ALLOCATED(DPXYZ_C_COMB)) DEALLOCATE(DPXYZ_C_COMB)
      IF (ALLOCATED(DNRM_COMB)) DEALLOCATE(DNRM_COMB)
      IF (ALLOCATED(DDS_COMB)) DEALLOCATE(DDS_COMB)
      IF (ALLOCATED(DPNSZ_COMB)) DEALLOCATE(DPNSZ_COMB)
      ALLOCATE(DPXYZ_C_COMB(DNELEM_TOTAL,3),DNRM_COMB(DNELEM_TOTAL,6),DDS_COMB(DNELEM_TOTAL),DPNSZ_COMB(DNELEM_TOTAL))

      CALL CalPanelCentre(DPXYZ_COMB,DNTND_TOTAL,DNELEM_TOTAL,DNCN_COMB,DNCON_COMB,DPXYZ_C_COMB)
      CALL CalPanelArea(DPXYZ_COMB,DNTND_TOTAL,DNELEM_TOTAL,DNCN_COMB,DNCON_COMB,DDS_COMB)
      CALL CalPanelChartLength(DPXYZ_COMB,DNTND_TOTAL,DNELEM_TOTAL,DNCN_COMB,DNCON_COMB,DPNSZ_COMB)
      CALL CalTransNormals(DPXYZ_COMB,DNTND_TOTAL,DNELEM_TOTAL,DNCN_COMB,DNCON_COMB,DNRM_COMB)
      DO B=1,NBODY
       N=DNELEM_MULTI(B)
       IF (N.GT.0) THEN
        I0=DSTART(B)
        CALL CalRotNormals(XR_MULTI(B,:),DPXYZ_C_COMB(I0:I0+N-1,:),DNRM_COMB(I0:I0+N-1,:),N)
       ENDIF
      ENDDO

      CALL DIPOLE_CHECK_GEOMETRY(DPXYZ_COMB,DNTND_TOTAL,DPXYZ_C_COMB,DPNSZ_COMB,DNELEM_TOTAL,XYZ_GLOBAL_MULTI_COMB_P,NH,  &
                                 OK,MSG)
      IF (.NOT.OK) RETURN
      CALL DIPOLE_CHECK_CLEARANCE(DPXYZ_COMB,DNCN_COMB,DNCON_COMB,DPNSZ_COMB,DNELEM_TOTAL,OK,MSG)
      IF (.NOT.OK) RETURN
      CALL DIPOLE_CHECK_COPLANAR(DPXYZ_COMB,DNCN_COMB,DNCON_COMB,DPXYZ_C_COMB,DNRM_COMB(:,1:3),DPNSZ_COMB,DNELEM_TOTAL, &
                                 XYZ_GLOBAL_MULTI_COMB,NCN_MULTI_COMB,NCON_MULTI_COMB,XYZ_GLOBAL_MULTI_COMB_P,             &
                                 GEONRM_MULTI_COMB,PNSZ_MULTI_COMB,NH,OK,MSG)
      IF (.NOT.OK) RETURN

      DO B=1,NBODY
       IF (DNELEM_MULTI(B).GT.0) PRINT *,' Dipole panels read from HullMesh_',B,'.pnl:',DNELEM_MULTI(B)
      ENDDO
      RETURN
      END SUBROUTINE DIPOLE_SETUP_MULTI

!-------------------------------------------------------------------------------
!     The matrix on NTOT = NH + ND unknowns and its LU factors (ZGETRF), and F_ij for the right sides. Needs
!     CALGREEN_MULTI at this frequency.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ASSB_LEFT_MULTI(AMAT,IPIV,NTOT)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NTOT
      COMPLEX*16,INTENT(OUT):: AMAT(NTOT,NTOT)
      INTEGER,INTENT(OUT):: IPIV(NTOT)
      INTEGER:: INFO

      CALL FILL_SQUARE(AMAT,NTOT,NTOT,.FALSE.)
      CALL ZGETRF(NTOT,NTOT,AMAT,NTOT,IPIV,INFO)
      IF (INFO.NE.0) THEN
       PRINT *,' Error: the matrix with dipole panels of several bodies is singular, ZGETRF INFO =',INFO
       STOP 1
      ENDIF

      RETURN
      END SUBROUTINE DIPOLE_ASSB_LEFT_MULTI

!-------------------------------------------------------------------------------
!     Rows 1 to NTOT (hull rows of all bodies, then dipole rows) into AMAT(1:NTOT, 1:NTOT), AMAT dimensioned (LDA, NTOT),
!     and F_ij. IRRP: the hull block from CALGREEN_MULTI_IRR (BODINT_LEFT_IRR_MULTI, IRR = 1).
!-------------------------------------------------------------------------------
      SUBROUTINE FILL_SQUARE(AMAT,LDA,NTOT,IRRP)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: LDA,NTOT
      LOGICAL,INTENT(IN):: IRRP
      COMPLEX*16,INTENT(INOUT):: AMAT(LDA,NTOT)
      INTEGER:: IEL,JEL,FLAG,NN,ID,NH,ND
      REAL*8:: DIST,SV,XF(3),XI(3),NF(3),A(4,3),ENN,RK(4),RKI(4)
      COMPLEX*16:: TINDP(4),GNN,GW(4)

      NH=DIPOLE_NH_MULTI()
      ND=DNELEM_TOTAL
      IF (ALLOCATED(FDIPM)) DEALLOCATE(FDIPM)
      ALLOCATE(FDIPM(ND,NH))
      AMAT(1:NTOT,:)=CMPLX(0.D0,0.D0)
      SV=IMAGE_SIGN()

!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(IEL,JEL,FLAG,DIST,TINDP,XF) SCHEDULE(DYNAMIC,16)
      DO IEL=1,NH
       ! hull columns, as ASSB_LEFT_MULTI with NSYS = 1
       AMAT(IEL,IEL)=2.D0*PI
       DO JEL=1,NH
        DIST=SQRT((XYZ_GLOBAL_MULTI_COMB_P(IEL,1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2                                &
                 +(XYZ_GLOBAL_MULTI_COMB_P(IEL,2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2                                &
                 +(XYZ_GLOBAL_MULTI_COMB_P(IEL,3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
        TINDP=CMPLX(0.D0,0.D0)
        IF (IRRP) THEN
         CALL BODINT_LEFT_IRR_MULTI(1,IEL,JEL,TINDP,1,FLAG)
        ELSE
         CALL BODINT_MULTI_LEFT(1,IEL,JEL,TINDP,FLAG)
        ENDIF
        AMAT(IEL,JEL)=AMAT(IEL,JEL)+RXY_MULTI(1,1)*TINDP(1)
       ENDDO
       ! dipole columns
       XF(:)=XYZ_GLOBAL_MULTI_COMB_P(IEL,:)
       DO JEL=1,ND
        AMAT(IEL,NH+JEL)=-SHEET_COLUMN(XF,JEL)
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(ID,JEL,XF,XI,NF,A,NN,ENN,GNN,RK,RKI,GW) SCHEDULE(DYNAMIC,4)
      DO ID=1,ND
       XF(:)=DPXYZ_C_COMB(ID,:)
       XI(1)=XF(1)
       XI(2)=XF(2)
       XI(3)=-XF(3)
       NF(:)=DNRM_COMB(ID,1:3)
       ! hull columns: H_ij, and F_ij for the right sides
       DO JEL=1,NH
        CALL HULL_NODES(JEL,A,NN)
        CALL DIPOLE_RANKINE_NN(XF,NF,A,NN,GEONRM_MULTI_COMB(JEL,:),SV,ENN)
        CALL GREEN_WAVE_NN(XF,NF,XYZ_GLOBAL_MULTI_COMB_P(JEL,:),GEONRM_MULTI_COMB(JEL,:),V,WVN,NK,H,GNN)
        AMAT(NH+ID,JEL)=DPNSZ_COMB(ID)*(ENN+GNN*DS_MULTI_COMB(JEL))
        IF (NN.EQ.3) THEN
         CALL SGLINTBD_MULTI_TRI2(XF,JEL,RK,1)
         CALL SGLINTBD_MULTI_TRI2(XI,JEL,RKI,1)
        ELSE
         CALL SGLINTBD_MULTI_QUAD2(XF,JEL,RK,1)
         CALL SGLINTBD_MULTI_QUAD2(XI,JEL,RKI,1)
        ENDIF
        CALL GREEN_WAVE_FIELD_GRAD(XF,XYZ_GLOBAL_MULTI_COMB_P(JEL,:),V,WVN,NK,H,GW)
        ! grad_x int 1/r = -RK(2:4); grad_x int 1/r1 = -R RKI(2:4), R = diag(1, 1, -1) (as DipoleAssb)
        FDIPM(ID,JEL)=-(NF(1)*RK(2)+NF(2)*RK(3)+NF(3)*RK(4))-SV*(NF(1)*RKI(2)+NF(2)*RKI(3)-NF(3)*RKI(4))        &
                      +(GW(2)*NF(1)+GW(3)*NF(2)+GW(4)*NF(3))*DS_MULTI_COMB(JEL)
       ENDDO
       ! dipole columns: -E_ij
       DO JEL=1,ND
        CALL SHEET_NODES(JEL,A,NN)
        CALL DIPOLE_RANKINE_NN(XF,NF,A,NN,DNRM_COMB(JEL,1:3),SV,ENN)
        CALL GREEN_WAVE_NN(XF,NF,DPXYZ_C_COMB(JEL,:),DNRM_COMB(JEL,1:3),V,WVN,NK,H,GNN)
        AMAT(NH+ID,NH+JEL)=-DPNSZ_COMB(ID)*(ENN+GNN*DDS_COMB(JEL))
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      RETURN
      END SUBROUTINE FILL_SQUARE

!-------------------------------------------------------------------------------
!     Radiation right sides of body IB, its modes 1-6, into BRMAT(1:NTOT, 1:6): hull rows of all bodies as ASSB_RBC_MULTI
!     on IB's hull panels (NSYS = 1); dipole rows as in the module header.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ASSB_RBC_MULTI(BRMAT,NTOT,IB)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NTOT,IB
      COMPLEX*16,INTENT(OUT):: BRMAT(NTOT,6)
      INTEGER:: IEL,JEL,FLAG,NH,J0,J1
      REAL*8:: DIST
      COMPLEX*16:: TINRD(4,6,4),BTMP(6)

      NH=DIPOLE_NH_MULTI()
      CALL HULL_RANGE(IB,J0,J1)
      BRMAT=CMPLX(0.D0,0.D0)
!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(IEL,JEL,FLAG,DIST,TINRD,BTMP) SCHEDULE(DYNAMIC,16)
      DO IEL=1,NH
       BTMP=CMPLX(0.D0,0.D0)
       DO JEL=J0,J1
        DIST=SQRT((XYZ_GLOBAL_MULTI_COMB_P(IEL,1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2                                &
                 +(XYZ_GLOBAL_MULTI_COMB_P(IEL,2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2                                &
                 +(XYZ_GLOBAL_MULTI_COMB_P(IEL,3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
        TINRD=CMPLX(0.D0,0.D0)
        CALL RBC_MULTI_RIGHT(1,IEL,JEL,TINRD,FLAG)
        BTMP(:)=BTMP(:)+TINRD(1,:,1)
       ENDDO
       BRMAT(IEL,:)=RXY_MULTI(1,1)*BTMP(:)
      ENDDO
!$OMP END PARALLEL DO
      CALL DIPOLE_ROWS_RBC(BRMAT,NH,IB)

      RETURN
      END SUBROUTINE DIPOLE_ASSB_RBC_MULTI

!-------------------------------------------------------------------------------
!     The dipole rows NH+1..NH+ND of the radiation right sides of body IB, times DPNSZ, from FDIPM: the F sum runs over
!     IB's hull panels; -4 pi nu_D only on IB's own sheet panels.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ROWS_RBC(BRMAT,NH,IB)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NH,IB
      COMPLEX*16,INTENT(INOUT):: BRMAT(NH+DNELEM_TOTAL,6)
      INTEGER:: ID,JEL,MD,J0,J1
      COMPLEX*16:: BTMP

      CALL HULL_RANGE(IB,J0,J1)
      DO ID=1,DNELEM_TOTAL
       DO MD=1,6
        BTMP=CMPLX(0.D0,0.D0)
        DO JEL=J0,J1
         BTMP=BTMP+FDIPM(ID,JEL)*DXYZ_MULTI_COMB(JEL,MD)
        ENDDO
        IF (DOWNER(ID).EQ.IB) BTMP=BTMP-FOURPI*DNRM_COMB(ID,MD)
        BRMAT(NH+ID,MD)=DPNSZ_COMB(ID)*BTMP
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DIPOLE_ROWS_RBC

!-------------------------------------------------------------------------------
!     Diffraction right side (scattered potential): hull rows as ASSB_DBC_MULTI (NSYS = 1); dipole rows with
!     v = -dphi_I/dn on all hulls and sheets, times DPNSZ.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ASSB_DBC_MULTI(BDMAT,NTOT)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NTOT
      COMPLEX*16,INTENT(OUT):: BDMAT(NTOT)
      INTEGER:: IEL,JEL,FLAG,NH
      REAL*8:: DIST
      COMPLEX*16:: TINRD(4,4),BTMP

      NH=DIPOLE_NH_MULTI()
      BDMAT=CMPLX(0.D0,0.D0)
!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(IEL,JEL,FLAG,DIST,TINRD,BTMP) SCHEDULE(DYNAMIC,16)
      DO IEL=1,NH
       BTMP=CMPLX(0.D0,0.D0)
       DO JEL=1,NH
        DIST=SQRT((XYZ_GLOBAL_MULTI_COMB_P(IEL,1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2                                &
                 +(XYZ_GLOBAL_MULTI_COMB_P(IEL,2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2                                &
                 +(XYZ_GLOBAL_MULTI_COMB_P(IEL,3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
        TINRD=CMPLX(0.D0,0.D0)
        CALL DBC_MULTI_RIGHT(1,IEL,JEL,TINRD,FLAG)
        BTMP=BTMP+TINRD(1,1)
       ENDDO
       BDMAT(IEL)=RXY_MULTI(1,1)*BTMP
      ENDDO
!$OMP END PARALLEL DO
      CALL DIPOLE_ROWS_DBC(BDMAT,NH)

      RETURN
      END SUBROUTINE DIPOLE_ASSB_DBC_MULTI

!-------------------------------------------------------------------------------
!     The dipole rows of the diffraction right side, v = -dphi_I/dn, times DPNSZ, from FDIPM.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ROWS_DBC(BDMAT,NH)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NH
      COMPLEX*16,INTENT(INOUT):: BDMAT(NH+DNELEM_TOTAL)
      INTEGER:: ID,JEL
      COMPLEX*16:: BTMP,DPOX,DPOY,DPOZ,VD
      COMPLEX*16,ALLOCATABLE:: VH(:)

      ALLOCATE(VH(NH))
      DO JEL=1,NH
       CALL DINP(XYZ_GLOBAL_MULTI_COMB_P(JEL,1),XYZ_GLOBAL_MULTI_COMB_P(JEL,2),XYZ_GLOBAL_MULTI_COMB_P(JEL,3),         &
                 XW(1),XW(2),BETA,DPOX,DPOY,DPOZ)
       VH(JEL)=-(DPOX*GEONRM_MULTI_COMB(JEL,1)+DPOY*GEONRM_MULTI_COMB(JEL,2)+DPOZ*GEONRM_MULTI_COMB(JEL,3))
      ENDDO
      DO ID=1,DNELEM_TOTAL
       BTMP=CMPLX(0.D0,0.D0)
       DO JEL=1,NH
        BTMP=BTMP+FDIPM(ID,JEL)*VH(JEL)
       ENDDO
       CALL DINP(DPXYZ_C_COMB(ID,1),DPXYZ_C_COMB(ID,2),DPXYZ_C_COMB(ID,3),XW(1),XW(2),BETA,DPOX,DPOY,DPOZ)
       VD=-(DPOX*DNRM_COMB(ID,1)+DPOY*DNRM_COMB(ID,2)+DPOZ*DNRM_COMB(ID,3))
       BDMAT(NH+ID)=DPNSZ_COMB(ID)*(BTMP-FOURPI*VD)
      ENDDO
      DEALLOCATE(VH)

      RETURN
      END SUBROUTINE DIPOLE_ROWS_DBC

!-------------------------------------------------------------------------------
!     Solve with the LU factors (DIPOLE_ASSB_LEFT_MULTI, or the normal matrix of DIPOLE_ASSB_LEFT_IRR_MULTI) for NRHS
!     right sides into MXPOT_MULTI_COMB(1:NTOT, MDFIRST:MDFIRST+NRHS-1, 1): radiation columns 6(b-1)+k, diffraction
!     column 6 NBODY + 1.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_SOLVE_MULTI(AMAT,IPIV,B,NTOT,NRHS,MDFIRST)
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
      MXPOT_MULTI_COMB(1:NTOT,MDFIRST:MDFIRST+NRHS-1,1)=X(:,:)
      DEALLOCATE(X)

      RETURN
      END SUBROUTINE DIPOLE_SOLVE_MULTI

!-------------------------------------------------------------------------------
!     The sheet's share of the force sums on body IB: F(k) = sum over IB's sheet panels m of (-mu_m) nu_Dk,m dS_m, for the
!     potential in MXPOT_MULTI_COMB(:, MD, 1). Zero without dipole panels on IB.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_FORCE_SUM_MULTI(IB,MD,F)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: IB,MD
      COMPLEX*16,INTENT(OUT):: F(6)
      INTEGER:: M,NH

      F=CMPLX(0.D0,0.D0)
      IF (DNELEM_TOTAL.LE.0) RETURN
      IF (DNELEM_MULTI(IB).LE.0) RETURN
      NH=DIPOLE_NH_MULTI()
      DO M=DSTART(IB),DEND(IB)
       F(:)=F(:)-MXPOT_MULTI_COMB(NH+M,MD,1)*DDS_COMB(M)*DNRM_COMB(M,1:6)
      ENDDO

      RETURN
      END SUBROUTINE DIPOLE_FORCE_SUM_MULTI

!-------------------------------------------------------------------------------
!     IRSP = 1 with dipole panels: the NTOT + NL rows (hull rows of all bodies, dipole rows, then the lid rows of all
!     bodies) into AFULL, and the normal matrix C = A^H A (ZGEMM) with its LU factors. The lid rows are those of
!     ASSB_LEFT_IRR_MULTI (BODINT_LEFT_IRR_MULTI, IRR = 3) plus -W_lj for every sheet panel. Needs CALGREEN_MULTI_IRR.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ASSB_LEFT_IRR_MULTI(AFULL,CMAT,IPIV,NTOT,NL)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NTOT,NL
      COMPLEX*16,INTENT(OUT):: AFULL(NTOT+NL,NTOT),CMAT(NTOT,NTOT)
      INTEGER,INTENT(OUT):: IPIV(NTOT)
      INTEGER:: L,JEL,FLAG,INFO,NH
      REAL*8:: DIST,XF(3)
      COMPLEX*16:: TINDP(4)

      NH=DIPOLE_NH_MULTI()
      CALL FILL_SQUARE(AFULL,NTOT+NL,NTOT,.TRUE.)
!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(L,JEL,FLAG,DIST,TINDP,XF) SCHEDULE(DYNAMIC,16)
      DO L=1,NL
       XF(:)=iXYZ_GLOBAL_MULTI_COMB_P(L,:)
       DO JEL=1,NH
        DIST=SQRT((XF(1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2+(XF(2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2                  &
                 +(XF(3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
        TINDP=CMPLX(0.D0,0.D0)
        CALL BODINT_LEFT_IRR_MULTI(1,L,JEL,TINDP,3,FLAG)
        AFULL(NTOT+L,JEL)=RXY_MULTI(1,1)*TINDP(1)
       ENDDO
       DO JEL=1,DNELEM_TOTAL
        AFULL(NTOT+L,NH+JEL)=-SHEET_COLUMN(XF,JEL)
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL ZGEMM('C','N',NTOT,NTOT,NTOT+NL,(1.D0,0.D0),AFULL,NTOT+NL,AFULL,NTOT+NL,(0.D0,0.D0),CMAT,NTOT)
      CALL ZGETRF(NTOT,NTOT,CMAT,NTOT,IPIV,INFO)
      IF (INFO.NE.0) THEN
       PRINT *,' Error: the normal matrix with dipole panels of several bodies is singular, ZGETRF INFO =',INFO
       STOP 1
      ENDIF

      RETURN
      END SUBROUTINE DIPOLE_ASSB_LEFT_IRR_MULTI

!-------------------------------------------------------------------------------
!     IRSP = 1 right sides B(1:NTOT+NL, NRHS): radiation of body IB (KIND = 1, NRHS = 6) or diffraction (KIND = 2,
!     NRHS = 1, IB not used). Hull rows as RBC_IRR_MULTI / DBC_IRR_MULTI (IRR = 1), dipole rows as DIPOLE_ROWS_RBC /
!     _DBC, lid rows as RBC_IRR_MULTI / DBC_IRR_MULTI (IRR = 3); then D = A^H B (ZGEMM).
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_ASSB_RHS_IRR_MULTI(AFULL,B,D,NTOT,NL,NRHS,KIND,IB)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NTOT,NL,NRHS,KIND,IB
      COMPLEX*16,INTENT(IN):: AFULL(NTOT+NL,NTOT)
      COMPLEX*16,INTENT(OUT):: B(NTOT+NL,NRHS),D(NTOT,NRHS)
      INTEGER:: IEL,JEL,FLAG,IRR,IROW,NH,J0,J1
      REAL*8:: DIST,XF(3)
      COMPLEX*16:: TINRD(4,6,4),TINRDD(4,4),BTMP(6)
      COMPLEX*16,ALLOCATABLE:: BSQ(:,:)

      NH=DIPOLE_NH_MULTI()
      IF (KIND.EQ.1) THEN
       CALL HULL_RANGE(IB,J0,J1)
      ELSE
       J0=1
       J1=NH
      ENDIF
      B=CMPLX(0.D0,0.D0)
!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(IROW,IEL,JEL,FLAG,IRR,DIST,XF,TINRD,TINRDD,BTMP) SCHEDULE(DYNAMIC,16)
      DO IROW=1,NH+NL
       IF (IROW.LE.NH) THEN
        IEL=IROW
        IRR=1
        XF(:)=XYZ_GLOBAL_MULTI_COMB_P(IEL,:)
       ELSE
        IEL=IROW-NH
        IRR=3
        XF(:)=iXYZ_GLOBAL_MULTI_COMB_P(IEL,:)
       ENDIF
       BTMP=CMPLX(0.D0,0.D0)
       DO JEL=J0,J1
        DIST=SQRT((XF(1)-XYZ_GLOBAL_MULTI_COMB_P(JEL,1))**2+(XF(2)-XYZ_GLOBAL_MULTI_COMB_P(JEL,2))**2                  &
                 +(XF(3)-XYZ_GLOBAL_MULTI_COMB_P(JEL,3))**2)
        IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
         FLAG=1
        ELSE
         FLAG=0
        ENDIF
        IF (KIND.EQ.1) THEN
         TINRD=CMPLX(0.D0,0.D0)
         CALL RBC_IRR_MULTI(1,IEL,JEL,TINRD,IRR,FLAG)
         BTMP(:)=BTMP(:)+TINRD(1,:,1)
        ELSE
         TINRDD=CMPLX(0.D0,0.D0)
         CALL DBC_IRR_MULTI(1,IEL,JEL,TINRDD,IRR,FLAG)
         BTMP(1)=BTMP(1)+TINRDD(1,1)
        ENDIF
       ENDDO
       IF (IROW.LE.NH) THEN
        B(IEL,1:NRHS)=RXY_MULTI(1,1)*BTMP(1:NRHS)
       ELSE
        B(NTOT+IEL,1:NRHS)=RXY_MULTI(1,1)*BTMP(1:NRHS)
       ENDIF
      ENDDO
!$OMP END PARALLEL DO
      ALLOCATE(BSQ(NTOT,NRHS))
      BSQ=CMPLX(0.D0,0.D0)
      IF (KIND.EQ.1) THEN
       CALL DIPOLE_ROWS_RBC(BSQ,NH,IB)
      ELSE
       CALL DIPOLE_ROWS_DBC(BSQ(:,1),NH)
      ENDIF
      B(NH+1:NTOT,:)=BSQ(NH+1:NTOT,:)
      DEALLOCATE(BSQ)

      CALL ZGEMM('C','N',NTOT,NRHS,NTOT+NL,(1.D0,0.D0),AFULL,NTOT+NL,B,NTOT+NL,(0.D0,0.D0),D,NTOT)

      RETURN
      END SUBROUTINE DIPOLE_ASSB_RHS_IRR_MULTI

      SUBROUTINE DIPOLE_RELEASE_MULTI
      IF (ALLOCATED(FDIPM)) DEALLOCATE(FDIPM)
      END SUBROUTINE DIPOLE_RELEASE_MULTI
!-------------------------------------------------------------------------------
END MODULE DipoleMulti
!*******************************************************************************
