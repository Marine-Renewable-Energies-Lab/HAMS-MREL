!
!  ------------------------------------------------------------------------------------------------------
!
!    Damping lid on the free surface for several bodies: the solve on the combined hull arrays. The lid input, its geometry,
!    its checks and its panel integrals are those of the single body, in DampingLid.f90, since the lid belongs to no body.
!
!  License:
!
!    This routine is part of HAMS-MREL.
!
!    HAMS-MREL is a free software framework: you can redistribute it and/or modify it
!    under the terms of the Apache License, Version 2.0 (the "License"); you may
!    not use this subroutine except in compliance with the License. The software is
!    distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND.
!
!    You should have received a copy of the Apache License, Version 2.0, along with
!    HAMS. If not, see <http://www.apache.org/licenses/LICENSE-2.0>.
!
!  ---------------------------------------------------------------------------------------------------------
!
!  DampingLidMulti must not USE PressureElevation or PressureElevationMulti, and nothing in SingularIntgr* or AssbMatx* may
!  USE it.
!
MODULE DampingLidMulti

   USE HAMS_mod,       ONLY: NTHREAD
   USE Const_mod
   USE PanelMesh_mod
   USE Lid_mod
   USE SingularIntgrMulti
   USE WaveDyn_mod,     ONLY: V,H,WVN,NK,BETA,NBODY
   USE Body_mod,        ONLY: XW
   USE Potentials_mod,  ONLY: MXPOT_MULTI_COMB
   USE INFG3D_Open,     ONLY: INFGREEN3D
   USE FinGreen3D_Open, ONLY: FINGREEN3D
   USE PatcVelct,       ONLY: VINP,DINP
   USE AssbMatxMulti,   ONLY: FILL_HULL_LEFT_MULTI,ASSB_RBC_MULTI,ASSB_DBC_MULTI
   USE DampingLid,      ONLY: LID_SINGLE_LAYER,SYM_SIGNS,REFLECTED

   IMPLICIT NONE

   PRIVATE

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: LID_ASSEMBLE_MULTI
   PUBLIC :: LID_RADIATION_MULTI
   PUBLIC :: LID_DIFFRACTION_MULTI

   ! ..... The multi-body system with a lid, built by LID_ASSEMBLE_MULTI once per frequency .....................................
   !   Unknowns [phi(1:NH), q(1:LNELEM)], NH the hull panels of all bodies in the order of the combined arrays; rows: hull, then
   !   lid. Square, so LU-factorised. CHECK_LID stops a run with IRSP = 1, which would add interior rows.
   COMPLEX*16,ALLOCATABLE:: ALIDM(:,:,:)                            ! LU factors
   INTEGER,ALLOCATABLE:: IPIVM(:,:)
   COMPLEX*16,ALLOCATABLE:: SROWM(:,:,:)                            ! single layer of hull panel JEL, symmetry image IS, at lid centroid K

CONTAINS
!=======================================================================
!SUBROUTINE HULL_LAYERS_MULTI

!   !   HULL_LAYERS of DampingLid.f90 on the combined arrays: the single layer SH and the double layer DH of hull panel JEL at
!   !   field point XF, already reflected. Exact 1/r panel integrals and IMAGE_TERM_EXACT_MULTI within 50 panel sizes, one point
!   !   beyond, with the wave part of the Green function at the centroid. DH is that of BODINT_MULTI_LEFT and SH that of
!   !   RBC_MULTI_RIGHT. Full precision: CGRN_MULTI_COMB and RKBN_MULTI_COMB are not used.

      SUBROUTINE HULL_LAYERS_MULTI(XF,JEL,SH,DH)
      IMPLICIT   NONE

      REAL*8,INTENT(IN):: XF(3)
      INTEGER,INTENT(IN):: JEL
      COMPLEX*16,INTENT(OUT):: SH,DH
      INTEGER:: FLAG
      REAL*8:: XQ(3),XC(3),RK(4),DIST,EN(3)
      COMPLEX*16:: GRN(4)

      XQ=XYZ_GLOBAL_MULTI_COMB_P(JEL,:)
      EN=GEONRM_MULTI_COMB(JEL,:)
      DIST=SQRT(SUM((XF-XQ)**2))
      IF (DIST.LE.50.D0*PNSZ_MULTI_COMB(JEL)) THEN
       FLAG=1
      ELSE
       FLAG=0
      ENDIF

      IF (H.LT.0.D0) THEN
       CALL INFGREEN3D(XQ(1),XF(1),XQ(2),XF(2),XQ(3),XF(3),V,GRN,FLAG)
      ELSE
       CALL FINGREEN3D(XQ(1),XF(1),XQ(2),XF(2),XQ(3),XF(3),V,WVN,NK,H,GRN,FLAG)
      ENDIF

      IF (FLAG.EQ.1) THEN
       XC=XF
       IF (NCN_MULTI_COMB(JEL).EQ.3) THEN
        CALL SGLINTBD_MULTI_TRI2(XC,JEL,RK,1)
       ELSE
        CALL SGLINTBD_MULTI_QUAD2(XC,JEL,RK,1)
       ENDIF
       CALL IMAGE_TERM_EXACT_MULTI(XC,XQ,V,JEL,1,GRN,RK(1),RK(2:4))
       SH=RK(1)+GRN(1)*DS_MULTI_COMB(JEL)
       DH=(RK(2)+GRN(2)*DS_MULTI_COMB(JEL))*EN(1)+(RK(3)+GRN(3)*DS_MULTI_COMB(JEL))*EN(2)+(RK(4)+GRN(4)*DS_MULTI_COMB(JEL))*EN(3)
      ELSE
       SH=GRN(1)*DS_MULTI_COMB(JEL)
       DH=(GRN(2)*EN(1)+GRN(3)*EN(2)+GRN(4)*EN(3))*DS_MULTI_COMB(JEL)
      ENDIF

      RETURN
      END SUBROUTINE HULL_LAYERS_MULTI

!=======================================================================
!SUBROUTINE LID_ASSEMBLE_MULTI

!   !   The left-hand side with a lid for several bodies, for each symmetry class IP, once per frequency. The equations are those
!   !   of LID_ASSEMBLE in DampingLid.f90, with the hull of every body in the combined arrays:
!   !   - hull: FILL_HULL_LEFT_MULTI, and -S of the lid panels at the hull centroids.
!   !   - lid: phi = 0, so no free term; D of the hull panels and -S of the lid panels at the lid centroids.
!   !   Symmetry images enter as in ASSB_LEFT_MULTI. Needs CALGREEN_MULTI.

      SUBROUTINE LID_ASSEMBLE_MULTI(NH,NSYS)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NH,NSYS
      INTEGER:: NC,IEL,JEL,K,L,IS,IP,INFO
      REAL*8:: XF(3)
      COMPLEX*16:: SL,SH,DH
      COMPLEX*16,ALLOCATABLE:: AH(:,:,:)

      NC=NH+LNELEM
      IF (.NOT.ALLOCATED(ALIDM)) ALLOCATE(ALIDM(NC,NC,NSYS),IPIVM(NC,NSYS),SROWM(LNELEM,NH,NSYS))

      ALIDM=(0.D0,0.D0)
      ALLOCATE(AH(NH,NH,NSYS))
      CALL FILL_HULL_LEFT_MULTI(AH,NH,NSYS)
      ALIDM(1:NH,1:NH,:)=AH
      DEALLOCATE(AH)

!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(IEL,IS,IP,L,XF,SL) SCHEDULE(DYNAMIC,16)
      DO IEL=1,NH
       DO IS=1,NSYS
        XF=REFLECTED(XYZ_GLOBAL_MULTI_COMB_P(IEL,:),IS)
        DO L=1,LNELEM
         CALL LID_SINGLE_LAYER(XF,L,SL)
         DO IP=1,NSYS
          ALIDM(IEL,NH+L,IP)=ALIDM(IEL,NH+L,IP)-RXY_MULTI(IS,IP)*SL
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(K,IS,IP,JEL,L,XF,SL,SH,DH) SCHEDULE(DYNAMIC,16)
      DO K=1,LNELEM
       DO IS=1,NSYS
        XF=REFLECTED(LXYZ_P(K,:),IS)
        DO JEL=1,NH
         CALL HULL_LAYERS_MULTI(XF,JEL,SH,DH)
         SROWM(K,JEL,IS)=SH
         DO IP=1,NSYS
          ALIDM(NH+K,JEL,IP)=ALIDM(NH+K,JEL,IP)+RXY_MULTI(IS,IP)*DH
         ENDDO
        ENDDO
        DO L=1,LNELEM
         CALL LID_SINGLE_LAYER(XF,L,SL)
         DO IP=1,NSYS
          ALIDM(NH+K,NH+L,IP)=ALIDM(NH+K,NH+L,IP)-RXY_MULTI(IS,IP)*SL
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      DO IP=1,NSYS
       CALL ZGETRF(NC,NC,ALIDM(:,:,IP),NC,IPIVM(:,IP),INFO)
       IF (INFO.NE.0) THEN
        PRINT *,' Error: the system with the damping lid is singular, ZGETRF INFO =',INFO
        STOP 1
       ENDIF
      ENDDO

      RETURN
      END SUBROUTINE LID_ASSEMBLE_MULTI

!=======================================================================
!SUBROUTINE LID_SOLVE_UNPACK_MULTI

!   !   Solves for the NRHS columns of B(NH+LNELEM,NRHS,NSYS) by the LU factors of LID_ASSEMBLE_MULTI, then puts the hull phi of
!   !   each symmetry image into MXPOT_MULTI_COMB(:,MD0+1:MD0+NRHS,:) and the lid q into MXQ, by the RXY_MULTI recombination and
!   !   the /NSYS of RADIATION_SOLVER_MULTI.

      SUBROUTINE LID_SOLVE_UNPACK_MULTI(B,NRHS,MD0,NH,NSYS)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NRHS,MD0,NH,NSYS
      COMPLEX*16,INTENT(INOUT):: B(NH+LNELEM,NRHS,NSYS)
      INTEGER:: NC,M,IP,IS,INFO

      NC=NH+LNELEM
      DO IP=1,NSYS
       CALL ZGETRS('N',NC,NRHS,ALIDM(:,:,IP),NC,IPIVM(:,IP),B(:,:,IP),NC,INFO)
      ENDDO

      DO M=1,NRHS
       DO IP=1,NSYS
        MXPOT_MULTI_COMB(:,MD0+M,IP)=(0.D0,0.D0)
        MXQ(:,MD0+M,IP)=(0.D0,0.D0)
        DO IS=1,NSYS
         MXPOT_MULTI_COMB(:,MD0+M,IP)=MXPOT_MULTI_COMB(:,MD0+M,IP)+B(1:NH,M,IS)*RXY_MULTI(IP,IS)
         MXQ(:,MD0+M,IP)=MXQ(:,MD0+M,IP)+B(NH+1:NC,M,IS)*RXY_MULTI(IP,IS)
        ENDDO
        MXPOT_MULTI_COMB(:,MD0+M,IP)=MXPOT_MULTI_COMB(:,MD0+M,IP)/NSYS
        MXQ(:,MD0+M,IP)=MXQ(:,MD0+M,IP)/NSYS
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE LID_SOLVE_UNPACK_MULTI

!=======================================================================
!SUBROUTINE LID_RADIATION_MULTI

!   !   The 6*NBODY radiation problems with the lid. Mode MD of body M moves the panels of body M only. Hull rows from
!   !   ASSB_RBC_MULTI over the panels of each body, as in the multi-body loop of HAMS_Prog.f90. Lid rows: sum_j S_xj V_j over
!   !   the panels j of that body, with S from LID_ASSEMBLE_MULTI and the mode signs of RBC_MULTI_RIGHT; on the lid phi = 0, so
!   !   nothing else. Fills MXPOT_MULTI_COMB(:,1:6*NBODY,:) and MXQ(:,1:6*NBODY,:).

      SUBROUTINE LID_RADIATION_MULTI(NH,NSYS)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NH,NSYS
      INTEGER:: K,JEL,IS,IP,MD,M,J0,J1
      INTEGER:: JFIRST(NBODY),JLAST(NBODY)
      REAL*8:: SGX,SGY,MS(6)
      COMPLEX*16:: BTMP(6,2)
      COMPLEX*16,ALLOCATABLE:: BH(:,:,:),B(:,:,:)

      ! The panels of each body in the combined arrays
      J1=0
      DO M=1,NBODY
       JFIRST(M)=J1+1
       J1=J1+NELEM_MULTI(M)
       JLAST(M)=J1
      ENDDO

      ALLOCATE(BH(NH,6*NBODY,NSYS),B(NH+LNELEM,6*NBODY,NSYS))
      BH=(0.D0,0.D0)
      DO M=1,NBODY
       CALL ASSB_RBC_MULTI(BH(:,6*(M-1)+1:6*M,:),NH,NSYS,JFIRST(M),JLAST(M))
      ENDDO
      B=(0.D0,0.D0)
      B(1:NH,:,:)=BH

!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(K,M,J0,J1,JEL,IS,IP,MD,SGX,SGY,MS,BTMP) SCHEDULE(DYNAMIC,4)
      DO K=1,LNELEM
       DO M=1,NBODY
        J0=JFIRST(M)
        J1=JLAST(M)
        BTMP=(0.D0,0.D0)
        DO IS=1,NSYS
         DO IP=1,NSYS
          CALL SYM_SIGNS(IS,IP,SGX,SGY)
          MS=(/SGX,SGY,1.D0,SGY,SGX,SGX*SGY/)
          DO JEL=J0,J1
           BTMP(:,IP)=BTMP(:,IP)+SROWM(K,JEL,IS)*MS*DXYZ_MULTI_COMB(JEL,:)
          ENDDO
         ENDDO
        ENDDO
        DO MD=1,6
         DO IP=1,NSYS
          DO IS=1,NSYS
           B(NH+K,6*(M-1)+MD,IP)=B(NH+K,6*(M-1)+MD,IP)+RXY_MULTI(IP,IS)*BTMP(MD,IS)
          ENDDO
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL LID_SOLVE_UNPACK_MULTI(B,6*NBODY,0,NH,NSYS)
      DEALLOCATE(BH,B)

      RETURN
      END SUBROUTINE LID_RADIATION_MULTI

!=======================================================================
!SUBROUTINE LID_DIFFRACTION_MULTI

!   !   The diffraction problem with the lid for several bodies, for the heading BETA, as LID_DIFFRACTION in DampingLid.f90:
!   !   hull rows from ASSB_DBC_MULTI; lid rows -sum_j S_xj dphi_I/dn_j over the panels of all bodies, plus 4 pi phi_I(x), since
!   !   the scattered potential is -phi_I on the lid. At zero frequency the incident term is left out. Fills
!   !   MXPOT_MULTI_COMB(:,6*NBODY+1,:) and MXQ(:,6*NBODY+1,:).

      SUBROUTINE LID_DIFFRACTION_MULTI(NH,NSYS)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NH,NSYS
      INTEGER:: K,JEL,IS,IP
      REAL*8:: SGX,SGY,XF(3)
      COMPLEX*16:: BTMP(2),DPOX,DPOY,DPOZ
      COMPLEX*16,ALLOCATABLE:: BH(:,:),B(:,:,:),DINC(:,:,:)

      ALLOCATE(BH(NH,NSYS),B(NH+LNELEM,1,NSYS),DINC(NH,NSYS,NSYS))
      CALL ASSB_DBC_MULTI(BH,NH,NSYS)
      B=(0.D0,0.D0)
      B(1:NH,1,:)=BH

      ! dphi_I/dn at hull panel JEL reflected between images IS and IP, as DBC_MULTI_RIGHT
      DO IS=1,NSYS
       DO IP=1,NSYS
        CALL SYM_SIGNS(IS,IP,SGX,SGY)
        DO JEL=1,NH
         CALL DINP(SGX*XYZ_GLOBAL_MULTI_COMB_P(JEL,1),SGY*XYZ_GLOBAL_MULTI_COMB_P(JEL,2),XYZ_GLOBAL_MULTI_COMB_P(JEL,3), &
                   XW(1),XW(2),BETA,DPOX,DPOY,DPOZ)
         DINC(JEL,IS,IP)=SGX*DPOX*GEONRM_MULTI_COMB(JEL,1)+SGY*DPOY*GEONRM_MULTI_COMB(JEL,2)+DPOZ*GEONRM_MULTI_COMB(JEL,3)
        ENDDO
       ENDDO
      ENDDO

!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(K,JEL,IS,IP,XF,BTMP) SCHEDULE(DYNAMIC,16)
      DO K=1,LNELEM
       BTMP=(0.D0,0.D0)
       DO IS=1,NSYS
        DO IP=1,NSYS
         DO JEL=1,NH
          BTMP(IP)=BTMP(IP)-SROWM(K,JEL,IS)*DINC(JEL,IS,IP)
         ENDDO
        ENDDO
       ENDDO
       IF (ABS(V).GE.1.D-8) THEN
        DO IS=1,NSYS
         XF=REFLECTED(LXYZ_P(K,:),IS)
         BTMP(IS)=BTMP(IS)+4.D0*PI*VINP(XF(1),XF(2),XF(3),XW(1),XW(2),BETA)
        ENDDO
       ENDIF
       DO IP=1,NSYS
        DO IS=1,NSYS
         B(NH+K,1,IP)=B(NH+K,1,IP)+RXY_MULTI(IP,IS)*BTMP(IS)
        ENDDO
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL LID_SOLVE_UNPACK_MULTI(B,1,6*NBODY,NH,NSYS)
      DEALLOCATE(BH,B,DINC)

      RETURN
      END SUBROUTINE LID_DIFFRACTION_MULTI

END MODULE DampingLidMulti
