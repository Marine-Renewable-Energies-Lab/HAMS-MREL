!
!  ------------------------------------------------------------------------------------------------------
!
!    Damping lid on the free surface: input, geometry and the checks that stop a run on a bad lid.
!    The lid is a list of panels on z = 0, in global coordinates, where the potential is zero (docs/input_files.md).
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
!  DampingLid must not USE PressureElevation, and nothing in SingularIntgr* or AssbMatx* may USE DampingLid.
!
MODULE DampingLid

   USE HAMS_mod,       ONLY: ISOLV
   USE Const_mod
   USE PanelMesh_mod
   USE Lid_mod
   USE NormalProcess
   USE SingularIntgr
   USE SingularIntgrMulti

   IMPLICIT NONE

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: READ_LID_MESH
   PUBLIC :: LID_GEOMETRY
   PUBLIC :: CHECK_LID
   PUBLIC :: ALLOCATE_LID_UNKNOWNS
   PUBLIC :: HULL_SOLID_ANGLE

   INTEGER,PRIVATE:: NBODY_LID=1                                    ! set by CHECK_LID: selects the single or multi-body hull arrays

CONTAINS
!=======================================================================
!SUBROUTINE READ_LID_MESH

!   !   Reads DampingLid.pnl, in the format of WaterPlaneMesh.pnl. Unlike the waterplane read in HAMS_Prog, the symmetry
!   !   flags of the file are compared with the hull's and never stored over them.

      SUBROUTINE READ_LID_MESH(INPUTDIR)
      IMPLICIT   NONE

      CHARACTER(LEN=*),INTENT(IN):: INPUTDIR
      INTEGER:: IU,IND,IEL,J,M,LISX,LISY

      OPEN(NEWUNIT=IU,FILE=TRIM(INPUTDIR)//"/DampingLid.pnl",STATUS="OLD",ACTION="READ")

      DO J=1,3
       READ(IU,*)
      ENDDO
      READ(IU,*) LNELEM,LNTND,LISX,LISY
      IF (LISX.NE.ISX.OR.LISY.NE.ISY) THEN
       PRINT *,' Error: lid symmetry differs from the hull.'
       PRINT *,' ISX, ISY of DampingLid.pnl:',LISX,LISY,'  of the hull:',ISX,ISY
       STOP 1
      ENDIF

      ALLOCATE(LXYZ(LNTND,3),LNCN(LNELEM),LNCON(LNELEM,4))
      LNCON=0

      DO J=1,2
       READ(IU,*)
      ENDDO
      DO IND=1,LNTND
       READ(IU,*) M,LXYZ(IND,1),LXYZ(IND,2),LXYZ(IND,3)
       IF (ABS(LXYZ(IND,3)).GT.1.D-10) THEN
        PRINT *,' Error: lid node is not on z = 0, node No.',IND
        STOP 1
       ENDIF
      ENDDO

      DO J=1,3
       READ(IU,*)
      ENDDO
      DO IEL=1,LNELEM
       READ(IU,*) M,LNCN(IEL),(LNCON(IEL,J),J=1,LNCN(IEL))
      ENDDO

      CLOSE(IU)

      RETURN
      END SUBROUTINE READ_LID_MESH

!=======================================================================
!SUBROUTINE LID_GEOMETRY

!   !   Centroids, areas and chart lengths of the lid panels. No normal: the lid has none, so the vertex order of a panel
!   !   does not matter.

      SUBROUTINE LID_GEOMETRY()
      IMPLICIT   NONE

      ALLOCATE(LXYZ_P(LNELEM,3),LDS(LNELEM),LPNSZ(LNELEM))

      CALL CalPanelCentre(LXYZ,LNTND,LNELEM,LNCN,LNCON,LXYZ_P)
      CALL CalPanelArea(LXYZ,LNTND,LNELEM,LNCN,LNCON,LDS)
      CALL CalPanelChartLength(LXYZ,LNTND,LNELEM,LNCN,LNCON,LPNSZ)

      RETURN
      END SUBROUTINE LID_GEOMETRY

!=======================================================================
!SUBROUTINE ALLOCATE_LID_UNKNOWNS

!   !   q = dphi/dz - nu*phi on the lid: the radiation modes of every body and the diffraction problem, for each symmetry image.

      SUBROUTINE ALLOCATE_LID_UNKNOWNS(NBODYIN)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NBODYIN

      IF (NBODYIN.EQ.1) THEN
       ALLOCATE(MXQ(LNELEM,7,NSYS))
      ELSE
       ALLOCATE(MXQ(LNELEM,6*NBODYIN+1,NSYS))
      ENDIF
      MXQ=(0.D0,0.D0)

      RETURN
      END SUBROUTINE ALLOCATE_LID_UNKNOWNS

!=======================================================================
!FUNCTION HULL_SOLID_ANGLE

!   !   The sum over the hull panels, over the symmetry images of the hull, of the double-layer Rankine integral of 1/r,
!   !   n.grad(1/r) over the panel, at a field point XF on z = 0, plus its image in z = 0. The image of the hull is
!   !   accounted for by evaluating at XF with z mirrored, and at z = 0 that equals the direct term, so the total is twice it.
!   !   Over a closed surface this is a solid angle: about 4 pi in magnitude for XF inside and 0 outside. The hull is open
!   !   at z = 0 and its image closes it. The sign follows the panel normals. The field point is reflected for the symmetry
!   !   images as CALGREEN does.

      REAL*8 FUNCTION HULL_SOLID_ANGLE(XF)
      IMPLICIT   NONE

      REAL*8,INTENT(IN):: XF(3)
      INTEGER:: IS,JEL,IMG,NH
      REAL*8:: XP(3),XQ(3),RK(4),SUMA

      SUMA=0.D0
      IF (NBODY_LID.EQ.1) THEN
       NH=NELEM
      ELSE
       NH=SUM(NELEM_MULTI)
      ENDIF

      DO IS=1,NSYS
       IF (NBODY_LID.EQ.1) THEN
        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XP(1)=SY(IS,1)*XF(1)
         XP(2)=SX(IS,1)*XF(2)
        ELSE
         XP(1)=SX(IS,1)*XF(1)
         XP(2)=SY(IS,1)*XF(2)
        ENDIF
       ELSE
        IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
         XP(1)=SY_MULTI(IS,1)*XF(1)
         XP(2)=SX_MULTI(IS,1)*XF(2)
        ELSE
         XP(1)=SX_MULTI(IS,1)*XF(1)
         XP(2)=SY_MULTI(IS,1)*XF(2)
        ENDIF
       ENDIF
       XP(3)=XF(3)

       DO IMG=1,2
        XQ=XP
        IF (IMG.EQ.2) XQ(3)=-XP(3)                                  ! the image of the hull in z = 0
        DO JEL=1,NH
         IF (NBODY_LID.EQ.1) THEN
          IF (NCN(JEL).EQ.3) THEN
           CALL SGLINTBD_TRI2(XQ,JEL,RK,1)
          ELSE
           CALL SGLINTBD_QUAD2(XQ,JEL,RK,1)
          ENDIF
          SUMA=SUMA+RK(2)*GEONRM_P(JEL,1)+RK(3)*GEONRM_P(JEL,2)+RK(4)*GEONRM_P(JEL,3)
         ELSE
          IF (NCN_MULTI_COMB(JEL).EQ.3) THEN
           CALL SGLINTBD_MULTI_TRI2(XQ,JEL,RK,1)
          ELSE
           CALL SGLINTBD_MULTI_QUAD2(XQ,JEL,RK,1)
          ENDIF
          SUMA=SUMA+RK(2)*GEONRM_MULTI_COMB(JEL,1)+RK(3)*GEONRM_MULTI_COMB(JEL,2)+RK(4)*GEONRM_MULTI_COMB(JEL,3)
         ENDIF
        ENDDO
       ENDDO
      ENDDO

      HULL_SOLID_ANGLE=SUMA

      RETURN
      END FUNCTION HULL_SOLID_ANGLE

!=======================================================================
!SUBROUTINE CHECK_LID

!   !   The stops that need the hull: the solver type, and a lid panel inside a hull. The file and node checks are made
!   !   in the control file read and in READ_LID_MESH. Call after the hull normals and LID_GEOMETRY.

      SUBROUTINE CHECK_LID(NBODYIN)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NBODYIN
      INTEGER:: IEL,IFIRST
      REAL*8:: SANG,SFIRST

      NBODY_LID=NBODYIN

      IF (ISOLV.EQ.2) THEN
       PRINT *,' Error: damping lid needs Solver_type 1 (direct LU). The iterative solver has no lid rows.'
       STOP 1
      ENDIF

      IFIRST=LNELEM+1
      SFIRST=0.D0
!$OMP PARALLEL DO PRIVATE(IEL,SANG) REDUCTION(MIN:IFIRST) SCHEDULE(DYNAMIC,16)
      DO IEL=1,LNELEM
       SANG=HULL_SOLID_ANGLE(LXYZ_P(IEL,:))
       IF (ABS(SANG).GT.2.D0*PI) IFIRST=MIN(IFIRST,IEL)
      ENDDO
!$OMP END PARALLEL DO

      IF (IFIRST.LE.LNELEM) THEN
       SFIRST=HULL_SOLID_ANGLE(LXYZ_P(IFIRST,:))
       PRINT *,' Error: lid panel inside a hull, panel No.',IFIRST
       PRINT *,' The solid angle of the hull at its centroid is',SFIRST,' (4 pi inside a hull, 0 outside).'
       PRINT *,' The lid belongs on open water. Use the water plane mesh with If_remove_irr_freq for an interior free surface.'
       STOP 1
      ENDIF

      RETURN
      END SUBROUTINE CHECK_LID

END MODULE DampingLid
