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
MODULE ReadPanelMesh

   USE Body_mod
   USE PanelMesh_mod
   USE Inerfs_mod
   USE NormalProcess
   USE DipoleMesh_mod
   
   IMPLICIT NONE

      ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: ReadBodyMesh
   PUBLIC :: ReadWTPLMesh
   
   PUBLIC :: CalNormals
   PUBLIC :: DIPOLE_SETUP
   PUBLIC :: DIPOLE_CHECK_OPTIONS
   PUBLIC :: DIPOLE_CHECK_GEOMETRY
   PUBLIC :: DIPOLE_RELEASE
   PUBLIC :: DIPOLE_CHECK_CLEARANCE
   PUBLIC :: DIPOLE_CHECK_COPLANAR

CONTAINS
!=======================================================================
!SUBROUTINE ReadPanelMesh

!   !   ReadPanelMesh is used to read nodal and elemental data on the 
!   !   body surface of an arbitrary body

      SUBROUTINE ReadBodyMesh
      IMPLICIT   NONE

      INTEGER M,N,IND,IEL,J
      REAL*8::  DX,DY,DZ

        DO 10 IND=1,NTND                                                              ! The 10 is a label. Using this, statements to be continued if the loop does not progress can be mentioned with the CONTINUE command. Reference: http://www.personal.psu.edu/jhm/f90/lectures/15.html
          READ(2,*) M, XYZ(IND,1), XYZ(IND,2), XYZ(IND,3)
10      CONTINUE

         DO J=1,3
          READ(2,*)
         ENDDO
         
        DO 30  IEL=1, NELEM
	    READ(2, *) M, NCN(IEL), (NCON(IEL,J), J=1, NCN(IEL))
30      CONTINUE

      RETURN
      END SUBROUTINE ReadBodyMesh

!=======================================================================
!SUBROUTINE ReadWTPLMesh


!   !   ReadWTPLMesh is used to read nodal and elemental data on the                     
!   !   inner water plane

      SUBROUTINE ReadWTPLMesh
      IMPLICIT   NONE  

      INTEGER M,N,IND,IEL,J
      REAL*8::  DX,DY,DZ
!
! -------------------------------------------------------------------------
! 
      DO 10 IND=1,iNTND
        READ(5,*) M, iXYZ(IND,1), iXYZ(IND,2), iXYZ(IND,3)
        IF (ABS(iXYZ(IND,3)).GT.1.E-10) THEN
         Print *,' Error: Z Coordinate is not zero at Node No.',IND
         STOP
        ENDIF
10    CONTINUE

         DO J=1,3
          READ(5,*)
         ENDDO
         
      DO 30  IEL=1, iNELEM
	    READ(5, *) M, iNCN(IEL),(iNCON(IEL,J), J=1, iNCN(IEL))    
30    CONTINUE

       RETURN
       END SUBROUTINE ReadWTPLMesh

        
!=======================================================================
!SUBROUTINE CalNormals


!   !   CalNormals is used to calculate normals on the immersed body surface 
!   !   and on the inner water plane


      SUBROUTINE CalNormals(IFLAG,INPUTDIR)
      IMPLICIT   NONE

      INTEGER IEL,IND
      INTEGER,INTENT(IN):: IFLAG
      CHARACTER(LEN=*),INTENT(IN):: INPUTDIR

      ! Local variables for generalized modes
      INTEGER:: GEN_DOF, GEN_PID, GEN_IEL, GEN_IOST
      REAL*8:: GEN_NORMAL
      CHARACTER(LEN=512):: GEN_FILENAME
      CHARACTER(LEN=10):: DOF_STR
      LOGICAL:: GEN_FILE_EXISTS

! -------------------------------------------------------------------------
! 
!      Calculate some panel properties
!
      CALL CalPanelCentre( XYZ, NTND, NELEM, NCN, NCON, XYZ_P)

      IEL=FIRST_PANEL_ON_Z0(XYZ_P(:,3),NELEM)
      IF (IEL.GT.0) THEN
         PRINT *,' Error: hull panel',IEL,' has its centroid on the free surface, z = 0.'
         PRINT *,' The image term of the Green function is singular at it. Keep hull panels below z = 0,'
         PRINT *,' and put an interior free surface in the water plane mesh.'
         STOP 1
      ENDIF

      CALL CalPanelArea( XYZ, NTND, NELEM, NCN, NCON, DS)

      CALL CalPanelChartLength( XYZ, NTND, NELEM, NCN, NCON, PNSZ)

      CALL CalTransNormals( XYZ, NTND, NELEM, NCN, NCON, DXYZ_P)

      CALL CalRotNormals(XR, XYZ_P, DXYZ_P, NELEM)

      GEONRM_P(1:NELEM,1:3)=DXYZ_P(1:NELEM,1:3)                                       ! Before a generalized mode overwrites DXYZ_P

      IF (IFLAG.NE.0) THEN                                                            ! IFLAG is used as an indicator for irregular frequency removal. If 1, these will be removed

         CALL CalPanelCentre(iXYZ,INTND,INELEM,INCN,INCON,iXYZ_P)

         CALL CalPanelArea(iXYZ,INTND,INELEM,INCN,INCON,IDS)

         CALL CalPanelChartLength(iXYZ,INTND,INELEM,INCN,INCON,IPNSZ)

         CALL CalTransNormals(iXYZ,INTND,INELEM,INCN,INCON,IDXYZ_P)

         CALL CalRotNormals(XR,iXYZ_P,IDXYZ_P,INELEM)
      ENDIF
!
! -------------------------------------------------------------------------
!     Read generalized mode normals from gen_mod files if present
!     File naming convention: gen_mod_1_{dof}.txt  (body number is always 1)
!     File format per line:   panel_id    normal_value
! -------------------------------------------------------------------------

      DO GEN_DOF = 1, 6
        WRITE(DOF_STR,'(I0)') GEN_DOF
        GEN_FILENAME = TRIM(INPUTDIR)//'/gen_mod_1_'// &
                        TRIM(ADJUSTL(DOF_STR))//'.txt'

        INQUIRE(FILE=TRIM(GEN_FILENAME), EXIST=GEN_FILE_EXISTS)

        IF (GEN_FILE_EXISTS) THEN
          PRINT *, ' Found generalized mode file: ', TRIM(GEN_FILENAME)
          PRINT *, '   Overriding DOF ', GEN_DOF, ' normals for body 1'

          OPEN(99, FILE=TRIM(GEN_FILENAME), STATUS='OLD', ACTION='READ')

          DO GEN_IEL = 1, NELEM
            READ(99, *, IOSTAT=GEN_IOST) GEN_PID, GEN_NORMAL
            IF (GEN_IOST .NE. 0) THEN
              PRINT *, '  Warning: Could not read line ', GEN_IEL, &
                       ' from ', TRIM(GEN_FILENAME)
              EXIT
            ENDIF
            DXYZ_P(GEN_PID, GEN_DOF) = GEN_NORMAL
          ENDDO

          CLOSE(99)
          PRINT *, '   Loaded ', NELEM, ' generalized normals.'
        ENDIF
      ENDDO

      Print *,' Calculating panel normals is finished...'
      Print * 
!

      RETURN
      END SUBROUTINE CalNormals

!=======================================================================
!SUBROUTINE DIPOLE_SETUP
!   !   Reads the optional dipole (zero-thickness) panel mesh INPUTDIR/DipoleMesh.pnl, in the format of HullMesh.pnl,
!   !   and sets its geometry in DipoleMesh_mod (notes/specs/2026-11-24-dipole-panels-design.md, Sec. 6). Call it after
!   !   CalNormals. Without the file DNELEM stays 0, nothing is allocated and OK is true. With it, the v1 limits
!   !   (DIPOLE_CHECK_OPTIONS) and the geometry (DIPOLE_CHECK_GEOMETRY) are checked; OK false and MSG say why.
!   !   The panel normal follows the node order as for hull panels (CalTransNormals); rotations about XR.
      SUBROUTINE DIPOLE_SETUP(INPUTDIR,NSYS_IN,ISOLV_IN,NBODY_IN,SO_ACTIVE,OK,MSG)
      IMPLICIT   NONE
      CHARACTER(LEN=*),INTENT(IN):: INPUTDIR
      INTEGER,INTENT(IN):: NSYS_IN,ISOLV_IN,NBODY_IN
      LOGICAL,INTENT(IN):: SO_ACTIVE
      LOGICAL,INTENT(OUT):: OK
      CHARACTER(LEN=*),INTENT(OUT):: MSG
      CHARACTER(LEN=1024):: FNAME
      LOGICAL:: EXISTS
      INTEGER:: U,J,IND,IEL,M,KSX,KSY,IOS

      OK=.TRUE.
      MSG=''
      FNAME=TRIM(INPUTDIR)//'/DipoleMesh.pnl'
      INQUIRE(FILE=TRIM(FNAME),EXIST=EXISTS)
      IF (.NOT.EXISTS) RETURN

      CALL DIPOLE_CHECK_OPTIONS(NSYS_IN,ISOLV_IN,NBODY_IN,SO_ACTIVE,OK,MSG)
      IF (.NOT.OK) RETURN

      OPEN(NEWUNIT=U,FILE=TRIM(FNAME),STATUS='OLD',ACTION='READ',IOSTAT=IOS)
      IF (IOS.NE.0) THEN
       OK=.FALSE.
       MSG=' Error: cannot open '//TRIM(FNAME)
       RETURN
      ENDIF
      DO J=1,3
       READ(U,*)
      ENDDO
      READ(U,*,IOSTAT=IOS) DNELEM,DNTND,KSX,KSY
      IF (IOS.NE.0.OR.DNELEM.LE.0.OR.DNTND.LE.0) THEN
       OK=.FALSE.
       MSG=' Error: DipoleMesh.pnl: cannot read the numbers of panels and nodes.'
       DNELEM=0
       DNTND=0
       CLOSE(U)
       RETURN
      ENDIF
      IF (KSX.NE.0.OR.KSY.NE.0) THEN
       OK=.FALSE.
       MSG=' Error: dipole panels do not support symmetry: set both symmetry flags of DipoleMesh.pnl to 0.'
       DNELEM=0
       DNTND=0
       CLOSE(U)
       RETURN
      ENDIF
      DO J=1,2
       READ(U,*)
      ENDDO
      ALLOCATE(DPXYZ(DNTND,3),DNCN(DNELEM),DNCON(DNELEM,4))
      ALLOCATE(DPXYZ_C(DNELEM,3),DNRM(DNELEM,6),DDS(DNELEM),DPNSZ(DNELEM))
      DNCON=0
      DO IND=1,DNTND
       READ(U,*) M,DPXYZ(IND,1),DPXYZ(IND,2),DPXYZ(IND,3)
      ENDDO
      DO J=1,3
       READ(U,*)
      ENDDO
      DO IEL=1,DNELEM
       READ(U,*) M,DNCN(IEL),(DNCON(IEL,J),J=1,DNCN(IEL))
      ENDDO
      CLOSE(U)

      CALL CalPanelCentre(DPXYZ,DNTND,DNELEM,DNCN,DNCON,DPXYZ_C)
      CALL CalPanelArea(DPXYZ,DNTND,DNELEM,DNCN,DNCON,DDS)
      CALL CalPanelChartLength(DPXYZ,DNTND,DNELEM,DNCN,DNCON,DPNSZ)
      CALL CalTransNormals(DPXYZ,DNTND,DNELEM,DNCN,DNCON,DNRM)
      CALL CalRotNormals(XR,DPXYZ_C,DNRM,DNELEM)

      IF (ALLOCATED(XYZ_P)) THEN
       CALL DIPOLE_CHECK_GEOMETRY(DPXYZ,DNTND,DPXYZ_C,DPNSZ,DNELEM,XYZ_P,MIN(NELEM,SIZE(XYZ_P,1)),OK,MSG)
      ELSE
       CALL DIPOLE_CHECK_GEOMETRY(DPXYZ,DNTND,DPXYZ_C,DPNSZ,DNELEM,DPXYZ_C,0,OK,MSG)        ! no hull panels to compare
      ENDIF
      IF (.NOT.OK) RETURN
      CALL DIPOLE_CHECK_CLEARANCE(DPXYZ,DNCN,DNCON,DPNSZ,DNELEM,OK,MSG)
      IF (.NOT.OK) RETURN
      ! Lid panels (IRSP = 1) lie on z = 0, which the clearance check keeps every sheet away from.
      IF (ALLOCATED(XYZ_P).AND.ALLOCATED(GEONRM_P).AND.ALLOCATED(NCON)) THEN
       CALL DIPOLE_CHECK_COPLANAR(DPXYZ,DNCN,DNCON,DPXYZ_C,DNRM(:,1:3),DPNSZ,DNELEM,                               &
                                  XYZ,NCN,NCON,XYZ_P,GEONRM_P,PNSZ,MIN(NELEM,SIZE(XYZ_P,1)),OK,MSG)
       IF (.NOT.OK) RETURN
      ENDIF

      PRINT *,' Dipole panels read from DipoleMesh.pnl:',DNELEM
      RETURN
      END SUBROUTINE DIPOLE_SETUP

!=======================================================================
!SUBROUTINE DIPOLE_CHECK_OPTIONS
!   !   The v1 limits of dipole panels (plan 2026-11-24-dipole-panels-first-order, Global Constraints): one body, no
!   !   symmetry, direct LU, first order only. OK false and a one-line MSG naming the option otherwise.
      SUBROUTINE DIPOLE_CHECK_OPTIONS(NSYS_IN,ISOLV_IN,NBODY_IN,SO_ACTIVE,OK,MSG)
      IMPLICIT   NONE
      INTEGER,INTENT(IN):: NSYS_IN,ISOLV_IN,NBODY_IN
      LOGICAL,INTENT(IN):: SO_ACTIVE
      LOGICAL,INTENT(OUT):: OK
      CHARACTER(LEN=*),INTENT(OUT):: MSG

      OK=.FALSE.
      IF (NBODY_IN.NE.1) THEN
       MSG=' Error: dipole panels (DipoleMesh.pnl) support one body only.'
      ELSEIF (NSYS_IN.NE.1) THEN
       MSG=' Error: dipole panels (DipoleMesh.pnl) do not support symmetry: set the hull symmetry flags to 0.'
      ELSEIF (ISOLV_IN.NE.1) THEN
       MSG=' Error: dipole panels (DipoleMesh.pnl) need the direct solver, ISOLV = 1.'
      ELSEIF (SO_ACTIVE) THEN
       MSG=' Error: dipole panels (DipoleMesh.pnl) are first order only: remove SecondOrder.in or switch it off.'
      ELSE
       OK=.TRUE.
       MSG=''
      ENDIF

      RETURN
      END SUBROUTINE DIPOLE_CHECK_OPTIONS

!=======================================================================
!SUBROUTINE DIPOLE_CHECK_GEOMETRY
!   !   Every dipole node below z = 0, and no dipole centroid within 1e-6 panel sizes of a hull centroid (a sheet must
!   !   not coincide with hull panels). OK false and MSG naming the first offender otherwise.
      SUBROUTINE DIPOLE_CHECK_GEOMETRY(PXYZ,NND,PC,PSZ,ND,HC,NH,OK,MSG)
      IMPLICIT   NONE
      INTEGER,INTENT(IN):: NND,ND,NH
      REAL*8,INTENT(IN):: PXYZ(:,:),PC(:,:),PSZ(:),HC(:,:)              ! assumed shape: HC may have more rows than NH
      LOGICAL,INTENT(OUT):: OK
      CHARACTER(LEN=*),INTENT(OUT):: MSG
      INTEGER:: I,J

      OK=.TRUE.
      MSG=''
      DO I=1,NND
       IF (PXYZ(I,3).GE.0.D0) THEN
        OK=.FALSE.
        WRITE(MSG,'(A,I0,A)') ' Error: dipole node ',I,' is not below z = 0; dipole panels must be fully submerged.'
        RETURN
       ENDIF
      ENDDO
      DO I=1,ND
       DO J=1,NH
        IF (SQRT((PC(I,1)-HC(J,1))**2+(PC(I,2)-HC(J,2))**2+(PC(I,3)-HC(J,3))**2).LT.1.D-6*PSZ(I)) THEN
         OK=.FALSE.
         WRITE(MSG,'(A,I0,A,I0,A)') ' Error: dipole panel ',I,' coincides with hull panel ',J,'.'
         RETURN
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DIPOLE_CHECK_GEOMETRY

!=======================================================================
!SUBROUTINE DIPOLE_CHECK_CLEARANCE
!   !   Every node of a dipole panel at least one panel size below z = 0. The wave part of the dipole rows is tested only
!   !   for deeper sheets (design note Sec. 3.2); next to the free surface the image Hessian and the logarithmic part of
!   !   G_w are not. OK false and MSG naming the first panel otherwise.
      SUBROUTINE DIPOLE_CHECK_CLEARANCE(PXYZ,PNCN,PNCON,PSZ,ND,OK,MSG)
      IMPLICIT   NONE
      INTEGER,INTENT(IN):: ND,PNCN(:),PNCON(:,:)
      REAL*8,INTENT(IN):: PXYZ(:,:),PSZ(:)
      LOGICAL,INTENT(OUT):: OK
      CHARACTER(LEN=*),INTENT(OUT):: MSG
      INTEGER:: I,K

      OK=.TRUE.
      MSG=''
      DO I=1,ND
       DO K=1,PNCN(I)
        IF (PXYZ(PNCON(I,K),3).GT.-PSZ(I)) THEN
         OK=.FALSE.
         WRITE(MSG,'(A,I0,A)') ' Error: dipole panel ',I,' has a node less than one panel size below z = 0;'// &
                               ' sheets this close to the free surface are not supported.'
         RETURN
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DIPOLE_CHECK_CLEARANCE

!=======================================================================
!SUBROUTINE DIPOLE_CHECK_COPLANAR
!   !   No dipole panel may overlap a coplanar hull panel: a sheet cannot close a body or lie on its surface (design
!   !   note Sec. 6). Coplanar: parallel normals (|n_D x n_H| < 1e-6, either direction) and the other centroid within
!   !   1e-6 panel sizes of the plane. Overlapping: the centroid of either panel inside the other's outline, projected
!   !   on the plane (edges included, to 1e-9 panel sizes). The centroid test does not see every partial overlap.
      SUBROUTINE DIPOLE_CHECK_COPLANAR(PXYZ,PNCN,PNCON,PC,PN,PSZ,ND,HXYZ,HNCN,HNCON,HC,HN,HSZ,NH,OK,MSG)
      IMPLICIT   NONE
      INTEGER,INTENT(IN):: ND,NH,PNCN(:),PNCON(:,:),HNCN(:),HNCON(:,:)
      REAL*8,INTENT(IN):: PXYZ(:,:),PC(:,:),PN(:,:),PSZ(:),HXYZ(:,:),HC(:,:),HN(:,:),HSZ(:)
      LOGICAL,INTENT(OUT):: OK
      CHARACTER(LEN=*),INTENT(OUT):: MSG
      INTEGER:: I,J
      REAL*8:: CR(3),SZ

      OK=.TRUE.
      MSG=''
      DO I=1,ND
       DO J=1,NH
        CR(1)=PN(I,2)*HN(J,3)-PN(I,3)*HN(J,2)
        CR(2)=PN(I,3)*HN(J,1)-PN(I,1)*HN(J,3)
        CR(3)=PN(I,1)*HN(J,2)-PN(I,2)*HN(J,1)
        IF (SQRT(CR(1)**2+CR(2)**2+CR(3)**2).GE.1.D-6) CYCLE
        SZ=MAX(PSZ(I),HSZ(J))
        IF (ABS(DOT_PRODUCT(HC(J,:)-PC(I,:),PN(I,:))).GE.1.D-6*SZ) CYCLE
        IF (INSIDE_PANEL(HC(J,:),PXYZ,PNCON(I,:),PNCN(I),PN(I,:),SZ).OR.                                          &
            INSIDE_PANEL(PC(I,:),HXYZ,HNCON(J,:),HNCN(J),PN(I,:),SZ)) THEN
         OK=.FALSE.
         WRITE(MSG,'(A,I0,A,I0,A)') ' Error: dipole panel ',I,' overlaps hull panel ',J,                         &
                                   ' in the same plane; a sheet cannot close a body or lie on it.'
         RETURN
        ENDIF
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE DIPOLE_CHECK_COPLANAR

!=======================================================================
!FUNCTION INSIDE_PANEL
!   !   Whether X, in the plane of the panel with nodes XYZ(NCON(1:NN),:) and normal N, lies inside its outline or on it
!   !   (to 1e-9 SZ): on the same side of every edge, for either node order.
      LOGICAL FUNCTION INSIDE_PANEL(X,XYZ,NCON,NN,N,SZ)
      IMPLICIT   NONE
      INTEGER,INTENT(IN):: NN,NCON(:)
      REAL*8,INTENT(IN):: X(3),XYZ(:,:),N(3),SZ
      INTEGER:: K,KN
      REAL*8:: E(3),W(3),S
      LOGICAL:: POS,NEG

      POS=.FALSE.
      NEG=.FALSE.
      DO K=1,NN
       KN=MOD(K,NN)+1
       E(:)=XYZ(NCON(KN),:)-XYZ(NCON(K),:)
       W(:)=X(:)-XYZ(NCON(K),:)
       S=((E(2)*W(3)-E(3)*W(2))*N(1)+(E(3)*W(1)-E(1)*W(3))*N(2)+(E(1)*W(2)-E(2)*W(1))*N(3))
       S=S/MAX(SQRT(E(1)**2+E(2)**2+E(3)**2),1.D-300)
       IF (S.GT.1.D-9*SZ) POS=.TRUE.
       IF (S.LT.-1.D-9*SZ) NEG=.TRUE.
      ENDDO
      INSIDE_PANEL=.NOT.(POS.AND.NEG)

      RETURN
      END FUNCTION INSIDE_PANEL

!=======================================================================
!SUBROUTINE DIPOLE_RELEASE
!   !   Frees the dipole arrays and sets DNELEM = 0 (end of run, and between tests).
      SUBROUTINE DIPOLE_RELEASE
      IMPLICIT   NONE

      IF (ALLOCATED(DPXYZ)) DEALLOCATE(DPXYZ)
      IF (ALLOCATED(DNCN)) DEALLOCATE(DNCN)
      IF (ALLOCATED(DNCON)) DEALLOCATE(DNCON)
      IF (ALLOCATED(DPXYZ_C)) DEALLOCATE(DPXYZ_C)
      IF (ALLOCATED(DNRM)) DEALLOCATE(DNRM)
      IF (ALLOCATED(DDS)) DEALLOCATE(DDS)
      IF (ALLOCATED(DPNSZ)) DEALLOCATE(DPNSZ)
      DNELEM=0
      DNTND=0

      RETURN
      END SUBROUTINE DIPOLE_RELEASE
      
!-------------------------------------------------------------------------------
    END MODULE ReadPanelMesh
!*******************************************************************************
