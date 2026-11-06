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

   USE HAMS_mod,       ONLY: ISOLV,IRSP,NTHREAD
   USE Const_mod
   USE PanelMesh_mod
   USE Inerfs_mod,      ONLY: iNELEM,iXYZ_P
   USE Lid_mod
   USE NormalProcess
   USE SingularIntgr
   USE SingularIntgrMulti
   USE WaveDyn_mod,     ONLY: V,H,WVN,NK,BETA
   USE Body_mod,        ONLY: XW
   USE Potentials_mod,  ONLY: MXPOT
   USE INFG3D_Open,     ONLY: INFGREEN3D,HavelockGF
   USE FinGreen3D_Open, ONLY: FINGREEN3D
   USE PatcVelct,       ONLY: VINP,DINP
   USE AssbMatx,        ONLY: FILL_HULL_LEFT,ASSB_RBC,ASSB_DBC

   IMPLICIT NONE

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: READ_LID_MESH
   PUBLIC :: LID_GEOMETRY
   PUBLIC :: CHECK_LID
   PUBLIC :: ALLOCATE_LID_UNKNOWNS
   PUBLIC :: HULL_SOLID_ANGLE
   PUBLIC :: LID_SINGLE_LAYER
   PUBLIC :: LID_DEEP_WAVE
   PUBLIC :: LID_ASSEMBLE
   PUBLIC :: LID_RADIATION
   PUBLIC :: LID_DIFFRACTION
   PUBLIC :: SYM_SIGNS
   PUBLIC :: REFLECTED

   INTEGER,PRIVATE:: NBODY_LID=1                                    ! set by CHECK_LID: selects the single or multi-body hull arrays
   INTEGER,PARAMETER,PRIVATE:: NGAUSS=6                             ! Gauss-Legendre points per direction of the wave term rules
   REAL*8,PARAMETER,PRIVATE:: RFAR=3.D0                             ! one point beyond RFAR panel sizes
   REAL*8,PARAMETER,PRIVATE:: RHMIN=0.03D0                          ! FINGREEN3D is sampled at R/h of at least RHMIN

   ! ..... The single-body system with a lid, built by LID_ASSEMBLE once per frequency ...........................................
   !   Unknowns [phi(1:NELEM), q(1:LNELEM)]; rows: hull, lid, then the INELEM interior rows when IRSP = 1. NROW of the rows
   !   after the hull have field points XROW: the lid centroids, then the interior waterplane centroids.
   INTEGER,PRIVATE:: NROW=0
   REAL*8,ALLOCATABLE,PRIVATE:: XROW(:,:)
   COMPLEX*16,ALLOCATABLE,PRIVATE:: ALID(:,:,:)                     ! LU factors when square (IRSP = 0), the matrix itself if not
   INTEGER,ALLOCATABLE,PRIVATE:: IPIVL(:,:)
   COMPLEX*16,ALLOCATABLE,PRIVATE:: SROW(:,:,:)                     ! single layer of hull panel JEL, symmetry image IS, at XROW(K)

CONTAINS
!=======================================================================
!SUBROUTINE READ_LID_MESH

!   !   Reads DampingLid.pnl, in the format of WaterPlaneMesh.pnl. Unlike the waterplane read in HAMS_Prog, the symmetry
!   !   flags of the file are compared with the hull's and never stored over them. Stops on a read error, on a panel with other
!   !   than 3 or 4 vertices, and on a vertex number outside 1 to the number of nodes.

      SUBROUTINE READ_LID_MESH(INPUTDIR)
      IMPLICIT   NONE

      CHARACTER(LEN=*),INTENT(IN):: INPUTDIR
      INTEGER:: IU,IND,IEL,J,M,LISX,LISY,IOS
      CHARACTER(LEN=512):: LINE

      OPEN(NEWUNIT=IU,FILE=TRIM(INPUTDIR)//"/DampingLid.pnl",STATUS="OLD",ACTION="READ",IOSTAT=IOS)
      IF (IOS.NE.0) CALL LID_READ_ERROR('the file cannot be opened',0)

      DO J=1,3
       READ(IU,*,IOSTAT=IOS)
       IF (IOS.NE.0) CALL LID_READ_ERROR('header line',J)
      ENDDO
      READ(IU,*,IOSTAT=IOS) LNELEM,LNTND,LISX,LISY
      IF (IOS.NE.0) CALL LID_READ_ERROR('line of the numbers of panels and nodes and the symmetry flags',0)
      IF (LNELEM.LT.1.OR.LNTND.LT.3) CALL LID_READ_ERROR('line of the numbers of panels and nodes (need at least 1 and 3)',0)
      IF (LISX.NE.ISX.OR.LISY.NE.ISY) THEN
       PRINT *,' Error: lid symmetry differs from the hull.'
       PRINT *,' ISX, ISY of DampingLid.pnl:',LISX,LISY,'  of the hull:',ISX,ISY
       STOP 1
      ENDIF

      ALLOCATE(LXYZ(LNTND,3),LNCN(LNELEM),LNCON(LNELEM,4))
      LNCON=0

      DO J=1,2
       READ(IU,*,IOSTAT=IOS)
       IF (IOS.NE.0) CALL LID_READ_ERROR('header line before the nodes',J)
      ENDDO
      DO IND=1,LNTND
       READ(IU,*,IOSTAT=IOS) M,LXYZ(IND,1),LXYZ(IND,2),LXYZ(IND,3)
       IF (IOS.NE.0) CALL LID_READ_ERROR('node line, node No.',IND)
       IF (ABS(LXYZ(IND,3)).GT.1.D-10) THEN
        PRINT *,' Error: lid node is not on z = 0, node No.',IND
        STOP 1
       ENDIF
      ENDDO

      DO J=1,3
       READ(IU,*,IOSTAT=IOS)
       IF (IOS.NE.0) CALL LID_READ_ERROR('header line before the panels',J)
      ENDDO
      DO IEL=1,LNELEM
       READ(IU,'(A)',IOSTAT=IOS) LINE                               ! the vertex count is checked before the vertices are read
       IF (IOS.EQ.0) READ(LINE,*,IOSTAT=IOS) M,LNCN(IEL)
       IF (IOS.NE.0) CALL LID_READ_ERROR('panel line, panel No.',IEL)
       IF (LNCN(IEL).NE.3.AND.LNCN(IEL).NE.4) THEN
        PRINT *,' Error: a lid panel in DampingLid.pnl must have 3 or 4 vertices, panel No.',IEL,' has',LNCN(IEL)
        STOP 1
       ENDIF
       READ(LINE,*,IOSTAT=IOS) M,LNCN(IEL),(LNCON(IEL,J),J=1,LNCN(IEL))
       IF (IOS.NE.0) CALL LID_READ_ERROR('panel line, panel No.',IEL)
       IF (ANY(LNCON(IEL,1:LNCN(IEL)).LT.1).OR.ANY(LNCON(IEL,1:LNCN(IEL)).GT.LNTND)) THEN
        PRINT *,' Error: a vertex of a lid panel in DampingLid.pnl is not a node number, panel No.',IEL
        STOP 1
       ENDIF
      ENDDO

      CLOSE(IU)

      RETURN
      END SUBROUTINE READ_LID_MESH

!=======================================================================
!   !   Stops on a read error in DampingLid.pnl, naming the kind of line and its number (0: none).

      SUBROUTINE LID_READ_ERROR(WHAT,N)
      IMPLICIT   NONE

      CHARACTER(LEN=*),INTENT(IN):: WHAT
      INTEGER,INTENT(IN):: N

      IF (N.GT.0) THEN
       PRINT *,' Error: cannot read DampingLid.pnl, ',WHAT,N
      ELSE
       PRINT *,' Error: cannot read DampingLid.pnl, ',WHAT
      ENDIF
      STOP 1
      END SUBROUTINE LID_READ_ERROR

!=======================================================================
!SUBROUTINE LID_GEOMETRY

!   !   Centroids, areas and chart lengths of the lid panels. No normal: the lid has none, so the vertex order of a panel
!   !   does not matter. Stops on a panel with a repeated vertex or with zero area.

      SUBROUTINE LID_GEOMETRY()
      IMPLICIT   NONE

      INTEGER:: IEL,I,J

      ALLOCATE(LXYZ_P(LNELEM,3),LDS(LNELEM),LPNSZ(LNELEM))

      CALL CalPanelCentre(LXYZ,LNTND,LNELEM,LNCN,LNCON,LXYZ_P)
      CALL CalPanelArea(LXYZ,LNTND,LNELEM,LNCN,LNCON,LDS)
      CALL CalPanelChartLength(LXYZ,LNTND,LNELEM,LNCN,LNCON,LPNSZ)

      DO IEL=1,LNELEM
       DO I=1,LNCN(IEL)-1
        DO J=I+1,LNCN(IEL)
         IF (LNCON(IEL,I).EQ.LNCON(IEL,J)) THEN
          PRINT *,' Error: a lid panel in DampingLid.pnl has a repeated vertex, panel No.',IEL
          STOP 1
         ENDIF
        ENDDO
       ENDDO
       IF (LDS(IEL).LT.1.D-12*LPNSZ(IEL)**2) THEN
        PRINT *,' Error: a lid panel in DampingLid.pnl has zero area, panel No.',IEL
        STOP 1
       ENDIF
      ENDDO

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

!   !   The stops that need the hull: the solver type, IRSP with several bodies, and a lid panel inside a hull. The file and
!   !   node checks are made in the control file read and in READ_LID_MESH. Call after the hull normals and LID_GEOMETRY.

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

      IF (NBODYIN.GT.1.AND.IRSP.EQ.1) THEN
       PRINT *,' Error: damping lid with several bodies needs If_remove_irr_freq 0.'
       PRINT *,' The interior rows of the multi-body irregular frequency removal are not combined with the lid'
       PRINT *,' (test/data/damping_lid/README.md). Use one body, or no irregular frequency removal.'
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

!=======================================================================
!SUBROUTINE LID_SINGLE_LAYER

!   !   SL, the integral of the Green function G over lid panel LEL, for a field point XF with z <= 0. The parts:
!   !   - 1/r and its image 1/r1, exact panel integrals (IRR = 5). The image takes the sign of IMAGE_TERM_EXACT: -1 at
!   !     infinite frequency, where G = 1/r - 1/r1.
!   !   - the deep-water wave term -V*GF of HavelockGF, by LID_DEEP_WAVE.
!   !   - in finite depth (H > 0), the area times the remainder G_finite - G_deep between the centroid and the depth of XF,
!   !     at the horizontal distance of the two but at least RHMIN*H, below which FINGREEN3D scatters
!   !     (test/data/damping_lid/README.md). The remainder has no log singularity, since the log coefficients of the two
!   !     match, as the README measures.
!   !   At V = 0 and at infinite frequency there is no wave term. V, WVN, NK and H are those of WaveDyn_mod.

      SUBROUTINE LID_SINGLE_LAYER(XF,LEL,SL)
      IMPLICIT   NONE

      REAL*8,INTENT(IN):: XF(3)
      INTEGER,INTENT(IN):: LEL
      COMPLEX*16,INTENT(OUT):: SL
      REAL*8:: XI(3),RK(4),RKI(4),SGN,RH
      COMPLEX*16:: GFN(4),GDP(4)

      CALL LID_RANKINE(XF,LEL,RK)
      IF (XF(3).EQ.0.D0) THEN
       RKI=RK                                                        ! on z = 0 the image integral is the direct one
      ELSE
       XI(1)=XF(1)
       XI(2)=XF(2)
       XI(3)=-XF(3)
       CALL LID_RANKINE(XI,LEL,RKI)
      ENDIF

      SGN=1.D0
      IF (ABS(V+1.D0).LT.1.E-8) SGN=-1.D0                           ! the limit of infinite frequency, as in IMAGE_TERM_EXACT

      SL=CMPLX(RK(1)+SGN*RKI(1),0.D0,KIND=8)

      IF (ABS(V).LT.1.E-8.OR.ABS(V+1.D0).LT.1.E-8) RETURN

      SL=SL+LID_DEEP_WAVE(XF,LEL,NGAUSS)

      IF (H.GT.0.D0) THEN
       RH=SQRT((XF(1)-LXYZ_P(LEL,1))**2+(XF(2)-LXYZ_P(LEL,2))**2)
       RH=MAX(RH,RHMIN*H)
       ! The same 1/R1 is added to both with FLAG = 1, so the difference is that of the wave parts
       CALL FINGREEN3D(RH,0.D0,0.D0,0.D0,LXYZ_P(LEL,3),XF(3),V,WVN,NK,H,GFN,1)
       CALL INFGREEN3D(RH,0.D0,0.D0,0.D0,LXYZ_P(LEL,3),XF(3),V,GDP,1)
       SL=SL+(GFN(1)-GDP(1))*LDS(LEL)
      ENDIF

      RETURN
      END SUBROUTINE LID_SINGLE_LAYER

!=======================================================================
!FUNCTION LID_DEEP_WAVE

!   !   The integral of the deep-water wave term -V*GF over lid panel LEL, at field point XF, with NG x NG Gauss-Legendre
!   !   points per rule. GF has a log singularity at zero distance on z = 0, so the rule depends on where XF is:
!   !   - farther than RFAR panel sizes from the centroid: one point at the centroid, as elsewhere in the code.
!   !   - projection of XF inside the panel or on its edge: the panel is split into triangles at the projection, and each
!   !     triangle is integrated by a Duffy transform, which removes the singularity at its vertex. Two refinements make
!   !     6 x 6 points enough (test duffy_converges): the radial coordinate u of the transform is graded as u = s**2, since
!   !     the log leaves u*log(u) in u and s**3*log(s) in s; and a triangle is split at the foot of the perpendicular from
!   !     the projection to its edge, where the log of the distance along the edge varies fastest. A triangle of zero area,
!   !     next to a node or edge that holds the projection, is skipped.
!   !   - otherwise: Gauss on the panel.
!   !   Valid only for V neither 0 nor -1.

      COMPLEX*16 FUNCTION LID_DEEP_WAVE(XF,LEL,NG)
      IMPLICIT   NONE

      REAL*8,INTENT(IN):: XF(3)
      INTEGER,INTENT(IN):: LEL,NG
      INTEGER:: K,N
      REAL*8:: XN(4,3),P(3),A(3),B(3),F(3),DIST,CRS(4),TOL,AREA,ALFA
      COMPLEX*16:: SUMW

      N=LNCN(LEL)
      DO K=1,N
       XN(K,:)=LXYZ(LNCON(LEL,K),:)
      ENDDO

      DIST=SQRT(SUM((XF-LXYZ_P(LEL,:))**2))
      IF (DIST.GT.RFAR*LPNSZ(LEL)) THEN
       LID_DEEP_WAVE=DEEP_WAVE_AT(XF,LXYZ_P(LEL,:))*LDS(LEL)
       RETURN
      ENDIF

      ! Projection of XF on the panel plane, and on which side of each edge it lies
      P(1)=XF(1)
      P(2)=XF(2)
      P(3)=LXYZ_P(LEL,3)
      DO K=1,N
       CRS(K)=EDGE_CROSS(XN(K,:),XN(MOD(K,N)+1,:),P)
      ENDDO
      TOL=1.D-10*LPNSZ(LEL)**2

      IF (ALL(CRS(1:N).GE.-TOL).OR.ALL(CRS(1:N).LE.TOL)) THEN
       SUMW=(0.D0,0.D0)
       DO K=1,N
        AREA=0.5D0*ABS(CRS(K))
        IF (AREA.LE.1.D-10*LDS(LEL)) CYCLE
        A=XN(K,:)
        B=XN(MOD(K,N)+1,:)
        ALFA=SUM((P-A)*(B-A))/SUM((B-A)**2)
        IF (ALFA.GT.0.D0.AND.ALFA.LT.1.D0) THEN                      ! split at the foot of the perpendicular from P
         F=A+ALFA*(B-A)
         SUMW=SUMW+MAPPED_RULE(XF,P,A,F,P,NG,.TRUE.)+MAPPED_RULE(XF,P,F,B,P,NG,.TRUE.)
        ELSE
         SUMW=SUMW+MAPPED_RULE(XF,P,A,B,P,NG,.TRUE.)
        ENDIF
       ENDDO
       LID_DEEP_WAVE=SUMW
      ELSEIF (N.EQ.4) THEN
       LID_DEEP_WAVE=MAPPED_RULE(XF,XN(1,:),XN(2,:),XN(3,:),XN(4,:),NG,.FALSE.)
      ELSE
       LID_DEEP_WAVE=MAPPED_RULE(XF,XN(1,:),XN(2,:),XN(3,:),XN(3,:),NG,.FALSE.)
      ENDIF

      RETURN
      END FUNCTION LID_DEEP_WAVE

!=======================================================================
!   !   z-component of (B - A) x (P - A): positive when P is to the left of the edge A -> B, seen from above.

      REAL*8 FUNCTION EDGE_CROSS(A,B,P)
      IMPLICIT   NONE

      REAL*8,INTENT(IN):: A(3),B(3),P(3)

      EDGE_CROSS=(B(1)-A(1))*(P(2)-A(2))-(B(2)-A(2))*(P(1)-A(1))

      RETURN
      END FUNCTION EDGE_CROSS

!=======================================================================
!   !   -V*GF of HavelockGF between field point XF and source point XQ, without the Rankine parts.

      COMPLEX*16 FUNCTION DEEP_WAVE_AT(XF,XQ)
      IMPLICIT   NONE

      REAL*8,INTENT(IN):: XF(3),XQ(3)
      REAL*8:: RR
      COMPLEX*16:: GF,GFH

      RR=SQRT((XF(1)-XQ(1))**2+(XF(2)-XQ(2))**2)
      CALL HavelockGF(V*RR,V*(XF(3)+XQ(3)),GF,GFH)
      DEEP_WAVE_AT=-V*GF

      RETURN
      END FUNCTION DEEP_WAVE_AT

!=======================================================================
!   !   NG x NG Gauss-Legendre rule for the deep-water wave term over the bilinear map of the unit square onto X1 X2 X3 X4.
!   !   With X4 = X1 the map is a triangle collapsed at X1, whose Jacobian vanishes there: the Duffy transform. GRADE puts
!   !   s = sigma**2, with the Gauss points in sigma, to resolve the log singularity at X1.

      COMPLEX*16 FUNCTION MAPPED_RULE(XF,X1,X2,X3,X4,NG,GRADE)
      IMPLICIT   NONE

      REAL*8,INTENT(IN):: XF(3),X1(3),X2(3),X3(3),X4(3)
      INTEGER,INTENT(IN):: NG
      LOGICAL,INTENT(IN):: GRADE
      INTEGER:: I,J
      REAL*8:: XG(NG),WG(NG),S,WS,T,X(3),DXS(3),DXT(3),JAC
      COMPLEX*16:: SUMW

      CALL GAUSS_LEGENDRE01(NG,XG,WG)
      SUMW=(0.D0,0.D0)
      DO I=1,NG
       IF (GRADE) THEN
        S=XG(I)**2
        WS=2.D0*XG(I)*WG(I)
       ELSE
        S=XG(I)
        WS=WG(I)
       ENDIF
       DO J=1,NG
        T=XG(J)
        X=(1.D0-S)*(1.D0-T)*X1+S*(1.D0-T)*X2+S*T*X3+(1.D0-S)*T*X4
        DXS=(1.D0-T)*(X2-X1)+T*(X3-X4)
        DXT=(1.D0-S)*(X4-X1)+S*(X3-X2)
        JAC=SQRT((DXS(2)*DXT(3)-DXS(3)*DXT(2))**2+(DXS(3)*DXT(1)-DXS(1)*DXT(3))**2+(DXS(1)*DXT(2)-DXS(2)*DXT(1))**2)
        SUMW=SUMW+DEEP_WAVE_AT(XF,X)*JAC*WS*WG(J)
       ENDDO
      ENDDO
      MAPPED_RULE=SUMW

      RETURN
      END FUNCTION MAPPED_RULE

!=======================================================================
!   !   Gauss-Legendre points and weights on [0, 1], by Newton iteration on the Legendre polynomial of degree N.

      SUBROUTINE GAUSS_LEGENDRE01(N,X,W)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: N
      REAL*8,INTENT(OUT):: X(N),W(N)
      INTEGER:: I,J,IT
      REAL*8:: Z,P0,P1,P2,DP,DZ

      DO I=1,N
       Z=COS(PI*(I-0.25D0)/(N+0.5D0))
       DO IT=1,100
        P0=1.D0
        P1=Z
        DO J=2,N
         P2=((2*J-1)*Z*P1-(J-1)*P0)/J
         P0=P1
         P1=P2
        ENDDO
        DP=N*(Z*P1-P0)/(Z*Z-1.D0)
        DZ=P1/DP
        Z=Z-DZ
        IF (ABS(DZ).LT.1.D-15) EXIT
       ENDDO
       X(I)=0.5D0*(1.D0-Z)
       W(I)=1.D0/((1.D0-Z*Z)*DP*DP)
      ENDDO

      RETURN
      END SUBROUTINE GAUSS_LEGENDRE01

!=======================================================================
!   !   The 1/r integral over lid panel LEL at XF, IRR = 5.

      SUBROUTINE LID_RANKINE(XF,LEL,RK)
      IMPLICIT   NONE

      REAL*8,INTENT(IN):: XF(3)
      INTEGER,INTENT(IN):: LEL
      REAL*8,INTENT(OUT):: RK(4)
      REAL*8:: XC(3)

      XC=XF
      IF (LNCN(LEL).EQ.3) THEN
       CALL SGLINTBD_TRI2(XC,LEL,RK,5)
      ELSE
       CALL SGLINTBD_QUAD2(XC,LEL,RK,5)
      ENDIF

      RETURN
      END SUBROUTINE LID_RANKINE

!=======================================================================
!SUBROUTINE SYM_SIGNS

!   !   The signs of x and y between symmetry images IS and IP, in the order of the hull routines: with ISX = 1 and ISY = 0
!   !   the x sign is SY and the y sign SX, otherwise the other way round. IP = 1 reflects a field point as CALGREEN does.

      SUBROUTINE SYM_SIGNS(IS,IP,SGX,SGY)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: IS,IP
      REAL*8,INTENT(OUT):: SGX,SGY

      IF (ISX.EQ.1.AND.ISY.EQ.0) THEN
       SGX=SY(IS,IP)
       SGY=SX(IS,IP)
      ELSE
       SGX=SX(IS,IP)
       SGY=SY(IS,IP)
      ENDIF

      RETURN
      END SUBROUTINE SYM_SIGNS

!=======================================================================
!   !   Point X reflected for symmetry image IS, as CALGREEN reflects the field point.

      FUNCTION REFLECTED(X,IS) RESULT(XP)
      IMPLICIT   NONE

      REAL*8,INTENT(IN):: X(3)
      INTEGER,INTENT(IN):: IS
      REAL*8:: XP(3),SGX,SGY

      CALL SYM_SIGNS(IS,1,SGX,SGY)
      XP(1)=SGX*X(1)
      XP(2)=SGY*X(2)
      XP(3)=X(3)

      END FUNCTION REFLECTED

!=======================================================================
!SUBROUTINE HULL_LAYERS

!   !   The single layer SH and the double layer DH of hull panel JEL at field point XF, already reflected, built as the IRR = 3
!   !   loop of CALGREEN_IRR builds the interior rows: exact 1/r panel integrals and IMAGE_TERM_EXACT within 50 panel sizes,
!   !   one point beyond, with the wave part of the Green function at the centroid. DH is that of BODINT_LEFT, n.grad G at the
!   !   source with the panel normal GEONRM_P, and SH that of RBC_RIGHT. Full precision: CGRN and RKBN are not used.

      SUBROUTINE HULL_LAYERS(XF,JEL,SH,DH)
      IMPLICIT   NONE

      REAL*8,INTENT(IN):: XF(3)
      INTEGER,INTENT(IN):: JEL
      COMPLEX*16,INTENT(OUT):: SH,DH
      INTEGER:: FLAG
      REAL*8:: XQ(3),XC(3),RK(4),DIST
      COMPLEX*16:: GRN(4)

      XQ=XYZ_P(JEL,:)
      DIST=SQRT(SUM((XF-XQ)**2))                                     ! FLAG uses the distance to the reflected field point XF;
                                                                     ! CALGREEN uses the distance to the unreflected one
      IF (DIST.LE.50.D0*PNSZ(JEL)) THEN
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
       IF (NCN(JEL).EQ.3) THEN
        CALL SGLINTBD_TRI2(XC,JEL,RK,1)
       ELSE
        CALL SGLINTBD_QUAD2(XC,JEL,RK,1)
       ENDIF
       CALL IMAGE_TERM_EXACT(XC,XQ,V,JEL,1,GRN,RK(1),RK(2:4))
       SH=RK(1)+GRN(1)*DS(JEL)
       DH=(RK(2)+GRN(2)*DS(JEL))*GEONRM_P(JEL,1)+(RK(3)+GRN(3)*DS(JEL))*GEONRM_P(JEL,2)+(RK(4)+GRN(4)*DS(JEL))*GEONRM_P(JEL,3)
      ELSE
       SH=GRN(1)*DS(JEL)
       DH=(GRN(2)*GEONRM_P(JEL,1)+GRN(3)*GEONRM_P(JEL,2)+GRN(4)*GEONRM_P(JEL,3))*DS(JEL)
      ENDIF

      RETURN
      END SUBROUTINE HULL_LAYERS

!=======================================================================
!SUBROUTINE LID_ASSEMBLE

!   !   The left-hand side of the single-body problem with a lid, for each symmetry class IP, once per frequency. Unknowns are
!   !   the hull phi and the lid q = dphi/dz - nu*phi. Green's identity in the fluid, with the normal out of the fluid, gives
!   !   for a field point x (docs/superpowers/specs/2026-11-05-damping-lid-design.md, Formulation):
!   !     c phi(x) + sum_j D_xj phi_j - sum_l S_xl q_l = sum_j S_xj V_j,    c = 2 pi on the hull, 4 pi on the lid, 0 inside
!   !   since dG/dz = nu G on z = 0 leaves the free-surface integral -int G q on the lid only. The rows:
!   !   - hull: FILL_HULL_LEFT, and -S of the lid panels at the hull centroids.
!   !   - lid: phi = 0, so no free term; D of the hull panels and -S of the lid panels at the lid centroids.
!   !   - interior (IRSP = 1): the same at the interior waterplane centroids, where x is outside the fluid.
!   !   Symmetry images enter as in ASSB_LEFT: field points reflected for image IS, weighted by RXY(IS,IP). Square systems are
!   !   LU-factorised here; with interior rows the matrix is kept for the least-squares solve. Needs CALGREEN.

      SUBROUTINE LID_ASSEMBLE(NELEM,NSYS)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NELEM,NSYS
      INTEGER:: NC,IEL,JEL,K,L,IS,IP,INFO
      REAL*8:: XF(3)
      COMPLEX*16:: SL,SH,DH
      COMPLEX*16,ALLOCATABLE:: AH(:,:,:)

      NC=NELEM+LNELEM
      IF (.NOT.ALLOCATED(ALID)) THEN
       NROW=LNELEM
       IF (IRSP.EQ.1) NROW=LNELEM+iNELEM
       ALLOCATE(XROW(NROW,3),ALID(NELEM+NROW,NC,NSYS),IPIVL(NC,NSYS),SROW(NROW,NELEM,NSYS))
       XROW(1:LNELEM,:)=LXYZ_P
       IF (IRSP.EQ.1) XROW(LNELEM+1:NROW,:)=iXYZ_P
      ENDIF

      ALID=(0.D0,0.D0)
      ALLOCATE(AH(NELEM,NELEM,NSYS))
      CALL FILL_HULL_LEFT(AH,NELEM,NSYS)
      ALID(1:NELEM,1:NELEM,:)=AH
      DEALLOCATE(AH)

!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(IEL,IS,IP,L,XF,SL) SCHEDULE(DYNAMIC,16)
      DO IEL=1,NELEM
       DO IS=1,NSYS
        XF=REFLECTED(XYZ_P(IEL,:),IS)
        DO L=1,LNELEM
         CALL LID_SINGLE_LAYER(XF,L,SL)
         DO IP=1,NSYS
          ALID(IEL,NELEM+L,IP)=ALID(IEL,NELEM+L,IP)-RXY(IS,IP)*SL
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(K,IS,IP,JEL,L,XF,SL,SH,DH) SCHEDULE(DYNAMIC,16)
      DO K=1,NROW
       DO IS=1,NSYS
        XF=REFLECTED(XROW(K,:),IS)
        DO JEL=1,NELEM
         CALL HULL_LAYERS(XF,JEL,SH,DH)
         SROW(K,JEL,IS)=SH
         DO IP=1,NSYS
          ALID(NELEM+K,JEL,IP)=ALID(NELEM+K,JEL,IP)+RXY(IS,IP)*DH
         ENDDO
        ENDDO
        DO L=1,LNELEM
         CALL LID_SINGLE_LAYER(XF,L,SL)
         DO IP=1,NSYS
          ALID(NELEM+K,NELEM+L,IP)=ALID(NELEM+K,NELEM+L,IP)-RXY(IS,IP)*SL
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      IF (NROW.EQ.LNELEM) THEN
       DO IP=1,NSYS
        CALL ZGETRF(NC,NC,ALID(:,:,IP),NC,IPIVL(:,IP),INFO)
        IF (INFO.NE.0) THEN
         PRINT *,' Error: the system with the damping lid is singular, ZGETRF INFO =',INFO
         STOP 1
        ENDIF
       ENDDO
      ENDIF

      RETURN
      END SUBROUTINE LID_ASSEMBLE

!=======================================================================
!SUBROUTINE LID_SOLVE

!   !   Solves ALID X = B for the NRHS columns of B(NELEM+NROW,NRHS,NSYS), in place: by the LU factors of LID_ASSEMBLE when the
!   !   system is square, by ZGELS least squares with the interior rows. X is in the first NELEM+LNELEM rows of B.

      SUBROUTINE LID_SOLVE(B,NRHS,NELEM,NSYS)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NRHS,NELEM,NSYS
      COMPLEX*16,INTENT(INOUT):: B(NELEM+NROW,NRHS,NSYS)
      INTEGER:: NR,NC,IP,INFO,LWORK
      COMPLEX*16:: WQ(1)
      COMPLEX*16,ALLOCATABLE:: AW(:,:),WORK(:)

      NR=NELEM+NROW
      NC=NELEM+LNELEM
      IF (NR.EQ.NC) THEN
       DO IP=1,NSYS
        CALL ZGETRS('N',NC,NRHS,ALID(:,:,IP),NC,IPIVL(:,IP),B(:,:,IP),NR,INFO)
       ENDDO
      ELSE
       ALLOCATE(AW(NR,NC))
       CALL ZGELS('N',NR,NC,NRHS,AW,NR,B(:,:,1),NR,WQ,-1,INFO)
       LWORK=INT(REAL(WQ(1)))
       ALLOCATE(WORK(LWORK))
       DO IP=1,NSYS
        AW=ALID(:,:,IP)
        CALL ZGELS('N',NR,NC,NRHS,AW,NR,B(:,:,IP),NR,WORK,LWORK,INFO)
        IF (INFO.NE.0) THEN
         PRINT *,' Error: the least-squares system with the damping lid is rank deficient, ZGELS INFO =',INFO
         STOP 1
        ENDIF
       ENDDO
       DEALLOCATE(AW,WORK)
      ENDIF

      RETURN
      END SUBROUTINE LID_SOLVE

!=======================================================================
!SUBROUTINE LID_UNPACK

!   !   From the solution of each symmetry class in B, the hull phi of each symmetry image into MXPOT(:,MD0+1:MD0+NRHS,:) and
!   !   the lid q into MXQ, by the RXY recombination and the /NSYS of RADIATION_SOLVER.

      SUBROUTINE LID_UNPACK(B,NRHS,MD0,NELEM,NSYS)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NRHS,MD0,NELEM,NSYS
      COMPLEX*16,INTENT(IN):: B(NELEM+NROW,NRHS,NSYS)
      INTEGER:: M,IP,IS

      DO M=1,NRHS
       DO IP=1,NSYS
        MXPOT(:,MD0+M,IP)=(0.D0,0.D0)
        MXQ(:,MD0+M,IP)=(0.D0,0.D0)
        DO IS=1,NSYS
         MXPOT(:,MD0+M,IP)=MXPOT(:,MD0+M,IP)+B(1:NELEM,M,IS)*RXY(IP,IS)
         MXQ(:,MD0+M,IP)=MXQ(:,MD0+M,IP)+B(NELEM+1:NELEM+LNELEM,M,IS)*RXY(IP,IS)
        ENDDO
        MXPOT(:,MD0+M,IP)=MXPOT(:,MD0+M,IP)/NSYS
        MXQ(:,MD0+M,IP)=MXQ(:,MD0+M,IP)/NSYS
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE LID_UNPACK

!=======================================================================
!SUBROUTINE LID_RADIATION

!   !   The six radiation problems with the lid. Hull rows from ASSB_RBC. Lid and interior rows: sum_j S_xj V_j with S from
!   !   LID_ASSEMBLE and the mode signs of RBC_IRR; on the lid phi = 0, so nothing else. Fills MXPOT(:,1:6,:) and MXQ(:,1:6,:).

      SUBROUTINE LID_RADIATION(NELEM,NSYS)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NELEM,NSYS
      INTEGER:: K,JEL,IS,IP,MD
      REAL*8:: SGX,SGY,MS(6)
      COMPLEX*16:: BTMP(6,2)
      COMPLEX*16,ALLOCATABLE:: BH(:,:,:),B(:,:,:)

      ALLOCATE(BH(NELEM,6,NSYS),B(NELEM+NROW,6,NSYS))
      CALL ASSB_RBC(BH,NELEM,NSYS)
      B=(0.D0,0.D0)
      B(1:NELEM,:,:)=BH

!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(K,JEL,IS,IP,MD,SGX,SGY,MS,BTMP) SCHEDULE(DYNAMIC,16)
      DO K=1,NROW
       BTMP=(0.D0,0.D0)
       DO IS=1,NSYS
        DO IP=1,NSYS
         CALL SYM_SIGNS(IS,IP,SGX,SGY)
         MS=(/SGX,SGY,1.D0,SGY,SGX,SGX*SGY/)
         DO JEL=1,NELEM
          BTMP(:,IP)=BTMP(:,IP)+SROW(K,JEL,IS)*MS*DXYZ_P(JEL,:)
         ENDDO
        ENDDO
       ENDDO
       DO MD=1,6
        DO IP=1,NSYS
         DO IS=1,NSYS
          B(NELEM+K,MD,IP)=B(NELEM+K,MD,IP)+RXY(IP,IS)*BTMP(MD,IS)
         ENDDO
        ENDDO
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL LID_SOLVE(B,6,NELEM,NSYS)
      CALL LID_UNPACK(B,6,0,NELEM,NSYS)
      DEALLOCATE(BH,B)

      RETURN
      END SUBROUTINE LID_RADIATION

!=======================================================================
!SUBROUTINE LID_DIFFRACTION

!   !   The diffraction problem with the lid, for the heading BETA. The code solves for the scattered potential. Hull rows from
!   !   ASSB_DBC. Lid and interior rows: -sum_j S_xj dphi_I/dn_j, as DBC_IRR. On the lid the total potential is zero, so the
!   !   scattered one is -phi_I, and the free term 4 pi phi_S moves to the right as +4 pi phi_I(x), at x reflected for each
!   !   image. Interior rows have no free term. q is that of the total potential as well, since phi_I satisfies the free-surface
!   !   condition. At zero frequency the incident term is left out, as CalPotential leaves out the incident potential there.
!   !   Fills MXPOT(:,7,:) and MXQ(:,7,:).

      SUBROUTINE LID_DIFFRACTION(NELEM,NSYS)
      IMPLICIT   NONE

      INTEGER,INTENT(IN):: NELEM,NSYS
      INTEGER:: K,JEL,IS,IP
      REAL*8:: SGX,SGY,XF(3)
      COMPLEX*16:: BTMP(2),DPOX,DPOY,DPOZ
      COMPLEX*16,ALLOCATABLE:: BH(:,:),B(:,:,:),DINC(:,:,:)

      ALLOCATE(BH(NELEM,NSYS),B(NELEM+NROW,1,NSYS),DINC(NELEM,NSYS,NSYS))
      CALL ASSB_DBC(BH,NELEM,NSYS)
      B=(0.D0,0.D0)
      B(1:NELEM,1,:)=BH

      ! dphi_I/dn at hull panel JEL reflected between images IS and IP, as DBC_RIGHT
      DO IS=1,NSYS
       DO IP=1,NSYS
        CALL SYM_SIGNS(IS,IP,SGX,SGY)
        DO JEL=1,NELEM
         CALL DINP(SGX*XYZ_P(JEL,1),SGY*XYZ_P(JEL,2),XYZ_P(JEL,3),XW(1),XW(2),BETA,DPOX,DPOY,DPOZ)
         DINC(JEL,IS,IP)=SGX*DPOX*GEONRM_P(JEL,1)+SGY*DPOY*GEONRM_P(JEL,2)+DPOZ*GEONRM_P(JEL,3)
        ENDDO
       ENDDO
      ENDDO

!$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(K,JEL,IS,IP,XF,BTMP) SCHEDULE(DYNAMIC,16)
      DO K=1,NROW
       BTMP=(0.D0,0.D0)
       DO IS=1,NSYS
        DO IP=1,NSYS
         DO JEL=1,NELEM
          BTMP(IP)=BTMP(IP)-SROW(K,JEL,IS)*DINC(JEL,IS,IP)
         ENDDO
        ENDDO
       ENDDO
       IF (K.LE.LNELEM.AND.ABS(V).GE.1.D-8) THEN
        DO IS=1,NSYS
         XF=REFLECTED(XROW(K,:),IS)
         BTMP(IS)=BTMP(IS)+4.D0*PI*VINP(XF(1),XF(2),XF(3),XW(1),XW(2),BETA)
        ENDDO
       ENDIF
       DO IP=1,NSYS
        DO IS=1,NSYS
         B(NELEM+K,1,IP)=B(NELEM+K,1,IP)+RXY(IP,IS)*BTMP(IS)
        ENDDO
       ENDDO
      ENDDO
!$OMP END PARALLEL DO

      CALL LID_SOLVE(B,1,NELEM,NSYS)
      CALL LID_UNPACK(B,1,6,NELEM,NSYS)
      DEALLOCATE(BH,B,DINC)

      RETURN
      END SUBROUTINE LID_DIFFRACTION

END MODULE DampingLid
