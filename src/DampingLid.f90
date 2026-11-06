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
   USE WaveDyn_mod,     ONLY: V,H,WVN,NK
   USE INFG3D_Open,     ONLY: INFGREEN3D,HavelockGF
   USE FinGreen3D_Open, ONLY: FINGREEN3D

   IMPLICIT NONE

   ! ..... Public Subroutines ...................................................................................................

   PUBLIC :: READ_LID_MESH
   PUBLIC :: LID_GEOMETRY
   PUBLIC :: CHECK_LID
   PUBLIC :: ALLOCATE_LID_UNKNOWNS
   PUBLIC :: HULL_SOLID_ANGLE
   PUBLIC :: LID_SINGLE_LAYER
   PUBLIC :: LID_DEEP_WAVE

   INTEGER,PRIVATE:: NBODY_LID=1                                    ! set by CHECK_LID: selects the single or multi-body hull arrays
   INTEGER,PARAMETER,PRIVATE:: NGAUSS=6                             ! Gauss-Legendre points per direction of the wave term rules
   REAL*8,PARAMETER,PRIVATE:: RFAR=3.D0                             ! one point beyond RFAR panel sizes
   REAL*8,PARAMETER,PRIVATE:: RHMIN=0.03D0                          ! FINGREEN3D is sampled at R/h of at least RHMIN

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

END MODULE DampingLid
