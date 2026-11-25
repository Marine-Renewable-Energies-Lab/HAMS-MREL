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
!  Rankine integrals of a flat panel with a constant normal dipole density (zero-thickness panels, first order).
!  Design and sources: docs/developer_documentation.md, "Dipole (zero-thickness) panels".
!
!  The panel is given by its nodes A(1:NNODE,3), NNODE = 3 or 4, and a unit normal NQ. Both routines are oriented about
!  NQ, whatever the node order: a panel whose node-order normal (right-hand rule) points against NQ is traversed backwards.
!
!    DIPOLE_SOLID_ANGLE  Omega(x) = int (x - xi).NQ / |x - xi|^3 dS, the signed solid angle of the panel, +2 pi just on
!                        its +NQ side (Van Oosterom and Strackee 1983, per triangle; a quadrilateral as two triangles).
!    DIPOLE_GRAD_OMEGA   grad_x Omega, the Biot-Savart integral of a unit vortex ring on the panel edges (Hess and Smith
!                        1964, Sec. 4.2; Katz and Plotkin 2001, Ch. 10), Eqs. (7)-(8) of the design note.
!    DIPOLE_RANKINE_NN   NF . grad Omega = the finite-part integral of d2(1/r)/dn_F dn_Q over the panel, plus, when SV
!                        is not zero, SV times the same for the image 1/r1 (field point mirrored in z = 0).
!
!  On the plane of a triangle the solid angle is set to 0 (exact outside it, the principal value on it). On the line of an
!  edge the edge is skipped (its term is 0 on the line outside the edge, the principal value on it).
!-------------------------------------------------------------------------------
MODULE DipoleIntgr

   USE INFG3D_Open, ONLY: INFGREEN3D
   USE FinGreen3D_Open, ONLY: FINGREEN3D
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: DIPOLE_SOLID_ANGLE
   PUBLIC :: DIPOLE_GRAD_OMEGA
   PUBLIC :: DIPOLE_RANKINE_NN
   PUBLIC :: GREEN_WAVE_GRAD
   PUBLIC :: GREEN_WAVE_FIELD_GRAD
   PUBLIC :: GREEN_WAVE_NN

   REAL*8,PARAMETER:: TOL_PLANE=1.D-12          ! |triple product| / size^3 below which the point is in the triangle's plane
   REAL*8,PARAMETER:: TOL_EDGE=1.D-9            ! distance to an edge line / size below which the edge is skipped

CONTAINS

!-------------------------------------------------------------------------------
!     +1 if the node-order normal of the panel points along NQ, -1 otherwise; SZ is the panel size (longest diagonal or
!     side), used to scale the tolerances.
!-------------------------------------------------------------------------------
      SUBROUTINE PANEL_ORIENT(A,NNODE,NQ,SGN,SZ)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NNODE
      REAL*8,INTENT(IN):: A(4,3),NQ(3)
      REAL*8,INTENT(OUT):: SGN,SZ
      REAL*8:: M(3),U(3),W(3)
      INTEGER:: I,J

      IF (NNODE.EQ.4) THEN
       U(:)=A(3,:)-A(1,:)
       W(:)=A(4,:)-A(2,:)
      ELSE
       U(:)=A(2,:)-A(1,:)
       W(:)=A(3,:)-A(1,:)
      ENDIF
      M(1)=U(2)*W(3)-U(3)*W(2)
      M(2)=U(3)*W(1)-U(1)*W(3)
      M(3)=U(1)*W(2)-U(2)*W(1)
      SGN=1.D0
      IF (M(1)*NQ(1)+M(2)*NQ(2)+M(3)*NQ(3).LT.0.D0) SGN=-1.D0
      SZ=0.D0
      DO I=1,NNODE
       DO J=I+1,NNODE
        SZ=MAX(SZ,SQRT((A(I,1)-A(J,1))**2+(A(I,2)-A(J,2))**2+(A(I,3)-A(J,3))**2))
       ENDDO
      ENDDO

      RETURN
      END SUBROUTINE PANEL_ORIENT

!-------------------------------------------------------------------------------
!     Signed solid angle of the triangle P1, P2, P3 seen from XF, positive on the side its right-hand normal points to.
!     tan(Omega/2) = -R1.(R2 x R3) / (r1 r2 r3 + (R1.R2) r3 + (R1.R3) r2 + (R2.R3) r1), Ri = Pi - XF
!     (Van Oosterom and Strackee 1983; the minus sign gives the orientation above [derived, design note Sec. 3.1]).
!-------------------------------------------------------------------------------
      REAL*8 FUNCTION TRI_SOLID_ANGLE(XF,P1,P2,P3,SZ)
      IMPLICIT NONE
      REAL*8,INTENT(IN):: XF(3),P1(3),P2(3),P3(3),SZ
      REAL*8:: R1(3),R2(3),R3(3),A1,A2,A3,TRIP,DEN

      R1(:)=P1(:)-XF(:)
      R2(:)=P2(:)-XF(:)
      R3(:)=P3(:)-XF(:)
      A1=SQRT(R1(1)**2+R1(2)**2+R1(3)**2)
      A2=SQRT(R2(1)**2+R2(2)**2+R2(3)**2)
      A3=SQRT(R3(1)**2+R3(2)**2+R3(3)**2)
      TRIP=R1(1)*(R2(2)*R3(3)-R2(3)*R3(2))+R1(2)*(R2(3)*R3(1)-R2(1)*R3(3))+R1(3)*(R2(1)*R3(2)-R2(2)*R3(1))
      IF (ABS(TRIP).LT.TOL_PLANE*SZ**3) THEN
       TRI_SOLID_ANGLE=0.D0
       RETURN
      ENDIF
      DEN=A1*A2*A3+DOT_PRODUCT(R1,R2)*A3+DOT_PRODUCT(R1,R3)*A2+DOT_PRODUCT(R2,R3)*A1
      TRI_SOLID_ANGLE=-2.D0*ATAN2(TRIP,DEN)

      RETURN
      END FUNCTION TRI_SOLID_ANGLE

!-------------------------------------------------------------------------------
!     Omega(XF) of the panel A(1:NNODE,:) oriented about NQ.
!-------------------------------------------------------------------------------
      REAL*8 FUNCTION DIPOLE_SOLID_ANGLE(XF,A,NNODE,NQ)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NNODE
      REAL*8,INTENT(IN):: XF(3),A(4,3),NQ(3)
      REAL*8:: SGN,SZ,OM

      CALL PANEL_ORIENT(A,NNODE,NQ,SGN,SZ)
      OM=TRI_SOLID_ANGLE(XF,A(1,:),A(2,:),A(3,:),SZ)
      IF (NNODE.EQ.4) OM=OM+TRI_SOLID_ANGLE(XF,A(1,:),A(3,:),A(4,:),SZ)
      DIPOLE_SOLID_ANGLE=SGN*OM

      RETURN
      END FUNCTION DIPOLE_SOLID_ANGLE

!-------------------------------------------------------------------------------
!     grad_x Omega at XF for the panel A(1:NNODE,:) oriented about NQ. Per edge from P to Q (design note Eq. 8):
!       -(r1 x r2) / |r1 x r2|^2 * r0.(r1/|r1| - r2/|r2|),  r1 = XF - P, r2 = XF - Q, r0 = Q - P,
!     summed over the edges in node order, which run right-handed about the node-order normal.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_GRAD_OMEGA(XF,A,NNODE,NQ,GRAD)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NNODE
      REAL*8,INTENT(IN):: XF(3),A(4,3),NQ(3)
      REAL*8,INTENT(OUT):: GRAD(3)
      INTEGER:: K,KN
      REAL*8:: SGN,SZ,R1(3),R2(3),R0(3),C(3),C2,L0,A1,A2,F

      CALL PANEL_ORIENT(A,NNODE,NQ,SGN,SZ)
      GRAD=0.D0
      DO K=1,NNODE
       KN=MOD(K,NNODE)+1
       R0(:)=A(KN,:)-A(K,:)
       L0=SQRT(R0(1)**2+R0(2)**2+R0(3)**2)
       IF (L0.LT.1.D-12*SZ) CYCLE                                   ! repeated node
       R1(:)=XF(:)-A(K,:)
       R2(:)=XF(:)-A(KN,:)
       C(1)=R1(2)*R2(3)-R1(3)*R2(2)
       C(2)=R1(3)*R2(1)-R1(1)*R2(3)
       C(3)=R1(1)*R2(2)-R1(2)*R2(1)
       C2=C(1)**2+C(2)**2+C(3)**2
       IF (C2.LT.(TOL_EDGE*SZ*L0)**2) CYCLE                         ! on the edge line: distance = |C| / L0
       A1=SQRT(R1(1)**2+R1(2)**2+R1(3)**2)
       A2=SQRT(R2(1)**2+R2(2)**2+R2(3)**2)
       F=(R0(1)*(R1(1)/A1-R2(1)/A2)+R0(2)*(R1(2)/A1-R2(2)/A2)+R0(3)*(R1(3)/A1-R2(3)/A2))/C2
       GRAD(:)=GRAD(:)-F*C(:)
      ENDDO
      GRAD(:)=SGN*GRAD(:)

      RETURN
      END SUBROUTINE DIPOLE_GRAD_OMEGA

!-------------------------------------------------------------------------------
!     ENN = NF . grad Omega(XF) + SV * (R NF) . grad Omega(XF mirrored), R = diag(1, 1, -1): the finite-part integral of
!     d2(1/r)/dn_F dn_Q over the panel, and SV times that of d2(1/r1)/dn_F dn_Q (design note Sec. 3.2). SV = 0 gives the
!     1/r part alone; SV = 1 the image of G (HAMS's 1/r1), SV = -1 at the infinite-frequency limit.
!-------------------------------------------------------------------------------
      SUBROUTINE DIPOLE_RANKINE_NN(XF,NF,A,NNODE,NQ,SV,ENN)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NNODE
      REAL*8,INTENT(IN):: XF(3),NF(3),A(4,3),NQ(3),SV
      REAL*8,INTENT(OUT):: ENN
      REAL*8:: GRAD(3),XI(3)

      CALL DIPOLE_GRAD_OMEGA(XF,A,NNODE,NQ,GRAD)
      ENN=NF(1)*GRAD(1)+NF(2)*GRAD(2)+NF(3)*GRAD(3)
      IF (SV.NE.0.D0) THEN
       XI(1)=XF(1)
       XI(2)=XF(2)
       XI(3)=-XF(3)
       CALL DIPOLE_GRAD_OMEGA(XI,A,NNODE,NQ,GRAD)
       ENN=ENN+SV*(NF(1)*GRAD(1)+NF(2)*GRAD(2)-NF(3)*GRAD(3))
      ENDIF

      RETURN
      END SUBROUTINE DIPOLE_RANKINE_NN

!-------------------------------------------------------------------------------
!     The wave part G_w of the Green function and its gradient with respect to the source point XQ, at the field point
!     XF: G = 1/r + 1/r1 + G_w (design note Sec. 1.2). INFGREEN3D (H < 0) or FINGREEN3D with FLAG = 1 return G_w + 1/r1
!     (-1/r1 at V = -1); the image is removed here as IMAGE_TERM_EXACT removes it. At V = 0 and V = -1, G_w = 0.
!     RRD and FIELD are passed to FINGREEN3D (the series choice of a difference stencil, see FINGREEN3D_R).
!-------------------------------------------------------------------------------
      SUBROUTINE GREEN_WAVE_GRAD(XF,XQ,V,WVN,NK,H,GW,RRD,FIELD)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NK
      REAL*8,INTENT(IN):: XF(3),XQ(3),V,WVN(NK),H
      COMPLEX*16,INTENT(OUT):: GW(4)
      REAL*8,INTENT(IN),OPTIONAL:: RRD
      LOGICAL,INTENT(IN),OPTIONAL:: FIELD
      REAL*8:: R1,SGN

      IF (ABS(V).LT.1.E-8.OR.ABS(V+1.D0).LT.1.E-8) THEN
       GW=CMPLX(0.D0,0.D0)
       RETURN
      ENDIF
      IF (H.LT.0.D0) THEN
       CALL INFGREEN3D(XQ(1),XF(1),XQ(2),XF(2),XQ(3),XF(3),V,GW,1)
      ELSE
       CALL FINGREEN3D(XQ(1),XF(1),XQ(2),XF(2),XQ(3),XF(3),V,WVN,NK,H,GW,1,RR_DECIDE=RRD,FIELD=FIELD)
      ENDIF
      SGN=1.D0
      R1=SQRT((XQ(1)-XF(1))**2+(XQ(2)-XF(2))**2+(XQ(3)+XF(3))**2)
      GW(1)=GW(1)-SGN/R1
      GW(2)=GW(2)+SGN*(XQ(1)-XF(1))/R1**3
      GW(3)=GW(3)+SGN*(XQ(2)-XF(2))/R1**3
      GW(4)=GW(4)+SGN*(XQ(3)+XF(3))/R1**3

      RETURN
      END SUBROUTINE GREEN_WAVE_GRAD

!-------------------------------------------------------------------------------
!     G_w and its gradient with respect to the field point XF: by the symmetry G(x, xi) = G(xi, x) (Wehausen and Laitone
!     1960, Sec. 13) this is the source gradient with the two points swapped.
!-------------------------------------------------------------------------------
      SUBROUTINE GREEN_WAVE_FIELD_GRAD(XF,XQ,V,WVN,NK,H,GW)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NK
      REAL*8,INTENT(IN):: XF(3),XQ(3),V,WVN(NK),H
      COMPLEX*16,INTENT(OUT):: GW(4)

      CALL GREEN_WAVE_GRAD(XQ,XF,V,WVN,NK,H,GW)

      RETURN
      END SUBROUTINE GREEN_WAVE_FIELD_GRAD

!-------------------------------------------------------------------------------
!     GNN = d2 G_w / dn_F dn_Q, the mixed second derivative of the wave part along the field normal NF at XF and the
!     source normal NQ at XQ (design note Sec. 3.2, plan Task 2): central differences in XF along NF of NQ . grad_Q G_w,
!     step DELTA = 1e-4 min(R1, 1/k0), at least 1e-5 m (R1 the image distance, k0 = V in deep water, WVN(1) otherwise).
!     In finite depth the stencil takes the series of its centre (RR_DECIDE) and the series without Limes (FIELD).
!-------------------------------------------------------------------------------
      SUBROUTINE GREEN_WAVE_NN(XF,NF,XQ,NQ,V,WVN,NK,H,GNN)
      IMPLICIT NONE
      INTEGER,INTENT(IN):: NK
      REAL*8,INTENT(IN):: XF(3),NF(3),XQ(3),NQ(3),V,WVN(NK),H
      COMPLEX*16,INTENT(OUT):: GNN
      REAL*8:: R1,RR,K0,DELTA,XP(3),XM(3)
      COMPLEX*16:: GP(4),GM(4)

      IF (ABS(V).LT.1.E-8.OR.ABS(V+1.D0).LT.1.E-8) THEN
       GNN=CMPLX(0.D0,0.D0)
       RETURN
      ENDIF
      R1=SQRT((XQ(1)-XF(1))**2+(XQ(2)-XF(2))**2+(XQ(3)+XF(3))**2)
      RR=SQRT((XQ(1)-XF(1))**2+(XQ(2)-XF(2))**2)
      IF (H.LT.0.D0) THEN
       K0=V
      ELSE
       K0=WVN(1)
      ENDIF
      DELTA=MAX(1.D-4*MIN(R1,1.D0/K0),1.D-5)
      XP(:)=XF(:)+DELTA*NF(:)
      XM(:)=XF(:)-DELTA*NF(:)
      IF (H.LT.0.D0) THEN
       CALL GREEN_WAVE_GRAD(XP,XQ,V,WVN,NK,H,GP)
       CALL GREEN_WAVE_GRAD(XM,XQ,V,WVN,NK,H,GM)
      ELSE
       CALL GREEN_WAVE_GRAD(XP,XQ,V,WVN,NK,H,GP,RRD=RR,FIELD=.TRUE.)
       CALL GREEN_WAVE_GRAD(XM,XQ,V,WVN,NK,H,GM,RRD=RR,FIELD=.TRUE.)
      ENDIF
      GNN=((GP(2)-GM(2))*NQ(1)+(GP(3)-GM(3))*NQ(2)+(GP(4)-GM(4))*NQ(3))/(2.D0*DELTA)

      RETURN
      END SUBROUTINE GREEN_WAVE_NN
!-------------------------------------------------------------------------------
END MODULE DipoleIntgr
!*******************************************************************************
