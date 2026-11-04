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
!  Code Original Author:
!
!    Vaibhav Raghavan
!
!  ------------------------------------------------------------------------------------------------------
MODULE SurfaceVelocity
   ! Fluid velocity on the body panels from the potential at the panel centroids. HAMS solves for the potential, not
   ! its gradient, so the tangential part is obtained by differencing over the panel neighbours and the normal part is
   ! the known normal derivative.

   USE SecondOrderTypes, ONLY: panel_mesh, surface_geom, first_order_state
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: velocity_differencing

   ! Normalised determinant of the 2x2 normal equations below which the neighbours are treated as collinear
   REAL(8), PARAMETER :: DEGENERATE_DET = 1D-6

CONTAINS

   SUBROUTINE velocity_differencing(mesh, geom, s, n_fallback)
      ! Fills s%grad_p, the gradient of the body-disturbance potential phi_p at the panel centroids:
      !   grad_p = n dphi_p_dn + t1 d(phi_p)/dt1 + t2 d(phi_p)/dt2
      ! with orthonormal tangents t1, t2. The tangential derivatives are the weighted least-squares fit (weight
      ! 1/|d|^2) of phi_j - phi_i over the neighbours j, with d the centroid offset projected on the tangents. The total
      ! velocity is s%grad_i + s%grad_p. A panel with fewer than two non-collinear neighbours keeps only its normal
      ! part, and is counted in n_fallback. s%grad_p must be allocated, as build_state does.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(first_order_state), INTENT(INOUT) :: s
      INTEGER, INTENT(OUT) :: n_fallback

      INTEGER :: i, m, j
      REAL(8) :: nv(3), t1(3), t2(3), a(3), d(3), d1, d2, r2, w
      REAL(8) :: a11, a12, a22, det
      COMPLEX(8) :: b1, b2, dphi, g1, g2

      n_fallback = 0
      DO i = 1, mesh%n
         nv = mesh%nrm(:, i)
         s%grad_p(:, i) = nv*s%dphi_p_dn(i)

         ! Tangents from the coordinate axis least aligned with the normal
         a = 0D0
         a(MINLOC(ABS(nv), 1)) = 1D0
         t1 = a - DOT_PRODUCT(a, nv)*nv
         t1 = t1/SQRT(DOT_PRODUCT(t1, t1))
         t2 = [nv(2)*t1(3) - nv(3)*t1(2), nv(3)*t1(1) - nv(1)*t1(3), nv(1)*t1(2) - nv(2)*t1(1)]

         a11 = 0D0; a12 = 0D0; a22 = 0D0
         b1 = (0D0, 0D0); b2 = (0D0, 0D0)
         DO m = geom%nbr_start(i), geom%nbr_start(i + 1) - 1
            j = geom%nbr(m)
            d = mesh%xc(:, j) - mesh%xc(:, i)
            d1 = DOT_PRODUCT(d, t1)
            d2 = DOT_PRODUCT(d, t2)
            r2 = d1*d1 + d2*d2
            IF (r2 <= 0D0) CYCLE
            w = 1D0/r2
            dphi = s%phi_p(j) - s%phi_p(i)
            a11 = a11 + w*d1*d1
            a12 = a12 + w*d1*d2
            a22 = a22 + w*d2*d2
            b1 = b1 + w*d1*dphi
            b2 = b2 + w*d2*dphi
         END DO

         det = a11*a22 - a12*a12
         IF (a11*a22 <= 0D0 .OR. det < DEGENERATE_DET*a11*a22) THEN
            n_fallback = n_fallback + 1
            CYCLE
         END IF
         g1 = (a22*b1 - a12*b2)/det
         g2 = (a11*b2 - a12*b1)/det
         s%grad_p(:, i) = s%grad_p(:, i) + t1*g1 + t2*g2
      END DO
   END SUBROUTINE velocity_differencing

END MODULE SurfaceVelocity
