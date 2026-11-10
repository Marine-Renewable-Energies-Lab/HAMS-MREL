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
MODULE QTFField
   ! Second-order pressure and free-surface elevation of a pair of first-order waves at a field point: the quadratic part
   ! (products of first-order quantities) plus the second-order incident potential. The second-order diffracted and radiated
   ! potentials are omitted, so near a body the result is partial.
   !
   ! Bernoulli's equation expanded about z = 0 (Dean and Dalrymple 1991), time factor e^{-i omega t}:
   !   pressure at a fixed point below z = 0:  p2 = -rho [phi2_t + 1/2 |grad phi1|^2]
   !   elevation on z = 0:  eta2 = -(1/g) [phi2_t + 1/2 |grad phi1|^2 + eta1 phi1_zt],
   !   with eta1 = (i omega/g) phi1 and phi1_zt = -i omega phi1_z.
   ! For the ordered pair (i, j) and the kernel P of NearFieldDrift (pair_minus: difference, pair_plus: sum), per unit
   ! amplitudes:
   !   quad_p   = -(rho/2) sum_m P(u_im, u_im, u_jm, u_jm)
   !   quad_eta = -(1/g) [ 1/2 sum_m P(u_im, u_im, u_jm, u_jm) + P(eta_i, -i omega_i u_i3, eta_j, -i omega_j u_j3) ]
   !   pot_p    = i rho Omega phi2,   pot_eta = (i Omega/g) phi2,
   ! with Omega = omega_i -/+ omega_j and phi2 the second-order incident potential of the pair (SecondOrderIncident). The
   ! physical e^{-i Omega t} amplitude of two waves a_i /= a_j is 2 a_i a_j (quad + pot); for the sum with i = j it is
   ! a_i^2 (quad + pot).

   USE SecondOrderIncident, ONLY: incident2, incident2_phi
   USE NearFieldDrift, ONLY: pair_product
   USE SecondOrderTypes, ONLY: panel_mesh, cs_points, CS_FREE
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: field_pair, build_field_points

CONTAINS

   PURE SUBROUTINE field_pair(kernel, w1, phi1, grad1, w2, phi2, grad2, inc, x, surface, rho, g, quad, pot)
      ! Quadratic part quad and second-order incident part pot of the pressure (surface = .FALSE.) or the elevation
      ! (surface = .TRUE., x(3) = 0) of the pair of unit waves with frequencies w1, w2, first-order potentials phi and
      ! gradients grad at x. inc is the second-order incident potential of the pair, made with the same kernel.
      INTEGER, INTENT(IN) :: kernel
      REAL(8), INTENT(IN) :: w1, w2, x(3), rho, g
      COMPLEX(8), INTENT(IN) :: phi1, grad1(3), phi2, grad2(3)
      TYPE(incident2), INTENT(IN) :: inc
      LOGICAL, INTENT(IN) :: surface
      COMPLEX(8), INTENT(OUT) :: quad, pot
      COMPLEX(8) :: speed, p2, eta1, eta2
      COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
      INTEGER :: m

      speed = (0D0, 0D0)
      DO m = 1, 3
         speed = speed + pair_product(kernel, grad1(m), grad1(m), grad2(m), grad2(m))
      END DO
      p2 = incident2_phi(inc, x)
      IF (surface) THEN
         eta1 = CI*w1/g*phi1
         eta2 = CI*w2/g*phi2
         quad = -(0.5D0*speed + pair_product(kernel, eta1, -CI*w1*grad1(3), eta2, -CI*w2*grad2(3)))/g
         pot = CI*inc%omega/g*p2
      ELSE
         quad = -0.5D0*rho*speed
         pot = CI*rho*inc%omega*p2
      END IF
   END SUBROUTINE field_pair

   SUBROUTINE build_field_points(mesh, xfp, pts, surface, ok, msg)
      ! Field points xfp(3,:) as a cs_points set (kind CS_FREE, nrm = 0, w = 0). Type as PressureElevation: |z| <= 1e-6 is
      ! a surface (elevation) point and gets z = 0, z < -1e-6 is a pressure point. psize and body are those of the nearest
      ! hull panel (centroids of the full panel set, so symmetry images are included). Rejected, checked in point order,
      ! the first failing point giving msg: above the free surface; within one local panel size of a hull panel centroid
      ! (the rule of check_control_surface); inside a body, decided by the nearest centroid c as (x - c).n > 0, with n the
      ! panel normal, which points into the body.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: xfp(:,:)
      TYPE(cs_points), INTENT(OUT) :: pts
      LOGICAL, ALLOCATABLE, INTENT(OUT) :: surface(:)
      LOGICAL, INTENT(OUT) :: ok
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: msg
      REAL(8), PARAMETER :: ZTOL = 1D-6
      REAL(8) :: x(3), d2, d2min
      INTEGER :: ip, i, imin, ihull, n
      CHARACTER(LEN=300) :: line

      n = SIZE(xfp, 2)
      pts%n = n
      ALLOCATE (pts%x(3, n), pts%nrm(3, n), pts%w(n), pts%psize(n), pts%kind(n), pts%body(n), surface(n))
      pts%x = xfp
      pts%nrm = 0D0
      pts%w = 0D0
      pts%psize = 0D0
      pts%kind = CS_FREE
      pts%body = 0
      surface = .FALSE.
      ok = .TRUE.
      msg = ''

      DO ip = 1, n
         x = xfp(:, ip)
         IF (x(3) > ZTOL) THEN
            WRITE (line, '(a,i0,a,3(1x,es11.3),a)') 'field point ', ip, ' at', x, ' is above the free surface (z > 0)'
            msg = TRIM(line)
            ok = .FALSE.
            RETURN
         END IF
         IF (ABS(x(3)) <= ZTOL) THEN
            surface(ip) = .TRUE.
            x(3) = 0D0
            pts%x(3, ip) = 0D0
         END IF

         imin = 0
         ihull = 0
         d2min = HUGE(1D0)
         DO i = 1, mesh%n
            d2 = SUM((x - mesh%xc(:, i))**2)
            IF (d2 < d2min) THEN
               d2min = d2
               imin = i
            END IF
            IF (ihull == 0 .AND. d2 < mesh%psize(i)**2) ihull = i
         END DO
         IF (imin == 0) CYCLE
         pts%psize(ip) = mesh%psize(imin)
         pts%body(ip) = mesh%body(imin)

         IF (ihull > 0) THEN
            WRITE (line, '(a,i0,a,3(1x,es11.3),a,i0,a,i0,a)') 'field point ', ip, ' at', x, ' lies on the hull of body ', &
               mesh%body(ihull), ' (within one panel size of panel ', ihull, ')'
            msg = TRIM(line)
            ok = .FALSE.
            RETURN
         END IF
         IF (SUM((x - mesh%xc(:, imin))*mesh%nrm(:, imin)) > 0D0) THEN
            WRITE (line, '(a,i0,a,3(1x,es11.3),a,i0)') 'field point ', ip, ' at', x, ' is inside body ', mesh%body(imin)
            msg = TRIM(line)
            ok = .FALSE.
            RETURN
         END IF
      END DO
   END SUBROUTINE build_field_points

END MODULE QTFField
