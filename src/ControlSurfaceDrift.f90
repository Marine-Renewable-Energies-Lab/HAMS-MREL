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
MODULE ControlSurfaceDrift
   ! Control-surface (middle-field) mean drift of each body, from the total first-order field on a closed surface C around the
   ! body and on the line Gamma_c where C meets z = 0: Chen, X.-B. (2007), J. Eng. Math. 59:61-82, eq. (39) to (41), the C and
   ! Gamma_c parts of eq. (37) and (38). Design spec docs/superpowers/specs/2026-11-06-control-surface-drift-design.md,
   ! section 3.2. Phase A: F_x, F_y and M_z.
   !
   ! Conventions, as in Chen: the normal n of C points into the fluid domain, towards the body; n-bar on Gamma_c is horizontal
   ! and points into the free surface F inside C. Time factor e^{-i omega t}. Moments are about xr of the body, a space-fixed
   ! point. Quadrature as in WAMIT's low-order method (spec section 3.5): the centroid value times the panel area on C, the
   ! segment midpoint value times the segment length on Gamma_c. Every product of two first-order quantities goes through
   ! pair_minus, so one state pair gives the mean load for s1 = s2 and the pair term of bidirectional drift otherwise.
   !
   ! The C and Gamma_c terms fill all six components, but only I = 1, 2, 6 are complete: the hull waterline, F and
   ! hydrostatic terms of F_z, M_x and M_y come in Phase B (spec section 3.3).

   USE SecondOrderTypes, ONLY: first_order_state, cs_points, cs_fields, CS_WALL, CS_LINE
   USE NearFieldDrift, ONLY: pair_minus, elevation, force_and_moment
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: control_surface_pair
   PUBLIC :: CS_T_LINE, CS_T_FLUX, CS_T_PRESSURE, CS_T_WATERLINE, CS_T_FREE, CS_T_HYDROSTATIC

   ! Term index of control_surface_pair. Phase A fills the first three.
   INTEGER, PARAMETER :: CS_T_LINE = 1, CS_T_FLUX = 2, CS_T_PRESSURE = 3, CS_T_WATERLINE = 4, CS_T_FREE = 5, &
                         CS_T_HYDROSTATIC = 6

CONTAINS

   SUBROUTINE control_surface_pair(pts, xr, s1, s2, f1, f2, terms)
      ! Control-surface load of the state pair (s1, s2), which may be the same state.
      !   pts     points of C (CS_WALL) and Gamma_c (CS_LINE) of all bodies; CS_FREE points are not used in Phase A
      !   xr      (3, nbody) reference point of each body, about which the moments are taken
      !   f1, f2  total phi and grad phi of s1 and s2 at the points; phi is read on CS_LINE points only
      !   terms   (6, 6, nbody) (component, term, body). Components: force, then moment about xr. Terms:
      !           CS_T_LINE      (rho g / 2) int_Gamma_c P(eta, eta) r' dl,               r' = n-bar or x x n-bar
      !           CS_T_FLUX      rho int_C r(P(phi_n, grad phi)) dS,                      r(v) = v or x x v
      !           CS_T_PRESSURE  -(rho / 2) int_C sum_i P(d_i phi, d_i phi) r(n) dS
      TYPE(cs_points), INTENT(IN) :: pts
      REAL(8), INTENT(IN) :: xr(:, :)
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      TYPE(cs_fields), INTENT(IN) :: f1, f2
      COMPLEX(8), INTENT(OUT) :: terms(:, :, :)
      REAL(8), PARAMETER :: E(3, 3) = RESHAPE([1D0, 0D0, 0D0, 0D0, 1D0, 0D0, 0D0, 0D0, 1D0], [3, 3])
      REAL(8) :: x(3), n(3), rho, rg
      COMPLEX(8) :: z1, z2, g1(3), g2(3), gn1, gn2, q, flux(6)
      INTEGER :: ip, b, i

      terms = (0D0, 0D0)
      rho = s1%rho
      rg = s1%rho*s1%g
      DO ip = 1, pts%n
         b = pts%body(ip)
         x = pts%x(:, ip) - xr(:, b)
         n = pts%nrm(:, ip)
         SELECT CASE (pts%kind(ip))
         CASE (CS_LINE)
            z1 = elevation(s1, f1%phi(ip))
            z2 = elevation(s2, f2%phi(ip))
            terms(:, CS_T_LINE, b) = terms(:, CS_T_LINE, b) &
                                     + 0.5D0*rg*pair_minus(z1, z1, z2, z2)*pts%w(ip)*force_and_moment(x, n)
         CASE (CS_WALL)
            g1 = f1%grad(:, ip)
            g2 = f2%grad(:, ip)
            gn1 = SUM(n*g1)
            gn2 = SUM(n*g2)
            q = (0D0, 0D0)
            flux = (0D0, 0D0)
            DO i = 1, 3
               q = q + pair_minus(g1(i), g1(i), g2(i), g2(i))
               ! force_and_moment is linear in its direction, so the complex vector goes in by components
               flux = flux + pair_minus(gn1, g1(i), gn2, g2(i))*force_and_moment(x, E(:, i))
            END DO
            terms(:, CS_T_FLUX, b) = terms(:, CS_T_FLUX, b) + rho*pts%w(ip)*flux
            terms(:, CS_T_PRESSURE, b) = terms(:, CS_T_PRESSURE, b) - 0.5D0*rho*q*pts%w(ip)*force_and_moment(x, n)
         END SELECT
      END DO
   END SUBROUTINE control_surface_pair

END MODULE ControlSurfaceDrift
