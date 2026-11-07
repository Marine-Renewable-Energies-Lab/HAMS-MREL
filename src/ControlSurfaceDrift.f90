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
   ! body, on the line Gamma_c where C meets z = 0, on the free surface F between the hull waterline Gamma and Gamma_c, and on
   ! Gamma itself: Chen, X.-B. (2007), J. Eng. Math. 59:61-82, eqs. (37) and (38). Phase A: F_x, F_y and M_z, from C and
   ! Gamma_c only, Chen (2007) eqs. (39) to (41). Phase B: F_z, M_x and M_y, which also take the Gamma, F and hydrostatic terms.
   ! Those terms give nothing to F_x, F_y and M_z (Chen 2007, p. 69), so Phase B leaves the Phase A rows unchanged.
   !
   ! Conventions, as in Chen: the normal n of C points into the fluid domain, towards the body; n-bar on Gamma_c is horizontal
   ! and points into the free surface F inside C. Time factor e^{-i omega t}. Moments are about xr of the body, a space-fixed
   ! point. Quadrature as in WAMIT's low-order method: the centroid value times the panel area on C and F, the segment midpoint
   ! value times the segment length on Gamma_c and Gamma. Every product of two first-order quantities goes through pair_minus,
   ! so one state pair gives the mean load for s1 = s2 and the pair term of bidirectional drift otherwise.
   !
   ! The hydrostatic term (F_S, M_S). Chen takes it from Molin and Hairault (1983), which was not available. Lee, C.-H. (2006),
   ! On the evaluation of quadratic forces on stationary bodies (WAMIT Inc.), derives the control-surface loads by momentum,
   ! eqs. (9) and (10), and states that they equal Chen's. His F_S and M_S, eqs. (16) and (21), equal WAMIT v7.5 eqs. 15.56 and
   ! 15.57 term by term [derived], which is what NearFieldDrift's add_hydrostatic_term evaluates. Lee assumes the hull is
   ! wall-sided at Gamma. For a flared hull, Chen's near-field waterline term, eqs. (17) and (18), is eta (2 zeta_3 - eta) n',
   ! while NearFieldDrift's terms I and Ia (WAMIT v7.5 eq. 15.54) square the relative elevation, (zeta - X_3)^2, in the flare
   ! term. The hull integrals of the two near-field forms are the same, so their totals agree only if Chen's F_S holds the
   ! difference, (rho g / 2) int_Gamma P(X_3, X_3) n_z / sqrt(1 - n_z^2) r'' dl with n out of the fluid [derived]. So
   ! CS_T_HYDROSTATIC is term V plus that flare part. The flare part is derived here from Chen (2007) eqs. (17), (18) and
   ! WAMIT v7.5 eq. 15.54, and is untested: it is zero on a wall-sided hull, and no validation hull is flared at z = 0.

   USE SecondOrderTypes, ONLY: first_order_state, cs_points, cs_fields, panel_mesh, surface_geom, CS_WALL, CS_LINE, CS_FREE
   USE NearFieldDrift, ONLY: pair_minus, elevation, force_and_moment, waterline_point, vertical_motion, rigid_displacement, &
                             add_hydrostatic_term, T_HYDROSTATIC
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: control_surface_pair
   PUBLIC :: CS_T_LINE, CS_T_FLUX, CS_T_PRESSURE, CS_T_WATERLINE, CS_T_FREE, CS_T_HYDROSTATIC

   ! Term index of control_surface_pair
   INTEGER, PARAMETER :: CS_T_LINE = 1, CS_T_FLUX = 2, CS_T_PRESSURE = 3, CS_T_WATERLINE = 4, CS_T_FREE = 5, &
                         CS_T_HYDROSTATIC = 6

   COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
   REAL(8), PARAMETER :: K_HAT(3) = [0D0, 0D0, 1D0]

CONTAINS

   SUBROUTINE control_surface_pair(pts, xr, s1, s2, f1, f2, terms, mesh, geom, phi_wl1, phi_wl2)
      ! Control-surface load of the state pair (s1, s2), which may be the same state.
      !   pts     points of C (CS_WALL), Gamma_c (CS_LINE) and F (CS_FREE) of all bodies
      !   xr      (3, nbody) reference point of each body, about which the moments are taken
      !   f1, f2  total phi and grad phi of s1 and s2 at the points; phi is read on CS_LINE and CS_FREE points only
      !   mesh, geom, phi_wl1, phi_wl2
      !           optional, all four or none: the hull, and the total potential of s1 and s2 at the midpoints of its
      !           waterline segments. They give the Gamma and hydrostatic terms, which are zero without them.
      !   terms   (6, 6, nbody) (component, term, body). Components: force, then moment about xr. Terms:
      !           CS_T_LINE         (rho g / 2) int_Gamma_c P(eta, eta) r' dl,              r' = n-bar or x x n-bar
      !           CS_T_FLUX         rho int_C r(P(phi_n, grad phi)) dS,                     r(v) = v or x x v
      !           CS_T_PRESSURE     -(rho / 2) int_C sum_i P(d_i phi, d_i phi) r(n) dS
      !           CS_T_WATERLINE    -(rho g / 2) int_Gamma <Xi> r'' dl,                      r'' = k or x x k
      !           CS_T_FREE         rho int_F [P(eta, Phi_zt) + (1/2) sum_i P(d_i phi, d_i phi)] r'' dS
      !           CS_T_HYDROSTATIC  (F_S, M_S), see the module header
      !           with Xi = (eta^2 - 2 eta zeta_3) n_3 / sqrt(1 - n_3^2) - 2 eta (X . n-bar) on Gamma (Chen 2007, p. 67).
      !           The last three fill F_z, M_x and M_y only.
      TYPE(cs_points), INTENT(IN) :: pts
      REAL(8), INTENT(IN) :: xr(:, :)
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      TYPE(cs_fields), INTENT(IN) :: f1, f2
      COMPLEX(8), INTENT(OUT) :: terms(:, :, :)
      TYPE(panel_mesh), INTENT(IN), OPTIONAL :: mesh
      TYPE(surface_geom), INTENT(IN), OPTIONAL :: geom
      COMPLEX(8), INTENT(IN), OPTIONAL :: phi_wl1(:), phi_wl2(:)
      REAL(8), PARAMETER :: E(3, 3) = RESHAPE([1D0, 0D0, 0D0, 0D0, 1D0, 0D0, 0D0, 0D0, 1D0], [3, 3])
      REAL(8) :: x(3), n(3), rho, rg
      COMPLEX(8) :: z1, z2, g1(3), g2(3), gn1, gn2, q, flux(6), zt1, zt2
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
         CASE (CS_FREE)
            ! Phi_zt has the amplitude -i omega dphi/dz. On z = 0 the fields give dphi/dz by the free-surface condition, nu
            ! phi_P for the perturbation and the exact incident gradient for the incident wave (FieldVelocity), so this is
            ! -i omega nu phi wherever the two parts have the same nu. The condition dphi/dz = nu phi_P on z = 0 holds exactly
            ! for cs_velocity = green (the fields come from the Green function, which satisfies it); with cs_velocity = source
            ! it holds only to the accuracy of the Green function.
            z1 = elevation(s1, f1%phi(ip))
            z2 = elevation(s2, f2%phi(ip))
            g1 = f1%grad(:, ip)
            g2 = f2%grad(:, ip)
            zt1 = -CI*s1%omega*g1(3)
            zt2 = -CI*s2%omega*g2(3)
            q = pair_minus(z1, zt1, z2, zt2)
            DO i = 1, 3
               q = q + 0.5D0*pair_minus(g1(i), g1(i), g2(i), g2(i))
            END DO
            terms(:, CS_T_FREE, b) = terms(:, CS_T_FREE, b) + rho*q*pts%w(ip)*force_and_moment(x, K_HAT)
         END SELECT
      END DO
      IF (PRESENT(mesh) .AND. PRESENT(geom) .AND. PRESENT(phi_wl1) .AND. PRESENT(phi_wl2)) THEN
         CALL add_hull_terms(mesh, geom, xr, s1, s2, phi_wl1, phi_wl2, terms)
      END IF
   END SUBROUTINE control_surface_pair

   SUBROUTINE add_hull_terms(mesh, geom, xr, s1, s2, phi_wl1, phi_wl2, terms)
      ! The Gamma term CS_T_WATERLINE and the hydrostatic term CS_T_HYDROSTATIC. Chen's normals point into the fluid,
      ! NearFieldDrift's (wl_nz, wl_nh) out of it, so n_3 = -wl_nz and n-bar = -wl_nh. X is the body displacement at the
      ! waterline point, xi + alpha x x, and zeta_3 its vertical part. add_hydrostatic_term takes its moments about mesh%xr,
      ! so xr must be mesh%xr.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      REAL(8), INTENT(IN) :: xr(:, :)
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      COMPLEX(8), INTENT(IN) :: phi_wl1(:), phi_wl2(:)
      COMPLEX(8), INTENT(INOUT) :: terms(:, :, :)
      COMPLEX(8) :: hydro(6, 6, SIZE(terms, 3)), z1, z2, x3_1, x3_2, xn1, xn2, xi_mean
      REAL(8) :: x(3), nbar(3), n3, flare, rg
      INTEGER :: iw, b

      rg = s1%rho*s1%g
      DO iw = 1, geom%nwl
         b = geom%wl_body(iw)
         x = waterline_point(geom, iw) - xr(:, b)
         nbar = -geom%wl_nh(:, iw)
         n3 = -geom%wl_nz(iw)
         flare = n3/SQRT(1D0 - n3**2)
         z1 = elevation(s1, phi_wl1(iw))
         z2 = elevation(s2, phi_wl2(iw))
         x3_1 = vertical_motion(s1%xi(:, b), x)
         x3_2 = vertical_motion(s2%xi(:, b), x)
         xn1 = SUM(rigid_displacement(s1%xi(:, b), x)*nbar)
         xn2 = SUM(rigid_displacement(s2%xi(:, b), x)*nbar)
         xi_mean = (pair_minus(z1, z1, z2, z2) - 2D0*pair_minus(z1, x3_1, z2, x3_2))*flare - 2D0*pair_minus(z1, xn1, z2, xn2)
         terms(:, CS_T_WATERLINE, b) = terms(:, CS_T_WATERLINE, b) &
                                       - 0.5D0*rg*xi_mean*geom%wl_len(iw)*force_and_moment(x, K_HAT)
         ! The flare part of (F_S, M_S), derived from Chen (2007) eqs. (17)-(18), untested (zero on wall-sided hulls):
         ! (rho g / 2) P(X_3, X_3) wl_nz / sqrt(1 - wl_nz^2), and that factor is -flare
         terms(:, CS_T_HYDROSTATIC, b) = terms(:, CS_T_HYDROSTATIC, b) &
                                         - 0.5D0*rg*flare*pair_minus(x3_1, x3_1, x3_2, x3_2)*geom%wl_len(iw) &
                                         *force_and_moment(x, K_HAT)
      END DO
      hydro = (0D0, 0D0)
      CALL add_hydrostatic_term(mesh, geom, s1, s2, hydro)
      terms(:, CS_T_HYDROSTATIC, :) = terms(:, CS_T_HYDROSTATIC, :) + hydro(:, T_HYDROSTATIC, :)
   END SUBROUTINE add_hull_terms

END MODULE ControlSurfaceDrift
