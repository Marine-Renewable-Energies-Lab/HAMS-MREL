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
MODULE NearFieldDrift
   ! Near-field mean drift force and moment on each body, from pressure integration over the mean wetted surface:
   ! WAMIT v7.5 eq. 15.54 to 15.57, as listed in docs/superpowers/specs/2026-11-04-second-order-waves-design.md,
   ! section 5.3. Time factor e^{-i omega t}, normals out of the fluid, moments about the body reference point xr.
   !
   ! Every product of two first-order quantities goes through pair_minus, section 5.4, so one state pair (s1, s2) gives
   ! the mean load for s1 = s2 and the pair term of bidirectional drift otherwise. Terms with products inside an integral
   ! are evaluated panel by panel; only term IV multiplies integrated quantities, because G and H are themselves factors.

   USE SecondOrderTypes, ONLY: panel_mesh, first_order_state, surface_geom
   USE FirstOrderState, ONLY: incident_phi, wave_gravity
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: pair_minus, waterline_point, waterline_extrapolate, near_field_pair, elevation, force_and_moment
   PUBLIC :: T_WATERLINE, T_FLARE, T_VELOCITY, T_DISPLACEMENT, T_ROTATION, T_HYDROSTATIC

   ! Term index of near_field_pair
   INTEGER, PARAMETER :: T_WATERLINE = 1, T_FLARE = 2, T_VELOCITY = 3, T_DISPLACEMENT = 4, T_ROTATION = 5, T_HYDROSTATIC = 6

   COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)

CONTAINS

   PURE COMPLEX(8) FUNCTION pair_minus(f1, g1, f2, g2)
      ! Difference-frequency coefficient of the product f g for components 1 and 2. For f1 = f2 and g1 = g2 it is the
      ! time average of Re(f e^{-i omega t}) Re(g e^{-i omega t}), which is Re(f g*)/2.
      COMPLEX(8), INTENT(IN) :: f1, g1, f2, g2

      pair_minus = 0.25D0*(f1*CONJG(g2) + CONJG(f2)*g1)
   END FUNCTION pair_minus

   SUBROUTINE waterline_extrapolate(mesh, geom, s, phi_wl)
      ! Total potential at the midpoint of each waterline segment, at z = 0: the incident wave exactly, plus the
      ! perturbation potential extrapolated from the adjacent panel with its gradient.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(first_order_state), INTENT(IN) :: s
      COMPLEX(8), INTENT(OUT) :: phi_wl(:)
      REAL(8) :: x(3), g_wave
      INTEGER :: iw, ip

      g_wave = wave_gravity(s%omega, s%k, s%h)
      DO iw = 1, geom%nwl
         ip = geom%wl_panel(iw)
         x = waterline_point(geom, iw)
         phi_wl(iw) = incident_phi(x, s%omega, s%k, s%h, s%amp, s%beta, g_wave) &
                      + s%phi_p(ip) + SUM(s%grad_p(:, ip)*(x - mesh%xc(:, ip)))
      END DO
   END SUBROUTINE waterline_extrapolate

   SUBROUTINE near_field_pair(mesh, geom, s1, s2, phi_wl1, phi_wl2, terms)
      ! Near-field load of the state pair (s1, s2), which may be the same state.
      !   phi_wl1, phi_wl2   total potential at the waterline segments (nwl)
      !   terms(6, 6, nbody) (component, term, body). Components: force, then moment about xr of the body. Terms: 1 I
      !                      waterline, 2 Ia flare, 3 II velocity, 4 III displacement, 5 IV rotation, 6 V hydrostatic.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      COMPLEX(8), INTENT(IN) :: phi_wl1(:), phi_wl2(:)
      COMPLEX(8), INTENT(OUT) :: terms(:,:,:)

      terms = (0D0, 0D0)
      CALL add_waterline_term(mesh, geom, s1, s2, phi_wl1, phi_wl2, terms)
      CALL add_flare_term(mesh, geom, s1, s2, phi_wl1, phi_wl2, terms)
      CALL add_velocity_term(mesh, s1, s2, terms)
      CALL add_displacement_term(mesh, s1, s2, terms)
      CALL add_rotation_term(mesh, s1, s2, terms)
      CALL add_hydrostatic_term(mesh, geom, s1, s2, terms)
   END SUBROUTINE near_field_pair

   SUBROUTINE add_waterline_term(mesh, geom, s1, s2, phi_wl1, phi_wl2, terms)
      ! I: (rho g / 2) int r' [ P(zeta, zeta) - 2 P(zeta, X3) ] dl,   r' = n' or x x n'
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      COMPLEX(8), INTENT(IN) :: phi_wl1(:), phi_wl2(:)
      COMPLEX(8), INTENT(INOUT) :: terms(:,:,:)
      REAL(8) :: x(3)
      COMPLEX(8) :: z1, z2, x3_1, x3_2, bracket
      INTEGER :: iw, b

      DO iw = 1, geom%nwl
         b = geom%wl_body(iw)
         x = waterline_point(geom, iw) - mesh%xr(:, b)
         z1 = elevation(s1, phi_wl1(iw))
         z2 = elevation(s2, phi_wl2(iw))
         x3_1 = vertical_motion(s1%xi(:, b), x)
         x3_2 = vertical_motion(s2%xi(:, b), x)
         bracket = pair_minus(z1, z1, z2, z2) - 2D0*pair_minus(z1, x3_1, z2, x3_2)
         terms(:, T_WATERLINE, b) = terms(:, T_WATERLINE, b) &
                                    + 0.5D0*s1%rho*s1%g*bracket*geom%wl_len(iw)*force_and_moment(x, geom%wl_nh(:, iw))
      END DO
   END SUBROUTINE add_waterline_term

   SUBROUTINE add_flare_term(mesh, geom, s1, s2, phi_wl1, phi_wl2, terms)
      ! Ia: (rho g / 2) int r'' [n_z / sqrt(1 - n_z^2)] P(zeta - X3, zeta - X3) dl,   r'' = k or x x k
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      COMPLEX(8), INTENT(IN) :: phi_wl1(:), phi_wl2(:)
      COMPLEX(8), INTENT(INOUT) :: terms(:,:,:)
      REAL(8), PARAMETER :: K_HAT(3) = [0D0, 0D0, 1D0]
      REAL(8) :: x(3), flare
      COMPLEX(8) :: r1, r2
      INTEGER :: iw, b

      DO iw = 1, geom%nwl
         b = geom%wl_body(iw)
         x = waterline_point(geom, iw) - mesh%xr(:, b)
         r1 = elevation(s1, phi_wl1(iw)) - vertical_motion(s1%xi(:, b), x)
         r2 = elevation(s2, phi_wl2(iw)) - vertical_motion(s2%xi(:, b), x)
         flare = geom%wl_nz(iw)/SQRT(1D0 - geom%wl_nz(iw)**2)
         terms(:, T_FLARE, b) = terms(:, T_FLARE, b) + 0.5D0*s1%rho*s1%g*flare*pair_minus(r1, r1, r2, r2) &
                                *geom%wl_len(iw)*force_and_moment(x, K_HAT)
      END DO
   END SUBROUTINE add_flare_term

   SUBROUTINE add_velocity_term(mesh, s1, s2, terms)
      ! II: -rho int r (1/2) sum_i P(d_i phi, d_i phi) dS,   r = n or x x n
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      COMPLEX(8), INTENT(INOUT) :: terms(:,:,:)
      COMPLEX(8) :: g1(3), g2(3), q
      INTEGER :: ip, b, i

      DO ip = 1, mesh%n
         b = mesh%body(ip)
         g1 = s1%grad_i(:, ip) + s1%grad_p(:, ip)
         g2 = s2%grad_i(:, ip) + s2%grad_p(:, ip)
         q = (0D0, 0D0)
         DO i = 1, 3
            q = q + pair_minus(g1(i), g1(i), g2(i), g2(i))
         END DO
         terms(:, T_VELOCITY, b) = terms(:, T_VELOCITY, b) - s1%rho*0.5D0*q*mesh%area(ip) &
                                   *force_and_moment(mesh%xc(:, ip) - mesh%xr(:, b), mesh%nrm(:, ip))
      END DO
   END SUBROUTINE add_velocity_term

   SUBROUTINE add_displacement_term(mesh, s1, s2, terms)
      ! III: -rho int r sum_i P(X_i, d_i phi_t) dS,   X = xi + alpha x x,   phi_t amplitude -i omega phi
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      COMPLEX(8), INTENT(INOUT) :: terms(:,:,:)
      REAL(8) :: x(3)
      COMPLEX(8) :: x1(3), x2(3), gt1(3), gt2(3), q
      INTEGER :: ip, b, i

      DO ip = 1, mesh%n
         b = mesh%body(ip)
         x = mesh%xc(:, ip) - mesh%xr(:, b)
         x1 = rigid_displacement(s1%xi(:, b), x)
         x2 = rigid_displacement(s2%xi(:, b), x)
         gt1 = -CI*s1%omega*(s1%grad_i(:, ip) + s1%grad_p(:, ip))
         gt2 = -CI*s2%omega*(s2%grad_i(:, ip) + s2%grad_p(:, ip))
         q = (0D0, 0D0)
         DO i = 1, 3
            q = q + pair_minus(x1(i), gt1(i), x2(i), gt2(i))
         END DO
         terms(:, T_DISPLACEMENT, b) = terms(:, T_DISPLACEMENT, b) - s1%rho*q*mesh%area(ip) &
                                       *force_and_moment(x, mesh%nrm(:, ip))
      END DO
   END SUBROUTINE add_displacement_term

   SUBROUTINE add_rotation_term(mesh, s1, s2, terms)
      ! IV force: alpha x G. IV moment: alpha x H + xi x G. Each cross product goes through P. On the panels of each
      ! body: G = -rho int n phi_t dS, H = -rho int (x x n) phi_t dS.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      COMPLEX(8), INTENT(INOUT) :: terms(:,:,:)
      COMPLEX(8) :: load1(6, mesh%nbody), load2(6, mesh%nbody)
      INTEGER :: b

      CALL dynamic_pressure_load(mesh, s1, load1)
      CALL dynamic_pressure_load(mesh, s2, load2)
      DO b = 1, mesh%nbody
         terms(1:3, T_ROTATION, b) = pair_cross(s1%xi(4:6, b), load1(1:3, b), s2%xi(4:6, b), load2(1:3, b))
         terms(4:6, T_ROTATION, b) = pair_cross(s1%xi(4:6, b), load1(4:6, b), s2%xi(4:6, b), load2(4:6, b)) &
                                     + pair_cross(s1%xi(1:3, b), load1(1:3, b), s2%xi(1:3, b), load2(1:3, b))
      END DO
   END SUBROUTINE add_rotation_term

   SUBROUTINE dynamic_pressure_load(mesh, s, load)
      ! load(1:3, b) = G and load(4:6, b) = H of body b.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(first_order_state), INTENT(IN) :: s
      COMPLEX(8), INTENT(OUT) :: load(:,:)
      INTEGER :: ip, b

      load = (0D0, 0D0)
      DO ip = 1, mesh%n
         b = mesh%body(ip)
         load(:, b) = load(:, b) + s%rho*CI*s%omega*(s%phi_i(ip) + s%phi_p(ip))*mesh%area(ip) &
                      *force_and_moment(mesh%xc(:, ip) - mesh%xr(:, b), mesh%nrm(:, ip))
      END DO
   END SUBROUTINE dynamic_pressure_load

   SUBROUTINE add_hydrostatic_term(mesh, geom, s1, s2, terms)
      ! V: F_S and M_S of eq. 15.56 and 15.57, with waterplane and volume properties about xr and Z_o = xr(3).
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      COMPLEX(8), INTENT(INOUT) :: terms(:,:,:)
      INTEGER, PARAMETER :: X1 = 1, X2 = 2, X3 = 3, A1 = 4, A2 = 5, A3 = 6
      REAL(8) :: rg, a, axf, ayf, zo, l1m2, l12, vxb, vyb
      INTEGER :: b

      rg = s1%rho*s1%g
      DO b = 1, mesh%nbody
         a = geom%awp(b)
         axf = a*geom%xf(b)
         ayf = a*geom%yf(b)
         zo = mesh%xr(3, b)
         l1m2 = geom%l11(b) - geom%l22(b)
         l12 = geom%l12(b)
         vxb = geom%vol(b)*geom%xb(b)
         vyb = geom%vol(b)*geom%yb(b)

         terms(3, T_HYDROSTATIC, b) = -rg*(p(A1, A3)*axf + p(A2, A3)*ayf + 0.5D0*(p(A1, A1) + p(A2, A2))*a*zo)
         terms(4, T_HYDROSTATIC, b) = rg*((p(X2, A2) - p(X3, A3))*axf - p(X2, A1)*ayf + p(A1, A2)*axf*zo &
                                          - (1.5D0*p(A1, A1) + 0.5D0*p(A2, A2))*ayf*zo - p(X2, X3)*a - p(X3, A1)*a*zo &
                                          + p(A2, A3)*l1m2 - 2D0*p(A1, A3)*l12 + p(A1, A2)*vxb &
                                          - 0.5D0*(p(A1, A1) + p(A3, A3))*vyb)
         terms(5, T_HYDROSTATIC, b) = rg*(-p(X1, A2)*axf + (p(X1, A1) - p(X3, A3))*ayf &
                                          + (0.5D0*p(A1, A1) + 1.5D0*p(A2, A2))*axf*zo - p(A1, A2)*ayf*zo &
                                          + p(X1, X3)*a - p(X3, A2)*a*zo &
                                          + p(A1, A3)*l1m2 + 2D0*p(A2, A3)*l12 + 0.5D0*(p(A2, A2) + p(A3, A3))*vxb)
      END DO

   CONTAINS

      PURE COMPLEX(8) FUNCTION p(i, j)
         INTEGER, INTENT(IN) :: i, j

         p = pair_minus(s1%xi(i, b), s1%xi(j, b), s2%xi(i, b), s2%xi(j, b))
      END FUNCTION p

   END SUBROUTINE add_hydrostatic_term

   PURE FUNCTION waterline_point(geom, iw) RESULT(x)
      ! Midpoint of a waterline segment, exactly on z = 0.
      TYPE(surface_geom), INTENT(IN) :: geom
      INTEGER, INTENT(IN) :: iw
      REAL(8) :: x(3)

      x = [geom%wl_xm(1, iw), geom%wl_xm(2, iw), 0D0]
   END FUNCTION waterline_point

   PURE COMPLEX(8) FUNCTION elevation(s, phi)
      ! Wave elevation zeta = (i omega / g) phi from the total potential at z = 0.
      TYPE(first_order_state), INTENT(IN) :: s
      COMPLEX(8), INTENT(IN) :: phi

      elevation = CI*s%omega/s%g*phi
   END FUNCTION elevation

   PURE FUNCTION force_and_moment(x, d) RESULT(r)
      ! The six weights of a force along d at arm x: d, then x x d.
      REAL(8), INTENT(IN) :: x(3), d(3)
      REAL(8) :: r(6)

      r(1:3) = d
      r(4) = x(2)*d(3) - x(3)*d(2)
      r(5) = x(3)*d(1) - x(1)*d(3)
      r(6) = x(1)*d(2) - x(2)*d(1)
   END FUNCTION force_and_moment

   PURE COMPLEX(8) FUNCTION vertical_motion(q, x)
      ! X3 = xi_3 + alpha_1 y - alpha_2 x, for x relative to xr.
      COMPLEX(8), INTENT(IN) :: q(6)
      REAL(8), INTENT(IN) :: x(3)

      vertical_motion = q(3) + q(4)*x(2) - q(5)*x(1)
   END FUNCTION vertical_motion

   PURE FUNCTION rigid_displacement(q, x) RESULT(d)
      ! X = xi + alpha x x, for x relative to xr.
      COMPLEX(8), INTENT(IN) :: q(6)
      REAL(8), INTENT(IN) :: x(3)
      COMPLEX(8) :: d(3)

      d(1) = q(1) + q(5)*x(3) - q(6)*x(2)
      d(2) = q(2) + q(6)*x(1) - q(4)*x(3)
      d(3) = q(3) + q(4)*x(2) - q(5)*x(1)
   END FUNCTION rigid_displacement

   PURE FUNCTION pair_cross(a1, b1, a2, b2) RESULT(c)
      ! sum_jk eps_ijk P(a_j, b_k)
      COMPLEX(8), INTENT(IN) :: a1(3), b1(3), a2(3), b2(3)
      COMPLEX(8) :: c(3)

      c(1) = pair_minus(a1(2), b1(3), a2(2), b2(3)) - pair_minus(a1(3), b1(2), a2(3), b2(2))
      c(2) = pair_minus(a1(3), b1(1), a2(3), b2(1)) - pair_minus(a1(1), b1(3), a2(1), b2(3))
      c(3) = pair_minus(a1(1), b1(2), a2(1), b2(2)) - pair_minus(a1(2), b1(1), a2(2), b2(1))
   END FUNCTION pair_cross

END MODULE NearFieldDrift
