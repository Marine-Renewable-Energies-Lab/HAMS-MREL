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
MODULE QTFPotential
   ! Body-forcing terms of the potential part of the difference-frequency QTF (indirect method: Molin 1979). The
   ! second-order body condition on the mean surface (Ogilvie 1983) is integrated against the radiation potentials psi_k at
   ! the difference frequency, one column per mode k = 1..6 nbody. psi_k satisfies dpsi_k/dn = n_k with the generalized
   ! normals, so moments come from psi_4..6 and no moment arm enters. Time factor e^{-i omega t}, normals out of the fluid.
   ! Every product of two first-order quantities goes through pair_product or pair_cross.
   !
   ! Second-derivative term B2, int_S psi n.[(X.grad)u] dS, with X = xi + alpha x (x - xr) and u the total velocity.
   ! For irrotational u and rigid X, curl(X x u) = alpha x u - (X.grad)u, so Stokes' theorem on S, bounded by the
   ! waterline Gamma, gives
   !
   !    int_S psi n.[(X.grad)u] dS = int_S [ psi n.(alpha x u) + (n x grad psi).(X x u) ] dS - oint_Gamma psi (X x u).t dl
   !
   ! with t the counter-clockwise direction of the waterline segments (body on the left), as Stokes' theorem orients it for
   ! n into the body; the unit test stokes_identity_smooth fixes this sign against the direct form. The identity holds for
   ! each part of the kernel pair, since each is bilinear in a rigid field and an irrotational, divergence-free one. Only
   ! first derivatives of phi and the tangential gradient of psi enter. On a sharp edge the direct integrand is not
   ! integrable (u ~ r^{-1/3}); this form is (r^{-2/3}). Only this form is implemented.
   !
   ! u at the waterline: without u1_wl, u2_wl it is extrapolated along the panel from its centroid value. With a
   ! differenced u the waterline row is fitted one-sided, which limits the line term to O(h) (measured in
   ! stokes_identity_differenced_velocity). potential_pair therefore passes u at the waterline from waterline_velocity
   ! (ruling R7): the incident wave exactly; for the disturbance, u_z from the free-surface condition, the normal part from
   ! the body condition and the tangential part by differencing phi along the waterline. That is O(h^2)
   ! (stokes_identity_waterline_velocity).
   !
   ! Control-surface form of the surface piece (qtf_b2 = control). On a sharp edge grad psi and u are both singular
   ! (rho^{-1/3}), so the surface piece int_S (n x grad psi).(X x u) dS, a rho^{-2/3} product integrated with one point per
   ! panel and a differenced grad psi, converges slowly. For a rigid field A = a + alpha x (x - xr) and harmonic phi and psi
   ! (v = grad phi), T_ab = v_a d_b psi + v_b d_a psi - delta_ab (v.grad psi) is divergence-free, and div(T A) = T : grad A = 0
   ! because grad A is antisymmetric: the two-potential, rigid-weighted extension of the Gauss identity of Chen (2007),
   ! J. Eng. Math. 59:61-82, Appendix C (eq. 63-65), as Cong, Teng, Chen and Gou (2020), Ocean Eng. 204, 107278 (eq. 35-39)
   ! use it for the second-order body condition. With (n x grad psi).(A x v) = (n.A)(v.grad psi) - (n.v)(A.grad psi) and
   ! Gauss on the fluid D between the hull S, a closed control surface C around the body (wall and bottom) and the free
   ! surface F between the waterline and C:
   !
   !    int_S (n x grad psi).(A x v) dS = int_S (A.v) psi_n dS
   !                                      + int_{C u F} [ (A.v) psi_no + (A.grad psi) v_no - (A.n_o)(v.grad psi) ] dS
   !
   ! n on S is the hull normal (out of the fluid, so out of D) and psi_n = n.grad psi = n_k (generalized normals) on the
   ! body of mode k, 0 on the others. n_o on C and F points out of D, the opposite of cs_points%nrm, which points into D
   ! (towards the body on C, (0, 0, -1) on F). No line integral enters. The hull term needs no grad psi and is singular
   ! like u only (rho^{-1/3}); grad psi on C and F is smooth. Under the kernel each term is bilinear in (A, v), so the
   ! per-point tensor M_ab = pair_product(A_a, v_b) of the two states carries every product, and psi enters linearly:
   !    surface piece       n^T M g - g^T M n                    (g = grad psi)
   !    hull term           psi_n tr M                           (b2_body_term)
   !    C and F term        psi_no tr M + g^T M n_o - n_o^T M g  (b2_control_term)
   ! One rigid A per control volume: the C and F of body b enclose body b alone and take its motion, so C must have a
   ! bottom (a C that reaches the sea bed leaves the sea bed between body and C as a boundary of D, not implemented).
   ! rot and line of the Stokes form stay on the hull.
   !
   ! Sharp edges (qtf_edge_weights = on). The hull term is linear in u, which goes like rho^p near an edge (p = -1/3 at a
   ! 90 deg edge; SharpEdge states the corner-flow result). Its one-point rule then misses the singular part on the panels
   ! next to the edge, so b2_body_term takes the panel weights of SharpEdge (edge_w) that make the rule exact for rho^p.
   ! The rot piece psi n.(alpha x u) is linear in u too but is not weighted: on the sharp-cylinder fixture against the
   ! exact reference the weight did not reduce its error consistently (F_z 27% to 26%, F_x 9.4% to 9.6%).
   !
   ! Assembly (potential_pair). At Omega = omega_1 - omega_2, per unit amplitudes, with the free-surface integral omitted
   ! (Phase 6):
   !
   !    T_pot,k = i Omega rho [ int_S phi_I2 n_k dS + int_S psi_k (B - dphi_I2/dn) dS ]
   !    B = (Xdot - u).(alpha x n)  [B1]  - n.[(X.grad)u]  [B2]  + n.d/dt(H x), H x = (1/2) alpha x (alpha x x)  [B3]
   !
   ! with X = xi + alpha x (x - xr), Xdot = -i omega X, d/dt = -i Omega on the pair, phi_I2 of SecondOrderIncident. The
   ! terms are kept apart in the order of the TP_ indices. n_k of the FK term are the generalized normals mesh%nj of the
   ! panels of the body of mode k.

   USE SecondOrderTypes, ONLY: panel_mesh, surface_geom, first_order_state, cs_points, CS_WALL, CS_FREE
   USE NearFieldDrift, ONLY: pair_cross, pair_product, rigid_displacement, waterline_point, KERNEL_DIFF
   USE SurfaceVelocity, ONLY: tangential_fit
   USE FirstOrderState, ONLY: incident_phi, incident_grad, wave_gravity
   USE SecondOrderIncident, ONLY: incident2, incident2_pair, incident2_phi, incident2_grad
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hessian_term_converted, psi_waterline, potential_pair, psi_surface_gradient, waterline_velocity
   PUBLIC :: b2_body_term, b2_control_term, b2_control_check

   ! Terms of potential_pair: Froude-Krylov of phi_I2, incident-diffraction (-dphi_I2/dn), B1, B2, B3
   INTEGER, PARAMETER, PUBLIC :: NPOT_TERMS = 5
   INTEGER, PARAMETER, PUBLIC :: TP_FK = 1, TP_DIFFRACTION = 2, TP_B1 = 3, TP_B2 = 4, TP_B3 = 5

   COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)

CONTAINS

   SUBROUTINE potential_pair(mesh, geom, s1, s2, wl1, wl2, psi, grad_psi, psi_wl, omega_psi, rho, g, terms, ok, &
                             phi2_constant, cs, grad_psi_cs, u1_cs, u2_cs, edge_w)
      ! Body-forcing terms of the difference-frequency potential QTF of the state pair (s1, s2), per unit amplitudes.
      !   wl1, wl2      total potential at the waterline segments (nwl), as for near_field_pair
      !   psi           (n, 6 nbody) radiation potentials at |Omega|, Omega = s1%omega - s2%omega
      !   grad_psi      (3, n, 6 nbody) their gradients (psi_surface_gradient)
      !   psi_wl        (nwl, 6 nbody) psi at the waterline (psi_waterline with the same grad_psi)
      !   omega_psi     the frequency psi was solved at; must be |Omega| within 1e-6 max(omega_1, omega_2)
      !   rho           water density
      !   g             gravity of the incident waves, for phi_I2 (SecondOrderIncident)
      !   terms         (6 nbody, NPOT_TERMS), T_pot,k per term
      !   ok            optional: .false. instead of a stop when psi does not belong to Omega (terms are then zero)
      !   phi2_constant optional test hook: phi_I2 replaced by this constant (dphi_I2/dn = 0)
      !   cs, grad_psi_cs, u1_cs, u2_cs   optional, all or none (qtf_b2 = control): the points of C and F of every body,
      !                 grad psi_k there at |Omega| (3, cs%n, 6 nbody), and the total velocity of s1 and s2 there (3, cs%n).
      !                 B2 then takes the surface piece from b2_body_term and b2_control_term; rot and line are unchanged.
      !   edge_w        optional (n), with cs only: panel weights of the hull term at sharp edges (SharpEdge)
      ! Omega = 0 gives zeros. For Omega < 0, psi(-|Omega|) = conj psi(|Omega|) is used (the radiation potential of the
      ! negative frequency), so terms(2, 1) = conj terms(1, 2). Exact equality s1%omega == s2%omega is the diagonal: the
      ! caller passes the same grid value for both.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      COMPLEX(8), INTENT(IN) :: wl1(:), wl2(:), psi(:, :), grad_psi(:, :, :), psi_wl(:, :)
      REAL(8), INTENT(IN) :: omega_psi, rho, g
      COMPLEX(8), INTENT(OUT) :: terms(:, :)
      LOGICAL, INTENT(OUT), OPTIONAL :: ok
      COMPLEX(8), INTENT(IN), OPTIONAL :: phi2_constant
      TYPE(cs_points), INTENT(IN), OPTIONAL :: cs
      COMPLEX(8), INTENT(IN), OPTIONAL :: grad_psi_cs(:, :, :), u1_cs(:, :), u2_cs(:, :)
      REAL(8), INTENT(IN), OPTIONAL :: edge_w(:)

      INTEGER :: nmode
      LOGICAL :: control
      REAL(8) :: omega
      CHARACTER(LEN=200) :: msg

      terms = (0D0, 0D0)
      IF (PRESENT(ok)) ok = .TRUE.
      omega = s1%omega - s2%omega
      IF (omega == 0D0) RETURN

      nmode = 6*mesh%nbody
      IF (SIZE(psi, 1) /= mesh%n .OR. SIZE(psi, 2) /= nmode .OR. SIZE(grad_psi, 1) /= 3 .OR. &
          SIZE(grad_psi, 2) /= mesh%n .OR. SIZE(grad_psi, 3) /= nmode .OR. SIZE(psi_wl, 1) /= geom%nwl .OR. &
          SIZE(psi_wl, 2) /= nmode .OR. SIZE(terms, 1) /= nmode .OR. SIZE(terms, 2) /= NPOT_TERMS .OR. &
          SIZE(wl1) /= geom%nwl .OR. SIZE(wl2) /= geom%nwl) THEN
         WRITE (*, '(a)') ' potential_pair: psi, grad_psi, psi_wl, wl or terms does not have the size of the mesh.'
         ERROR STOP 1
      END IF
      control = PRESENT(cs)
      IF ((PRESENT(grad_psi_cs) .NEQV. control) .OR. (PRESENT(u1_cs) .NEQV. control) .OR. &
          (PRESENT(u2_cs) .NEQV. control)) THEN
         WRITE (*, '(a)') ' potential_pair: give cs, grad_psi_cs, u1_cs and u2_cs together.'
         ERROR STOP 1
      END IF
      IF (control) THEN
         IF (SIZE(grad_psi_cs, 1) /= 3 .OR. SIZE(grad_psi_cs, 2) /= cs%n .OR. SIZE(grad_psi_cs, 3) /= nmode .OR. &
             SIZE(u1_cs, 1) /= 3 .OR. SIZE(u1_cs, 2) /= cs%n .OR. SIZE(u2_cs, 1) /= 3 .OR. SIZE(u2_cs, 2) /= cs%n) THEN
            WRITE (*, '(a)') ' potential_pair: grad_psi_cs, u1_cs or u2_cs does not have the size of the control surface.'
            ERROR STOP 1
         END IF
      END IF
      IF (PRESENT(edge_w)) THEN
         IF (.NOT. control .OR. SIZE(edge_w) /= mesh%n) THEN
            WRITE (*, '(a)') ' potential_pair: edge_w needs the control surface and one weight per panel.'
            ERROR STOP 1
         END IF
      END IF
      IF (ABS(omega_psi - ABS(omega)) > 1D-6*MAX(s1%omega, s2%omega)) THEN
         WRITE (msg, '(a,ES14.6,a,ES14.6,a)') ' potential_pair: psi was solved at omega =', omega_psi, &
            ' rad/s, but the pair needs Omega =', ABS(omega), ' rad/s.'
         CALL refuse()
         RETURN
      END IF
      IF (ALL(psi == (0D0, 0D0))) THEN
         WRITE (msg, '(a,ES14.6,a)') ' potential_pair: psi at Omega =', ABS(omega), &
            ' rad/s is zero: it was not stored.'
         CALL refuse()
         RETURN
      END IF

      IF (omega > 0D0) THEN
         IF (control) THEN
            CALL assemble(psi, grad_psi, psi_wl, grad_psi_cs)
         ELSE
            CALL assemble(psi, grad_psi, psi_wl)
         END IF
      ELSE
         IF (control) THEN
            CALL assemble(CONJG(psi), CONJG(grad_psi), CONJG(psi_wl), CONJG(grad_psi_cs))
         ELSE
            CALL assemble(CONJG(psi), CONJG(grad_psi), CONJG(psi_wl))
         END IF
      END IF

   CONTAINS

      SUBROUTINE refuse()
         IF (PRESENT(ok)) THEN
            ok = .FALSE.
            RETURN
         END IF
         WRITE (*, '(a)') TRIM(msg)
         ERROR STOP 1
      END SUBROUTINE refuse

      SUBROUTINE assemble(ps, gps, pw, gcs)
         ! The terms with ps = psi at Omega; gcs, grad psi on C and F at Omega, for qtf_b2 = control
         COMPLEX(8), INTENT(IN) :: ps(:, :), gps(:, :, :), pw(:, :)
         COMPLEX(8), INTENT(IN), OPTIONAL :: gcs(:, :, :)
         TYPE(incident2) :: inc
         COMPLEX(8), ALLOCATABLE :: u1(:, :), u2(:, :), uw1(:, :), uw2(:, :), hs(:), pieces(:, :), hb(:), hc(:)
         COMPLEX(8) :: p2, gr2(3), dn2, d1(3), d2(3), v1(3), v2(3), r1(3), r2(3), hx(3), b1v, b3v
         COMPLEX(8) :: zero3(3)
         REAL(8) :: nv(3), x(3)
         INTEGER :: i, b, j, k

         zero3 = (0D0, 0D0)
         IF (.NOT. PRESENT(phi2_constant)) THEN
            CALL incident2_pair(s1%omega, s1%k, s1%beta, s2%omega, s2%k, s2%beta, s1%h, g, inc)
         END IF
         ALLOCATE (u1(3, mesh%n), u2(3, mesh%n), hs(nmode))
         u1 = s1%grad_i + s1%grad_p
         u2 = s2%grad_i + s2%grad_p

         DO i = 1, mesh%n
            b = mesh%body(i)
            nv = mesh%nrm(:, i)
            x = mesh%xc(:, i) - mesh%xr(:, b)
            IF (PRESENT(phi2_constant)) THEN
               p2 = phi2_constant
               dn2 = (0D0, 0D0)
            ELSE
               p2 = incident2_phi(inc, mesh%xc(:, i))
               CALL incident2_grad(inc, mesh%xc(:, i), gr2)
               dn2 = SUM(nv*gr2)
            END IF
            DO j = 1, 6
               terms(6*(b - 1) + j, TP_FK) = terms(6*(b - 1) + j, TP_FK) + p2*mesh%nj(j, i)*mesh%area(i)
            END DO

            ! B1 = (Xdot - u).(alpha x n) = n.[(Xdot - u) x alpha]
            d1 = rigid_displacement(s1%xi(:, b), x)
            d2 = rigid_displacement(s2%xi(:, b), x)
            v1 = -CI*s1%omega*d1 - u1(:, i)
            v2 = -CI*s2%omega*d2 - u2(:, i)
            b1v = SUM(nv*pair_cross(KERNEL_DIFF, v1, s1%xi(4:6, b), v2, s2%xi(4:6, b)))

            ! B3 = n.d/dt(H x), H x = (1/2) alpha x (alpha x x)
            r1 = rigid_displacement([zero3, s1%xi(4:6, b)], x)
            r2 = rigid_displacement([zero3, s2%xi(4:6, b)], x)
            hx = 0.5D0*pair_cross(KERNEL_DIFF, s1%xi(4:6, b), r1, s2%xi(4:6, b), r2)
            b3v = -CI*omega*SUM(nv*hx)

            DO k = 1, nmode
               terms(k, TP_DIFFRACTION) = terms(k, TP_DIFFRACTION) - ps(i, k)*dn2*mesh%area(i)
               terms(k, TP_B1) = terms(k, TP_B1) + ps(i, k)*b1v*mesh%area(i)
               terms(k, TP_B3) = terms(k, TP_B3) + ps(i, k)*b3v*mesh%area(i)
            END DO
         END DO

         ! B2 = -n.[(X.grad)u], with u at the waterline from the free-surface and body conditions
         ALLOCATE (pieces(nmode, 3))
         IF (geom%nwl > 0) THEN
            ALLOCATE (uw1(3, geom%nwl), uw2(3, geom%nwl))
            CALL waterline_velocity(mesh, geom, s1, wl1, body_normal_velocity(s1), uw1)
            CALL waterline_velocity(mesh, geom, s2, wl2, body_normal_velocity(s2), uw2)
            CALL hessian_term_converted(mesh, geom, ps, gps, pw, s1%xi, u1, s2%xi, u2, KERNEL_DIFF, hs, uw1, uw2, &
                                        pieces=pieces)
         ELSE
            CALL hessian_term_converted(mesh, geom, ps, gps, pw, s1%xi, u1, s2%xi, u2, KERNEL_DIFF, hs, pieces=pieces)
         END IF
         IF (PRESENT(gcs)) THEN
            ! qtf_b2 = control: rot and line as above, the surface piece from the hull term and the C and F term
            ALLOCATE (hb(nmode), hc(nmode))
            CALL b2_body_term(mesh, gps, s1%xi, u1, s2%xi, u2, KERNEL_DIFF, hb, edge_w)   ! edge_w may be absent
            CALL b2_control_term(cs, mesh%xr, gcs, s1%xi, u1_cs, s2%xi, u2_cs, KERNEL_DIFF, hc)
            hs = pieces(:, 1) + pieces(:, 3) + hb + hc
         END IF
         terms(:, TP_B2) = -hs

         terms = CI*omega*rho*terms
      END SUBROUTINE assemble

      FUNCTION body_normal_velocity(s) RESULT(vn)
         ! Normal velocity of the body at the waterline segment midpoints, -i omega X.n with the normal of the panel
         TYPE(first_order_state), INTENT(IN) :: s
         COMPLEX(8) :: vn(geom%nwl)
         INTEGER :: iw, b

         DO iw = 1, geom%nwl
            b = geom%wl_body(iw)
            vn(iw) = -CI*s%omega*SUM(mesh%nrm(:, geom%wl_panel(iw)) &
                                     *rigid_displacement(s%xi(:, b), waterline_point(geom, iw) - mesh%xr(:, b)))
         END DO
      END FUNCTION body_normal_velocity

   END SUBROUTINE potential_pair

   SUBROUTINE psi_surface_gradient(mesh, geom, psi, grad_psi)
      ! grad psi_k at the panel centroids, as velocity_differencing builds grad phi: the normal part dpsi_k/dn = n_k (mesh%nj)
      ! on the panels of the body of mode k and 0 on the others, the tangential part from tangential_fit of the panel
      ! values. A panel without a fit keeps the normal part only.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      COMPLEX(8), INTENT(IN) :: psi(:, :)
      COMPLEX(8), INTENT(OUT) :: grad_psi(:, :, :)
      INTEGER :: i, k, kb, j
      REAL(8) :: t1(3), t2(3), dn
      COMPLEX(8) :: g1, g2
      LOGICAL :: ok

      DO k = 1, SIZE(psi, 2)
         kb = (k - 1)/6 + 1
         j = k - 6*(kb - 1)
         DO i = 1, mesh%n
            dn = 0D0
            IF (mesh%body(i) == kb) dn = mesh%nj(j, i)
            grad_psi(:, i, k) = mesh%nrm(:, i)*dn
            CALL tangential_fit(mesh, geom, i, psi(:, k), t1, t2, g1, g2, ok)
            IF (ok) grad_psi(:, i, k) = grad_psi(:, i, k) + t1*g1 + t2*g2
         END DO
      END DO
   END SUBROUTINE psi_surface_gradient

   SUBROUTINE waterline_velocity(mesh, geom, s, phi_wl, vn_wl, u_wl)
      ! Total first-order velocity at the waterline segment midpoints (z = 0), without the one-sided fit of the waterline
      ! row (ruling R7). The incident wave is exact. The disturbance phi_p = phi_wl - phi_I has, in the frame of the
      ! segment direction t, the horizontal normal nh (wl_nh) and the vertical e_z:
      !   u_z = nu phi_p                                  free-surface condition -omega^2 phi + g phi_z = 0, nu = s%nu
      !   u_t = d(phi_p)/ds along the waterline           3-point differencing of the midpoint values of each loop
      !   u_nh = (vn - n.u_I - n_z u_z) / (n.nh)          body condition n.u = vn, with n the panel normal (n.t = 0)
      ! phi_wl from waterline_extrapolate is O(h^2), and so are these. Neighbours on the loop are geom%wl_prev, wl_next. A
      ! segment with one neighbour takes the one-sided difference, one with none the tangential part of the panel velocity.
      !   phi_wl   (nwl) total potential at the midpoints        vn_wl  (nwl) normal velocity of the body there
      !   u_wl     (3, nwl)
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(first_order_state), INTENT(IN) :: s
      COMPLEX(8), INTENT(IN) :: phi_wl(:), vn_wl(:)
      COMPLEX(8), INTENT(OUT) :: u_wl(:, :)
      COMPLEX(8), ALLOCATABLE :: php(:)
      COMPLEX(8) :: ui(3), upz, upt, upn
      REAL(8) :: g_wave, x(3), t(3), nv(3), nh(3), h1, h2
      INTEGER :: iw, ip, jp, jn

      IF (geom%nwl == 0) RETURN
      IF (.NOT. ALLOCATED(geom%wl_next)) THEN
         WRITE (*, '(a)') ' waterline_velocity: geom has no waterline loop neighbours (build it with build_surface_geom).'
         ERROR STOP 1
      END IF
      g_wave = wave_gravity(s%omega, s%k, s%h)
      ALLOCATE (php(geom%nwl))
      DO iw = 1, geom%nwl
         php(iw) = phi_wl(iw) - incident_phi(waterline_point(geom, iw), s%omega, s%k, s%h, s%amp, s%beta, g_wave)
      END DO

      DO iw = 1, geom%nwl
         ip = geom%wl_panel(iw)
         x = waterline_point(geom, iw)
         t = (geom%wl_x2(:, iw) - geom%wl_x1(:, iw))/geom%wl_len(iw)
         nv = mesh%nrm(:, ip)
         nh = geom%wl_nh(:, iw)
         CALL incident_grad(x, s%omega, s%k, s%h, s%amp, s%beta, g_wave, ui)
         upz = s%nu*php(iw)

         jp = geom%wl_prev(iw)
         jn = geom%wl_next(iw)
         h1 = 0D0
         h2 = 0D0
         IF (jp > 0) h1 = -SUM((waterline_point(geom, jp) - x)*t)
         IF (jn > 0) h2 = SUM((waterline_point(geom, jn) - x)*t)
         IF (h1 > 0D0 .AND. h2 > 0D0) THEN
            upt = -h2/(h1*(h1 + h2))*php(jp) + (h2 - h1)/(h1*h2)*php(iw) + h1/(h2*(h1 + h2))*php(jn)
         ELSE IF (h2 > 0D0) THEN
            upt = (php(jn) - php(iw))/h2
         ELSE IF (h1 > 0D0) THEN
            upt = (php(iw) - php(jp))/h1
         ELSE
            upt = SUM(s%grad_p(:, ip)*t)
         END IF

         upn = (vn_wl(iw) - SUM(nv*ui) - nv(3)*upz)/SUM(nv*nh)
         u_wl(:, iw) = ui + t*upt + nh*upn
         u_wl(3, iw) = u_wl(3, iw) + upz
      END DO
   END SUBROUTINE waterline_velocity

   SUBROUTINE hessian_term_converted(mesh, geom, psi, grad_psi, psi_wl, x1, u1, x2, u2, kernel, out, u1_wl, u2_wl, &
                                     pieces)
      ! Kernel pair (1, 2) of int_S psi_k n.[(X.grad)u] dS for each mode column k, in the converted form above.
      !   psi        (n, nmode) radiation potentials at the panel centroids, nmode = 6 nbody
      !   grad_psi   (3, n, nmode) their gradients; only the tangential part enters
      !   psi_wl     (nwl, nmode) psi at the waterline segments (psi_waterline)
      !   x1, x2     (6, nbody) rigid motions, translations then rotations about xr of the body
      !   u1, u2     (3, n) total first-order velocities at the panel centroids
      !   kernel     KERNEL_DIFF or KERNEL_SUM of NearFieldDrift
      !   out        (nmode)
      !   u1_wl, u2_wl  optional (3, nwl): u at the waterline segment midpoints (waterline_velocity), in place of the
      !              extrapolation from the panel centroid. Both or neither.
      !   pieces     optional (nmode, 3): the three pieces of out apart, rot int psi n.(alpha x u), surface
      !              int (n x grad psi).(X x u) and line -oint psi (X x u).t; out is their sum (to rounding)
      ! The sum runs over the panels and waterline segments of all bodies: psi_k of one body acts on every surface, with
      ! the motion of the body that owns the panel.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      COMPLEX(8), INTENT(IN) :: psi(:, :), grad_psi(:, :, :), psi_wl(:, :)
      COMPLEX(8), INTENT(IN) :: x1(:, :), u1(:, :), x2(:, :), u2(:, :)
      INTEGER, INTENT(IN) :: kernel
      COMPLEX(8), INTENT(OUT) :: out(:)
      COMPLEX(8), INTENT(IN), OPTIONAL :: u1_wl(:, :), u2_wl(:, :)
      COMPLEX(8), INTENT(OUT), OPTIONAL :: pieces(:, :)

      INTEGER :: i, k, b, iw, ip
      REAL(8) :: nv(3), x(3), t(3)
      COMPLEX(8) :: d1(3), d2(3), rot(3), xu(3), uw1(3), uw2(3)
      COMPLEX(8), ALLOCATABLE :: uc1(:, :), uc2(:, :)

      out = (0D0, 0D0)
      IF (PRESENT(pieces)) pieces = (0D0, 0D0)
      DO i = 1, mesh%n
         b = mesh%body(i)
         nv = mesh%nrm(:, i)
         d1 = rigid_displacement(x1(:, b), mesh%xc(:, i) - mesh%xr(:, b))
         d2 = rigid_displacement(x2(:, b), mesh%xc(:, i) - mesh%xr(:, b))
         rot = pair_cross(kernel, x1(4:6, b), u1(:, i), x2(4:6, b), u2(:, i))   ! alpha x u
         xu = pair_cross(kernel, d1, u1(:, i), d2, u2(:, i))                   ! X x u
         DO k = 1, SIZE(out)
            out(k) = out(k) + mesh%area(i)*(psi(i, k)*SUM(nv*rot) + SUM(cross_rc(nv, grad_psi(:, i, k))*xu))
         END DO
         IF (PRESENT(pieces)) THEN
            DO k = 1, SIZE(out)
               pieces(k, 1) = pieces(k, 1) + mesh%area(i)*psi(i, k)*SUM(nv*rot)
               pieces(k, 2) = pieces(k, 2) + mesh%area(i)*SUM(cross_rc(nv, grad_psi(:, i, k))*xu)
            END DO
         END IF
      END DO

      IF (geom%nwl == 0) RETURN
      IF (PRESENT(u1_wl) .NEQV. PRESENT(u2_wl)) THEN
         WRITE (*, '(a)') ' hessian_term_converted: give u1_wl and u2_wl together.'
         ERROR STOP 1
      END IF
      IF (.NOT. PRESENT(u1_wl)) THEN
         uc1 = TRANSPOSE(u1)   ! (n, 3): each component a contiguous panel scalar for tangential_fit
         uc2 = TRANSPOSE(u2)
      END IF
      DO iw = 1, geom%nwl
         ip = geom%wl_panel(iw)
         b = geom%wl_body(iw)
         x = waterline_point(geom, iw)
         t = (geom%wl_x2(:, iw) - geom%wl_x1(:, iw))/geom%wl_len(iw)
         d1 = rigid_displacement(x1(:, b), x - mesh%xr(:, b))
         d2 = rigid_displacement(x2(:, b), x - mesh%xr(:, b))
         IF (PRESENT(u1_wl)) THEN
            uw1 = u1_wl(:, iw)
            uw2 = u2_wl(:, iw)
         ELSE
            uw1 = extrapolated_velocity(uc1, ip, x)
            uw2 = extrapolated_velocity(uc2, ip, x)
         END IF
         xu = pair_cross(kernel, d1, uw1, d2, uw2)
         DO k = 1, SIZE(out)
            out(k) = out(k) - psi_wl(iw, k)*SUM(xu*t)*geom%wl_len(iw)
         END DO
         IF (PRESENT(pieces)) pieces(:, 3) = pieces(:, 3) - psi_wl(iw, :)*SUM(xu*t)*geom%wl_len(iw)
      END DO

   CONTAINS

      FUNCTION extrapolated_velocity(uc, ip, xw) RESULT(uw)
         ! u at the waterline point xw, extrapolated from the centroid of panel ip along the panel. xw lies in the panel
         ! plane, so the tangential fit of each component is all the extrapolation needs. The centroid value alone would
         ! leave an error of order the panel height in the line integral. A panel without a fit keeps its centroid value.
         COMPLEX(8), INTENT(IN) :: uc(:, :)
         INTEGER, INTENT(IN) :: ip
         REAL(8), INTENT(IN) :: xw(3)
         COMPLEX(8) :: uw(3)
         REAL(8) :: t1(3), t2(3), dx(3)
         COMPLEX(8) :: g1, g2
         LOGICAL :: ok
         INTEGER :: m

         dx = xw - mesh%xc(:, ip)
         DO m = 1, 3
            CALL tangential_fit(mesh, geom, ip, uc(:, m), t1, t2, g1, g2, ok)
            uw(m) = uc(ip, m) + g1*DOT_PRODUCT(t1, dx) + g2*DOT_PRODUCT(t2, dx)
         END DO
      END FUNCTION extrapolated_velocity

   END SUBROUTINE hessian_term_converted

   SUBROUTINE b2_body_term(mesh, grad_psi, x1, u1, x2, u2, kernel, out, edge_w)
      ! Hull term of the control-surface form of the B2 surface piece: out(k) = int_S (X.u) psi_n dS, kernel pair (1, 2),
      ! with psi_n = n.grad_psi(:, i, k) (n_k on the body of mode k when grad_psi comes from psi_surface_gradient) and the
      ! motion of the body that owns the panel. Arguments as hessian_term_converted.
      !   edge_w   optional (n): panel weights for the velocity singularity at sharp edges (SharpEdge, qtf_edge_weights);
      !            the panel area is multiplied by them. Without it the plain one-point rule.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      COMPLEX(8), INTENT(IN) :: grad_psi(:, :, :), x1(:, :), u1(:, :), x2(:, :), u2(:, :)
      INTEGER, INTENT(IN) :: kernel
      COMPLEX(8), INTENT(OUT) :: out(:)
      REAL(8), INTENT(IN), OPTIONAL :: edge_w(:)
      INTEGER :: i, b, k
      COMPLEX(8) :: d1(3), d2(3), m(3, 3), tr
      REAL(8) :: nv(3), w

      out = (0D0, 0D0)
      DO i = 1, mesh%n
         b = mesh%body(i)
         nv = mesh%nrm(:, i)
         d1 = rigid_displacement(x1(:, b), mesh%xc(:, i) - mesh%xr(:, b))
         d2 = rigid_displacement(x2(:, b), mesh%xc(:, i) - mesh%xr(:, b))
         m = pair_tensor(kernel, d1, u1(:, i), d2, u2(:, i))
         tr = m(1, 1) + m(2, 2) + m(3, 3)
         w = mesh%area(i)
         IF (PRESENT(edge_w)) w = w*edge_w(i)
         DO k = 1, SIZE(out)
            out(k) = out(k) + w*SUM(nv*grad_psi(:, i, k))*tr
         END DO
      END DO
   END SUBROUTINE b2_body_term

   SUBROUTINE b2_control_term(pts, xr, grad_psi_cs, x1, u1_cs, x2, u2_cs, kernel, out)
      ! C and F term of the control-surface form of the B2 surface piece, kernel pair (1, 2):
      !    out(k) = sum over the CS_WALL and CS_FREE points of w [ (A.v) psi_no + (A.grad psi) v_no - (A.n_o)(v.grad psi) ]
      ! with n_o = -pts%nrm (out of the fluid between hull, C and F), A the rigid displacement of the body of the point
      ! (about xr(:, b)), v = u1_cs, u2_cs and grad psi = grad_psi_cs(:, p, k). CS_LINE points are not used.
      !   xr           (3, nbody) reference points of the bodies, as mesh%xr
      !   grad_psi_cs  (3, pts%n, nmode)       u1_cs, u2_cs  (3, pts%n)       x1, x2  (6, nbody)
      TYPE(cs_points), INTENT(IN) :: pts
      REAL(8), INTENT(IN) :: xr(:, :)
      COMPLEX(8), INTENT(IN) :: grad_psi_cs(:, :, :), x1(:, :), u1_cs(:, :), x2(:, :), u2_cs(:, :)
      INTEGER, INTENT(IN) :: kernel
      COMPLEX(8), INTENT(OUT) :: out(:)
      INTEGER :: p, b, k
      COMPLEX(8) :: d1(3), d2(3), m(3, 3), tr, gr(3)
      REAL(8) :: no(3)

      out = (0D0, 0D0)
      DO p = 1, pts%n
         IF (pts%kind(p) /= CS_WALL .AND. pts%kind(p) /= CS_FREE) CYCLE
         b = pts%body(p)
         no = -pts%nrm(:, p)
         d1 = rigid_displacement(x1(:, b), pts%x(:, p) - xr(:, b))
         d2 = rigid_displacement(x2(:, b), pts%x(:, p) - xr(:, b))
         m = pair_tensor(kernel, d1, u1_cs(:, p), d2, u2_cs(:, p))
         tr = m(1, 1) + m(2, 2) + m(3, 3)
         DO k = 1, SIZE(out)
            gr = grad_psi_cs(:, p, k)
            ! (A.v) psi_no + (A.grad psi) v_no - (A.n_o)(v.grad psi) = psi_no tr M + g^T M n_o - n_o^T M g
            out(k) = out(k) + pts%w(p)*(SUM(no*gr)*tr + SUM(gr*MATMUL(m, no)) - SUM(no*MATMUL(m, gr)))
         END DO
      END DO
   END SUBROUTINE b2_control_term

   PURE FUNCTION pair_tensor(kernel, a1, v1, a2, v2) RESULT(m)
      ! m(a, b) = pair_product(kernel, a1(a), v1(b), a2(a), v2(b)): every bilinear product of a rigid field and a velocity
      INTEGER, INTENT(IN) :: kernel
      COMPLEX(8), INTENT(IN) :: a1(3), v1(3), a2(3), v2(3)
      COMPLEX(8) :: m(3, 3)
      INTEGER :: ia, ib

      DO ib = 1, 3
         DO ia = 1, 3
            m(ia, ib) = pair_product(kernel, a1(ia), v1(ib), a2(ia), v2(ib))
         END DO
      END DO
   END FUNCTION pair_tensor

   SUBROUTINE b2_control_check(nbody, pts, ok, msg)
      ! qtf_b2 = control needs, for every body, F (CS_FREE points) and a C with a bottom (a CS_WALL point whose normal points
      ! up, towards the body). A C that reaches the sea bed leaves the sea bed between body and C as a boundary of the
      ! fluid, whose term is not implemented. ok = .false. with the reason in msg otherwise.
      ! The bottom test is one CS_WALL point with n_z > 0.5. That is enough for a generated C, whose bottom is a whole disc
      ! or rectangle when it exists. A C read from a file (cs_shape = file) with only part of a bottom also passes, and its
      ! gap is then an open boundary of the fluid that the identity does not account for. Nothing checks that a C from a
      ! file is closed.
      INTEGER, INTENT(IN) :: nbody
      TYPE(cs_points), INTENT(IN) :: pts
      LOGICAL, INTENT(OUT) :: ok
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: msg
      CHARACTER(LEN=300) :: line
      INTEGER :: b

      ok = .TRUE.
      msg = ''
      DO b = 1, nbody
         IF (.NOT. ANY(pts%kind == CS_FREE .AND. pts%body == b)) THEN
            WRITE (line, '(a,I0,a)') 'qtf_b2 = control: body ', b, ' has no free surface F inside its control surface'// &
               ' (see the warning above). Give FSMesh_<b>.pnl, change cs_shape or cs_margin, or set qtf_b2 = body.'
            ok = .FALSE.
            msg = TRIM(line)
            RETURN
         END IF
         IF (.NOT. ANY(pts%kind == CS_WALL .AND. pts%body == b .AND. pts%nrm(3, :) > 0.5D0)) THEN
            WRITE (line, '(a,I0,a)') 'qtf_b2 = control: the control surface of body ', b, ' has no bottom (it reaches the'// &
               ' sea bed, or the file has none). For a generated C, set cs_bottom_margin below h minus the draft to keep'// &
               ' a bottom, or set qtf_b2 = body.'
            ok = .FALSE.
            msg = TRIM(line)
            RETURN
         END IF
      END DO
   END SUBROUTINE b2_control_check

   SUBROUTINE psi_waterline(mesh, geom, psi, grad_psi, psi_wl)
      ! psi at the midpoint of each waterline segment, at z = 0, extrapolated from the adjacent panel with its gradient, as
      ! waterline_extrapolate does for the disturbance potential.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      COMPLEX(8), INTENT(IN) :: psi(:, :), grad_psi(:, :, :)
      COMPLEX(8), INTENT(OUT) :: psi_wl(:, :)
      INTEGER :: iw, ip, k
      REAL(8) :: dx(3)

      DO iw = 1, geom%nwl
         ip = geom%wl_panel(iw)
         dx = waterline_point(geom, iw) - mesh%xc(:, ip)
         DO k = 1, SIZE(psi, 2)
            psi_wl(iw, k) = psi(ip, k) + SUM(grad_psi(:, ip, k)*dx)
         END DO
      END DO
   END SUBROUTINE psi_waterline

   PURE FUNCTION cross_rc(a, b) RESULT(c)
      ! a x b for real a and complex b
      REAL(8), INTENT(IN) :: a(3)
      COMPLEX(8), INTENT(IN) :: b(3)
      COMPLEX(8) :: c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   END FUNCTION cross_rc

END MODULE QTFPotential
