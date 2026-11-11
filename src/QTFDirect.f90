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
MODULE QTFDirect
   ! Direct method of the difference-frequency second-order problem (Phase 7a): the right side of the boundary integral
   ! equation for phi_D at Omega. HAMS solves 2 pi phi + int phi dG/dn dS = int G dphi/dn dS (ASSB_LEFT, ASSB_RBC), so with
   ! dphi_D/dn = b on S the right side is r_i = int_S G(x_i, xi) b(xi) dS_xi, b = B - dphi_I2/dn (QTFPotential gives B).
   ! The parts of b known at the panels go through ASSB_RHS_FIELD. This module gives the B2 part.
   !
   ! B2 = -n.[(X.grad)u] contains grad u, which is not known pointwise at a sharp edge. Ruling R7-1: the Stokes form of
   ! QTFPotential (hessian_term_converted) with psi replaced by G(x_i, .):
   !
   !    int_S G n.[(X.grad)u] dS = int_S [ G n.(alpha x u) + (n x grad_xi G).(X x u) ] dS - oint_Gamma G (X x u).t dl
   !
   ! with t the counter-clockwise direction of the waterline segments, as there. G(x_i, .) is singular at x_i on S. Its
   ! surface gradient goes like 1/rho^2 there, with a leading part odd in the direction, so the surface term is a principal
   ! value: a disc of radius eps about x_i cut out of S adds a line integral around the disc that is O(eps), so the identity
   ! keeps its form in the limit. The panel integrals give that principal value: on the own panel the 1/r part of
   ! int grad G dS is the in-plane part, finite, and its normal part (the jump) is left out (SGLINTBD_*2 skip the normal
   ! part in the plane). Only n x grad G enters, so the normal part does not matter.
   !
   ! Near panels, linear density. With one point per panel for the density c = (X x u) x n, the surface term converges
   ! like h^0.65 at collocation points on S (measured, Phase 7a Task 2): grad_xi G ~ 1/rho^2 near x_i, so the linear
   ! variation of c over the own panel and its neighbours carries an O(h) part that the centroid value drops. So on the
   ! panels that CALGREEN takes as near (within 50 panel sizes, the flag of the panel integrals) the Rankine part
   ! 1/r + 1/r1 is integrated against c linear on the panel, c(xi) = c_p + C_p (xi - xi_p): the zeroth moment as before,
   ! plus
   !    sum_ab C_p(a, b) int_p d_a(1/r) (xi - xi_p)_b dS     (and the same for 1/r1)
   ! with the first moments in closed form (rankine_moments). C_p comes from tangential_fit of the components of w = X x u
   ! over the panel neighbours, crossed with n_p, so that c stays normal to n_p and only the in-plane gradient enters. The
   ! density a = n.(alpha x u) of the term int G a dS is taken linear in the same way, with the first moment
   ! int_p (1/r)(xi - xi_p) dS. The wave part and the far panels keep the one-point rule.
   !
   ! Per panel p and field point x the integrals are int_p G dS and int_p grad_xi G dS, with grad_xi the gradient with
   ! respect to the source point xi, written as G = 1/r + 1/r1 + rest:
   !   1/r     panel integral of 1/|x - xi| and of its gradient in xi (rankine_panel_integral)
   !   1/r1    the same at the mirrored field point (x, y, -z): 1/r1 = 1/|(x, y, -z) - xi|
   !   rest    one point: the Green function with the panel centroid first, so that GRN(2:4) is the gradient in xi. G is
   !           symmetric, so this is the same function. With flag = 1 GRN holds 1/r1, which is taken out at the centroid.
   ! This is source_panel_integrals of SurfaceVelocity with the roles of the points swapped. source_panel_integrals is not
   ! used because it gives the gradient in the field point, which the source solve needs; the Stokes form needs it in the
   ! source point. For 1/r the two differ in sign only, for 1/r1 in the sign of the horizontal part, but the wave part in
   ! finite depth is not a function of z_x + z_xi alone (cosh k(z_x + h) cosh k(z_xi + h)), so its gradient in xi cannot
   ! be had from the one in x: it needs a call of the Green function with the points swapped.
   ! The line integral of G over a waterline segment takes 1/r and 1/r1 in closed form along the segment and the rest at the
   ! midpoint.
   !
   ! Symmetry images. The full panel p is image b of base panel j. With M_b the reflection, int_{M_b q} G(x, xi) dS =
   ! int_q G(M_b x, xi) dS and int_{M_b q} grad G = M_b int_q grad G(M_b x, .), as build_source_operator has it. The rows
   ! are the full panels; the right side of block ip at base panel i is sum_a RXY(ip, a) r(image a of i), the block form of
   ! ASSB_RBC.
   !
   ! Stored integrals. In production (direct_rhs_b2, direct_rhs_b2_block with stored, the default) int_p G and
   ! int_p grad_xi G come from the arrays of the left side at Omega: CGRN, RKBN and DS (CGRN_MULTI_COMB, RKBN_MULTI_COMB
   ! and DS_MULTI_COMB for several bodies), with the flag rule of CALGREEN (base-to-base distance within 50 panel sizes).
   ! CALGREEN puts the source panel first, so GRN(2:4) is grad_xi G, and IMAGE_TERM_EXACT treats 1/r1 as
   ! green_panel_integrals does. Row image a against panel image b is base row against stored image a xor b, with the
   ! gradient mirrored by M_b, so no value is computed twice and B2 uses the discrete G of the left side and of
   ! ASSB_RHS_FIELD. The stored values are single precision (GRN_KIND = 4). green_panel_integrals stays for the tests,
   ! which set no stored arrays; the two agree to the single-precision storage (rhs_b2_stored_matches_direct).
   ! The near-panel moments and the waterline line integrals are computed as before.
   !
   ! Sharp edges. The linear density on the near panels uses a fitted gradient of w = X x u and alpha x u. Next to a sharp
   ! edge u is singular (rho^p) and the fit of its gradient is not meaningful, which is what the Stokes form avoids. So a
   ! panel that shares a node with a panel of the same body whose normal differs by more than edge_deg (DEFAULT_EDGE_DEG of
   ! SharpEdge, or qtf_edge_angle) keeps the one-point rule: cl = 0 there (sharp_edge_panels). build_neighbours already
   ! leaves out neighbours across such an edge, so the fit elsewhere is over the same face.
   !
   ! The Green function is taken at the wave of the solver globals (V, WVN, NK, H): the caller sets them to Omega, and
   ! with stored integrals CALGREEN (CALGREEN_MULTI) must have been run at Omega.

   USE Const_mod, ONLY: RXY
   USE HAMS_mod, ONLY: NTHREAD
   USE WaveDyn_mod, ONLY: H, NK, V, WVN
   USE PanelMesh_mod, ONLY: ISX, ISY, DS, PNSZ, XYZ_P, DS_MULTI_COMB, PNSZ_MULTI_COMB, XYZ_GLOBAL_MULTI_COMB_P
   USE Potentials_mod, ONLY: CGRN, RKBN, CGRN_MULTI_COMB, RKBN_MULTI_COMB
   USE SurfaceGeometry, ONLY: body_sizes, merge_nodes
   USE SharpEdge, ONLY: DEFAULT_EDGE_DEG
   USE INFG3D_Open, ONLY: INFGREEN3D
   USE FinGreen3D_Open, ONLY: FINGREEN3D
   USE SecondOrderTypes, ONLY: panel_mesh, surface_geom, first_order_state
   USE FirstOrderState, ONLY: image_signs
   USE NearFieldDrift, ONLY: pair_cross, rigid_displacement, waterline_point, KERNEL_DIFF
   USE WaterlineField, ONLY: rankine_panel_integral
   USE QTFPotential, ONLY: waterline_velocity
   USE SurfaceVelocity, ONLY: tangential_fit
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: direct_rhs_b2, direct_rhs_b2_block, b2_stokes_densities, b2_stokes_at_point
   PUBLIC :: green_panel_integrals, stored_panel_integrals, green_line_integral, rankine_moments, sharp_edge_panels

   COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
   ! As CALGREEN: a panel closer than this many panel sizes gets panel integrals of 1/r and of its image
   REAL(8), PARAMETER :: NEAR_PANEL_SIZES = 50D0

CONTAINS

   SUBROUTINE direct_rhs_b2(mesh, geom, s1, s2, wl1, wl2, rhs, is_multi, stored, edge_deg)
      ! The B2 part of the right side for the pair (s1, s2), difference-frequency kernel:
      !   rhs(i, ip) = sum_a RXY(ip, a) r(x of image a of base panel i),  r(x) = -int_S G(x, .) n.[(X.grad)u] dS
      ! in the Stokes form above, with u at the waterline from waterline_velocity, as potential_pair has it.
      !   wl1, wl2   total potential at the waterline segments (nwl), as for potential_pair
      !   rhs        (nbase, nsys), the layout of one column of ASSB_RHS_FIELD
      !   is_multi   optional: the panel integrals read the combined multi-body arrays; default mesh%nbody > 1
      !   stored     optional: the panel integrals of G from the stored arrays (default .true.), or fresh
      !              (green_panel_integrals, for tests without stored arrays)
      !   edge_deg   optional: the sharp-edge angle of sharp_edge_panels; default DEFAULT_EDGE_DEG
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      COMPLEX(8), INTENT(IN) :: wl1(:), wl2(:)
      COMPLEX(8), INTENT(OUT) :: rhs(:, :)
      LOGICAL, INTENT(IN), OPTIONAL :: is_multi, stored
      REAL(8), INTENT(IN), OPTIONAL :: edge_deg
      COMPLEX(8), ALLOCATABLE :: a(:, :), c(:, :, :), cl(:, :, :, :), e(:, :), r(:, :, :)
      LOGICAL, ALLOCATABLE :: one_point(:)
      REAL(8) :: deg

      IF (SIZE(rhs, 1) /= mesh%nbase .OR. SIZE(rhs, 2) /= mesh%nsys) THEN
         WRITE (*, '(a)') ' direct_rhs_b2: rhs must be (nbase, nsys).'
         ERROR STOP 1
      END IF
      deg = DEFAULT_EDGE_DEG
      IF (PRESENT(edge_deg)) deg = edge_deg
      CALL sharp_edge_panels(mesh, deg, one_point)
      ALLOCATE (a(mesh%n, 1), c(3, mesh%n, 1), cl(3, 4, mesh%n, 1), e(geom%nwl, 1), r(mesh%nbase, 1, mesh%nsys))
      CALL b2_stokes_densities(mesh, geom, s1, s2, wl1, wl2, a(:, 1), c(:, :, 1), e(:, 1), cl(:, :, :, 1), one_point)
      CALL direct_rhs_b2_block(mesh, geom, a, c, e, r, is_multi, cl, stored)
      rhs = r(:, 1, :)
   END SUBROUTINE direct_rhs_b2

   SUBROUTINE b2_stokes_densities(mesh, geom, s1, s2, wl1, wl2, a, c, e, cl, one_point)
      ! The pair data of the Stokes form, so that for any potential psi
      !   int_S psi n.[(X.grad)u] dS = sum_p [ psi_p a_p + (int_p grad psi).c_p ] - sum_w psi_w e_w
      ! with psi_p and int_p grad psi panel integrals and psi_w the line integral over segment w:
      !   a(p) = n.(alpha x u)       c(:, p) = (X x u) x n   [(n x g).w = g.(w x n)]       e(w) = (X x u).t
      ! kernel KERNEL_DIFF, the motion of the body that owns the panel or segment.
      !   cl  optional (3, 4, n): the linear parts on each panel, from tangential_fit over the panel neighbours (zero on a
      !       panel without a fit): c(xi) = c_p + cl(:, 1:3, p) (xi - xi_p), from the components of w = X x u crossed
      !       with n_p, and a(xi) = a_p + cl(:, 4, p).(xi - xi_p), from the components of alpha x u dotted with n_p.
      !   one_point  optional (n), with cl: panels that keep cl = 0, the panels at sharp edges (sharp_edge_panels)
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      COMPLEX(8), INTENT(IN) :: wl1(:), wl2(:)
      COMPLEX(8), INTENT(OUT) :: a(:), c(:, :), e(:)
      COMPLEX(8), INTENT(OUT), OPTIONAL :: cl(:, :, :)
      LOGICAL, INTENT(IN), OPTIONAL :: one_point(:)
      COMPLEX(8), ALLOCATABLE :: uw1(:, :), uw2(:, :), w(:, :)
      COMPLEX(8) :: u1(3), u2(3), d1(3), d2(3), rot(3), xu(3), dw1(6), dw2(6), g1, g2
      REAL(8) :: nv(3), x(3), t(3), t1(3), t2(3)
      INTEGER :: i, b, iw, m
      LOGICAL :: ok

      IF (SIZE(a) /= mesh%n .OR. SIZE(c, 1) /= 3 .OR. SIZE(c, 2) /= mesh%n .OR. SIZE(e) /= geom%nwl .OR. &
          SIZE(wl1) /= geom%nwl .OR. SIZE(wl2) /= geom%nwl) THEN
         WRITE (*, '(a)') ' b2_stokes_densities: a, c, e or wl does not have the size of the mesh.'
         ERROR STOP 1
      END IF
      DO i = 1, mesh%n
         b = mesh%body(i)
         nv = mesh%nrm(:, i)
         u1 = s1%grad_i(:, i) + s1%grad_p(:, i)
         u2 = s2%grad_i(:, i) + s2%grad_p(:, i)
         d1 = rigid_displacement(s1%xi(:, b), mesh%xc(:, i) - mesh%xr(:, b))
         d2 = rigid_displacement(s2%xi(:, b), mesh%xc(:, i) - mesh%xr(:, b))
         rot = pair_cross(KERNEL_DIFF, s1%xi(4:6, b), u1, s2%xi(4:6, b), u2)
         xu = pair_cross(KERNEL_DIFF, d1, u1, d2, u2)
         a(i) = SUM(nv*rot)
         c(:, i) = cross_cr(xu, nv)
         IF (PRESENT(cl)) THEN
            IF (.NOT. ALLOCATED(w)) ALLOCATE (w(mesh%n, 6))
            w(i, 1:3) = xu
            w(i, 4:6) = rot
         END IF
      END DO

      IF (PRESENT(cl)) THEN
         IF (SIZE(cl, 1) /= 3 .OR. SIZE(cl, 2) /= 4 .OR. SIZE(cl, 3) /= mesh%n) THEN
            WRITE (*, '(a)') ' b2_stokes_densities: cl does not have the size of the mesh.'
            ERROR STOP 1
         END IF
         DO i = 1, mesh%n
            nv = mesh%nrm(:, i)
            cl(:, :, i) = (0D0, 0D0)
            IF (PRESENT(one_point)) THEN
               IF (one_point(i)) CYCLE
            END IF
            ok = .TRUE.
            DO m = 1, 6
               CALL tangential_fit(mesh, geom, i, w(:, m), t1, t2, g1, g2, ok)
               IF (.NOT. ok) EXIT
               dw1(m) = g1
               dw2(m) = g2
            END DO
            cl(:, :, i) = (0D0, 0D0)
            IF (.NOT. ok) CYCLE
            ! a = n_p.(alpha x u) on the flat panel: fit the vector alpha x u, not a, whose panel values also carry the
            ! change of the normal from panel to panel
            cl(:, 4, i) = SUM(nv*dw1(4:6))*t1 + SUM(nv*dw2(4:6))*t2
            ! c = w x n, so dc/dt_m = (dw/dt_m) x n, and C = sum_m (dc/dt_m) t_m^T
            dw1(1:3) = cross_cr(dw1(1:3), nv)
            dw2(1:3) = cross_cr(dw2(1:3), nv)
            DO m = 1, 3
               cl(:, m, i) = dw1(1:3)*t1(m) + dw2(1:3)*t2(m)
            END DO
         END DO
      END IF

      IF (geom%nwl == 0) RETURN
      ALLOCATE (uw1(3, geom%nwl), uw2(3, geom%nwl))
      CALL waterline_velocity(mesh, geom, s1, wl1, body_normal_velocity(s1), uw1)
      CALL waterline_velocity(mesh, geom, s2, wl2, body_normal_velocity(s2), uw2)
      DO iw = 1, geom%nwl
         b = geom%wl_body(iw)
         x = waterline_point(geom, iw)
         t = (geom%wl_x2(:, iw) - geom%wl_x1(:, iw))/geom%wl_len(iw)
         d1 = rigid_displacement(s1%xi(:, b), x - mesh%xr(:, b))
         d2 = rigid_displacement(s2%xi(:, b), x - mesh%xr(:, b))
         xu = pair_cross(KERNEL_DIFF, d1, uw1(:, iw), d2, uw2(:, iw))
         e(iw) = SUM(xu*t)
      END DO

   CONTAINS

      FUNCTION body_normal_velocity(s) RESULT(vn)
         ! Normal velocity of the body at the waterline segment midpoints, -i omega X.n, as in potential_pair
         TYPE(first_order_state), INTENT(IN) :: s
         COMPLEX(8) :: vn(geom%nwl)
         INTEGER :: jw, jb

         DO jw = 1, geom%nwl
            jb = geom%wl_body(jw)
            vn(jw) = -CI*s%omega*SUM(mesh%nrm(:, geom%wl_panel(jw)) &
                                     *rigid_displacement(s%xi(:, jb), waterline_point(geom, jw) - mesh%xr(:, jb)))
         END DO
      END FUNCTION body_normal_velocity

   END SUBROUTINE b2_stokes_densities

   SUBROUTINE direct_rhs_b2_block(mesh, geom, a, c, e, rhs, is_multi, cl, stored)
      ! The B2 part of the right side for nrhs pairs at once, from their densities (b2_stokes_densities):
      !   a (n, nrhs), c (3, n, nrhs), e (nwl, nrhs)  ->  rhs (nbase, nrhs, nsys), the layout of ASSB_RHS_FIELD
      ! The panel integrals do not depend on the pair, so one pass over the field points serves all of them.
      !   cl  optional (3, 4, n, nrhs): the linear parts of c and a (b2_stokes_densities); without it they are constant
      !       per panel
      !   stored  optional: the panel integrals of G from the stored arrays at Omega (default .true.; see the header)
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      COMPLEX(8), INTENT(IN) :: a(:, :), c(:, :, :), e(:, :)
      COMPLEX(8), INTENT(OUT) :: rhs(:, :, :)
      LOGICAL, INTENT(IN), OPTIONAL :: is_multi
      COMPLEX(8), INTENT(IN), OPTIONAL :: cl(:, :, :, :)
      LOGICAL, INTENT(IN), OPTIONAL :: stored
      COMPLEX(8), ALLOCATABLE :: r(:, :)
      INTEGER :: nb, ns, nrhs, i, ip, im
      LOGICAL :: multi, from_store

      nb = mesh%nbase
      ns = mesh%nsys
      nrhs = SIZE(a, 2)
      IF (SIZE(a, 1) /= mesh%n .OR. SIZE(c, 1) /= 3 .OR. SIZE(c, 2) /= mesh%n .OR. SIZE(c, 3) /= nrhs .OR. &
          SIZE(e, 1) /= geom%nwl .OR. SIZE(e, 2) /= nrhs .OR. SIZE(rhs, 1) /= nb .OR. SIZE(rhs, 2) /= nrhs .OR. &
          SIZE(rhs, 3) /= ns) THEN
         WRITE (*, '(a)') ' direct_rhs_b2_block: a, c, e or rhs does not have the size of the mesh.'
         ERROR STOP 1
      END IF
      multi = mesh%nbody > 1
      IF (PRESENT(is_multi)) multi = is_multi
      from_store = .TRUE.
      IF (PRESENT(stored)) from_store = stored

      ALLOCATE (r(nrhs, mesh%n))
      !$OMP PARALLEL DO NUM_THREADS(NTHREAD) SCHEDULE(DYNAMIC, 4)
      DO i = 1, mesh%n
         IF (from_store .AND. PRESENT(cl)) THEN
            CALL b2_stokes_at_point(mesh, geom, mesh%xc(:, i), a, c, e, multi, r(:, i), cl, i)
         ELSE IF (from_store) THEN
            CALL b2_stokes_at_point(mesh, geom, mesh%xc(:, i), a, c, e, multi, r(:, i), row=i)
         ELSE IF (PRESENT(cl)) THEN
            CALL b2_stokes_at_point(mesh, geom, mesh%xc(:, i), a, c, e, multi, r(:, i), cl)
         ELSE
            CALL b2_stokes_at_point(mesh, geom, mesh%xc(:, i), a, c, e, multi, r(:, i))
         END IF
      END DO
      !$OMP END PARALLEL DO

      ! r is -int G n.[(X.grad)u] dS at the full panels; the block form as ASSB_RBC
      DO ip = 1, ns
         DO i = 1, nb
            rhs(i, :, ip) = (0D0, 0D0)
            DO im = 1, ns
               rhs(i, :, ip) = rhs(i, :, ip) - RXY(ip, im)*r(:, (im - 1)*nb + i)
            END DO
         END DO
      END DO
   END SUBROUTINE direct_rhs_b2_block

   SUBROUTINE b2_stokes_at_point(mesh, geom, x, a, c, e, is_multi, out, cl, row)
      ! int_S G(x, .) n.[(X.grad)u] dS for each column of the densities, at the field point x (full coordinates):
      !   out = sum_p [ int_p G a_p + (int_p grad_xi G).c_p ] - sum_w (int_w G dl) e_w
      ! With cl (3, 4, n, ncol), on the near panels also sum_ab cl(a, b, p) int_p d_a(1/r + 1/r1) (xi - xi_p)_b dS and
      ! cl(:, 4, p).int_p (1/r + 1/r1)(xi - xi_p) dS.
      ! With row (a full panel index, whose centroid x must be), int_p G and int_p grad G come from the stored arrays
      ! (stored_panel_integrals), otherwise from green_panel_integrals.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      REAL(8), INTENT(IN) :: x(3)
      COMPLEX(8), INTENT(IN) :: a(:, :), c(:, :, :), e(:, :)
      LOGICAL, INTENT(IN) :: is_multi
      COMPLEX(8), INTENT(OUT) :: out(:)
      COMPLEX(8), INTENT(IN), OPTIONAL :: cl(:, :, :, :)
      INTEGER, INTENT(IN), OPTIONAL :: row
      COMPLEX(8) :: g_int, grad_int(3), g_line
      REAL(8) :: xv(3, 4), m0, g0(3), m1(3, 3), m1i(3, 3), mg(3), mgi(3)
      INTEGER :: p, iw, k, nv

      out = (0D0, 0D0)
      DO p = 1, mesh%n
         IF (PRESENT(row)) THEN
            CALL stored_panel_integrals(mesh, row, p, is_multi, g_int, grad_int)
         ELSE
            CALL green_panel_integrals(mesh, x, p, is_multi, g_int, grad_int)
         END IF
         DO k = 1, SIZE(out)
            out(k) = out(k) + g_int*a(p, k) + SUM(grad_int*c(:, p, k))
         END DO
         IF (.NOT. PRESENT(cl)) CYCLE
         IF (NORM2(x - mesh%xc(:, p)) > NEAR_PANEL_SIZES*mesh%psize(p)) CYCLE
         nv = mesh%ncn(p)
         xv(:, 1:nv) = mesh%xnode(:, mesh%ncon(1:nv, p))
         CALL rankine_moments(xv, nv, mesh%nrm(:, p), x, mesh%xc(:, p), m0, g0, m1, mg)
         CALL rankine_moments(xv, nv, mesh%nrm(:, p), [x(1), x(2), -x(3)], mesh%xc(:, p), m0, g0, m1i, mgi)
         m1 = m1 + m1i
         mg = mg + mgi
         DO k = 1, SIZE(out)
            out(k) = out(k) + SUM(m1*cl(:, 1:3, p, k)) + SUM(mg*cl(:, 4, p, k))
         END DO
      END DO
      DO iw = 1, geom%nwl
         CALL green_line_integral(x, geom%wl_x1(:, iw), geom%wl_x2(:, iw), mesh%psize(geom%wl_panel(iw)), g_line)
         out = out - g_line*e(iw, :)
      END DO
   END SUBROUTINE b2_stokes_at_point

   SUBROUTINE stored_panel_integrals(mesh, row, p, is_multi, g_int, grad_int)
      ! int_p G and int_p grad_xi G over full panel p at the centroid of full panel row, from the stored arrays at Omega
      ! (header): base row i, base panel j, stored image a xor b, gradient mirrored by M_b. The flag rule of CALGREEN.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      INTEGER, INTENT(IN) :: row, p
      LOGICAL, INTENT(IN) :: is_multi
      COMPLEX(8), INTENT(OUT) :: g_int, grad_int(3)
      COMPLEX(8) :: grn(4)
      REAL(8) :: rk(4), f(3), area, dist, psz
      INTEGER :: i, j, is

      i = mesh%base(row)
      j = mesh%base(p)
      is = MERGE(1, 2, mesh%image(row) == mesh%image(p))
      f = image_signs(mesh%image(p), ISX, ISY)
      IF (is_multi) THEN
         grn = CGRN_MULTI_COMB(i, j, is, :)
         rk = RKBN_MULTI_COMB(i, j, is, :)
         area = DS_MULTI_COMB(j)
         psz = PNSZ_MULTI_COMB(j)
         dist = NORM2(XYZ_GLOBAL_MULTI_COMB_P(i, :) - XYZ_GLOBAL_MULTI_COMB_P(j, :))
      ELSE
         grn = CGRN(i, j, is, :)
         rk = RKBN(i, j, is, :)
         area = DS(j)
         psz = PNSZ(j)
         dist = NORM2(XYZ_P(i, :) - XYZ_P(j, :))
      END IF
      IF (dist <= NEAR_PANEL_SIZES*psz) THEN
         g_int = rk(1) + grn(1)*area
         grad_int = f*(rk(2:4) + grn(2:4)*area)
      ELSE
         g_int = grn(1)*area
         grad_int = f*grn(2:4)*area
      END IF
   END SUBROUTINE stored_panel_integrals

   SUBROUTINE sharp_edge_panels(mesh, edge_deg, mask)
      ! mask(i) is .true. when panel i shares a node (after merging coincident nodes, so across the symmetry plane too)
      ! with a panel of the same body whose normal differs by more than edge_deg degrees: the panels at sharp edges,
      ! convex or concave, where the linear density is not used (header).
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: edge_deg
      LOGICAL, ALLOCATABLE, INTENT(OUT) :: mask(:)
      REAL(8), ALLOCATABLE :: bsize(:)
      INTEGER, ALLOCATABLE :: rep(:), nstart(:), nlist(:), fill(:)
      REAL(8) :: cosmin
      INTEGER :: i, m, k, r, r2, pi_, pj

      ALLOCATE (mask(mesh%n))
      mask = .FALSE.
      IF (mesh%n == 0) RETURN
      cosmin = COS(edge_deg*3.141592653589793D0/180D0)
      CALL body_sizes(mesh, bsize)
      CALL merge_nodes(mesh, bsize, rep)
      ALLOCATE (nstart(mesh%nnode + 1), fill(mesh%nnode))
      nstart = 0
      DO i = 1, mesh%n
         DO m = 1, mesh%ncn(i)
            k = rep(mesh%ncon(m, i))
            nstart(k + 1) = nstart(k + 1) + 1
         END DO
      END DO
      nstart(1) = 1
      DO k = 1, mesh%nnode
         nstart(k + 1) = nstart(k) + nstart(k + 1)
      END DO
      fill = nstart(1:mesh%nnode)
      ALLOCATE (nlist(nstart(mesh%nnode + 1) - 1))
      DO i = 1, mesh%n
         DO m = 1, mesh%ncn(i)
            k = rep(mesh%ncon(m, i))
            nlist(fill(k)) = i
            fill(k) = fill(k) + 1
         END DO
      END DO
      DO k = 1, mesh%nnode
         DO r = nstart(k), nstart(k + 1) - 1
            pi_ = nlist(r)
            DO r2 = r + 1, nstart(k + 1) - 1
               pj = nlist(r2)
               IF (mesh%body(pi_) /= mesh%body(pj)) CYCLE
               IF (DOT_PRODUCT(mesh%nrm(:, pi_), mesh%nrm(:, pj)) < cosmin) THEN
                  mask(pi_) = .TRUE.
                  mask(pj) = .TRUE.
               END IF
            END DO
         END DO
      END DO
   END SUBROUTINE sharp_edge_panels

   SUBROUTINE green_panel_integrals(mesh, x, p, is_multi, g_int, grad_int)
      ! int_p G(x, xi) dS_xi and int_p grad_xi G(x, xi) dS_xi over full panel p (image b of base panel j), x in full
      ! coordinates. Near panels (within NEAR_PANEL_SIZES panel sizes, as CALGREEN) take 1/r and 1/r1 as panel integrals,
      ! far ones the centroid value. On the own panel the 1/r part of the gradient is the principal value (in-plane part).
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: x(3)
      INTEGER, INTENT(IN) :: p
      LOGICAL, INTENT(IN) :: is_multi
      COMPLEX(8), INTENT(OUT) :: g_int, grad_int(3)
      REAL(8) :: f(3), xf(3), xq(3), d(3), r1, area, rk(4), rk_image(4)
      COMPLEX(8) :: grn(4)
      INTEGER :: j, flag

      j = mesh%base(p)
      f = image_signs(mesh%image(p), ISX, ISY)
      xf = f*x                 ! the field point in the frame where panel p is base panel j
      xq = f*mesh%xc(:, p)     ! the centroid of base panel j
      area = mesh%area(p)
      flag = MERGE(1, 0, NORM2(x - mesh%xc(:, p)) <= NEAR_PANEL_SIZES*mesh%psize(p))
      ! Centroid first: GRN(2:4) is the gradient in the source point
      IF (H < 0D0) THEN
         CALL INFGREEN3D(xq(1), xf(1), xq(2), xf(2), xq(3), xf(3), V, grn, flag)
      ELSE
         CALL FINGREEN3D(xq(1), xf(1), xq(2), xf(2), xq(3), xf(3), V, WVN, NK, H, grn, flag)
      END IF
      IF (flag /= 1) THEN
         g_int = grn(1)*area
         grad_int = f*grn(2:4)*area
         RETURN
      END IF

      ! As source_panel_integrals: the V = -1 branch of IMAGE_TERM_EXACT is left out (limit frequencies return before).
      ! rk(2:4) and rk_image(2:4) are the panel integrals of the gradient in xi of 1/r and 1/r1.
      CALL rankine_panel_integral(xf, j, mesh%ncn(p), is_multi, rk)
      CALL rankine_panel_integral([xf(1), xf(2), -xf(3)], j, mesh%ncn(p), is_multi, rk_image)
      d = [xq(1) - xf(1), xq(2) - xf(2), xq(3) + xf(3)]   ! grad_xi (1/r1) = -d/r1^3 at the centroid
      r1 = NORM2(d)
      g_int = rk(1) + rk_image(1) + (grn(1) - 1D0/r1)*area
      grad_int = f*(rk(2:4) + rk_image(2:4) + (grn(2:4) + d/r1**3)*area)
   END SUBROUTINE green_panel_integrals

   SUBROUTINE green_line_integral(x, x1, x2, psize, g_line)
      ! int G(x, xi) dl over the straight segment from x1 to x2 on z = 0 (the z of x1 and x2 is set to 0, where the
      ! waterline midpoints are). 1/r and 1/r1 in closed form, int dl/r = ln((R1 + R2 + L)/(R1 + R2 - L)), the rest at the
      ! midpoint. psize decides near and far as for the panels. x must not lie on the segment.
      REAL(8), INTENT(IN) :: x(3), x1(3), x2(3), psize
      COMPLEX(8), INTENT(OUT) :: g_line
      REAL(8) :: a1(3), a2(3), xm(3), xi(3), d(3), len, r1
      COMPLEX(8) :: grn(4)
      INTEGER :: flag

      a1 = [x1(1), x1(2), 0D0]
      a2 = [x2(1), x2(2), 0D0]
      xm = 0.5D0*(a1 + a2)
      len = NORM2(a2 - a1)
      flag = MERGE(1, 0, NORM2(x - xm) <= NEAR_PANEL_SIZES*psize)
      IF (H < 0D0) THEN
         CALL INFGREEN3D(xm(1), x(1), xm(2), x(2), xm(3), x(3), V, grn, flag)
      ELSE
         CALL FINGREEN3D(xm(1), x(1), xm(2), x(2), xm(3), x(3), V, WVN, NK, H, grn, flag)
      END IF
      IF (flag /= 1) THEN
         g_line = grn(1)*len
         RETURN
      END IF
      xi = [x(1), x(2), -x(3)]
      d = [xm(1) - x(1), xm(2) - x(2), xm(3) + x(3)]
      r1 = NORM2(d)
      g_line = segment_inverse_r(x) + segment_inverse_r(xi) + (grn(1) - 1D0/r1)*len

   CONTAINS

      REAL(8) FUNCTION segment_inverse_r(y)
         ! int dl/|y - xi| over the segment
         REAL(8), INTENT(IN) :: y(3)
         REAL(8) :: q1, q2
         q1 = NORM2(y - a1)
         q2 = NORM2(y - a2)
         segment_inverse_r = LOG((q1 + q2 + len)/(q1 + q2 - len))
      END FUNCTION segment_inverse_r

   END SUBROUTINE green_line_integral

   PURE SUBROUTINE rankine_moments(xv, nv, nrm, x, xc, m0, g0, m1, mg)
      ! Moments of 1/r = 1/|x - xi| over the flat panel with vertices xv(:, 1:nv) (3 or 4, projected on the plane through
      ! their average with normal nrm), in closed form:
      !   m0 = int 1/r dS
      !   g0 = int grad_xi(1/r) dS, in-plane part
      !   m1 = int grad_xi(1/r) (xi - xc)^T dS, both indices in-plane (xc projected on the plane)
      !   mg = int (1/r) (xi - xc) dS, in-plane part (optional)
      ! With X the projection of x, Z its height, rho = xi - X in the plane, and for each edge its outward in-plane normal nu,
      ! direction tau, length L, end distances R1, R2 (3D), foot parameter s0 and h = (P1 - X).nu:
      !   int_edge dl/r = I0 = ln((R1 + R2 + L)/(R1 + R2 - L))
      !   m0 = sum_e [ h I0 + |Z| (atan(|Z| t2/(h R2)) - atan(|Z| t1/(h R1)) - atan(t2/h) + atan(t1/h)) ],
      !        t1 = -s0, t2 = L - s0, from div[rho (r - |Z|)/rho^2] = 1/r in the plane
      !   g0 = sum_e nu I0
      !   int d_a(1/r) rho_b dS = sum_e nu_a [ h nu_b I0 + tau_b (R2 - R1) ] - delta_ab m0      (by parts in the plane)
      ! and m1 = that + g0 (X - xc)^T. With grad_xi r = (xi - x)/r, int (1/r) rho dS = sum_e nu oint_e r dl, with
      ! int r dt = [t r + d^2 asinh(t/d)]/2, d^2 = h^2 + Z^2, and mg = that + m0 (X - xc). With x on the panel (Z = 0) the integrand of m1 is O(1/rho), integrable, and g0 is the
      ! principal value: the circle about x adds nothing as its radius goes to 0, the convention of SGLINTBD_*2. x must not
      ! lie on an edge. Tested against adaptive quadrature (test_qtf_direct, rankine_moments_quadrature).
      REAL(8), INTENT(IN) :: xv(:, :), nrm(3), x(3), xc(3)
      INTEGER, INTENT(IN) :: nv
      REAL(8), INTENT(OUT) :: m0, g0(3), m1(3, 3)
      REAL(8), INTENT(OUT), OPTIONAL :: mg(3)
      REAL(8) :: e1(3), e2(3), o(3), q(2, 4), xp(2), qc(2), cen(2), z, az, p1(2), p2(2), tau(2), nu(2), el, s0, h, &
                 r1, r2, i0, t1, t2, g2(2), m2(2, 2), tiny, mr(2), dd
      INTEGER :: k

      o = SUM(xv(:, 1:nv), 2)/nv
      e1 = xv(:, 1) - o
      e1 = e1 - DOT_PRODUCT(e1, nrm)*nrm
      e1 = e1/NORM2(e1)
      e2 = [nrm(2)*e1(3) - nrm(3)*e1(2), nrm(3)*e1(1) - nrm(1)*e1(3), nrm(1)*e1(2) - nrm(2)*e1(1)]
      DO k = 1, nv
         q(:, k) = [DOT_PRODUCT(xv(:, k) - o, e1), DOT_PRODUCT(xv(:, k) - o, e2)]
      END DO
      cen = SUM(q(:, 1:nv), 2)/nv
      xp = [DOT_PRODUCT(x - o, e1), DOT_PRODUCT(x - o, e2)]
      qc = [DOT_PRODUCT(xc - o, e1), DOT_PRODUCT(xc - o, e2)]
      z = DOT_PRODUCT(x - o, nrm)
      az = ABS(z)

      m0 = 0D0
      g2 = 0D0
      m2 = 0D0
      mr = 0D0
      DO k = 1, nv
         p1 = q(:, k)
         p2 = q(:, MOD(k, nv) + 1)
         el = NORM2(p2 - p1)
         IF (el <= 0D0) CYCLE
         tiny = 1D-12*el
         tau = (p2 - p1)/el
         nu = [tau(2), -tau(1)]
         IF (DOT_PRODUCT(nu, 0.5D0*(p1 + p2) - cen) < 0D0) nu = -nu
         s0 = DOT_PRODUCT(xp - p1, tau)
         h = DOT_PRODUCT(p1 - xp, nu)
         r1 = SQRT(SUM((p1 - xp)**2) + z*z)
         r2 = SQRT(SUM((p2 - xp)**2) + z*z)
         i0 = 0D0
         IF (r1 + r2 - el > tiny) i0 = LOG((r1 + r2 + el)/(r1 + r2 - el))
         IF (ABS(h) > tiny) THEN
            t1 = -s0
            t2 = el - s0
            m0 = m0 + h*i0
            IF (az > 0D0) m0 = m0 + az*(ATAN(az*t2/(h*r2)) - ATAN(az*t1/(h*r1)) - ATAN(t2/h) + ATAN(t1/h))
         END IF
         g2 = g2 + nu*i0
         dd = SQRT(h*h + z*z)
         mr = mr + nu*(r_line(el - s0, r2) - r_line(-s0, r1))
         m2(:, 1) = m2(:, 1) + nu*(h*nu(1)*i0 + tau(1)*(r2 - r1))
         m2(:, 2) = m2(:, 2) + nu*(h*nu(2)*i0 + tau(2)*(r2 - r1))
      END DO
      m2(1, 1) = m2(1, 1) - m0
      m2(2, 2) = m2(2, 2) - m0
      m2(:, 1) = m2(:, 1) + g2*(xp(1) - qc(1))
      m2(:, 2) = m2(:, 2) + g2*(xp(2) - qc(2))

      g0 = g2(1)*e1 + g2(2)*e2
      DO k = 1, 3
         m1(k, :) = (m2(1, 1)*e1(k) + m2(2, 1)*e2(k))*e1 + (m2(1, 2)*e1(k) + m2(2, 2)*e2(k))*e2
      END DO
      IF (PRESENT(mg)) THEN
         mr = mr + m0*(xp - qc)
         mg = mr(1)*e1 + mr(2)*e2
      END IF

   CONTAINS

      PURE REAL(8) FUNCTION r_line(t, rt)
         ! The primitive of r = sqrt(t^2 + dd^2) at t, with rt = r there
         REAL(8), INTENT(IN) :: t, rt
         r_line = 0.5D0*t*rt
         IF (dd > 0D0) r_line = r_line + 0.5D0*dd*dd*ASINH(t/dd)
      END FUNCTION r_line

   END SUBROUTINE rankine_moments

   PURE FUNCTION cross_cr(a, b) RESULT(c)
      ! a x b for complex a and real b
      COMPLEX(8), INTENT(IN) :: a(3)
      REAL(8), INTENT(IN) :: b(3)
      COMPLEX(8) :: c(3)

      c(1) = a(2)*b(3) - a(3)*b(2)
      c(2) = a(3)*b(1) - a(1)*b(3)
      c(3) = a(1)*b(2) - a(2)*b(1)
   END FUNCTION cross_cr

END MODULE QTFDirect
