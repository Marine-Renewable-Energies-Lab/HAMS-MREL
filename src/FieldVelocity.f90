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
MODULE FieldVelocity
   ! The total first-order potential and its gradient at points in the fluid, for the control-surface drift (Chen
   ! 2007, eq. 39 to 41). The incident part is exact. The
   ! perturbation part phi_P comes from the boundary integral formula of WaterlineField, and its gradient from central
   ! differences of phi_P with step delta:
   !   - below z = -delta: the 6 points x +- delta e_c, second order in delta [derived: the error of a central difference is
   !     O(delta^2)];
   !   - from z = -delta up: the 4 horizontal points, with dphi_P/dz = nu phi_P, the linear free-surface condition, and the
   !     exact incident dphi_I/dz. No stencil point is above z = 0. Below z = 0 but within -delta, nu phi_P for dphi_P/dz is an
   !     O(delta) approximation, exact only at z = 0 (the linear condition holds on z = 0); the generated Gamma_c and F points sit
   !     at z = 0.
   ! Differences are used because the analytic gradient of the formula needs the gradient of dG/dn, a second derivative of G,
   ! which FINGREEN3D and INFGREEN3D do not give. All stencil points of a point take the near-or-far decision of its centre
   ! (x_decide of field_influence_row), so that no jump between the two panel integrals enters a difference. They also take
   ! the centre's choice of the FINGREEN3D series, which do not agree exactly where FINGREEN3D_R switches between them.
   ! And they take field_mode: FINGREEN3D then sums its series directly, since Limes in EIGENE and PIDCOCK returns
   ! isolated outliers, which a difference takes over 2 delta (FinGreen3D.f90, FINGREEN3D_R). The exception is PIDCOCK
   ! for nu*h > 50, which keeps Limes (FIELD_SERIES in FinGreen3D.f90).
   ! The rows of one stencil point serve every state before the next point, and no row is stored: the rows of all points
   ! would take (stencil points x N_C) x N_hull x nsys x 2 complex numbers, about 5 GB at 7168 hull panels and 4000 C
   ! panels [derived].
   !
   ! The second route, fields_source, takes phi_P and its gradient from the source strengths sigma of the auxiliary source
   ! solve of SurfaceVelocity: phi_P = (1/4 pi) sum sigma int G dS, and the gradient with the analytic gradient of G, so no
   ! stencil is needed.

   USE SecondOrderTypes, ONLY: panel_mesh, first_order_state, cs_points, cs_fields, CS_WALL
   USE FirstOrderState, ONLY: incident_phi, incident_grad, wave_gravity, image_signs
   USE HAMS_mod, ONLY: NTHREAD
   USE PanelMesh_mod, ONLY: ISX, ISY
   USE SurfaceVelocity, ONLY: source_panel_integrals
   USE WaterlineField, ONLY: field_influence_row, field_potential
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: fields_green, fields_source

   REAL(8), PARAMETER :: PI = 3.141592653589793D0
   ! As build_source_operator: a source panel closer than this many panel sizes gets panel integrals of 1/r and of its image
   REAL(8), PARAMETER :: NEAR_PANEL_SIZES = 50D0

CONTAINS

   SUBROUTINE fields_green(mesh, pts, states, step, is_multi, fields)
      ! phi and grad phi of every state at every point of pts. step is delta as a fraction of the point size pts%psize.
      ! phi is computed at the centre for points from z = -delta up, which need it for dphi/dz, and for CS_LINE and CS_FREE
      ! points. At a CS_WALL point further down it is not needed, and it is left zero.
      ! All states must be at the current global frequency: field_influence_row reads the module globals V, WVN and H.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(cs_points), INTENT(IN) :: pts
      TYPE(first_order_state), INTENT(IN) :: states(:)
      REAL(8), INTENT(IN) :: step
      LOGICAL, INTENT(IN) :: is_multi
      TYPE(cs_fields), INTENT(OUT) :: fields(SIZE(states))
      COMPLEX(8), ALLOCATABLE :: row_g(:, :), row_dg(:, :)
      COMPLEX(8) :: phi_c(SIZE(states)), phi_plus(SIZE(states)), dphi(3, SIZE(states)), grad_i(3), phi_i
      REAL(8) :: x0(3), xs(3), delta
      LOGICAL :: at_surface, with_centre
      INTEGER :: i, c, ncomp, ist, nst

      nst = SIZE(states)
      ALLOCATE (row_g(mesh%nbase, mesh%nsys), row_dg(mesh%nbase, mesh%nsys))
      DO ist = 1, nst
         ALLOCATE (fields(ist)%phi(pts%n), fields(ist)%grad(3, pts%n))
         fields(ist)%phi = (0D0, 0D0)
         fields(ist)%grad = (0D0, 0D0)
      END DO

      DO i = 1, pts%n
         x0 = pts%x(:, i)
         delta = step*pts%psize(i)
         at_surface = x0(3) >= -delta
         with_centre = at_surface .OR. pts%kind(i) /= CS_WALL
         ncomp = MERGE(2, 3, at_surface)

         phi_c = (0D0, 0D0)
         IF (with_centre) THEN
            CALL field_influence_row(mesh, x0, is_multi, row_g, row_dg, x0, field_mode=.TRUE.)
            DO ist = 1, nst
               phi_c(ist) = field_potential(mesh, states(ist), row_g, row_dg)
            END DO
         END IF

         DO c = 1, ncomp
            xs = x0
            xs(c) = x0(c) + delta
            CALL field_influence_row(mesh, xs, is_multi, row_g, row_dg, x0, field_mode=.TRUE.)
            DO ist = 1, nst
               phi_plus(ist) = field_potential(mesh, states(ist), row_g, row_dg)
            END DO
            xs(c) = x0(c) - delta
            CALL field_influence_row(mesh, xs, is_multi, row_g, row_dg, x0, field_mode=.TRUE.)
            DO ist = 1, nst
               dphi(c, ist) = (phi_plus(ist) - field_potential(mesh, states(ist), row_g, row_dg))/(2D0*delta)
            END DO
         END DO

         DO ist = 1, nst
            ASSOCIATE (s => states(ist))
               CALL incident_at(s, x0, phi_i, grad_i)
               IF (with_centre) fields(ist)%phi(i) = phi_i + phi_c(ist)
               fields(ist)%grad(1:ncomp, i) = grad_i(1:ncomp) + dphi(1:ncomp, ist)
               ! Each part exactly: s%nu is V = omega^2/g of Const_mod, which phi_P satisfies through the Green function, while the
               ! incident wave has omega^2/g_w = k tanh kh
               IF (at_surface) fields(ist)%grad(3, i) = grad_i(3) + s%nu*phi_c(ist)
            END ASSOCIATE
         END DO
      END DO
   END SUBROUTINE fields_green

   SUBROUTINE fields_source(mesh, pts, states, sigmas, is_multi, fields)
      ! phi and grad phi of every state at every point of pts, from the source strengths of velocity_source. With the image of
      ! base panel j under symmetry element is at full index (is-1)*nbase + j, and G(x, M q) = G(M x, q):
      !   phi_P(x)      = (1/4 pi) sum_is sum_j sigma g_int(M x, j)
      !   grad phi_P(x) = (1/4 pi) sum_is sum_j sigma M grad_int(M x, j)
      ! with g_int and grad_int from source_panel_integrals. The field point is not on a panel, so there is no jump term.
      ! sigmas is (n, nstate). A source panel within 50 panel sizes of the image of the point takes panel integrals of 1/r and of
      ! its image (flag = 1), as in build_source_operator. phi is computed at every point, CS_WALL points included, since here
      ! it costs nothing extra. All states must be at the current global frequency: source_panel_integrals reads the module
      ! globals V, WVN and H.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(cs_points), INTENT(IN) :: pts
      TYPE(first_order_state), INTENT(IN) :: states(:)
      COMPLEX(8), INTENT(IN) :: sigmas(:, :)
      LOGICAL, INTENT(IN) :: is_multi
      TYPE(cs_fields), INTENT(OUT) :: fields(SIZE(states))
      COMPLEX(8) :: g_int, grad_int(3), phi_p(SIZE(states)), grad_p(3, SIZE(states)), grad_i(3), phi_i, sig
      REAL(8) :: x0(3), xf(3), f(3)
      INTEGER :: i, is, j, ist, nst, nb, ns, flag

      nst = SIZE(states)
      nb = mesh%nbase
      ns = mesh%nsys
      DO ist = 1, nst
         ALLOCATE (fields(ist)%phi(pts%n), fields(ist)%grad(3, pts%n))
      END DO

      !$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(x0, xf, f, is, j, ist, flag, g_int, grad_int, phi_p, grad_p, grad_i, phi_i, sig) &
      !$OMP SCHEDULE(DYNAMIC, 4)
      DO i = 1, pts%n
         x0 = pts%x(:, i)
         phi_p = (0D0, 0D0)
         grad_p = (0D0, 0D0)
         DO is = 1, ns
            f = image_signs(is, ISX, ISY)
            xf = f*x0
            DO j = 1, nb
               flag = MERGE(1, 0, NORM2(xf - mesh%xc(:, j)) <= NEAR_PANEL_SIZES*mesh%psize(j))
               CALL source_panel_integrals(mesh, xf, j, flag, is_multi, g_int, grad_int)
               DO ist = 1, nst
                  sig = sigmas((is - 1)*nb + j, ist)
                  phi_p(ist) = phi_p(ist) + sig*g_int
                  grad_p(:, ist) = grad_p(:, ist) + sig*f*grad_int
               END DO
            END DO
         END DO
         DO ist = 1, nst
            CALL incident_at(states(ist), x0, phi_i, grad_i)
            fields(ist)%phi(i) = phi_i + phi_p(ist)/(4D0*PI)
            fields(ist)%grad(:, i) = grad_i + grad_p(:, ist)/(4D0*PI)
         END DO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE fields_source

   SUBROUTINE incident_at(s, x, phi_i, grad_i)
      ! The incident wave the solver used, as in build_state, with g_w = wave_gravity(omega, k, h)
      TYPE(first_order_state), INTENT(IN) :: s
      REAL(8), INTENT(IN) :: x(3)
      COMPLEX(8), INTENT(OUT) :: phi_i, grad_i(3)
      REAL(8) :: g_w

      g_w = wave_gravity(s%omega, s%k, s%h)
      phi_i = incident_phi(x, s%omega, s%k, s%h, s%amp, s%beta, g_w)
      CALL incident_grad(x, s%omega, s%k, s%h, s%amp, s%beta, g_w, grad_i)
   END SUBROUTINE incident_at

END MODULE FieldVelocity
