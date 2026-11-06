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
   ! The total first-order potential and its gradient at points in the fluid, for the control-surface drift (design spec
   ! docs/superpowers/specs/2026-11-06-control-surface-drift-design.md, section 4.1). The incident part is exact. The
   ! perturbation part phi_P comes from the boundary integral formula of WaterlineField, and its gradient from central
   ! differences of phi_P with step delta:
   !   - below z = -delta: the 6 points x +- delta e_c, second order in delta;
   !   - from z = -delta up: the 4 horizontal points, with dphi_P/dz = nu phi_P, the linear free-surface condition, and the
   !     exact incident dphi_I/dz. No stencil point is above z = 0.
   ! Differences are used because the analytic gradient of the formula needs the gradient of dG/dn, a second derivative of G,
   ! which FINGREEN3D and INFGREEN3D do not give. All stencil points of a point take the near-or-far decision of its centre
   ! (x_decide of field_influence_row), so that no jump between the two panel integrals enters a difference.
   ! The rows of one stencil point serve every state before the next point, and no row is stored: the rows of all points
   ! would take gigabytes (section 4.3).

   USE SecondOrderTypes, ONLY: panel_mesh, first_order_state, cs_points, cs_fields, CS_WALL
   USE FirstOrderState, ONLY: incident_phi, incident_grad, wave_gravity
   USE WaterlineField, ONLY: field_influence_row, field_potential
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: fields_green

CONTAINS

   SUBROUTINE fields_green(mesh, pts, states, step, is_multi, fields)
      ! phi and grad phi of every state at every point of pts. step is delta as a fraction of the point size pts%psize.
      ! phi is computed at the centre for points from z = -delta up, which need it for dphi/dz, and for CS_LINE and CS_FREE
      ! points. At a CS_WALL point further down it is not needed, and it is left zero.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(cs_points), INTENT(IN) :: pts
      TYPE(first_order_state), INTENT(IN) :: states(:)
      REAL(8), INTENT(IN) :: step
      LOGICAL, INTENT(IN) :: is_multi
      TYPE(cs_fields), INTENT(OUT) :: fields(SIZE(states))
      COMPLEX(8), ALLOCATABLE :: row_g(:, :), row_dg(:, :)
      COMPLEX(8) :: phi_c(SIZE(states)), phi_plus(SIZE(states)), dphi(3, SIZE(states)), grad_i(3)
      REAL(8) :: x0(3), xs(3), delta, g_w
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
            CALL field_influence_row(mesh, x0, is_multi, row_g, row_dg, x0)
            DO ist = 1, nst
               phi_c(ist) = field_potential(mesh, states(ist), row_g, row_dg)
            END DO
         END IF

         DO c = 1, ncomp
            xs = x0
            xs(c) = x0(c) + delta
            CALL field_influence_row(mesh, xs, is_multi, row_g, row_dg, x0)
            DO ist = 1, nst
               phi_plus(ist) = field_potential(mesh, states(ist), row_g, row_dg)
            END DO
            xs(c) = x0(c) - delta
            CALL field_influence_row(mesh, xs, is_multi, row_g, row_dg, x0)
            DO ist = 1, nst
               dphi(c, ist) = (phi_plus(ist) - field_potential(mesh, states(ist), row_g, row_dg))/(2D0*delta)
            END DO
         END DO

         DO ist = 1, nst
            ASSOCIATE (s => states(ist))
               ! The incident wave the solver used, as in build_state
               g_w = wave_gravity(s%omega, s%k, s%h)
               CALL incident_grad(x0, s%omega, s%k, s%h, s%amp, s%beta, g_w, grad_i)
               IF (with_centre) fields(ist)%phi(i) = incident_phi(x0, s%omega, s%k, s%h, s%amp, s%beta, g_w) + phi_c(ist)
               fields(ist)%grad(1:ncomp, i) = grad_i(1:ncomp) + dphi(1:ncomp, ist)
               ! Each part exactly: s%nu is V = omega^2/g of Const_mod, which phi_P satisfies through the Green function, while the
               ! incident wave has omega^2/g_w = k tanh kh
               IF (at_surface) fields(ist)%grad(3, i) = grad_i(3) + s%nu*phi_c(ist)
            END ASSOCIATE
         END DO
      END DO
   END SUBROUTINE fields_green

END MODULE FieldVelocity
