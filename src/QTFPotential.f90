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
   ! u at the waterline is extrapolated along the panel from its centroid value. With a differenced u the waterline row is
   ! fitted one-sided, which limits the line term to O(h) (measured in stokes_identity_differenced_velocity).

   USE SecondOrderTypes, ONLY: panel_mesh, surface_geom
   USE NearFieldDrift, ONLY: pair_cross, rigid_displacement, waterline_point
   USE SurfaceVelocity, ONLY: tangential_fit
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: hessian_term_converted, psi_waterline

CONTAINS

   SUBROUTINE hessian_term_converted(mesh, geom, psi, grad_psi, psi_wl, x1, u1, x2, u2, kernel, out)
      ! Kernel pair (1, 2) of int_S psi_k n.[(X.grad)u] dS for each mode column k, in the converted form above.
      !   psi        (n, nmode) radiation potentials at the panel centroids, nmode = 6 nbody
      !   grad_psi   (3, n, nmode) their gradients; only the tangential part enters
      !   psi_wl     (nwl, nmode) psi at the waterline segments (psi_waterline)
      !   x1, x2     (6, nbody) rigid motions, translations then rotations about xr of the body
      !   u1, u2     (3, n) total first-order velocities at the panel centroids
      !   kernel     KERNEL_DIFF or KERNEL_SUM of NearFieldDrift
      !   out        (nmode)
      ! The sum runs over the panels and waterline segments of all bodies: psi_k of one body acts on every surface, with
      ! the motion of the body that owns the panel.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      COMPLEX(8), INTENT(IN) :: psi(:, :), grad_psi(:, :, :), psi_wl(:, :)
      COMPLEX(8), INTENT(IN) :: x1(:, :), u1(:, :), x2(:, :), u2(:, :)
      INTEGER, INTENT(IN) :: kernel
      COMPLEX(8), INTENT(OUT) :: out(:)

      INTEGER :: i, k, b, iw, ip
      REAL(8) :: nv(3), x(3), t(3)
      COMPLEX(8) :: d1(3), d2(3), rot(3), xu(3), uw1(3), uw2(3)
      COMPLEX(8), ALLOCATABLE :: uc1(:, :), uc2(:, :)

      out = (0D0, 0D0)
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
      END DO

      IF (geom%nwl == 0) RETURN
      uc1 = TRANSPOSE(u1)   ! (n, 3): each component a contiguous panel scalar for tangential_fit
      uc2 = TRANSPOSE(u2)
      DO iw = 1, geom%nwl
         ip = geom%wl_panel(iw)
         b = geom%wl_body(iw)
         x = waterline_point(geom, iw)
         t = (geom%wl_x2(:, iw) - geom%wl_x1(:, iw))/geom%wl_len(iw)
         d1 = rigid_displacement(x1(:, b), x - mesh%xr(:, b))
         d2 = rigid_displacement(x2(:, b), x - mesh%xr(:, b))
         uw1 = waterline_velocity(uc1, ip, x)
         uw2 = waterline_velocity(uc2, ip, x)
         xu = pair_cross(kernel, d1, uw1, d2, uw2)
         DO k = 1, SIZE(out)
            out(k) = out(k) - psi_wl(iw, k)*SUM(xu*t)*geom%wl_len(iw)
         END DO
      END DO

   CONTAINS

      FUNCTION waterline_velocity(uc, ip, xw) RESULT(uw)
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
      END FUNCTION waterline_velocity

   END SUBROUTINE hessian_term_converted

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
