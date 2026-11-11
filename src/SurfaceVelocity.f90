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
   ! Fluid velocity on the body panels. HAMS solves for the potential, not its gradient, so there are two ways to it.
   !
   ! By differencing: the tangential part from the potential over the panel neighbours, and the normal part from the known
   ! normal derivative.
   !
   ! By source solve: phi_p = (1/4 pi) int sigma G dS, with normals out of the fluid and sigma from
   !
   !    (2 pi I + D) sigma = 4 pi dphi_p/dn,   D_ij = n_i . int_j grad_p G dS
   !    grad phi_p(i) = (1/4 pi) [ 2 pi sigma_i n_i + sum_j sigma_j int_j grad_p G dS ]
   !
   ! with grad_p the gradient with respect to the field point. The stored arrays CGRN and RKBN are not used. For the entry (i, j)
   ! they hold derivatives with respect to the source point, and on near panels RKBN adds the z parts of 1/r and of its image,
   ! whose z gradients with respect to the field point have opposite signs. In finite depth the wave part is not a function of
   ! z_p + z_q alone, so its z gradient with respect to the field point is a different number. The transposed entry
   ! CGRN(j, i, is, 2:4) does hold the wave part of the gradient, from FINGREEN3D with xq = x_i and xp = M_is x_j. So a stored
   ! route exists for the most part. The Rankine z part would still need one fresh panel integral for each near pair, and the
   ! stored values are single precision. What it would save has not been measured. This route calls the Green function again
   ! in double precision, with the field point in the first place, in source_panel_integrals.

   USE Const_mod, ONLY: RXY
   USE HAMS_mod, ONLY: NTHREAD
   USE WaveDyn_mod, ONLY: H, NK, V, WVN
   USE PanelMesh_mod, ONLY: ISX, ISY
   USE INFG3D_Open, ONLY: INFGREEN3D
   USE FinGreen3D_Open, ONLY: FINGREEN3D
   USE SecondOrderTypes, ONLY: panel_mesh, surface_geom, first_order_state
   USE FirstOrderState, ONLY: image_signs
   USE WaterlineField, ONLY: rankine_panel_integral
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: velocity_differencing, tangential_fit
   PUBLIC :: source_operator, source_panel_integrals, build_source_operator, velocity_source

   REAL(8), PARAMETER :: PI = 3.141592653589793D0
   ! Normalised determinant of the 2x2 normal equations below which the neighbours are treated as collinear
   REAL(8), PARAMETER :: DEGENERATE_DET = 1D-6
   ! As CALGREEN: a source panel closer than this many panel sizes gets panel integrals of 1/r and of its image
   REAL(8), PARAMETER :: NEAR_PANEL_SIZES = 50D0

   TYPE :: source_operator
      ! The source solve of one frequency, per symmetry block ip: 1 alone without symmetry; 1 even and 2 odd in the image.
      COMPLEX(8), ALLOCATABLE :: lu(:, :, :)            ! (nbase, nbase, nsys), factors of 2 pi I + D
      INTEGER, ALLOCATABLE :: ipiv(:, :)                ! (nbase, nsys)
      COMPLEX(8), ALLOCATABLE :: gradient(:, :, :, :)   ! (3, j, i, ip), int_j grad_p G dS at panel i
      COMPLEX(8), ALLOCATABLE :: potential(:, :, :)     ! (i, j, ip), int_j G dS at panel i. Only for diagnostics
   END TYPE source_operator

CONTAINS

   SUBROUTINE velocity_differencing(mesh, geom, s, n_fallback)
      ! Fills s%grad_p, the gradient of the body-disturbance potential phi_p at the panel centroids:
      !   grad_p = n dphi_p_dn + t1 d(phi_p)/dt1 + t2 d(phi_p)/dt2
      ! with the tangential derivatives from tangential_fit. The total velocity is s%grad_i + s%grad_p. A panel with fewer
      ! than two non-collinear neighbours keeps only its normal part, and is counted in n_fallback. s%grad_p must be
      ! allocated, as build_state does.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(first_order_state), INTENT(INOUT) :: s
      INTEGER, INTENT(OUT) :: n_fallback

      INTEGER :: i
      REAL(8) :: t1(3), t2(3)
      COMPLEX(8) :: g1, g2
      LOGICAL :: ok

      n_fallback = 0
      DO i = 1, mesh%n
         s%grad_p(:, i) = mesh%nrm(:, i)*s%dphi_p_dn(i)
         CALL tangential_fit(mesh, geom, i, s%phi_p, t1, t2, g1, g2, ok)
         IF (.NOT. ok) THEN
            n_fallback = n_fallback + 1
            CYCLE
         END IF
         s%grad_p(:, i) = s%grad_p(:, i) + t1*g1 + t2*g2
      END DO
   END SUBROUTINE velocity_differencing

   SUBROUTINE tangential_fit(mesh, geom, i, f, t1, t2, g1, g2, ok)
      ! Tangential derivatives g1 = df/dt1, g2 = df/dt2 of a panel scalar f at panel i, with orthonormal tangents t1, t2:
      ! the weighted least-squares fit (weight 1/|d|^2) of f_j - f_i over the neighbours j, with d the centroid offset
      ! projected on the tangents. ok is .false. when the panel has fewer than two non-collinear neighbours; g1 and g2 are
      ! then zero.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      INTEGER, INTENT(IN) :: i
      COMPLEX(8), INTENT(IN) :: f(:)
      REAL(8), INTENT(OUT) :: t1(3), t2(3)
      COMPLEX(8), INTENT(OUT) :: g1, g2
      LOGICAL, INTENT(OUT) :: ok

      INTEGER :: m, j
      REAL(8) :: nv(3), a(3), d(3), d1, d2, r2, w
      REAL(8) :: a11, a12, a22, det
      COMPLEX(8) :: b1, b2, df

      nv = mesh%nrm(:, i)

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
         df = f(j) - f(i)
         a11 = a11 + w*d1*d1
         a12 = a12 + w*d1*d2
         a22 = a22 + w*d2*d2
         b1 = b1 + w*d1*df
         b2 = b2 + w*d2*df
      END DO

      g1 = (0D0, 0D0)
      g2 = (0D0, 0D0)
      det = a11*a22 - a12*a12
      ok = .NOT. (a11*a22 <= 0D0 .OR. det < DEGENERATE_DET*a11*a22)
      IF (.NOT. ok) RETURN
      g1 = (a22*b1 - a12*b2)/det
      g2 = (a11*b2 - a12*b1)/det
   END SUBROUTINE tangential_fit

   SUBROUTINE source_panel_integrals(mesh, xf, jel, flag, is_multi, g_int, grad_int)
      ! G and its gradient with respect to the field point, integrated over base panel jel, at xf, which is already reflected for
      ! the symmetry image. Written as the terms 1/r, its image 1/r1 and the rest, with q on the panel:
      !   1/r    a function of xf - q, so the gradient is minus the panel integral of the gradient with respect to q
      !   1/r1   of the horizontal part of xf - q and of z_f + z_q, so the z part keeps its sign
      !   rest   from the Green function, one point on the panel. With the field point first, GRN(2:4) is its gradient there
      ! flag = 1 takes 1/r and 1/r1 as panel integrals, like IMAGE_TERM_EXACT, and flag = 0 all at one point. The Green function
      ! holds 1/r1 at flag = 1, so it is taken out at the centroid.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: xf(3)
      INTEGER, INTENT(IN) :: jel, flag
      LOGICAL, INTENT(IN) :: is_multi
      COMPLEX(8), INTENT(OUT) :: g_int, grad_int(3)
      REAL(8) :: xq(3), d(3), r1, area, rk(4), rk_image(4)
      COMPLEX(8) :: grn(4)

      xq = mesh%xc(:, jel)
      area = mesh%area(jel)
      IF (H < 0D0) THEN
         CALL INFGREEN3D(xf(1), xq(1), xf(2), xq(2), xf(3), xq(3), V, grn, flag)
      ELSE
         CALL FINGREEN3D(xf(1), xq(1), xf(2), xq(2), xf(3), xq(3), V, WVN, NK, H, grn, flag)
      END IF
      IF (flag /= 1) THEN
         g_int = grn(1)*area
         grad_int = grn(2:4)*area
         RETURN
      END IF

      ! The V = -1 branch of IMAGE_TERM_EXACT is left out, since limit frequencies return before this code.
      CALL rankine_panel_integral(xf, jel, mesh%ncn(jel), is_multi, rk)
      CALL rankine_panel_integral([xf(1), xf(2), -xf(3)], jel, mesh%ncn(jel), is_multi, rk_image)
      d = [xf(1) - xq(1), xf(2) - xq(2), xf(3) + xq(3)]
      r1 = NORM2(d)
      g_int = rk(1) + rk_image(1) + (grn(1) - 1D0/r1)*area
      grad_int(1:2) = -rk(2:3) - rk_image(2:3) + (grn(2:3) + d(1:2)/r1**3)*area
      grad_int(3) = -rk(4) + rk_image(4) + (grn(4) + d(3)/r1**3)*area
   END SUBROUTINE source_panel_integrals

   SUBROUTINE build_source_operator(mesh, is_multi, with_potential, op)
      ! Assembles and factorizes the source solve for the current frequency. The panels stay in place and the field point is
      ! mirrored, as in ASSB_LEFT: G(x, M q) = G(M x, q), whose gradient in x is M grad_p G(M x, q). is_multi selects the panel
      ! integrals that read the combined multi-body arrays. with_potential also keeps int_j G dS.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      LOGICAL, INTENT(IN) :: is_multi, with_potential
      TYPE(source_operator), INTENT(INOUT) :: op
      INTEGER :: nb, ns, i, j, is, ip, flag, info
      REAL(8) :: f(3)
      COMPLEX(8) :: g_int, grad_int(3)

      nb = mesh%nbase
      ns = mesh%nsys
      IF (.NOT. ALLOCATED(op%lu)) THEN
         ALLOCATE (op%lu(nb, nb, ns), op%ipiv(nb, ns), op%gradient(3, nb, nb, ns))
         IF (with_potential) ALLOCATE (op%potential(nb, nb, ns))
      END IF
      op%gradient = (0D0, 0D0)
      IF (with_potential) op%potential = (0D0, 0D0)

      !$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(j, is, ip, flag, f, g_int, grad_int) SCHEDULE(DYNAMIC, 16)
      DO i = 1, nb
         DO j = 1, nb
            flag = MERGE(1, 0, NORM2(mesh%xc(:, i) - mesh%xc(:, j)) <= NEAR_PANEL_SIZES*mesh%psize(j))
            DO is = 1, ns
               f = image_signs(is, ISX, ISY)
               CALL source_panel_integrals(mesh, f*mesh%xc(:, i), j, flag, is_multi, g_int, grad_int)
               DO ip = 1, ns
                  op%gradient(:, j, i, ip) = op%gradient(:, j, i, ip) + RXY(is, ip)*f*grad_int
                  IF (with_potential) op%potential(i, j, ip) = op%potential(i, j, ip) + RXY(is, ip)*g_int
               END DO
            END DO
            DO ip = 1, ns
               op%lu(i, j, ip) = SUM(mesh%nrm(:, i)*op%gradient(:, j, i, ip))
            END DO
         END DO
         op%lu(i, i, :) = op%lu(i, i, :) + 2D0*PI
      END DO
      !$OMP END PARALLEL DO

      DO ip = 1, ns
         CALL ZGETRF(nb, nb, op%lu(:, :, ip), nb, op%ipiv(:, ip), info)
         IF (info /= 0) THEN
            WRITE (*, *) 'Mean drift: the source solve could not factorize its matrix, ZGETRF info = ', info
            ERROR STOP 1
         END IF
      END DO
   END SUBROUTINE build_source_operator

   SUBROUTINE velocity_source(mesh, op, s, phi_source, sigma_full)
      ! Fills s%grad_p, the same quantity as velocity_differencing, from s%dphi_p_dn. The even and odd parts of a state are solved
      ! apart. A panel and its image take their sum and difference, and the gradient at the image is mirrored.
      ! phi_source, if given, is (1/4 pi) int sigma G dS at the panels, for which op must hold the potential.
      ! sigma_full, if given, is sigma on every panel, rebuilt from the blocks as grad_p is, without the image sign (sigma is a
      ! scalar): sigma at image im of base panel i is sum_ip RXY(im, ip) sigma(i, ip)/ns, at full index (im-1) nb + i.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(source_operator), INTENT(IN) :: op
      TYPE(first_order_state), INTENT(INOUT) :: s
      COMPLEX(8), INTENT(OUT), OPTIONAL :: phi_source(:), sigma_full(:)
      COMPLEX(8) :: sigma(mesh%nbase, mesh%nsys), block_grad(3, mesh%nbase, mesh%nsys)
      COMPLEX(8) :: block_phi(mesh%nbase, mesh%nsys)
      INTEGER :: nb, ns, i, j, im, ip, info

      nb = mesh%nbase
      ns = mesh%nsys
      DO ip = 1, ns
         DO i = 1, nb
            sigma(i, ip) = (0D0, 0D0)
            DO im = 1, ns
               sigma(i, ip) = sigma(i, ip) + RXY(im, ip)*4D0*PI*s%dphi_p_dn((im - 1)*nb + i)
            END DO
         END DO
         CALL ZGETRS('No transpose', nb, 1, op%lu(:, :, ip), nb, op%ipiv(:, ip), sigma(:, ip), nb, info)
         IF (info /= 0) THEN
            WRITE (*, *) 'Mean drift: the source solve failed, ZGETRS info = ', info
            ERROR STOP 1
         END IF
      END DO

      !$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(ip, j) COLLAPSE(2)
      DO i = 1, nb
         DO ip = 1, ns
            block_grad(:, i, ip) = 2D0*PI*sigma(i, ip)*mesh%nrm(:, i)
            DO j = 1, nb
               block_grad(:, i, ip) = block_grad(:, i, ip) + op%gradient(:, j, i, ip)*sigma(j, ip)
            END DO
            IF (PRESENT(phi_source)) block_phi(i, ip) = SUM(op%potential(i, :, ip)*sigma(:, ip))
         END DO
      END DO
      !$OMP END PARALLEL DO

      DO im = 1, ns
         DO i = 1, nb
            s%grad_p(:, (im - 1)*nb + i) = (0D0, 0D0)
            DO ip = 1, ns
               s%grad_p(:, (im - 1)*nb + i) = s%grad_p(:, (im - 1)*nb + i) + RXY(im, ip)*block_grad(:, i, ip)
            END DO
            s%grad_p(:, (im - 1)*nb + i) = image_signs(im, ISX, ISY)*s%grad_p(:, (im - 1)*nb + i)/(4D0*PI*ns)
            IF (PRESENT(phi_source)) phi_source((im - 1)*nb + i) = SUM(RXY(im, 1:ns)*block_phi(i, 1:ns))/(4D0*PI*ns)
            IF (PRESENT(sigma_full)) sigma_full((im - 1)*nb + i) = SUM(RXY(im, 1:ns)*sigma(i, 1:ns))/ns
         END DO
      END DO
   END SUBROUTINE velocity_source

END MODULE SurfaceVelocity
