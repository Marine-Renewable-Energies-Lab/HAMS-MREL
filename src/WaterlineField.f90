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
MODULE WaterlineField
   ! The waterline potential by field evaluation. The perturbation potential comes from the boundary integral formula at a point
   ! in the fluid close to the waterline, and the incident potential is added exactly at z = 0:
   !
   !    phi_P(x) = (1 / 4 pi) sum_images sum_j [ G dphi_P/dn - (dG/dn) phi_P ]
   !
   ! with G and dG/dn integrated over panel j, as in CalPotential and CalPotentialMulti, the free-surface image term included. The
   ! point lies in the fluid, so the free term is 4 pi, if the wall at the waterline widens downward by no more than 45 degrees
   ! from the vertical (see field_point_inside_body) and no other body is within offset*wl_len of the point. Only hull panels are sources: with IRSP = 1 the lid acts through the
   ! hull potentials. The integrals depend on the frequency and not on the state, so one set of rows serves all the headings.

   USE HAMS_mod, ONLY: NTHREAD
   USE WaveDyn_mod, ONLY: H, NK, V, WVN
   USE PanelMesh_mod, ONLY: ISX, ISY
   USE INFG3D_Open, ONLY: INFGREEN3D
   USE FinGreen3D_Open, ONLY: FINGREEN3D
   USE SingularIntgr, ONLY: IMAGE_TERM_EXACT, SGLINTBD_QUAD2, SGLINTBD_TRI2
   USE SingularIntgrMulti, ONLY: IMAGE_TERM_EXACT_MULTI, SGLINTBD_MULTI_QUAD2, SGLINTBD_MULTI_TRI2
   USE SecondOrderTypes, ONLY: panel_mesh, first_order_state, surface_geom
   USE FirstOrderState, ONLY: image_signs, incident_phi, wave_gravity
   USE NearFieldDrift, ONLY: waterline_point
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: field_influence_row, field_potential, waterline_rows, waterline_field, field_point_inside_body
   PUBLIC :: rankine_panel_integral

   REAL(8), PARAMETER :: PI = 3.141592653589793D0
   INTEGER, PARAMETER :: IRR_HULL = 1          ! source panel on the hull, in SGLINTBD_*2

CONTAINS

   SUBROUTINE field_influence_row(mesh, xp, is_multi, row_g, row_dg, x_decide)
      ! G and dG/dn at xp, integrated over each base panel for each image. The panels stay in place and xp is mirrored, as in
      ! CalPotential. is_multi selects the panel integrals that read the combined multi-body arrays.
      ! A panel within 50 panel sizes of the point gets panel integrals of 1/r and of its image, and one further away the
      ! centroid value. x_decide, if present, is the point of that decision instead of xp. The points of one difference
      ! stencil pass its centre, so that they all take the same branch: otherwise the jump between the two would enter the
      ! difference, divided by the step (design spec, section 4.1). For the same reason they also take the series that
      ! FINGREEN3D picks at the centre (RR_DECIDE): the series do not agree exactly at the switches (FINGREEN3D_R).
      !   row_g, row_dg   (nbase, nsys)
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: xp(3)
      LOGICAL, INTENT(IN) :: is_multi
      COMPLEX(8), INTENT(OUT) :: row_g(:, :), row_dg(:, :)
      REAL(8), INTENT(IN), OPTIONAL :: x_decide(3)
      REAL(8) :: xq(3), xf(3), xd(3), xdi(3), nrm(3), area, rkn(4), rr_decide
      COMPLEX(8) :: grn(4)
      INTEGER :: jel, is, flag

      xd = xp
      IF (PRESENT(x_decide)) xd = x_decide

      !$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(jel, is, flag, xq, xf, xdi, nrm, area, rkn, grn, rr_decide) &
      !$OMP SCHEDULE(DYNAMIC, 16)
      DO jel = 1, mesh%nbase
         xq = mesh%xc(:, jel)
         nrm = mesh%nrm(:, jel)
         area = mesh%area(jel)
         DO is = 1, mesh%nsys
            xf = image_signs(is, ISX, ISY)*xp
            xdi = image_signs(is, ISX, ISY)*xd
            flag = MERGE(1, 0, NORM2(xdi - xq) <= 50D0*mesh%psize(jel))

            IF (H < 0D0) THEN
               CALL INFGREEN3D(xq(1), xf(1), xq(2), xf(2), xq(3), xf(3), V, grn, flag)
            ELSE IF (PRESENT(x_decide)) THEN
               rr_decide = SQRT((xq(1) - xdi(1))**2 + (xq(2) - xdi(2))**2)   ! as RR of FINGREEN3D for xf
               CALL FINGREEN3D(xq(1), xf(1), xq(2), xf(2), xq(3), xf(3), V, WVN, NK, H, grn, flag, RR_DECIDE=rr_decide)
            ELSE
               CALL FINGREEN3D(xq(1), xf(1), xq(2), xf(2), xq(3), xf(3), V, WVN, NK, H, grn, flag)
            END IF

            rkn = 0D0
            IF (flag == 1) THEN
               CALL rankine_panel_integral(xf, jel, mesh%ncn(jel), is_multi, rkn)
               IF (is_multi) THEN
                  CALL IMAGE_TERM_EXACT_MULTI(xf, xq, V, jel, IRR_HULL, grn, rkn(1), rkn(2:4))
               ELSE
                  CALL IMAGE_TERM_EXACT(xf, xq, V, jel, IRR_HULL, grn, rkn(1), rkn(2:4))
               END IF
            END IF

            row_g(jel, is) = rkn(1) + grn(1)*area
            row_dg(jel, is) = SUM((rkn(2:4) + grn(2:4)*area)*nrm)
         END DO
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE field_influence_row

   SUBROUTINE rankine_panel_integral(xf, jel, ncn, is_multi, rk)
      ! The integral of 1/r over hull panel jel at the point xf in rk(1), and of its gradient with respect to the panel in rk(2:4).
      ! is_multi selects the routines that read the combined multi-body arrays.
      REAL(8), INTENT(IN) :: xf(3)
      INTEGER, INTENT(IN) :: jel, ncn
      LOGICAL, INTENT(IN) :: is_multi
      REAL(8), INTENT(OUT) :: rk(4)

      IF (is_multi) THEN
         IF (ncn == 3) THEN
            CALL SGLINTBD_MULTI_TRI2(xf, jel, rk, IRR_HULL)
         ELSE
            CALL SGLINTBD_MULTI_QUAD2(xf, jel, rk, IRR_HULL)
         END IF
      ELSE
         IF (ncn == 3) THEN
            CALL SGLINTBD_TRI2(xf, jel, rk, IRR_HULL)
         ELSE
            CALL SGLINTBD_QUAD2(xf, jel, rk, IRR_HULL)
         END IF
      END IF
   END SUBROUTINE rankine_panel_integral

   PURE ELEMENTAL LOGICAL FUNCTION field_point_inside_body(wl_nz)
      ! True if the field point of a segment lies inside the body. The point is d out of the body along the horizontal normal and
      ! d down. The unit panel normal (n_h, n_z) points into the body, so the point is a distance -d*(|n_h| + n_z) along it. That
      ! is negative, which is the fluid side, when n_z > -|n_h| = -sqrt(1 - n_z**2).
      REAL(8), INTENT(IN) :: wl_nz   ! n_z of the unit normal of the waterline panel
      field_point_inside_body = wl_nz < -SQRT(MAX(1D0 - wl_nz**2, 0D0))
   END FUNCTION field_point_inside_body

   SUBROUTINE waterline_rows(mesh, geom, offset, is_multi, rows_g, rows_dg)
      ! The rows of the point of each waterline segment, offset by a fraction of the segment length: out into the fluid along the
      ! horizontal normal, and down.
      !   rows_g, rows_dg   (nbase, nsys, nwl)
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      REAL(8), INTENT(IN) :: offset
      LOGICAL, INTENT(IN) :: is_multi
      COMPLEX(8), INTENT(OUT) :: rows_g(:, :, :), rows_dg(:, :, :)
      REAL(8) :: x(3), distance
      INTEGER :: iw

      DO iw = 1, geom%nwl
         distance = offset*geom%wl_len(iw)
         x = waterline_point(geom, iw) - distance*geom%wl_nh(:, iw)
         x(3) = -distance
         CALL field_influence_row(mesh, x, is_multi, rows_g(:, :, iw), rows_dg(:, :, iw))
      END DO
   END SUBROUTINE waterline_rows

   PURE COMPLEX(8) FUNCTION field_potential(mesh, s, row_g, row_dg)
      ! Perturbation potential of a state at the point of a row.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(first_order_state), INTENT(IN) :: s
      COMPLEX(8), INTENT(IN) :: row_g(:, :), row_dg(:, :)
      INTEGER :: i

      field_potential = (0D0, 0D0)
      DO i = 1, mesh%n
         field_potential = field_potential + row_g(mesh%base(i), mesh%image(i))*s%dphi_p_dn(i) &
                           - row_dg(mesh%base(i), mesh%image(i))*s%phi_p(i)
      END DO
      field_potential = field_potential/(4D0*PI)
   END FUNCTION field_potential

   SUBROUTINE waterline_field(mesh, geom, s, rows_g, rows_dg, phi_wl)
      ! Total potential at the waterline segments, from the rows of waterline_rows.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(first_order_state), INTENT(IN) :: s
      COMPLEX(8), INTENT(IN) :: rows_g(:, :, :), rows_dg(:, :, :)
      COMPLEX(8), INTENT(OUT) :: phi_wl(:)
      REAL(8) :: g_wave
      INTEGER :: iw

      g_wave = wave_gravity(s%omega, s%k, s%h)
      DO iw = 1, geom%nwl
         phi_wl(iw) = incident_phi(waterline_point(geom, iw), s%omega, s%k, s%h, s%amp, s%beta, g_wave) &
                      + field_potential(mesh, s, rows_g(:, :, iw), rows_dg(:, :, iw))
      END DO
   END SUBROUTINE waterline_field

END MODULE WaterlineField
