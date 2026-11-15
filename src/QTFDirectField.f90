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
SUBMODULE (QTFDirect) QTFDirectField_sm
   ! The second-order potential phi_D of the direct solution (Phase 7a, QTFDirect) at field points off the hull (Phase 7c-1),
   ! from its solved panel values by Green's formula.
   !
   ! The 4 pi. Green's second identity over the fluid volume, with G(x, .) of the solver globals at Omega and the normal n
   ! out of the fluid (into the body, as HAMS has it on S):
   !    int_fluid (phi lap G - G lap phi) dV = int_boundary (phi dG/dn - G dphi/dn) dS
   ! lap phi = 0, and lap G = -4 pi delta_x (G = 1/r + ..., lap(1/r) = -4 pi delta), so the left side is -4 pi phi(x) at a
   ! point x in the fluid. On the free surface G and phi_D satisfy the same condition at Omega (the free-surface forcing of
   ! phi_D is not included, as in direct_at_omega), and on the bottom and at infinity they also satisfy the same conditions,
   ! so only S remains:
   !    4 pi phi_D(x) = int_S [ G dphi_D/dn - phi_D dG/dn ] dS
   ! At a collocation point on S the same identity gives 2 pi (the integral equation of HAMS, ASSB_LEFT and ASSB_RBC, see
   ! the header of QTFDirect), and at a point inside the body 0. A point on z = 0 off the hull keeps 4 pi: G holds the
   ! image 1/r1, which coincides with 1/r there, so the singularity 2/r over the half space below z = 0 gives 4 pi too
   ! (WaterlineField has the same free term for the first-order field).
   !
   ! The right side. dphi_D/dn = b on S splits as in direct_at_omega: the panel part vf = B1 + B3 - dphi_I2/dn
   ! (pair_forcing) and B2 = -n.[(X.grad)u], which enters in the Stokes form of QTFDirect (ruling R7-1):
   !    int_S G B2 dS = -r_B2,   r_B2 = int_S G n.[(X.grad)u] dS = b2_stokes_at_point(x)
   ! with fresh panel integrals at x (green_panel_integrals; no row, no lid_row: x is not a stored collocation or lid
   ! point), the same sign as direct_rhs_b2_block, whose right side is -r. G is smooth on S for x off S, so the Stokes form
   ! is no principal value here. So
   !    4 pi phi_D(x) = sum_f [ int_f G vf_f - (int_f dG/dn) phi_D,f ] - r_B2(x)
   ! over the full panels f, with vf and phi_D constant per panel, as the solve has them.
   !
   ! The image sums. Full panel f is image is of base panel j (mesh%base, mesh%image). With M_is the reflection of image
   ! is (image_signs), G(M x, M xi) = G(x, xi), so int_{M q} G(x, xi) dS = int_q G(M x, xi) dS, and the normal of the image
   ! panel is M n, so int_{M q} dG/dn = int_q dG/dn (M x, .) too. field_influence_row gives these as row_g(j, is) and
   ! row_dg(j, is), with the panels in place and x mirrored, as CalPotential does. vf has the layout of pair_forcing,
   ! vf(j, col, is) on image is of base panel j, so
   !    4 pi phi_D(x) = sum_is sum_j [ row_g(j, is) vf(j, col, is) - row_dg(j, is) phi_D((is - 1) nbase + j) ] - r_B2(x)
   ! written below over the full panels with mesh%base and mesh%image, as field_potential of WaterlineField.
   !
   ! Both parts are taken in field mode: the rows (field_mode of field_influence_row) and the B2 part (field_mode of
   ! b2_stokes_at_point, which passes it to green_panel_integrals and green_line_integral). In finite depth FINGREEN3D then
   ! sums its series without Limes, whose isolated outliers (FinGreen3D.f90, FINGREEN3D_R) would otherwise land on single
   ! field points. The first-order field at the same points (fields_green, as SecondOrderDriver takes it for the field
   ! points) uses the same route, so phi_D and the first-order field share one Green function. The hull solve (CALGREEN,
   ! and the B2 right side of direct_at_omega) keeps the default route; at depth 10 on the truncated cylinder the routes
   ! differ by 1.1e-5 relative in the rows and 2.1e-7 in the B2 part (measured, tests direct_field_radiation_identity and
   ! direct_field_b2_field_route), far below the mesh error. In deep water the routes are the same.
   !
   ! The lid of irregular frequency removal (IRSP = 1) is inside the body: Green's formula over S alone holds at a point
   ! outside, so the lid takes no part here.
   !
   ! Omega < 0. direct_at_omega solves the pairs with Omega < 0 conjugated (header of QTFDirect): vf, a, c, e, cl and the
   ! solution are those of |Omega| before the final conjugation. This routine takes the columns as solved and returns the
   ! field at |Omega|; the caller conjugates the columns of Omega < 0 afterwards, as direct_at_omega does for phi_D.

   ! A submodule of QTFDirect: direct_at_omega calls it (the field_points hook) and it calls b2_stokes_at_point of QTFDirect,
   ! so as a separate module the two would use each other. QTFDirect holds the interface and exports direct_field_at_points;
   ! NTHREAD, panel_mesh, surface_geom and b2_stokes_at_point come by host association.
   USE WaterlineField, ONLY: field_influence_row
   IMPLICIT NONE

   REAL(8), PARAMETER :: PI = 3.141592653589793D0

CONTAINS

   MODULE SUBROUTINE direct_field_at_points(mesh, geom, xp, phi_d, vf, a, c, e, cl, is_multi, phi_field)
      ! phi_D at the points xp off the hull for each column (pair) of the solution (header):
      !   xp         (3, npt), points in the fluid, full coordinates
      !   phi_d      (n, nc), the solved phi_D on the full panels, as direct_solve returns it
      !   vf         (nbase, nc, nsys), the panel forcing B1 + B3 - dphi_I2/dn in the layout of pair_forcing
      !   a, c, e    (n, nc), (3, n, nc), (nwl, nc), the B2 densities of b2_stokes_densities
      !   cl         (3, 4, n, nc), their linear parts
      !   phi_field  (npt, nc), out
      ! All at the solver globals of Omega, as direct_at_omega requires; no conjugation here (header).
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      REAL(8), INTENT(IN) :: xp(:, :)
      COMPLEX(8), INTENT(IN) :: phi_d(:, :), vf(:, :, :), a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :)
      LOGICAL, INTENT(IN) :: is_multi
      COMPLEX(8), INTENT(OUT) :: phi_field(:, :)
      COMPLEX(8), ALLOCATABLE :: row_g(:, :), row_dg(:, :), r(:)
      COMPLEX(8) :: acc
      INTEGER :: npt, nc, ip, col, f, j, is

      npt = SIZE(xp, 2)
      nc = SIZE(phi_d, 2)
      IF (SIZE(xp, 1) /= 3 .OR. SIZE(phi_d, 1) /= mesh%n .OR. SIZE(vf, 1) /= mesh%nbase .OR. SIZE(vf, 2) /= nc .OR. &
          SIZE(vf, 3) /= mesh%nsys .OR. SIZE(a, 1) /= mesh%n .OR. SIZE(a, 2) /= nc .OR. SIZE(c, 1) /= 3 .OR. &
          SIZE(c, 2) /= mesh%n .OR. SIZE(c, 3) /= nc .OR. SIZE(e, 1) /= geom%nwl .OR. SIZE(e, 2) /= nc .OR. &
          SIZE(cl, 1) /= 3 .OR. SIZE(cl, 2) /= 4 .OR. SIZE(cl, 3) /= mesh%n .OR. SIZE(cl, 4) /= nc .OR. &
          SIZE(phi_field, 1) /= npt .OR. SIZE(phi_field, 2) /= nc) THEN
         WRITE (*, '(a)') ' direct_field_at_points: xp, phi_d, vf, a, c, e, cl or phi_field does not have the size of the'// &
            ' mesh and of each other.'
         ERROR STOP 1
      END IF
      IF (npt == 0 .OR. nc == 0) RETURN

      ! The points are independent and write only their own row of phi_field, so the result does not depend on the thread
      ! count. The loop over the panels in field_influence_row is itself an OpenMP loop, which runs on one thread inside
      ! this region (nested parallelism is off by default), as in fields_green and direct_rhs_b2_block.
      !$OMP PARALLEL NUM_THREADS(NTHREAD) PRIVATE(row_g, row_dg, r, acc, ip, col, f, j, is)
      ALLOCATE (row_g(mesh%nbase, mesh%nsys), row_dg(mesh%nbase, mesh%nsys), r(nc))
      !$OMP DO SCHEDULE(DYNAMIC, 1)
      DO ip = 1, npt
         CALL field_influence_row(mesh, xp(:, ip), is_multi, row_g, row_dg, field_mode=.TRUE.)
         CALL b2_stokes_at_point(mesh, geom, xp(:, ip), a, c, e, is_multi, r, cl, field_mode=.TRUE.)
         DO col = 1, nc
            acc = (0D0, 0D0)
            DO f = 1, mesh%n
               j = mesh%base(f)
               is = mesh%image(f)
               acc = acc + row_g(j, is)*vf(j, col, is) - row_dg(j, is)*phi_d(f, col)
            END DO
            phi_field(ip, col) = (acc - r(col))/(4D0*PI)
         END DO
      END DO
      !$OMP END DO
      DEALLOCATE (row_g, row_dg, r)
      !$OMP END PARALLEL
   END SUBROUTINE direct_field_at_points

END SUBMODULE QTFDirectField_sm
