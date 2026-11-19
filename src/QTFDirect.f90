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
   !
   ! The solve at Omega (direct_at_omega). The matrix of the phi_D problem is AMAT of the radiation problem at Omega
   ! (AMAT_MULTI for several bodies). With ISOLV = 1 ASSB_LEFT leaves it as its LU factors in place (ZGETRF), and
   ! RADIATION_SOLVER and DIFFRACTION_SOLVER call ZGETRS on them. With ISOLV = 2 ASSB_LEFT leaves AMAT unfactored for
   ! GMRES, and direct_factor factors it in place once per Omega: the direct solve is the last use of AMAT at that Omega.
   ! The right sides of all the pairs of one Omega are then one block: ZGETRS with NRHS = number of pairs, per symmetry
   ! block, as RADIATION_SOLVER. For each pair (s1, s2), Omega = omega_1 - omega_2:
   !    b = B1 + B3 - dphi_I2/dn       panel values, through ASSB_RHS_FIELD
   !      + B2                          through direct_rhs_b2_block (the Stokes form with G, ruling R7-1)
   !    T_dir,k = i Omega rho int_S (phi_I2 + phi_D) n_k dS
   ! with B1, B3, phi_I2 and the n_k (generalized normals of the body of mode k) taken as potential_pair takes them, so
   ! that T_dir equals T_pot of potential_pair (qtf_b2 = body) by Green's theorem, to the mesh error. A pair with Omega < 0
   ! is treated as potential_pair treats it: its phi_D solves the problem at -|Omega|, whose matrix and Green function are
   ! the conjugates of those at |Omega|, so phi_D = conj(A^-1 r[conj b]), with r[.] the right side at |Omega|.
   !
   ! The free-surface forcing (qtf_terms = full, Phase 7b). With the optional fs of direct_at_omega the right side gets
   ! r_FS = (s_D / g) int_{z=0} G Q_D dS (QTFDirectFS: s_D, the inner Green matrix and the Graf amplitudes of G outside the
   ! circle R) at every row, hull and lid. Per pair: Q_D at the inner points (fs_forcing_points) and the outer moments
   ! (fs_outer_moments) of QTFFreeSurface, from the 6a store that fs points to (fields at the inner points, incident fields
   ! there, partial waves of the scattered potentials outside R). Omega < 0: the right side at |Omega| is r[conj b], and
   ! conj Q_D(1, 2) = Q_D(2, 1) (pair_product of KERNEL_DIFF; fs_pair_terms of QTFFreeSurface has the derivation), so the
   ! pair is evaluated as (2, 1), which is ordered (omega_2 > omega_1) as fs_outer_moments needs. The moments of (2, 1)
   ! are the integrals of the B_mn at |Omega| against Q_D(2, 1) = conj Q_D(1, 2): no conjugate of the moments, since the
   ! B_mn are those of G at |Omega|, not at -|Omega|. direct_fs_block (submodule QTFDirectFSBlock_sm: it calls
   ! QTFDirectFS, which uses this module) does it. The term T_dir,k = i Omega rho int_S phi_D n_k dS then holds the
   ! free-surface part, which the indirect method has as FS_inner + FS_outer (fs_pair_terms) by Green's theorem: test
   ! direct_fs_equals_indirect_full pins s_D with it. Without fs nothing of this runs.
   !
   ! Irregular frequency removal (IRSP = 1, Phase 7a Task 3b). HAMS adds rows at the lid points x_l on the interior
   ! waterplane, where the same integral equation holds with 0 in place of 2 pi phi (ASSB_LEFT_IRR, ASSB_RBC_IRR): the
   ! system A x = r has nbase + nlid rows and nbase unknowns, and is solved in the least-squares form of HAMS,
   ! (A^T A) x = A^T r, with CMAT = A^T A factored by ASSB_LEFT_IRR (ISOLV = 1) and A^T r formed as ASSB_RBC_IRR forms
   ! DRMAT (plain transpose, no conjugate). The direct right side follows the same construction with b in place of n_k:
   !    r_l = int_S G(x_l, xi) b(xi) dS_xi
   ! the field integral from the lid points, with the stored lid integrals DGRN and PKBN (ASSB_RHS_FIELD_IRR), and the
   ! B2 part by the same Stokes form with G(x_l, .), which is smooth on S since x_l is off the hull: the stored lid
   ! integrals in place of CGRN and RKBN, and the near-panel linear density and the waterline line integral as on the
   ! hull rows. AMAT keeps the full system; it is not factored with the lid. The Omega < 0 rule above holds unchanged,
   ! since A^T uses no conjugate: conj((A^T A)^-1 A^T r[conj b]) solves the problem at -|Omega|.
   !
   ! Sum frequency (Phase 7c-2). With kernels(p) = KERNEL_SUM, pair p of direct_at_omega is the sum pair at
   ! Omega+ = omega_1 + omega_2: B1, B2 (b2_stokes_densities) and B3 by pair_plus, phi_I2 of the sum, B3 with -i Omega+
   ! (QTFPotential has the derivation), and T_dir,k = i Omega+ rho int_S (phi_I2 + phi_D) n_k dS. Omega+ > 0, so the
   ! forcing is never conjugated: phi_D solves the problem at Omega+ with the matrix and Green function of the radiation
   ! solve there. Sum and difference pairs of the same Omega share one factorisation and one block of right sides; each
   ! column carries its own signed Omega (dws). By the Green's identity of QTFPotential at Omega+, T_dir = T_pot of
   ! potential_pair with KERNEL_SUM to the mesh error (test direct_sum_matches_indirect_pair). With fs, a sum column gets
   ! the free-surface forcing Q_S of pair_plus (direct_fs_block: the pair ordered by frequency, never conjugated), and
   ! T_dir holds the free-surface part that the indirect method has as FS_inner + FS_outer of fs_pair_terms with
   ! KERNEL_SUM (test direct_fs_sum_equals_indirect); a block may mix sum and difference columns (direct_fs_mixed_block).

   USE Const_mod, ONLY: RXY, RXY_MULTI
   USE HAMS_mod, ONLY: NTHREAD, ISOLV, IRSP
   USE WaveDyn_mod, ONLY: H, NK, V, WVN, W1
   USE PanelMesh_mod, ONLY: ISX, ISY, DS, PNSZ, XYZ_P, DS_MULTI_COMB, PNSZ_MULTI_COMB, XYZ_GLOBAL_MULTI_COMB_P, IPIV, &
                            IPIV_MULTI_COMB, iNELEM_TOTAL, iXYZ_GLOBAL_MULTI_COMB_P
   USE Inerfs_mod, ONLY: iNELEM, iXYZ_P
   USE SingularIntgr, ONLY: WAVE_PART_PANEL
   USE SingularIntgrMulti, ONLY: WAVE_PART_PANEL_MULTI
   USE LinearMatrix_mod, ONLY: AMAT, AMAT_MULTI, CMAT, CMAT_MULTI
   USE Potentials_mod, ONLY: CGRN, RKBN, CGRN_MULTI_COMB, RKBN_MULTI_COMB, DGRN, PKBN, DGRN_MULTI_COMB, PKBN_MULTI_COMB
   USE AssbMatx, ONLY: ASSB_RHS_FIELD_IRR
   USE AssbMatxMulti, ONLY: ASSB_RHS_FIELD_MULTI_IRR
   USE SurfaceGeometry, ONLY: body_sizes, merge_nodes
   USE SharpEdge, ONLY: DEFAULT_EDGE_DEG
   USE INFG3D_Open, ONLY: INFGREEN3D
   USE FinGreen3D_Open, ONLY: FINGREEN3D
   USE SecondOrderTypes, ONLY: panel_mesh, surface_geom, first_order_state
   USE FirstOrderState, ONLY: image_signs, wave_gravity
   USE NearFieldDrift, ONLY: pair_cross, rigid_displacement, waterline_point, KERNEL_DIFF, KERNEL_SUM
   USE WaterlineField, ONLY: rankine_panel_integral
   USE QTFPotential, ONLY: waterline_velocity
   USE SecondOrderIncident, ONLY: incident2, incident2_pair, incident2_phi, incident2_grad
   USE SurfaceVelocity, ONLY: tangential_fit
   USE PartialWaves, ONLY: pw_set
   USE QTFFreeSurface, ONLY: fs_point_field
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: direct_rhs_b2, direct_rhs_b2_block, b2_stokes_densities, b2_stokes_at_point
   PUBLIC :: green_panel_integrals, stored_panel_integrals, green_line_integral, rankine_moments, sharp_edge_panels
   PUBLIC :: direct_setup, direct_edge_mask, direct_factor, direct_solve, direct_at_omega, direct_lid_points, &
             direct_field_at_points
   PUBLIC :: direct_fs_data

   TYPE :: direct_fs_data
      ! The free-surface part of the direct right side at one Omega (QTFDirectFS: direct_fs_setup fills the first block,
      ! the caller of direct_at_omega the pointers). Defined here, not in QTFDirectFS, so that direct_at_omega can take it:
      ! QTFDirectFS uses this module.
      !   gfs(nrow, nfs)          the inner Green matrix, (1 / w_j) int_{panel j} G(x_i, .) dS; rows the full hull panels,
      !                           then the full lid points
      !   a(-mg:mg, 0:nev, nrow)  the Graf amplitudes of G(x_i, .) outside R in the contraction form (QTFDirectFS header);
      !                           zero for |m| > mcut
      !   w(nfs)                  the inner panel areas
      !   kev(nev)                the evanescent wave numbers at Omega, k0 the propagating one
      !   cut_error               q^(mcut + 1) / ((mcut + 1) (1 - q)), q = r_max / R, when mcut < mg, else 0
      ! The 6a store of the pairs (QTFDriver), indexed like the store of direct_at_omega (frequency, heading):
      !   f(nfs, nfreq, nbeta)    the total first-order fields at the inner points (fs_inner_fields)
      !   inc(nfs, nfreq, nbeta)  the incident fields there (incident_at, phi_zz = k^2 phi)
      !   pw(nfreq, nbeta)        the partial waves of the scattered potentials outside R (fs_circle_modes)
      !   np                      the chords of the polygon of a generated inner mesh (the sliver), 0 for a file mesh
      INTEGER :: nrow = 0, nfs = 0, nbase = 0, nl = 0, nsys = 1, n = 0, mg = 0, nev = 0, mcut = -1
      LOGICAL :: multi = .FALSE.
      REAL(8) :: omega = 0D0, g = 0D0, h = 0D0, radius = 0D0, centre(2) = 0D0, k0 = 0D0, r_max = 0D0, cut_error = 0D0
      REAL(8), ALLOCATABLE :: kev(:), w(:)
      COMPLEX(8), ALLOCATABLE :: gfs(:, :), a(:, :, :)
      TYPE(fs_point_field), POINTER :: f(:, :, :) => NULL(), inc(:, :, :) => NULL()
      TYPE(pw_set), POINTER :: pw(:, :) => NULL()
      INTEGER :: np = 0
   END TYPE direct_fs_data

   INTERFACE
      MODULE SUBROUTINE direct_fs_block(fs, store, pairs, dws, kernels, rhs, seconds)
         ! The free-surface part of the right sides of a block of pairs (QTFDirectFSBlock.f90, a submodule of this module,
         ! which it needs for direct_fs_rhs of QTFDirectFS): added to rhs (nbase + nl, npair, nsys), the layout of
         ! ASSB_RHS_FIELD_IRR. pairs (4, npair), dws (npair) and kernels (npair) as direct_at_omega has them: dws is
         ! omega_1 - omega_2 of a difference column and omega_1 + omega_2 of a sum column. The difference columns of
         ! dws < 0 get the right side at |Omega| of the conjugated forcing (header); a sum column is never conjugated and
         ! its pair is ordered by frequency, not by the sign of dws (QTFDirectFSBlock header). seconds(2): the time of
         ! Q_D and the outer moments, and of the contraction (direct_fs_rhs), added to.
         TYPE(direct_fs_data), INTENT(IN) :: fs
         TYPE(first_order_state), INTENT(IN) :: store(:, :)
         INTEGER, INTENT(IN) :: pairs(:, :)
         REAL(8), INTENT(IN) :: dws(:)
         INTEGER, INTENT(IN) :: kernels(:)
         COMPLEX(8), INTENT(INOUT) :: rhs(:, :, :)
         REAL(8), INTENT(INOUT) :: seconds(2)
      END SUBROUTINE direct_fs_block
   END INTERFACE

   COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
   ! As CALGREEN: a panel closer than this many panel sizes gets panel integrals of 1/r and of its image
   REAL(8), PARAMETER :: NEAR_PANEL_SIZES = 50D0
   ! Lid points: the wave part is integrated over source panels closer than this many panel sizes (WNEAR of CALGREEN_IRR)
   REAL(8), PARAMETER :: LID_WAVE_NEAR = 2D0
   ! direct_at_omega takes the pairs of one Omega in blocks of at most this many bytes of work arrays (about 25 complex
   ! values per panel and pair with a lid: the fields, the right sides, the B2 densities and their linear parts, phi_D)
   REAL(8), PARAMETER :: DIRECT_BLOCK_BYTES = 256D6
   ! The sharp-edge mask of the run (direct_setup): panels that keep the one-point rule in the B2 right side
   LOGICAL, ALLOCATABLE :: edge_mask(:)
   ! With ISOLV /= 1, the state of the matrix that direct_factor factored last: W1 and the last diagonal entry of the
   ! factors (direct_factor). A second factorisation of the same factors would give a wrong phi_D with no error.
   REAL(8) :: factored_omega = -1D0
   COMPLEX(8) :: factored_entry = (0D0, 0D0)

   INTERFACE
      MODULE SUBROUTINE direct_field_at_points(mesh, geom, xp, phi_d, vf, a, c, e, cl, is_multi, phi_field)
         ! phi_D at field points off the hull by Green's formula (QTFDirectField.f90, a submodule of this module, which it
         ! needs for b2_stokes_at_point); the hook of direct_at_omega
         TYPE(panel_mesh), INTENT(IN) :: mesh
         TYPE(surface_geom), INTENT(IN) :: geom
         REAL(8), INTENT(IN) :: xp(:, :)
         COMPLEX(8), INTENT(IN) :: phi_d(:, :), vf(:, :, :), a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :)
         LOGICAL, INTENT(IN) :: is_multi
         COMPLEX(8), INTENT(OUT) :: phi_field(:, :)
      END SUBROUTINE direct_field_at_points
   END INTERFACE

CONTAINS

   SUBROUTINE direct_setup(mesh, edge_deg, nmask)
      ! Once per run, before direct_at_omega: the sharp-edge mask of sharp_edge_panels for the B2 right side, with edge_deg
      ! (qtf_edge_angle). nmask, optional: the number of panels that keep the one-point rule.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: edge_deg
      INTEGER, INTENT(OUT), OPTIONAL :: nmask

      ! With IRSP = 1 the lid rows of the right side need lid points: stop before the first-order solves if there are none.
      ! The counts are known here: iNELEM, and iNELEM_TOTAL for several bodies (HAMS_Prog sums it while it reads the
      ! meshes, before so_setup_multi). direct_at_omega checks again, for callers without the set-up.
      CALL check_lid(mesh%nbody > 1)
      IF (ALLOCATED(edge_mask)) DEALLOCATE (edge_mask)
      CALL sharp_edge_panels(mesh, edge_deg, edge_mask)
      IF (PRESENT(nmask)) nmask = COUNT(edge_mask)
   END SUBROUTINE direct_setup

   SUBROUTINE check_lid(is_multi)
      ! Stops with a message if IRSP = 1 and there are no lid points: the solve would run without the removal
      LOGICAL, INTENT(IN) :: is_multi

      IF (IRSP == 1 .AND. direct_lid_points(is_multi) == 0) THEN
         WRITE (*, '(a)') ' QTF direct: If_remove_irr_freq = 1, but the mesh has no lid points (WaterPlaneMesh.pnl is'// &
            ' empty or has no panels). Add the lid panels, or set If_remove_irr_freq = 0.'
         ERROR STOP 1
      END IF
   END SUBROUTINE check_lid

   FUNCTION direct_edge_mask() RESULT(mask)
      ! The mask of direct_setup (empty before it)
      LOGICAL, ALLOCATABLE :: mask(:)

      IF (ALLOCATED(edge_mask)) THEN
         mask = edge_mask
      ELSE
         ALLOCATE (mask(0))
      END IF
   END FUNCTION direct_edge_mask

   INTEGER FUNCTION direct_lid_points(is_multi)
      ! The number of lid points of irregular frequency removal: iNELEM (iNELEM_TOTAL for several bodies) with IRSP = 1,
      ! else 0. The right sides of direct_solve have this many rows after the nbase hull rows.
      LOGICAL, INTENT(IN) :: is_multi

      direct_lid_points = 0
      IF (IRSP /= 1) RETURN
      IF (is_multi) THEN
         IF (ALLOCATED(iNELEM_TOTAL)) direct_lid_points = iNELEM_TOTAL
      ELSE
         direct_lid_points = iNELEM
      END IF
   END FUNCTION direct_lid_points

   SUBROUTINE direct_factor(mesh, is_multi)
      ! With ISOLV /= 1, the LU factors in place of the matrix of the direct solve, with the pivots in IPIV
      ! (IPIV_MULTI_COMB): AMAT (AMAT_MULTI), or with the lid (IRSP = 1) the least-squares matrix CMAT (CMAT_MULTI). With
      ! ISOLV = 1 ASSB_LEFT (ASSB_LEFT_IRR) has factored it already, and nothing is done. Once per Omega, after the
      ! radiation solve there. Stops if the matrix is the one it factored last (same W1 and the same last diagonal entry
      ! as its factors left it): its factors would be factored again. Stops on a failed or singular factorisation.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      LOGICAL, INTENT(IN), OPTIONAL :: is_multi
      INTEGER :: ip, nb, info
      LOGICAL :: multi, lid

      IF (ISOLV == 1) RETURN
      multi = mesh%nbody > 1
      IF (PRESENT(is_multi)) multi = is_multi
      lid = direct_lid_points(multi) > 0
      nb = mesh%nbase
      IF (ABS(W1 - factored_omega) <= 1D-12*W1) THEN
         IF (last_entry() == factored_entry) THEN
            WRITE (*, '(a,ES14.6,a)') ' QTF direct: the matrix at Omega =', W1, &
               ' rad/s is factored already; direct_factor must run once per radiation solve.'
            ERROR STOP 1
         END IF
      END IF
      DO ip = 1, mesh%nsys
         IF (multi .AND. lid) THEN
            CALL ZGETRF(nb, nb, CMAT_MULTI(:, :, ip), nb, IPIV_MULTI_COMB(:, ip), info)
         ELSE IF (multi) THEN
            CALL ZGETRF(nb, nb, AMAT_MULTI(:, :, ip), nb, IPIV_MULTI_COMB(:, ip), info)
         ELSE IF (lid) THEN
            CALL ZGETRF(nb, nb, CMAT(:, :, ip), nb, IPIV(:, ip), info)
         ELSE
            CALL ZGETRF(nb, nb, AMAT(:, :, ip), nb, IPIV(:, ip), info)
         END IF
         IF (info /= 0) THEN
            WRITE (*, '(a,I0,a,ES14.6,a,I0,a)') ' QTF direct: the factorisation of block ', ip, ' at Omega =', W1, &
               ' rad/s failed, ZGETRF INFO = ', info, MERGE(' (singular matrix) ', ' (illegal argument)', info > 0)
            ERROR STOP 1
         END IF
      END DO
      factored_omega = W1
      factored_entry = last_entry()

   CONTAINS

      COMPLEX(8) FUNCTION last_entry()
         ! The last diagonal entry of the last block of the matrix that is factored
         INTEGER :: ns
         ns = mesh%nsys
         IF (multi .AND. lid) THEN
            last_entry = CMAT_MULTI(nb, nb, ns)
         ELSE IF (multi) THEN
            last_entry = AMAT_MULTI(nb, nb, ns)
         ELSE IF (lid) THEN
            last_entry = CMAT(nb, nb, ns)
         ELSE
            last_entry = AMAT(nb, nb, ns)
         END IF
      END FUNCTION last_entry

   END SUBROUTINE direct_factor

   SUBROUTINE direct_solve(mesh, rhs, phi, is_multi)
      ! phi_D for nrhs right sides at once, with the LU factors at Omega:
      !   rhs  (nbase + nlid, nrhs, nsys), the block layout of ASSB_RHS_FIELD_IRR, nlid = direct_lid_points (0 without the
      !        lid); overwritten: rows 1 to nbase hold the solution of each block
      !   phi  (n, nrhs), phi_D on the full panels: image a of base panel i is sum_ip RXY(a, ip) x(i, :, ip)/nsys, as
      !        RADIATION_SOLVER builds MXPOT and radiation_psi expands it (RXY_MULTI for several bodies)
      ! Without the lid, ZGETRS on the factors of AMAT (AMAT_MULTI). With the lid, the least-squares form of HAMS: A^H r with
      ! A = AMAT (AMAT_MULTI) of nbase + nlid rows, as NORMAL_RHS forms DRMAT, then ZGETRS on the factors of CMAT
      ! (CMAT_MULTI) = A^H A (NORMAL_MATRIX). Stops if ZGETRS reports an error.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      COMPLEX(8), INTENT(INOUT) :: rhs(:, :, :)
      COMPLEX(8), INTENT(OUT) :: phi(:, :)
      LOGICAL, INTENT(IN), OPTIONAL :: is_multi
      COMPLEX(8), ALLOCATABLE :: d(:, :)
      INTEGER :: ip, nb, nt, nrhs, info, f, im
      LOGICAL :: multi

      multi = mesh%nbody > 1
      IF (PRESENT(is_multi)) multi = is_multi
      nb = mesh%nbase
      nt = nb + direct_lid_points(multi)
      nrhs = SIZE(rhs, 2)
      IF (SIZE(rhs, 1) /= nt .OR. SIZE(rhs, 3) /= mesh%nsys .OR. SIZE(phi, 1) /= mesh%n .OR. SIZE(phi, 2) /= nrhs) THEN
         WRITE (*, '(a)') ' direct_solve: rhs or phi does not have the size of the mesh (and the lid).'
         ERROR STOP 1
      END IF
      IF (nrhs == 0) RETURN
      IF (nt > nb) THEN
         IF ((multi .AND. SIZE(AMAT_MULTI, 1) /= nt) .OR. (.NOT. multi .AND. SIZE(AMAT, 1) /= nt)) THEN
            WRITE (*, '(a)') ' direct_solve: AMAT does not have the hull and lid rows of the right side.'
            ERROR STOP 1
         END IF
         ALLOCATE (d(nb, nrhs))
      END IF
      DO ip = 1, mesh%nsys
         IF (nt > nb) THEN
            ! d = A^H r over the hull and lid rows (the first nb columns of AMAT), as DRMAT of NORMAL_RHS
            IF (multi) THEN
               CALL ZGEMM('C', 'N', nb, nrhs, nt, (1D0, 0D0), AMAT_MULTI(:, :, ip), nt, rhs(:, :, ip), nt, (0D0, 0D0), d, nb)
               CALL ZGETRS('No transpose', nb, nrhs, CMAT_MULTI(:, :, ip), nb, IPIV_MULTI_COMB(:, ip), d, nb, info)
            ELSE
               CALL ZGEMM('C', 'N', nb, nrhs, nt, (1D0, 0D0), AMAT(:, :, ip), nt, rhs(:, :, ip), nt, (0D0, 0D0), d, nb)
               CALL ZGETRS('No transpose', nb, nrhs, CMAT(:, :, ip), nb, IPIV(:, ip), d, nb, info)
            END IF
            rhs(1:nb, :, ip) = d
         ELSE IF (multi) THEN
            CALL ZGETRS('No transpose', nb, nrhs, AMAT_MULTI(:, :, ip), nb, IPIV_MULTI_COMB(:, ip), rhs(:, :, ip), nb, info)
         ELSE
            CALL ZGETRS('No transpose', nb, nrhs, AMAT(:, :, ip), nb, IPIV(:, ip), rhs(:, :, ip), nb, info)
         END IF
         IF (info /= 0) THEN
            WRITE (*, '(a,I0,a,ES14.6,a,I0)') ' QTF direct: the solve of block ', ip, ' at Omega =', W1, &
               ' rad/s failed, ZGETRS INFO = ', info
            ERROR STOP 1
         END IF
      END DO
      DO f = 1, mesh%n
         im = mesh%image(f)
         phi(f, :) = (0D0, 0D0)
         DO ip = 1, mesh%nsys
            IF (multi) THEN
               phi(f, :) = phi(f, :) + RXY_MULTI(im, ip)*rhs(mesh%base(f), :, ip)
            ELSE
               phi(f, :) = phi(f, :) + RXY(im, ip)*rhs(mesh%base(f), :, ip)
            END IF
         END DO
         phi(f, :) = phi(f, :)/mesh%nsys
      END DO
   END SUBROUTINE direct_solve

   SUBROUTINE direct_at_omega(mesh, geom, omega, store, store_wl, pairs, tdir, is_multi, omega_tol, phi2_constant, seconds, parts, &
                              fs, field_points, phi_field, kernels, factored)
      ! The direct load T_dir of every pair of the second-order frequency omega (> 0), per unit amplitudes (header):
      !   store, store_wl  the first-order states (nfreq, nbeta) and their waterline potentials (nwl, nfreq, nbeta), as
      !                    QTFDriver keeps them
      !   pairs            (4, npair): frequency and heading of state 1, then of state 2, in store; each difference
      !                    pair with |omega_1 - omega_2| = omega, each sum pair with omega_1 + omega_2 = omega, within
      !                    omega_tol (default 1e-6 max(omega_1, omega_2))
      !   kernels          optional (npair): KERNEL_DIFF or KERNEL_SUM per pair (default all KERNEL_DIFF, bit-identical
      !                    to the code without it). A sum pair is not conjugated (header), with fs as without.
      !   tdir             (6 nbody, npair), T_dir,k of each pair, e^{-i omega t}, as the sum of the terms of potential_pair
      !   phi2_constant    optional test hook, as potential_pair: phi_I2 replaced by this constant (dphi_I2/dn = 0)
      !   parts            optional (6 nbody, npair, 2): the part of tdir from phi_I2 (:, :, 1) and from phi_D (:, :, 2);
      !                    their sum is tdir. (6 nbody, npair, 3): the part from phi_D is split into the part of the body
      !                    forcing (:, :, 2) and that of the free-surface forcing (:, :, 3, zero without fs); the three
      !                    sum to tdir. For that split the free-surface right sides are solved as their own columns.
      !   seconds          optional (4): wall time of the pair forcing with ASSB_RHS_FIELD_IRR, of the B2 part (the
      !                    densities of b2_stokes_densities and direct_rhs_b2_block), of the factorisation, solve and
      !                    load, and of the field evaluation (the hook below; 0 without it), for the log
      !   fs               optional: the free-surface forcing (qtf_terms = full, header), direct_fs_setup at omega with
      !                    the pointers to the 6a store set; its time (Q_D and the outer moments, the contraction) is
      !                    written to the log per call
      !   field_points     optional (3, npt), with phi_field (npt, npair) out, both or neither: phi_D of every pair at the
      !                    points off the hull by Green's formula (direct_field_at_points), e^{-i omega t}. phi_field is
      !                    the phi_D of the body forcing only: with fs the free-surface columns are then solved apart
      !                    (as for the third slice of parts), since Green's formula over the hull alone holds only for
      !                    the body-forced phi_D. tdir and parts still have the whole load. A sum column's phi_field
      !                    is phi_D at Omega+ as solved; only the difference columns with Omega < 0 are conjugated.
      ! The solver globals must be at omega, with CALGREEN (CALGREEN_MULTI) and ASSB_LEFT (ASSB_LEFT_MULTI) run there, as
      ! RadiationOnlySolves leaves them after the radiation solve at omega; W1 is checked. direct_setup must have been
      ! called for this mesh. With IRSP = 1 the right sides have the lid rows too (header), and the radiation solve must
      ! have been the one with the lid (CALGREEN_IRR, ASSB_LEFT_IRR). The matrix is factored in place when ISOLV /= 1
      ! (direct_factor), so call this once per radiation solve; or with factored = .TRUE. (optional, default .FALSE.) for
      ! a further call at the same Omega after one that factored the matrix (with pairs): direct_factor is then skipped.
      ! QTFDriver solves the sum columns of an Omega in such a call of their own, so that the difference columns keep the
      ! block sizes, and so the bits, of a run without the sum (MKL does not promise column p of ZGETRS or ZGEMM to be
      ! the same for another number of right sides).
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      REAL(8), INTENT(IN) :: omega
      TYPE(first_order_state), INTENT(IN) :: store(:, :)
      COMPLEX(8), INTENT(IN) :: store_wl(:, :, :)
      INTEGER, INTENT(IN) :: pairs(:, :)
      COMPLEX(8), INTENT(OUT) :: tdir(:, :)
      LOGICAL, INTENT(IN), OPTIONAL :: is_multi
      REAL(8), INTENT(IN), OPTIONAL :: omega_tol
      COMPLEX(8), INTENT(IN), OPTIONAL :: phi2_constant
      REAL(8), INTENT(OUT), OPTIONAL :: seconds(4)
      COMPLEX(8), INTENT(OUT), OPTIONAL :: parts(:, :, :)
      TYPE(direct_fs_data), INTENT(IN), OPTIONAL :: fs
      REAL(8), INTENT(IN), OPTIONAL :: field_points(:, :)
      COMPLEX(8), INTENT(OUT), OPTIONAL :: phi_field(:, :)
      INTEGER, INTENT(IN), OPTIONAL :: kernels(:)
      LOGICAL, INTENT(IN), OPTIONAL :: factored
      COMPLEX(8), ALLOCATABLE :: vf(:, :, :), rf(:, :, :), r2(:, :, :), a(:, :), c(:, :, :), cl(:, :, :, :), e(:, :), &
                                 phi(:, :), p2(:, :), rfs(:, :, :), phis(:, :)
      REAL(8), ALLOCATABLE :: dws(:)
      INTEGER, ALLOCATABLE :: kns(:)
      INTEGER :: npair, nb, nt, nl, ns, n, nmode, first, last, nc, col, ncol_max
      INTEGER(8) :: c0, c1, c2, rate
      REAL(8) :: sec(4), fsec(2), per_col
      LOGICAL :: multi, fs_apart, fs_slice

      npair = SIZE(pairs, 2)
      nb = mesh%nbase
      ns = mesh%nsys
      n = mesh%n
      nmode = 6*mesh%nbody
      IF (SIZE(pairs, 1) /= 4 .OR. SIZE(tdir, 1) /= nmode .OR. SIZE(tdir, 2) /= npair .OR. &
          SIZE(store_wl, 1) /= geom%nwl) THEN
         WRITE (*, '(a)') ' direct_at_omega: pairs, tdir or store_wl does not have the size of the mesh.'
         ERROR STOP 1
      END IF
      IF (.NOT. ALLOCATED(edge_mask)) THEN
         WRITE (*, '(a)') ' direct_at_omega: direct_setup has not been called.'
         ERROR STOP 1
      END IF
      IF (SIZE(edge_mask) /= n) THEN
         WRITE (*, '(a)') ' direct_at_omega: direct_setup was called for another mesh.'
         ERROR STOP 1
      END IF
      IF (.NOT. (omega > 0D0) .OR. ABS(W1 - omega) > 1D-12*omega) THEN
         WRITE (*, '(a,ES14.6,a,ES14.6,a)') ' direct_at_omega: Omega =', omega, ' rad/s, but the solver is at', W1, &
            ' rad/s; the matrix and the Green function must be those of Omega.'
         ERROR STOP 1
      END IF
      multi = mesh%nbody > 1
      IF (PRESENT(is_multi)) multi = is_multi
      CALL check_lid(multi)
      nl = direct_lid_points(multi)
      nt = nb + nl
      tdir = (0D0, 0D0)
      sec = 0D0
      IF (PRESENT(seconds)) seconds = 0D0
      fs_apart = .FALSE.
      fs_slice = .FALSE.
      IF (PRESENT(parts)) THEN
         IF (SIZE(parts, 1) /= nmode .OR. SIZE(parts, 2) /= npair .OR. SIZE(parts, 3) < 2 .OR. SIZE(parts, 3) > 3) THEN
            WRITE (*, '(a)') ' direct_at_omega: parts must be (6 nbody, npair, 2) or (6 nbody, npair, 3).'
            ERROR STOP 1
         END IF
         parts = (0D0, 0D0)
         fs_slice = SIZE(parts, 3) == 3
         fs_apart = PRESENT(fs) .AND. fs_slice
      END IF
      ALLOCATE (kns(npair))
      kns = KERNEL_DIFF
      IF (PRESENT(kernels)) THEN
         IF (SIZE(kernels) /= npair) THEN
            WRITE (*, '(a)') ' direct_at_omega: kernels must have one entry per pair.'
            ERROR STOP 1
         END IF
         IF (ANY(kernels /= KERNEL_DIFF .AND. kernels /= KERNEL_SUM)) THEN
            WRITE (*, '(a)') ' direct_at_omega: kernels must be KERNEL_DIFF or KERNEL_SUM.'
            ERROR STOP 1
         END IF
         kns = kernels
      END IF
      IF (PRESENT(fs)) CALL check_fs()
      fsec = 0D0
      IF (PRESENT(field_points) .NEQV. PRESENT(phi_field)) THEN
         WRITE (*, '(a)') ' direct_at_omega: field_points and phi_field must be given together.'
         ERROR STOP 1
      END IF
      IF (PRESENT(phi_field)) THEN
         IF (SIZE(field_points, 1) /= 3 .OR. SIZE(phi_field, 1) /= SIZE(field_points, 2) .OR. SIZE(phi_field, 2) /= npair) THEN
            WRITE (*, '(a)') ' direct_at_omega: field_points must be (3, npt) and phi_field (npt, npair).'
            ERROR STOP 1
         END IF
         phi_field = (0D0, 0D0)
         ! The field hook takes the body-forced phi_D only, so the free-surface columns are solved apart (header)
         IF (PRESENT(fs)) fs_apart = .TRUE.
      END IF
      IF (npair == 0) RETURN

      CALL SYSTEM_CLOCK(c0, rate)
      IF (PRESENT(factored)) THEN
         IF (.NOT. factored) CALL direct_factor(mesh, multi)
      ELSE
         CALL direct_factor(mesh, multi)
      END IF
      CALL SYSTEM_CLOCK(c1)
      sec(3) = sec(3) + REAL(c1 - c0, 8)/REAL(rate, 8)
      ncol_max = MAX(1, INT(DIRECT_BLOCK_BYTES/(16D0*25D0*REAL(n, 8))))
      IF (PRESENT(fs)) THEN
         ! The free-surface columns add Q_D at the inner points, the moments and, apart, a right side and phi_D per pair
         per_col = 25D0*REAL(n, 8) + REAL(fs%nfs, 8) + REAL((2*fs%mg + 1)*(fs%nev + 1), 8)
         IF (fs_apart) per_col = per_col + REAL(nt*ns + n, 8)
         ncol_max = MAX(1, INT(DIRECT_BLOCK_BYTES/(16D0*per_col)))
      END IF
      DO first = 1, npair, ncol_max
         last = MIN(npair, first + ncol_max - 1)
         nc = last - first + 1
         ALLOCATE (vf(nb, nc, ns), rf(nt, nc, ns), r2(nt, nc, ns), a(n, nc), c(3, n, nc), cl(3, 4, n, nc), &
                   e(geom%nwl, nc), phi(n, nc), p2(n, nc), dws(nc))
         IF (fs_apart) ALLOCATE (rfs(nt, nc, ns), phis(n, nc))
         CALL SYSTEM_CLOCK(c0)
         DO col = 1, nc
            CALL pair_forcing(first + col - 1, col)
         END DO
         IF (multi) THEN
            CALL ASSB_RHS_FIELD_MULTI_IRR(vf, rf, nb, nl, ns, nc)
         ELSE
            CALL ASSB_RHS_FIELD_IRR(vf, rf, nb, nl, ns, nc)
         END IF
         CALL SYSTEM_CLOCK(c1)
         sec(1) = sec(1) + REAL(c1 - c0, 8)/REAL(rate, 8)
         DO col = 1, nc
            CALL pair_densities(first + col - 1, col)
         END DO
         CALL direct_rhs_b2_block(mesh, geom, a, c, e, r2, multi, cl)
         CALL SYSTEM_CLOCK(c0)
         sec(2) = sec(2) + REAL(c0 - c1, 8)/REAL(rate, 8)
         rf = rf + r2
         IF (fs_apart) THEN
            ! The free-surface right sides as their own columns, so that their part of the load is known (parts slice 3)
            rfs = (0D0, 0D0)
            CALL direct_fs_block(fs, store, pairs(:, first:last), dws, kns(first:last), rfs, fsec)
         ELSE IF (PRESENT(fs)) THEN
            CALL direct_fs_block(fs, store, pairs(:, first:last), dws, kns(first:last), rf, fsec)
         END IF
         CALL SYSTEM_CLOCK(c0)
         CALL direct_solve(mesh, rf, phi, multi)
         IF (fs_apart) CALL direct_solve(mesh, rfs, phis, multi)
         IF (PRESENT(phi_field)) THEN
            ! phi_D at the field points from the columns as solved (vf, a, c, e, cl, phi are those of |Omega|). For
            ! Omega < 0 the solved phi is conj(phi_D), since G(-|Omega|) = conj(G(|Omega|)) and the data are conjugated.
            ! Green's formula is linear in (b, phi_D) with the kernel rows of the solver at |Omega|, so the field of the
            ! solved columns is conj(phi_D at the point); the columns with dws < 0 are conjugated, as phi_D below.
            CALL SYSTEM_CLOCK(c1)
            CALL direct_field_at_points(mesh, geom, field_points, phi, vf, a, c, e, cl, multi, phi_field(:, first:last))
            DO col = 1, nc
               IF (dws(col) < 0D0) phi_field(:, first + col - 1) = CONJG(phi_field(:, first + col - 1))
            END DO
            ! The field time is logged apart (sec(4)) and is not in the solve and load time (sec(3)): c0 moves past it.
            CALL SYSTEM_CLOCK(c2)
            sec(4) = sec(4) + REAL(c2 - c1, 8)/REAL(rate, 8)
            c0 = c0 + (c2 - c1)
         END IF
         DO col = 1, nc
            IF (dws(col) < 0D0) phi(:, col) = CONJG(phi(:, col))
            IF (fs_apart) THEN
               IF (dws(col) < 0D0) phis(:, col) = CONJG(phis(:, col))
            END IF
            CALL pair_load(first + col - 1, col)
         END DO
         CALL SYSTEM_CLOCK(c1)
         sec(3) = sec(3) + REAL(c1 - c0, 8)/REAL(rate, 8)
         DEALLOCATE (vf, rf, r2, a, c, cl, e, phi, p2, dws)
         IF (fs_apart) DEALLOCATE (rfs, phis)
      END DO
      IF (PRESENT(seconds)) seconds = sec
      IF (PRESENT(fs)) THEN
         ! The time of the free-surface forcing per Omega, for the cost of the moments on large runs (ledger ruling F12)
         WRITE (*, '(a,F10.5,a,I0,a,F9.2,a,F9.2,a)') ' QTF direct free surface: Omega =', omega, ' rad/s, ', npair, &
            ' pairs, Q_D and outer moments', fsec(1), ' s, contraction', fsec(2), ' s'
      END IF

   CONTAINS

      SUBROUTINE pair_forcing(p, col)
         ! Column col of the block for pair p: the panel field b - B2 (vf) and phi_I2 at the panels (p2), under the kernel
         ! of the pair. The forcing is conjugated for a difference pair with Omega < 0 (header); p2 is not. A sum pair has
         ! dws = omega_1 + omega_2 > 0 and is never conjugated.
         INTEGER, INTENT(IN) :: p, col
         TYPE(incident2) :: inc
         COMPLEX(8) :: gr2(3), dn2, d1(3), d2(3), v1(3), v2(3), q1(3), q2(3), hx(3), b1v, b3v, zero3(3), field
         REAL(8) :: nv(3), x(3), dw, tol, g_wave
         INTEGER :: f, b

         ASSOCIATE (s1 => store(pairs(1, p), pairs(2, p)), s2 => store(pairs(3, p), pairs(4, p)))
            tol = 1D-6*MAX(s1%omega, s2%omega)
            IF (PRESENT(omega_tol)) tol = omega_tol
            IF (kns(p) == KERNEL_SUM) THEN
               dw = s1%omega + s2%omega
               IF (ABS(dw - omega) > tol) THEN
                  WRITE (*, '(a,I0,a,ES14.6,a,ES14.6,a)') ' direct_at_omega: sum pair ', p, ' has Omega =', dw, &
                     ' rad/s, not', omega, ' rad/s.'
                  ERROR STOP 1
               END IF
            ELSE
               dw = s1%omega - s2%omega
               IF (dw == 0D0 .OR. ABS(ABS(dw) - omega) > tol) THEN
                  WRITE (*, '(a,I0,a,ES14.6,a,ES14.6,a)') ' direct_at_omega: pair ', p, ' has Omega =', dw, &
                     ' rad/s, not +-', omega, ' rad/s.'
                  ERROR STOP 1
               END IF
            END IF
            dws(col) = dw
            zero3 = (0D0, 0D0)
            ! phi_I2 with the gravity of the first-order wave, as potential_part
            g_wave = wave_gravity(s1%omega, s1%k, s1%h)
            IF (.NOT. PRESENT(phi2_constant)) THEN
               CALL incident2_pair(s1%omega, s1%k, s1%beta, s2%omega, s2%k, s2%beta, s1%h, g_wave, inc, kns(p))
            END IF
            DO f = 1, n
               b = mesh%body(f)
               nv = mesh%nrm(:, f)
               x = mesh%xc(:, f) - mesh%xr(:, b)
               IF (PRESENT(phi2_constant)) THEN
                  p2(f, col) = phi2_constant
                  dn2 = (0D0, 0D0)
               ELSE
                  p2(f, col) = incident2_phi(inc, mesh%xc(:, f))
                  CALL incident2_grad(inc, mesh%xc(:, f), gr2)
                  dn2 = SUM(nv*gr2)
               END IF
               ! B1 = (Xdot - u).(alpha x n) = n.[(Xdot - u) x alpha], as potential_pair
               d1 = rigid_displacement(s1%xi(:, b), x)
               d2 = rigid_displacement(s2%xi(:, b), x)
               v1 = -CI*s1%omega*d1 - (s1%grad_i(:, f) + s1%grad_p(:, f))
               v2 = -CI*s2%omega*d2 - (s2%grad_i(:, f) + s2%grad_p(:, f))
               b1v = SUM(nv*pair_cross(kns(p), v1, s1%xi(4:6, b), v2, s2%xi(4:6, b)))
               ! B3 = n.d/dt(H x), H x = (1/2) alpha x (alpha x x), d/dt = -i dw (dw = Omega+ for the sum)
               q1 = rigid_displacement([zero3, s1%xi(4:6, b)], x)
               q2 = rigid_displacement([zero3, s2%xi(4:6, b)], x)
               hx = 0.5D0*pair_cross(kns(p), s1%xi(4:6, b), q1, s2%xi(4:6, b), q2)
               b3v = -CI*dw*SUM(nv*hx)
               field = b1v + b3v - dn2
               IF (dw < 0D0) field = CONJG(field)
               vf(mesh%base(f), col, mesh%image(f)) = field
            END DO
         END ASSOCIATE
      END SUBROUTINE pair_forcing

      SUBROUTINE pair_densities(p, col)
         ! Column col of the B2 densities for pair p under its kernel, with the sharp-edge mask of the run (direct_setup),
         ! conjugated for Omega < 0 as the forcing
         INTEGER, INTENT(IN) :: p, col

         ASSOCIATE (s1 => store(pairs(1, p), pairs(2, p)), s2 => store(pairs(3, p), pairs(4, p)))
            CALL b2_stokes_densities(mesh, geom, s1, s2, store_wl(:, pairs(1, p), pairs(2, p)), &
                                     store_wl(:, pairs(3, p), pairs(4, p)), a(:, col), c(:, :, col), e(:, col), &
                                     cl(:, :, :, col), edge_mask, kns(p))
            IF (dws(col) < 0D0) THEN
               a(:, col) = CONJG(a(:, col))
               c(:, :, col) = CONJG(c(:, :, col))
               e(:, col) = CONJG(e(:, col))
               cl(:, :, :, col) = CONJG(cl(:, :, :, col))
            END IF
         END ASSOCIATE
      END SUBROUTINE pair_densities

      SUBROUTINE pair_load(p, col)
         ! T_dir,k = i Omega rho int_S (phi_I2 + phi_D) n_k dS over the panels of the body of mode k, Omega = dws(col)
         ! (signed for the difference, Omega+ for the sum); with parts, the
         ! two integrals are kept apart, and with three slices the free-surface part of phi_D (phis) too
         INTEGER, INTENT(IN) :: p, col
         COMPLEX(8) :: t2(6*mesh%nbody), td(6*mesh%nbody), tf(6*mesh%nbody), fac
         INTEGER :: f, b, j

         t2 = (0D0, 0D0)
         td = (0D0, 0D0)
         DO f = 1, n
            b = mesh%body(f)
            DO j = 1, 6
               t2(6*(b - 1) + j) = t2(6*(b - 1) + j) + p2(f, col)*mesh%nj(j, f)*mesh%area(f)
               td(6*(b - 1) + j) = td(6*(b - 1) + j) + phi(f, col)*mesh%nj(j, f)*mesh%area(f)
            END DO
         END DO
         fac = CI*dws(col)*store(pairs(1, p), pairs(2, p))%rho
         IF (fs_apart) THEN
            tf = (0D0, 0D0)
            DO f = 1, n
               b = mesh%body(f)
               DO j = 1, 6
                  tf(6*(b - 1) + j) = tf(6*(b - 1) + j) + phis(f, col)*mesh%nj(j, f)*mesh%area(f)
               END DO
            END DO
            tdir(:, p) = fac*(t2 + td + tf)
            IF (fs_slice) parts(:, p, 3) = fac*tf     ! parts has the third slice only with write_terms >= 1
         ELSE
            tdir(:, p) = fac*(t2 + td)
         END IF
         IF (PRESENT(parts)) THEN
            parts(:, p, 1) = fac*t2
            parts(:, p, 2) = fac*td
         END IF
      END SUBROUTINE pair_load

      SUBROUTINE check_fs()
         ! fs must be set up at omega for this mesh and lid, with the 6a store of the pairs: fields and partial waves for
         ! every state of store
         LOGICAL :: ok

         ok = ABS(fs%omega - omega) <= 1D-12*omega .AND. fs%n == n .AND. fs%nbase == nb .AND. fs%nsys == ns .AND. &
              fs%nl == nl .AND. (fs%multi .EQV. multi)
         IF (.NOT. ok) THEN
            WRITE (*, '(a)') ' direct_at_omega: the free-surface data are not those of this Omega, mesh and lid'// &
               ' (direct_fs_setup).'
            ERROR STOP 1
         END IF
         ok = ASSOCIATED(fs%f) .AND. ASSOCIATED(fs%inc) .AND. ASSOCIATED(fs%pw)
         IF (ok) ok = SIZE(fs%f, 1) == fs%nfs .AND. SIZE(fs%f, 2) == SIZE(store, 1) .AND. SIZE(fs%f, 3) == SIZE(store, 2) &
                      .AND. SIZE(fs%inc, 1) == fs%nfs .AND. SIZE(fs%inc, 2) == SIZE(store, 1) .AND. &
                      SIZE(fs%inc, 3) == SIZE(store, 2) .AND. SIZE(fs%pw, 1) == SIZE(store, 1) .AND. &
                      SIZE(fs%pw, 2) == SIZE(store, 2)
         IF (.NOT. ok) THEN
            WRITE (*, '(a)') ' direct_at_omega: the free-surface data need the fields at the inner points and the partial'// &
               ' waves of every state of the store.'
            ERROR STOP 1
         END IF
      END SUBROUTINE check_fs

   END SUBROUTINE direct_at_omega


   SUBROUTINE direct_rhs_b2(mesh, geom, s1, s2, wl1, wl2, rhs, is_multi, stored, edge_deg)
      ! The B2 part of the right side for the pair (s1, s2), difference-frequency kernel:
      !   rhs(i, ip) = sum_a RXY(ip, a) r(x of image a of base panel i),  r(x) = -int_S G(x, .) n.[(X.grad)u] dS
      ! in the Stokes form above, with u at the waterline from waterline_velocity, as potential_pair has it.
      !   wl1, wl2   total potential at the waterline segments (nwl), as for potential_pair
      !   rhs        (nbase, nsys), the layout of one column of ASSB_RHS_FIELD (hull rows only)
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

   SUBROUTINE b2_stokes_densities(mesh, geom, s1, s2, wl1, wl2, a, c, e, cl, one_point, kernel)
      ! The pair data of the Stokes form, so that for any potential psi
      !   int_S psi n.[(X.grad)u] dS = sum_p [ psi_p a_p + (int_p grad psi).c_p ] - sum_w psi_w e_w
      ! with psi_p and int_p grad psi panel integrals and psi_w the line integral over segment w:
      !   a(p) = n.(alpha x u)       c(:, p) = (X x u) x n   [(n x g).w = g.(w x n)]       e(w) = (X x u).t
      ! under kernel (optional, KERNEL_DIFF by default, or KERNEL_SUM), the motion of the body that owns the panel or
      ! segment. The waterline velocity (waterline_velocity) is first order and does not depend on the kernel.
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
      INTEGER, INTENT(IN), OPTIONAL :: kernel
      COMPLEX(8), ALLOCATABLE :: uw1(:, :), uw2(:, :), w(:, :)
      COMPLEX(8) :: u1(3), u2(3), d1(3), d2(3), rot(3), xu(3), dw1(6), dw2(6), g1, g2
      REAL(8) :: nv(3), x(3), t(3), t1(3), t2(3)
      INTEGER :: i, b, iw, m, kn
      LOGICAL :: ok

      kn = KERNEL_DIFF
      IF (PRESENT(kernel)) kn = kernel
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
         rot = pair_cross(kn, s1%xi(4:6, b), u1, s2%xi(4:6, b), u2)
         xu = pair_cross(kn, d1, u1, d2, u2)
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
         xu = pair_cross(kn, d1, uw1(:, iw), d2, uw2(:, iw))
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
      !   a (n, nrhs), c (3, n, nrhs), e (nwl, nrhs)  ->  rhs (nbase + nl, nrhs, nsys), the layout of ASSB_RHS_FIELD_IRR,
      !   with nl = 0 (hull rows only) or nl = direct_lid_points (the lid rows too, IRSP = 1)
      ! The panel integrals do not depend on the pair, so one pass over the field points serves all of them. The field
      ! points are the full panels, then the full lid points: image a of base lid point l is row n + (a - 1) nl + l, at
      ! M_a x_l, with the stored lid integrals (stored_lid_integrals) or fresh ones (green_panel_integrals at that point).
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
      REAL(8) :: x(3), rs
      INTEGER :: nb, ns, nrhs, nl, i, ip, im, row
      LOGICAL :: multi, from_store

      nb = mesh%nbase
      ns = mesh%nsys
      nrhs = SIZE(a, 2)
      multi = mesh%nbody > 1
      IF (PRESENT(is_multi)) multi = is_multi
      nl = SIZE(rhs, 1) - nb
      IF (SIZE(a, 1) /= mesh%n .OR. SIZE(c, 1) /= 3 .OR. SIZE(c, 2) /= mesh%n .OR. SIZE(c, 3) /= nrhs .OR. &
          SIZE(e, 1) /= geom%nwl .OR. SIZE(e, 2) /= nrhs .OR. SIZE(rhs, 2) /= nrhs .OR. SIZE(rhs, 3) /= ns .OR. &
          (nl /= 0 .AND. nl /= direct_lid_points(multi))) THEN
         WRITE (*, '(a)') ' direct_rhs_b2_block: a, c, e or rhs does not have the size of the mesh (and the lid).'
         ERROR STOP 1
      END IF
      from_store = .TRUE.
      IF (PRESENT(stored)) from_store = stored

      ALLOCATE (r(nrhs, mesh%n + ns*nl))
      !$OMP PARALLEL DO NUM_THREADS(NTHREAD) SCHEDULE(DYNAMIC, 4) PRIVATE(x, row)
      DO i = 1, mesh%n + ns*nl
         IF (i <= mesh%n) THEN
            x = mesh%xc(:, i)
         ELSE
            row = i - mesh%n
            x = lid_point(row)
         END IF
         IF (i > mesh%n .AND. from_store) THEN
            IF (PRESENT(cl)) THEN
               CALL b2_stokes_at_point(mesh, geom, x, a, c, e, multi, r(:, i), cl, lid_row=row)
            ELSE
               CALL b2_stokes_at_point(mesh, geom, x, a, c, e, multi, r(:, i), lid_row=row)
            END IF
         ELSE IF (from_store .AND. PRESENT(cl)) THEN
            CALL b2_stokes_at_point(mesh, geom, x, a, c, e, multi, r(:, i), cl, i)
         ELSE IF (from_store) THEN
            CALL b2_stokes_at_point(mesh, geom, x, a, c, e, multi, r(:, i), row=i)
         ELSE IF (PRESENT(cl)) THEN
            CALL b2_stokes_at_point(mesh, geom, x, a, c, e, multi, r(:, i), cl, lid_point=i > mesh%n)
         ELSE
            CALL b2_stokes_at_point(mesh, geom, x, a, c, e, multi, r(:, i), lid_point=i > mesh%n)
         END IF
      END DO
      !$OMP END PARALLEL DO

      ! r is -int G n.[(X.grad)u] dS at the full panels and lid points; the block form as ASSB_RBC (ASSB_RBC_IRR)
      DO ip = 1, ns
         DO i = 1, nb + nl
            rhs(i, :, ip) = (0D0, 0D0)
            DO im = 1, ns
               rs = RXY(ip, im)
               IF (multi) rs = RXY_MULTI(ip, im)
               IF (i <= nb) THEN
                  rhs(i, :, ip) = rhs(i, :, ip) - rs*r(:, (im - 1)*nb + i)
               ELSE
                  rhs(i, :, ip) = rhs(i, :, ip) - rs*r(:, mesh%n + (im - 1)*nl + i - nb)
               END IF
            END DO
         END DO
      END DO

   CONTAINS

      FUNCTION lid_point(lrow) RESULT(xl)
         ! Full lid row lrow = (a - 1) nl + l: base lid point l mirrored by image a
         INTEGER, INTENT(IN) :: lrow
         REAL(8) :: xl(3)
         INTEGER :: l, ia

         ia = (lrow - 1)/nl + 1
         l = lrow - (ia - 1)*nl
         IF (multi) THEN
            xl = iXYZ_GLOBAL_MULTI_COMB_P(l, :)
         ELSE
            xl = iXYZ_P(l, :)
         END IF
         IF (ia == 2) THEN
            IF (ISX == 1 .AND. ISY == 0) THEN
               xl(1) = -xl(1)
            ELSE
               xl(2) = -xl(2)
            END IF
         END IF
      END FUNCTION lid_point

   END SUBROUTINE direct_rhs_b2_block

   SUBROUTINE b2_stokes_at_point(mesh, geom, x, a, c, e, is_multi, out, cl, row, lid_row, field_mode, lid_point)
      ! int_S G(x, .) n.[(X.grad)u] dS for each column of the densities, at the field point x (full coordinates):
      !   out = sum_p [ int_p G a_p + (int_p grad_xi G).c_p ] - sum_w (int_w G dl) e_w
      ! With cl (3, 4, n, ncol), on the near panels also sum_ab cl(a, b, p) int_p d_a(1/r + 1/r1) (xi - xi_p)_b dS and
      ! cl(:, 4, p).int_p (1/r + 1/r1)(xi - xi_p) dS.
      ! With row (a full panel index, whose centroid x must be), int_p G and int_p grad G come from the stored arrays
      ! (stored_panel_integrals); with lid_row (a full lid row of direct_rhs_b2_block, whose point x must be) from the
      ! stored lid integrals (stored_lid_integrals); otherwise from green_panel_integrals.
      ! field_mode, optional (default .false.): the fresh integrals (green_panel_integrals, green_line_integral) pass FIELD
      ! to FINGREEN3D, the route of the field rows (field_mode of field_influence_row), for field points off the hull
      ! (QTFDirectField). Absent, every caller of Phase 7a keeps the default route bit for bit.
      ! lid_point, optional (default .false.): x is a lid point; passed on to green_panel_integrals (fresh integrals only).
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      REAL(8), INTENT(IN) :: x(3)
      COMPLEX(8), INTENT(IN) :: a(:, :), c(:, :, :), e(:, :)
      LOGICAL, INTENT(IN) :: is_multi
      COMPLEX(8), INTENT(OUT) :: out(:)
      COMPLEX(8), INTENT(IN), OPTIONAL :: cl(:, :, :, :)
      INTEGER, INTENT(IN), OPTIONAL :: row, lid_row
      LOGICAL, INTENT(IN), OPTIONAL :: field_mode, lid_point
      COMPLEX(8) :: g_int, grad_int(3), g_line, acc
      REAL(8) :: xv(3, 4), m0, g0(3), m1(3, 3), m1i(3, 3), mg(3), mgi(3), xim(3)
      INTEGER :: p, iw, k, nv, m, ia, ib
      ! The per-panel loop below has no array temporaries: under /heap-arrays0 an array temporary (an array
      ! constructor, a vector subscript, SUM of an array expression) is a heap allocation, and the threads of
      ! direct_rhs_b2_block then wait on each other in the allocator. With this loop, rankine_moments and
      ! stored_panel_integrals written so, the B2 right side on 480 panels went from 0.9 s at both 1 and 4 threads to
      ! 0.68 s and 0.34 s (Phase 7a Task 3). Outside it, xim (once per point) is an array constructor.

      xim = [x(1), x(2), -x(3)]
      out = (0D0, 0D0)
      DO p = 1, mesh%n
         IF (PRESENT(row)) THEN
            CALL stored_panel_integrals(mesh, row, p, is_multi, g_int, grad_int)
         ELSE IF (PRESENT(lid_row)) THEN
            CALL stored_lid_integrals(mesh, lid_row, p, is_multi, g_int, grad_int)
         ELSE
            CALL green_panel_integrals(mesh, x, p, is_multi, g_int, grad_int, field_mode, lid_point)
         END IF
         DO k = 1, SIZE(out)
            out(k) = out(k) + g_int*a(p, k) + (grad_int(1)*c(1, p, k) + grad_int(2)*c(2, p, k) + grad_int(3)*c(3, p, k))
         END DO
         IF (.NOT. PRESENT(cl)) CYCLE
         IF ((x(1) - mesh%xc(1, p))**2 + (x(2) - mesh%xc(2, p))**2 + (x(3) - mesh%xc(3, p))**2 > &
             (NEAR_PANEL_SIZES*mesh%psize(p))**2) CYCLE
         nv = mesh%ncn(p)
         DO m = 1, nv
            xv(:, m) = mesh%xnode(:, mesh%ncon(m, p))
         END DO
         CALL rankine_moments(xv, nv, mesh%nrm(:, p), x, mesh%xc(:, p), m0, g0, m1, mg)
         CALL rankine_moments(xv, nv, mesh%nrm(:, p), xim, mesh%xc(:, p), m0, g0, m1i, mgi)
         m1 = m1 + m1i
         mg = mg + mgi
         DO k = 1, SIZE(out)
            acc = (0D0, 0D0)
            DO ib = 1, 3
               DO ia = 1, 3
                  acc = acc + m1(ia, ib)*cl(ia, ib, p, k)
               END DO
               acc = acc + mg(ib)*cl(ib, 4, p, k)
            END DO
            out(k) = out(k) + acc
         END DO
      END DO
      DO iw = 1, geom%nwl
         CALL green_line_integral(x, geom%wl_x1(:, iw), geom%wl_x2(:, iw), mesh%psize(geom%wl_panel(iw)), g_line, &
                                  field_mode)
         DO k = 1, SIZE(out)
            out(k) = out(k) - g_line*e(iw, k)
         END DO
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
      REAL(8) :: rk(4), f(3), area, dist2, psz
      INTEGER :: i, j, is, k
      ! Written out, with no array expressions (see b2_stokes_at_point): this runs in its per-panel loop, n^2 times per
      ! Omega

      i = mesh%base(row)
      j = mesh%base(p)
      is = MERGE(1, 2, mesh%image(row) == mesh%image(p))
      ! image_signs(mesh%image(p), ISX, ISY)
      f = 1D0
      IF (mesh%image(p) == 2) THEN
         IF (ISX == 1 .AND. ISY == 0) THEN
            f(1) = -1D0
         ELSE
            f(2) = -1D0
         END IF
      END IF
      dist2 = 0D0
      IF (is_multi) THEN
         DO k = 1, 4
            grn(k) = CGRN_MULTI_COMB(i, j, is, k)
            rk(k) = RKBN_MULTI_COMB(i, j, is, k)
         END DO
         area = DS_MULTI_COMB(j)
         psz = PNSZ_MULTI_COMB(j)
         DO k = 1, 3
            dist2 = dist2 + (XYZ_GLOBAL_MULTI_COMB_P(i, k) - XYZ_GLOBAL_MULTI_COMB_P(j, k))**2
         END DO
      ELSE
         DO k = 1, 4
            grn(k) = CGRN(i, j, is, k)
            rk(k) = RKBN(i, j, is, k)
         END DO
         area = DS(j)
         psz = PNSZ(j)
         DO k = 1, 3
            dist2 = dist2 + (XYZ_P(i, k) - XYZ_P(j, k))**2
         END DO
      END IF
      IF (SQRT(dist2) <= NEAR_PANEL_SIZES*psz) THEN
         g_int = rk(1) + grn(1)*area
         DO k = 1, 3
            grad_int(k) = f(k)*(rk(k + 1) + grn(k + 1)*area)
         END DO
      ELSE
         g_int = grn(1)*area
         DO k = 1, 3
            grad_int(k) = f(k)*grn(k + 1)*area
         END DO
      END IF
   END SUBROUTINE stored_panel_integrals

   SUBROUTINE stored_lid_integrals(mesh, lid_row, p, is_multi, g_int, grad_int)
      ! As stored_panel_integrals at the full lid point lid_row = (a - 1) nl + l (image a of base lid point l), from the
      ! stored lid integrals of CALGREEN_IRR (DGRN, PKBN; the MULTI_COMB forms for several bodies): stored image a xor b
      ! for full panel p (image b of base panel j), gradient mirrored by M_b, and the flag rule of ASSB_RBC_IRR (base lid
      ! point to base panel centroid within 50 panel sizes), so that the B2 part uses the discrete G of the lid rows.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      INTEGER, INTENT(IN) :: lid_row, p
      LOGICAL, INTENT(IN) :: is_multi
      COMPLEX(8), INTENT(OUT) :: g_int, grad_int(3)
      COMPLEX(8) :: grn(4)
      REAL(8) :: rk(4), f(3), area, dist2, psz
      INTEGER :: l, ia, nl, j, is, k

      IF (is_multi) THEN
         nl = SIZE(DGRN_MULTI_COMB, 1)
      ELSE
         nl = SIZE(DGRN, 1)
      END IF
      ia = (lid_row - 1)/nl + 1
      l = lid_row - (ia - 1)*nl
      j = mesh%base(p)
      is = MERGE(1, 2, ia == mesh%image(p))
      f = 1D0
      IF (mesh%image(p) == 2) THEN
         IF (ISX == 1 .AND. ISY == 0) THEN
            f(1) = -1D0
         ELSE
            f(2) = -1D0
         END IF
      END IF
      dist2 = 0D0
      IF (is_multi) THEN
         DO k = 1, 4
            grn(k) = DGRN_MULTI_COMB(l, j, is, k)
            rk(k) = PKBN_MULTI_COMB(l, j, is, k)
         END DO
         area = DS_MULTI_COMB(j)
         psz = PNSZ_MULTI_COMB(j)
         DO k = 1, 3
            dist2 = dist2 + (iXYZ_GLOBAL_MULTI_COMB_P(l, k) - XYZ_GLOBAL_MULTI_COMB_P(j, k))**2
         END DO
      ELSE
         DO k = 1, 4
            grn(k) = DGRN(l, j, is, k)
            rk(k) = PKBN(l, j, is, k)
         END DO
         area = DS(j)
         psz = PNSZ(j)
         DO k = 1, 3
            dist2 = dist2 + (iXYZ_P(l, k) - XYZ_P(j, k))**2
         END DO
      END IF
      IF (SQRT(dist2) <= NEAR_PANEL_SIZES*psz) THEN
         g_int = rk(1) + grn(1)*area
         DO k = 1, 3
            grad_int(k) = f(k)*(rk(k + 1) + grn(k + 1)*area)
         END DO
      ELSE
         g_int = grn(1)*area
         DO k = 1, 3
            grad_int(k) = f(k)*grn(k + 1)*area
         END DO
      END IF
   END SUBROUTINE stored_lid_integrals

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

   SUBROUTINE green_panel_integrals(mesh, x, p, is_multi, g_int, grad_int, field_mode, lid_point)
      ! int_p G(x, xi) dS_xi and int_p grad_xi G(x, xi) dS_xi over full panel p (image b of base panel j), x in full
      ! coordinates. Near panels (within NEAR_PANEL_SIZES panel sizes, as CALGREEN) take 1/r and 1/r1 as panel integrals,
      ! far ones the centroid value. On the own panel the 1/r part of the gradient is the principal value (in-plane part).
      ! field_mode, optional (default .false.): FIELD of FINGREEN3D (see b2_stokes_at_point); no effect in deep water.
      ! lid_point, optional (default .false.): x is a lid point of irregular frequency removal. As CALGREEN_IRR
      ! (CALGREEN_MULTI_IRR) does for the stored lid integrals, the wave part is then integrated over the panel
      ! (WAVE_PART_PANEL) when x is within LID_WAVE_NEAR panel sizes of its centroid, at finite frequency.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: x(3)
      INTEGER, INTENT(IN) :: p
      LOGICAL, INTENT(IN) :: is_multi
      COMPLEX(8), INTENT(OUT) :: g_int, grad_int(3)
      LOGICAL, INTENT(IN), OPTIONAL :: field_mode, lid_point
      REAL(8) :: f(3), xf(3), xq(3), d(3), r1, area, rk(4), rk_image(4), pnsz_j
      COMPLEX(8) :: grn(4), gint(4)
      INTEGER :: j, flag
      LOGICAL :: lid

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
         ! An absent field_mode passes on as an absent FIELD: the default route
         CALL FINGREEN3D(xq(1), xf(1), xq(2), xf(2), xq(3), xf(3), V, WVN, NK, H, grn, flag, FIELD=field_mode)
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
      lid = .FALSE.
      IF (PRESENT(lid_point)) lid = lid_point
      IF (lid .AND. ABS(V) > 1D-8 .AND. ABS(V + 1D0) > 1D-8) THEN
         ! Not MERGE: it evaluates both arguments, and PNSZ_MULTI_COMB is not allocated for a single body
         IF (is_multi) THEN
            pnsz_j = PNSZ_MULTI_COMB(j)
         ELSE
            pnsz_j = PNSZ(j)
         END IF
         IF (NORM2(xf - xq) <= LID_WAVE_NEAR*pnsz_j) THEN
            IF (is_multi) THEN
               CALL WAVE_PART_PANEL_MULTI(xf, j, V, WVN, NK, H, gint)
            ELSE
               CALL WAVE_PART_PANEL(xf, j, V, WVN, NK, H, gint)
            END IF
            g_int = rk(1) + rk_image(1) + gint(1)
            grad_int = f*(rk(2:4) + rk_image(2:4) + gint(2:4))
            RETURN
         END IF
      END IF
      g_int = rk(1) + rk_image(1) + (grn(1) - 1D0/r1)*area
      grad_int = f*(rk(2:4) + rk_image(2:4) + (grn(2:4) + d/r1**3)*area)
   END SUBROUTINE green_panel_integrals

   SUBROUTINE green_line_integral(x, x1, x2, psize, g_line, field_mode)
      ! int G(x, xi) dl over the straight segment from x1 to x2 on z = 0 (the z of x1 and x2 is set to 0, where the
      ! waterline midpoints are). 1/r and 1/r1 in closed form, int dl/r = ln((R1 + R2 + L)/(R1 + R2 - L)), the rest at the
      ! midpoint. psize decides near and far as for the panels. x must not lie on the segment.
      ! field_mode, optional (default .false.): FIELD of FINGREEN3D, as green_panel_integrals.
      REAL(8), INTENT(IN) :: x(3), x1(3), x2(3), psize
      COMPLEX(8), INTENT(OUT) :: g_line
      LOGICAL, INTENT(IN), OPTIONAL :: field_mode
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
         CALL FINGREEN3D(xm(1), x(1), xm(2), x(2), xm(3), x(3), V, WVN, NK, H, grn, flag, FIELD=field_mode)
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
                 r1, r2, i0, t1, t2, g2(2), m2(2, 2), tiny, mr(2), dd, dv(3), dp(2)
      INTEGER :: k
      ! Written without array expressions as arguments or constructors: under /heap-arrays0 each is a heap allocation, and
      ! the threads of direct_rhs_b2_block then wait on each other in the allocator (this routine runs twice per near pair)

      o = 0D0
      DO k = 1, nv
         o = o + xv(:, k)
      END DO
      o = o/nv
      e1 = xv(:, 1) - o
      e1 = e1 - (e1(1)*nrm(1) + e1(2)*nrm(2) + e1(3)*nrm(3))*nrm
      e1 = e1/SQRT(e1(1)**2 + e1(2)**2 + e1(3)**2)
      e2(1) = nrm(2)*e1(3) - nrm(3)*e1(2)
      e2(2) = nrm(3)*e1(1) - nrm(1)*e1(3)
      e2(3) = nrm(1)*e1(2) - nrm(2)*e1(1)
      cen = 0D0
      DO k = 1, nv
         dv = xv(:, k) - o
         q(1, k) = dv(1)*e1(1) + dv(2)*e1(2) + dv(3)*e1(3)
         q(2, k) = dv(1)*e2(1) + dv(2)*e2(2) + dv(3)*e2(3)
         cen = cen + q(:, k)
      END DO
      cen = cen/nv
      dv = x - o
      xp(1) = dv(1)*e1(1) + dv(2)*e1(2) + dv(3)*e1(3)
      xp(2) = dv(1)*e2(1) + dv(2)*e2(2) + dv(3)*e2(3)
      z = dv(1)*nrm(1) + dv(2)*nrm(2) + dv(3)*nrm(3)
      dv = xc - o
      qc(1) = dv(1)*e1(1) + dv(2)*e1(2) + dv(3)*e1(3)
      qc(2) = dv(1)*e2(1) + dv(2)*e2(2) + dv(3)*e2(3)
      az = ABS(z)

      m0 = 0D0
      g2 = 0D0
      m2 = 0D0
      mr = 0D0
      DO k = 1, nv
         p1 = q(:, k)
         p2 = q(:, MOD(k, nv) + 1)
         dp = p2 - p1
         el = SQRT(dp(1)**2 + dp(2)**2)
         IF (el <= 0D0) CYCLE
         tiny = 1D-12*el
         tau = dp/el
         nu(1) = tau(2)
         nu(2) = -tau(1)
         dp = 0.5D0*(p1 + p2) - cen
         IF (nu(1)*dp(1) + nu(2)*dp(2) < 0D0) nu = -nu
         dp = p1 - xp
         s0 = -(dp(1)*tau(1) + dp(2)*tau(2))
         h = dp(1)*nu(1) + dp(2)*nu(2)
         r1 = SQRT(dp(1)**2 + dp(2)**2 + z*z)
         dp = p2 - xp
         r2 = SQRT(dp(1)**2 + dp(2)**2 + z*z)
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
