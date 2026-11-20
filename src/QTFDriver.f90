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
MODULE QTFDriver
   ! Quadratic part of the difference- and sum-frequency QTFs (qtf in SecondOrder.in), from the near-field pressure
   ! integration of NearFieldDrift with two states at different frequencies.
   !
   ! The QTFs need the first-order solution of every frequency at once, so the converted state and the waterline potential
   ! of each frequency and heading are kept in memory (qtf_keep, called by the driver after the hull velocity is filled).
   ! Only what the near-field terms read is kept: phi_i, phi_p, grad_i, grad_p, xi. After the frequency loop, qtf_finish
   ! evaluates every pair and writes QTF_Quadratic.10d and .10s.
   !
   ! With components a_i of amplitude and the time factor e^{-i omega t}, the second-order load is
   ! (docs/output_files.md, there conjugated)
   !    F(t) = sum_i sum_j Re[ T-_ij a_i a_j* e^{-i(w_i - w_j)t} + T+_ij a_i a_j e^{-i(w_i + w_j)t} ],
   ! over all i and j. T-_ij is near_field_pair(state i, state j, KERNEL_DIFF) and T+_ij the same with KERNEL_SUM. The files
   ! hold their conjugates, scaled as MeanDriftPressure.9 (SecondOrderOutput), on the full square of frequency pairs.
   ! Heading pairs: beta1 = beta2 for mean_drift = 1, every pair for mean_drift = 2; state i has heading beta1 and state
   ! j heading beta2.
   !
   ! With qtf_terms = body the potential part of the difference QTF comes from the indirect method (Molin 1979), which needs
   ! the radiation potentials psi_k at every difference frequency Omega = omega_i - omega_j > 0 of the grid within
   ! qtf_max_dw. qtf_radiation_list gives these Omega, with no interpolation: a non-uniform grid can have N(N-1)/2 of them.
   ! Two Omega are one only to rounding (qtf_psi_tolerance, 1e-9 of the largest frequency), and the lookup of an Omega
   ! (qtf_psi_index) and the check of potential_pair use that same tolerance. qtf_init checks before the solves that every
   ! pair finds its Omega, and refuses a grid with two frequencies equal to within that tolerance but not exactly.
   ! An Omega on the main grid takes psi from the radiation solve of that frequency; the others are solved after the
   ! frequency loop (qtf_radiation_omegas). psi_store keeps psi of every Omega, filled by qtf_keep_psi.
   !
   ! qtf_finish then writes QTF_Indirect.11d = QTF_Quadratic.10d + T_pot, with T_pot the body-forcing terms of potential_pair
   ! (QTFPotential), and with qtf_terms = full also the free-surface integral (below). Same rows, order and
   ! scaling as .10d, qtf_max_dw included. On the exact diagonal (equal frequencies, also with two headings) T_pot = 0, so
   ! there .11d is .10d. T_pot is computed Omega by Omega (potential_part), so that grad psi and psi at the waterline are
   ! built once per difference frequency, and kept per pair in tpot until the files are written. With write_terms >= 1 the
   ! five terms of each pair (seven with qtf_terms = full) go to QTF_PotentialTerms.txt as they are computed: its rows are
   ! in ascending |Omega|, not in the order of .11d, so a reader matches them by their keys PER1 PER2 BETA1 BETA2 I (header
   ! in the file).
   !
   ! With qtf_b2 = control the surface piece of B2 is moved to the control surface C and the free surface F of each body
   ! (QTFPotential). The driver then hands over the points of C and F (qtf_set_control_surface), the total velocity there of
   ! every kept frequency and heading (qtf_keep_cs_velocity, indexed like store) and grad psi_k there at every difference
   ! frequency (qtf_keep_cs_psi, indexed like psi_store). With qtf_edge_weights the panel weights of the hull term at sharp
   ! edges come with qtf_set_edge_weights.
   !
   ! With qtf_terms = full the store also holds what the free-surface integral of the second-order potential needs
   ! (QTFFreeSurface, PartialWaves), allocated by qtf_init: per kept frequency and heading the total first-order fields at
   ! the n_fs inner free-surface points (fs_point_field: phi, grad phi, phi_zz; qtf_keep_fs_state, indexed like store) and
   ! the partial waves of the scattered potential outside the circle r = R; per difference frequency psi_k at the inner
   ! points (value only) and the partial waves of each psi_k (qtf_keep_fs_psi, indexed like psi_store). Each pw_set holds
   ! (2M + 1)(N_ev + 1) complex coefficients, so the free-surface share is 16 (5 n_fs + n_pw) bytes per frequency and
   ! heading plus 16 x 6 nbody (n_fs + n_pw) bytes per difference frequency (qtf_stored_mb with nfs and npw).
   ! potential_part then adds the free-surface terms FS_inner and FS_outer of each pair (fs_pair_terms) to T_pot, Omega by
   ! Omega: the inner panel sum, the outer partial-wave integral and the sliver between the polygon of a generated inner
   ! mesh and the circle. Omega = 0 (the diagonal) has none. A pair with omega_i < omega_j is evaluated as (j, i) and
   ! conjugated, so that the outer integral always has omega_1 > omega_2. The pairs run one after the other; the outer
   ! integral is parallel over its radial intervals (one OpenMP level, the same bits for any number of threads).
   !
   ! With qtf_direct = 1 (Phase 7a) the potential part is also found by the direct method: phi_D is solved on the body at
   ! each difference frequency (QTFDirect), in the loop of radiation-only solves after the frequency loop, where the store
   ! is complete and the solver is at Omega. qtf_direct_at_omega puts T_dir of every pair of that Omega into tdir, which
   ! has the slots of tpot. With qtf_terms = full (Phase 7b) the caller passes the free-surface data of direct_fs_setup,
   ! which qtf_direct_at_omega points to the free-surface store, so that T_dir has the free-surface forcing too.
   !
   ! The sum frequency (qtf = 2 or 3 with qtf_terms = body or full; Phase 7c-2 with qtf_direct = 1, Phase 6b also with
   ! qtf_direct = 0, where only the indirect T+_pot is computed): the list psi_need is
   ! the union of the difference and the sum frequencies Omega = omega_i + omega_j (qtf_omega_list), and psi, grad psi on
   ! C and F and psi on the inner free surface are kept at every one of them, indexed like psi_store. qtf_direct_at_omega
   ! solves every difference pair and every sum pair of an Omega with one factorisation, the sum pairs in a
   ! direct_at_omega call of their own, so that the difference results keep the bits of a run without the sum.
   ! T+ is symmetric, T+(i, b1, j, b2) = T+(j, b2, i, b1) (pair_plus; ruling R6 of the plan), so only the canonical half of
   ! the sum pairs is solved (qtf_direct_sum_pairs: i < j, or i = j with b1 <= b2) and each result fills both slots of
   ! tdir_sum. qtf_finish then writes QTF_Direct.12s: the rows of .10s whose sum is kept (sum_pair_index > 0), in its
   ! order, each .10s + T+_dir, the diagonal included (Omega = 2 omega_i is a regular pair of the sum). It also computes
   ! the indirect sum load T+_pot of the same pairs from the same psi (potential_part_sum, into tpot_sum; ruling R1), and
   ! with write_terms >= 1 writes it beside the parts of T+_dir in QTF_DirectTermsSum.txt, the in-run check of .12s.
   ! Phase 6b: the same T+_pot is computed with qtf_direct = 0, where no direct solve runs (the radiation-only solves
   ! give the same psi at the sum frequencies), and qtf_finish writes QTF_Indirect.11s, the rows of .10s whose sum is kept
   ! (sum_pair_index > 0), each .10s + T+_pot, the diagonal included, and with write_terms >= 1 the five body terms (and
   ! with qtf_terms = full FS_inner and FS_outer) of every row in QTF_PotentialTermsSum.txt. The two routes give one .11s.
   ! The difference stores (tpot, tdir and their parts) are allocated only with the difference QTF (qtf = 1 or 3).
   !
   ! Phase 6c, qtf_motion = 1: the second-order motion xi2 of the free bodies at every row of .11d and .11s (and with
   ! qtf_direct = 1 of .12d and .12s), from the force of that row: per Omega of the list, the added mass and damping from
   ! psi (radiation_coefficients of QTFMotion), the matrix of the equation of motion with the body matrices that qtf_init
   ! is given (mass, hydrostatic restoring, external stiffness and linear damper: B(Omega) + BLNR only, ruling R5) and one
   ! solve over every pair of that Omega, per kernel and route (motion_part, into xi2 and xi2_sum). On the difference
   ! diagonal Omega = 0 the equation is static, (C + K_ext) xi2 = F (R6). An Omega whose matrix has rcond < 1e-12 is logged
   ! with the modes without restoring, and its rows are left out of the motion files (R6). qtf_finish writes
   ! QTF_Motion.11d, .11s (and .12d, .12s): the rows of the force file in its order, scaled by A^2/L (translations) and
   ! A^2/L^2 (rotations) and conjugated as .10d (R4), after a header that names the route and what is left out. A fixed
   ! body (motion_state = fixed) gets none: qtf_init logs it and allocates nothing for the motion.
   ! With qtf_field = 1 and qtf_direct = 1 the radiated wave of that motion also goes into QTF_Field.txt (R7): at each
   ! Omega qtf_direct_at_omega takes psi_j at the field points (psi_at_points, into mot_psi_pts), and qtf_finish, after the
   ! motion, adds phi_R = sum_j (-i Omega xi_j) psi_j of the direct route to the direct part of every row whose Omega has a
   ! regular matrix (motion_field_part; QTFField has the forms and the conjugation of Omega < 0).
   !
   ! A psi with NaN or Inf (a radiation solve that failed) is kept, but counted and logged by qtf_keep_psi, and every pair
   ! whose terms are not finite is counted by potential_part. The files are still written (.10d does not use psi); the main
   ! program ends with exit code 1 when qtf_nonfinite_psi() or qtf_nonfinite_pairs() is not zero, as for first-order NaN.

   USE SecondOrderTypes, ONLY: so_options, panel_mesh, first_order_state, surface_geom, cs_points, QTF_TERMS_BODY, &
                               QTF_TERMS_FULL
   USE QTFFreeSurface, ONLY: fs_point_field, fs_pair_terms, fs_wl_state, fs_waterline_state, fs_edge_points, fs_edge_state, &
                             fs_edge_psi
   USE FieldVelocity, ONLY: incident_at
   USE PartialWaves, ONLY: pw_set
   USE NearFieldDrift, ONLY: near_field_pair, KERNEL_DIFF, KERNEL_SUM
   USE SecondOrderOutput, ONLY: write_qtf_row, write_qtf_terms_row, write_qtf_terms_header, write_direct_terms_header, &
                                write_direct_sum_terms_header, write_potential_sum_terms_header, write_motion_row, &
                                write_motion_header
   USE QTFMotion, ONLY: radiation_coefficients, motion_matrix, motion_solve_omega
   USE QTFPotential, ONLY: potential_pair, psi_surface_gradient, psi_waterline, NPOT_TERMS
   USE QTFDirect, ONLY: direct_at_omega, direct_fs_data, direct_field_at_points
   USE QTFField, ONLY: qfield_direct_active, qfield_direct_points, qfield_keep_direct, qfield_radiated_active, &
                       qfield_add_radiated, radiated_potential
   USE WaveDyn_mod, ONLY: W1
   USE FirstOrderState, ONLY: wave_gravity
   USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qtf_init, qtf_keep, qtf_pair_load, qtf_finish, qtf_stored_mb, qtf_kept
   PUBLIC :: qtf_omega_list, irregular_frequency_estimate, sum_pair_index, qtf_sum_limits
   PUBLIC :: qtf_radiation_list, qtf_radiation_omegas, qtf_direct_omegas, qtf_omega_on_grid, qtf_keep_psi, qtf_psi_index, qtf_psi_count, qtf_psi_tolerance
   PUBLIC :: qtf_nonfinite_psi, qtf_nonfinite_pairs
   PUBLIC :: qtf_direct_at_omega, qtf_direct_pairs, qtf_direct_missing, qtf_diff_omega, qtf_direct_sum_pairs, qtf_sum_omega
   PUBLIC :: qtf_direct_sum_load, qtf_direct_diff_load, qtf_indirect_sum_load, potential_part_sum, psi_bad_files
   PUBLIC :: qtf_set_control_surface, qtf_keep_cs_velocity, qtf_keep_cs_psi, qtf_set_edge_weights
   PUBLIC :: qtf_keep_fs_state, qtf_keep_fs_psi, qtf_fs_state, qtf_fs_psi, qtf_fs_active, qtf_set_fs_inner
   PUBLIC :: potential_part, qtf_motion_part, qtf_motion_load, psi_at_points

   LOGICAL :: active = .FALSE.
   INTEGER :: kind_qtf = 0, nbeta = 0, nkept = 0, all_headings_pairs = 0
   REAL(8) :: max_dw = 0D0
   CHARACTER(LEN=:), ALLOCATABLE :: outdir
   TYPE(first_order_state), ALLOCATABLE :: store(:, :)   ! (nfreq, nbeta)
   COMPLEX(8), ALLOCATABLE :: store_wl(:, :, :)          ! (nwl, nfreq, nbeta)
   REAL(8), ALLOCATABLE :: store_oufr(:), beta_deg(:)    ! (nfreq), (nbeta)
   ! Radiation potentials at the difference frequencies, for qtf_terms = body
   REAL(8), ALLOCATABLE :: psi_need(:)                   ! (nneed), the Omega, ascending
   INTEGER, ALLOCATABLE :: psi_reuse(:)                  ! (nneed), main-grid index of each Omega, or 0
   REAL(8) :: psi_tol = 0D0                              ! two Omega closer than this are the same (qtf_psi_tolerance)
   ! Relative tolerance of the difference frequencies: rounding only. Differences of distinct grid values that are closer
   ! than this are not merged into one psi (a period grid has many that differ by 1e-6 and more)
   REAL(8), PARAMETER :: PSI_REL_TOL = 1D-9
   ! Phase 7c-2: the list is the union of the difference frequencies (is_diff) and, with the sum on, of the sum frequencies
   ! omega_i + omega_j (is_sum). sum_cut_irr is the estimate of the lowest irregular frequency (HUGE: none applies): sums at
   ! or above it are left out. max_sw is qtf_max_sw (0: all).
   LOGICAL, ALLOCATABLE :: psi_is_diff(:), psi_is_sum(:) ! (nneed)
   LOGICAL :: diff_on = .FALSE., sum_on = .FALSE.        ! the difference and the sum kernel of the direct solve
   REAL(8) :: sum_cut_irr = HUGE(1D0), max_sw = 0D0
   COMPLEX(8), ALLOCATABLE :: psi_store(:, :, :)         ! (n_panels, 6*nbody, nneed)
   LOGICAL, ALLOCATABLE :: psi_kept(:)                   ! (nneed)
   ! Potential part of the difference QTF (qtf_terms = body)
   LOGICAL :: body_active = .FALSE.                      ! QTF_Indirect.11d is written
   LOGICAL :: pot_terms = .FALSE.                        ! QTF_PotentialTerms.txt is written (write_terms >= 1)
   LOGICAL :: pot_sum_terms = .FALSE.                    ! QTF_PotentialTermsSum.txt is written (write_terms >= 1, sum on)
   COMPLEX(8), ALLOCATABLE :: tpot(:, :, :, :)           ! (6*nbody, nfreq, nfreq, heading pairs), e^{-i omega t}
   ! Direct method (qtf_direct = 1, Phase 7a): T_dir of each pair, in the slots of tpot, filled Omega by Omega by
   ! qtf_direct_at_omega during the radiation-only loop after the frequency loop, when the store is complete
   LOGICAL :: direct_active = .FALSE.
   COMPLEX(8), ALLOCATABLE :: tdir(:, :, :, :)           ! (6*nbody, nfreq, nfreq, heading pairs), e^{-i omega t}
   LOGICAL, ALLOCATABLE :: tdir_kept(:)                  ! (nneed), the direct solve has run at psi_need(m)
   ! With write_terms >= 1 (direct_terms), the two parts of tdir, for QTF_DirectTerms.txt: from phi_I2 and from phi_D
   COMPLEX(8), ALLOCATABLE :: tdir_i2(:, :, :, :), tdir_d(:, :, :, :)   ! as tdir
   ! With qtf_terms = full, the part from phi_D splits into that of the body forcing (tdir_d) and of the free-surface
   ! forcing (tdir_fs), the third column pair of QTF_DirectTerms.txt
   COMPLEX(8), ALLOCATABLE :: tdir_fs(:, :, :, :)                       ! as tdir
   LOGICAL :: direct_terms = .FALSE.                     ! QTF_DirectTerms.txt is written
   ! The sum frequency (Phase 7c-2): T+_dir and T+_pot of each pair, both slots of a symmetric pair filled (R6), and with
   ! write_terms >= 1 the parts of T+_dir, as those of tdir, for QTF_DirectTermsSum.txt
   COMPLEX(8), ALLOCATABLE :: tdir_sum(:, :, :, :), tpot_sum(:, :, :, :)          ! (6*nbody, nfreq, nfreq, heading pairs)
   COMPLEX(8), ALLOCATABLE :: tdir_sum_i2(:, :, :, :), tdir_sum_d(:, :, :, :), tdir_sum_fs(:, :, :, :)   ! as tdir_sum
   LOGICAL :: sum_terms = .FALSE.                        ! QTF_DirectTermsSum.txt is written (write_terms >= 1)
   ! qtf_b2 = control: points of C and F of every body, and the fields there
   LOGICAL :: b2_control = .FALSE.
   TYPE(cs_points) :: cs
   COMPLEX(8), ALLOCATABLE :: cs_u(:, :, :, :)           ! (3, cs%n, nfreq, nbeta), total velocity, as store
   LOGICAL, ALLOCATABLE :: cs_u_kept(:, :)               ! (nfreq, nbeta)
   COMPLEX(8), ALLOCATABLE :: cs_gpsi(:, :, :, :)        ! (3, cs%n, 6*nbody, nneed), grad psi_k, as psi_store
   LOGICAL, ALLOCATABLE :: cs_psi_kept(:)                ! (nneed)
   REAL(8), ALLOCATABLE :: edge_w(:)                     ! (n_panels), weights of the B2 hull term at sharp edges, if set
   ! qtf_terms = full: the inner free surface, the circle and the fields there
   LOGICAL :: fs_active = .FALSE.
   REAL(8), PARAMETER :: FS_GEOM_TOL = 1D-12                ! relative, for the circle and depth of a kept pw_set
   TYPE(cs_points) :: fs_store_pts                       ! the inner free-surface points (build_fs_inner)
   REAL(8) :: fs_store_centre(2) = 0D0, fs_store_radius = 0D0   ! the circle r = R of the partial waves
   REAL(8) :: fs_store_depth = 0D0                       ! the water depth h of the partial waves
   INTEGER :: fs_store_mmax = 0, fs_store_nev = 0        ! M and N_ev of every pw_set of the store
   INTEGER :: fs_store_np = 0                            ! chords of the polygon of a generated inner mesh, 0 for a file mesh
   ! fs_f, fs_pw and fs_inc_store are TARGETs: the free-surface data of the direct solve point to them (qtf_direct_at_omega)
   TYPE(fs_point_field), ALLOCATABLE, TARGET :: fs_f(:, :, :)    ! (n_fs, nfreq, nbeta), total first-order fields, as store
   TYPE(pw_set), ALLOCATABLE, TARGET :: fs_pw(:, :)              ! (nfreq, nbeta), partial waves of the scattered potential
   ! (n_fs, nfreq, nbeta), the incident fields at the inner points, for the direct solve with qtf_terms = full
   TYPE(fs_point_field), ALLOCATABLE, TARGET :: fs_inc_store(:, :, :)
   LOGICAL, ALLOCATABLE :: fs_kept(:, :)                 ! (nfreq, nbeta)
   COMPLEX(8), ALLOCATABLE :: fs_psi(:, :, :)            ! (n_fs, 6*nbody, nneed), psi_k at the inner points, as psi_store
   TYPE(pw_set), ALLOCATABLE :: fs_pw_psi(:, :)          ! (6*nbody, nneed), partial waves of each psi_k
   LOGICAL, ALLOCATABLE :: fs_psi_kept(:)                ! (nneed)
   ! Second-order motion of the free bodies (qtf_motion = 1, Phase 6c); route 1 is the indirect force (.11d, .11s), route 2
   ! the direct one (.12d, .12s, qtf_direct = 1)
   LOGICAL :: mot_active = .FALSE.                       ! the QTF_Motion files are written
   INTEGER :: mot_nroute = 0                             ! 1, or 2 with qtf_direct = 1
   REAL(8), ALLOCATABLE :: mot_mass(:, :), mot_c(:, :), mot_k(:, :), mot_b(:, :)   ! (6*nbody, 6*nbody): M, C, K_ext, BLNR
   COMPLEX(8), ALLOCATABLE :: xi2(:, :, :, :, :)         ! (6*nbody, nfreq, nfreq, heading pairs, mot_nroute), difference
   COMPLEX(8), ALLOCATABLE :: xi2_sum(:, :, :, :, :)     ! as xi2, sum; both e^{-i Omega t}, unscaled
   ! Per Omega of the list, index 0 for Omega = 0 (the difference diagonal): rcond >= 1e-12 (its rows are written), the
   ! rcond, and the modes without restoring when it is not
   LOGICAL, ALLOCATABLE :: mot_ok(:)                     ! (0:nneed)
   REAL(8), ALLOCATABLE :: mot_rcond(:)                  ! (0:nneed)
   LOGICAL, ALLOCATABLE :: mot_free(:, :)                ! (6*nbody, 0:nneed)
   LOGICAL, ALLOCATABLE :: mot_solved(:)                 ! (0:nneed), the Omega has rows and was solved
   ! psi_j at the field points of QTFField per Omega of the list (qtf_field = 1 with qtf_direct = 1), for the radiated wave
   COMPLEX(8), ALLOCATABLE :: mot_psi_pts(:, :, :)       ! (npt, 6*nbody, nneed), e^{-i Omega t} at Omega > 0
   LOGICAL, ALLOCATABLE :: mot_psi_kept(:)               ! (nneed)
   ! qtf_fs_inner = source (item 2 of the field-pass speed-ups): the inner integral without phi_zz, which needs grad_h psi_k
   ! at the inner points too, (n_fs, 2, 6*nbody, nneed), kept by qtf_keep_fs_psi
   LOGICAL :: fs_source = .FALSE.
   COMPLEX(8), ALLOCATABLE :: fs_gpsi(:, :, :, :)
   ! qtf_stored_mb counts fs_gpsi: set by qtf_set_fs_inner before the estimate (the store is allocated later, by qtf_init)
   LOGICAL :: estimate_gpsi = .FALSE.
   ! Non-finite results of the potential part, over the run (reset by qtf_init, kept after qtf_finish for the exit code)
   INTEGER :: n_bad_psi = 0                              ! difference frequencies whose psi has NaN or Inf
   INTEGER :: n_bad_pairs = 0                            ! frequency and heading pairs whose T_pot has NaN or Inf

CONTAINS

   SUBROUTINE qtf_init(opts, nfreq, beta_in_deg, nwl, n_panels, nbody, outputdir, omegas_main, fs_pts, fs_centre, fs_radius, &
                       fs_depth, fs_mmax, fs_nev, fs_np_circle, geom, g, irsp, mass, c_hyd, k_ext, b_ext)
      ! Allocates the whole store for up to nfreq frequencies, the headings beta_in_deg, n_panels panels and nbody bodies,
      ! so that a run that does not fit stops here, before the first-order solves. outputdir holds Wamit_format.
      ! omegas_main (rad/s), the nfreq frequencies of the QTF grid, is needed for qtf_terms = body: the store then also holds
      ! psi at every difference frequency. qtf_terms = full needs as well the inner free-surface points fs_pts, the circle
      ! (fs_centre, fs_radius), the depth fs_depth and the mode counts fs_mmax, fs_nev of the partial waves, for the
      ! free-surface store, and fs_np_circle, the chords of the polygon a generated inner mesh ends on (build_fs_inner; 0 or
      ! absent for a mesh from FSInnerMesh.pnl, which gets no sliver).
      ! geom, g and irsp (Phase 7c-2, the sum frequency of the direct solve): with irsp = 0 the sum pairs at or above the
      ! estimate of the lowest irregular frequency (irregular_frequency_estimate) are left out of the list and of the pairs,
      ! with a warning. Without them no such cutoff applies.
      ! mass, c_hyd, k_ext, b_ext (Phase 6c, needed with qtf_motion = 1 and a free body): the body mass matrix, the
      ! hydrostatic restoring, the external stiffness and the external linear damping, each (6 nbody, 6 nbody) and
      ! block-diagonal over the bodies (MATX, CRS, KSTF, BLNR, or their _MULTI forms)
      TYPE(so_options), INTENT(IN) :: opts
      INTEGER, INTENT(IN) :: nfreq, nwl, n_panels, nbody
      REAL(8), INTENT(IN) :: beta_in_deg(:)
      CHARACTER(LEN=*), INTENT(IN) :: outputdir
      REAL(8), INTENT(IN), OPTIONAL :: omegas_main(:)
      TYPE(cs_points), INTENT(IN), OPTIONAL :: fs_pts
      REAL(8), INTENT(IN), OPTIONAL :: fs_centre(2), fs_radius, fs_depth
      INTEGER, INTENT(IN), OPTIONAL :: fs_mmax, fs_nev, fs_np_circle
      TYPE(surface_geom), INTENT(IN), OPTIONAL :: geom
      REAL(8), INTENT(IN), OPTIONAL :: g
      INTEGER, INTENT(IN), OPTIONAL :: irsp
      REAL(8), INTENT(IN), OPTIONAL :: mass(:, :), c_hyd(:, :), k_ext(:, :), b_ext(:, :)
      INTEGER :: st, i, b
      REAL(8) :: sum_top

      active = opts%qtf > 0
      IF (.NOT. active) RETURN
      ! The sum is on for qtf = 2 or 3 with qtf_terms = body or full, with or without the direct solve (Phase 6b: with
      ! qtf_direct = 0 the sum QTF is the indirect one, .10s + T+_pot, and no direct solve runs)
      diff_on = opts%qtf == 1 .OR. opts%qtf == 3
      sum_on = opts%qtf >= 2 .AND. opts%qtf_terms >= QTF_TERMS_BODY
      body_active = opts%qtf_terms >= QTF_TERMS_BODY .AND. (diff_on .OR. sum_on)
      ! The motion needs the potential part (validate_second_order checks qtf_terms and qtf) and a free body: a fixed one
      ! does not move, so nothing is computed for it (plan review focus 2), and this is the one log line that says so
      IF (opts%qtf_motion == 1 .AND. opts%motion_fixed) &
         WRITE (*, '(a)') ' QTF note: qtf_motion = 1 with motion_state = fixed: no second-order motion is computed.'
      ! The stores of a previous run first (free_motion_stores clears mot_active), then the one assignment of this run,
      ! which sum_warnings below reads already
      CALL free_motion_stores()
      mot_active = opts%qtf_motion == 1 .AND. .NOT. opts%motion_fixed .AND. body_active
      IF (body_active) THEN
         IF (.NOT. PRESENT(omegas_main)) THEN
            WRITE (*, '(a)') ' QTF: qtf_init needs the main-grid frequencies for qtf_terms = body.'
            ERROR STOP 1
         END IF
         ! The difference frequencies first, so that a grid they do not fit stops before anything is allocated
         IF (ALLOCATED(psi_need)) DEALLOCATE (psi_need)
         IF (ALLOCATED(psi_reuse)) DEALLOCATE (psi_reuse)
         psi_tol = qtf_psi_tolerance(omegas_main)
         CALL check_grid(omegas_main)
         max_sw = opts%qtf_max_sw
         sum_cut_irr = HUGE(1D0)
         sum_top = 0D0
         IF (sum_on) THEN
            CALL qtf_sum_limits(opts, geom, g, irsp, sum_cut_irr, sum_top)
            CALL sum_warnings(omegas_main)
         END IF
         IF (ALLOCATED(psi_is_diff)) DEALLOCATE (psi_is_diff, psi_is_sum)
         CALL qtf_omega_list(omegas_main, opts%qtf_max_dw, sum_top, psi_tol, diff_on, sum_on, psi_need, psi_reuse, &
                             psi_is_diff, psi_is_sum)
         CALL check_pairs(omegas_main)
      END IF
      kind_qtf = opts%qtf
      max_dw = opts%qtf_max_dw
      all_headings_pairs = MERGE(1, 0, opts%mean_drift == 2)
      nbeta = SIZE(beta_in_deg)
      nkept = 0
      outdir = outputdir
      IF (ALLOCATED(store)) DEALLOCATE (store)
      IF (ALLOCATED(store_wl)) DEALLOCATE (store_wl)
      IF (ALLOCATED(store_oufr)) DEALLOCATE (store_oufr)
      ALLOCATE (store(nfreq, nbeta), store_wl(nwl, nfreq, nbeta), store_oufr(nfreq), STAT=st)
      IF (st /= 0) CALL out_of_memory()
      DO b = 1, nbeta
         DO i = 1, nfreq
            ASSOCIATE (c => store(i, b))
               ALLOCATE (c%phi_i(n_panels), c%phi_p(n_panels), c%grad_i(3, n_panels), c%grad_p(3, n_panels), &
                         c%xi(6, nbody), STAT=st)
               IF (st /= 0) CALL out_of_memory()
               ! Touch the arrays, so that memory that is only promised (overcommit) is taken now
               c%phi_i = (0D0, 0D0)
               c%phi_p = (0D0, 0D0)
               c%grad_i = (0D0, 0D0)
               c%grad_p = (0D0, 0D0)
               c%xi = (0D0, 0D0)
            END ASSOCIATE
         END DO
      END DO
      store_wl = (0D0, 0D0)
      beta_deg = beta_in_deg

      IF (.NOT. body_active) THEN
         IF (ALLOCATED(psi_need)) DEALLOCATE (psi_need)
         IF (ALLOCATED(psi_reuse)) DEALLOCATE (psi_reuse)
         IF (ALLOCATED(psi_is_diff)) DEALLOCATE (psi_is_diff, psi_is_sum)
      END IF
      IF (ALLOCATED(psi_store)) DEALLOCATE (psi_store)
      IF (ALLOCATED(psi_kept)) DEALLOCATE (psi_kept)
      IF (ALLOCATED(tpot)) DEALLOCATE (tpot)
      IF (ALLOCATED(tdir)) DEALLOCATE (tdir)
      IF (ALLOCATED(tdir_kept)) DEALLOCATE (tdir_kept)
      IF (ALLOCATED(tdir_i2)) DEALLOCATE (tdir_i2, tdir_d)
      IF (ALLOCATED(tdir_fs)) DEALLOCATE (tdir_fs)
      CALL free_sum_stores()
      CALL free_control_surface()
      CALL free_free_surface()
      n_bad_psi = 0
      n_bad_pairs = 0
      pot_terms = body_active .AND. diff_on .AND. opts%write_terms >= 1
      pot_sum_terms = body_active .AND. sum_on .AND. opts%write_terms >= 1
      IF (body_active) THEN
         ALLOCATE (psi_store(n_panels, 6*nbody, SIZE(psi_need)), psi_kept(SIZE(psi_need)), STAT=st)
         IF (st /= 0) CALL out_of_memory(body=.TRUE.)
         psi_store = (0D0, 0D0)
         psi_kept = .FALSE.
      END IF
      IF (body_active .AND. diff_on) THEN
         ALLOCATE (tpot(6*nbody, nfreq, nfreq, heading_pairs(nbeta)), STAT=st)
         IF (st /= 0) CALL out_of_memory(body=.TRUE.)
         tpot = (0D0, 0D0)
      END IF
      IF (body_active .AND. opts%qtf_terms == QTF_TERMS_FULL) CALL init_free_surface()
      direct_active = body_active .AND. opts%qtf_direct == 1
      IF (direct_active) THEN
         ALLOCATE (tdir_kept(SIZE(psi_need)), STAT=st)
         IF (st /= 0) CALL out_of_memory(body=.TRUE.)
         tdir_kept = .FALSE.
      END IF
      IF (direct_active .AND. diff_on) THEN
         ALLOCATE (tdir(6*nbody, nfreq, nfreq, heading_pairs(nbeta)), STAT=st)
         IF (st /= 0) CALL out_of_memory(body=.TRUE.)
         tdir = (0D0, 0D0)
      END IF
      ! The indirect sum load tpot_sum (R1) comes with the sum, with or without the direct solve (Phase 6b); the direct
      ! one, tdir_sum (R6), only with qtf_direct = 1; with write_terms >= 1 also the parts of tdir_sum for
      ! QTF_DirectTermsSum.txt
      IF (body_active .AND. sum_on) THEN
         ALLOCATE (tpot_sum(6*nbody, nfreq, nfreq, heading_pairs(nbeta)), STAT=st)
         IF (st /= 0) CALL out_of_memory(body=.TRUE.)
         tpot_sum = (0D0, 0D0)
      END IF
      IF (direct_active .AND. sum_on) THEN
         ALLOCATE (tdir_sum(6*nbody, nfreq, nfreq, heading_pairs(nbeta)), STAT=st)
         IF (st /= 0) CALL out_of_memory(body=.TRUE.)
         tdir_sum = (0D0, 0D0)
      END IF
      sum_terms = direct_active .AND. sum_on .AND. opts%write_terms >= 1
      IF (sum_terms) THEN
         ALLOCATE (tdir_sum_i2(6*nbody, nfreq, nfreq, heading_pairs(nbeta)), tdir_sum_d(6*nbody, nfreq, nfreq, &
                   heading_pairs(nbeta)), STAT=st)
         IF (st /= 0) CALL out_of_memory(body=.TRUE.)
         tdir_sum_i2 = (0D0, 0D0)
         tdir_sum_d = (0D0, 0D0)
         IF (fs_active) THEN
            ALLOCATE (tdir_sum_fs(6*nbody, nfreq, nfreq, heading_pairs(nbeta)), STAT=st)
            IF (st /= 0) CALL out_of_memory(body=.TRUE., fs=.TRUE.)
            tdir_sum_fs = (0D0, 0D0)
         END IF
      END IF
      direct_terms = direct_active .AND. diff_on .AND. opts%write_terms >= 1
      IF (direct_terms) THEN
         ALLOCATE (tdir_i2(6*nbody, nfreq, nfreq, heading_pairs(nbeta)), tdir_d(6*nbody, nfreq, nfreq, heading_pairs(nbeta)), &
                   STAT=st)
         IF (st /= 0) CALL out_of_memory(body=.TRUE.)
         tdir_i2 = (0D0, 0D0)
         tdir_d = (0D0, 0D0)
         IF (fs_active) THEN
            ALLOCATE (tdir_fs(6*nbody, nfreq, nfreq, heading_pairs(nbeta)), STAT=st)
            IF (st /= 0) CALL out_of_memory(body=.TRUE., fs=.TRUE.)
            tdir_fs = (0D0, 0D0)
         END IF
      END IF
      IF (mot_active) CALL init_motion()

   CONTAINS

      SUBROUTINE init_motion()
         ! The body matrices and the stores of the second-order motion (qtf_motion = 1, Phase 6c): xi2 with the difference
         ! QTF, xi2_sum with the sum, one slice per heading pair and route
         INTEGER :: nmode

         nmode = 6*nbody
         IF (.NOT. (PRESENT(mass) .AND. PRESENT(c_hyd) .AND. PRESENT(k_ext) .AND. PRESENT(b_ext))) THEN
            WRITE (*, '(a)') ' QTF: qtf_init needs the mass, restoring, stiffness and damping matrices for qtf_motion = 1.'
            ERROR STOP 1
         END IF
         IF (ANY(SHAPE(mass) /= [nmode, nmode]) .OR. ANY(SHAPE(c_hyd) /= [nmode, nmode]) .OR. &
             ANY(SHAPE(k_ext) /= [nmode, nmode]) .OR. ANY(SHAPE(b_ext) /= [nmode, nmode])) THEN
            WRITE (*, '(a)') ' QTF: the body matrices of qtf_motion = 1 must be 6 nbody x 6 nbody.'
            ERROR STOP 1
         END IF
         mot_mass = mass
         mot_c = c_hyd
         mot_k = k_ext
         mot_b = b_ext
         mot_nroute = MERGE(2, 1, direct_active)
         IF (diff_on) THEN
            ALLOCATE (xi2(nmode, nfreq, nfreq, heading_pairs(nbeta), mot_nroute), STAT=st)
            IF (st /= 0) CALL out_of_memory(body=.TRUE.)
            xi2 = (0D0, 0D0)
         END IF
         IF (sum_on) THEN
            ALLOCATE (xi2_sum(nmode, nfreq, nfreq, heading_pairs(nbeta), mot_nroute), STAT=st)
            IF (st /= 0) CALL out_of_memory(body=.TRUE.)
            xi2_sum = (0D0, 0D0)
         END IF
         ALLOCATE (mot_ok(0:SIZE(psi_need)), mot_rcond(0:SIZE(psi_need)), mot_free(nmode, 0:SIZE(psi_need)), &
                   mot_solved(0:SIZE(psi_need)))
         mot_ok = .FALSE.
         mot_rcond = 0D0
         mot_free = .FALSE.
         mot_solved = .FALSE.
      END SUBROUTINE init_motion

      SUBROUTINE init_free_surface()
         ! The free-surface store of qtf_terms = full, touched (default initialisation and zeros) so that it is taken now
         INTEGER :: m, k

         IF (.NOT. (PRESENT(fs_pts) .AND. PRESENT(fs_centre) .AND. PRESENT(fs_radius) .AND. PRESENT(fs_depth) .AND. &
                    PRESENT(fs_mmax) .AND. &
                    PRESENT(fs_nev))) THEN
            WRITE (*, '(a)') ' QTF: qtf_init needs the inner free surface, the circle and the mode counts for qtf_terms = full.'
            ERROR STOP 1
         END IF
         fs_store_pts = fs_pts
         fs_store_centre = fs_centre
         fs_store_radius = fs_radius
         fs_store_depth = fs_depth
         fs_store_mmax = fs_mmax
         fs_store_nev = fs_nev
         fs_store_np = 0
         IF (PRESENT(fs_np_circle)) fs_store_np = fs_np_circle
         ALLOCATE (fs_f(fs_pts%n, nfreq, nbeta), fs_pw(nfreq, nbeta), fs_kept(nfreq, nbeta), &
                   fs_psi(fs_pts%n, 6*nbody, SIZE(psi_need)), fs_pw_psi(6*nbody, SIZE(psi_need)), &
                   fs_psi_kept(SIZE(psi_need)), STAT=st)
         IF (st /= 0) CALL out_of_memory(body=.TRUE., fs=.TRUE.)
         fs_f = fs_point_field()
         fs_psi = (0D0, 0D0)
         DO b = 1, nbeta
            DO i = 1, nfreq
               CALL alloc_pw(fs_pw(i, b))
            END DO
         END DO
         DO m = 1, SIZE(psi_need)
            DO k = 1, 6*nbody
               CALL alloc_pw(fs_pw_psi(k, m))
            END DO
         END DO
         fs_kept = .FALSE.
         fs_psi_kept = .FALSE.
         fs_source = opts%qtf_fs_inner == 2
         IF (fs_source) THEN
            ALLOCATE (fs_gpsi(fs_pts%n, 2, 6*nbody, SIZE(psi_need)), STAT=st)
            IF (st /= 0) CALL out_of_memory(body=.TRUE., fs=.TRUE.)
            fs_gpsi = (0D0, 0D0)
         END IF
         fs_active = .TRUE.
      END SUBROUTINE init_free_surface

      SUBROUTINE alloc_pw(pw)
         TYPE(pw_set), INTENT(OUT) :: pw

         pw%mmax = fs_store_mmax
         pw%nev = fs_store_nev
         pw%h = fs_store_depth
         pw%radius = fs_store_radius
         pw%centre = fs_store_centre
         ALLOCATE (pw%kev(fs_store_nev), pw%a(-fs_store_mmax:fs_store_mmax, 0:fs_store_nev), STAT=st)
         IF (st /= 0) CALL out_of_memory(body=.TRUE., fs=.TRUE.)
         pw%kev = 0D0
         pw%a = (0D0, 0D0)
      END SUBROUTINE alloc_pw

      SUBROUTINE check_grid(w)
         ! Two frequencies of the grid that are equal to within psi_tol but not exactly would give an Omega that is neither
         ! the diagonal nor a difference frequency of psi
         REAL(8), INTENT(IN) :: w(:)
         INTEGER :: i, j

         DO i = 1, SIZE(w)
            DO j = i + 1, SIZE(w)
               IF (w(i) /= w(j) .AND. ABS(w(i) - w(j)) <= psi_tol) THEN
                  WRITE (*, '(a,I0,a,I0,a,ES22.14,a,ES22.14,a)') ' QTF: frequencies ', i, ' and ', j, ' of the grid,', &
                     w(i), ' and', w(j), ' rad/s, are equal to within rounding but not exactly.'
                  WRITE (*, '(a)') ' Remove the duplicate from the frequency list (qtf_terms = body needs distinct frequencies).'
                  ERROR STOP 1
               END IF
            END DO
         END DO
      END SUBROUTINE check_grid

      SUBROUTINE check_pairs(w)
         ! Every pair of the difference QTF finds its Omega among the difference frequencies within psi_tol, the tolerance
         ! potential_pair is given: a stop here costs no solve
         REAL(8), INTENT(IN) :: w(:)
         REAL(8) :: dw
         INTEGER :: i, j

         DO i = 1, SIZE(w)
            DO j = 1, SIZE(w)
               IF (sum_on .AND. j >= i) THEN
                  IF (sum_pair_kept(w(i) + w(j))) THEN
                     IF (qtf_psi_index(w(i) + w(j)) == 0) THEN
                        WRITE (*, '(a,ES22.14,a,ES14.6,a,ES14.6,a)') ' QTF: Omega =', w(i) + w(j), ' rad/s of the pair', w(i), &
                           ',', w(j), ' rad/s is not among the frequencies of psi.'
                        WRITE (*, '(a)') ' This is a bug: please report it with the frequency list.'
                        ERROR STOP 1
                     END IF
                  END IF
               END IF
               IF (.NOT. diff_on) CYCLE
               dw = w(i) - w(j)
               IF (dw <= 0D0) CYCLE
               IF (opts%qtf_max_dw > 0D0 .AND. dw > opts%qtf_max_dw*(1D0 + 1D-9)) CYCLE
               IF (qtf_psi_index(dw) == 0) THEN
                  WRITE (*, '(a,ES22.14,a,ES14.6,a,ES14.6,a)') ' QTF: Omega =', dw, ' rad/s of the pair', w(i), ',', w(j), &
                     ' rad/s is not among the difference frequencies of psi.'
                  WRITE (*, '(a)') ' This is a bug: please report it with the frequency list.'
                  ERROR STOP 1
               END IF
            END DO
         END DO
      END SUBROUTINE check_pairs

      SUBROUTINE sum_warnings(w)
         ! The log of what the sum leaves out (rulings R2 and R5, docs/developer_documentation.md, Phase 7c-2)
         REAL(8), INTENT(IN) :: w(:)
         INTEGER :: i, j, nleft, ntot
         CHARACTER(LEN=40) :: left_out, also   ! the files of the sum this run writes: the first, and the one after it

         left_out = 'QTF_Indirect.11s'
         also = ''
         IF (opts%qtf_direct == 1) THEN
            left_out = 'QTF_Direct.12s'
            also = ', and of QTF_Indirect.11s'
         END IF
         IF (opts%qtf_terms == QTF_TERMS_BODY) THEN
            WRITE (*, '(a)') ' QTF warning: qtf_terms = body with the sum frequency leaves out the free-surface integral. The'
            WRITE (*, '(a)') ' sum QTF without it is not usable (WAMIT advises against it); use qtf_terms = full.'
            IF (mot_active) WRITE (*, '(a)') ' The same holds for QTF_Motion.11s, the motion from that force.'
         END IF
         IF (sum_cut_irr < HUGE(1D0)) THEN
            nleft = 0
            ntot = 0
            DO i = 1, SIZE(w)
               DO j = i, SIZE(w)
                  ntot = ntot + 1
                  IF (w(i) + w(j) >= sum_cut_irr) nleft = nleft + 1
               END DO
            END DO
            IF (nleft > 0) THEN
               WRITE (*, '(a,F8.4,a,I0,a,I0,a)') ' QTF warning: the sum frequencies at or above', sum_cut_irr, &
                  ' rad/s, the estimate of the lowest irregular frequency, are left out of '//TRIM(left_out)//': ', nleft, ' of ', &
                  ntot, ' pairs (i <= j)'//TRIM(also)//'.'
               WRITE (*, '(a)') ' Set If_remove_irr_freq = 1 to keep them. The estimate is a lower bound for a wall-sided hull'// &
                  ' only; near the true irregular frequency a row can still show a spike.'
               IF (mot_active) WRITE (*, '(a)') ' QTF_Motion.11s leaves out the same pairs.'
            END IF
         END IF
      END SUBROUTINE sum_warnings

   END SUBROUTINE qtf_init

   LOGICAL FUNCTION sum_pair_kept(omega)
      ! The sum Omega = omega_i + omega_j is inside the cutoff of the sum: not above qtf_max_sw (inclusive, to rounding, as
      ! qtf_max_dw) and below the estimate of the lowest irregular frequency
      REAL(8), INTENT(IN) :: omega

      sum_pair_kept = omega < sum_cut_irr
      IF (max_sw > 0D0 .AND. omega > max_sw*(1D0 + 1D-9)) sum_pair_kept = .FALSE.
   END FUNCTION sum_pair_kept

   SUBROUTINE qtf_sum_limits(opts, geom, g, irsp, cut_irr, sum_top)
      ! The cutoffs of the sum frequencies (Phase 7c-2): cut_irr, the estimate of the lowest irregular frequency, when
      ! irsp = 0 and geom and g are given (HUGE otherwise: no cutoff); sum_top, the largest sum Omega that has a psi, for
      ! qtf_omega_list: the smaller of qtf_max_sw (when above 0) and cut_irr (a hair below it, as sums at it are left out),
      ! 0 when neither applies. qtf_init and the driver log call it, so that the list is one.
      TYPE(so_options), INTENT(IN) :: opts
      TYPE(surface_geom), INTENT(IN), OPTIONAL :: geom
      REAL(8), INTENT(IN), OPTIONAL :: g
      INTEGER, INTENT(IN), OPTIONAL :: irsp
      REAL(8), INTENT(OUT) :: cut_irr, sum_top

      cut_irr = HUGE(1D0)
      sum_top = 0D0
      IF (PRESENT(geom) .AND. PRESENT(g) .AND. PRESENT(irsp)) THEN
         IF (irsp == 0) cut_irr = irregular_frequency_estimate(geom, g)
      END IF
      IF (opts%qtf_max_sw > 0D0) sum_top = opts%qtf_max_sw
      IF (cut_irr < HUGE(1D0)) THEN
         IF (sum_top > 0D0) THEN
            sum_top = MIN(sum_top, cut_irr*(1D0 - 1D-9))
         ELSE
            sum_top = cut_irr*(1D0 - 1D-9)
         END IF
      END IF
   END SUBROUTINE qtf_sum_limits

   REAL(8) FUNCTION irregular_frequency_estimate(geom, g)
      ! Estimate of the lowest irregular frequency (rad/s) of the closed waterline loops of geom, for the cutoff of the sum
      ! frequencies (ruling R2): the minimum over the loops of every body of sqrt(g j01 / a_eq), with j01 = 2.404826 the first
      ! zero of J_0 and a_eq = sqrt(A / pi), A the shoelace area of the loop. Why: for a wall-sided hull the interior
      ! eigenfunction J_m(kr) sinh k(z + d) has omega^2 = g k coth(k d) with coth >= 1 [derived], so omega^2 >= g k;
      ! by Faber-Krahn the disc of equal area has the smallest first Dirichlet eigenvalue (Polya and Szego, 1951), so
      ! k >= j01 / a_eq and the value is a lower bound. A hull that is not wall-sided (flare, rounded bilge) is not covered.
      ! It gives 1.536 rad/s for a cylinder of radius 10 m (exact 1.548) and 1.983 for the OC4 hull. HUGE(1D0) if there is
      ! no closed loop. A loop is closed when wl_next leads back to its first segment; open chains are skipped.
      TYPE(surface_geom), INTENT(IN) :: geom
      REAL(8), INTENT(IN) :: g
      REAL(8), PARAMETER :: J01 = 2.404826D0
      LOGICAL, ALLOCATABLE :: seen(:)
      REAL(8) :: area, a_eq
      INTEGER :: iw, j, nseg

      irregular_frequency_estimate = HUGE(1D0)
      IF (.NOT. ALLOCATED(geom%wl_next) .OR. geom%nwl < 1) RETURN
      ALLOCATE (seen(geom%nwl))
      seen = .FALSE.
      DO iw = 1, geom%nwl
         IF (seen(iw)) CYCLE
         ! Walk wl_next from iw: a closed loop returns to iw, a chain ends at 0 or runs into a segment already seen
         area = 0D0
         nseg = 0
         j = iw
         DO
            seen(j) = .TRUE.
            area = area + geom%wl_x1(1, j)*geom%wl_x2(2, j) - geom%wl_x2(1, j)*geom%wl_x1(2, j)
            nseg = nseg + 1
            j = geom%wl_next(j)
            IF (j == 0 .OR. nseg > geom%nwl) EXIT
            IF (j == iw) EXIT
            IF (seen(j)) THEN
               j = 0
               EXIT
            END IF
         END DO
         IF (j /= iw) CYCLE
         area = 0.5D0*ABS(area)
         IF (area <= 0D0) CYCLE
         a_eq = SQRT(area/ACOS(-1D0))
         irregular_frequency_estimate = MIN(irregular_frequency_estimate, SQRT(g*J01/a_eq))
      END DO
   END FUNCTION irregular_frequency_estimate

   INTEGER FUNCTION sum_pair_index(i, j)
      ! The pair rule of the sum QTF, as difference_pair_index is that of the difference QTF: 0 when omega_i + omega_j is
      ! beyond the cutoff of the sum (qtf_max_sw, inclusive to rounding, or the estimate of the lowest irregular frequency
      ! with IRSP = 0, exclusive: sums at or above it are left out), else the index of omega_i + omega_j in psi_need
      ! (qtf_psi_index), or -1 if it is not there. 0 as well when the sum is not on.
      INTEGER, INTENT(IN) :: i, j
      REAL(8) :: w

      sum_pair_index = 0
      IF (.NOT. sum_on) RETURN
      w = store(i, 1)%omega + store(j, 1)%omega
      IF (.NOT. sum_pair_kept(w)) RETURN
      sum_pair_index = qtf_psi_index(w)
      IF (sum_pair_index == 0) sum_pair_index = -1
   END FUNCTION sum_pair_index

   REAL(8) PURE FUNCTION qtf_psi_tolerance(omegas_main)
      ! The tolerance of the difference frequencies of the grid omegas_main (rad/s): qtf_radiation_list puts every Omega
      ! within it of its psi, and an Omega is found (qtf_psi_index) and accepted (potential_pair) within it. Rounding only
      REAL(8), INTENT(IN) :: omegas_main(:)

      qtf_psi_tolerance = PSI_REL_TOL*MAXVAL(omegas_main)
   END FUNCTION qtf_psi_tolerance

   SUBROUTINE qtf_set_control_surface(pts)
      ! qtf_b2 = control: keeps the points of C and F and allocates the fields there for the store of qtf_init. Call after
      ! qtf_init, before the frequency loop, only with qtf_terms = body or full and qtf = 1, 2 or 3.
      TYPE(cs_points), INTENT(IN) :: pts
      INTEGER :: st

      IF (.NOT. body_active) THEN
         WRITE (*, '(a)') ' QTF: qtf_set_control_surface needs qtf_terms = body or full and qtf = 1, 2 or 3.'
         ERROR STOP 1
      END IF
      CALL free_control_surface()
      cs = pts
      ALLOCATE (cs_u(3, cs%n, SIZE(store, 1), nbeta), cs_u_kept(SIZE(store, 1), nbeta), &
                cs_gpsi(3, cs%n, SIZE(psi_store, 2), SIZE(psi_need)), cs_psi_kept(SIZE(psi_need)), STAT=st)
      IF (st /= 0) CALL out_of_memory(body=.TRUE.)
      cs_u = (0D0, 0D0)
      cs_gpsi = (0D0, 0D0)
      cs_u_kept = .FALSE.
      cs_psi_kept = .FALSE.
      b2_control = .TRUE.
   END SUBROUTINE qtf_set_control_surface

   SUBROUTINE qtf_set_edge_weights(w)
      ! qtf_b2 = control with qtf_edge_weights: the panel weights w(n_panels) of the B2 hull term (SharpEdge). Call after
      ! qtf_set_control_surface.
      REAL(8), INTENT(IN) :: w(:)

      IF (.NOT. b2_control .OR. SIZE(w) /= SIZE(psi_store, 1)) THEN
         WRITE (*, '(a)') ' QTF: qtf_set_edge_weights needs qtf_b2 = control and one weight per panel.'
         ERROR STOP 1
      END IF
      edge_w = w
   END SUBROUTINE qtf_set_edge_weights

   SUBROUTINE free_control_surface()
      b2_control = .FALSE.
      IF (ALLOCATED(edge_w)) DEALLOCATE (edge_w)
      IF (ALLOCATED(cs_u)) DEALLOCATE (cs_u)
      IF (ALLOCATED(cs_u_kept)) DEALLOCATE (cs_u_kept)
      IF (ALLOCATED(cs_gpsi)) DEALLOCATE (cs_gpsi)
      IF (ALLOCATED(cs_psi_kept)) DEALLOCATE (cs_psi_kept)
      cs%n = 0
   END SUBROUTINE free_control_surface

   SUBROUTINE qtf_keep_cs_velocity(ib, u)
      ! Keeps u(3, cs%n), the total velocity at the points of C and F of heading ib of the frequency kept last by qtf_keep
      INTEGER, INTENT(IN) :: ib
      COMPLEX(8), INTENT(IN) :: u(:, :)

      IF (.NOT. b2_control) RETURN
      IF (nkept == 0 .OR. ib < 1 .OR. ib > nbeta .OR. SIZE(u, 1) /= 3 .OR. SIZE(u, 2) /= cs%n) THEN
         WRITE (*, '(a,I0)') ' QTF: qtf_keep_cs_velocity called out of order or with another size, heading ', ib
         ERROR STOP 1
      END IF
      cs_u(:, :, nkept, ib) = u
      cs_u_kept(nkept, ib) = .TRUE.
   END SUBROUTINE qtf_keep_cs_velocity

   SUBROUTINE qtf_keep_cs_psi(omega, grad)
      ! Keeps grad(3, cs%n, 6*nbody), grad psi_k at the points of C and F, of the difference frequency omega
      REAL(8), INTENT(IN) :: omega
      COMPLEX(8), INTENT(IN) :: grad(:, :, :)
      INTEGER :: m

      IF (.NOT. b2_control) RETURN
      m = qtf_psi_index(omega)
      IF (m == 0 .OR. SIZE(grad, 1) /= 3 .OR. SIZE(grad, 2) /= cs%n .OR. SIZE(grad, 3) /= SIZE(cs_gpsi, 3)) THEN
         WRITE (*, '(a,ES14.6,a)') ' QTF: grad psi on the control surface at omega =', omega, &
            ' rad/s is not of a difference frequency, or has another size.'
         ERROR STOP 1
      END IF
      cs_gpsi(:, :, :, m) = grad
      cs_psi_kept(m) = .TRUE.
   END SUBROUTINE qtf_keep_cs_psi

   SUBROUTINE free_free_surface()
      fs_active = .FALSE.
      IF (ALLOCATED(fs_f)) DEALLOCATE (fs_f)
      IF (ALLOCATED(fs_pw)) DEALLOCATE (fs_pw)
      IF (ALLOCATED(fs_inc_store)) DEALLOCATE (fs_inc_store)
      IF (ALLOCATED(fs_kept)) DEALLOCATE (fs_kept)
      IF (ALLOCATED(fs_psi)) DEALLOCATE (fs_psi)
      IF (ALLOCATED(fs_pw_psi)) DEALLOCATE (fs_pw_psi)
      IF (ALLOCATED(fs_psi_kept)) DEALLOCATE (fs_psi_kept)
      IF (ALLOCATED(fs_gpsi)) DEALLOCATE (fs_gpsi)
      fs_source = .FALSE.
      IF (ALLOCATED(fs_store_pts%x)) DEALLOCATE (fs_store_pts%x)
      IF (ALLOCATED(fs_store_pts%nrm)) DEALLOCATE (fs_store_pts%nrm)
      IF (ALLOCATED(fs_store_pts%w)) DEALLOCATE (fs_store_pts%w)
      IF (ALLOCATED(fs_store_pts%psize)) DEALLOCATE (fs_store_pts%psize)
      IF (ALLOCATED(fs_store_pts%kind)) DEALLOCATE (fs_store_pts%kind)
      IF (ALLOCATED(fs_store_pts%body)) DEALLOCATE (fs_store_pts%body)
      fs_store_pts%n = 0
   END SUBROUTINE free_free_surface

   LOGICAL FUNCTION qtf_fs_active()
      ! The free-surface store of qtf_terms = full is allocated
      qtf_fs_active = fs_active
   END FUNCTION qtf_fs_active

   LOGICAL FUNCTION pw_fits(pw)
      ! pw has the mode counts, the circle and the depth of the free-surface store
      TYPE(pw_set), INTENT(IN) :: pw

      pw_fits = ALLOCATED(pw%a) .AND. ALLOCATED(pw%kev)
      IF (pw_fits) pw_fits = pw%mmax == fs_store_mmax .AND. pw%nev == fs_store_nev
      ! The circle and the depth of the store, to rounding
      IF (pw_fits) pw_fits = near_equal(pw%radius, fs_store_radius) .AND. near_equal(pw%h, fs_store_depth) .AND. &
                             NORM2(pw%centre - fs_store_centre) <= FS_GEOM_TOL*MAX(fs_store_radius, NORM2(fs_store_centre))
   CONTAINS
      LOGICAL FUNCTION near_equal(a, b)
         REAL(8), INTENT(IN) :: a, b
         near_equal = ABS(a - b) <= FS_GEOM_TOL*MAX(ABS(a), ABS(b))
      END FUNCTION near_equal
   END FUNCTION pw_fits

   SUBROUTINE qtf_keep_fs_state(ib, f, pw)
      ! qtf_terms = full: keeps f(n_fs), the total first-order fields at the inner free-surface points, and pw, the partial
      ! waves of the scattered potential outside R, of heading ib of the frequency kept last by qtf_keep
      INTEGER, INTENT(IN) :: ib
      TYPE(fs_point_field), INTENT(IN) :: f(:)
      TYPE(pw_set), INTENT(IN) :: pw

      IF (.NOT. fs_active) RETURN
      IF (nkept == 0 .OR. ib < 1 .OR. ib > nbeta .OR. SIZE(f) /= fs_store_pts%n .OR. .NOT. pw_fits(pw)) THEN
         WRITE (*, '(a,I0)') ' QTF: qtf_keep_fs_state called out of order or with another size, heading ', ib
         ERROR STOP 1
      END IF
      fs_f(:, nkept, ib) = f
      fs_pw(nkept, ib) = pw
      fs_kept(nkept, ib) = .TRUE.
   END SUBROUTINE qtf_keep_fs_state

   SUBROUTINE qtf_keep_fs_psi(omega, psi_pts, pw_psi, ok, gpsi)
      ! qtf_terms = full: keeps psi_pts(n_fs, 6*nbody), psi_k at the inner free-surface points, and pw_psi(6*nbody), the
      ! partial waves of each psi_k, of the difference frequency omega. An omega that is not a difference frequency, or
      ! values of another size, stop the run; with ok present they are not kept and ok is false instead.
      ! gpsi(n_fs, 2, 6*nbody), grad_h psi_k at the points: needed with qtf_fs_inner = source (missing counts as another
      ! size), not kept otherwise.
      REAL(8), INTENT(IN) :: omega
      COMPLEX(8), INTENT(IN) :: psi_pts(:, :)
      TYPE(pw_set), INTENT(IN) :: pw_psi(:)
      LOGICAL, INTENT(OUT), OPTIONAL :: ok
      COMPLEX(8), INTENT(IN), OPTIONAL :: gpsi(:, :, :)
      INTEGER :: m, k
      LOGICAL :: fits

      IF (PRESENT(ok)) ok = .FALSE.
      IF (.NOT. fs_active) RETURN
      m = qtf_psi_index(omega)
      fits = m > 0 .AND. SIZE(psi_pts, 1) == SIZE(fs_psi, 1) .AND. SIZE(psi_pts, 2) == SIZE(fs_psi, 2) .AND. &
             SIZE(pw_psi) == SIZE(fs_pw_psi, 1)
      IF (fits) THEN
         DO k = 1, SIZE(pw_psi)
            fits = fits .AND. pw_fits(pw_psi(k))
         END DO
      END IF
      IF (fits .AND. fs_source) THEN
         fits = PRESENT(gpsi)
         IF (fits) fits = SIZE(gpsi, 1) == SIZE(fs_gpsi, 1) .AND. SIZE(gpsi, 2) == 2 .AND. SIZE(gpsi, 3) == SIZE(fs_gpsi, 3)
      END IF
      IF (.NOT. fits) THEN
         IF (PRESENT(ok)) RETURN
         WRITE (*, '(a,ES14.6,a)') ' QTF: psi on the inner free surface at omega =', omega, &
            ' rad/s is not of a difference frequency, or has another size.'
         ERROR STOP 1
      END IF
      fs_psi(:, :, m) = psi_pts
      fs_pw_psi(:, m) = pw_psi
      IF (fs_source) fs_gpsi(:, :, :, m) = gpsi
      fs_psi_kept(m) = .TRUE.
      IF (PRESENT(ok)) ok = .TRUE.
   END SUBROUTINE qtf_keep_fs_psi

   SUBROUTINE qtf_fs_state(i, ib, f, pw)
      ! The free-surface fields and partial waves kept for frequency i (in the order of qtf_keep) and heading ib. Stops if
      ! they were not kept.
      INTEGER, INTENT(IN) :: i, ib
      TYPE(fs_point_field), ALLOCATABLE, INTENT(OUT) :: f(:)
      TYPE(pw_set), INTENT(OUT) :: pw
      LOGICAL :: kept

      kept = fs_active
      IF (kept) kept = i >= 1 .AND. i <= nkept .AND. ib >= 1 .AND. ib <= nbeta
      IF (kept) kept = fs_kept(i, ib)
      IF (.NOT. kept) THEN
         WRITE (*, '(a,I0,a,I0,a)') ' QTF: the free-surface fields of frequency ', i, ' and heading ', ib, ' were not stored.'
         ERROR STOP 1
      END IF
      f = fs_f(:, i, ib)
      pw = fs_pw(i, ib)
   END SUBROUTINE qtf_fs_state

   SUBROUTINE qtf_fs_psi(omega, psi_pts, pw_psi, ok)
      ! psi_k at the inner free-surface points and the partial waves of each psi_k kept for the difference frequency omega;
      ! ok is false (and nothing is returned) if they were not kept.
      REAL(8), INTENT(IN) :: omega
      COMPLEX(8), ALLOCATABLE, INTENT(OUT) :: psi_pts(:, :)
      TYPE(pw_set), ALLOCATABLE, INTENT(OUT) :: pw_psi(:)
      LOGICAL, INTENT(OUT) :: ok
      INTEGER :: m

      ok = .FALSE.
      IF (.NOT. fs_active) RETURN
      m = qtf_psi_index(omega)
      IF (m == 0) RETURN
      IF (.NOT. fs_psi_kept(m)) RETURN
      psi_pts = fs_psi(:, :, m)
      pw_psi = fs_pw_psi(:, m)
      ok = .TRUE.
   END SUBROUTINE qtf_fs_psi

   INTEGER FUNCTION heading_pairs(nb)
      ! Number of heading pairs of the files for nb headings
      INTEGER, INTENT(IN) :: nb

      heading_pairs = MERGE(nb*nb, nb, all_headings_pairs == 1)
   END FUNCTION heading_pairs

   INTEGER FUNCTION heading_pair(b1, b2)
      ! Index of the heading pair (b1, b2) in tpot; for mean_drift = 1 only b1 = b2 is used
      INTEGER, INTENT(IN) :: b1, b2

      heading_pair = MERGE((b1 - 1)*nbeta + b2, b1, all_headings_pairs == 1)
   END FUNCTION heading_pair

   SUBROUTINE qtf_radiation_list(omegas_main, max_dw_in, tol, need, reuse_index)
      ! The difference frequencies only: qtf_omega_list without the sum, which returns exactly this
      REAL(8), INTENT(IN) :: omegas_main(:), max_dw_in, tol
      REAL(8), ALLOCATABLE, INTENT(OUT) :: need(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: reuse_index(:)
      LOGICAL, ALLOCATABLE :: is_diff(:), is_sum(:)

      CALL qtf_omega_list(omegas_main, max_dw_in, 0D0, tol, .TRUE., .FALSE., need, reuse_index, is_diff, is_sum)
   END SUBROUTINE qtf_radiation_list

   SUBROUTINE qtf_omega_list(omegas_main, max_dw_in, sum_top, tol, want_diff, want_sum, need, reuse_index, is_diff, is_sum)
      ! The distinct frequencies Omega at which psi is needed, ascending: with want_diff the differences
      ! Omega = omega_i - omega_j > 0 of the grid omegas_main with Omega <= max_dw_in (all of them for max_dw_in = 0), as in
      ! qtf_finish; with want_sum the sums omega_i + omega_j (i <= j) with Omega <= sum_top (all of them for sum_top = 0).
      ! Every one is within tol of its need(m): the values within tol/2 of the smallest of a cluster are one, and
      ! reuse_index(m) is the main-grid index whose frequency is need(m) within tol/2, else 0; need(m) is then that frequency
      ! exactly, so that the radiation solve of the main grid finds it. is_diff(m) and is_sum(m) say which kernels use
      ! need(m): a frequency that is both a difference and a sum (0.6 from 0.9 - 0.3 and 0.3 + 0.3) has one solve.
      ! qtf_init passes qtf_psi_tolerance, which merges rounding only. With want_sum = .FALSE. the list is the one of the
      ! Phases 6 and 7a to the bit.
      REAL(8), INTENT(IN) :: omegas_main(:), max_dw_in, sum_top, tol
      LOGICAL, INTENT(IN) :: want_diff, want_sum
      REAL(8), ALLOCATABLE, INTENT(OUT) :: need(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: reuse_index(:)
      LOGICAL, ALLOCATABLE, INTENT(OUT) :: is_diff(:), is_sum(:)
      REAL(8), ALLOCATABLE :: d(:), first(:), raw(:)
      REAL(8) :: dw
      INTEGER :: n, i, j, m, nd, ns, lo, hi, mid

      n = SIZE(omegas_main)
      ALLOCATE (d(n*n + n*(n + 1)/2))
      nd = 0
      IF (want_diff) THEN
         DO i = 1, n
            DO j = 1, n
               dw = omegas_main(i) - omegas_main(j)
               IF (dw <= tol) CYCLE
               IF (max_dw_in > 0D0 .AND. dw > max_dw_in*(1D0 + 1D-9)) CYCLE
               nd = nd + 1
               d(nd) = dw
            END DO
         END DO
      END IF
      ns = nd
      IF (want_sum) THEN
         DO i = 1, n
            DO j = i, n
               dw = omegas_main(i) + omegas_main(j)
               IF (sum_top > 0D0 .AND. dw > sum_top*(1D0 + 1D-9)) CYCLE
               nd = nd + 1
               d(nd) = dw
            END DO
         END DO
      END IF
      ! raw keeps the values in the order of d, the differences first (d(1:ns)) and the sums after them; the flags of a
      ! cluster come from the values that fall in it
      ALLOCATE (raw(nd))
      raw = d(1:nd)
      CALL heap_sort(d(1:nd))
      ALLOCATE (need(nd))
      m = 0
      DO i = 1, nd
         IF (m > 0) THEN
            IF (d(i) - need(m) <= 0.5D0*tol) CYCLE
         END IF
         m = m + 1
         need(m) = d(i)
      END DO
      need = need(1:m)
      first = need
      ALLOCATE (is_diff(m), is_sum(m))
      is_diff = .FALSE.
      is_sum = .FALSE.
      DO i = 1, nd
         ! The cluster of a value is the last one whose first value is at or below it (bisection)
         lo = 1
         hi = m
         DO WHILE (lo < hi)
            mid = (lo + hi + 1)/2
            IF (first(mid) <= raw(i)) THEN
               lo = mid
            ELSE
               hi = mid - 1
            END IF
         END DO
         IF (i <= ns) THEN
            is_diff(lo) = .TRUE.
         ELSE
            is_sum(lo) = .TRUE.
         END IF
      END DO
      ALLOCATE (reuse_index(m))
      reuse_index = 0
      DO i = 1, m
         DO j = 1, n
            IF (ABS(omegas_main(j) - need(i)) <= 0.5D0*tol) THEN
               reuse_index(i) = j
               need(i) = omegas_main(j)
               EXIT
            END IF
         END DO
      END DO
   END SUBROUTINE qtf_omega_list

   PURE SUBROUTINE heap_sort(a)
      ! Ascending, in place
      REAL(8), INTENT(INOUT) :: a(:)
      INTEGER :: n, i
      REAL(8) :: t

      n = SIZE(a)
      DO i = n/2, 1, -1
         CALL sift(i, n)
      END DO
      DO i = n, 2, -1
         t = a(1)
         a(1) = a(i)
         a(i) = t
         CALL sift(1, i - 1)
      END DO

   CONTAINS

      PURE SUBROUTINE sift(first, last)
         INTEGER, INTENT(IN) :: first, last
         INTEGER :: root, child
         REAL(8) :: v

         root = first
         v = a(root)
         DO
            child = 2*root
            IF (child > last) EXIT
            IF (child < last) THEN
               IF (a(child + 1) > a(child)) child = child + 1
            END IF
            IF (a(child) <= v) EXIT
            a(root) = a(child)
            root = child
         END DO
         a(root) = v
      END SUBROUTINE sift

   END SUBROUTINE heap_sort

   FUNCTION qtf_radiation_omegas() RESULT(w)
      ! The off-grid entries of the union list: each needs a radiation-only solve. These are the difference frequencies
      ! (qtf = 1 or 3) and, with the sum on, the sum frequencies Omega = omega_i + omega_j. RadiationOnlySolves uses the
      ! list with qtf_direct = 0. Empty unless qtf_terms = body or full.
      REAL(8), ALLOCATABLE :: w(:)

      IF (ALLOCATED(psi_need)) THEN
         w = PACK(psi_need, psi_reuse == 0)
      ELSE
         ALLOCATE (w(0))
      END IF
   END FUNCTION qtf_radiation_omegas

   FUNCTION qtf_direct_omegas() RESULT(w)
      ! Every difference frequency Omega > 0 within qtf_max_dw, those on the main grid too: the direct solve of Phase 7a runs
      ! at each of them. Empty unless qtf_terms = body.
      REAL(8), ALLOCATABLE :: w(:)

      IF (ALLOCATED(psi_need)) THEN
         w = psi_need
      ELSE
         ALLOCATE (w(0))
      END IF
   END FUNCTION qtf_direct_omegas

   LOGICAL FUNCTION qtf_omega_on_grid(omega)
      ! The difference frequency omega is a main-grid frequency, whose psi the radiation solve of the main grid kept
      REAL(8), INTENT(IN) :: omega
      INTEGER :: m

      qtf_omega_on_grid = .FALSE.
      m = qtf_psi_index(omega)
      IF (m > 0) qtf_omega_on_grid = psi_reuse(m) > 0
   END FUNCTION qtf_omega_on_grid

   INTEGER FUNCTION qtf_psi_count(reused)
      ! Number of difference frequencies, or with reused = .TRUE. of those on the main grid
      LOGICAL, INTENT(IN), OPTIONAL :: reused

      qtf_psi_count = 0
      IF (.NOT. ALLOCATED(psi_need)) RETURN
      qtf_psi_count = SIZE(psi_need)
      IF (PRESENT(reused)) THEN
         IF (reused) qtf_psi_count = COUNT(psi_reuse > 0)
      END IF
   END FUNCTION qtf_psi_count

   INTEGER FUNCTION qtf_psi_index(omega)
      ! Index of the difference frequency nearest to omega, if within psi_tol, or 0. psi_need is ascending: bisection to the
      ! first need >= omega - psi_tol, then the nearest up to omega + psi_tol
      REAL(8), INTENT(IN) :: omega
      INTEGER :: m, lo, hi, mid

      qtf_psi_index = 0
      IF (.NOT. ALLOCATED(psi_need)) RETURN
      lo = 1
      hi = SIZE(psi_need) + 1
      DO WHILE (lo < hi)
         mid = (lo + hi)/2
         IF (psi_need(mid) < omega - psi_tol) THEN
            lo = mid + 1
         ELSE
            hi = mid
         END IF
      END DO
      DO m = lo, SIZE(psi_need)
         IF (psi_need(m) > omega + psi_tol) EXIT
         IF (qtf_psi_index == 0) THEN
            qtf_psi_index = m
         ELSE IF (ABS(psi_need(m) - omega) < ABS(psi_need(qtf_psi_index) - omega)) THEN
            qtf_psi_index = m
         END IF
      END DO
   END FUNCTION qtf_psi_index

   PURE FUNCTION psi_bad_files(is_sum, is_diff, terms_diff, terms_sum) RESULT(names)
      ! The output files a non-finite psi feeds, for the message of qtf_keep_psi: QTF_Indirect.11d if the psi is of a
      ! difference frequency (is_diff, or neither flag), QTF_Indirect.11s if of a sum (is_sum), each with its terms file
      ! (QTF_PotentialTerms.txt, QTF_PotentialTermsSum.txt) only when this run writes it (terms_diff, terms_sum: write_terms >= 1)
      LOGICAL, INTENT(IN) :: is_sum, is_diff, terms_diff, terms_sum
      CHARACTER(LEN=120) :: names
      CHARACTER(LEN=60) :: first, second
      LOGICAL :: has_diff, has_sum

      has_sum = is_sum
      has_diff = is_diff .OR. .NOT. is_sum
      first = ''
      second = ''
      IF (has_diff) THEN
         first = 'QTF_Indirect.11d'
         IF (terms_diff) first = TRIM(first)//', QTF_PotentialTerms.txt'
      END IF
      IF (has_sum) THEN
         second = 'QTF_Indirect.11s'
         IF (terms_sum) second = TRIM(second)//', QTF_PotentialTermsSum.txt'
      END IF
      IF (has_diff .AND. has_sum) THEN
         names = TRIM(first)//', '//TRIM(second)
      ELSE
         names = TRIM(first)//TRIM(second)
      END IF
      ! Join the last two names with "and", as the single-route messages do
      names = and_last(names)
   END FUNCTION psi_bad_files

   PURE FUNCTION and_last(list) RESULT(text)
      ! "a, b, c" -> "a, b and c"; a list of one name is unchanged
      CHARACTER(LEN=*), INTENT(IN) :: list
      CHARACTER(LEN=LEN(list)+4) :: text
      INTEGER :: k

      k = INDEX(TRIM(list), ',', BACK=.TRUE.)
      IF (k == 0) THEN
         text = list
      ELSE
         text = list(1:k-1)//' and'//list(k+1:)
      END IF
   END FUNCTION and_last

   SUBROUTINE qtf_keep_psi(omega, psi, ok)
      ! Keeps psi(n_panels, 6*nbody) of the difference frequency omega. An omega that is not a difference frequency, or a psi
      ! of another size, stops the run; with ok present, psi is not kept and ok is false instead. A psi with NaN or Inf is
      ! kept, logged with omega and the first mode and body concerned, and counted in qtf_nonfinite_psi().
      REAL(8), INTENT(IN) :: omega
      COMPLEX(8), INTENT(IN) :: psi(:, :)
      LOGICAL, INTENT(OUT), OPTIONAL :: ok
      INTEGER :: m, nbad, k
      CHARACTER(LEN=120) :: bad_files   ! the files this psi feeds: those of the difference and of the sum

      IF (PRESENT(ok)) ok = .FALSE.
      m = qtf_psi_index(omega)
      IF (m == 0) THEN
         IF (PRESENT(ok)) RETURN
         WRITE (*, '(a,ES14.6,a)') ' QTF: psi at omega =', omega, ' rad/s, which is not a difference frequency of the grid.'
         ERROR STOP 1
      END IF
      IF (SIZE(psi, 1) /= SIZE(psi_store, 1) .OR. SIZE(psi, 2) /= SIZE(psi_store, 2)) THEN
         IF (PRESENT(ok)) RETURN
         WRITE (*, '(a)') ' QTF: psi does not have the size the store was set up for.'
         ERROR STOP 1
      END IF
      nbad = COUNT(.NOT. (IEEE_IS_FINITE(REAL(psi, 8)) .AND. IEEE_IS_FINITE(AIMAG(psi))))
      IF (nbad > 0) THEN
         n_bad_psi = n_bad_psi + 1
         DO k = 1, SIZE(psi, 2)
            IF (.NOT. ALL(IEEE_IS_FINITE(REAL(psi(:, k), 8)) .AND. IEEE_IS_FINITE(AIMAG(psi(:, k))))) EXIT
         END DO
         WRITE (*, '(a,ES14.6,a,I0,a,I0,a,I0,a,a)') ' QTF error: psi at Omega =', omega, ' rad/s has ', nbad, &
            ' values with NaN or Inf, the first in mode ', k - 6*((k - 1)/6), ' of body ', (k - 1)/6 + 1, ' (', &
            TRIM(MERGE('radiation solve of the main grid', 'radiation-only solve            ', psi_reuse(m) > 0))//').'
         bad_files = psi_bad_files(psi_is_sum(m), psi_is_diff(m), pot_terms, pot_sum_terms)
         ! "is" for one file, "are" for a list (psi_bad_files joins a list with a comma or "and")
         IF (INDEX(bad_files, ',') > 0 .OR. INDEX(bad_files, ' and ') > 0) THEN
            WRITE (*, '(a)') ' '//TRIM(bad_files)//' are not valid at that Omega; the run will end with exit code 1.'
         ELSE
            WRITE (*, '(a)') ' '//TRIM(bad_files)//' is not valid at that Omega; the run will end with exit code 1.'
         END IF
      END IF
      psi_store(:, :, m) = psi
      psi_kept(m) = .TRUE.
      IF (PRESENT(ok)) ok = .TRUE.
   END SUBROUTINE qtf_keep_psi

   INTEGER FUNCTION qtf_nonfinite_psi()
      ! Number of difference frequencies whose psi had NaN or Inf in this run
      qtf_nonfinite_psi = n_bad_psi
   END FUNCTION qtf_nonfinite_psi

   INTEGER FUNCTION qtf_nonfinite_pairs()
      ! Number of frequency and heading pairs of QTF_Indirect.11d and .11s whose T_pot or T+_pot had NaN or Inf in this run
      qtf_nonfinite_pairs = n_bad_pairs
   END FUNCTION qtf_nonfinite_pairs

   SUBROUTINE qtf_direct_pairs(omega, pairs)
      ! The pairs of the difference QTF whose |omega_i - omega_j| is the difference frequency omega (qtf_psi_index), with
      ! omega_i /= omega_j and within qtf_max_dw, as potential_part takes them: pairs(4, npair) = frequency and heading of
      ! state 1, then of state 2, in the store. Empty unless qtf_direct = 1, or before the frequencies are kept.
      REAL(8), INTENT(IN) :: omega
      INTEGER, ALLOCATABLE, INTENT(OUT) :: pairs(:, :)
      INTEGER, ALLOCATABLE :: tmp(:, :)
      INTEGER :: m, i, j, b1, b2, np

      ALLOCATE (pairs(4, 0))
      IF (.NOT. direct_active .OR. .NOT. diff_on) RETURN
      m = qtf_psi_index(omega)
      IF (m == 0) RETURN
      ALLOCATE (tmp(4, nkept*nkept*heading_pairs(nbeta)))
      np = 0
      DO b1 = 1, nbeta
         DO b2 = 1, nbeta
            IF (all_headings_pairs == 0 .AND. b1 /= b2) CYCLE
            DO i = 1, nkept
               DO j = 1, nkept
                  IF (difference_pair_index(i, j) /= m) CYCLE
                  np = np + 1
                  tmp(:, np) = [i, b1, j, b2]
               END DO
            END DO
         END DO
      END DO
      pairs = tmp(:, 1:np)
   END SUBROUTINE qtf_direct_pairs

   INTEGER FUNCTION difference_pair_index(i, j)
      ! The pair rule of the difference QTF, shared by potential_part (T_pot) and qtf_direct_pairs (T_dir) so that both
      ! cover the same pairs: 0 when kept frequencies i and j are not a pair (omega_i = omega_j, or |omega_i - omega_j|
      ! beyond qtf_max_dw), else the index of |omega_i - omega_j| in psi_need (qtf_psi_index), or -1 if it is not there.
      ! Heading 1 stands for all: qtf_keep stores every heading of a frequency from the same solve, so store(i, b)%omega
      ! is the same for every b (the grid is shared by construction).
      INTEGER, INTENT(IN) :: i, j
      REAL(8) :: dw

      difference_pair_index = 0
      dw = store(i, 1)%omega - store(j, 1)%omega
      IF (dw == 0D0) RETURN
      IF (max_dw > 0D0 .AND. ABS(dw) > max_dw*(1D0 + 1D-9)) RETURN
      difference_pair_index = qtf_psi_index(ABS(dw))
      IF (difference_pair_index == 0) difference_pair_index = -1
   END FUNCTION difference_pair_index

   LOGICAL FUNCTION qtf_diff_omega(omega)
      ! The frequency omega of the list is a difference frequency that the difference kernel of the direct solve uses. False
      ! for a sum-only frequency (Phase 7c-2), and when qtf_init has not built the list.
      REAL(8), INTENT(IN) :: omega
      INTEGER :: m

      qtf_diff_omega = .FALSE.
      IF (.NOT. ALLOCATED(psi_is_diff)) RETURN
      m = qtf_psi_index(omega)
      IF (m > 0) qtf_diff_omega = diff_on .AND. psi_is_diff(m)
   END FUNCTION qtf_diff_omega

   LOGICAL FUNCTION qtf_sum_omega(omega)
      ! The frequency omega of the list is a sum frequency that the sum kernel of the direct solve uses (Phase 7c-2). False
      ! when the sum is off, and when qtf_init has not built the list.
      REAL(8), INTENT(IN) :: omega
      INTEGER :: m

      qtf_sum_omega = .FALSE.
      IF (.NOT. ALLOCATED(psi_is_sum)) RETURN
      m = qtf_psi_index(omega)
      IF (m > 0) qtf_sum_omega = sum_on .AND. psi_is_sum(m)
   END FUNCTION qtf_sum_omega

   SUBROUTINE qtf_direct_sum_pairs(omega, pairs)
      ! The canonical half of the sum pairs whose omega_i + omega_j is the sum frequency omega (sum_pair_index), as
      ! pairs(4, npair) = frequency and heading of state 1, then of state 2, in the store: i < j, or i = j with b1 <= b2
      ! (ruling R6). T+ is symmetric (pair_plus = (f1 g2 + f2 g1)/4 in every term, and the sum Omega and psi do not depend
      ! on the order), so the other half is the same pair and is not solved again. The diagonal i = j, b1 = b2 is a pair
      ! of its own. Empty unless the sum is on (with or without qtf_direct, Phase 6b), or before the frequencies are kept.
      REAL(8), INTENT(IN) :: omega
      INTEGER, ALLOCATABLE, INTENT(OUT) :: pairs(:, :)
      INTEGER, ALLOCATABLE :: tmp(:, :)
      INTEGER :: m, i, j, b1, b2, np

      ALLOCATE (pairs(4, 0))
      IF (.NOT. sum_on) RETURN
      m = qtf_psi_index(omega)
      IF (m == 0) RETURN
      ALLOCATE (tmp(4, nkept*nkept*heading_pairs(nbeta)))
      np = 0
      DO b1 = 1, nbeta
         DO b2 = 1, nbeta
            IF (all_headings_pairs == 0 .AND. b1 /= b2) CYCLE
            DO i = 1, nkept
               DO j = i, nkept
                  IF (i == j .AND. b1 > b2) CYCLE
                  IF (sum_pair_index(i, j) /= m) CYCLE
                  np = np + 1
                  tmp(:, np) = [i, b1, j, b2]
               END DO
            END DO
         END DO
      END DO
      pairs = tmp(:, 1:np)
   END SUBROUTINE qtf_direct_sum_pairs

   SUBROUTINE qtf_direct_sum_load(i, b1, j, b2, t)
      ! T+_dir(6 nbody) of the stored pair (frequency i, heading b1; frequency j, heading b2) as tdir_sum holds it, for the
      ! tests of the slots (R6); zero before the direct solve at omega_i + omega_j, and when the sum is off
      INTEGER, INTENT(IN) :: i, b1, j, b2
      COMPLEX(8), INTENT(OUT) :: t(:)

      t = (0D0, 0D0)
      IF (ALLOCATED(tdir_sum)) t = tdir_sum(:, i, j, heading_pair(b1, b2))
   END SUBROUTINE qtf_direct_sum_load

   SUBROUTINE qtf_indirect_sum_load(i, b1, j, b2, t)
      ! T+_pot(6 nbody) of the stored pair (frequency i, heading b1; frequency j, heading b2) as tpot_sum holds it, for the
      ! tests of QTF_Indirect.11s; zero before qtf_finish computes it, and when the sum is off
      INTEGER, INTENT(IN) :: i, b1, j, b2
      COMPLEX(8), INTENT(OUT) :: t(:)

      t = (0D0, 0D0)
      IF (ALLOCATED(tpot_sum)) t = tpot_sum(:, i, j, heading_pair(b1, b2))
   END SUBROUTINE qtf_indirect_sum_load

   SUBROUTINE qtf_direct_diff_load(i, b1, j, b2, t, p_i2, p_d)
      ! T_dir(6 nbody) of the stored difference pair as tdir holds it, and its parts from phi_I2 and phi_D (write_terms
      ! >= 1, zero otherwise), for the test that qtf = 3 leaves the difference results of qtf = 1 unchanged
      INTEGER, INTENT(IN) :: i, b1, j, b2
      COMPLEX(8), INTENT(OUT) :: t(:), p_i2(:), p_d(:)

      t = (0D0, 0D0)
      p_i2 = (0D0, 0D0)
      p_d = (0D0, 0D0)
      IF (ALLOCATED(tdir)) t = tdir(:, i, j, heading_pair(b1, b2))
      IF (ALLOCATED(tdir_i2)) p_i2 = tdir_i2(:, i, j, heading_pair(b1, b2))
      IF (ALLOCATED(tdir_d)) p_d = tdir_d(:, i, j, heading_pair(b1, b2))
   END SUBROUTINE qtf_direct_diff_load

   SUBROUTINE free_sum_stores()
      IF (ALLOCATED(tdir_sum)) DEALLOCATE (tdir_sum)
      IF (ALLOCATED(tpot_sum)) DEALLOCATE (tpot_sum)
      IF (ALLOCATED(tdir_sum_i2)) DEALLOCATE (tdir_sum_i2)
      IF (ALLOCATED(tdir_sum_d)) DEALLOCATE (tdir_sum_d)
      IF (ALLOCATED(tdir_sum_fs)) DEALLOCATE (tdir_sum_fs)
   END SUBROUTINE free_sum_stores

   SUBROUTINE qtf_direct_at_omega(mesh, geom, omega, npair, seconds, fs, nsum)
      ! The direct method at the frequency omega of the list (qtf_direct = 1): T_dir of every difference pair with that
      ! Omega (qtf_direct_pairs) and, with the sum on (Phase 7c-2), T+_dir of the canonical half of its sum pairs
      ! (qtf_direct_sum_pairs), with one factorisation: two direct_at_omega calls, the difference columns as without the
      ! sum, then the sum columns. The difference results go into tdir, each sum result into both slots of tdir_sum
      ! (R6). Call after the frequency loop, when every frequency is kept, with the solver at omega after its radiation solve
      ! (RadiationOnlySolves). Pairs with NaN or Inf are counted in qtf_nonfinite_pairs() and logged. npair, nsum and
      ! seconds (optional): the number of difference and of sum pairs solved, and the times of direct_at_omega (4: forcing,
      ! B2, solve and load, field evaluation), for the log. A frequency with neither kind of pair is only marked as done.
      ! fs (optional, qtf_terms = full): the free-surface data of direct_fs_setup at omega; its pointers are set here to
      ! the free-surface store (fields and incident fields at the inner points, partial waves) and the sliver chords, so
      ! that T_dir has the free-surface forcing. Stops if a free-surface field was not stored.
      ! With qtf_field = 1 the field points get phi_D of the body forcing only (direct_at_omega solves the free-surface
      ! part apart when both are present), as the header of QTF_Field.txt says; only the difference columns are kept there
      ! (the sum rows of QTF_Field.txt come from the second call below).
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      REAL(8), INTENT(IN) :: omega
      INTEGER, INTENT(OUT), OPTIONAL :: npair
      REAL(8), INTENT(OUT), OPTIONAL :: seconds(4)
      TYPE(direct_fs_data), INTENT(INOUT), OPTIONAL :: fs
      INTEGER, INTENT(OUT), OPTIONAL :: nsum
      REAL(8) :: sec(4), sec2(4)
      INTEGER, ALLOCATABLE :: pd(:, :), ps(:, :), kernels(:)
      COMPLEX(8), ALLOCATABLE :: t(:, :), parts(:, :, :), phi_field(:, :)
      REAL(8), ALLOCATABLE :: xfield(:, :)
      INTEGER :: m, p, nbad, h, hm, nd, ns

      IF (PRESENT(npair)) npair = 0
      IF (PRESENT(nsum)) nsum = 0
      IF (PRESENT(seconds)) seconds = 0D0
      IF (.NOT. direct_active) RETURN
      m = qtf_psi_index(omega)
      IF (m == 0) THEN
         WRITE (*, '(a,ES14.6,a)') ' QTF direct: omega =', omega, ' rad/s is not a difference or sum frequency of the grid.'
         ERROR STOP 1
      END IF
      IF (.NOT. (qtf_diff_omega(omega) .OR. qtf_sum_omega(omega))) THEN
         ! Neither kernel has a pair here (a sum frequency of a run without the sum cannot occur: the list is built
         ! with the kernels that are on), so there is nothing to solve
         tdir_kept(m) = .TRUE.
         RETURN
      END IF
      IF (nkept /= SIZE(store, 1)) THEN
         WRITE (*, '(a,I0,a,I0,a)') ' QTF direct: ', nkept, ' of ', SIZE(store, 1), &
            ' frequencies are kept; the direct solve runs after the frequency loop.'
         ERROR STOP 1
      END IF
      IF (qtf_diff_omega(omega)) THEN
         CALL qtf_direct_pairs(omega, pd)
      ELSE
         ALLOCATE (pd(4, 0))
      END IF
      IF (qtf_sum_omega(omega)) THEN
         CALL qtf_direct_sum_pairs(omega, ps)
      ELSE
         ALLOCATE (ps(4, 0))
      END IF
      nd = SIZE(pd, 2)
      ns = SIZE(ps, 2)
      ! Two calls of direct_at_omega with one factorisation: first the difference columns,
      ! exactly as without the sum (the call of 3f00f4d, so the same column blocks, ZGETRS and ZGEMM sizes and bits; MKL
      ! does not promise column p to be the same for another number of right sides), then the sum columns with
      ! factored = .TRUE. (the factors of the first call; with no difference pair the second call factors)
      sec = 0D0
      sec2 = 0D0
      IF (PRESENT(fs)) CALL direct_fs_store(fs)
      nbad = 0
      IF (nd > 0) THEN
         ALLOCATE (t(6*mesh%nbody, nd))
         IF (PRESENT(fs)) THEN
            ! The third slice keeps the free-surface part apart, for QTF_DirectTerms.txt
            ALLOCATE (parts(6*mesh%nbody, nd, MERGE(3, 2, direct_terms)))
         ELSE
            ALLOCATE (parts(6*mesh%nbody, nd, 2))
         END IF
         ! An absent fs is passed on as absent: the direct solve of Phase 7a
         IF (qfield_direct_active()) THEN
            ! qtf_field = 1 with qtf_direct = 1: phi_D of every pair at the field points, for QTF_Field.txt (QTFField)
            CALL qfield_direct_points(xfield)
            ALLOCATE (phi_field(SIZE(xfield, 2), nd))
            CALL direct_at_omega(mesh, geom, omega, store, store_wl, pd, t, omega_tol=psi_tol, seconds=sec, parts=parts, &
                                 fs=fs, field_points=xfield, phi_field=phi_field)
            DO p = 1, nd
               CALL qfield_keep_direct(pd(1, p), pd(2, p), pd(3, p), pd(4, p), phi_field(:, p))
            END DO
            DEALLOCATE (phi_field, xfield)
         ELSE
            CALL direct_at_omega(mesh, geom, omega, store, store_wl, pd, t, omega_tol=psi_tol, seconds=sec, parts=parts, &
                                 fs=fs)
         END IF
         DO p = 1, nd
            h = heading_pair(pd(2, p), pd(4, p))
            tdir(:, pd(1, p), pd(3, p), h) = t(:, p)
            IF (direct_terms) THEN
               tdir_i2(:, pd(1, p), pd(3, p), h) = parts(:, p, 1)
               tdir_d(:, pd(1, p), pd(3, p), h) = parts(:, p, 2)
               IF (SIZE(parts, 3) == 3) tdir_fs(:, pd(1, p), pd(3, p), h) = parts(:, p, 3)
            END IF
            IF (.NOT. ALL(IEEE_IS_FINITE(REAL(t(:, p), 8)) .AND. IEEE_IS_FINITE(AIMAG(t(:, p))))) nbad = nbad + 1
         END DO
         DEALLOCATE (t, parts)
      END IF
      IF (ns > 0) THEN
         ALLOCATE (t(6*mesh%nbody, ns), kernels(ns))
         kernels = KERNEL_SUM
         IF (PRESENT(fs)) THEN
            ALLOCATE (parts(6*mesh%nbody, ns, MERGE(3, 2, sum_terms)))
         ELSE
            ALLOCATE (parts(6*mesh%nbody, ns, 2))
         END IF
         IF (qfield_direct_active()) THEN
            ! qtf_field = 1 with qtf_direct = 1: phi_D of every sum pair at the field points for the sum rows of
            ! QTF_Field.txt (QTFField), body-forced as for the difference columns (R4). Not conjugated: Omega+ > 0.
            CALL qfield_direct_points(xfield)
            ALLOCATE (phi_field(SIZE(xfield, 2), ns))
            CALL direct_at_omega(mesh, geom, omega, store, store_wl, ps, t, omega_tol=psi_tol, seconds=sec2, parts=parts, &
                                 fs=fs, kernels=kernels, factored=nd > 0, field_points=xfield, phi_field=phi_field)
            DO p = 1, ns
               ! One call fills both slots of the symmetric pair (R6)
               CALL qfield_keep_direct(ps(1, p), ps(2, p), ps(3, p), ps(4, p), phi_field(:, p), KERNEL_SUM)
            END DO
            DEALLOCATE (phi_field, xfield)
         ELSE
            CALL direct_at_omega(mesh, geom, omega, store, store_wl, ps, t, omega_tol=psi_tol, seconds=sec2, parts=parts, &
                                 fs=fs, kernels=kernels, factored=nd > 0)
         END IF
         sec = sec + sec2
         DO p = 1, ns
            ! Both slots of the symmetric pair (R6): (i, b1, j, b2) and its mirror (j, b2, i, b1), the same values
            h = heading_pair(ps(2, p), ps(4, p))
            hm = heading_pair(ps(4, p), ps(2, p))
            tdir_sum(:, ps(1, p), ps(3, p), h) = t(:, p)
            tdir_sum(:, ps(3, p), ps(1, p), hm) = t(:, p)
            IF (sum_terms) THEN
               tdir_sum_i2(:, ps(1, p), ps(3, p), h) = parts(:, p, 1)
               tdir_sum_i2(:, ps(3, p), ps(1, p), hm) = parts(:, p, 1)
               tdir_sum_d(:, ps(1, p), ps(3, p), h) = parts(:, p, 2)
               tdir_sum_d(:, ps(3, p), ps(1, p), hm) = parts(:, p, 2)
               IF (SIZE(parts, 3) == 3 .AND. ALLOCATED(tdir_sum_fs)) THEN
                  tdir_sum_fs(:, ps(1, p), ps(3, p), h) = parts(:, p, 3)
                  tdir_sum_fs(:, ps(3, p), ps(1, p), hm) = parts(:, p, 3)
               END IF
            END IF
            IF (.NOT. ALL(IEEE_IS_FINITE(REAL(t(:, p), 8)) .AND. IEEE_IS_FINITE(AIMAG(t(:, p))))) nbad = nbad + 1
         END DO
      END IF
      ! Phase 6c: psi_j at the field points at this Omega, for the radiated wave of the motion in QTF_Field.txt (R7); after
      ! the direct solves, whose results it does not touch. Its time goes with the field evaluation (sec(4)).
      IF (mot_active .AND. qfield_radiated_active()) CALL keep_psi_points()
      IF (PRESENT(npair)) npair = nd
      IF (PRESENT(nsum)) nsum = ns
      IF (PRESENT(seconds)) seconds = sec
      tdir_kept(m) = .TRUE.
      IF (nbad > 0) THEN
         n_bad_pairs = n_bad_pairs + nbad
         WRITE (*, '(a,I0,a,ES14.6,a)') ' QTF error: T_dir has NaN or Inf at ', nbad, ' frequency and heading pairs with'// &
            ' Omega =', psi_need(m), ' rad/s.'
      END IF

   CONTAINS

      SUBROUTINE direct_fs_store(d)
         ! The pointers of d to the free-surface store, indexed like store; the incident fields at the inner points are
         ! built once, at the first Omega (fs_incident_all)
         TYPE(direct_fs_data), INTENT(INOUT) :: d

         IF (.NOT. fs_active) THEN
            WRITE (*, '(a)') ' QTF direct: the free-surface forcing needs the free-surface store (qtf_terms = full).'
            ERROR STOP 1
         END IF
         IF (.NOT. ALL(fs_kept)) THEN
            WRITE (*, '(a)') ' QTF direct: the free-surface fields of a kept frequency were not stored.'
            ERROR STOP 1
         END IF
         IF (.NOT. ALLOCATED(fs_inc_store)) CALL fs_incident_all(fs_inc_store)
         d%f => fs_f
         d%inc => fs_inc_store
         d%pw => fs_pw
         d%np = fs_store_np
      END SUBROUTINE direct_fs_store

      SUBROUTINE keep_psi_points()
         ! psi of this Omega (kept by qtf_keep_psi before the direct solve) at the field points, into mot_psi_pts(:, :, m)
         REAL(8), ALLOCATABLE :: xp(:, :)
         INTEGER(8) :: c0, c1, rate

         IF (.NOT. psi_kept(m)) THEN
            WRITE (*, '(a,ES14.6,a)') ' QTF motion: psi at Omega =', omega, ' rad/s is not kept before the direct solve.'
            ERROR STOP 1
         END IF
         CALL SYSTEM_CLOCK(c0, rate)
         CALL qfield_direct_points(xp)
         IF (.NOT. ALLOCATED(mot_psi_pts)) THEN
            ALLOCATE (mot_psi_pts(SIZE(xp, 2), SIZE(psi_store, 2), SIZE(psi_need)), mot_psi_kept(SIZE(psi_need)))
            mot_psi_pts = (0D0, 0D0)
            mot_psi_kept = .FALSE.
         END IF
         CALL psi_at_points(mesh, psi_store(:, :, m), omega, xp, mot_psi_pts(:, :, m))
         mot_psi_kept(m) = .TRUE.
         CALL SYSTEM_CLOCK(c1)
         sec(4) = sec(4) + REAL(c1 - c0, 8)/REAL(rate, 8)
      END SUBROUTINE keep_psi_points

   END SUBROUTINE qtf_direct_at_omega

   SUBROUTINE psi_at_points(mesh, psi, omega, points, psi_pts)
      ! The radiation potentials psi_j (6 nbody columns) at the field points, e^{-i omega t}, by Green's formula with their
      ! body values and dpsi_j/dn = n_j: direct_field_at_points of QTFDirect on the 6 nbody columns (ruling R7 of Phase 6c).
      !   psi      (n, 6 nbody)    radiation_psi at omega on the full panels, as psi_store holds it
      !   points   (3, npt)        points in the fluid off the hull, full coordinates
      !   psi_pts  (npt, 6 nbody)  out
      ! The solver globals must be those of omega (> 0), as for direct_at_omega: W1 is checked.
      !
      ! Signs and images. psi_k is the unit-velocity potential of mode k with dpsi_k/dn = n_k on the body of mode k and 0
      ! on the others, n out of the fluid (into the body): the normal of radiation_psi, of radiation_coefficients (QTFMotion
      ! header) and of Green's formula in direct_field_at_points, 4 pi phi(x) = int_S [G dphi/dn - phi dG/dn] dS
      ! (QTFDirectField header). So the panel forcing vf of column k is n_k, in the layout of pair_forcing:
      ! vf(base(f), k, image(f)) on the full panel f. mesh%nj of an image panel holds the image signs RX, RY already (as
      ! radiation_coefficients uses it with no factor), so vf takes it as it is; the images of the Green function are those
      ! of direct_field_at_points (the panels in place, x mirrored). The radiation identity of 7c-1
      ! (direct_field_radiation_identity: psi_3 through direct_field_at_points equals fields_green of the state with
      ! phi_p = psi_3, dphi_p/dn = n_3) pins this sign. psi has no B2 part: a, c, e and cl are zero, and an empty waterline
      ! (a surface_geom with nwl = 0) skips the line integrals of b2_stokes_at_point, which would add zero.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      COMPLEX(8), INTENT(IN) :: psi(:, :)
      REAL(8), INTENT(IN) :: omega, points(:, :)
      COMPLEX(8), INTENT(OUT) :: psi_pts(:, :)
      TYPE(surface_geom) :: no_waterline
      COMPLEX(8), ALLOCATABLE :: vf(:, :, :), a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :)
      INTEGER :: nmode, k, kb, km, f

      nmode = 6*mesh%nbody
      IF (.NOT. (omega > 0D0) .OR. ABS(W1 - omega) > 1D-12*omega) THEN
         WRITE (*, '(a,ES14.6,a,ES14.6,a)') ' psi_at_points: Omega =', omega, ' rad/s, but the solver is at', W1, ' rad/s.'
         ERROR STOP 1
      END IF
      IF (SIZE(psi, 1) /= mesh%n .OR. SIZE(psi, 2) /= nmode .OR. SIZE(points, 1) /= 3 .OR. &
          SIZE(psi_pts, 1) /= SIZE(points, 2) .OR. SIZE(psi_pts, 2) /= nmode) THEN
         WRITE (*, '(a)') ' psi_at_points: psi, points or psi_pts does not have the size of the mesh and of each other.'
         ERROR STOP 1
      END IF
      ALLOCATE (vf(mesh%nbase, nmode, mesh%nsys), a(mesh%n, nmode), c(3, mesh%n, nmode), e(0, nmode), &
                cl(3, 4, mesh%n, nmode))
      vf = (0D0, 0D0)
      DO k = 1, nmode
         kb = (k - 1)/6 + 1
         km = k - 6*(kb - 1)
         DO f = 1, mesh%n
            IF (mesh%body(f) == kb) vf(mesh%base(f), k, mesh%image(f)) = CMPLX(mesh%nj(km, f), 0D0, KIND=8)
         END DO
      END DO
      a = (0D0, 0D0)
      c = (0D0, 0D0)
      cl = (0D0, 0D0)
      ! The multi-body arrays for several bodies, as direct_at_omega decides by default
      CALL direct_field_at_points(mesh, no_waterline, points, psi, vf, a, c, e, cl, mesh%nbody > 1, psi_pts)
   END SUBROUTINE psi_at_points

   SUBROUTINE fs_incident_all(inc)
      ! The incident fields of every kept state at the inner free-surface points, as fs_inner_fields adds them: phi_I,
      ! grad phi_I and phi_I,zz = k^2 phi_I (incident_at); inc(n_fs, nkept, nbeta)
      TYPE(fs_point_field), ALLOCATABLE, INTENT(OUT) :: inc(:, :, :)
      TYPE(fs_point_field) :: f
      INTEGER :: ip, is, ib

      ALLOCATE (inc(fs_store_pts%n, nkept, nbeta))
      DO ib = 1, nbeta
         DO is = 1, nkept
            DO ip = 1, fs_store_pts%n
               CALL incident_at(store(is, ib), fs_store_pts%x(:, ip), f%phi, f%grad)
               f%phizz = store(is, ib)%k**2*f%phi
               inc(ip, is, ib) = f
            END DO
         END DO
      END DO
   END SUBROUTINE fs_incident_all

   SUBROUTINE fs_boundary_all(mesh, geom, bnd, xe, wb)
      ! qtf_fs_inner = source (item 2 of the field-pass speed-ups): every kept state on the boundary of the inner free
      ! surface for the line integral (fs_wl_state), bnd(nkept, nbeta): first the waterline segments of geom
      ! (fs_waterline_state, with the waterline potential of the store, as the near-field drift takes it), then the points
      ! xe(2, ne) of the outer edge (fs_edge_points of the store's circle, M and polygon; fs_edge_state from the state's
      ! partial waves and incident wave). wb(nwl + ne): the segment lengths, then the edge weights.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      TYPE(fs_wl_state), ALLOCATABLE, INTENT(OUT) :: bnd(:, :)
      REAL(8), ALLOCATABLE, INTENT(OUT) :: xe(:, :), wb(:)
      TYPE(fs_wl_state) :: wl, ed
      REAL(8), ALLOCATABLE :: nhe(:, :), we(:)
      INTEGER :: is, ib, nw

      CALL fs_edge_points(fs_store_np, fs_store_centre, fs_store_radius, fs_store_mmax, xe, nhe, we)
      nw = geom%nwl
      wb = [geom%wl_len(1:nw), we]
      ALLOCATE (bnd(nkept, nbeta))
      DO ib = 1, nbeta
         DO is = 1, nkept
            CALL fs_waterline_state(mesh, geom, store(is, ib), store_wl(:, is, ib), wl)
            CALL fs_edge_state(fs_pw(is, ib), store(is, ib)%omega, store(is, ib)%k, store(is, ib)%beta, xe, nhe, ed)
            bnd(is, ib)%phi = [wl%phi, ed%phi]
            bnd(is, ib)%dn = [wl%dn, ed%dn]
            bnd(is, ib)%phi_i = [wl%phi_i, ed%phi_i]
            bnd(is, ib)%dn_i = [wl%dn_i, ed%dn_i]
         END DO
      END DO
   END SUBROUTINE fs_boundary_all

   SUBROUTINE fs_boundary_psi(nwl, pwl, pw_psi, xe, psi_b)
      ! psi_k on the boundary of fs_boundary_all: pwl(nwl, 6*nbody) at the waterline segments (psi_waterline, as the body
      ! terms take it), then on the outer edge points xe from the partial waves pw_psi (fs_edge_psi)
      INTEGER, INTENT(IN) :: nwl
      COMPLEX(8), INTENT(IN) :: pwl(:, :)
      TYPE(pw_set), INTENT(IN) :: pw_psi(:)
      REAL(8), INTENT(IN) :: xe(:, :)
      COMPLEX(8), INTENT(OUT) :: psi_b(:, :)

      psi_b(1:nwl, :) = pwl
      CALL fs_edge_psi(pw_psi, xe, psi_b(nwl + 1:, :))
   END SUBROUTINE fs_boundary_psi

   REAL(8) FUNCTION qtf_direct_missing()
      ! The smallest difference frequency of the grid (rad/s) at which the direct solve has not run, or 0 if there is none
      ! or qtf_direct is off (or qtf_init has not set up its store). qtf_finish stops on a non-zero value before it writes
      ! QTF_Direct.12d.
      INTEGER :: m

      qtf_direct_missing = 0D0
      IF (.NOT. direct_active) RETURN
      IF (.NOT. ALLOCATED(tdir_kept)) RETURN
      DO m = 1, SIZE(tdir_kept)
         IF (.NOT. tdir_kept(m)) THEN
            qtf_direct_missing = psi_need(m)
            RETURN
         END IF
      END DO
   END FUNCTION qtf_direct_missing

   SUBROUTINE qtf_keep(ib, s, phi_wl, oufr)
      ! Keeps the state of heading ib of the current frequency. Heading 1 starts a new frequency, so the driver must call
      ! this for every heading of a frequency, in order, and only for the frequencies of the QTF grid.
      INTEGER, INTENT(IN) :: ib
      TYPE(first_order_state), INTENT(IN) :: s
      COMPLEX(8), INTENT(IN) :: phi_wl(:)
      REAL(8), INTENT(IN) :: oufr
      IF (.NOT. active) RETURN
      IF (ib < 1 .OR. ib > nbeta .OR. (ib /= 1 .AND. nkept == 0)) THEN
         WRITE (*, '(a,I0)') ' QTF: qtf_keep called out of order, heading ', ib
         ERROR STOP 1
      END IF
      IF (ib == 1) THEN
         nkept = nkept + 1
         IF (nkept > SIZE(store, 1)) THEN
            WRITE (*, '(a)') ' QTF: more frequencies than the store was set up for.'
            ERROR STOP 1
         END IF
         store_oufr(nkept) = oufr
      END IF
      IF (SIZE(s%phi_p) /= SIZE(store(nkept, ib)%phi_p) .OR. SIZE(s%xi, 2) /= SIZE(store(nkept, ib)%xi, 2)) THEN
         WRITE (*, '(a)') ' QTF: the state does not have the size the store was set up for.'
         ERROR STOP 1
      END IF
      ASSOCIATE (c => store(nkept, ib))
         c%omega = s%omega
         c%k = s%k
         c%nu = s%nu
         c%h = s%h
         c%amp = s%amp
         c%beta = s%beta
         c%g = s%g
         c%rho = s%rho
         c%phi_i = s%phi_i
         c%phi_p = s%phi_p
         c%grad_i = s%grad_i
         c%grad_p = s%grad_p
         c%xi = s%xi
      END ASSOCIATE
      store_wl(:, nkept, ib) = phi_wl
   END SUBROUTINE qtf_keep

   INTEGER FUNCTION qtf_kept()
      ! Number of frequencies kept so far
      qtf_kept = nkept
   END FUNCTION qtf_kept

   SUBROUTINE qtf_set_fs_inner(mode)
      ! The route of the inner free-surface fields for qtf_stored_mb: qtf_fs_inner, 1 stencil, 2 source (grad psi counted)
      INTEGER, INTENT(IN) :: mode

      estimate_gpsi = mode == 2
   END SUBROUTINE qtf_set_fs_inner

   REAL(8) FUNCTION qtf_stored_mb(n_panels, nwl, nbody, nfreq, nb, npsi, npairs, ncs, nfs, npw, nsum_pairs, nmotion)
      ! Size of the store in MB for nfreq frequencies and nb headings: 8 complex values per panel, 6 per body, 1 per segment.
      ! With npsi, also psi at npsi difference frequencies: 6*nbody complex values per panel. With npairs heading pairs, also
      ! T_pot: 6*nbody complex values per frequency pair and heading pair. With ncs points of C and F (qtf_b2 = control, with
      ! npsi), also the velocity there (3 complex values per point, frequency and heading) and grad psi_k (3 per point, mode
      ! and difference frequency). With nfs inner free-surface points and npw = (2M + 1)(N_ev + 1) partial-wave coefficients
      ! (qtf_terms = full), also the fields there (5 complex values per point) and a pw_set per frequency and heading, and, with
      ! npsi, psi_k at the points (1 per point) and a pw_set per mode and difference frequency.
      ! With nsum_pairs (Phase 7c-2), the sum stores tdir_sum and tpot_sum: 6*nbody complex values per frequency pair and
      ! per slice, nsum_pairs = 2 x the heading pairs (Phase 6b: 1 x with qtf_direct = 0, where tpot_sum is the only
      ! store). npsi then counts the sum frequencies of the list too, and nfs the points of the inner mesh built for the
      ! largest sum frequency.
      ! With nmotion (Phase 6c, qtf_motion = 1), the motion stores xi2 and xi2_sum: 6*nbody complex values per frequency
      ! pair and per slice, nmotion = the routes x the heading pairs, for each of the difference and the sum that is on.
      INTEGER, INTENT(IN) :: n_panels, nwl, nbody, nfreq, nb
      INTEGER, INTENT(IN), OPTIONAL :: npsi, npairs, ncs, nfs, npw, nsum_pairs, nmotion

      qtf_stored_mb = 16D0*REAL(8*n_panels + 6*nbody + nwl, 8)*REAL(nfreq, 8)*REAL(nb, 8)/1D6
      IF (PRESENT(npsi)) qtf_stored_mb = qtf_stored_mb + 16D0*REAL(n_panels, 8)*REAL(6*nbody, 8)*REAL(npsi, 8)/1D6
      IF (PRESENT(npairs)) qtf_stored_mb = qtf_stored_mb + 16D0*REAL(6*nbody, 8)*REAL(nfreq, 8)**2*REAL(npairs, 8)/1D6
      IF (PRESENT(nsum_pairs)) qtf_stored_mb = qtf_stored_mb + 16D0*REAL(6*nbody, 8)*REAL(nfreq, 8)**2*REAL(nsum_pairs, 8)/1D6
      IF (PRESENT(nmotion)) qtf_stored_mb = qtf_stored_mb + 16D0*REAL(6*nbody, 8)*REAL(nfreq, 8)**2*REAL(nmotion, 8)/1D6
      IF (PRESENT(ncs)) THEN
         qtf_stored_mb = qtf_stored_mb + 16D0*3D0*REAL(ncs, 8)*REAL(nfreq, 8)*REAL(nb, 8)/1D6
         IF (PRESENT(npsi)) qtf_stored_mb = qtf_stored_mb + 16D0*3D0*REAL(ncs, 8)*REAL(6*nbody, 8)*REAL(npsi, 8)/1D6
      END IF
      IF (PRESENT(nfs) .AND. PRESENT(npw)) THEN
         qtf_stored_mb = qtf_stored_mb + 16D0*REAL(5*nfs + npw, 8)*REAL(nfreq, 8)*REAL(nb, 8)/1D6
         IF (PRESENT(npsi)) qtf_stored_mb = qtf_stored_mb + 16D0*REAL(6*nbody, 8)*REAL(nfs + npw, 8)*REAL(npsi, 8)/1D6
         ! qtf_fs_inner = source (qtf_set_fs_inner): grad_h psi_k too, 2 per point, mode and psi frequency
         IF (PRESENT(npsi) .AND. estimate_gpsi) qtf_stored_mb = qtf_stored_mb + &
                                                                16D0*REAL(6*nbody, 8)*2D0*REAL(nfs, 8)*REAL(npsi, 8)/1D6
      END IF
   END FUNCTION qtf_stored_mb

   SUBROUTINE qtf_pair_load(mesh, geom, kernel, i1, b1, i2, b2, total)
      ! Near-field load of the stored pair (frequency i1, heading b1; frequency i2, heading b2) under kernel, summed over
      ! the terms: total(6, nbody), force then moment about xr of each body, e^{-i omega t}, per unit amplitudes.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      INTEGER, INTENT(IN) :: kernel, i1, b1, i2, b2
      COMPLEX(8), INTENT(OUT) :: total(:, :)
      COMPLEX(8) :: terms(6, 6, mesh%nbody)

      CALL near_field_pair(mesh, geom, store(i1, b1), store(i2, b2), store_wl(:, i1, b1), store_wl(:, i2, b2), terms, kernel)
      total = SUM(terms, DIM=2)
   END SUBROUTINE qtf_pair_load

   SUBROUTINE qtf_finish(mesh, geom, refl)
      ! Writes QTF_Quadratic.10d (qtf = 1 or 3) and .10s (qtf = 2 or 3) and frees the store.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      REAL(8), INTENT(IN) :: refl

      IF (.NOT. active) RETURN
      IF (body_active .AND. diff_on) CALL potential_part(mesh, geom, refl)
      IF (direct_active) THEN
         CALL check_direct_kept()
         IF (direct_terms) CALL write_direct_terms(refl)
      END IF
      ! T+_pot with or without the direct solves (Phase 6b), once: with qtf_direct = 1 after them, as QTF_DirectTermsSum.txt
      ! has the parts of tdir_sum beside it; before the writers, which add it to .10s
      IF (body_active .AND. sum_on) CALL potential_part_sum(mesh, geom, refl)
      ! The motion from the forces of the files, after both potential parts and before the writers
      IF (mot_active) CALL qtf_motion_part(mesh, geom)
      ! Its radiated wave at the field points, into the store of QTFField before qfield_finish writes QTF_Field.txt
      IF (mot_active) CALL motion_field_part()
      IF (kind_qtf == 1 .OR. kind_qtf == 3) CALL write_file('/Wamit_format/QTF_Quadratic.10d', KERNEL_DIFF)
      IF (kind_qtf == 2 .OR. kind_qtf == 3) CALL write_file('/Wamit_format/QTF_Quadratic.10s', KERNEL_SUM)
      WRITE (*, '(a,I0,a,I0,a)') ' QTF quadratic part: ', nkept, ' frequencies and ', nbeta, ' headings'
      DEALLOCATE (store, store_wl, store_oufr)
      IF (ALLOCATED(psi_store)) DEALLOCATE (psi_store, psi_kept, psi_need, psi_reuse)
      IF (ALLOCATED(psi_is_diff)) DEALLOCATE (psi_is_diff, psi_is_sum)
      IF (ALLOCATED(tpot)) DEALLOCATE (tpot)
      IF (ALLOCATED(tdir)) DEALLOCATE (tdir)
      IF (ALLOCATED(tdir_kept)) DEALLOCATE (tdir_kept)
      IF (ALLOCATED(tdir_i2)) DEALLOCATE (tdir_i2, tdir_d)
      IF (ALLOCATED(tdir_fs)) DEALLOCATE (tdir_fs)
      CALL free_sum_stores()
      CALL free_control_surface()
      CALL free_free_surface()
      CALL free_motion_stores()
      active = .FALSE.
      body_active = .FALSE.
      direct_active = .FALSE.
      direct_terms = .FALSE.
      sum_terms = .FALSE.
      diff_on = .FALSE.
      sum_on = .FALSE.
      pot_sum_terms = .FALSE.
      sum_cut_irr = HUGE(1D0)
      max_sw = 0D0

   CONTAINS

      SUBROUTINE write_file(name, kernel)
         ! With qtf_terms = body, the difference kernel also writes QTF_Indirect.11d, row for row with .10d, and with
         ! qtf_direct = 1 QTF_Direct.12d, the same with T_dir in place of T_pot. With the sum on (Phase 7c-2) the sum kernel
         ! writes QTF_Direct.12s: the rows of .10s whose sum is kept (sum_pair_index > 0), in the order of .10s, each
         ! .10s + T+_dir, the diagonal included. Phase 6b: with the sum on, with or without qtf_direct, it also writes
         ! QTF_Indirect.11s, the same rows, each .10s + T+_pot (tpot_sum, from potential_part_sum), the diagonal included.
         ! Phase 6c, qtf_motion = 1: QTF_Motion.11d (.11s) beside .11d (.11s), and with qtf_direct = 1 QTF_Motion.12d (.12s)
         ! beside .12d (.12s), each the rows of its force file whose Omega has a regular matrix (mot_ok), in the same order
         CHARACTER(LEN=*), INTENT(IN) :: name
         INTEGER, INTENT(IN) :: kernel
         COMPLEX(8) :: total(6, mesh%nbody)
         INTEGER :: unit, unit_indirect, unit_direct, unit_dsum, unit_isum, b1, b2, i, j, b, c, k
         INTEGER :: unit_mot(2), r, mrow, h
         LOGICAL :: sum_row, mot_row

         OPEN (NEWUNIT=unit, FILE=outdir//name, STATUS='REPLACE', ACTION='WRITE')
         unit_indirect = 0
         IF (kernel == KERNEL_DIFF .AND. body_active) THEN
            OPEN (NEWUNIT=unit_indirect, FILE=outdir//'/Wamit_format/QTF_Indirect.11d', STATUS='REPLACE', ACTION='WRITE')
         END IF
         unit_direct = 0
         IF (kernel == KERNEL_DIFF .AND. direct_active) THEN
            OPEN (NEWUNIT=unit_direct, FILE=outdir//'/Wamit_format/QTF_Direct.12d', STATUS='REPLACE', ACTION='WRITE')
         END IF
         unit_dsum = 0
         IF (kernel == KERNEL_SUM .AND. direct_active .AND. sum_on) THEN
            OPEN (NEWUNIT=unit_dsum, FILE=outdir//'/Wamit_format/QTF_Direct.12s', STATUS='REPLACE', ACTION='WRITE')
         END IF
         unit_isum = 0
         IF (kernel == KERNEL_SUM .AND. body_active .AND. sum_on) THEN
            OPEN (NEWUNIT=unit_isum, FILE=outdir//'/Wamit_format/QTF_Indirect.11s', STATUS='REPLACE', ACTION='WRITE')
         END IF
         unit_mot = 0
         IF (mot_active .AND. ((kernel == KERNEL_DIFF .AND. diff_on) .OR. (kernel == KERNEL_SUM .AND. sum_on))) THEN
            DO r = 1, mot_nroute
               CALL open_motion_file(kernel, r, unit_mot(r))
            END DO
         END IF
         DO b1 = 1, nbeta
            DO b2 = 1, nbeta
               IF (all_headings_pairs == 0 .AND. b1 /= b2) CYCLE
               DO i = 1, nkept
                  DO j = 1, nkept
                     IF (kernel == KERNEL_DIFF .AND. max_dw > 0D0) THEN
                        IF (ABS(store(i, b1)%omega - store(j, b2)%omega) > max_dw*(1D0 + 1D-9)) CYCLE
                     END IF
                     CALL qtf_pair_load(mesh, geom, kernel, i, b1, j, b2, total)
                     sum_row = .FALSE.
                     IF (unit_dsum /= 0 .OR. unit_isum /= 0) sum_row = sum_pair_index(i, j) > 0
                     ! The row of the motion files: its Omega (0 on the difference diagonal) has a regular matrix
                     mot_row = .FALSE.
                     IF (unit_mot(1) /= 0) THEN
                        IF (kernel == KERNEL_DIFF) THEN
                           ! 0 is the diagonal; off it, 0 (beyond qtf_max_dw, filtered above) or -1 is no row
                           mrow = 0
                           IF (store(i, 1)%omega /= store(j, 1)%omega) THEN
                              mrow = difference_pair_index(i, j)
                              IF (mrow == 0) mrow = -1
                           END IF
                           mot_row = mrow >= 0
                        ELSE
                           mrow = sum_pair_index(i, j)
                           mot_row = sum_row .AND. mrow > 0
                        END IF
                        IF (mot_row) mot_row = mot_ok(mrow)
                     END IF
                     h = heading_pair(b1, b2)
                     DO b = 1, mesh%nbody
                        DO c = 1, 6
                           k = 6*(b - 1) + c
                           CALL write_qtf_row(unit, store_oufr(i), store_oufr(j), beta_deg(b1), beta_deg(b2), &
                                              k, total(c, b), refl, c > 3)
                           IF (sum_row .AND. unit_dsum /= 0) CALL write_qtf_row(unit_dsum, store_oufr(i), store_oufr(j), &
                              beta_deg(b1), beta_deg(b2), k, total(c, b) + tdir_sum(k, i, j, heading_pair(b1, b2)), refl, c > 3)
                           IF (sum_row .AND. unit_isum /= 0) CALL write_qtf_row(unit_isum, store_oufr(i), store_oufr(j), &
                              beta_deg(b1), beta_deg(b2), k, total(c, b) + tpot_sum(k, i, j, heading_pair(b1, b2)), refl, c > 3)
                           ! On the diagonal T_pot = T_dir = 0 and .11d, .12d are .10d to the bit (adding a zero could
                           ! flip a -0)
                           IF (unit_indirect /= 0) CALL write_sum(unit_indirect, tpot, total(c, b), i, j, b1, b2, 6*(b - 1) + c, c > 3)
                           IF (unit_direct /= 0) CALL write_sum(unit_direct, tdir, total(c, b), i, j, b1, b2, 6*(b - 1) + c, c > 3)
                           IF (.NOT. mot_row) CYCLE
                           DO r = 1, mot_nroute
                              IF (kernel == KERNEL_DIFF) THEN
                                 CALL write_motion_row(unit_mot(r), store_oufr(i), store_oufr(j), beta_deg(b1), beta_deg(b2), &
                                                       k, xi2(k, i, j, h, r), refl, c > 3)
                              ELSE
                                 CALL write_motion_row(unit_mot(r), store_oufr(i), store_oufr(j), beta_deg(b1), beta_deg(b2), &
                                                       k, xi2_sum(k, i, j, h, r), refl, c > 3)
                              END IF
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
         CLOSE (unit)
         IF (unit_indirect /= 0) CLOSE (unit_indirect)
         IF (unit_direct /= 0) CLOSE (unit_direct)
         IF (unit_dsum /= 0) CLOSE (unit_dsum)
         IF (unit_isum /= 0) CLOSE (unit_isum)
         DO r = 1, 2
            IF (unit_mot(r) /= 0) CLOSE (unit_mot(r))
         END DO

      END SUBROUTINE write_file

      SUBROUTINE open_motion_file(kernel, route, u)
         ! Opens the motion file of kernel and route and writes its header: the force file it follows and every Omega of
         ! that kernel whose matrix is singular (mot_ok false), with its rcond and the modes without restoring
         INTEGER, INTENT(IN) :: kernel, route
         INTEGER, INTENT(OUT) :: u
         CHARACTER(LEN=16), PARAMETER :: motion_names(2, 2) = RESHAPE([CHARACTER(LEN=16) :: 'QTF_Motion.11d', &
                                         'QTF_Motion.12d', 'QTF_Motion.11s', 'QTF_Motion.12s'], [2, 2])
         CHARACTER(LEN=16), PARAMETER :: force_names(2, 2) = RESHAPE([CHARACTER(LEN=16) :: 'QTF_Indirect.11d', &
                                         'QTF_Direct.12d', 'QTF_Indirect.11s', 'QTF_Direct.12s'], [2, 2])
         CHARACTER(LEN=200), ALLOCATABLE :: left(:)
         INTEGER :: m, kc

         kc = MERGE(1, 2, kernel == KERNEL_DIFF)
         OPEN (NEWUNIT=u, FILE=outdir//'/Wamit_format/'//TRIM(motion_names(route, kc)), STATUS='REPLACE', ACTION='WRITE')
         ALLOCATE (left(0))
         DO m = 0, SIZE(psi_need)
            IF (m == 0 .AND. kernel /= KERNEL_DIFF) CYCLE
            IF (m > 0) THEN
               IF (kernel == KERNEL_DIFF .AND. .NOT. psi_is_diff(m)) CYCLE
               IF (kernel == KERNEL_SUM .AND. .NOT. psi_is_sum(m)) CYCLE
            END IF
            IF (mot_ok(m) .OR. .NOT. mot_solved(m)) CYCLE
            left = [left, omega_left_text(m)]
         END DO
         CALL write_motion_header(u, TRIM(force_names(route, kc)), MERGE('indirect', 'direct  ', route == 1), left)
      END SUBROUTINE open_motion_file

      SUBROUTINE write_sum(u, tpart, base, i, j, b1, b2, k, moment)
         ! One row of .11d or .12d into unit u: the row base of .10d plus the potential part tpart of the same pair and mode
         ! k, which is left out on the diagonal (zero there; adding it could flip a -0)
         INTEGER, INTENT(IN) :: u, i, j, b1, b2, k
         COMPLEX(8), INTENT(IN) :: tpart(:, :, :, :), base
         LOGICAL, INTENT(IN) :: moment
         COMPLEX(8) :: value

         value = base
         IF (store(i, b1)%omega /= store(j, b2)%omega) value = value + tpart(k, i, j, heading_pair(b1, b2))
         CALL write_qtf_row(u, store_oufr(i), store_oufr(j), beta_deg(b1), beta_deg(b2), k, value, refl, moment)
      END SUBROUTINE write_sum

      SUBROUTINE check_direct_kept()
         ! Every difference frequency of the pairs must have been solved by the direct method before .12d is written
         REAL(8) :: missing

         missing = qtf_direct_missing()
         IF (missing > 0D0) THEN
            WRITE (*, '(a,ES14.6,a)') ' QTF direct: qtf_direct = 1, but the direct solve did not run at Omega =', &
               missing, ' rad/s; QTF_Direct.12d or .12s would be incomplete.'
            ERROR STOP 1
         END IF
      END SUBROUTINE check_direct_kept

      SUBROUTINE write_direct_terms(refl)
         ! QTF_DirectTerms.txt: for every pair with omega_i /= omega_j, by ascending |Omega| as QTF_PotentialTerms.txt, the
         ! parts of T_dir from phi_I2 and from phi_D, scaled and conjugated as .10d; with qtf_terms = full the part from
         ! phi_D of the body forcing (D_D) and that of the free-surface forcing (D_FS), a third column pair
         REAL(8), INTENT(IN) :: refl
         INTEGER :: unit, m, b1, b2, i, j, b, c, k, hp

         OPEN (NEWUNIT=unit, FILE=outdir//'/Wamit_format/QTF_DirectTerms.txt', STATUS='REPLACE', ACTION='WRITE')
         CALL write_direct_terms_header(unit, ALLOCATED(tdir_fs))
         DO m = 1, SIZE(psi_need)
            DO b1 = 1, nbeta
               DO b2 = 1, nbeta
                  IF (all_headings_pairs == 0 .AND. b1 /= b2) CYCLE
                  DO i = 1, nkept
                     DO j = 1, nkept
                        IF (difference_pair_index(i, j) /= m) CYCLE
                        DO b = 1, mesh%nbody
                           DO c = 1, 6
                              k = 6*(b - 1) + c
                              hp = heading_pair(b1, b2)
                              IF (ALLOCATED(tdir_fs)) THEN
                                 CALL write_qtf_terms_row(unit, store_oufr(i), store_oufr(j), beta_deg(b1), beta_deg(b2), k, &
                                    [tdir_i2(k, i, j, hp), tdir_d(k, i, j, hp), tdir_fs(k, i, j, hp)], refl, c > 3)
                              ELSE
                                 CALL write_qtf_terms_row(unit, store_oufr(i), store_oufr(j), beta_deg(b1), beta_deg(b2), k, &
                                    [tdir_i2(k, i, j, hp), tdir_d(k, i, j, hp)], refl, c > 3)
                              END IF
                           END DO
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
         CLOSE (unit)
      END SUBROUTINE write_direct_terms

   END SUBROUTINE qtf_finish

   SUBROUTINE potential_part(mesh, geom, refl)
      ! T_pot of every pair of the difference QTF with omega_i /= omega_j within qtf_max_dw, into tpot, Omega by Omega:
      ! grad psi and psi at the waterline are built once for each difference frequency. With qtf_terms = full the
      ! free-surface terms FS_inner and FS_outer (fs_pair) follow the five of potential_pair, and the time of each Omega goes
      ! to the log. With write_terms >= 1 the terms go to QTF_PotentialTerms.txt, in that order (ascending |Omega|, not the
      ! .11d order). Stops if psi of a needed Omega, or a free-surface field, was not stored. Pairs with NaN or Inf terms
      ! are written as they are, counted in qtf_nonfinite_pairs() and logged per Omega.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      REAL(8), INTENT(IN) :: refl
      INTEGER, ALLOCATABLE :: pair_m(:, :)                 ! (nkept, nkept), index of |omega_i - omega_j| in psi_need, or 0
      COMPLEX(8), ALLOCATABLE :: grad(:, :, :), pwl(:, :), terms(:, :)
      TYPE(fs_point_field), ALLOCATABLE :: fs_inc(:, :, :)   ! (n_fs, nkept, nbeta), incident fields at the inner points
      ! qtf_fs_inner = source: the states on the boundary of the inner free surface, psi_k there, the line weights
      TYPE(fs_wl_state), ALLOCATABLE :: fs_bnd(:, :)
      COMPLEX(8), ALLOCATABLE :: psi_b(:, :)
      REAL(8), ALLOCATABLE :: xe(:, :), wbnd(:)
      REAL(8) :: dw, g_wave, rf, rf_max, fs_seconds
      INTEGER(8) :: t0, t1, rate
      INTEGER :: unit, i, j, m, b1, b2, b, c, k, nmode, nbad, nterms, npairs_fs

      nmode = 6*mesh%nbody
      ALLOCATE (pair_m(nkept, nkept))
      pair_m = 0
      DO i = 1, nkept
         DO j = 1, nkept
            m = difference_pair_index(i, j)
            IF (m == 0) CYCLE
            IF (m < 0) THEN
               WRITE (*, '(a,ES14.6,a)') ' QTF: Omega =', ABS(store(i, 1)%omega - store(j, 1)%omega), &
                  ' rad/s is not among the difference frequencies of psi.'
               ERROR STOP 1
            END IF
            IF (.NOT. psi_kept(m)) THEN
               WRITE (*, '(a,ES14.6,a)') ' QTF: psi at Omega =', psi_need(m), ' rad/s was not stored.'
               ERROR STOP 1
            END IF
            IF (b2_control) THEN
               IF (.NOT. cs_psi_kept(m)) THEN
                  WRITE (*, '(a,ES14.6,a)') ' QTF: grad psi on the control surface at Omega =', psi_need(m), &
                     ' rad/s was not stored.'
                  ERROR STOP 1
               END IF
               IF (.NOT. (ALL(cs_u_kept(i, :)) .AND. ALL(cs_u_kept(j, :)))) THEN
                  WRITE (*, '(a)') ' QTF: the velocity on the control surface of a kept frequency was not stored.'
                  ERROR STOP 1
               END IF
            END IF
            IF (fs_active) THEN
               IF (.NOT. fs_psi_kept(m)) THEN
                  WRITE (*, '(a,ES14.6,a)') ' QTF: psi on the inner free surface at Omega =', psi_need(m), &
                     ' rad/s was not stored.'
                  ERROR STOP 1
               END IF
               IF (.NOT. (ALL(fs_kept(i, :)) .AND. ALL(fs_kept(j, :)))) THEN
                  WRITE (*, '(a)') ' QTF: the free-surface fields of a kept frequency were not stored.'
                  ERROR STOP 1
               END IF
            END IF
            pair_m(i, j) = m
         END DO
      END DO

      unit = 0
      IF (pot_terms) THEN
         OPEN (NEWUNIT=unit, FILE=outdir//'/Wamit_format/QTF_PotentialTerms.txt', STATUS='REPLACE', ACTION='WRITE')
         CALL write_qtf_terms_header(unit, fs_active)
      END IF
      nterms = NPOT_TERMS
      IF (fs_active) THEN
         nterms = NPOT_TERMS + 2
         CALL fs_incident_all(fs_inc)
         IF (fs_source) THEN
            CALL fs_boundary_all(mesh, geom, fs_bnd, xe, wbnd)
            ALLOCATE (psi_b(SIZE(wbnd), nmode))
         END IF
      END IF
      ALLOCATE (grad(3, mesh%n, nmode), pwl(geom%nwl, nmode), terms(nmode, nterms))
      fs_seconds = 0D0
      DO m = 1, SIZE(psi_need)
         IF (.NOT. ANY(pair_m == m)) CYCLE
         CALL psi_surface_gradient(mesh, geom, psi_store(:, :, m), grad)
         CALL psi_waterline(mesh, geom, psi_store(:, :, m), grad, pwl)
         IF (fs_active .AND. fs_source) CALL fs_boundary_psi(geom%nwl, pwl, fs_pw_psi(:, m), xe, psi_b)
         nbad = 0
         npairs_fs = 0
         rf_max = 0D0
         IF (fs_active) CALL SYSTEM_CLOCK(t0, rate)
         DO b1 = 1, nbeta
            DO b2 = 1, nbeta
               IF (all_headings_pairs == 0 .AND. b1 /= b2) CYCLE
               DO i = 1, nkept
                  DO j = 1, nkept
                     IF (pair_m(i, j) /= m) CYCLE
                     ASSOCIATE (s1 => store(i, b1), s2 => store(j, b2))
                        ! phi_I2 with the gravity of the first-order incident wave, which has the k of DISPERSION
                        ! (wave_gravity): both orders of the incident wave then have one dispersion relation (ruling R8)
                        g_wave = wave_gravity(s1%omega, s1%k, s1%h)
                        IF (b2_control .AND. ALLOCATED(edge_w)) THEN
                           CALL potential_pair(mesh, geom, s1, s2, store_wl(:, i, b1), store_wl(:, j, b2), &
                                               psi_store(:, :, m), grad, pwl, psi_need(m), s1%rho, g_wave, &
                                               terms(:, 1:NPOT_TERMS), &
                                               cs=cs, grad_psi_cs=cs_gpsi(:, :, :, m), u1_cs=cs_u(:, :, i, b1), &
                                               u2_cs=cs_u(:, :, j, b2), edge_w=edge_w, omega_tol=psi_tol)
                        ELSE IF (b2_control) THEN
                           CALL potential_pair(mesh, geom, s1, s2, store_wl(:, i, b1), store_wl(:, j, b2), &
                                               psi_store(:, :, m), grad, pwl, psi_need(m), s1%rho, g_wave, &
                                               terms(:, 1:NPOT_TERMS), &
                                               cs=cs, grad_psi_cs=cs_gpsi(:, :, :, m), u1_cs=cs_u(:, :, i, b1), &
                                               u2_cs=cs_u(:, :, j, b2), omega_tol=psi_tol)
                        ELSE
                           CALL potential_pair(mesh, geom, s1, s2, store_wl(:, i, b1), store_wl(:, j, b2), &
                                               psi_store(:, :, m), grad, pwl, psi_need(m), s1%rho, g_wave, &
                                               terms(:, 1:NPOT_TERMS), omega_tol=psi_tol)
                        END IF
                     END ASSOCIATE
                     IF (fs_active) THEN
                        CALL fs_pair(i, b1, j, b2, m, terms(:, NPOT_TERMS + 1), terms(:, NPOT_TERMS + 2), rf)
                        npairs_fs = npairs_fs + 1
                        rf_max = MAX(rf_max, rf)
                     END IF
                     tpot(:, i, j, heading_pair(b1, b2)) = SUM(terms, DIM=2)
                     IF (.NOT. ALL(IEEE_IS_FINITE(REAL(terms, 8)) .AND. IEEE_IS_FINITE(AIMAG(terms)))) nbad = nbad + 1
                     IF (unit == 0) CYCLE
                     DO b = 1, mesh%nbody
                        DO c = 1, 6
                           k = 6*(b - 1) + c
                           CALL write_qtf_terms_row(unit, store_oufr(i), store_oufr(j), beta_deg(b1), beta_deg(b2), k, &
                                                    terms(k, :), refl, c > 3)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
         IF (fs_active) THEN
            CALL SYSTEM_CLOCK(t1)
            fs_seconds = fs_seconds + REAL(t1 - t0, 8)/REAL(rate, 8)
            WRITE (*, '(a,F10.5,a,I0,a,F10.1,a,F10.2,a)') ' QTF free surface integral: Omega', psi_need(m), ' rad/s, ', &
               npairs_fs, ' pairs, r_far up to', rf_max, ' m,', REAL(t1 - t0, 8)/REAL(rate, 8), ' s'
         END IF
         IF (nbad > 0) THEN
            n_bad_pairs = n_bad_pairs + nbad
            WRITE (*, '(a,I0,a,ES14.6,a)') ' QTF error: T_pot has NaN or Inf at ', nbad, ' frequency and heading pairs with'// &
               ' Omega =', psi_need(m), ' rad/s; those rows of QTF_Indirect.11d are not valid.'
         END IF
      END DO
      IF (unit /= 0) CLOSE (unit)
      IF (fs_active) WRITE (*, '(a,F10.2,a)') ' QTF free surface integral: ', fs_seconds, ' s in all (body terms included)'

   CONTAINS

      SUBROUTINE fs_pair(ia, ba, ja, bb, mm, t_in, t_out, r_far)
         ! FS_inner and FS_outer of the stored pair (frequency ia, heading ba; frequency ja, heading bb) with psi_k of the
         ! difference frequency mm (fs_pair_terms: Omega = 0 gives none, omega_ia < omega_ja is evaluated as (ja, ia) and
         ! conjugated), with the gravity of the first-order states
         INTEGER, INTENT(IN) :: ia, ba, ja, bb, mm
         COMPLEX(8), INTENT(OUT) :: t_in(:), t_out(:)
         REAL(8), INTENT(OUT) :: r_far

         ASSOCIATE (s1 => store(ia, ba), s2 => store(ja, bb))
            IF (fs_source) THEN
               CALL fs_pair_terms(s1%omega, fs_f(:, ia, ba), fs_inc(:, ia, ba), fs_pw(ia, ba), s1%beta, s1%k, &
                                  s2%omega, fs_f(:, ja, bb), fs_inc(:, ja, bb), fs_pw(ja, bb), s2%beta, s2%k, &
                                  fs_psi(:, :, mm), fs_pw_psi(:, mm), fs_store_pts%w, fs_store_np, s1%g, s1%rho, t_in, &
                                  t_out, r_far, gpsi=fs_gpsi(:, :, :, mm), b1=fs_bnd(ia, ba), b2=fs_bnd(ja, bb), &
                                  psi_b=psi_b, wb=wbnd)
            ELSE
               CALL fs_pair_terms(s1%omega, fs_f(:, ia, ba), fs_inc(:, ia, ba), fs_pw(ia, ba), s1%beta, s1%k, &
                                  s2%omega, fs_f(:, ja, bb), fs_inc(:, ja, bb), fs_pw(ja, bb), s2%beta, s2%k, &
                                  fs_psi(:, :, mm), fs_pw_psi(:, mm), fs_store_pts%w, fs_store_np, s1%g, s1%rho, t_in, &
                                  t_out, r_far)
            END IF
         END ASSOCIATE
      END SUBROUTINE fs_pair

   END SUBROUTINE potential_part

   SUBROUTINE potential_part_sum(mesh, geom, refl)
      ! The indirect sum load T+_pot (Phase 7c-2, ruling R1 of its plan) of the canonical half of the sum pairs
      ! (qtf_direct_sum_pairs), Omega by Omega over the sum frequencies of the list, into both slots of tpot_sum (R6):
      ! potential_pair with KERNEL_SUM from the same psi, grad psi, psi at the waterline, psi on C and F and psi on the
      ! inner free surface as potential_part, indexed by the union list, and with qtf_terms = full the free-surface terms
      ! of fs_pair_terms with KERNEL_SUM (ordered w1 >= w2 there, no conjugate on the swap). With write_terms >= 1 every
      ! row of .12s goes to QTF_DirectTermsSum.txt: the parts of T+_dir (D_I2, D_D and with full D_FS) beside P_BODY, the
      ! sum of the five body-forcing terms, and with full P_FS = FS_inner + FS_outer, in ascending Omega; a pair and its
      ! mirror (j, b2, i, b1) both have a row, with the same values. With write_terms >= 1 the five body terms and with
      ! full FS_inner and FS_outer of every such row also go to QTF_PotentialTermsSum.txt (Phase 6b), in the same order,
      ! with or without the direct solve. Nothing here reads a direct-solve store (tdir_sum_* only for QTF_DirectTermsSum.txt),
      ! so the two routes fill the same slots of tpot_sum from the same psi (R2). Stops if psi of a needed Omega, or a
      ! free-surface field, was not stored. Pairs with NaN or Inf are counted in qtf_nonfinite_pairs() and logged per Omega.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      REAL(8), INTENT(IN) :: refl
      INTEGER, ALLOCATABLE :: pairs(:, :)
      COMPLEX(8), ALLOCATABLE :: grad(:, :, :), pwl(:, :), terms(:, :), tfs_in(:), tfs_out(:)
      TYPE(fs_point_field), ALLOCATABLE :: fs_inc(:, :, :)   ! (n_fs, nkept, nbeta), incident fields at the inner points
      ! qtf_fs_inner = source: as potential_part
      TYPE(fs_wl_state), ALLOCATABLE :: fs_bnd(:, :)
      COMPLEX(8), ALLOCATABLE :: psi_b(:, :)
      REAL(8), ALLOCATABLE :: xe(:, :), wbnd(:)
      REAL(8) :: g_wave, rf, rf_max, fs_seconds
      INTEGER(8) :: t0, t1, rate
      INTEGER :: unit, unit_ps, m, p, i, j, b1, b2, h, hm, nmode, nbad
      CHARACTER(LEN=120) :: bad_files   ! the files the rows of a non-finite pair are in

      nmode = 6*mesh%nbody
      bad_files = 'QTF_Indirect.11s'
      IF (pot_sum_terms) bad_files = TRIM(bad_files)//', QTF_PotentialTermsSum.txt'
      IF (sum_terms) bad_files = TRIM(bad_files)//', QTF_DirectTermsSum.txt'
      unit_ps = 0
      IF (pot_sum_terms) THEN
         OPEN (NEWUNIT=unit_ps, FILE=outdir//'/Wamit_format/QTF_PotentialTermsSum.txt', STATUS='REPLACE', ACTION='WRITE')
         CALL write_potential_sum_terms_header(unit_ps, fs_active)
      END IF
      unit = 0
      IF (sum_terms) THEN
         OPEN (NEWUNIT=unit, FILE=outdir//'/Wamit_format/QTF_DirectTermsSum.txt', STATUS='REPLACE', ACTION='WRITE')
         CALL write_direct_sum_terms_header(unit, fs_active)
      END IF
      IF (fs_active) CALL fs_incident_all(fs_inc)
      IF (fs_active .AND. fs_source) THEN
         CALL fs_boundary_all(mesh, geom, fs_bnd, xe, wbnd)
         ALLOCATE (psi_b(SIZE(wbnd), nmode))
      END IF
      ALLOCATE (grad(3, mesh%n, nmode), pwl(geom%nwl, nmode), terms(nmode, NPOT_TERMS), tfs_in(nmode), tfs_out(nmode))
      tfs_in = (0D0, 0D0)
      tfs_out = (0D0, 0D0)
      fs_seconds = 0D0
      DO m = 1, SIZE(psi_need)
         IF (.NOT. psi_is_sum(m)) CYCLE
         CALL qtf_direct_sum_pairs(psi_need(m), pairs)
         IF (SIZE(pairs, 2) == 0) CYCLE
         CALL check_stored(m)
         CALL psi_surface_gradient(mesh, geom, psi_store(:, :, m), grad)
         CALL psi_waterline(mesh, geom, psi_store(:, :, m), grad, pwl)
         IF (fs_active .AND. fs_source) CALL fs_boundary_psi(geom%nwl, pwl, fs_pw_psi(:, m), xe, psi_b)
         nbad = 0
         rf_max = 0D0
         IF (fs_active) CALL SYSTEM_CLOCK(t0, rate)
         DO p = 1, SIZE(pairs, 2)
            i = pairs(1, p)
            b1 = pairs(2, p)
            j = pairs(3, p)
            b2 = pairs(4, p)
            ASSOCIATE (s1 => store(i, b1), s2 => store(j, b2))
               ! phi_I2 with the gravity of the first-order incident wave, as potential_part (ruling R8 of Phase 6)
               g_wave = wave_gravity(s1%omega, s1%k, s1%h)
               IF (b2_control .AND. ALLOCATED(edge_w)) THEN
                  CALL potential_pair(mesh, geom, s1, s2, store_wl(:, i, b1), store_wl(:, j, b2), psi_store(:, :, m), grad, &
                                      pwl, psi_need(m), s1%rho, g_wave, terms, cs=cs, grad_psi_cs=cs_gpsi(:, :, :, m), &
                                      u1_cs=cs_u(:, :, i, b1), u2_cs=cs_u(:, :, j, b2), edge_w=edge_w, omega_tol=psi_tol, &
                                      kernel=KERNEL_SUM)
               ELSE IF (b2_control) THEN
                  CALL potential_pair(mesh, geom, s1, s2, store_wl(:, i, b1), store_wl(:, j, b2), psi_store(:, :, m), grad, &
                                      pwl, psi_need(m), s1%rho, g_wave, terms, cs=cs, grad_psi_cs=cs_gpsi(:, :, :, m), &
                                      u1_cs=cs_u(:, :, i, b1), u2_cs=cs_u(:, :, j, b2), omega_tol=psi_tol, kernel=KERNEL_SUM)
               ELSE
                  CALL potential_pair(mesh, geom, s1, s2, store_wl(:, i, b1), store_wl(:, j, b2), psi_store(:, :, m), grad, &
                                      pwl, psi_need(m), s1%rho, g_wave, terms, omega_tol=psi_tol, kernel=KERNEL_SUM)
               END IF
               IF (fs_active .AND. fs_source) THEN
                  CALL fs_pair_terms(s1%omega, fs_f(:, i, b1), fs_inc(:, i, b1), fs_pw(i, b1), s1%beta, s1%k, &
                                     s2%omega, fs_f(:, j, b2), fs_inc(:, j, b2), fs_pw(j, b2), s2%beta, s2%k, &
                                     fs_psi(:, :, m), fs_pw_psi(:, m), fs_store_pts%w, fs_store_np, s1%g, s1%rho, tfs_in, &
                                     tfs_out, rf, kernel=KERNEL_SUM, gpsi=fs_gpsi(:, :, :, m), b1=fs_bnd(i, b1), &
                                     b2=fs_bnd(j, b2), psi_b=psi_b, wb=wbnd)
                  rf_max = MAX(rf_max, rf)
               ELSE IF (fs_active) THEN
                  CALL fs_pair_terms(s1%omega, fs_f(:, i, b1), fs_inc(:, i, b1), fs_pw(i, b1), s1%beta, s1%k, &
                                     s2%omega, fs_f(:, j, b2), fs_inc(:, j, b2), fs_pw(j, b2), s2%beta, s2%k, &
                                     fs_psi(:, :, m), fs_pw_psi(:, m), fs_store_pts%w, fs_store_np, s1%g, s1%rho, tfs_in, &
                                     tfs_out, rf, kernel=KERNEL_SUM)
                  rf_max = MAX(rf_max, rf)
               END IF
            END ASSOCIATE
            h = heading_pair(b1, b2)
            hm = heading_pair(b2, b1)
            tpot_sum(:, i, j, h) = SUM(terms, DIM=2) + tfs_in + tfs_out
            tpot_sum(:, j, i, hm) = tpot_sum(:, i, j, h)
            IF (.NOT. (ALL(IEEE_IS_FINITE(REAL(terms, 8)) .AND. IEEE_IS_FINITE(AIMAG(terms))) .AND. &
                       ALL(IEEE_IS_FINITE(REAL(tpot_sum(:, i, j, h), 8)) .AND. IEEE_IS_FINITE(AIMAG(tpot_sum(:, i, j, h)))))) &
               nbad = nbad + 1
            IF (unit_ps /= 0) THEN
               CALL write_pot_rows(i, b1, j, b2)
               IF (i /= j .OR. b1 /= b2) CALL write_pot_rows(j, b2, i, b1)
            END IF
            IF (unit == 0) CYCLE
            CALL write_rows(i, b1, j, b2)
            IF (i /= j .OR. b1 /= b2) CALL write_rows(j, b2, i, b1)
         END DO
         IF (fs_active) THEN
            CALL SYSTEM_CLOCK(t1)
            fs_seconds = fs_seconds + REAL(t1 - t0, 8)/REAL(rate, 8)
            WRITE (*, '(a,F10.5,a,I0,a,F10.1,a,F10.2,a)') ' QTF sum indirect check: Omega', psi_need(m), ' rad/s, ', &
               SIZE(pairs, 2), ' pairs, r_far up to', rf_max, ' m,', REAL(t1 - t0, 8)/REAL(rate, 8), ' s'
         END IF
         IF (nbad > 0) THEN
            n_bad_pairs = n_bad_pairs + nbad
            WRITE (*, '(a,I0,a,ES14.6,a)') ' QTF error: T+_pot has NaN or Inf at ', nbad, ' frequency and heading pairs'// &
               ' with Omega =', psi_need(m), ' rad/s; those rows of '//TRIM(bad_files)//' are not valid.'
         END IF
      END DO
      IF (unit /= 0) CLOSE (unit)
      IF (unit_ps /= 0) CLOSE (unit_ps)
      IF (fs_active) WRITE (*, '(a,F10.2,a)') ' QTF sum indirect check: ', fs_seconds, ' s in all (body terms included)'

   CONTAINS

      SUBROUTINE write_pot_rows(ia, ba, ja, bb)
         ! The rows of the pair (ia, ba, ja, bb) in QTF_PotentialTermsSum.txt, one per mode: the five body terms and with
         ! full FS_inner and FS_outer, the values of the pair (i, b1, j, b2) being solved (the mirror has the same)
         INTEGER, INTENT(IN) :: ia, ba, ja, bb
         INTEGER :: b, c, k

         DO b = 1, mesh%nbody
            DO c = 1, 6
               k = 6*(b - 1) + c
               IF (fs_active) THEN
                  CALL write_qtf_terms_row(unit_ps, store_oufr(ia), store_oufr(ja), beta_deg(ba), beta_deg(bb), k, &
                     [terms(k, :), tfs_in(k), tfs_out(k)], refl, c > 3)
               ELSE
                  CALL write_qtf_terms_row(unit_ps, store_oufr(ia), store_oufr(ja), beta_deg(ba), beta_deg(bb), k, &
                     terms(k, :), refl, c > 3)
               END IF
            END DO
         END DO
      END SUBROUTINE write_pot_rows

      SUBROUTINE write_rows(ia, ba, ja, bb)
         ! The rows of the pair (ia, ba, ja, bb) in QTF_DirectTermsSum.txt, one per mode, the values of slot h
         INTEGER, INTENT(IN) :: ia, ba, ja, bb
         INTEGER :: b, c, k

         DO b = 1, mesh%nbody
            DO c = 1, 6
               k = 6*(b - 1) + c
               IF (fs_active) THEN
                  CALL write_qtf_terms_row(unit, store_oufr(ia), store_oufr(ja), beta_deg(ba), beta_deg(bb), k, &
                     [tdir_sum_i2(k, i, j, h), tdir_sum_d(k, i, j, h), tdir_sum_fs(k, i, j, h), SUM(terms(k, :)), &
                      tfs_in(k) + tfs_out(k)], refl, c > 3)
               ELSE
                  CALL write_qtf_terms_row(unit, store_oufr(ia), store_oufr(ja), beta_deg(ba), beta_deg(bb), k, &
                     [tdir_sum_i2(k, i, j, h), tdir_sum_d(k, i, j, h), SUM(terms(k, :))], refl, c > 3)
               END IF
            END DO
         END DO
      END SUBROUTINE write_rows

      SUBROUTINE check_stored(mm)
         ! psi and the fields the pairs of the sum frequency mm read were stored
         INTEGER, INTENT(IN) :: mm

         IF (.NOT. psi_kept(mm)) THEN
            WRITE (*, '(a,ES14.6,a)') ' QTF: psi at the sum frequency Omega =', psi_need(mm), ' rad/s was not stored.'
            ERROR STOP 1
         END IF
         IF (b2_control) THEN
            IF (.NOT. (cs_psi_kept(mm) .AND. ALL(cs_u_kept))) THEN
               WRITE (*, '(a,ES14.6,a)') ' QTF: the fields on the control surface at the sum frequency Omega =', &
                  psi_need(mm), ' rad/s, or of a kept frequency, were not stored.'
               ERROR STOP 1
            END IF
         END IF
         IF (fs_active) THEN
            IF (.NOT. (fs_psi_kept(mm) .AND. ALL(fs_kept))) THEN
               WRITE (*, '(a,ES14.6,a)') ' QTF: the free-surface fields at the sum frequency Omega =', psi_need(mm), &
                  ' rad/s, or of a kept frequency, were not stored.'
               ERROR STOP 1
            END IF
         END IF
      END SUBROUTINE check_stored

   END SUBROUTINE potential_part_sum

   SUBROUTINE free_motion_stores()
      mot_active = .FALSE.
      mot_nroute = 0
      IF (ALLOCATED(mot_mass)) DEALLOCATE (mot_mass, mot_c, mot_k, mot_b)
      IF (ALLOCATED(xi2)) DEALLOCATE (xi2)
      IF (ALLOCATED(xi2_sum)) DEALLOCATE (xi2_sum)
      IF (ALLOCATED(mot_ok)) DEALLOCATE (mot_ok, mot_rcond, mot_free, mot_solved)
      IF (ALLOCATED(mot_psi_pts)) DEALLOCATE (mot_psi_pts, mot_psi_kept)
   END SUBROUTINE free_motion_stores

   SUBROUTINE motion_field_part()
      ! The radiated wave of the second-order motion at the field points (R7), after qtf_motion_part: for every row of
      ! QTF_Field.txt with a direct part (K = 1: omega_i /= omega_j within qtf_max_dw, the pairs of qtf_direct_pairs; K = 2:
      ! the kept sums, sum_pair_index > 0) whose Omega has a regular matrix (mot_ok), phi_R = radiated_potential(Omega,
      ! xi, psi at the points) with Omega signed and xi of the direct route (route 2, QTF_Motion.12d and .12s, whose force
      ! is the one of the direct part), handed to qfield_add_radiated. A row of an Omega left out of the motion (R6) gets
      ! none, so QTF_Field.txt marks it has_rad = 0. The difference diagonal has no direct part (and phi_R = 0 at
      ! Omega = 0). No-op unless QTFField takes the radiated wave.
      COMPLEX(8), ALLOCATABLE :: phi_r(:)
      INTEGER :: b1, b2, i, j, m, h
      REAL(8) :: om

      IF (.NOT. (mot_active .AND. qfield_radiated_active())) RETURN
      ! mot_psi_pts is allocated by the first direct solve; a list with no Omega (one frequency, or qtf_max_dw below the
      ! step) has none and no row with a direct part, so the loop below does nothing and the files are still written
      IF (mot_nroute /= 2) THEN
         WRITE (*, '(a)') ' QTF motion: the radiated wave in QTF_Field.txt needs qtf_direct = 1.'
         ERROR STOP 1
      END IF
      DO b1 = 1, nbeta
         DO b2 = 1, nbeta
            IF (all_headings_pairs == 0 .AND. b1 /= b2) CYCLE
            h = heading_pair(b1, b2)
            DO i = 1, nkept
               DO j = 1, nkept
                  IF (diff_on) THEN
                     m = difference_pair_index(i, j)
                     IF (m > 0) THEN
                        IF (mot_ok(m)) THEN
                           CALL check_psi(m)
                           ! Omega as QTFField takes it for the row (field_direct), signed
                           om = store(i, b1)%omega - store(j, b2)%omega
                           phi_r = radiated_potential(om, xi2(:, i, j, h, 2), mot_psi_pts(:, :, m))
                           CALL qfield_add_radiated(i, b1, j, b2, phi_r, KERNEL_DIFF)
                        END IF
                     END IF
                  END IF
                  IF (sum_on) THEN
                     m = sum_pair_index(i, j)
                     IF (m > 0) THEN
                        IF (mot_ok(m)) THEN
                           CALL check_psi(m)
                           om = store(i, b1)%omega + store(j, b2)%omega
                           phi_r = radiated_potential(om, xi2_sum(:, i, j, h, 2), mot_psi_pts(:, :, m))
                           CALL qfield_add_radiated(i, b1, j, b2, phi_r, KERNEL_SUM)
                        END IF
                     END IF
                  END IF
               END DO
            END DO
         END DO
      END DO
   CONTAINS
      SUBROUTINE check_psi(mm)
         INTEGER, INTENT(IN) :: mm
         LOGICAL :: kept

         ! No store at all (no direct solve ran) counts as missing; the test is split so mot_psi_kept is not indexed then
         kept = ALLOCATED(mot_psi_kept)
         IF (kept) kept = mot_psi_kept(mm)
         IF (.NOT. kept) THEN
            WRITE (*, '(a,ES14.6,a)') ' QTF motion: psi at the field points is missing at Omega =', psi_need(mm), ' rad/s.'
            ERROR STOP 1
         END IF
      END SUBROUTINE check_psi
   END SUBROUTINE motion_field_part

   SUBROUTINE qtf_motion_part(mesh, geom)
      ! The second-order motion of every row of the difference and the sum files (Phase 6c, rulings R4 to R6 of its plan),
      ! Omega by Omega, into xi2 and xi2_sum. Call after potential_part, potential_part_sum and the direct solves, which
      ! fill the forces; qtf_finish does. Per Omega of the list with rows: A and B from psi (radiation_coefficients, with
      ! the rho of the states), the matrix of (2) of QTFMotion (motion_matrix), and per kernel and route one solve over the
      ! force columns of all the pairs of that Omega (motion_solve_omega). Each Omega is logged with its rcond; one below
      ! 1e-12 (mot_ok false) also with the modes without restoring, and its rows are not written.
      !
      ! Omega = 0, the difference diagonal (R6): the equation is static, (C + K_ext) xi2 = F. motion_matrix at Omega = 0
      ! gives exactly that: the -Omega^2 (M + A) and -i Omega (B + B_ext) terms vanish, so A and B are passed as zero.
      !
      ! The sign of Omega (R4). The force of the difference row (i, j) is the coefficient of e^{-i Omega t} with
      ! Omega = omega_i - omega_j, negative when omega_i < omega_j. The matrix of (2) at Omega < 0 needs A and B at -|Omega|.
      ! The radiation force (Omega^2 A + i Omega B) xi is the transform of a real response, so its coefficient at -Omega is
      ! the conjugate of that at Omega: A and B are even in Omega. With M, C, K_ext and B_ext real, the matrix at -|Omega| is
      ! then conj(lhs), lhs the matrix at |Omega|. conj(lhs) xi = F is lhs conj(xi) = conj(F): the column conj(F) is solved
      ! with the factors of lhs, and xi is the conjugate of the result. So an Omega has one matrix for both signs.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      INTEGER, ALLOCATABLE :: diff_m(:, :), sum_m(:, :)   ! (nkept, nkept): Omega index of each row, -1 (diff) or 0 (sum): none
      REAL(8), ALLOCATABLE :: a(:, :), b(:, :)
      COMPLEX(8), ALLOCATABLE :: lhs(:, :)
      REAL(8) :: rho
      INTEGER :: nmode, i, j, m

      IF (.NOT. mot_active) RETURN
      nmode = 6*mesh%nbody
      IF (nmode /= SIZE(mot_mass, 1) .OR. nkept /= SIZE(store, 1)) THEN
         WRITE (*, '(a)') ' QTF motion: the mesh or the kept frequencies do not fit the store of qtf_init.'
         ERROR STOP 1
      END IF
      ! rho of the first-order states, which potential_pair uses for the forces too
      rho = store(1, 1)%rho
      ALLOCATE (diff_m(nkept, nkept), sum_m(nkept, nkept), a(nmode, nmode), b(nmode, nmode), lhs(nmode, nmode))
      diff_m = -1
      sum_m = 0
      DO i = 1, nkept
         DO j = 1, nkept
            IF (diff_on) THEN
               IF (store(i, 1)%omega == store(j, 1)%omega) THEN
                  diff_m(i, j) = 0
               ELSE
                  ! 0 (beyond qtf_max_dw: not a row) and -1 (not in the list: potential_part has stopped) are no row
                  diff_m(i, j) = difference_pair_index(i, j)
                  IF (diff_m(i, j) == 0) diff_m(i, j) = -1
               END IF
            END IF
            IF (sum_on) sum_m(i, j) = MAX(sum_pair_index(i, j), 0)
         END DO
      END DO
      IF (diff_on) THEN
         a = 0D0
         b = 0D0
         CALL motion_matrix(0D0, mot_mass, a, b, mot_b, mot_c, mot_k, lhs)
         CALL solve_kernel(KERNEL_DIFF, 0)
         CALL log_omega(0, 0D0)
      END IF
      DO m = 1, SIZE(psi_need)
         IF (.NOT. (ANY(diff_m == m) .OR. ANY(sum_m == m))) CYCLE
         ! An unkept psi is zero in the store and would give A = B = 0 without a word
         IF (.NOT. psi_kept(m)) THEN
            WRITE (*, '(a,ES14.6,a)') ' QTF motion: psi at Omega =', psi_need(m), ' rad/s is not kept for A and B.'
            ERROR STOP 1
         END IF
         CALL radiation_coefficients(mesh, psi_store(:, :, m), psi_need(m), rho, a, b)
         CALL motion_matrix(psi_need(m), mot_mass, a, b, mot_b, mot_c, mot_k, lhs)
         IF (ANY(diff_m == m)) CALL solve_kernel(KERNEL_DIFF, m)
         IF (ANY(sum_m == m)) CALL solve_kernel(KERNEL_SUM, m)
         CALL log_omega(m, psi_need(m))
      END DO

   CONTAINS

      SUBROUTINE solve_kernel(kernel, mm)
         ! Every row of kernel at the Omega of index mm, one motion_solve_omega per route on lhs
         INTEGER, INTENT(IN) :: kernel, mm
         INTEGER, ALLOCATABLE :: pairs(:, :)
         LOGICAL, ALLOCATABLE :: negative(:)
         COMPLEX(8), ALLOCATABLE :: forces(:, :), xs(:, :)
         COMPLEX(8) :: f(nmode)
         REAL(8) :: rcond
         LOGICAL :: ok
         INTEGER :: b1, b2, ii, jj, np, p, r, h

         ALLOCATE (pairs(4, nkept*nkept*heading_pairs(nbeta)), negative(nkept*nkept*heading_pairs(nbeta)))
         np = 0
         DO b1 = 1, nbeta
            DO b2 = 1, nbeta
               IF (all_headings_pairs == 0 .AND. b1 /= b2) CYCLE
               DO ii = 1, nkept
                  DO jj = 1, nkept
                     IF (kernel == KERNEL_DIFF) THEN
                        IF (diff_m(ii, jj) /= mm) CYCLE
                     ELSE
                        IF (sum_m(ii, jj) /= mm) CYCLE
                     END IF
                     np = np + 1
                     pairs(:, np) = [ii, b1, jj, b2]
                     negative(np) = kernel == KERNEL_DIFF .AND. store(ii, 1)%omega < store(jj, 1)%omega
                  END DO
               END DO
            END DO
         END DO
         ALLOCATE (forces(nmode, np), xs(nmode, np))
         DO r = 1, mot_nroute
            DO p = 1, np
               CALL pair_force(mesh, geom, kernel, r, pairs(1, p), pairs(2, p), pairs(3, p), pairs(4, p), f)
               forces(:, p) = MERGE(CONJG(f), f, negative(p))
            END DO
            CALL motion_solve_omega(lhs, forces, xs, ok, rcond)
            DO p = 1, np
               h = heading_pair(pairs(2, p), pairs(4, p))
               IF (kernel == KERNEL_DIFF) THEN
                  xi2(:, pairs(1, p), pairs(3, p), h, r) = MERGE(CONJG(xs(:, p)), xs(:, p), negative(p))
               ELSE
                  xi2_sum(:, pairs(1, p), pairs(3, p), h, r) = xs(:, p)
               END IF
            END DO
         END DO
         ! The matrix is the same for every route and kernel of the Omega, and so are ok and rcond
         mot_ok(mm) = ok
         mot_rcond(mm) = rcond
         mot_solved(mm) = .TRUE.
         IF (.NOT. ok) mot_free(:, mm) = free_modes(lhs)
      END SUBROUTINE solve_kernel

      SUBROUTINE log_omega(mm, om)
         ! One log line per Omega: its rcond, and when it is below the bound the modes and that its rows are left out
         INTEGER, INTENT(IN) :: mm
         REAL(8), INTENT(IN) :: om

         IF (mot_ok(mm)) THEN
            WRITE (*, '(a,F10.5,a,ES10.3)') ' QTF motion: Omega =', om, ' rad/s, rcond', mot_rcond(mm)
         ELSE
            WRITE (*, '(a)') ' QTF motion: '//TRIM(omega_left_text(mm))//'; its pairs are left out of the QTF_Motion files.'
         END IF
      END SUBROUTINE log_omega

   END SUBROUTINE qtf_motion_part

   FUNCTION free_modes(lhs) RESULT(free)
      ! The modes without restoring in a singular matrix of the motion: those whose row or column is zero to rounding
      ! (1e-12 of the largest entry), as surge, sway and yaw of C + K_ext for an unmoored body (R6)
      COMPLEX(8), INTENT(IN) :: lhs(:, :)
      LOGICAL :: free(SIZE(lhs, 1))
      REAL(8) :: bound
      INTEGER :: k

      bound = 1D-12*MAXVAL(ABS(lhs))
      DO k = 1, SIZE(lhs, 1)
         free(k) = ALL(ABS(lhs(k, :)) <= bound) .OR. ALL(ABS(lhs(:, k)) <= bound)
      END DO
   END FUNCTION free_modes

   FUNCTION omega_left_text(m) RESULT(text)
      ! "Omega = ... rad/s, rcond ... < 1e-12 (modes 1 2 6 without restoring)" for the Omega of index m (0: Omega = 0),
      ! for the log and the header of the motion files. The modes are I of the files (6 (body - 1) + mode).
      INTEGER, INTENT(IN) :: m
      CHARACTER(LEN=200) :: text
      CHARACTER(LEN=12) :: num
      CHARACTER(LEN=:), ALLOCATABLE :: modes
      REAL(8) :: om
      INTEGER :: k

      om = 0D0
      IF (m > 0) om = psi_need(m)
      modes = ''
      DO k = 1, SIZE(mot_free, 1)
         IF (.NOT. mot_free(k, m)) CYCLE
         WRITE (num, '(I0)') k
         modes = modes//' '//TRIM(num)
      END DO
      IF (LEN(modes) > 0) THEN
         modes = '(modes'//modes//' without restoring)'
      ELSE
         modes = '(no single mode without restoring)'
      END IF
      WRITE (text, '(a,F10.5,a,ES10.3,a)') 'Omega =', om, ' rad/s, rcond', mot_rcond(m), ' < 1e-12 '//modes
   END FUNCTION omega_left_text

   SUBROUTINE pair_force(mesh, geom, kernel, route, i, b1, j, b2, f)
      ! The force f(6 nbody) of the row (i, b1, j, b2) of the file of kernel and route, unscaled, e^{-i Omega t}: the
      ! quadratic part (qtf_pair_load) plus T_pot (route 1) or T_dir (route 2), or for the sum T+_pot or T+_dir. The same
      ! sums as write_file, element by element, so that the force is the one of the file to the bit; on the difference
      ! diagonal the potential part is left out, as write_sum does.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      INTEGER, INTENT(IN) :: kernel, route, i, b1, j, b2
      COMPLEX(8), INTENT(OUT) :: f(:)
      COMPLEX(8) :: total(6, mesh%nbody)
      INTEGER :: h

      CALL qtf_pair_load(mesh, geom, kernel, i, b1, j, b2, total)
      ! Mode k = 6 (body - 1) + c, as the rows of the files
      f = RESHAPE(total, [6*mesh%nbody])
      h = heading_pair(b1, b2)
      IF (kernel == KERNEL_DIFF) THEN
         IF (store(i, b1)%omega == store(j, b2)%omega) RETURN
         IF (route == 1) THEN
            f = f + tpot(:, i, j, h)
         ELSE
            f = f + tdir(:, i, j, h)
         END IF
      ELSE
         IF (route == 1) THEN
            f = f + tpot_sum(:, i, j, h)
         ELSE
            f = f + tdir_sum(:, i, j, h)
         END IF
      END IF
   END SUBROUTINE pair_force

   SUBROUTINE qtf_motion_load(mesh, geom, kernel, route, i, b1, j, b2, xi, f, ok)
      ! For the tests: xi2 (6 nbody) of the stored row (frequency i, heading b1; frequency j, heading b2) of kernel and
      ! route as qtf_motion_part left it, the force f it was solved for (pair_force) and ok, the regularity of its Omega.
      ! Call after qtf_motion_part and before qtf_finish, which frees the stores.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      INTEGER, INTENT(IN) :: kernel, route, i, b1, j, b2
      COMPLEX(8), INTENT(OUT) :: xi(:), f(:)
      LOGICAL, INTENT(OUT) :: ok
      INTEGER :: m

      IF (.NOT. mot_active .OR. route < 1 .OR. route > mot_nroute) THEN
         WRITE (*, '(a)') ' QTF motion: qtf_motion_load needs qtf_motion = 1 and a route of the run.'
         ERROR STOP 1
      END IF
      CALL pair_force(mesh, geom, kernel, route, i, b1, j, b2, f)
      ok = .FALSE.
      IF (kernel == KERNEL_DIFF) THEN
         xi = xi2(:, i, j, heading_pair(b1, b2), route)
         ! 0 is the diagonal; a pair beyond qtf_max_dw (difference_pair_index 0) is no row
         m = 0
         IF (store(i, 1)%omega /= store(j, 1)%omega) m = MERGE(-1, difference_pair_index(i, j), &
                                                              difference_pair_index(i, j) == 0)
         IF (m >= 0) ok = mot_ok(m)
      ELSE
         xi = xi2_sum(:, i, j, heading_pair(b1, b2), route)
         m = sum_pair_index(i, j)
         IF (m > 0) ok = mot_ok(m)
      END IF
   END SUBROUTINE qtf_motion_load

   SUBROUTINE out_of_memory(body, fs)
      ! body: the store of qtf_terms = body (psi, T_pot, the fields on C and F), whose size the log gave before; fs: also
      ! the free-surface store of qtf_terms = full
      LOGICAL, INTENT(IN), OPTIONAL :: body, fs

      WRITE (*, '(a)') ' QTF: not enough memory to keep the first-order states of every frequency. Use fewer frequencies,'
      WRITE (*, '(a)') ' headings or panels, or qtf = 0.'
      IF (PRESENT(body)) THEN
         IF (body) THEN
            WRITE (*, '(a)') ' With qtf_terms = body the store also holds psi at every difference frequency and, with qtf = 2'
            WRITE (*, '(a)') ' or 3, every sum frequency (count in the log above): a smaller qtf_max_dw or qtf_max_sw, or a'
            WRITE (*, '(a)') ' uniform omega grid (whose differences coincide), needs fewer.'
            WRITE (*, '(a)') ' With qtf_b2 = control the fields on C and F add to it: qtf_b2 = body avoids them, a larger'
            WRITE (*, '(a)') ' cs_panel_size reduces them.'
         END IF
      END IF
      IF (PRESENT(fs)) THEN
         IF (fs) THEN
            WRITE (*, '(a)') ' With qtf_terms = full the fields at the inner free-surface points and the partial waves add to it'
            WRITE (*, '(a)') ' (size in the log above): qtf_terms = body avoids them, a smaller qtf_fs_radius or qtf_fs_modes'
            WRITE (*, '(a)') ' reduces them.'
         END IF
      END IF
      ERROR STOP 1
   END SUBROUTINE out_of_memory

END MODULE QTFDriver
