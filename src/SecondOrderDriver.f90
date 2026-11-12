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
MODULE SecondOrderDriver
   ! Connects the second-order modules to the solver: copies the first-order mesh and solution of the single-body or
   ! multi-body arrays into the NBODY-independent types, and writes the results.
   !
   ! All states are per unit wave amplitude. FarFieldDrift needs this, and HydroStatic and HydroStaticMulti always set
   ! AMP = 1, so MXPOT and DSPL are already per unit amplitude. The set-up stops if that ever changes.
   !
   ! PanelDiagnostics.txt (write_terms = 2) is test and debugging output, in tagged rows of plain columns. See
   ! docs/output_files.md. QTF_PsiCheck.txt (qtf_terms = body, write_terms >= 1) is test output of the same kind: the added
   ! mass and damping and the excitation by the indirect method, from psi (header in the file).
   !
   ! The control-surface method (drift_methods = control) needs the states of all the headings of a frequency at once: the
   ! fields on C are computed in one pass over the points at the end of the frequency, for every heading, since
   ! the rows of all points are not stored (FieldVelocity.f90). So with it the states of every heading are
   ! kept, as for mean_drift = 2.
   !
   ! qtf_b2 = control (with qtf_terms = body) moves the surface piece of B2 to C and F (QTFPotential). It builds C and F as
   ! the control-surface method does, with the same keys, also when that method is off, and needs F and a bottom of C for
   ! every body. The fields there are evaluated once for both: the total velocity of every heading at the end of each
   ! frequency, and grad psi_k (fields_green of one state per mode, no incident wave) when psi of a difference frequency is
   ! taken, while the solver globals are at that frequency. Both go to the QTF store.
   !
   ! Its F_z, M_x and M_y (Phase B) also need the free surface F inside C and the hull waterline potential. They are written for
   ! each body whose C has a bottom, when F could be built. The hull waterline potential then comes from waterline_method, with
   ! the hull velocity of velocity_method when it is extrapolated, as for the near-field method.
   !
   ! qtf_terms = full (Phase 6a) adds the free-surface integral of the second-order potential, split at a circle r = R about
   ! all bodies. At set-up: the inner free-surface mesh (build_fs_inner, with the largest frequency of the grid), the offsets
   ! of its points (fs_inner_steps), the mode counts M and N_ev, and the samples on the circle (pw_sample_points). At the
   ! end of each frequency, after the field points of qtf_field: the first-order fields of every heading at the inner
   ! points (fs_inner_fields, 5-point stencil, stock route) and the partial waves of the scattered potential outside R
   ! (fs_circle_modes, the samples by the centroid route), with the mode check. When psi of a difference frequency is
   ! taken: psi_k at the inner points (value only) and its partial waves, as keep_cs_psi does for C and F. All go to the
   ! QTF store, whose size the log gives before qtf_init allocates it.

   USE Const_mod, ONLY: G, RHO
   USE HAMS_mod, ONLY: SYBO
   USE WaveDyn_mod, ONLY: AMP, BETA, BLNR, H, NBETA, NBODY, NK, NPER, OUFR, TP, V, W1, WK, WVFQ, WVHD, WVNB, AMAS, BDMP, EXFC, &
                          DSPL, AMAS_MULTI_COMB, BDMP_MULTI_COMB, BLNR_MULTI, EXFC_MULTI_COMB, DSPL_MULTI_COMB
   USE Body_mod, ONLY: CRS, KSTF, MATX, REFL, XR, XR_MULTI, CRS_MULTI, KSTF_MULTI, MATX_MULTI
   USE PanelMesh_mod, ONLY: DS, DXYZ_P, GEONRM_P, GEONRM_MULTI_COMB, ISX, ISY, NCN, NCON, NELEM, NSYS, PNSZ, XYZ, XYZ_P, &
                            DS_MULTI_COMB, DXYZ_MULTI_COMB, NCN_MULTI_COMB, NCON_MULTI_COMB, NELEM_MULTI, NTND_MULTI, &
                            PNSZ_MULTI_COMB, XYZ_GLOBAL_MULTI_COMB, XYZ_GLOBAL_MULTI_COMB_P
   USE Inerfs_mod, ONLY: NELEM_PE
   USE FieldOutput_mod, ONLY: NFP, XFP
   USE Potentials_mod, ONLY: MXPOT, MXPOT_MULTI_COMB
   USE IO, ONLY: CreateDirectory
   USE ImplementSubs, ONLY: InputOmega
   USE SecondOrderTypes, ONLY: so_options, panel_mesh, first_order_state, surface_geom, cs_points, cs_fields, CS_WALL, CS_LINE, &
                               CS_FREE, QTF_TERMS_BODY, QTF_TERMS_FULL, QTF_B2_CONTROL
   USE FirstOrderState, ONLY: expand_images, build_state, radiation_psi
   USE SurfaceGeometry, ONLY: build_surface_geom, DEFAULT_SHARP_DEG
   USE SurfaceVelocity, ONLY: velocity_differencing, source_operator, build_source_operator, velocity_source
   USE FarFieldDrift, ONLY: far_amplitude, far_field_pair, far_field_pairs
   USE NearFieldDrift, ONLY: waterline_extrapolate, near_field_pair, T_WATERLINE, T_HYDROSTATIC
   USE WaterlineField, ONLY: waterline_rows, waterline_field, field_point_inside_body
   USE ControlSurfaceMesh, ONLY: build_control_surface, check_control_surface, write_control_surface, build_fs_inner
   USE FieldVelocity, ONLY: fields_green, fields_source
   USE ControlSurfaceDrift, ONLY: control_surface_pair, CS_T_LINE, CS_T_PRESSURE, CS_T_HYDROSTATIC
   USE SecondOrderOutput, ONLY: write_drift_row, write_terms_row
   USE QTFDriver, ONLY: qtf_init, qtf_keep, qtf_finish, qtf_stored_mb, qtf_keep_psi, qtf_psi_index, qtf_psi_count, &
                        qtf_radiation_list, qtf_psi_tolerance, &
                        qtf_radiation_omegas, qtf_set_control_surface, qtf_keep_cs_velocity, qtf_keep_cs_psi, &
                        qtf_set_edge_weights, qtf_keep_fs_state, qtf_keep_fs_psi
   USE QTFFreeSurface, ONLY: fs_point_field, fs_inner_fields, fs_inner_steps, fs_sample_phi, fs_circle_modes, fs_mode_warning, &
                             fs_auto_modes, fs_auto_evanescent, fs_radius_message
   USE PartialWaves, ONLY: pw_set, pw_sample_points, pw_evanescent_from_wvn, PW_KR_MAX
   USE QTFPotential, ONLY: b2_control_check
   USE SharpEdge, ONLY: sharp_edge_weights
   USE QTFField, ONLY: build_field_points, qfield_init, qfield_keep, qfield_finish, qfield_stored_mb, FIELD_FD_STEP
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: so_setup_single, so_setup_multi, so_after_motion_single, so_after_motion_multi, so_end_frequency, so_finish
   PUBLIC :: so_after_radiation, so_after_radiation_only

   REAL(8), PARAMETER :: PI = 3.141592653589793D0
   COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
   INTEGER, PARAMETER :: DIAG_N_ANGLES = 72            ! angles of the far-field amplitudes in the diagnostics

   TYPE(panel_mesh) :: mesh
   TYPE(surface_geom) :: geom
   TYPE(first_order_state), ALLOCATABLE :: states(:)   ! one per heading for mean_drift = 2 or control, otherwise the current
   COMPLEX(8), ALLOCATABLE :: phi_wl(:,:)              ! (nwl, size(states)), total potential on the waterline segments; allocated
                                                       ! for the near-field method and for the Phase B control-surface loads
   COMPLEX(8), ALLOCATABLE :: rows_g(:,:,:), rows_dg(:,:,:)   ! (nbase, nsys, nwl), of the current frequency, for waterline_method = field
   LOGICAL :: multi_arrays = .FALSE.                   ! the solver's mesh is in the combined multi-body arrays
   INTEGER :: unit_far = 0, unit_near = 0, unit_terms = 0, unit_diag = 0   ! 0: not open (NEWUNIT gives negative numbers)
   INTEGER :: n_velocity_fields = 0, n_fallback_total = 0                  ! over the run, logged by so_finish
   TYPE(source_operator) :: source_op                  ! of the current frequency, for velocity_method = source
   COMPLEX(8), ALLOCATABLE :: phi_source(:)            ! (n), potential of the sources of the current state, for the diagnostics
   INTEGER :: n_source_builds = 0, n_source_solves = 0 ! over the run, logged by so_finish with the two times
   REAL(8) :: source_build_seconds = 0D0, source_solve_seconds = 0D0
   LOGICAL :: hull_uses_source = .FALSE.               ! velocity_method = source for the hull velocity (phi_wl allocated)
   ! Control-surface method
   TYPE(cs_points) :: cs_pts                           ! points of C, Gamma_c and (for Phase B) F of all bodies
   LOGICAL :: have_free = .FALSE.                      ! F was built for a body: the Phase B loads are computed
   LOGICAL, ALLOCATABLE :: cs_vertical(:)              ! (nbody), F_z, M_x and M_y are written: F was built and C has a bottom
   COMPLEX(8), ALLOCATABLE :: cs_sigma(:,:)            ! (n, NBETA), source strengths of each heading, for cs_velocity = source
   COMPLEX(8), ALLOCATABLE :: far_diag(:,:)            ! (3, NBETA), far-field F_x, F_y, M_z of each heading, for the log
   INTEGER :: unit_control = 0
   ! Radiation potentials of the indirect method (qtf_terms = body)
   LOGICAL :: body_terms = .FALSE.                     ! psi is kept at the difference frequencies
   COMPLEX(8), ALLOCATABLE :: psi_main(:,:)            ! (n, 6*nbody), psi of the current main-grid frequency
   INTEGER :: unit_psi = 0                             ! QTF_PsiCheck.txt (write_terms >= 1)
   LOGICAL :: b2_control = .FALSE.                     ! qtf_b2 = control: fields on C and F go to the QTF store
   TYPE(cs_points) :: b2_pts                           ! the points of C and F of cs_pts (no Gamma_c), for B2
   INTEGER, ALLOCATABLE :: b2_idx(:)                   ! (b2_pts%n), their index in cs_pts
   ! Partial second-order field at the field points of ControlFile.in (qtf_field = 1)
   TYPE(cs_points) :: fp_pts                           ! the field points, kind CS_FREE
   ! Free-surface integral of the difference QTF (qtf_terms = full)
   LOGICAL :: fs_terms = .FALSE.                       ! the fields on the inner free surface and the partial waves are kept
   TYPE(cs_points) :: fs_pts                           ! the inner free-surface points (build_fs_inner)
   TYPE(cs_points) :: fs_circle                        ! the samples on the circle r = R (pw_sample_points)
   REAL(8) :: fs_centre(2) = 0D0, fs_radius = 0D0      ! the circle
   REAL(8) :: fs_step = 0D0                            ! offset of the inner points as a fraction of their size (cs_fd_step)
   INTEGER :: fs_mmax = 0, fs_nev = 0                  ! angular modes M and evanescent modes N_ev outside R
   REAL(8), ALLOCATABLE :: fs_delta(:)                 ! (fs_pts%n), offsets of fs_inner_fields, smaller next to a waterline

CONTAINS

   SUBROUTINE so_setup_single(opts, inputdir, outputdir)
      ! Call after CalNormals.
      TYPE(so_options), INTENT(IN) :: opts
      CHARACTER(LEN=*), INTENT(IN) :: inputdir, outputdir
      INTEGER, ALLOCATABLE :: body(:)

      ALLOCATE (body(NELEM))
      body = 1
      CALL expand_images(XYZ_P, GEONRM_P, DXYZ_P, DS, PNSZ, XYZ, NCN, NCON, body, RESHAPE(XR, [3, 1]), ISX, ISY, NSYS, mesh)
      multi_arrays = .FALSE.
      CALL common_setup(opts, inputdir, outputdir)
   END SUBROUTINE so_setup_single

   SUBROUTINE so_setup_multi(opts, inputdir, outputdir)
      ! Call after BuildCombinedArrays. MXPOT_MULTI_COMB has NELEM_PE rows: with IRSP = 1 the global TNELEM also counts
      ! the waterplane panels. Node numbers in NCON_MULTI_COMB already index XYZ_GLOBAL_MULTI_COMB.
      TYPE(so_options), INTENT(IN) :: opts
      CHARACTER(LEN=*), INTENT(IN) :: inputdir, outputdir
      INTEGER, ALLOCATABLE :: body(:)
      INTEGER :: b, first

      ALLOCATE (body(NELEM_PE))
      first = 1
      DO b = 1, NBODY
         body(first:first + NELEM_MULTI(b) - 1) = b
         first = first + NELEM_MULTI(b)
      END DO
      CALL expand_images(XYZ_GLOBAL_MULTI_COMB_P(1:NELEM_PE, :), GEONRM_MULTI_COMB(1:NELEM_PE, :), DXYZ_MULTI_COMB(1:NELEM_PE, :), &
                         DS_MULTI_COMB(1:NELEM_PE), PNSZ_MULTI_COMB(1:NELEM_PE), &
                         XYZ_GLOBAL_MULTI_COMB(1:SUM(NTND_MULTI), :), NCN_MULTI_COMB(1:NELEM_PE), &
                         NCON_MULTI_COMB(1:NELEM_PE, :), body, TRANSPOSE(XR_MULTI), ISX, ISY, NSYS, mesh)
      multi_arrays = .TRUE.
      CALL common_setup(opts, inputdir, outputdir)
   END SUBROUTINE so_setup_multi

   SUBROUTINE common_setup(opts, inputdir, outputdir)
      TYPE(so_options), INTENT(IN) :: opts
      CHARACTER(LEN=*), INTENT(IN) :: inputdir, outputdir
      CHARACTER(LEN=:), ALLOCATABLE :: msg
      LOGICAL :: ok
      INTEGER :: nqtf

      IF (ABS(AMP - 1D0) > 1D-12) THEN
         WRITE (*, *) 'Mean drift needs first-order results per unit wave amplitude, but AMP = ', AMP
         ERROR STOP 1
      END IF

      b2_control = uses_b2_control(opts)
      IF (opts%qtf_b2_given .AND. .NOT. (opts%qtf_terms >= QTF_TERMS_BODY .AND. (opts%qtf == 1 .OR. opts%qtf == 3))) THEN
         WRITE (*, '(a)') ' QTF note: qtf_b2 is ignored; it applies with qtf_terms = body or full and the difference QTF only.'
      END IF
      IF (opts%qtf_edge_given .AND. .NOT. b2_control) THEN
         WRITE (*, '(a)') ' QTF note: qtf_edge_weights and qtf_edge_angle are ignored; they apply with qtf_terms = body or full,'// &
            ' the difference QTF and qtf_b2 = control only.'
      END IF
      IF (opts%qtf_fs_given .AND. .NOT. (opts%qtf_terms == QTF_TERMS_FULL .AND. (opts%qtf == 1 .OR. opts%qtf == 3))) THEN
         WRITE (*, '(a)') ' QTF note: qtf_fs_radius, qtf_fs_rows, qtf_fs_modes and qtf_fs_evanescent are ignored; they apply'// &
            ' with qtf_terms = full and the difference QTF only.'
      END IF
      ALLOCATE (states(MERGE(NBETA, 1, keep_all_headings(opts))))
      ! Only the near-field and control-surface methods use the surface geometry (the control surface takes the waterplane
      ! extent from it), so only they can stop on the geometry's errors
      IF (opts%do_near .OR. opts%do_control) THEN
         CALL build_surface_geom(mesh, opts%wl_z_tolerance, DEFAULT_SHARP_DEG, geom, ok, msg)
         IF (.NOT. ok) THEN
            WRITE (*, *) 'Mean drift: ', msg
            ERROR STOP 1
         END IF
         IF (LEN_TRIM(msg) > 0) WRITE (*, *) 'Mean drift warning: ', msg
      END IF
      IF (opts%do_control) CALL control_setup(opts, inputdir, outputdir)
      IF (b2_control) CALL b2_control_setup(opts, inputdir, outputdir)
      ! The hull waterline potential, for the near-field method and the Phase B control-surface loads
      IF (opts%do_near .OR. have_free) THEN
         hull_uses_source = opts%velocity_method == 2
         ALLOCATE (phi_wl(geom%nwl, SIZE(states)))
         IF (opts%velocity_method == 2 .AND. opts%write_terms == 2) ALLOCATE (phi_source(mesh%n))
         IF (opts%waterline_method == 2) THEN
            IF (ANY(field_point_inside_body(geom%wl_nz))) THEN
               WRITE (*, '(a)') 'Mean drift: waterline_method = field puts its point inside the body where'
               WRITE (*, '(a)') 'the wall at the waterline widens downward by more than 45 degrees from the'
               WRITE (*, '(a)') 'vertical. Use waterline_method = extrapolate for this body.'
               ERROR STOP 1
            END IF
            ALLOCATE (rows_g(mesh%nbase, mesh%nsys, geom%nwl), rows_dg(mesh%nbase, mesh%nsys, geom%nwl))
            IF (H >= 0D0) THEN
               WRITE (*, '(a)') 'Mean drift warning: with waterline_method = field in finite depth, term I'
               WRITE (*, '(a)') 'has been seen to jump with the mesh or the offset. That was measured before'
               WRITE (*, '(a)') 'the Green function used Linton''s method at small R/h, and not since. Check'
               WRITE (*, '(a)') 'against extrapolate. See Limits under SecondOrder.in in docs/input_files.md.'
            END IF
         END IF
      END IF
      IF (opts%do_near) THEN
         OPEN (NEWUNIT=unit_near, FILE=outputdir//'/Wamit_format/MeanDriftPressure.9', STATUS='REPLACE', ACTION='WRITE')
      END IF
      IF (opts%qtf > 0) THEN
         ! validate_second_order makes sure the near-field method, and so phi_wl, is there
         ! The zero and infinite frequencies of SYBO = 1 are not kept
         nqtf = NPER - MERGE(2, 0, SYBO == 1)
         WRITE (*, '(a,F10.1,a,I0,a,I0,a)') ' QTF: the first-order states take', &
            qtf_stored_mb(mesh%n, geom%nwl, mesh%nbody, nqtf, NBETA), ' MB for ', nqtf, ' frequencies and ', NBETA, ' headings'
         body_terms = opts%qtf_terms >= QTF_TERMS_BODY .AND. (opts%qtf == 1 .OR. opts%qtf == 3)
         IF (body_terms) THEN
            CALL body_terms_setup(opts, nqtf, inputdir, outputdir)
         ELSE
            CALL qtf_init(opts, nqtf, WVHD(1:NBETA), geom%nwl, mesh%n, mesh%nbody, outputdir)
         END IF
      END IF
      IF (opts%qtf_field == 1) CALL field_setup(opts, nqtf, outputdir)
      IF ((opts%do_near .OR. opts%do_control) .AND. opts%write_terms >= 1) THEN
         OPEN (NEWUNIT=unit_terms, FILE=outputdir//'/Wamit_format/MeanDriftTerms.txt', STATUS='REPLACE', ACTION='WRITE')
         WRITE (unit_terms, '(a)') '# PER BETA1 BETA2 body I term Re Im, scaled and conjugated as MeanDriftPressure.9'
         IF (opts%do_near) THEN
            WRITE (unit_terms, '(a)') '# term: 1 waterline, 2 flare, 3 velocity, 4 displacement, 5 rotation, 6 hydrostatic'
         END IF
         IF (opts%do_control) THEN
            WRITE (unit_terms, '(a)') '# term: 11 control-surface line Gamma_c, 12 flux through C, 13 pressure on C (I = 1 to 6)'
            IF (have_free) THEN
               WRITE (unit_terms, '(a)') '# term: 14 hull waterline Gamma, 15 free surface F, 16 hydrostatic (I = 3, 4, 5)'
            END IF
         END IF
      END IF
      IF (opts%do_far) THEN
         OPEN (NEWUNIT=unit_far, FILE=outputdir//'/Wamit_format/MeanDriftMomentum.8', STATUS='REPLACE', ACTION='WRITE')
      END IF
      IF (opts%write_terms == 2) THEN
         CALL CreateDirectory(outputdir//'/SecondOrder', ok)
         IF (.NOT. ok) THEN
            WRITE (*, *) 'Mean drift: could not create the directory ', outputdir//'/SecondOrder'
            ERROR STOP 1
         END IF
         OPEN (NEWUNIT=unit_diag, FILE=outputdir//'/SecondOrder/PanelDiagnostics.txt', STATUS='REPLACE', ACTION='WRITE')
         WRITE (unit_diag, '(a)') '# Test and debugging output of the mean drift calculation. Time factor e^{-i omega t}.'
         WRITE (unit_diag, '(a)') '# WAVE omega k nu h'
         WRITE (unit_diag, '(a)') '# MODE omega j theta_deg Re(a_j) Im(a_j)      amplitude radiated by unit motion in mode j'
         WRITE (unit_diag, '(a)') '# RESIDUAL omega beta_deg residual            of the equation of motion, relative to |F|'
         WRITE (unit_diag, '(a)') '# FORCE omega beta_deg j Re(F_j) Im(F_j)      excitation force, F = i*EXFC'
         WRITE (unit_diag, '(a)') '# TOTAL omega beta_deg theta_deg Re(a) Im(a)  amplitude of the perturbation wave'
         WRITE (unit_diag, '(a)') '# PANEL omega beta_deg panel x y z Re(phi_D) Im(phi_D) Re(phi_total) Im(phi_total)'
         IF (ALLOCATED(phi_source)) THEN
            WRITE (unit_diag, '(a)') '# SOURCE omega beta_deg panel phi_P phi_source dphi_P_dn grad_P.n grad_P(1:3)   each Re Im'
         END IF
      END IF
   END SUBROUTINE common_setup

   SUBROUTINE field_setup(opts, nqtf, outputdir)
      ! The field points of ControlFile.in as cs_points, and the store of the partial second-order field. A point above the
      ! free surface, on a hull or inside a body stops the run, with the point and the body in the message. Call after the
      ! QTF set-up, with nqtf the number of frequencies of the QTF grid.
      TYPE(so_options), INTENT(IN) :: opts
      INTEGER, INTENT(IN) :: nqtf
      CHARACTER(LEN=*), INTENT(IN) :: outputdir
      CHARACTER(LEN=:), ALLOCATABLE :: msg
      LOGICAL, ALLOCATABLE :: surface(:)
      LOGICAL :: ok

      IF (NFP < 1) THEN
         WRITE (*, '(a)') ' qtf_field = 1 needs field points: Number_of_field_points in ControlFile.in'
         ERROR STOP 1
      END IF
      CALL build_field_points(mesh, TRANSPOSE(XFP), fp_pts, surface, ok, msg)
      IF (.NOT. ok) THEN
         WRITE (*, '(a,a)') ' QTF field points: ', msg
         ERROR STOP 1
      END IF
      WRITE (*, '(a,I0,a,I0,a,I0,a,F10.1,a)') ' QTF field points: ', COUNT(.NOT. surface), ' pressure points and ', &
         COUNT(surface), ' elevation points, ', nqtf, ' frequencies; the store takes', &
         qfield_stored_mb(fp_pts%n, nqtf, NBETA), ' MB'
      CALL qfield_init(opts, nqtf, WVHD(1:NBETA), fp_pts, surface, outputdir)
   END SUBROUTINE field_setup

   SUBROUTINE body_terms_setup(opts, nqtf, inputdir, outputdir)
      ! The QTF store with psi at the difference frequencies of the main grid, and QTF_PsiCheck.txt (write_terms >= 1).
      ! The main-grid frequencies come from the input values (WVNB, before the frequency loop converts them). With
      ! qtf_terms = full also the inner free surface and the circle (fs_setup), and their share of the store.
      TYPE(so_options), INTENT(IN) :: opts
      INTEGER, INTENT(IN) :: nqtf
      CHARACTER(LEN=*), INTENT(IN) :: inputdir, outputdir
      REAL(8) :: omegas(nqtf), mb, mb_cs, mb_fs
      REAL(8), ALLOCATABLE :: need(:)
      INTEGER, ALLOCATABLE :: reuse(:)
      INTEGER :: i, first, nneed, nreuse, npairs

      first = MERGE(2, 0, SYBO == 1)
      DO i = 1, nqtf
         omegas(i) = InputOmega(WVNB(first + i))
      END DO
      ! The counts and the size go to the log before qtf_init allocates, so that a run that does not fit has them; qtf_init
      ! builds the same list
      CALL qtf_radiation_list(omegas, opts%qtf_max_dw, qtf_psi_tolerance(omegas), need, reuse)
      nneed = SIZE(need)
      nreuse = COUNT(reuse > 0)
      WRITE (*, '(a,I0,a,I0,a,I0,a)') ' QTF body forcing: psi at ', nneed, ' difference frequencies, ', nreuse, &
         ' from the main grid and ', nneed - nreuse, ' radiation-only solves (each costs about one radiation solve)'
      npairs = MERGE(NBETA*NBETA, NBETA, opts%mean_drift == 2)
      mb = qtf_stored_mb(mesh%n, geom%nwl, mesh%nbody, nqtf, NBETA, npsi=nneed, npairs=npairs)
      fs_terms = uses_fs(opts)
      IF (fs_terms) THEN
         CALL fs_setup(opts, inputdir, omegas)
         mb_fs = qtf_stored_mb(mesh%n, geom%nwl, mesh%nbody, nqtf, NBETA, npsi=nneed, npairs=npairs, nfs=fs_pts%n, &
                               npw=(2*fs_mmax + 1)*(fs_nev + 1)) - mb
         WRITE (*, '(a,F10.1,a,I0,a)') ' QTF free surface: the store takes', mb_fs, ' MB more for the fields at the ', &
            fs_pts%n, ' inner points and the partial waves (qtf_terms = full)'
      END IF
      IF (b2_control) THEN
         ! The share of C and F is what the points add to the store, so qtf_stored_mb is the one formula
         mb_cs = qtf_stored_mb(mesh%n, geom%nwl, mesh%nbody, nqtf, NBETA, npsi=nneed, npairs=npairs, ncs=b2_pts%n) - mb
         WRITE (*, '(a,F10.1,a,F10.1,a,I0,a)') ' QTF: with psi and the potential part the store takes', mb + mb_cs, &
            ' MB, of which', mb_cs, ' MB for the fields at the ', b2_pts%n, ' points of C and F (qtf_b2 = control)'
      ELSE
         WRITE (*, '(a,F10.1,a)') ' QTF: with psi and the potential part the store takes', mb, ' MB'
      END IF
      IF (fs_terms) THEN
         CALL qtf_init(opts, nqtf, WVHD(1:NBETA), geom%nwl, mesh%n, mesh%nbody, outputdir, omegas_main=omegas, fs_pts=fs_pts, &
                       fs_centre=fs_centre, fs_radius=fs_radius, fs_depth=H, fs_mmax=fs_mmax, fs_nev=fs_nev)
      ELSE
         CALL qtf_init(opts, nqtf, WVHD(1:NBETA), geom%nwl, mesh%n, mesh%nbody, outputdir, omegas_main=omegas)
      END IF
      IF (qtf_psi_count() /= nneed .OR. qtf_psi_count(reused=.TRUE.) /= nreuse) THEN
         WRITE (*, '(a)') ' QTF: qtf_init found other difference frequencies than the log above.'
         ERROR STOP 1
      END IF
      IF (b2_control) THEN
         CALL qtf_set_control_surface(b2_pts)
         IF (opts%qtf_edge_weights) CALL b2_edge_weights(opts)
      END IF
      IF (fs_terms) THEN
         WRITE (*, '(a)') ' QTF_Indirect.11d: body forcing only for now; the free-surface fields are kept (qtf_terms = full),'// &
            ' their integral is not yet added.'
      ELSE
         WRITE (*, '(a)') ' QTF_Indirect.11d: body forcing only, the free-surface integral is omitted.'
      END IF
      ALLOCATE (psi_main(mesh%n, 6*mesh%nbody))
      psi_main = (0D0, 0D0)
      IF (opts%write_terms >= 1) THEN
         OPEN (NEWUNIT=unit_psi, FILE=outputdir//'/Wamit_format/QTF_PsiCheck.txt', STATUS='REPLACE', ACTION='WRITE')
         WRITE (unit_psi, '(a)') '# Checks of the radiation potentials psi_k of the indirect method (qtf_terms = body).'
         WRITE (unit_psi, '(a)') '# SI units, time factor e^{-i omega t}, per unit wave amplitude, moments about XR of each body.'
         WRITE (unit_psi, '(a)') '# PSI omega source k l Re Im   rho int_{S_l} psi_k n_l dS = A_lk + i B_lk/omega;'// &
            ' source 1 main grid, 2 radiation-only solve'
         WRITE (unit_psi, '(a)') '# EXC omega beta_deg k Re Im   i omega rho [int_{S_k} phi_I n_k dS - int_S psi_k dphi_I/dn dS]'
      END IF
   END SUBROUTINE body_terms_setup

   SUBROUTINE fs_setup(opts, inputdir, omegas)
      ! qtf_terms = full: the inner free surface (build_fs_inner, panel size from the largest frequency of the grid omegas),
      ! the offsets of its points (fs_inner_steps, once: they do not depend on the frequency), the mode counts of the
      ! partial waves (M = qtf_fs_modes, or ceil(k0_max R) + 10 with k0_max the wave number of the largest frequency;
      ! N_ev = qtf_fs_evanescent, or max(8, ceil(9 h / (pi (R - r_wl)))) with r_wl the largest distance from the centre to
      ! a waterline node, fs_auto_evanescent) and the samples on the circle (pw_sample_points). Stops if the mesh cannot be
      ! built or the evanescent modes do not fit at R.
      TYPE(so_options), INTENT(IN) :: opts
      CHARACTER(LEN=*), INTENT(IN) :: inputdir
      REAL(8), INTENT(IN) :: omegas(:)
      CHARACTER(LEN=:), ALLOCATABLE :: msg
      REAL(8) :: omega_max, wvn1(1), k0_max, r_wl
      LOGICAL :: ok
      INTEGER :: nshrunk, nfloor, iw
      EXTERNAL :: DISPERSION

      omega_max = MAXVAL(omegas)
      CALL build_fs_inner(mesh, geom, opts, H, inputdir, fs_pts, fs_centre, fs_radius, ok, msg, omega_max=omega_max)
      IF (.NOT. ok) THEN
         WRITE (*, '(a,a)') ' QTF free surface: ', msg
         ERROR STOP 1
      END IF
      CALL DISPERSION(wvn1, 1, omega_max, H)
      k0_max = wvn1(1)
      fs_mmax = fs_auto_modes(opts%qtf_fs_modes, k0_max, fs_radius)
      ! r_wl: the largest distance from the centre to a waterline node, of all bodies (as build_fs_inner)
      r_wl = 0D0
      DO iw = 1, geom%nwl
         r_wl = MAX(r_wl, NORM2(geom%wl_x1(1:2, iw) - fs_centre), NORM2(geom%wl_x2(1:2, iw) - fs_centre))
      END DO
      ! R <= r_wl (the equality too) cuts a body; build_fs_inner refuses it for a given R, this holds for the file mesh too
      msg = fs_radius_message(fs_radius, r_wl)
      IF (LEN(msg) > 0) THEN
         WRITE (*, '(a,a)') ' QTF free surface: ', msg
         ERROR STOP 1
      END IF
      fs_nev = fs_auto_evanescent(opts%qtf_fs_evanescent, H, fs_radius, r_wl)
      IF (fs_nev > NK - 1) THEN
         WRITE (*, '(a,I0,a,I0,a)') ' QTF free surface: qtf_fs_evanescent = ', fs_nev, ' is more than the ', NK - 1, &
            ' evanescent wave numbers of the solver (WVN); lower it.'
         ERROR STOP 1
      END IF
      ! k_n < n pi / h at every frequency (DISPERSION): K_m(k_n R) stays clear of underflow (pw_decompose) below PW_KR_MAX
      IF (fs_nev*PI/H*fs_radius > PW_KR_MAX) THEN
         WRITE (*, '(a,I0,a,F10.2,a,F10.2,a,F6.0,a)') ' QTF free surface: qtf_fs_evanescent = ', fs_nev, ' at R =', fs_radius, &
            ' m in depth', H, ' m puts k_n R above ', PW_KR_MAX, ', where K_m(k_n R) underflows; lower qtf_fs_evanescent'// &
            ' or qtf_fs_radius.'
         ERROR STOP 1
      END IF
      CALL pw_sample_points(fs_centre, fs_radius, H, fs_mmax, fs_nev, fs_circle)
      ! Twice cs_fd_step: the 5-point Laplacian over delta is the difference of central differences over delta/2 (the
      ! 25-row stencil of 9fdd193 at cs_fd_step) to rounding, and a smaller delta lets the rounding of the Green function,
      ! divided by delta^2, take over (fs_laplace_residual: 1.6e-2 at delta = 0.02 psize against 7.2e-4 at 0.1)
      fs_step = 2D0*opts%cs_fd_step
      ALLOCATE (fs_delta(fs_pts%n))
      CALL fs_inner_steps(mesh, fs_pts, fs_step, fs_delta, nshrunk, nfloor)
      WRITE (*, '(a,I0,a,F10.3,a,F10.3,a,F10.3,a)') ' QTF free surface: ', fs_pts%n, ' inner points, circle R =', fs_radius, &
         ' m about (', fs_centre(1), ',', fs_centre(2), ')'
      WRITE (*, '(a,I0,a,I0,a,a,I0,a,F8.4,a,F8.3,a)') ' QTF free surface: partial waves outside R with M = ', fs_mmax, &
         ' and N_ev = ', fs_nev, TRIM(MERGE(' (automatic)', ' (given)    ', opts%qtf_fs_evanescent == 0)), &
         ' evanescent modes, from ', fs_circle%n, ' samples on the circle (k0_max R =', k0_max*fs_radius, ', r_wl =', r_wl, ' m)'
      IF (nshrunk > 0) THEN
         WRITE (*, '(a,I0,a)') ' QTF free surface: ', nshrunk, ' inner points next to a waterline take a smaller offset for'// &
            ' the horizontal Laplacian (the stencil must stay outside the waterplanes, about 16 offsets from the waterline).'
      END IF
      IF (nfloor > 0) THEN
         WRITE (*, '(a,I0,a)') ' QTF free surface warning: ', nfloor, ' inner points are so close to a waterline that their'// &
            ' offset is held at 1e-3 of step x size (their stencil is clear of the waterplane, but closer to it than 16'// &
            ' offsets): check the inner mesh there.'
      END IF
   END SUBROUTINE fs_setup

   SUBROUTINE b2_edge_weights(opts)
      ! The panel weights of the B2 hull term at the sharp edges of each body (SharpEdge), into the QTF store
      TYPE(so_options), INTENT(IN) :: opts
      REAL(8), ALLOCATABLE :: w(:)
      INTEGER, ALLOCATABLE :: nedge(:), nseg(:), nweighted(:)
      INTEGER :: b

      ALLOCATE (w(mesh%n), nedge(mesh%nbody), nseg(mesh%nbody), nweighted(mesh%nbody))
      CALL sharp_edge_weights(mesh, opts%qtf_edge_angle, w, nedge, nseg, nweighted)
      DO b = 1, mesh%nbody
         WRITE (*, '(a,I0,a,I0,a,I0,a,I0,a,F6.1,a)') ' QTF B2 edge weights: body ', b, ', ', nedge(b), &
            ' sharp edges (', nseg(b), ' panel edges), ', nweighted(b), ' panels weighted (qtf_edge_angle', &
            opts%qtf_edge_angle, ')'
      END DO
      CALL qtf_set_edge_weights(w)
   END SUBROUTINE b2_edge_weights

   SUBROUTINE write_psi_rows(omega, source, psi)
      ! PSI rows of QTF_PsiCheck.txt: rho int psi_k n_l dS over the panels of the body of mode l, for all k and l
      REAL(8), INTENT(IN) :: omega
      INTEGER, INTENT(IN) :: source
      COMPLEX(8), INTENT(IN) :: psi(:, :)
      COMPLEX(8) :: s
      INTEGER :: k, l, lb, i

      IF (unit_psi == 0) RETURN
      DO k = 1, 6*mesh%nbody
         DO l = 1, 6*mesh%nbody
            lb = (l - 1)/6 + 1
            s = (0D0, 0D0)
            DO i = 1, mesh%n
               IF (mesh%body(i) == lb) s = s + psi(i, k)*mesh%nj(l - 6*(lb - 1), i)*mesh%area(i)
            END DO
            WRITE (unit_psi, '(a,ES24.15,I3,2I6,2ES24.15)') 'PSI ', omega, source, k, l, RHO*REAL(s, 8), RHO*AIMAG(s)
         END DO
      END DO
   END SUBROUTINE write_psi_rows

   SUBROUTINE write_exc_rows(s, beta_deg)
      ! EXC rows of QTF_PsiCheck.txt: the excitation by the indirect method (Haskind), with psi of the main-grid frequency
      TYPE(first_order_state), INTENT(IN) :: s
      REAL(8), INTENT(IN) :: beta_deg
      COMPLEX(8) :: incident, scattered, f
      INTEGER :: k, kb, i

      IF (unit_psi == 0) RETURN
      DO k = 1, 6*mesh%nbody
         kb = (k - 1)/6 + 1
         incident = (0D0, 0D0)
         scattered = (0D0, 0D0)
         DO i = 1, mesh%n
            IF (mesh%body(i) == kb) incident = incident + s%phi_i(i)*mesh%nj(k - 6*(kb - 1), i)*mesh%area(i)
            scattered = scattered + psi_main(i, k)*SUM(s%grad_i(:, i)*mesh%nrm(:, i))*mesh%area(i)
         END DO
         f = CI*s%omega*RHO*(incident - scattered)
         WRITE (unit_psi, '(a,2ES24.15,I6,2ES24.15)') 'EXC ', s%omega, beta_deg, k, REAL(f, 8), AIMAG(f)
      END DO
   END SUBROUTINE write_exc_rows

   SUBROUTINE so_after_radiation(opts, kk)
      ! Call after the radiation solve of main-grid frequency kk, before the diffraction solves overwrite the scattered
      ! column of MXPOT. With qtf_terms = body, takes psi of this frequency: for the excitation check of the headings, and into
      ! the QTF store if the frequency is also a difference frequency.
      TYPE(so_options), INTENT(IN) :: opts
      INTEGER, INTENT(IN) :: kk

      IF (.NOT. body_terms .OR. opts%qtf_terms < QTF_TERMS_BODY .OR. is_limit_frequency()) RETURN
      IF (ABS(W1 - WVFQ(kk)) > 1D-12*W1) THEN
         WRITE (*, '(a,I0)') ' QTF: so_after_radiation called with the wave of another frequency than ', kk
         ERROR STOP 1
      END IF
      IF (multi_arrays) THEN
         CALL radiation_psi(mesh, W1, MXPOT_MULTI_COMB, psi_main)
      ELSE
         CALL radiation_psi(mesh, W1, MXPOT, psi_main)
      END IF
      IF (qtf_psi_index(W1) > 0) THEN
         CALL qtf_keep_psi(W1, psi_main)
         IF (b2_control) CALL keep_cs_psi(opts, W1, psi_main)
         IF (fs_terms) CALL keep_fs_psi(W1, psi_main)
      END IF
      CALL write_psi_rows(W1, 1, psi_main)
   END SUBROUTINE so_after_radiation

   SUBROUTINE so_after_radiation_only(opts, omega)
      ! Call after the radiation-only solve at the difference frequency omega (one of qtf_radiation_omegas()).
      TYPE(so_options), INTENT(IN) :: opts
      REAL(8), INTENT(IN) :: omega
      COMPLEX(8), ALLOCATABLE :: psi(:, :)

      IF (.NOT. body_terms .OR. opts%qtf_terms < QTF_TERMS_BODY) RETURN
      ALLOCATE (psi(mesh%n, 6*mesh%nbody))
      IF (multi_arrays) THEN
         CALL radiation_psi(mesh, omega, MXPOT_MULTI_COMB, psi)
      ELSE
         CALL radiation_psi(mesh, omega, MXPOT, psi)
      END IF
      CALL qtf_keep_psi(omega, psi)
      IF (b2_control) CALL keep_cs_psi(opts, omega, psi)
      IF (fs_terms) CALL keep_fs_psi(omega, psi)
      CALL write_psi_rows(omega, 2, psi)
   END SUBROUTINE so_after_radiation_only

   SUBROUTINE keep_cs_psi(opts, omega, psi)
      ! grad psi_k at the points of C and F, into the QTF store (qtf_b2 = control). psi_k is the perturbation potential of a
      ! state with omega, no incident wave and dpsi_k/dn = n_k (mesh%nj) on the panels of the body of mode k, 0 on the others,
      ! as radiation_psi defines it; fields_green takes its gradient from the boundary integral formula. The solver globals
      ! (V, WK, H and the Green function) must be at omega, as after the radiation solve there.
      TYPE(so_options), INTENT(IN) :: opts
      REAL(8), INTENT(IN) :: omega
      COMPLEX(8), INTENT(IN) :: psi(:, :)
      TYPE(first_order_state), ALLOCATABLE :: ps(:)
      TYPE(cs_fields), ALLOCATABLE :: f(:)
      COMPLEX(8), ALLOCATABLE :: grad(:, :, :)
      INTEGER(8) :: t0, t1, rate
      INTEGER :: k, nmode

      IF (ABS(W1 - omega) > 1D-12*omega) THEN
         WRITE (*, '(a,ES14.6,a)') ' QTF: grad psi on the control surface at omega =', omega, &
            ' rad/s, but the solver is at another frequency.'
         ERROR STOP 1
      END IF
      nmode = 6*mesh%nbody
      ALLOCATE (f(nmode), grad(3, b2_pts%n, nmode))
      CALL psi_states(omega, psi, ps)
      CALL SYSTEM_CLOCK(t0, rate)
      CALL fields_green(mesh, b2_pts, ps, opts%cs_fd_step, multi_arrays, f)
      CALL SYSTEM_CLOCK(t1)
      DO k = 1, nmode
         grad(:, :, k) = f(k)%grad
      END DO
      CALL qtf_keep_cs_psi(omega, grad)
      WRITE (*, '(a,F10.5,a,F10.2,a,I0,a,I0,a)') ' QTF B2 control surface grad psi: omega', omega, ',', &
         REAL(t1 - t0, 8)/REAL(rate, 8), ' s for ', b2_pts%n, ' points and ', nmode, ' modes'
   END SUBROUTINE keep_cs_psi

   SUBROUTINE psi_states(omega, psi, ps)
      ! One state per mode k for the field routines: psi_k as the perturbation potential, with omega, no incident wave and
      ! dpsi_k/dn = n_k (mesh%nj) on the panels of the body of mode k, 0 on the others, as radiation_psi defines it. The
      ! solver globals (V, WK, H) must be at omega.
      REAL(8), INTENT(IN) :: omega
      COMPLEX(8), INTENT(IN) :: psi(:, :)
      TYPE(first_order_state), ALLOCATABLE, INTENT(OUT) :: ps(:)
      INTEGER :: k, kb, i, nmode

      nmode = 6*mesh%nbody
      ALLOCATE (ps(nmode))
      DO k = 1, nmode
         kb = (k - 1)/6 + 1
         ps(k)%omega = omega
         ps(k)%k = WK
         ps(k)%nu = V
         ps(k)%h = H
         ps(k)%amp = 0D0
         ps(k)%beta = 0D0
         ps(k)%g = G
         ps(k)%rho = RHO
         ALLOCATE (ps(k)%phi_p(mesh%n), ps(k)%dphi_p_dn(mesh%n))
         ps(k)%phi_p = psi(:, k)
         DO i = 1, mesh%n
            ps(k)%dphi_p_dn(i) = (0D0, 0D0)
            IF (mesh%body(i) == kb) ps(k)%dphi_p_dn(i) = mesh%nj(k - 6*(kb - 1), i)
         END DO
      END DO
   END SUBROUTINE psi_states

   SUBROUTINE keep_fs_psi(omega, psi)
      ! qtf_terms = full: psi_k at the inner free-surface points (value only, stock route) and the partial waves of each
      ! psi_k outside R (samples by the centroid route, with the mode check), into the QTF store. The solver globals
      ! (V, WK, WVN, H and the Green function) must be at omega, as after the radiation solve there.
      REAL(8), INTENT(IN) :: omega
      COMPLEX(8), INTENT(IN) :: psi(:, :)
      TYPE(first_order_state), ALLOCATABLE :: ps(:)
      COMPLEX(8), ALLOCATABLE :: psi_pts(:, :)
      TYPE(pw_set), ALLOCATABLE :: pw(:)
      REAL(8), ALLOCATABLE :: miss(:)
      REAL(8) :: kev(fs_nev)
      INTEGER(8) :: t0, t1, rate
      INTEGER :: nmode

      IF (ABS(W1 - omega) > 1D-12*omega) THEN
         WRITE (*, '(a,ES14.6,a)') ' QTF: psi on the inner free surface at omega =', omega, &
            ' rad/s, but the solver is at another frequency.'
         ERROR STOP 1
      END IF
      nmode = 6*mesh%nbody
      CALL psi_states(omega, psi, ps)
      ALLOCATE (psi_pts(fs_pts%n, nmode), pw(nmode), miss(nmode))
      CALL SYSTEM_CLOCK(t0, rate)
      CALL fs_sample_phi(mesh, fs_pts, ps, multi_arrays, psi_pts)
      CALL pw_evanescent_from_wvn(fs_nev, kev)
      CALL fs_circle_modes(mesh, fs_circle, ps, multi_arrays, fs_centre, fs_radius, H, WK, kev, fs_mmax, fs_nev, pw, miss)
      CALL SYSTEM_CLOCK(t1)
      CALL log_mode_check(miss, omega, 'psi of mode')
      CALL qtf_keep_fs_psi(omega, psi_pts, pw)
      WRITE (*, '(a,F10.5,a,F10.2,a,I0,a,I0,a,I0,a)') ' QTF free surface psi: omega', omega, ',', &
         REAL(t1 - t0, 8)/REAL(rate, 8), ' s for ', fs_pts%n, ' inner points, ', fs_circle%n, ' samples on the circle and ', &
         nmode, ' modes'
   END SUBROUTINE keep_fs_psi

   SUBROUTINE log_mode_check(miss, omega, what)
      ! The mode check of fs_circle_modes: one warning for the worst of the states, naming it (what and its index)
      REAL(8), INTENT(IN) :: miss(:), omega
      CHARACTER(LEN=*), INTENT(IN) :: what
      CHARACTER(LEN=:), ALLOCATABLE :: msg

      IF (SIZE(miss) == 0) RETURN
      msg = fs_mode_warning(MAXVAL(miss))
      IF (LEN(msg) > 0) WRITE (*, '(a,a,a,F10.5,a,a,a,I0,a)') ' ', msg, ' (omega', omega, ' rad/s, ', what, ' ', &
         MAXLOC(miss, 1), ')'
   END SUBROUTINE log_mode_check

   LOGICAL FUNCTION keep_all_headings(opts)
      ! The states of every heading of a frequency are kept: the pairs of mean_drift = 2, the fields of the control surface,
      ! those at the field points of qtf_field and those on the free surface of qtf_terms = full
      TYPE(so_options), INTENT(IN) :: opts

      keep_all_headings = opts%mean_drift == 2 .OR. opts%do_control .OR. uses_b2_control(opts) .OR. &
                          opts%qtf_field == 1 .OR. uses_fs(opts)
   END FUNCTION keep_all_headings

   LOGICAL FUNCTION uses_fs(opts)
      ! qtf_terms = full applies: the free-surface integral of the difference QTF is computed
      TYPE(so_options), INTENT(IN) :: opts

      uses_fs = opts%qtf_terms == QTF_TERMS_FULL .AND. (opts%qtf == 1 .OR. opts%qtf == 3)
   END FUNCTION uses_fs

   LOGICAL FUNCTION uses_b2_control(opts)
      ! qtf_b2 = control applies: the potential part of the difference QTF is computed
      TYPE(so_options), INTENT(IN) :: opts

      uses_b2_control = opts%qtf_terms >= QTF_TERMS_BODY .AND. opts%qtf_b2 == QTF_B2_CONTROL .AND. &
                        (opts%qtf == 1 .OR. opts%qtf == 3)
   END FUNCTION uses_b2_control

   SUBROUTINE b2_control_setup(opts, inputdir, outputdir)
      ! C and F of every body for qtf_b2 = control. With drift_methods = control they are those of control_setup; otherwise
      ! they are built and checked here in the same way, and C is written as CSMesh_<b>.pnl to outputdir. Unlike the drift,
      ! B2 has no partial fallback: a body without F, or whose C has no bottom, stops the run.
      TYPE(so_options), INTENT(IN) :: opts
      CHARACTER(LEN=*), INTENT(IN) :: inputdir, outputdir
      TYPE(panel_mesh), ALLOCATABLE :: cms(:)
      CHARACTER(LEN=:), ALLOCATABLE :: msg
      LOGICAL :: ok
      LOGICAL, ALLOCATABLE :: free_built(:)
      INTEGER :: i

      IF (.NOT. opts%do_control) THEN
         CALL build_control_surface(mesh, geom, opts, H, inputdir, cs_pts, cms, ok, msg, with_free=.TRUE., &
                                    free_built=free_built)
         IF (.NOT. ok) THEN
            WRITE (*, '(a,a)') ' QTF B2 control surface: ', msg
            ERROR STOP 1
         END IF
         IF (LEN(msg) > 0) WRITE (*, '(a,a)') ' QTF B2 control surface: ', msg
         CALL check_control_surface(mesh, cs_pts, cms, opts%cs_shape /= 3, ok, msg)
         IF (.NOT. ok) THEN
            WRITE (*, '(a,a)') ' QTF B2 control surface: ', msg
            ERROR STOP 1
         END IF
         CALL write_control_surface(outputdir, cms)
         IF (opts%cs_velocity == 2) ALLOCATE (cs_sigma(mesh%n, NBETA))
      END IF
      CALL b2_control_check(mesh%nbody, cs_pts, ok, msg)
      IF (.NOT. ok) THEN
         WRITE (*, '(a,a)') ' ', msg
         ERROR STOP 1
      END IF
      ! B2 uses the points of C and F only: the fields at the segments of Gamma_c are neither stored nor, without the
      ! control-surface drift, evaluated (cs_pts then holds the points of C and F only)
      b2_idx = PACK([(i, i = 1, cs_pts%n)], cs_pts%kind /= CS_LINE)
      CALL subset_points(cs_pts, b2_idx, b2_pts)
      IF (.NOT. opts%do_control) THEN
         cs_pts = b2_pts
         b2_idx = [(i, i = 1, cs_pts%n)]
      END IF
      WRITE (*, '(a,I0,a,I0,a,I0,a)') ' QTF B2 control surface: ', COUNT(cs_pts%kind == CS_WALL), ' panels on C and ', &
         COUNT(cs_pts%kind == CS_FREE), ' panels on F, for ', mesh%nbody, ' bodies'
   END SUBROUTINE b2_control_setup

   SUBROUTINE subset_points(pts, idx, sub)
      ! The points idx of pts
      TYPE(cs_points), INTENT(IN) :: pts
      INTEGER, INTENT(IN) :: idx(:)
      TYPE(cs_points), INTENT(OUT) :: sub

      sub%n = SIZE(idx)
      sub%x = pts%x(:, idx)
      sub%nrm = pts%nrm(:, idx)
      sub%w = pts%w(idx)
      sub%psize = pts%psize(idx)
      sub%kind = pts%kind(idx)
      sub%body = pts%body(idx)
   END SUBROUTINE subset_points

   SUBROUTINE control_setup(opts, inputdir, outputdir)
      ! Builds and checks C and F of every body, writes C as CSMesh_<b>.pnl to outputdir, and opens MeanDriftControlSurface.7.
      ! If the generated F of a body cannot be built (for example a waterline that is not star-shaped about the body centre),
      ! the run goes on without F for that body: its F_x, F_y and M_z are written, and the log says why its F_z, M_x and M_y
      ! are not. The other bodies keep their F. A FSMesh_<b>.pnl that cannot be used stops the run, as a failing C does.
      TYPE(so_options), INTENT(IN) :: opts
      CHARACTER(LEN=*), INTENT(IN) :: inputdir, outputdir
      TYPE(panel_mesh), ALLOCATABLE :: cms(:)
      CHARACTER(LEN=:), ALLOCATABLE :: msg
      LOGICAL :: ok
      LOGICAL, ALLOCATABLE :: free_built(:)
      INTEGER :: b

      CALL build_control_surface(mesh, geom, opts, H, inputdir, cs_pts, cms, ok, msg, with_free=.TRUE., free_built=free_built)
      IF (.NOT. ok) THEN
         WRITE (*, '(a,a)') ' Mean drift control surface: ', msg
         ERROR STOP 1
      END IF
      have_free = ANY(free_built)
      ! msg holds the warnings: the waterline of a C read from a file, and each body without F
      IF (LEN(msg) > 0) WRITE (*, '(a,a)') ' Mean drift control surface: ', msg
      CALL check_control_surface(mesh, cs_pts, cms, opts%cs_shape /= 3, ok, msg)
      IF (.NOT. ok) THEN
         WRITE (*, '(a,a)') ' Mean drift control surface: ', msg
         ERROR STOP 1
      END IF
      CALL write_control_surface(outputdir, cms)
      WRITE (*, '(a,I0,a,I0,a,I0,a,I0,a)') ' Mean drift control surface: ', COUNT(cs_pts%kind == CS_WALL), ' panels on C, ', &
         COUNT(cs_pts%kind == CS_LINE), ' segments on Gamma_c and ', COUNT(cs_pts%kind == CS_FREE), ' panels on F, for ', &
         mesh%nbody, ' bodies'
      ! A C that reaches the sea bed has no bottom panel, whose normal points up towards the body. Its domain is then also
      ! bounded by the sea bed between the body and C, which the vertical loads would need (not implemented).
      ALLOCATE (cs_vertical(mesh%nbody))
      DO b = 1, mesh%nbody
         cs_vertical(b) = free_built(b) .AND. ANY(cs_pts%kind == CS_WALL .AND. cs_pts%body == b .AND. cs_pts%nrm(3, :) > 0.5D0)
         IF (free_built(b) .AND. .NOT. cs_vertical(b)) THEN
            WRITE (*, '(a,I0,a,a,a)') ' Mean drift: control surface of body ', b, ' has no bottom (it reaches the sea bed, or', &
               ' the file has none): F_z, M_x, M_y not computed. For a generated C, set cs_bottom_margin below h minus', &
               ' the draft to keep a bottom.'
         END IF
      END DO
      IF (opts%cs_velocity == 2) ALLOCATE (cs_sigma(mesh%n, NBETA))
      IF (opts%do_far) ALLOCATE (far_diag(3, NBETA))
      OPEN (NEWUNIT=unit_control, FILE=outputdir//'/Wamit_format/MeanDriftControlSurface.7', STATUS='REPLACE', ACTION='WRITE')
   END SUBROUTINE control_setup

   LOGICAL FUNCTION is_limit_frequency()
      ! The zero and infinite frequencies added by SYBO = 1
      is_limit_frequency = ABS(TP) < 1D-6 .OR. ABS(TP + 1D0) < 1D-6
   END FUNCTION is_limit_frequency

   SUBROUTINE so_after_motion_single(opts, kk, ib)
      ! Call inside the heading loop, after SolveMotion and before the next diffraction solve overwrites MXPOT.
      TYPE(so_options), INTENT(IN) :: opts
      INTEGER, INTENT(IN) :: kk, ib

      IF (is_limit_frequency()) RETURN
      CALL after_motion(opts, ib, MXPOT, DSPL(kk, ib, :), EXFC(kk, ib, :), MATX, AMAS(kk, :, :), BDMP(kk, :, :) + BLNR, &
                        CRS + KSTF)
   END SUBROUTINE so_after_motion_single

   SUBROUTINE so_after_motion_multi(opts, kk, ib)
      ! Call inside the heading loop, after SolveMotionMulti and before the next diffraction solve overwrites MXPOT.
      TYPE(so_options), INTENT(IN) :: opts
      INTEGER, INTENT(IN) :: kk, ib

      IF (is_limit_frequency()) RETURN
      CALL after_motion(opts, ib, MXPOT_MULTI_COMB, DSPL_MULTI_COMB(kk, ib, :), EXFC_MULTI_COMB(kk, ib, :), &
                        block_diagonal(MATX_MULTI), AMAS_MULTI_COMB(kk, :, :), &
                        BDMP_MULTI_COMB(kk, :, :) + block_diagonal(BLNR_MULTI), &
                        block_diagonal(CRS_MULTI) + block_diagonal(KSTF_MULTI))
   END SUBROUTINE so_after_motion_multi

   FUNCTION block_diagonal(per_body) RESULT(full)
      REAL(8), INTENT(IN) :: per_body(:, :, :)   ! (nbody, 6, 6)
      REAL(8) :: full(6*SIZE(per_body, 1), 6*SIZE(per_body, 1))
      INTEGER :: b, r0

      full = 0D0
      DO b = 1, SIZE(per_body, 1)
         r0 = 6*(b - 1)
         full(r0 + 1:r0 + 6, r0 + 1:r0 + 6) = per_body(b, :, :)
      END DO
   END FUNCTION block_diagonal

   SUBROUTINE after_motion(opts, ib, mxpot_in, dspl_in, exfc_in, mass, added, damping, stiffness)
      ! mass, added, damping, stiffness and exfc_in are used for the diagnostics only.
      TYPE(so_options), INTENT(IN) :: opts
      INTEGER, INTENT(IN) :: ib
      COMPLEX(8), INTENT(IN) :: mxpot_in(:, :, :), dspl_in(:), exfc_in(:)
      REAL(8), INTENT(IN) :: mass(:, :), added(:, :), damping(:, :), stiffness(:, :)
      COMPLEX(8) :: f(3)
      TYPE(first_order_state) :: scratch
      INTEGER :: slot, n_fallback
      LOGICAL :: control_source

      slot = MERGE(ib, 1, keep_all_headings(opts))
      control_source = (opts%do_control .OR. b2_control) .AND. opts%cs_velocity == 2
      CALL build_state(mesh, W1, WK, V, H, 1D0, BETA, G, RHO, mxpot_in, dspl_in, opts%motion_fixed, states(slot))
      ! As for the waterline rows below: the operator holds for all the headings of a frequency
      IF (ib == 1 .AND. (hull_uses_source .OR. control_source)) CALL timed_source_build()
      IF (control_source .AND. .NOT. hull_uses_source) THEN
         ! Only sigma is wanted: the solve also overwrites grad_p, which the near-field method may take from differencing
         scratch = states(slot)
         CALL timed_source_velocity(scratch, cs_sigma(:, slot))
      END IF
      ! The hull velocity and waterline potential, for the near-field method and the Phase B control-surface loads
      IF (ALLOCATED(phi_wl)) THEN
         IF (opts%velocity_method == 2) THEN
            IF (control_source) THEN
               CALL timed_source_velocity(states(slot), cs_sigma(:, slot))
            ELSE
               CALL timed_source_velocity(states(slot))
            END IF
         ELSE
            CALL velocity_differencing(mesh, geom, states(slot), n_fallback)
            n_fallback_total = n_fallback_total + n_fallback
         END IF
         n_velocity_fields = n_velocity_fields + 1
         IF (opts%waterline_method == 2) THEN
            ! The headings of a frequency come in order, and the rows hold for all of them
            IF (ib == 1) CALL waterline_rows(mesh, geom, opts%waterline_offset, multi_arrays, rows_g, rows_dg)
            CALL waterline_field(mesh, geom, states(slot), rows_g, rows_dg, phi_wl(:, slot))
         ELSE
            CALL waterline_extrapolate(mesh, geom, states(slot), phi_wl(:, slot))
         END IF
      END IF

      ! The QTFs need every frequency: keep the state with its hull velocity and waterline potential
      IF (opts%qtf > 0) CALL qtf_keep(ib, states(slot), phi_wl(:, slot), OUFR)
      ! psi of this frequency was taken by so_after_radiation, before the diffraction solves
      IF (body_terms) CALL write_exc_rows(states(slot), WVHD(ib))

      IF (opts%write_terms == 2) THEN
         IF (ib == 1) CALL write_wave_and_modes(states(slot), mxpot_in)
         CALL write_heading_diagnostics(opts, states(slot), mxpot_in, exfc_in, mass, added, damping, stiffness)
      END IF

      IF (opts%do_far .AND. opts%mean_drift == 1) THEN
         CALL far_field_pair(mesh, states(slot), states(slot), opts, f)
         CALL write_far_rows(WVHD(ib), WVHD(ib), f)
         IF (ALLOCATED(far_diag)) far_diag(:, ib) = f
      END IF
      IF (opts%do_near .AND. opts%mean_drift == 1) THEN
         IF (opts%do_far) THEN
            CALL near_field_output(slot, slot, WVHD(ib), WVHD(ib), f)
         ELSE
            CALL near_field_output(slot, slot, WVHD(ib), WVHD(ib))
         END IF
      END IF
   END SUBROUTINE after_motion

   SUBROUTINE timed_source_build()
      INTEGER(8) :: t0, t1, rate

      CALL SYSTEM_CLOCK(t0, rate)
      CALL build_source_operator(mesh, multi_arrays, ALLOCATED(phi_source), source_op)
      CALL SYSTEM_CLOCK(t1)
      source_build_seconds = source_build_seconds + REAL(t1 - t0, 8)/REAL(rate, 8)
      n_source_builds = n_source_builds + 1
   END SUBROUTINE timed_source_build

   SUBROUTINE timed_source_velocity(s, sigma)
      ! sigma, if given, receives the source strength on every panel, for fields_source
      TYPE(first_order_state), INTENT(INOUT) :: s
      COMPLEX(8), INTENT(OUT), OPTIONAL :: sigma(:)
      INTEGER(8) :: t0, t1, rate

      CALL SYSTEM_CLOCK(t0, rate)
      IF (ALLOCATED(phi_source) .AND. PRESENT(sigma)) THEN
         CALL velocity_source(mesh, source_op, s, phi_source, sigma_full=sigma)
      ELSE IF (ALLOCATED(phi_source)) THEN
         CALL velocity_source(mesh, source_op, s, phi_source)
      ELSE IF (PRESENT(sigma)) THEN
         CALL velocity_source(mesh, source_op, s, sigma_full=sigma)
      ELSE
         CALL velocity_source(mesh, source_op, s)
      END IF
      CALL SYSTEM_CLOCK(t1)
      source_solve_seconds = source_solve_seconds + REAL(t1 - t0, 8)/REAL(rate, 8)
      n_source_solves = n_source_solves + 1
   END SUBROUTINE timed_source_velocity

   SUBROUTINE so_end_frequency(opts, kk)
      ! Call after the heading loop of frequency kk. Evaluates the heading pairs when mean_drift = 2, the field points of
      ! qtf_field, the free-surface fields of qtf_terms = full and the control-surface method.
      TYPE(so_options), INTENT(IN) :: opts
      INTEGER, INTENT(IN) :: kk
      TYPE(cs_fields), ALLOCATABLE :: fields(:)
      INTEGER :: ib

      IF (is_limit_frequency()) RETURN
      IF (opts%mean_drift == 2) CALL heading_pairs(opts)
      IF (opts%qtf_field == 1) CALL field_frequency()
      IF (fs_terms) CALL fs_frequency()
      IF (.NOT. (opts%do_control .OR. b2_control)) RETURN
      ! The fields on C and F, once for the control-surface drift and for B2 (qtf_b2 = control)
      ALLOCATE (fields(SIZE(states)))
      CALL control_fields(opts, fields)
      IF (opts%do_control) CALL control_frequency(opts, fields)
      IF (b2_control) THEN
         DO ib = 1, NBETA
            CALL qtf_keep_cs_velocity(ib, fields(ib)%grad(:, b2_idx))
         END DO
      END IF
   END SUBROUTINE so_end_frequency

   SUBROUTINE field_frequency()
      ! The first-order potential and gradient of every heading at the field points, kept for the partial second-order field.
      ! The states are at the current global frequency, which fields_green reads. The route is always the Green function with
      ! the central-difference step FIELD_FD_STEP of the nearest hull panel, whatever velocity_method is. Above nu*h = 50 the
      ! field mode of FINGREEN3D keeps the Limes acceleration of PIDCOCK, as on the control surface, so the same warning.
      TYPE(cs_fields), ALLOCATABLE :: fields(:)
      REAL(8), PARAMETER :: LIMES_NUH = 50D0
      INTEGER(8) :: t0, t1, rate

      IF (H >= 0D0 .AND. states(1)%nu*H > LIMES_NUH) THEN
         WRITE (*, '(a,F10.5,a,F8.2,a)') ' QTF field points warning: omega', W1, ', nu*h =', states(1)%nu*H, &
            ' > 50: the finite-depth field values may scatter (the Green function series still uses Limes there).'
      END IF
      ALLOCATE (fields(SIZE(states)))
      CALL SYSTEM_CLOCK(t0, rate)
      CALL fields_green(mesh, fp_pts, states, FIELD_FD_STEP, multi_arrays, fields)
      CALL SYSTEM_CLOCK(t1)
      WRITE (*, '(a,F10.5,a,F10.2,a,I0,a,I0,a)') ' QTF field points: omega', W1, ',', REAL(t1 - t0, 8)/REAL(rate, 8), &
         ' s for ', fp_pts%n, ' points and ', SIZE(states), ' headings'
      CALL qfield_keep(states, fields, OUFR)
   END SUBROUTINE field_frequency

   SUBROUTINE fs_frequency()
      ! qtf_terms = full: the total first-order fields of every heading at the inner free-surface points (fs_inner_fields:
      ! stock route, horizontal Laplacian from the offsets fs_delta) and the partial waves of the scattered potential outside
      ! R (fs_circle_modes: samples by the centroid route, with the mode check), into the QTF store. The states are at the
      ! current global frequency, whose WVN gives the evanescent wave numbers. Above nu*h = 50 the field mode of FINGREEN3D
      ! keeps the Limes acceleration of PIDCOCK, as on the control surface, so the same warning.
      TYPE(fs_point_field), ALLOCATABLE :: f(:, :)
      TYPE(pw_set), ALLOCATABLE :: pw(:)
      REAL(8), ALLOCATABLE :: miss(:)
      REAL(8), PARAMETER :: LIMES_NUH = 50D0
      REAL(8) :: kev(fs_nev)
      INTEGER(8) :: t0, t1, rate
      INTEGER :: ib

      IF (states(1)%nu*H > LIMES_NUH) THEN
         WRITE (*, '(a,F10.5,a,F8.2,a)') ' QTF free surface warning: omega', W1, ', nu*h =', states(1)%nu*H, &
            ' > 50: the finite-depth field values may scatter (the Green function series still uses Limes there).'
      END IF
      ALLOCATE (f(fs_pts%n, SIZE(states)), pw(SIZE(states)), miss(SIZE(states)))
      CALL SYSTEM_CLOCK(t0, rate)
      CALL fs_inner_fields(mesh, fs_pts, states, fs_step, multi_arrays, f, fs_delta)
      CALL pw_evanescent_from_wvn(fs_nev, kev)
      CALL fs_circle_modes(mesh, fs_circle, states, multi_arrays, fs_centre, fs_radius, H, WK, kev, fs_mmax, fs_nev, pw, miss)
      CALL SYSTEM_CLOCK(t1)
      CALL log_mode_check(miss, W1, 'heading')
      DO ib = 1, SIZE(states)
         CALL qtf_keep_fs_state(ib, f(:, ib), pw(ib))
      END DO
      WRITE (*, '(a,F10.5,a,F10.2,a,I0,a,I0,a,I0,a)') ' QTF free surface fields: omega', W1, ',', &
         REAL(t1 - t0, 8)/REAL(rate, 8), ' s for ', fs_pts%n, ' inner points, ', fs_circle%n, ' samples on the circle and ', &
         SIZE(states), ' headings'
   END SUBROUTINE fs_frequency

   SUBROUTINE heading_pairs(opts)
      ! The far-field and near-field loads of every heading pair of the frequency
      TYPE(so_options), INTENT(IN) :: opts
      COMPLEX(8), ALLOCATABLE :: f(:, :, :)
      INTEGER :: i, j

      IF (opts%do_far) THEN
         ALLOCATE (f(3, NBETA, NBETA))
         CALL far_field_pairs(mesh, states, opts, f)
         DO j = 1, NBETA
            DO i = 1, NBETA
               CALL write_far_rows(WVHD(i), WVHD(j), f(:, i, j))
            END DO
            IF (ALLOCATED(far_diag)) far_diag(:, j) = f(:, j, j)
         END DO
      END IF
      IF (opts%do_near) THEN
         DO j = 1, NBETA
            DO i = 1, NBETA
               IF (opts%do_far .AND. i == j) THEN
                  CALL near_field_output(i, j, WVHD(i), WVHD(j), f(:, i, j))
               ELSE
                  CALL near_field_output(i, j, WVHD(i), WVHD(j))
               END IF
            END DO
         END DO
      END IF
   END SUBROUTINE heading_pairs

   SUBROUTINE control_fields(opts, fields)
      ! The fields on C (and F) of every heading of the frequency in one pass over the points. The states are at the current
      ! global frequency, which fields_green and fields_source read.
      ! Above nu*h = 50 (VHLIN in FinGreen3D.f90) the field mode of FINGREEN3D keeps the Limes acceleration of PIDCOCK
      ! (FIELD_SERIES), whose isolated outliers a central difference takes over (FieldVelocity.f90), so the green route warns.
      TYPE(so_options), INTENT(IN) :: opts
      TYPE(cs_fields), INTENT(OUT) :: fields(:)
      REAL(8), PARAMETER :: LIMES_NUH = 50D0
      INTEGER(8) :: t0, t1, rate
      CHARACTER(LEN=:), ALLOCATABLE :: who

      who = 'Mean drift control surface'
      IF (.NOT. opts%do_control) who = 'QTF B2 control surface'
      IF (opts%cs_velocity /= 2 .AND. H >= 0D0 .AND. states(1)%nu*H > LIMES_NUH) THEN
         WRITE (*, '(a,F10.5,a,F8.2,a)') ' '//who//' warning: omega', W1, ', nu*h =', states(1)%nu*H, &
            ' > 50: the finite-depth field values on C may scatter (the Green function series still uses Limes there).'
         IF (opts%do_control) THEN
            WRITE (*, '(a)') ' Check against drift_methods = far, or use deep water (H < 0) where kh is large.'
         ELSE
            WRITE (*, '(a)') ' B2 takes the velocity there: check against qtf_b2 = body, or use cs_velocity = source or deep'// &
               ' water (H < 0) where kh is large.'
         END IF
      END IF
      CALL SYSTEM_CLOCK(t0, rate)
      IF (opts%cs_velocity == 2) THEN
         CALL fields_source(mesh, cs_pts, states, cs_sigma, multi_arrays, fields)
      ELSE
         CALL fields_green(mesh, cs_pts, states, opts%cs_fd_step, multi_arrays, fields)
      END IF
      CALL SYSTEM_CLOCK(t1)
      WRITE (*, '(a,F10.5,a,F10.2,a,I0,a,I0,a)') ' '//who//' fields: omega', W1, ',', &
         REAL(t1 - t0, 8)/REAL(rate, 8), ' s for ', cs_pts%n, ' points and ', SIZE(states), ' headings'
   END SUBROUTINE control_fields

   SUBROUTINE control_frequency(opts, fields)
      ! The control-surface load of each pair of the frequency, from the fields of control_fields: (i, i) for mean_drift = 1,
      ! all pairs for mean_drift = 2.
      TYPE(so_options), INTENT(IN) :: opts
      TYPE(cs_fields), INTENT(IN) :: fields(:)
      INTEGER :: i, j

      DO j = 1, NBETA
         DO i = 1, NBETA
            IF (opts%mean_drift /= 2 .AND. i /= j) CYCLE
            IF (ALLOCATED(far_diag) .AND. i == j) THEN
               CALL control_output(fields, i, j, WVHD(i), WVHD(j), far_diag(:, i))
            ELSE
               CALL control_output(fields, i, j, WVHD(i), WVHD(j))
            END IF
         END DO
      END DO
   END SUBROUTINE control_frequency

   SUBROUTINE control_output(fields, i1, i2, beta1_deg, beta2_deg, far)
      ! Writes the control-surface load of the state pair (i1, i2), per body, about XR of each body: I = 1, 2, 6 (Phase A),
      ! and I = 3, 4, 5 (Phase B) for the bodies of cs_vertical. The terms file has terms 11 to 13 for I = 1, 2, 6, which
      ! the Phase B terms do not change, and 11 to 16 for I = 3, 4, 5. With far, logs the difference from it, with M_z moved
      ! to the origin as in near_field_output.
      TYPE(cs_fields), INTENT(IN) :: fields(:)
      INTEGER, INTENT(IN) :: i1, i2
      REAL(8), INTENT(IN) :: beta1_deg, beta2_deg
      COMPLEX(8), INTENT(IN), OPTIONAL :: far(3)
      COMPLEX(8) :: terms(6, 6, mesh%nbody), total(6), control(3)
      INTEGER :: b, i, term, last

      IF (have_free) THEN
         CALL control_surface_pair(cs_pts, mesh%xr, states(i1), states(i2), fields(i1), fields(i2), terms, mesh, geom, &
                                   phi_wl(:, i1), phi_wl(:, i2))
      ELSE
         CALL control_surface_pair(cs_pts, mesh%xr, states(i1), states(i2), fields(i1), fields(i2), terms)
      END IF
      control = (0D0, 0D0)
      DO b = 1, mesh%nbody
         total = SUM(terms(:, :, b), DIM=2)
         DO i = 1, 6
            IF (i >= 3 .AND. i <= 5 .AND. .NOT. cs_vertical(b)) CYCLE
            CALL write_drift_row(unit_control, OUFR, beta1_deg, beta2_deg, 6*(b - 1) + i, total(i), REFL, i > 3)
            IF (unit_terms == 0) CYCLE
            last = MERGE(CS_T_HYDROSTATIC, CS_T_PRESSURE, i >= 3 .AND. i <= 5)
            DO term = CS_T_LINE, last
               CALL write_terms_row(unit_terms, OUFR, beta1_deg, beta2_deg, b, i, 10 + term, terms(i, term, b), REFL, i > 3)
            END DO
         END DO
         control(1:2) = control(1:2) + total(1:2)
         control(3) = control(3) + total(6) + mesh%xr(1, b)*total(2) - mesh%xr(2, b)*total(1)
      END DO
      IF (PRESENT(far)) THEN
         CALL log_against_far('control', 'F_x', beta1_deg, far(1), control(1))
         CALL log_against_far('control', 'F_y', beta1_deg, far(2), control(2))
         CALL log_against_far('control', 'M_z', beta1_deg, far(3), control(3))
      END IF
   END SUBROUTINE control_output

   SUBROUTINE so_finish()
      ! Writes the QTFs, closes the result files and logs the velocity method: what differencing could not do, or the time of
      ! the source solve. Call after the frequency loop.
      CALL qtf_finish(mesh, geom, REFL)
      CALL qfield_finish(REFL)
      IF (n_source_builds > 0) THEN
         WRITE (*, '(a,F10.2,a,I0,a,F10.2,a,I0,a)') ' Mean drift source solve: ', source_build_seconds, &
            ' s to assemble and factorize at ', n_source_builds, ' frequencies, ', source_solve_seconds, &
            ' s to solve ', n_source_solves, ' first-order solutions'
      END IF
      IF (ALLOCATED(phi_wl) .AND. .NOT. hull_uses_source) THEN
         WRITE (*, '(a,I0,a,I0,a)') ' Mean drift velocity fallback: ', n_fallback_total, &
            ' panel values over ', n_velocity_fields, ' first-order solutions have only the normal part of the velocity'
      END IF
      IF (unit_far /= 0) CLOSE (unit_far)
      IF (unit_near /= 0) CLOSE (unit_near)
      IF (unit_control /= 0) CLOSE (unit_control)
      IF (unit_terms /= 0) CLOSE (unit_terms)
      IF (unit_diag /= 0) CLOSE (unit_diag)
      IF (unit_psi /= 0) CLOSE (unit_psi)
      IF (ALLOCATED(psi_main)) DEALLOCATE (psi_main)
   END SUBROUTINE so_finish

   SUBROUTINE near_field_output(i1, i2, beta1_deg, beta2_deg, far)
      ! Writes the near-field load of the state pair (i1, i2), per body, about XR of each body. With far (F_x, F_y, M_z
      ! about the origin), also logs the difference from it, which needs M_z about the origin: M_z + (xr x F)_z per body.
      INTEGER, INTENT(IN) :: i1, i2
      REAL(8), INTENT(IN) :: beta1_deg, beta2_deg
      COMPLEX(8), INTENT(IN), OPTIONAL :: far(3)
      COMPLEX(8) :: terms(6, 6, mesh%nbody), total(6), near(3)
      INTEGER :: b, i, term

      CALL near_field_pair(mesh, geom, states(i1), states(i2), phi_wl(:, i1), phi_wl(:, i2), terms)
      near = (0D0, 0D0)
      DO b = 1, mesh%nbody
         total = SUM(terms(:, :, b), DIM=2)
         DO i = 1, 6
            CALL write_drift_row(unit_near, OUFR, beta1_deg, beta2_deg, 6*(b - 1) + i, total(i), REFL, i > 3)
            IF (unit_terms == 0) CYCLE
            DO term = T_WATERLINE, T_HYDROSTATIC
               CALL write_terms_row(unit_terms, OUFR, beta1_deg, beta2_deg, b, i, term, terms(i, term, b), REFL, i > 3)
            END DO
         END DO
         near(1:2) = near(1:2) + total(1:2)
         near(3) = near(3) + total(6) + mesh%xr(1, b)*total(2) - mesh%xr(2, b)*total(1)
      END DO
      IF (PRESENT(far)) THEN
         CALL log_against_far('near', 'F_x', beta1_deg, far(1), near(1))
         CALL log_against_far('near', 'F_y', beta1_deg, far(2), near(2))
         CALL log_against_far('near', 'M_z', beta1_deg, far(3), near(3))
      END IF
   END SUBROUTINE near_field_output

   SUBROUTINE log_against_far(method, label, beta_deg, far, other)
      ! The relative difference has no meaning for a load that is zero by symmetry, so the values are logged with it.
      ! method is 'near' or 'control'.
      CHARACTER(LEN=*), INTENT(IN) :: method, label
      REAL(8), INTENT(IN) :: beta_deg
      COMPLEX(8), INTENT(IN) :: far, other

      WRITE (*, '(a,a,a,F10.5,a,F8.3,a,ES12.4,a,ES12.4,a,ES10.3)') ' Mean drift '//method//' vs far ', label, ': omega', W1, &
         ' beta', beta_deg, ' far', REAL(far, 8), ' '//method, REAL(other, 8), ' relative difference', &
         ABS(far - other)/MAX(ABS(far), TINY(1D0))
   END SUBROUTINE log_against_far

   SUBROUTINE write_far_rows(beta1_deg, beta2_deg, f)
      ! f = (F_x, F_y, M_z), written as I = 1, 2, 6
      REAL(8), INTENT(IN) :: beta1_deg, beta2_deg
      COMPLEX(8), INTENT(IN) :: f(3)

      CALL write_drift_row(unit_far, OUFR, beta1_deg, beta2_deg, 1, f(1), REFL, .FALSE.)
      CALL write_drift_row(unit_far, OUFR, beta1_deg, beta2_deg, 2, f(2), REFL, .FALSE.)
      CALL write_drift_row(unit_far, OUFR, beta1_deg, beta2_deg, 6, f(3), REFL, .TRUE.)
   END SUBROUTINE write_far_rows

   SUBROUTINE write_wave_and_modes(s, mxpot_in)
      ! The wave properties and, for each mode, the amplitude of the wave radiated by unit motion: phi = -i omega psi_j
      ! on all panels, d(phi)/dn = -i omega n_j on the panels of the moving body.
      TYPE(first_order_state), INTENT(IN) :: s
      COMPLEX(8), INTENT(IN) :: mxpot_in(:, :, :)
      TYPE(first_order_state) :: unit_motion
      COMPLEX(8) :: a, da
      REAL(8) :: theta_deg
      INTEGER :: j, jb, jl, i, it

      WRITE (unit_diag, '(a,4ES24.15)') 'WAVE ', s%omega, s%k, s%nu, s%h
      unit_motion%omega = s%omega
      unit_motion%k = s%k
      unit_motion%nu = s%nu
      unit_motion%h = s%h
      unit_motion%g = s%g
      ALLOCATE (unit_motion%phi_p(mesh%n), unit_motion%dphi_p_dn(mesh%n))
      DO j = 1, 6*mesh%nbody
         jb = (j - 1)/6 + 1
         jl = j - 6*(jb - 1)
         DO i = 1, mesh%n
            unit_motion%phi_p(i) = -CI*s%omega*mxpot_in(mesh%base(i), j, mesh%image(i))
            unit_motion%dphi_p_dn(i) = (0D0, 0D0)
            IF (mesh%body(i) == jb) unit_motion%dphi_p_dn(i) = -CI*s%omega*mesh%nj(jl, i)
         END DO
         DO it = 1, DIAG_N_ANGLES
            theta_deg = 360D0*REAL(it - 1, 8)/REAL(DIAG_N_ANGLES, 8)
            CALL far_amplitude(mesh, unit_motion, theta_deg*PI/180D0, a, da)
            WRITE (unit_diag, '(a,ES24.15,I6,3ES24.15)') 'MODE ', s%omega, j, theta_deg, REAL(a, 8), AIMAG(a)
         END DO
      END DO
   END SUBROUTINE write_wave_and_modes

   SUBROUTINE write_heading_diagnostics(opts, s, mxpot_in, exfc_in, mass, added, damping, stiffness)
      TYPE(so_options), INTENT(IN) :: opts
      TYPE(first_order_state), INTENT(IN) :: s
      COMPLEX(8), INTENT(IN) :: mxpot_in(:, :, :), exfc_in(:)
      REAL(8), INTENT(IN) :: mass(:, :), added(:, :), damping(:, :), stiffness(:, :)
      COMPLEX(8) :: force(SIZE(exfc_in)), xi(SIZE(exfc_in)), a, da, phi_d
      REAL(8) :: beta_deg, theta_deg
      INTEGER :: j, i, it

      beta_deg = s%beta*180D0/PI
      force = CI*exfc_in

      ! With motion_state = fixed the state has no motion, so the equation of motion does not apply to it.
      IF (.NOT. opts%motion_fixed) THEN
         xi = RESHAPE(s%xi, [SIZE(xi)])
         WRITE (unit_diag, '(a,3ES24.15)') 'RESIDUAL ', s%omega, beta_deg, &
            norm_c(MATMUL(-s%omega**2*(mass + added) - CI*s%omega*damping + stiffness, xi) - force)/norm_c(force)
      END IF
      DO j = 1, SIZE(force)
         WRITE (unit_diag, '(a,2ES24.15,I6,2ES24.15)') 'FORCE ', s%omega, beta_deg, j, REAL(force(j), 8), AIMAG(force(j))
      END DO

      DO it = 1, DIAG_N_ANGLES
         theta_deg = 360D0*REAL(it - 1, 8)/REAL(DIAG_N_ANGLES, 8)
         CALL far_amplitude(mesh, s, theta_deg*PI/180D0, a, da)
         WRITE (unit_diag, '(a,5ES24.15)') 'TOTAL ', s%omega, beta_deg, theta_deg, REAL(a, 8), AIMAG(a)
      END DO

      DO i = 1, mesh%n
         phi_d = s%phi_i(i) + CI*mxpot_in(mesh%base(i), SIZE(mxpot_in, 2), mesh%image(i))
         WRITE (unit_diag, '(a,2ES24.15,I8,7ES24.15)') 'PANEL ', s%omega, beta_deg, i, mesh%xc(:, i), &
            REAL(phi_d, 8), AIMAG(phi_d), REAL(s%phi_i(i) + s%phi_p(i), 8), AIMAG(s%phi_i(i) + s%phi_p(i))
      END DO

      IF (ALLOCATED(phi_source)) THEN
         DO i = 1, mesh%n
            WRITE (unit_diag, '(a,2ES24.15,I8,14ES24.15)') 'SOURCE ', s%omega, beta_deg, i, s%phi_p(i), phi_source(i), &
               s%dphi_p_dn(i), SUM(s%grad_p(:, i)*mesh%nrm(:, i)), s%grad_p(:, i)
         END DO
      END IF
   END SUBROUTINE write_heading_diagnostics

   PURE REAL(8) FUNCTION norm_c(z)
      COMPLEX(8), INTENT(IN) :: z(:)
      norm_c = SQRT(SUM(REAL(z, 8)**2 + AIMAG(z)**2))
   END FUNCTION norm_c

END MODULE SecondOrderDriver
