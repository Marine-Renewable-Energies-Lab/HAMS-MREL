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
   ! docs/output_files.md.

   USE Const_mod, ONLY: G, RHO
   USE HAMS_mod, ONLY: IRSP
   USE WaveDyn_mod, ONLY: AMP, BETA, BLNR, H, NBETA, NBODY, OUFR, TP, V, W1, WK, WVHD, AMAS, BDMP, EXFC, DSPL, &
                          AMAS_MULTI_COMB, BDMP_MULTI_COMB, BLNR_MULTI, EXFC_MULTI_COMB, DSPL_MULTI_COMB
   USE Body_mod, ONLY: CRS, KSTF, MATX, REFL, XR, XR_MULTI, CRS_MULTI, KSTF_MULTI, MATX_MULTI
   USE PanelMesh_mod, ONLY: DS, DXYZ_P, GEONRM_P, GEONRM_MULTI_COMB, ISX, ISY, NCN, NCON, NELEM, NSYS, PNSZ, XYZ, XYZ_P, &
                            DS_MULTI_COMB, DXYZ_MULTI_COMB, NCN_MULTI_COMB, NCON_MULTI_COMB, NELEM_MULTI, NTND_MULTI, &
                            PNSZ_MULTI_COMB, XYZ_GLOBAL_MULTI_COMB, XYZ_GLOBAL_MULTI_COMB_P, iXYZ_GLOBAL_MULTI_COMB, &
                            iNCN_MULTI_COMB, iNCON_MULTI_COMB
   USE Inerfs_mod, ONLY: NELEM_PE, iXYZ, iNCN, iNCON
   USE Potentials_mod, ONLY: MXPOT, MXPOT_MULTI_COMB
   USE IO, ONLY: CreateDirectory
   USE SecondOrderTypes, ONLY: so_options, panel_mesh, first_order_state, surface_geom
   USE FirstOrderState, ONLY: expand_images, build_state
   USE SurfaceGeometry, ONLY: build_surface_geom, DEFAULT_SHARP_DEG
   USE SurfaceVelocity, ONLY: velocity_differencing
   USE FarFieldDrift, ONLY: far_amplitude, far_field_pair, far_field_pairs
   USE NearFieldDrift, ONLY: waterline_extrapolate, near_field_pair, T_WATERLINE, T_HYDROSTATIC
   USE WaterlineField, ONLY: waterline_rows, waterline_field
   USE SecondOrderOutput, ONLY: write_drift_row, write_terms_row
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: so_setup_single, so_setup_multi, so_after_motion_single, so_after_motion_multi, so_end_frequency, so_finish

   REAL(8), PARAMETER :: PI = 3.141592653589793D0
   COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
   INTEGER, PARAMETER :: DIAG_N_ANGLES = 72            ! angles of the far-field amplitudes in the diagnostics

   TYPE(panel_mesh) :: mesh
   TYPE(surface_geom) :: geom
   TYPE(first_order_state), ALLOCATABLE :: states(:)   ! one per heading for mean_drift = 2, otherwise only the current
   COMPLEX(8), ALLOCATABLE :: phi_wl(:,:)              ! (nwl, size(states)), total potential on the waterline segments
   COMPLEX(8), ALLOCATABLE :: rows_g(:,:,:), rows_dg(:,:,:)   ! (nbase, nsys, nwl), of the current frequency, for waterline_method = field
   LOGICAL :: multi_arrays = .FALSE.                   ! the solver's mesh is in the combined multi-body arrays
   INTEGER :: unit_far = 0, unit_near = 0, unit_terms = 0, unit_diag = 0   ! 0: not open (NEWUNIT gives negative numbers)
   INTEGER :: n_velocity_fields = 0, n_fallback_total = 0                  ! over the run, logged by so_finish

CONTAINS

   SUBROUTINE so_setup_single(opts, outputdir)
      ! Call after CalNormals.
      TYPE(so_options), INTENT(IN) :: opts
      CHARACTER(LEN=*), INTENT(IN) :: outputdir
      INTEGER, ALLOCATABLE :: body(:)

      ALLOCATE (body(NELEM))
      body = 1
      CALL expand_images(XYZ_P, GEONRM_P, DXYZ_P, DS, PNSZ, XYZ, NCN, NCON, body, RESHAPE(XR, [3, 1]), ISX, ISY, NSYS, mesh)
      multi_arrays = .FALSE.
      CALL common_setup(opts, outputdir)
   END SUBROUTINE so_setup_single

   SUBROUTINE so_setup_multi(opts, outputdir)
      ! Call after BuildCombinedArrays. MXPOT_MULTI_COMB has NELEM_PE rows: with IRSP = 1 the global TNELEM also counts
      ! the waterplane panels. Node numbers in NCON_MULTI_COMB already index XYZ_GLOBAL_MULTI_COMB.
      TYPE(so_options), INTENT(IN) :: opts
      CHARACTER(LEN=*), INTENT(IN) :: outputdir
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
      CALL common_setup(opts, outputdir)
   END SUBROUTINE so_setup_multi

   SUBROUTINE common_setup(opts, outputdir)
      TYPE(so_options), INTENT(IN) :: opts
      CHARACTER(LEN=*), INTENT(IN) :: outputdir
      CHARACTER(LEN=:), ALLOCATABLE :: msg
      LOGICAL :: ok

      IF (ABS(AMP - 1D0) > 1D-12) THEN
         WRITE (*, *) 'Mean drift needs first-order results per unit wave amplitude, but AMP = ', AMP
         ERROR STOP 1
      END IF

      ALLOCATE (states(MERGE(NBETA, 1, opts%mean_drift == 2)))
      ! Only the near-field method uses the surface geometry, so only it can stop on the geometry's errors
      IF (opts%do_near) THEN
         CALL build_surface_geom(mesh, opts%wl_z_tolerance, DEFAULT_SHARP_DEG, geom, ok, msg)
         IF (.NOT. ok) THEN
            WRITE (*, *) 'Mean drift: ', msg
            ERROR STOP 1
         END IF
         IF (LEN_TRIM(msg) > 0) WRITE (*, *) 'Mean drift warning: ', msg
         ALLOCATE (phi_wl(geom%nwl, SIZE(states)))
         IF (opts%waterline_method == 2) THEN
            ALLOCATE (rows_g(mesh%nbase, mesh%nsys, geom%nwl), rows_dg(mesh%nbase, mesh%nsys, geom%nwl))
            IF (H >= 0D0) WRITE (*, *) 'Mean drift warning: with waterline_method = field in finite depth, the Green function is', &
               ' erratic near the free surface at a small horizontal distance, so term I can jump with the mesh or the offset.'
         END IF
         OPEN (NEWUNIT=unit_near, FILE=outputdir//'/Wamit_format/MeanDriftPressure.9', STATUS='REPLACE', ACTION='WRITE')
         IF (opts%write_terms >= 1) THEN
            OPEN (NEWUNIT=unit_terms, FILE=outputdir//'/Wamit_format/MeanDriftTerms.txt', STATUS='REPLACE', ACTION='WRITE')
            WRITE (unit_terms, '(a)') '# PER BETA1 BETA2 body I term Re Im, scaled and conjugated as MeanDriftPressure.9'
            WRITE (unit_terms, '(a)') '# term: 1 waterline, 2 flare, 3 velocity, 4 displacement, 5 rotation, 6 hydrostatic'
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
      END IF
   END SUBROUTINE common_setup

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
      INTEGER :: slot, n_fallback

      slot = MERGE(ib, 1, opts%mean_drift == 2)
      CALL build_state(mesh, W1, WK, V, H, 1D0, BETA, G, RHO, mxpot_in, dspl_in, opts%motion_fixed, states(slot))
      IF (opts%do_near) THEN
         CALL velocity_differencing(mesh, geom, states(slot), n_fallback)
         n_velocity_fields = n_velocity_fields + 1
         n_fallback_total = n_fallback_total + n_fallback
         IF (opts%waterline_method == 2) THEN
            ! The headings of a frequency come in order, and the rows hold for all of them
            IF (ib == 1) CALL waterline_rows(mesh, geom, opts%waterline_offset, multi_arrays, rows_g, rows_dg)
            CALL waterline_field(mesh, geom, states(slot), rows_g, rows_dg, phi_wl(:, slot))
         ELSE
            CALL waterline_extrapolate(mesh, geom, states(slot), phi_wl(:, slot))
         END IF
      END IF

      IF (opts%write_terms == 2) THEN
         IF (ib == 1) CALL write_wave_and_modes(states(slot), mxpot_in)
         CALL write_heading_diagnostics(opts, states(slot), mxpot_in, exfc_in, mass, added, damping, stiffness)
      END IF

      IF (opts%do_far .AND. opts%mean_drift == 1) THEN
         CALL far_field_pair(mesh, states(slot), states(slot), opts, f)
         CALL write_far_rows(WVHD(ib), WVHD(ib), f)
      END IF
      IF (opts%do_near .AND. opts%mean_drift == 1) THEN
         IF (opts%do_far) THEN
            CALL near_field_output(slot, slot, WVHD(ib), WVHD(ib), f)
         ELSE
            CALL near_field_output(slot, slot, WVHD(ib), WVHD(ib))
         END IF
      END IF
   END SUBROUTINE after_motion

   SUBROUTINE so_end_frequency(opts, kk)
      ! Call after the heading loop of frequency kk. Evaluates the heading pairs when mean_drift = 2.
      TYPE(so_options), INTENT(IN) :: opts
      INTEGER, INTENT(IN) :: kk
      COMPLEX(8), ALLOCATABLE :: f(:, :, :)
      INTEGER :: i, j

      IF (is_limit_frequency() .OR. opts%mean_drift /= 2) RETURN
      IF (opts%do_far) THEN
         ALLOCATE (f(3, NBETA, NBETA))
         CALL far_field_pairs(mesh, states, opts, f)
         DO j = 1, NBETA
            DO i = 1, NBETA
               CALL write_far_rows(WVHD(i), WVHD(j), f(:, i, j))
            END DO
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
   END SUBROUTINE so_end_frequency

   SUBROUTINE so_finish()
      ! Closes the result files and logs what the velocity differencing could not do. Call after the frequency loop.
      IF (unit_near /= 0) WRITE (*, '(a,I0,a,I0,a)') ' Mean drift velocity fallback: ', n_fallback_total, &
         ' panel values over ', n_velocity_fields, ' first-order solutions have only the normal part of the velocity'
      IF (unit_far /= 0) CLOSE (unit_far)
      IF (unit_near /= 0) CLOSE (unit_near)
      IF (unit_terms /= 0) CLOSE (unit_terms)
      IF (unit_diag /= 0) CLOSE (unit_diag)
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
         CALL log_near_far('F_x', beta1_deg, far(1), near(1))
         CALL log_near_far('F_y', beta1_deg, far(2), near(2))
         CALL log_near_far('M_z', beta1_deg, far(3), near(3))
      END IF
   END SUBROUTINE near_field_output

   SUBROUTINE log_near_far(label, beta_deg, far, near)
      ! The relative difference has no meaning for a load that is zero by symmetry, so the values are logged with it.
      CHARACTER(LEN=*), INTENT(IN) :: label
      REAL(8), INTENT(IN) :: beta_deg
      COMPLEX(8), INTENT(IN) :: far, near

      WRITE (*, '(a,a,a,F10.5,a,F8.3,a,ES12.4,a,ES12.4,a,ES10.3)') ' Mean drift near vs far ', label, ': omega', W1, &
         ' beta', beta_deg, ' far', REAL(far, 8), ' near', REAL(near, 8), ' relative difference', &
         ABS(far - near)/MAX(ABS(far), TINY(1D0))
   END SUBROUTINE log_near_far

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
   END SUBROUTINE write_heading_diagnostics

   PURE REAL(8) FUNCTION norm_c(z)
      COMPLEX(8), INTENT(IN) :: z(:)
      norm_c = SQRT(SUM(REAL(z, 8)**2 + AIMAG(z)**2))
   END FUNCTION norm_c

END MODULE SecondOrderDriver
