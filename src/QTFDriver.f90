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
   ! With components a_i of amplitude and the time factor e^{-i omega t}, the second-order load is (spec 5.4)
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
   ! (QTFPotential); the free-surface integral of the second-order potential is omitted (Phase 6). Same rows, order and
   ! scaling as .10d, qtf_max_dw included. On the exact diagonal (equal frequencies, also with two headings) T_pot = 0, so
   ! there .11d is .10d. T_pot is computed Omega by Omega (potential_part), so that grad psi and psi at the waterline are
   ! built once per difference frequency, and kept per pair in tpot until the files are written. With write_terms >= 1 the
   ! five terms of each pair go to QTF_PotentialTerms.txt as they are computed: its rows are in ascending |Omega|, not in
   ! the order of .11d, so a reader matches them by their keys PER1 PER2 BETA1 BETA2 I (header in the file).
   !
   ! With qtf_b2 = control the surface piece of B2 is moved to the control surface C and the free surface F of each body
   ! (QTFPotential). The driver then hands over the points of C and F (qtf_set_control_surface), the total velocity there of
   ! every kept frequency and heading (qtf_keep_cs_velocity, indexed like store) and grad psi_k there at every difference
   ! frequency (qtf_keep_cs_psi, indexed like psi_store). With qtf_edge_weights the panel weights of the hull term at sharp
   ! edges come with qtf_set_edge_weights.
   !
   ! With qtf_direct = 1 (Phase 7a) the potential part is also found by the direct method: phi_D is solved on the body at
   ! each difference frequency (QTFDirect), in the loop of radiation-only solves after the frequency loop, where the store
   ! is complete and the solver is at Omega. qtf_direct_at_omega puts T_dir of every pair of that Omega into tdir, which
   ! has the slots of tpot.
   !
   ! A psi with NaN or Inf (a radiation solve that failed) is kept, but counted and logged by qtf_keep_psi, and every pair
   ! whose terms are not finite is counted by potential_part. The files are still written (.10d does not use psi); the main
   ! program ends with exit code 1 when qtf_nonfinite_psi() or qtf_nonfinite_pairs() is not zero, as for first-order NaN.

   USE SecondOrderTypes, ONLY: so_options, panel_mesh, first_order_state, surface_geom, cs_points, QTF_TERMS_BODY
   USE NearFieldDrift, ONLY: near_field_pair, KERNEL_DIFF, KERNEL_SUM
   USE SecondOrderOutput, ONLY: write_qtf_row, write_qtf_terms_row
   USE QTFPotential, ONLY: potential_pair, psi_surface_gradient, psi_waterline, NPOT_TERMS
   USE QTFDirect, ONLY: direct_at_omega
   USE FirstOrderState, ONLY: wave_gravity
   USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qtf_init, qtf_keep, qtf_pair_load, qtf_finish, qtf_stored_mb, qtf_kept
   PUBLIC :: qtf_radiation_list, qtf_radiation_omegas, qtf_direct_omegas, qtf_omega_on_grid, qtf_keep_psi, qtf_psi_index, qtf_psi_count, qtf_psi_tolerance
   PUBLIC :: qtf_nonfinite_psi, qtf_nonfinite_pairs
   PUBLIC :: qtf_direct_at_omega, qtf_direct_pairs
   PUBLIC :: qtf_set_control_surface, qtf_keep_cs_velocity, qtf_keep_cs_psi, qtf_set_edge_weights

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
   COMPLEX(8), ALLOCATABLE :: psi_store(:, :, :)         ! (n_panels, 6*nbody, nneed)
   LOGICAL, ALLOCATABLE :: psi_kept(:)                   ! (nneed)
   ! Potential part of the difference QTF (qtf_terms = body)
   LOGICAL :: body_active = .FALSE.                      ! QTF_Indirect.11d is written
   LOGICAL :: pot_terms = .FALSE.                        ! QTF_PotentialTerms.txt is written (write_terms >= 1)
   COMPLEX(8), ALLOCATABLE :: tpot(:, :, :, :)           ! (6*nbody, nfreq, nfreq, heading pairs), e^{-i omega t}
   ! Direct method (qtf_direct = 1, Phase 7a): T_dir of each pair, in the slots of tpot, filled Omega by Omega by
   ! qtf_direct_at_omega during the radiation-only loop after the frequency loop, when the store is complete
   LOGICAL :: direct_active = .FALSE.
   COMPLEX(8), ALLOCATABLE :: tdir(:, :, :, :)           ! (6*nbody, nfreq, nfreq, heading pairs), e^{-i omega t}
   LOGICAL, ALLOCATABLE :: tdir_kept(:)                  ! (nneed), the direct solve has run at psi_need(m)
   ! qtf_b2 = control: points of C and F of every body, and the fields there
   LOGICAL :: b2_control = .FALSE.
   TYPE(cs_points) :: cs
   COMPLEX(8), ALLOCATABLE :: cs_u(:, :, :, :)           ! (3, cs%n, nfreq, nbeta), total velocity, as store
   LOGICAL, ALLOCATABLE :: cs_u_kept(:, :)               ! (nfreq, nbeta)
   COMPLEX(8), ALLOCATABLE :: cs_gpsi(:, :, :, :)        ! (3, cs%n, 6*nbody, nneed), grad psi_k, as psi_store
   LOGICAL, ALLOCATABLE :: cs_psi_kept(:)                ! (nneed)
   REAL(8), ALLOCATABLE :: edge_w(:)                     ! (n_panels), weights of the B2 hull term at sharp edges, if set
   ! Non-finite results of the potential part, over the run (reset by qtf_init, kept after qtf_finish for the exit code)
   INTEGER :: n_bad_psi = 0                              ! difference frequencies whose psi has NaN or Inf
   INTEGER :: n_bad_pairs = 0                            ! frequency and heading pairs whose T_pot has NaN or Inf

CONTAINS

   SUBROUTINE qtf_init(opts, nfreq, beta_in_deg, nwl, n_panels, nbody, outputdir, omegas_main)
      ! Allocates the whole store for up to nfreq frequencies, the headings beta_in_deg, n_panels panels and nbody bodies,
      ! so that a run that does not fit stops here, before the first-order solves. outputdir holds Wamit_format.
      ! omegas_main (rad/s), the nfreq frequencies of the QTF grid, is needed for qtf_terms = body: the store then also holds
      ! psi at every difference frequency.
      TYPE(so_options), INTENT(IN) :: opts
      INTEGER, INTENT(IN) :: nfreq, nwl, n_panels, nbody
      REAL(8), INTENT(IN) :: beta_in_deg(:)
      CHARACTER(LEN=*), INTENT(IN) :: outputdir
      REAL(8), INTENT(IN), OPTIONAL :: omegas_main(:)
      INTEGER :: st, i, b

      active = opts%qtf > 0
      IF (.NOT. active) RETURN
      body_active = opts%qtf_terms == QTF_TERMS_BODY .AND. (opts%qtf == 1 .OR. opts%qtf == 3)
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
         CALL qtf_radiation_list(omegas_main, opts%qtf_max_dw, psi_tol, psi_need, psi_reuse)
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
      END IF
      IF (ALLOCATED(psi_store)) DEALLOCATE (psi_store)
      IF (ALLOCATED(psi_kept)) DEALLOCATE (psi_kept)
      IF (ALLOCATED(tpot)) DEALLOCATE (tpot)
      IF (ALLOCATED(tdir)) DEALLOCATE (tdir)
      IF (ALLOCATED(tdir_kept)) DEALLOCATE (tdir_kept)
      CALL free_control_surface()
      n_bad_psi = 0
      n_bad_pairs = 0
      pot_terms = body_active .AND. opts%write_terms >= 1
      IF (body_active) THEN
         ALLOCATE (psi_store(n_panels, 6*nbody, SIZE(psi_need)), psi_kept(SIZE(psi_need)), STAT=st)
         IF (st /= 0) CALL out_of_memory(body=.TRUE.)
         psi_store = (0D0, 0D0)
         psi_kept = .FALSE.
         ALLOCATE (tpot(6*nbody, nfreq, nfreq, heading_pairs(nbeta)), STAT=st)
         IF (st /= 0) CALL out_of_memory(body=.TRUE.)
         tpot = (0D0, 0D0)
      END IF
      direct_active = body_active .AND. opts%qtf_direct == 1
      IF (direct_active) THEN
         ALLOCATE (tdir(6*nbody, nfreq, nfreq, heading_pairs(nbeta)), tdir_kept(SIZE(psi_need)), STAT=st)
         IF (st /= 0) CALL out_of_memory(body=.TRUE.)
         tdir = (0D0, 0D0)
         tdir_kept = .FALSE.
      END IF

   CONTAINS

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

   END SUBROUTINE qtf_init

   REAL(8) PURE FUNCTION qtf_psi_tolerance(omegas_main)
      ! The tolerance of the difference frequencies of the grid omegas_main (rad/s): qtf_radiation_list puts every Omega
      ! within it of its psi, and an Omega is found (qtf_psi_index) and accepted (potential_pair) within it. Rounding only
      REAL(8), INTENT(IN) :: omegas_main(:)

      qtf_psi_tolerance = PSI_REL_TOL*MAXVAL(omegas_main)
   END FUNCTION qtf_psi_tolerance

   SUBROUTINE qtf_set_control_surface(pts)
      ! qtf_b2 = control: keeps the points of C and F and allocates the fields there for the store of qtf_init. Call after
      ! qtf_init, before the frequency loop, only with qtf_terms = body.
      TYPE(cs_points), INTENT(IN) :: pts
      INTEGER :: st

      IF (.NOT. body_active) THEN
         WRITE (*, '(a)') ' QTF: qtf_set_control_surface needs qtf_terms = body and the difference QTF.'
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
      ! The distinct difference frequencies Omega = omega_i - omega_j > 0 of the grid omegas_main, ascending, with
      ! Omega <= max_dw_in (all of them for max_dw_in = 0), as in qtf_finish. Every difference is within tol of its need(m):
      ! the differences within tol/2 of the smallest of them are one, and reuse_index(m) is the main-grid index whose
      ! frequency is need(m) within tol/2, else 0; need(m) is then that frequency exactly, so that the radiation solve of the
      ! main grid finds it. qtf_init passes qtf_psi_tolerance, which merges rounding only.
      REAL(8), INTENT(IN) :: omegas_main(:), max_dw_in, tol
      REAL(8), ALLOCATABLE, INTENT(OUT) :: need(:)
      INTEGER, ALLOCATABLE, INTENT(OUT) :: reuse_index(:)
      REAL(8), ALLOCATABLE :: d(:)
      REAL(8) :: dw
      INTEGER :: n, i, j, m, nd

      n = SIZE(omegas_main)
      ALLOCATE (d(n*n))
      nd = 0
      DO i = 1, n
         DO j = 1, n
            dw = omegas_main(i) - omegas_main(j)
            IF (dw <= tol) CYCLE
            IF (max_dw_in > 0D0 .AND. dw > max_dw_in*(1D0 + 1D-9)) CYCLE
            nd = nd + 1
            d(nd) = dw
         END DO
      END DO
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
   END SUBROUTINE qtf_radiation_list

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
      ! The difference frequencies that are not on the main grid: each needs a radiation-only solve. Empty unless
      ! qtf_terms = body.
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

   SUBROUTINE qtf_keep_psi(omega, psi, ok)
      ! Keeps psi(n_panels, 6*nbody) of the difference frequency omega. An omega that is not a difference frequency, or a psi
      ! of another size, stops the run; with ok present, psi is not kept and ok is false instead. A psi with NaN or Inf is
      ! kept, logged with omega and the first mode and body concerned, and counted in qtf_nonfinite_psi().
      REAL(8), INTENT(IN) :: omega
      COMPLEX(8), INTENT(IN) :: psi(:, :)
      LOGICAL, INTENT(OUT), OPTIONAL :: ok
      INTEGER :: m, nbad, k

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
         WRITE (*, '(a)') ' QTF_Indirect.11d and QTF_PotentialTerms.txt are not valid at that Omega; the run will end with'// &
            ' exit code 1.'
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
      ! Number of frequency and heading pairs of QTF_Indirect.11d whose T_pot had NaN or Inf in this run
      qtf_nonfinite_pairs = n_bad_pairs
   END FUNCTION qtf_nonfinite_pairs

   SUBROUTINE qtf_direct_pairs(omega, pairs)
      ! The pairs of the difference QTF whose |omega_i - omega_j| is the difference frequency omega (qtf_psi_index), with
      ! omega_i /= omega_j and within qtf_max_dw, as potential_part takes them: pairs(4, npair) = frequency and heading of
      ! state 1, then of state 2, in the store. Empty unless qtf_direct = 1, or before the frequencies are kept.
      REAL(8), INTENT(IN) :: omega
      INTEGER, ALLOCATABLE, INTENT(OUT) :: pairs(:, :)
      INTEGER, ALLOCATABLE :: tmp(:, :)
      REAL(8) :: dw
      INTEGER :: m, i, j, b1, b2, np

      ALLOCATE (pairs(4, 0))
      IF (.NOT. direct_active) RETURN
      m = qtf_psi_index(omega)
      IF (m == 0) RETURN
      ALLOCATE (tmp(4, nkept*nkept*heading_pairs(nbeta)))
      np = 0
      DO b1 = 1, nbeta
         DO b2 = 1, nbeta
            IF (all_headings_pairs == 0 .AND. b1 /= b2) CYCLE
            DO i = 1, nkept
               DO j = 1, nkept
                  ! Heading 1 stands for all, as in potential_part
                  dw = store(i, 1)%omega - store(j, 1)%omega
                  IF (dw == 0D0) CYCLE
                  IF (max_dw > 0D0 .AND. ABS(dw) > max_dw*(1D0 + 1D-9)) CYCLE
                  IF (qtf_psi_index(ABS(dw)) /= m) CYCLE
                  np = np + 1
                  tmp(:, np) = [i, b1, j, b2]
               END DO
            END DO
         END DO
      END DO
      pairs = tmp(:, 1:np)
   END SUBROUTINE qtf_direct_pairs

   SUBROUTINE qtf_direct_at_omega(mesh, geom, omega, npair, seconds)
      ! The direct method at the difference frequency omega (qtf_direct = 1): T_dir of every pair with that Omega
      ! (qtf_direct_pairs), by direct_at_omega, into tdir. Call after the frequency loop, when every frequency is kept,
      ! with the solver at omega after its radiation solve (RadiationOnlySolves). Pairs with NaN or Inf are counted in
      ! qtf_nonfinite_pairs() and logged. npair and seconds (optional): the number of pairs, and the times of
      ! direct_at_omega, for the log.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      REAL(8), INTENT(IN) :: omega
      INTEGER, INTENT(OUT), OPTIONAL :: npair
      REAL(8), INTENT(OUT), OPTIONAL :: seconds(3)
      REAL(8) :: sec(3)
      INTEGER, ALLOCATABLE :: pairs(:, :)
      COMPLEX(8), ALLOCATABLE :: t(:, :)
      INTEGER :: m, p, nbad

      IF (PRESENT(npair)) npair = 0
      IF (PRESENT(seconds)) seconds = 0D0
      IF (.NOT. direct_active) RETURN
      m = qtf_psi_index(omega)
      IF (m == 0) THEN
         WRITE (*, '(a,ES14.6,a)') ' QTF direct: omega =', omega, ' rad/s is not a difference frequency of the grid.'
         ERROR STOP 1
      END IF
      IF (nkept /= SIZE(store, 1)) THEN
         WRITE (*, '(a,I0,a,I0,a)') ' QTF direct: ', nkept, ' of ', SIZE(store, 1), &
            ' frequencies are kept; the direct solve runs after the frequency loop.'
         ERROR STOP 1
      END IF
      CALL qtf_direct_pairs(omega, pairs)
      ALLOCATE (t(6*mesh%nbody, SIZE(pairs, 2)))
      CALL direct_at_omega(mesh, geom, omega, store, store_wl, pairs, t, omega_tol=psi_tol, seconds=sec)
      IF (PRESENT(npair)) npair = SIZE(pairs, 2)
      IF (PRESENT(seconds)) seconds = sec
      nbad = 0
      DO p = 1, SIZE(pairs, 2)
         tdir(:, pairs(1, p), pairs(3, p), heading_pair(pairs(2, p), pairs(4, p))) = t(:, p)
         IF (.NOT. ALL(IEEE_IS_FINITE(REAL(t(:, p), 8)) .AND. IEEE_IS_FINITE(AIMAG(t(:, p))))) nbad = nbad + 1
      END DO
      tdir_kept(m) = .TRUE.
      IF (nbad > 0) THEN
         n_bad_pairs = n_bad_pairs + nbad
         WRITE (*, '(a,I0,a,ES14.6,a)') ' QTF error: T_dir has NaN or Inf at ', nbad, ' frequency and heading pairs with'// &
            ' Omega =', psi_need(m), ' rad/s.'
      END IF
   END SUBROUTINE qtf_direct_at_omega

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

   REAL(8) FUNCTION qtf_stored_mb(n_panels, nwl, nbody, nfreq, nb, npsi, npairs, ncs)
      ! Size of the store in MB for nfreq frequencies and nb headings: 8 complex values per panel, 6 per body, 1 per segment.
      ! With npsi, also psi at npsi difference frequencies: 6*nbody complex values per panel. With npairs heading pairs, also
      ! T_pot: 6*nbody complex values per frequency pair and heading pair. With ncs points of C and F (qtf_b2 = control, with
      ! npsi), also the velocity there (3 complex values per point, frequency and heading) and grad psi_k (3 per point, mode
      ! and difference frequency).
      INTEGER, INTENT(IN) :: n_panels, nwl, nbody, nfreq, nb
      INTEGER, INTENT(IN), OPTIONAL :: npsi, npairs, ncs

      qtf_stored_mb = 16D0*REAL(8*n_panels + 6*nbody + nwl, 8)*REAL(nfreq, 8)*REAL(nb, 8)/1D6
      IF (PRESENT(npsi)) qtf_stored_mb = qtf_stored_mb + 16D0*REAL(n_panels, 8)*REAL(6*nbody, 8)*REAL(npsi, 8)/1D6
      IF (PRESENT(npairs)) qtf_stored_mb = qtf_stored_mb + 16D0*REAL(6*nbody, 8)*REAL(nfreq, 8)**2*REAL(npairs, 8)/1D6
      IF (PRESENT(ncs)) THEN
         qtf_stored_mb = qtf_stored_mb + 16D0*3D0*REAL(ncs, 8)*REAL(nfreq, 8)*REAL(nb, 8)/1D6
         IF (PRESENT(npsi)) qtf_stored_mb = qtf_stored_mb + 16D0*3D0*REAL(ncs, 8)*REAL(6*nbody, 8)*REAL(npsi, 8)/1D6
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
      IF (body_active) CALL potential_part(mesh, geom, refl)
      IF (kind_qtf == 1 .OR. kind_qtf == 3) CALL write_file('/Wamit_format/QTF_Quadratic.10d', KERNEL_DIFF)
      IF (kind_qtf == 2 .OR. kind_qtf == 3) CALL write_file('/Wamit_format/QTF_Quadratic.10s', KERNEL_SUM)
      WRITE (*, '(a,I0,a,I0,a)') ' QTF quadratic part: ', nkept, ' frequencies and ', nbeta, ' headings'
      DEALLOCATE (store, store_wl, store_oufr)
      IF (ALLOCATED(psi_store)) DEALLOCATE (psi_store, psi_kept, psi_need, psi_reuse)
      IF (ALLOCATED(tpot)) DEALLOCATE (tpot)
      IF (ALLOCATED(tdir)) DEALLOCATE (tdir, tdir_kept)
      CALL free_control_surface()
      active = .FALSE.
      body_active = .FALSE.
      direct_active = .FALSE.

   CONTAINS

      SUBROUTINE write_file(name, kernel)
         ! With qtf_terms = body, the difference kernel also writes QTF_Indirect.11d, row for row with .10d
         CHARACTER(LEN=*), INTENT(IN) :: name
         INTEGER, INTENT(IN) :: kernel
         COMPLEX(8) :: total(6, mesh%nbody), value
         INTEGER :: unit, unit_indirect, b1, b2, i, j, b, c

         OPEN (NEWUNIT=unit, FILE=outdir//name, STATUS='REPLACE', ACTION='WRITE')
         unit_indirect = 0
         IF (kernel == KERNEL_DIFF .AND. body_active) THEN
            OPEN (NEWUNIT=unit_indirect, FILE=outdir//'/Wamit_format/QTF_Indirect.11d', STATUS='REPLACE', ACTION='WRITE')
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
                     DO b = 1, mesh%nbody
                        DO c = 1, 6
                           CALL write_qtf_row(unit, store_oufr(i), store_oufr(j), beta_deg(b1), beta_deg(b2), &
                                              6*(b - 1) + c, total(c, b), refl, c > 3)
                           IF (unit_indirect == 0) CYCLE
                           ! On the diagonal T_pot = 0 and .11d is .10d to the bit (adding a zero could flip a -0)
                           value = total(c, b)
                           IF (store(i, b1)%omega /= store(j, b2)%omega) THEN
                              value = value + tpot(6*(b - 1) + c, i, j, heading_pair(b1, b2))
                           END IF
                           CALL write_qtf_row(unit_indirect, store_oufr(i), store_oufr(j), beta_deg(b1), beta_deg(b2), &
                                              6*(b - 1) + c, value, refl, c > 3)
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
         CLOSE (unit)
         IF (unit_indirect /= 0) CLOSE (unit_indirect)
      END SUBROUTINE write_file

   END SUBROUTINE qtf_finish

   SUBROUTINE potential_part(mesh, geom, refl)
      ! T_pot of every pair of the difference QTF with omega_i /= omega_j within qtf_max_dw, into tpot, Omega by Omega:
      ! grad psi and psi at the waterline are built once for each difference frequency. With write_terms >= 1 the terms go to
      ! QTF_PotentialTerms.txt, in that order (ascending |Omega|, not the .11d order). Stops if psi of a needed Omega was not
      ! stored. Pairs with NaN or Inf terms are written as they are, counted in qtf_nonfinite_pairs() and logged per Omega.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(surface_geom), INTENT(IN) :: geom
      REAL(8), INTENT(IN) :: refl
      INTEGER, ALLOCATABLE :: pair_m(:, :)                 ! (nkept, nkept), index of |omega_i - omega_j| in psi_need, or 0
      COMPLEX(8), ALLOCATABLE :: grad(:, :, :), pwl(:, :), terms(:, :)
      REAL(8) :: dw, g_wave
      INTEGER :: unit, i, j, m, b1, b2, b, c, k, nmode, nbad

      nmode = 6*mesh%nbody
      ALLOCATE (pair_m(nkept, nkept))
      pair_m = 0
      DO i = 1, nkept
         DO j = 1, nkept
            ! Heading 1 stands for all: qtf_keep stores every heading of a frequency from the same solve, so
            ! store(i, b)%omega is the same for every b (the grid is shared by construction)
            dw = store(i, 1)%omega - store(j, 1)%omega
            IF (dw == 0D0) CYCLE
            IF (max_dw > 0D0 .AND. ABS(dw) > max_dw*(1D0 + 1D-9)) CYCLE
            m = qtf_psi_index(ABS(dw))
            IF (m == 0) THEN
               WRITE (*, '(a,ES14.6,a)') ' QTF: Omega =', ABS(dw), ' rad/s is not among the difference frequencies of psi.'
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
            pair_m(i, j) = m
         END DO
      END DO

      unit = 0
      IF (pot_terms) THEN
         OPEN (NEWUNIT=unit, FILE=outdir//'/Wamit_format/QTF_PotentialTerms.txt', STATUS='REPLACE', ACTION='WRITE')
         WRITE (unit, '(a)') '# Body-forcing terms of the potential part of QTF_Indirect.11d (qtf_terms = body), scaled'// &
            ' and conjugated as QTF_Quadratic.10d.'
         WRITE (unit, '(a)') '# PER1 PER2 BETA1 BETA2 I, then Re Im of each term: 1 Froude-Krylov of phi_I2, 2 incident'// &
            ' diffraction, 3 B1, 4 B2, 5 B3.'
         WRITE (unit, '(a)') '# .11d - .10d is their sum; the free-surface integral is omitted. One row per pair with'// &
            ' omega_i /= omega_j, by ascending |omega_i - omega_j|.'
         WRITE (unit, '(a)') '# The rows are not in the order of .11d: match them by PER1 PER2 BETA1 BETA2 I.'
      END IF
      ALLOCATE (grad(3, mesh%n, nmode), pwl(geom%nwl, nmode), terms(nmode, NPOT_TERMS))
      DO m = 1, SIZE(psi_need)
         IF (.NOT. ANY(pair_m == m)) CYCLE
         CALL psi_surface_gradient(mesh, geom, psi_store(:, :, m), grad)
         CALL psi_waterline(mesh, geom, psi_store(:, :, m), grad, pwl)
         nbad = 0
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
                                               psi_store(:, :, m), grad, pwl, psi_need(m), s1%rho, g_wave, terms, &
                                               cs=cs, grad_psi_cs=cs_gpsi(:, :, :, m), u1_cs=cs_u(:, :, i, b1), &
                                               u2_cs=cs_u(:, :, j, b2), edge_w=edge_w, omega_tol=psi_tol)
                        ELSE IF (b2_control) THEN
                           CALL potential_pair(mesh, geom, s1, s2, store_wl(:, i, b1), store_wl(:, j, b2), &
                                               psi_store(:, :, m), grad, pwl, psi_need(m), s1%rho, g_wave, terms, &
                                               cs=cs, grad_psi_cs=cs_gpsi(:, :, :, m), u1_cs=cs_u(:, :, i, b1), &
                                               u2_cs=cs_u(:, :, j, b2), omega_tol=psi_tol)
                        ELSE
                           CALL potential_pair(mesh, geom, s1, s2, store_wl(:, i, b1), store_wl(:, j, b2), &
                                               psi_store(:, :, m), grad, pwl, psi_need(m), s1%rho, g_wave, terms, &
                                               omega_tol=psi_tol)
                        END IF
                     END ASSOCIATE
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
         IF (nbad > 0) THEN
            n_bad_pairs = n_bad_pairs + nbad
            WRITE (*, '(a,I0,a,ES14.6,a)') ' QTF error: T_pot has NaN or Inf at ', nbad, ' frequency and heading pairs with'// &
               ' Omega =', psi_need(m), ' rad/s; those rows of QTF_Indirect.11d are not valid.'
         END IF
      END DO
      IF (unit /= 0) CLOSE (unit)
   END SUBROUTINE potential_part

   SUBROUTINE out_of_memory(body)
      ! body: the store of qtf_terms = body (psi, T_pot, the fields on C and F), whose size the log gave before
      LOGICAL, INTENT(IN), OPTIONAL :: body

      WRITE (*, '(a)') ' QTF: not enough memory to keep the first-order states of every frequency. Use fewer frequencies,'
      WRITE (*, '(a)') ' headings or panels, or qtf = 0.'
      IF (PRESENT(body)) THEN
         IF (body) THEN
            WRITE (*, '(a)') ' With qtf_terms = body the store also holds psi at every difference frequency (count in the log'
            WRITE (*, '(a)') ' above): a smaller qtf_max_dw, or a uniform omega grid (whose differences coincide), needs fewer.'
            WRITE (*, '(a)') ' With qtf_b2 = control the fields on C and F add to it: qtf_b2 = body avoids them, a larger'
            WRITE (*, '(a)') ' cs_panel_size reduces them.'
         END IF
      END IF
      ERROR STOP 1
   END SUBROUTINE out_of_memory

END MODULE QTFDriver
