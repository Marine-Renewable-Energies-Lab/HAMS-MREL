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

   USE SecondOrderTypes, ONLY: so_options, panel_mesh, first_order_state, surface_geom
   USE NearFieldDrift, ONLY: near_field_pair, KERNEL_DIFF, KERNEL_SUM
   USE SecondOrderOutput, ONLY: write_qtf_row
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: qtf_init, qtf_keep, qtf_pair_load, qtf_finish, qtf_stored_mb, qtf_kept

   LOGICAL :: active = .FALSE.
   INTEGER :: kind_qtf = 0, nbeta = 0, nkept = 0, all_headings_pairs = 0
   REAL(8) :: max_dw = 0D0
   CHARACTER(LEN=:), ALLOCATABLE :: outdir
   TYPE(first_order_state), ALLOCATABLE :: store(:, :)   ! (nfreq, nbeta)
   COMPLEX(8), ALLOCATABLE :: store_wl(:, :, :)          ! (nwl, nfreq, nbeta)
   REAL(8), ALLOCATABLE :: store_oufr(:), beta_deg(:)    ! (nfreq), (nbeta)

CONTAINS

   SUBROUTINE qtf_init(opts, nfreq, beta_in_deg, nwl, n_panels, nbody, outputdir)
      ! Allocates the whole store for up to nfreq frequencies, the headings beta_in_deg, n_panels panels and nbody bodies,
      ! so that a run that does not fit stops here, before the first-order solves. outputdir holds Wamit_format.
      TYPE(so_options), INTENT(IN) :: opts
      INTEGER, INTENT(IN) :: nfreq, nwl, n_panels, nbody
      REAL(8), INTENT(IN) :: beta_in_deg(:)
      CHARACTER(LEN=*), INTENT(IN) :: outputdir
      INTEGER :: st, i, b

      active = opts%qtf > 0
      IF (.NOT. active) RETURN
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
   END SUBROUTINE qtf_init

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

   REAL(8) FUNCTION qtf_stored_mb(n_panels, nwl, nbody, nfreq, nb)
      ! Size of the store in MB for nfreq frequencies and nb headings: 8 complex values per panel, 6 per body, 1 per segment
      INTEGER, INTENT(IN) :: n_panels, nwl, nbody, nfreq, nb

      qtf_stored_mb = 16D0*REAL(8*n_panels + 6*nbody + nwl, 8)*REAL(nfreq, 8)*REAL(nb, 8)/1D6
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
      IF (kind_qtf == 1 .OR. kind_qtf == 3) CALL write_file('/Wamit_format/QTF_Quadratic.10d', KERNEL_DIFF)
      IF (kind_qtf == 2 .OR. kind_qtf == 3) CALL write_file('/Wamit_format/QTF_Quadratic.10s', KERNEL_SUM)
      WRITE (*, '(a,I0,a,I0,a)') ' QTF quadratic part: ', nkept, ' frequencies and ', nbeta, ' headings'
      DEALLOCATE (store, store_wl, store_oufr)
      active = .FALSE.

   CONTAINS

      SUBROUTINE write_file(name, kernel)
         CHARACTER(LEN=*), INTENT(IN) :: name
         INTEGER, INTENT(IN) :: kernel
         COMPLEX(8) :: total(6, mesh%nbody)
         INTEGER :: unit, b1, b2, i, j, b, c

         OPEN (NEWUNIT=unit, FILE=outdir//name, STATUS='REPLACE', ACTION='WRITE')
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
                        END DO
                     END DO
                  END DO
               END DO
            END DO
         END DO
         CLOSE (unit)
      END SUBROUTINE write_file

   END SUBROUTINE qtf_finish

   SUBROUTINE out_of_memory()
      WRITE (*, '(a)') ' QTF: not enough memory to keep the first-order states of every frequency. Use fewer frequencies,'
      WRITE (*, '(a)') ' headings or panels, or qtf = 0.'
      ERROR STOP 1
   END SUBROUTINE out_of_memory

END MODULE QTFDriver
