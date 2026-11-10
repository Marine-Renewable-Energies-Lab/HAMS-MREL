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
MODULE QTFField
   ! Second-order pressure and free-surface elevation of a pair of first-order waves at a field point: the quadratic part
   ! (products of first-order quantities) plus the second-order incident potential. The second-order diffracted and radiated
   ! potentials are omitted, so near a body the result is partial.
   !
   ! Bernoulli's equation expanded about z = 0 (Dean and Dalrymple 1991), time factor e^{-i omega t}:
   !   pressure at a fixed point below z = 0:  p2 = -rho [phi2_t + 1/2 |grad phi1|^2]
   !   elevation on z = 0:  eta2 = -(1/g) [phi2_t + 1/2 |grad phi1|^2 + eta1 phi1_zt],
   !   with eta1 = (i omega/g) phi1 and phi1_zt = -i omega phi1_z.
   ! For the ordered pair (i, j) and the kernel P of NearFieldDrift (pair_minus: difference, pair_plus: sum), per unit
   ! amplitudes:
   !   quad_p   = -(rho/2) sum_m P(u_im, u_im, u_jm, u_jm)
   !   quad_eta = -(1/g) [ 1/2 sum_m P(u_im, u_im, u_jm, u_jm) + P(eta_i, -i omega_i u_i3, eta_j, -i omega_j u_j3) ]
   !   pot_p    = i rho Omega phi2,   pot_eta = (i Omega/g) phi2,
   ! with Omega = omega_i -/+ omega_j and phi2 the second-order incident potential of the pair (SecondOrderIncident). The
   ! physical e^{-i Omega t} amplitude of two waves a_i /= a_j is 2 a_i a_j (quad + pot); for the sum with i = j it is
   ! a_i^2 (quad + pot).
   !
   ! The store and the file (qtf_field = 1 in SecondOrder.in). At the end of each frequency of the QTF grid the driver gives
   ! qfield_keep the first-order potential and gradient of every heading at the field points (fields_green of
   ! FieldVelocity); these are all that the pairs need, so the store is small: 4 complex values per point, frequency and
   ! heading. After the frequency loop qfield_finish forms every pair of the QTFs of QTFDriver, with the same kernels
   ! (qtf = 1 difference, 2 sum, 3 both), heading pairs (mean_drift = 2: all, else beta_1 = beta_2) and qtf_max_dw on the
   ! difference pairs, so that its rows are those of QTF_Quadratic.10d and .10s at each point.

   USE Const_mod, ONLY: G
   USE SecondOrderIncident, ONLY: incident2, incident2_pair, incident2_phi
   USE NearFieldDrift, ONLY: pair_product, KERNEL_DIFF, KERNEL_SUM
   USE SecondOrderTypes, ONLY: so_options, first_order_state, panel_mesh, cs_points, cs_fields, CS_FREE
   USE SecondOrderOutput, ONLY: write_field_row
   USE FirstOrderState, ONLY: wave_gravity
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: field_pair, build_field_points
   PUBLIC :: qfield_init, qfield_keep, qfield_finish, qfield_stored_mb, FIELD_FD_STEP

   ! Central-difference step of the first-order gradient at the field points, as a fraction of the size of the nearest hull
   ! panel: the default of cs_fd_step, which fields_green takes for the control surface
   REAL(8), PARAMETER :: FIELD_FD_STEP = 0.02D0

   LOGICAL :: active = .FALSE.
   INTEGER :: kind_qtf = 0, nbeta = 0, nkept = 0, npt = 0, all_headings_pairs = 0
   REAL(8) :: max_dw = 0D0
   CHARACTER(LEN=:), ALLOCATABLE :: outdir
   REAL(8), ALLOCATABLE :: px(:, :)                     ! (3, npt), z = 0 at the surface points
   LOGICAL, ALLOCATABLE :: is_surface(:)                ! (npt), elevation (z = 0) or pressure
   COMPLEX(8), ALLOCATABLE :: phi_store(:, :, :)        ! (npt, nfreq, nbeta), total first-order potential
   COMPLEX(8), ALLOCATABLE :: grad_store(:, :, :, :)    ! (3, npt, nfreq, nbeta), its gradient
   REAL(8), ALLOCATABLE :: w_store(:, :), k_store(:, :), h_store(:, :), b_store(:, :), rho_store(:, :)   ! (nfreq, nbeta)
   REAL(8), ALLOCATABLE :: oufr_store(:), beta_deg(:)   ! (nfreq), (nbeta)

CONTAINS

   PURE SUBROUTINE field_pair(kernel, w1, phi1, grad1, w2, phi2, grad2, inc, x, surface, rho, g, quad, pot)
      ! Quadratic part quad and second-order incident part pot of the pressure (surface = .FALSE.) or the elevation
      ! (surface = .TRUE., x(3) = 0) of the pair of unit waves with frequencies w1, w2, first-order potentials phi and
      ! gradients grad at x. inc is the second-order incident potential of the pair, made with the same kernel.
      INTEGER, INTENT(IN) :: kernel
      REAL(8), INTENT(IN) :: w1, w2, x(3), rho, g
      COMPLEX(8), INTENT(IN) :: phi1, grad1(3), phi2, grad2(3)
      TYPE(incident2), INTENT(IN) :: inc
      LOGICAL, INTENT(IN) :: surface
      COMPLEX(8), INTENT(OUT) :: quad, pot
      COMPLEX(8) :: speed, p2, eta1, eta2
      COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
      INTEGER :: m

      speed = (0D0, 0D0)
      DO m = 1, 3
         speed = speed + pair_product(kernel, grad1(m), grad1(m), grad2(m), grad2(m))
      END DO
      p2 = incident2_phi(inc, x)
      IF (surface) THEN
         eta1 = CI*w1/g*phi1
         eta2 = CI*w2/g*phi2
         quad = -(0.5D0*speed + pair_product(kernel, eta1, -CI*w1*grad1(3), eta2, -CI*w2*grad2(3)))/g
         pot = CI*inc%omega/g*p2
      ELSE
         quad = -0.5D0*rho*speed
         pot = CI*rho*inc%omega*p2
      END IF
   END SUBROUTINE field_pair

   SUBROUTINE build_field_points(mesh, xfp, pts, surface, ok, msg)
      ! Field points xfp(3,:) as a cs_points set (kind CS_FREE, nrm = 0, w = 0). Type as PressureElevation: |z| <= 1e-6 is
      ! a surface (elevation) point and gets z = 0, z < -1e-6 is a pressure point. psize and body are those of the nearest
      ! hull panel (centroids of the full panel set, so symmetry images are included). Rejected, checked in point order,
      ! the first failing point giving msg: above the free surface; within one local panel size of a hull panel centroid
      ! (a point has no panel size of its own, so the limit is the hull panel's size); inside a body, decided by the nearest centroid c as (x - c).n > 0, with n the
      ! panel normal, which points into the body.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: xfp(:,:)
      TYPE(cs_points), INTENT(OUT) :: pts
      LOGICAL, ALLOCATABLE, INTENT(OUT) :: surface(:)
      LOGICAL, INTENT(OUT) :: ok
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: msg
      REAL(8), PARAMETER :: ZTOL = 1D-6
      REAL(8) :: x(3), d2, d2min
      INTEGER :: ip, i, imin, ihull, n
      CHARACTER(LEN=300) :: line

      n = SIZE(xfp, 2)
      pts%n = n
      ALLOCATE (pts%x(3, n), pts%nrm(3, n), pts%w(n), pts%psize(n), pts%kind(n), pts%body(n), surface(n))
      pts%x = xfp
      pts%nrm = 0D0
      pts%w = 0D0
      pts%psize = 0D0
      pts%kind = CS_FREE
      pts%body = 0
      surface = .FALSE.
      ok = .TRUE.
      msg = ''

      DO ip = 1, n
         x = xfp(:, ip)
         IF (x(3) > ZTOL) THEN
            WRITE (line, '(a,i0,a,3(1x,es11.3),a)') 'field point ', ip, ' at', x, ' is above the free surface (z > 0)'
            msg = TRIM(line)
            ok = .FALSE.
            RETURN
         END IF
         IF (ABS(x(3)) <= ZTOL) THEN
            surface(ip) = .TRUE.
            x(3) = 0D0
            pts%x(3, ip) = 0D0
         END IF

         imin = 0
         ihull = 0
         d2min = HUGE(1D0)
         DO i = 1, mesh%n
            d2 = SUM((x - mesh%xc(:, i))**2)
            IF (d2 < d2min) THEN
               d2min = d2
               imin = i
            END IF
            IF (ihull == 0 .AND. d2 < mesh%psize(i)**2) ihull = i
         END DO
         IF (imin == 0) CYCLE
         pts%psize(ip) = mesh%psize(imin)
         pts%body(ip) = mesh%body(imin)

         IF (ihull > 0) THEN
            WRITE (line, '(a,i0,a,3(1x,es11.3),a,i0,a,i0,a)') 'field point ', ip, ' at', x, ' lies on the hull of body ', &
               mesh%body(ihull), ' (within one panel size of panel ', ihull, ')'
            msg = TRIM(line)
            ok = .FALSE.
            RETURN
         END IF
         IF (SUM((x - mesh%xc(:, imin))*mesh%nrm(:, imin)) > 0D0) THEN
            WRITE (line, '(a,i0,a,3(1x,es11.3),a,i0)') 'field point ', ip, ' at', x, ' is inside body ', mesh%body(imin)
            msg = TRIM(line)
            ok = .FALSE.
            RETURN
         END IF
      END DO
   END SUBROUTINE build_field_points

   SUBROUTINE qfield_init(opts, nfreq, beta_in_deg, pts, surface, outputdir)
      ! Allocates the store for up to nfreq frequencies, the headings beta_in_deg and the field points pts (of
      ! build_field_points, with their types in surface), so that a run that does not fit stops here, before the first-order
      ! solves. outputdir holds Wamit_format.
      TYPE(so_options), INTENT(IN) :: opts
      INTEGER, INTENT(IN) :: nfreq
      REAL(8), INTENT(IN) :: beta_in_deg(:)
      TYPE(cs_points), INTENT(IN) :: pts
      LOGICAL, INTENT(IN) :: surface(:)
      CHARACTER(LEN=*), INTENT(IN) :: outputdir
      INTEGER :: st

      active = opts%qtf_field == 1 .AND. opts%qtf > 0
      IF (.NOT. active) RETURN
      kind_qtf = opts%qtf
      max_dw = opts%qtf_max_dw
      all_headings_pairs = MERGE(1, 0, opts%mean_drift == 2)
      nbeta = SIZE(beta_in_deg)
      npt = pts%n
      nkept = 0
      outdir = outputdir
      px = pts%x
      is_surface = surface
      beta_deg = beta_in_deg
      IF (ALLOCATED(phi_store)) DEALLOCATE (phi_store, grad_store, w_store, k_store, h_store, b_store, rho_store, oufr_store)
      ALLOCATE (phi_store(npt, nfreq, nbeta), grad_store(3, npt, nfreq, nbeta), w_store(nfreq, nbeta), k_store(nfreq, nbeta), &
                h_store(nfreq, nbeta), b_store(nfreq, nbeta), rho_store(nfreq, nbeta), oufr_store(nfreq), STAT=st)
      IF (st /= 0) THEN
         WRITE (*, '(a)') ' QTF field: not enough memory for the first-order field at the field points. Use fewer frequencies,'
         WRITE (*, '(a)') ' headings or field points, or qtf_field = 0.'
         ERROR STOP 1
      END IF
      ! Touch the arrays, so that memory that is only promised (overcommit) is taken now
      phi_store = (0D0, 0D0)
      grad_store = (0D0, 0D0)
      w_store = 0D0
      k_store = 0D0
      h_store = -1D0
      b_store = 0D0
      rho_store = 0D0
      oufr_store = 0D0
   END SUBROUTINE qfield_init

   REAL(8) FUNCTION qfield_stored_mb(npoints, nfreq, nb)
      ! Size of the store in MB: 4 complex values per point, frequency and heading
      INTEGER, INTENT(IN) :: npoints, nfreq, nb

      qfield_stored_mb = 64D0*REAL(npoints, 8)*REAL(nfreq, 8)*REAL(nb, 8)/1D6
   END FUNCTION qfield_stored_mb

   SUBROUTINE qfield_keep(states, fields, oufr)
      ! Keeps the first-order fields at the points of every heading of the current frequency: states(:) and fields(:) in
      ! heading order, the fields of fields_green at the frequency of the states. One call per frequency of the QTF grid.
      TYPE(first_order_state), INTENT(IN) :: states(:)
      TYPE(cs_fields), INTENT(IN) :: fields(:)
      REAL(8), INTENT(IN) :: oufr
      INTEGER :: ib

      IF (.NOT. active) RETURN
      IF (SIZE(states) /= nbeta .OR. SIZE(fields) /= nbeta) THEN
         WRITE (*, '(a)') ' QTF field: qfield_keep needs the states and fields of every heading.'
         ERROR STOP 1
      END IF
      nkept = nkept + 1
      IF (nkept > SIZE(phi_store, 2)) THEN
         WRITE (*, '(a)') ' QTF field: more frequencies than the store was set up for.'
         ERROR STOP 1
      END IF
      oufr_store(nkept) = oufr
      DO ib = 1, nbeta
         w_store(nkept, ib) = states(ib)%omega
         k_store(nkept, ib) = states(ib)%k
         h_store(nkept, ib) = states(ib)%h
         b_store(nkept, ib) = states(ib)%beta
         rho_store(nkept, ib) = states(ib)%rho
         phi_store(:, nkept, ib) = fields(ib)%phi
         grad_store(:, :, nkept, ib) = fields(ib)%grad
      END DO
   END SUBROUTINE qfield_keep

   SUBROUTINE qfield_finish(refl)
      ! Writes Wamit_format/QTF_FieldPartial.txt and frees the store. The rows are ordered by kernel (K = 1 difference,
      ! then K = 2 sum), heading pair, frequency pair, then point. No-op when inactive.
      REAL(8), INTENT(IN) :: refl
      TYPE(incident2) :: inc
      COMPLEX(8) :: quad, pot
      REAL(8) :: g_wave
      INTEGER :: unit, kern, kn, b1, b2, i, j, ip, nrows

      IF (.NOT. active) RETURN
      OPEN (NEWUNIT=unit, FILE=outdir//'/Wamit_format/QTF_FieldPartial.txt', STATUS='REPLACE', ACTION='WRITE')
      CALL write_header(unit)
      nrows = 0
      DO kern = 1, 2
         IF (kern == 1 .AND. .NOT. (kind_qtf == 1 .OR. kind_qtf == 3)) CYCLE
         IF (kern == 2 .AND. .NOT. (kind_qtf == 2 .OR. kind_qtf == 3)) CYCLE
         kn = MERGE(KERNEL_DIFF, KERNEL_SUM, kern == 1)
         DO b1 = 1, nbeta
            DO b2 = 1, nbeta
               IF (all_headings_pairs == 0 .AND. b1 /= b2) CYCLE
               DO i = 1, nkept
                  DO j = 1, nkept
                     IF (kn == KERNEL_DIFF .AND. max_dw > 0D0) THEN
                        IF (ABS(w_store(i, b1) - w_store(j, b2)) > max_dw*(1D0 + 1D-9)) CYCLE
                     END IF
                     ! The gravity of the first-order incident wave, as for QTF_Indirect.11d (wave_gravity); Bernoulli
                     ! takes the g of Const_mod, the one of the first-order elevation
                     g_wave = wave_gravity(w_store(i, b1), k_store(i, b1), h_store(i, b1))
                     CALL incident2_pair(w_store(i, b1), k_store(i, b1), b_store(i, b1), w_store(j, b2), k_store(j, b2), &
                                         b_store(j, b2), h_store(i, b1), g_wave, inc, kn)
                     DO ip = 1, npt
                        CALL field_pair(kn, w_store(i, b1), phi_store(ip, i, b1), grad_store(:, ip, i, b1), &
                                        w_store(j, b2), phi_store(ip, j, b2), grad_store(:, ip, j, b2), inc, px(:, ip), &
                                        is_surface(ip), rho_store(i, b1), G, quad, pot)
                        CALL write_field_row(unit, oufr_store(i), oufr_store(j), beta_deg(b1), beta_deg(b2), kern, ip, &
                                             is_surface(ip), quad, pot, refl)
                        nrows = nrows + 1
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      CLOSE (unit)
      WRITE (*, '(a,I0,a,I0,a,I0,a)') ' QTF field points: ', nrows, ' rows at ', npt, ' points, ', nkept, ' frequencies'
      DEALLOCATE (phi_store, grad_store, w_store, k_store, h_store, b_store, rho_store, oufr_store, px, is_surface, beta_deg)
      active = .FALSE.
   END SUBROUTINE qfield_finish

   SUBROUTINE write_header(unit)
      INTEGER, INTENT(IN) :: unit
      INTEGER :: ip

      WRITE (unit, '(a)') '# QTF_FieldPartial.txt: PARTIAL second-order pressure and elevation at the field points of '// &
         'ControlFile.in'
      WRITE (unit, '(a)') '# OMITTED: the second-order diffracted and radiated waves. Near the bodies the omitted part is'// &
         ' not small;'
      WRITE (unit, '(a)') '#          use these values for the far field and for checks only.'
      WRITE (unit, '(a)') '# Included: the quadratic part (products of first-order quantities) and the second-order'// &
         ' incident wave.'
      WRITE (unit, '(a)') '# Columns: PER1 PER2 BETA1 BETA2 K IPT TYPE Mod Pha Re Im Re_quad Im_quad Re_inc Im_inc'
      WRITE (unit, '(a)') '# K: 1 difference frequency (T-), 2 sum frequency (T+). TYPE: 1 pressure, 2 elevation.'
      WRITE (unit, '(a)') '# e^{+i omega t} as the other WAMIT-format files; pressure / (rho g A^2 / L), elevation / (A^2 / L),'// &
         ' L = Reference_body_length.'
      WRITE (unit, '(a)') '# POINT IPT X Y Z TYPE'
      DO ip = 1, npt
         WRITE (unit, '(a,I6,3ES14.6,I6)') '# POINT', ip, px(:, ip), MERGE(2, 1, is_surface(ip))
      END DO
   END SUBROUTINE write_header

END MODULE QTFField
