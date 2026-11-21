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
   ! physical e^{-i Omega t} amplitude is, for i /= j, 2 a_i a_j* (quad + pot) for the difference and 2 a_i a_j (quad + pot)
   ! for the sum; for i = j it is a_i^2 (quad + pot) for the sum and |a_i|^2 (quad + pot) for the mean (difference).
   !
   ! The store and the file (qtf_field = 1 in SecondOrder.in). At the end of each frequency of the QTF grid the driver gives
   ! qfield_keep the first-order potential and gradient of every heading at the field points (fields_green of
   ! FieldVelocity); these are all that the pairs need, so the store is small: 4 complex values per point, frequency and
   ! heading. After the frequency loop qfield_finish forms every pair of the QTFs of QTFDriver, with the same kernels
   ! (qtf = 1 difference, 2 sum, 3 both), heading pairs (mean_drift = 2: all, else beta_1 = beta_2) and qtf_max_dw on the
   ! difference pairs, so that its rows are those of QTF_Quadratic.10d and .10s at each point.
   !
   ! The direct part (qtf_field = 1 with qtf_direct = 1, Phase 7c-1). The direct solve of QTFDriver (qtf_direct_at_omega)
   ! gives phi_D of every difference-frequency pair at the field points (direct_at_omega, field_points), and the driver hands
   ! them to qfield_keep_direct. The second-order potential enters the pressure and the elevation as the second-order incident
   ! potential does (pot_p, pot_eta above), so with e^{-i Omega t} and Omega = omega_i - omega_j signed:
   !   p_D = -rho d(phi_D)/dt = i rho Omega phi_D,   eta_D = -(1/g) d(phi_D)/dt = (i Omega/g) phi_D   (field_direct)
   ! per unit amplitudes, like quad and pot; the amplitude 2 a_i a_j* (i /= j) of the file header multiplies it as well.
   ! direct_at_omega returns phi_D of a pair with Omega < 0 conjugated, i.e. as the pair (j, i) with Omega > 0 would give it
   ! conjugated, and the signed Omega here is the one that goes with it: p_D of the pair (j, i) is the conjugate of that of
   ! (i, j). Rows the direct solve does not cover (mean rows, i = j, and the pairs it did not keep) have no direct part: 0 and
   ! has_direct = 0.
   !
   ! The sum frequency (Phase 7c-2, qtf = 2 or 3 with qtf_direct = 1). The sum columns of direct_at_omega (kernel = KERNEL_SUM,
   ! Omega = omega_i + omega_j > 0, never conjugated: the ordered pair w1 >= w2 of QTFDriver) go to qfield_keep_direct with the
   ! kernel. The sum pair is symmetric, phi_D(i, j) = phi_D(j, i), so one call fills both slots. The pressure and the
   ! elevation are the same forms with the sum Omega, i rho Omega phi_D and (i Omega/g) phi_D, per unit amplitudes; the
   ! physical amplitude is 2 a_i a_j (i /= j) or a_i^2 (i = j), the 5b factors above, applied by the reader as for quad and
   ! pot. A sum row the direct solve did not keep (outside qtf_max_sw, beyond the irregular-frequency cutoff of QTFDriver)
   ! has 0 and has_direct = 0. The direct part is the body-forced phi_D for both kernels: with qtf_terms = full the
   ! free-surface forcing is left out of the field, as for the difference.
   ! The value with the direct part is the sum of the Re Im of QTF_FieldPartial.txt and Re_dir Im_dir of QTF_Field.txt; the
   ! first file is written as without qtf_direct, to the byte.
   !
   ! The radiated wave of the second-order motion (Phase 6c, qtf_motion = 1 with qtf_direct = 1).
   ! A free body moves at Omega with xi (QTF_Motion.12d, .12s: the motion of the direct route, e^{-i Omega t}, unscaled),
   ! and radiates phi_R. With the velocity -i Omega xi_j of mode j and psi_j the unit-velocity radiation potential at Omega
   ! (header of QTFMotion: phi = psi_k (-i Omega xi_k)),
   !   phi_R = sum_j (-i Omega xi_j) psi_j                                                         (radiated_potential)
   ! psi_j at the points comes from Green's formula with the body values of psi_j and dpsi_j/dn = n_j (psi_at_points of
   ! QTFDriver: direct_field_at_points on the 6 nbody columns). phi_R is a second-order potential at Omega like phi_D, so
   ! it enters as phi_D does, p_R = i rho Omega phi_R and eta_R = (i Omega/g) phi_R (field_direct), and the driver adds
   ! it to the stored phi_D of the row (qfield_add_radiated): Re_dir Im_dir then hold phi_D + phi_R.
   ! The sign of Omega. psi is solved at |Omega| > 0. A difference row with Omega = omega_i - omega_j < 0 holds phi_D
   ! conjugated (above): the e^{-i Omega t} amplitude at the signed Omega, conj of the pair (j, i). Its xi is stored the
   ! same way (qtf_motion_part: conj(lhs) xi = F at Omega < 0). The radiation potential of unit velocity amplitude at
   ! -|Omega| is conj(psi(|Omega|)): the velocity Re(V e^{+i|Omega|t}) = Re(conj(V) e^{-i|Omega|t}) radiates
   ! Re(conj(V) psi e^{-i|Omega|t}) = Re(V conj(psi) e^{+i|Omega|t}). So phi_R of the row (i, j) at Omega < 0 is
   ! sum_j (-i Omega xi_j) conj(psi_j), which is conj(phi_R of (j, i)): the row is conj of its mirror, as phi_D, and
   ! field_direct with the signed Omega gives p_R(i, j) = conj(p_R(j, i)), as for the direct rows. The file conjugates
   ! once more (e^{+i omega t}, write_field_row), as for every column.
   ! Rows left out of the motion (an Omega whose matrix has rcond < 1e-12, as QTF_Motion leaves them out) get no
   ! phi_R. QTF_Field.txt then has the column has_rad after has_direct, as has_direct marks the rows without a direct part:
   ! 1 where Re_dir Im_dir include phi_R, 0 where they do not (no direct part, or the motion is left out). With
   ! qtf_motion = 0 neither the column nor the header lines change: the file is that of Phase 7c-2, to the byte.

   USE Const_mod, ONLY: G
   USE SecondOrderIncident, ONLY: incident2, incident2_pair, incident2_phi
   USE NearFieldDrift, ONLY: pair_product, KERNEL_DIFF, KERNEL_SUM
   USE SecondOrderTypes, ONLY: so_options, first_order_state, panel_mesh, cs_points, cs_fields, CS_FREE
   USE SecondOrderOutput, ONLY: write_field_row
   USE FirstOrderState, ONLY: wave_gravity
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: field_pair, build_field_points, body_winding
   PUBLIC :: qfield_init, qfield_keep, qfield_finish, qfield_stored_mb, FIELD_FD_STEP
   PUBLIC :: field_direct, qfield_keep_direct, qfield_direct_active, qfield_direct_points
   PUBLIC :: radiated_potential, qfield_add_radiated, qfield_radiated_active

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
   ! The direct part: phi_D of the pair (i, j) at the points, heading pair as in QTFDriver (b1 = b2 unless mean_drift = 2)
   LOGICAL :: direct_on = .FALSE., direct_diff = .FALSE., direct_sum = .FALSE.   ! any kernel, difference, sum
   COMPLEX(8), ALLOCATABLE :: dir_store(:, :, :, :)     ! (npt, nfreq, nfreq, nheading pairs), phi_D, 0 if not kept
   LOGICAL, ALLOCATABLE :: dir_have(:, :, :)            ! (nfreq, nfreq, nheading pairs), phi_D of the pair was kept
   COMPLEX(8), ALLOCATABLE :: dir_store_sum(:, :, :, :) ! the same for the sum kernel, both slots of a pair filled
   LOGICAL, ALLOCATABLE :: dir_have_sum(:, :, :)
   ! The radiated wave of the second-order motion (qtf_motion = 1, Phase 6c): phi_R is added to dir_store (dir_store_sum)
   ! of the row, and rad_have (rad_have_sum) marks the rows that have it, as dir_have
   LOGICAL :: radiated_on = .FALSE.
   LOGICAL, ALLOCATABLE :: rad_have(:, :, :), rad_have_sum(:, :, :)   ! as dir_have, dir_have_sum

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

   PURE COMPLEX(8) FUNCTION field_direct(omega, phi_d, rho, g, surface)
      ! The direct part of the pressure (surface = .FALSE.) or the elevation (surface = .TRUE.) of a pair with the signed
      ! difference frequency omega = omega_1 - omega_2, from phi_D of the pair at the point (header: the same forms as the
      ! second-order incident part, with phi_D for phi2). The time factor is e^{-i omega t}, so d(phi_D)/dt = -i omega phi_D
      ! and p = -rho d(phi_D)/dt = i rho omega phi_D, eta = -(1/g) d(phi_D)/dt = (i omega/g) phi_D. Per unit amplitudes.
      REAL(8), INTENT(IN) :: omega, rho, g
      COMPLEX(8), INTENT(IN) :: phi_d
      LOGICAL, INTENT(IN) :: surface
      COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)

      IF (surface) THEN
         field_direct = CI*omega/g*phi_d
      ELSE
         field_direct = CI*rho*omega*phi_d
      END IF
   END FUNCTION field_direct

   PURE FUNCTION radiated_potential(omega, xi, psi_pts) RESULT(phi_r)
      ! The radiated wave phi_R = sum_j (-i omega xi_j) psi_j at the points (header), e^{-i omega t}, of a row with the
      ! signed frequency omega (difference omega_1 - omega_2, or the sum) and the motion xi (6 nbody, as xi2 of QTFDriver
      ! stores it). psi_pts (npt, 6 nbody) is psi at the points at |omega| (psi_at_points); for omega < 0 the radiation
      ! potential at omega is its conjugate (header), so conj(psi_pts) is taken. omega = 0 gives 0.
      REAL(8), INTENT(IN) :: omega
      COMPLEX(8), INTENT(IN) :: xi(:), psi_pts(:, :)
      COMPLEX(8) :: phi_r(SIZE(psi_pts, 1))
      COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
      INTEGER :: k

      phi_r = (0D0, 0D0)
      DO k = 1, SIZE(xi)
         IF (omega < 0D0) THEN
            phi_r = phi_r + (-CI*omega*xi(k))*CONJG(psi_pts(:, k))
         ELSE
            phi_r = phi_r + (-CI*omega*xi(k))*psi_pts(:, k)
         END IF
      END DO
   END FUNCTION radiated_potential

   SUBROUTINE build_field_points(mesh, xfp, pts, surface, ok, msg)
      ! Field points xfp(3,:) as a cs_points set (kind CS_FREE, nrm = 0, w = 0). Type as PressureElevation: |z| <= 1e-6 is
      ! a surface (elevation) point and gets z = 0, z < -1e-6 is a pressure point. psize and body are those of the nearest
      ! hull panel (centroids of the full panel set, so symmetry images are included). Rejected, checked in point order,
      ! the first failing point giving msg: above the free surface; within one local panel size of a hull panel centroid
      ! (a point has no panel size of its own, so the limit is the hull panel's size); inside a body, decided by the
      ! winding number of body_winding. The nearest centroid alone decides wrongly near edges where panel sizes change.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: xfp(:,:)
      TYPE(cs_points), INTENT(OUT) :: pts
      LOGICAL, ALLOCATABLE, INTENT(OUT) :: surface(:)
      LOGICAL, INTENT(OUT) :: ok
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: msg
      REAL(8), PARAMETER :: ZTOL = 1D-6
      REAL(8) :: x(3), d2, d2min
      REAL(8), ALLOCATABLE :: wind(:)
      INTEGER :: ip, i, imin, ihull, n, ib
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
         CALL body_winding(mesh, x, wind)
         ib = MAXLOC(wind, 1)
         IF (wind(ib) > 0.5D0) THEN
            WRITE (line, '(a,i0,a,3(1x,es11.3),a,i0)') 'field point ', ip, ' at', x, ' is inside body ', ib
            msg = TRIM(line)
            ok = .FALSE.
            RETURN
         END IF
      END DO
   END SUBROUTINE build_field_points

   SUBROUTINE body_winding(mesh, x, wind)
      ! Winding number wind(b) of each body b about x: the solid angle of its hull panels plus their mirror in z = 0, over
      ! 4 pi. A surface-piercing hull and its mirror form a closed surface; a submerged hull is closed and its mirror lies
      ! above z = 0. So wind is 1 inside a body (or its waterplane, for x on z = 0) and 0 outside, for a watertight hull.
      ! Panels are split into triangles (1,2,3) and (1,3,4), each turned to face out of the body (against nrm, which points
      ! into it), and the solid angle of a triangle is that of Van Oosterom and Strackee (1983).
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: x(3)
      REAL(8), ALLOCATABLE, INTENT(OUT) :: wind(:)
      REAL(8), PARAMETER :: PI = 3.14159265358979323846D0
      REAL(8) :: v(3, 4), a(3), b(3), c(3), f(3), nin(3), la, lb, lc, num, den
      INTEGER :: i, is, it, m

      ALLOCATE (wind(mesh%nbody))
      wind = 0D0
      DO i = 1, mesh%n
         DO is = 1, 2
            f = [1D0, 1D0, MERGE(1D0, -1D0, is == 1)]
            nin = f*mesh%nrm(:, i)
            DO m = 1, mesh%ncn(i)
               v(:, m) = f*mesh%xnode(:, mesh%ncon(m, i)) - x
            END DO
            DO it = 1, mesh%ncn(i) - 2
               a = v(:, 1)
               b = v(:, it + 1)
               c = v(:, it + 2)
               IF (SUM(cross(b - a, c - a)*nin) > 0D0) THEN
                  b = v(:, it + 2)
                  c = v(:, it + 1)
               END IF
               la = NORM2(a)
               lb = NORM2(b)
               lc = NORM2(c)
               num = SUM(a*cross(b, c))
               den = la*lb*lc + SUM(a*b)*lc + SUM(a*c)*lb + SUM(b*c)*la
               wind(mesh%body(i)) = wind(mesh%body(i)) + 2D0*ATAN2(num, den)
            END DO
         END DO
      END DO
      wind = wind/(4D0*PI)
   CONTAINS
      PURE FUNCTION cross(p, q)
         REAL(8), INTENT(IN) :: p(3), q(3)
         REAL(8) :: cross(3)
         cross = [p(2)*q(3) - p(3)*q(2), p(3)*q(1) - p(1)*q(3), p(1)*q(2) - p(2)*q(1)]
      END FUNCTION cross
   END SUBROUTINE body_winding

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
      ! The direct solve gives the difference kernel with qtf = 1 or 3 and the sum kernel with qtf = 2 or 3 (Phase 7c-2)
      direct_diff = opts%qtf_direct == 1 .AND. (opts%qtf == 1 .OR. opts%qtf == 3)
      direct_sum = opts%qtf_direct == 1 .AND. (opts%qtf == 2 .OR. opts%qtf == 3)
      direct_on = direct_diff .OR. direct_sum
      ! The radiated wave needs the direct part and a free body: a fixed one has no motion (QTFDriver computes none)
      radiated_on = direct_on .AND. opts%qtf_motion == 1 .AND. .NOT. opts%motion_fixed
      nbeta = SIZE(beta_in_deg)
      npt = pts%n
      nkept = 0
      outdir = outputdir
      px = pts%x
      is_surface = surface
      beta_deg = beta_in_deg
      IF (ALLOCATED(phi_store)) DEALLOCATE (phi_store, grad_store, w_store, k_store, h_store, b_store, rho_store, oufr_store)
      IF (ALLOCATED(dir_store)) DEALLOCATE (dir_store, dir_have)
      IF (ALLOCATED(dir_store_sum)) DEALLOCATE (dir_store_sum, dir_have_sum)
      IF (ALLOCATED(rad_have)) DEALLOCATE (rad_have)
      IF (ALLOCATED(rad_have_sum)) DEALLOCATE (rad_have_sum)
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
      IF (direct_diff) THEN
         ALLOCATE (dir_store(npt, nfreq, nfreq, MERGE(nbeta*nbeta, nbeta, all_headings_pairs == 1)), STAT=st)
         IF (st == 0) ALLOCATE (dir_have(nfreq, nfreq, SIZE(dir_store, 4)), STAT=st)
         IF (st /= 0) CALL direct_memory_stop()
         dir_store = (0D0, 0D0)
         dir_have = .FALSE.
      END IF
      IF (direct_sum) THEN
         ALLOCATE (dir_store_sum(npt, nfreq, nfreq, MERGE(nbeta*nbeta, nbeta, all_headings_pairs == 1)), STAT=st)
         IF (st == 0) ALLOCATE (dir_have_sum(nfreq, nfreq, SIZE(dir_store_sum, 4)), STAT=st)
         IF (st /= 0) CALL direct_memory_stop()
         dir_store_sum = (0D0, 0D0)
         dir_have_sum = .FALSE.
      END IF
      IF (radiated_on .AND. direct_diff) THEN
         ALLOCATE (rad_have(nfreq, nfreq, SIZE(dir_store, 4)))
         rad_have = .FALSE.
      END IF
      IF (radiated_on .AND. direct_sum) THEN
         ALLOCATE (rad_have_sum(nfreq, nfreq, SIZE(dir_store_sum, 4)))
         rad_have_sum = .FALSE.
      END IF
   CONTAINS
      SUBROUTINE direct_memory_stop()
         WRITE (*, '(a)') ' QTF field: not enough memory for phi_D of the direct solve at the field points. Use fewer'
         WRITE (*, '(a)') ' frequencies, headings or field points, or qtf_direct = 0.'
         ERROR STOP 1
      END SUBROUTINE direct_memory_stop
   END SUBROUTINE qfield_init

   REAL(8) FUNCTION qfield_stored_mb(npoints, nfreq, nb, direct_pairs)
      ! Size of the store in MB: 4 complex values per point, frequency and heading, and with direct_pairs (the number of
      ! heading pairs of the direct part, when it is kept) one complex value per point, frequency pair and heading pair
      INTEGER, INTENT(IN) :: npoints, nfreq, nb
      INTEGER, INTENT(IN), OPTIONAL :: direct_pairs

      qfield_stored_mb = 64D0*REAL(npoints, 8)*REAL(nfreq, 8)*REAL(nb, 8)/1D6
      IF (PRESENT(direct_pairs)) qfield_stored_mb = qfield_stored_mb + &
         16D0*REAL(npoints, 8)*REAL(nfreq, 8)**2*REAL(direct_pairs, 8)/1D6
   END FUNCTION qfield_stored_mb

   LOGICAL FUNCTION qfield_direct_active()
      ! The store keeps the direct part: qtf_field = 1 with qtf_direct = 1 (and the difference QTF), set up by qfield_init
      qfield_direct_active = active .AND. direct_on
   END FUNCTION qfield_direct_active

   SUBROUTINE qfield_direct_points(x)
      ! The field points (3, npt) of the store, for direct_at_omega. Stops if the direct part is not kept.
      REAL(8), ALLOCATABLE, INTENT(OUT) :: x(:, :)

      IF (.NOT. qfield_direct_active()) THEN
         WRITE (*, '(a)') ' QTF field: qfield_direct_points called without qtf_field = 1 and qtf_direct = 1.'
         ERROR STOP 1
      END IF
      x = px
   END SUBROUTINE qfield_direct_points

   SUBROUTINE qfield_keep_direct(i, b1, j, b2, phi_d, kernel)
      ! Keeps phi_D at the points (npt) of the pair of frequencies i, j (in the order of qfield_keep) and headings b1, b2:
      ! as direct_at_omega returns it, e^{-i omega t}, conjugated for omega_i < omega_j (header). No-op unless the direct part
      ! is kept. The heading pair must be one of the QTFs (b1 = b2 unless mean_drift = 2). kernel (default KERNEL_DIFF, as
      ! before) = KERNEL_SUM keeps the sum column, not conjugated, in both slots (i, b1, j, b2) and (j, b2, i, b1) (header).
      INTEGER, INTENT(IN) :: i, b1, j, b2
      COMPLEX(8), INTENT(IN) :: phi_d(:)
      INTEGER, INTENT(IN), OPTIONAL :: kernel
      INTEGER :: h, hm
      LOGICAL :: is_sum

      IF (.NOT. qfield_direct_active()) RETURN
      is_sum = .FALSE.
      IF (PRESENT(kernel)) is_sum = kernel == KERNEL_SUM
      IF (is_sum .AND. .NOT. direct_sum) THEN
         WRITE (*, '(a)') ' QTF field: qfield_keep_direct called with the sum kernel for a run without the direct sum.'
         ERROR STOP 1
      END IF
      IF (.NOT. is_sum .AND. .NOT. direct_diff) THEN
         WRITE (*, '(a)') ' QTF field: qfield_keep_direct called with the difference kernel for a run without it.'
         ERROR STOP 1
      END IF
      IF (i < 1 .OR. j < 1 .OR. i > SIZE(phi_store, 2) .OR. j > SIZE(phi_store, 2) .OR. b1 < 1 .OR. b2 < 1 .OR. &
          b1 > nbeta .OR. b2 > nbeta .OR. SIZE(phi_d) /= npt) THEN
         WRITE (*, '(a)') ' QTF field: qfield_keep_direct called with a pair or a size that is not in the store.'
         ERROR STOP 1
      END IF
      IF (all_headings_pairs == 0 .AND. b1 /= b2) THEN
         WRITE (*, '(a)') ' QTF field: qfield_keep_direct called for a heading pair that mean_drift does not take.'
         ERROR STOP 1
      END IF
      h = heading_pair_index(b1, b2)
      IF (is_sum) THEN
         hm = heading_pair_index(b2, b1)
         dir_store_sum(:, i, j, h) = phi_d
         dir_have_sum(i, j, h) = .TRUE.
         dir_store_sum(:, j, i, hm) = phi_d
         dir_have_sum(j, i, hm) = .TRUE.
         RETURN
      END IF
      dir_store(:, i, j, h) = phi_d
      dir_have(i, j, h) = .TRUE.
   END SUBROUTINE qfield_keep_direct

   LOGICAL FUNCTION qfield_radiated_active()
      ! The store takes the radiated wave of the second-order motion: the direct part with qtf_motion = 1 and a free body
      ! (header), set up by qfield_init
      qfield_radiated_active = qfield_direct_active() .AND. radiated_on
   END FUNCTION qfield_radiated_active

   SUBROUTINE qfield_add_radiated(i1, b1, i2, b2, phi_r, kernel)
      ! Adds phi_R at the points (npt, radiated_potential) to the direct part of the row (frequency i1, heading b1; frequency
      ! i2, heading b2) of kernel (KERNEL_DIFF or KERNEL_SUM), and marks the row (has_rad = 1). phi_R is in the form of the
      ! stored phi_D of that row: e^{-i Omega t} at the signed Omega (header). One slot only: for the sum, whose two slots
      ! keep_direct fills with one call, the driver calls this for each slot, with the xi of each. No-op unless the store
      ! takes the radiated wave. Stops for a row without a direct part (the radiated wave goes only into rows that have
      ! one) or a row given twice.
      INTEGER, INTENT(IN) :: i1, b1, i2, b2, kernel
      COMPLEX(8), INTENT(IN) :: phi_r(:)
      INTEGER :: h

      IF (.NOT. qfield_radiated_active()) RETURN
      IF (i1 < 1 .OR. i2 < 1 .OR. i1 > SIZE(phi_store, 2) .OR. i2 > SIZE(phi_store, 2) .OR. b1 < 1 .OR. b2 < 1 .OR. &
          b1 > nbeta .OR. b2 > nbeta .OR. SIZE(phi_r) /= npt .OR. (all_headings_pairs == 0 .AND. b1 /= b2)) THEN
         WRITE (*, '(a)') ' QTF field: qfield_add_radiated called with a pair or a size that is not in the store.'
         ERROR STOP 1
      END IF
      h = heading_pair_index(b1, b2)
      IF (kernel == KERNEL_SUM .AND. direct_sum) THEN
         IF (.NOT. dir_have_sum(i1, i2, h) .OR. rad_have_sum(i1, i2, h)) CALL row_stop()
         dir_store_sum(:, i1, i2, h) = dir_store_sum(:, i1, i2, h) + phi_r
         rad_have_sum(i1, i2, h) = .TRUE.
      ELSE IF (kernel == KERNEL_DIFF .AND. direct_diff) THEN
         IF (.NOT. dir_have(i1, i2, h) .OR. rad_have(i1, i2, h)) CALL row_stop()
         dir_store(:, i1, i2, h) = dir_store(:, i1, i2, h) + phi_r
         rad_have(i1, i2, h) = .TRUE.
      ELSE
         WRITE (*, '(a)') ' QTF field: qfield_add_radiated called with a kernel the direct part does not have.'
         ERROR STOP 1
      END IF
   CONTAINS
      SUBROUTINE row_stop()
         WRITE (*, '(a,4(1x,I0))') ' QTF field: qfield_add_radiated: the row has no direct part or has phi_R already:', &
            i1, b1, i2, b2
         ERROR STOP 1
      END SUBROUTINE row_stop
   END SUBROUTINE qfield_add_radiated

   PURE INTEGER FUNCTION heading_pair_index(b1, b2)
      ! The heading pair of the direct store: all pairs with mean_drift = 2, else only b1 = b2
      INTEGER, INTENT(IN) :: b1, b2

      heading_pair_index = MERGE((b1 - 1)*nbeta + b2, b1, all_headings_pairs == 1)
   END FUNCTION heading_pair_index

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
      COMPLEX(8) :: quad, pot, d
      REAL(8) :: g_wave
      INTEGER :: unit, unit_full, kern, kn, b1, b2, i, j, ip, nrows, h
      LOGICAL :: have, rad

      IF (.NOT. active) RETURN
      OPEN (NEWUNIT=unit, FILE=outdir//'/Wamit_format/QTF_FieldPartial.txt', STATUS='REPLACE', ACTION='WRITE')
      CALL write_header(unit, .FALSE.)
      ! With the direct part the same rows go to QTF_Field.txt as well, with the columns of the direct part after theirs
      unit_full = 0
      IF (direct_on) THEN
         OPEN (NEWUNIT=unit_full, FILE=outdir//'/Wamit_format/QTF_Field.txt', STATUS='REPLACE', ACTION='WRITE')
         CALL write_header(unit_full, .TRUE.)
      END IF
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
                        IF (direct_on) THEN
                           ! Only the difference rows of pairs the direct solve gave; the mean rows (i = j), the sum rows
                           ! and the pairs outside qtf_max_dw have none (header)
                           h = heading_pair_index(b1, b2)
                           have = .FALSE.
                           IF (kn == KERNEL_DIFF .AND. direct_diff) have = i /= j .AND. dir_have(i, j, h)
                           IF (kn == KERNEL_SUM .AND. direct_sum) have = dir_have_sum(i, j, h)
                           d = (0D0, 0D0)
                           IF (have .AND. kn == KERNEL_DIFF) &
                              d = field_direct(w_store(i, b1) - w_store(j, b2), dir_store(ip, i, j, h), rho_store(i, b1), G, &
                                               is_surface(ip))
                           IF (have .AND. kn == KERNEL_SUM) &
                              d = field_direct(w_store(i, b1) + w_store(j, b2), dir_store_sum(ip, i, j, h), rho_store(i, b1), &
                                               G, is_surface(ip))
                           IF (radiated_on) THEN
                              ! has_rad: the row has phi_R in d (qfield_add_radiated); 0 where it has no direct part
                              rad = .FALSE.
                              IF (have .AND. kn == KERNEL_DIFF) rad = rad_have(i, j, h)
                              IF (have .AND. kn == KERNEL_SUM) rad = rad_have_sum(i, j, h)
                              CALL write_field_row(unit_full, oufr_store(i), oufr_store(j), beta_deg(b1), beta_deg(b2), &
                                                   kern, ip, is_surface(ip), quad, pot, refl, d, have, rad)
                           ELSE
                              CALL write_field_row(unit_full, oufr_store(i), oufr_store(j), beta_deg(b1), beta_deg(b2), &
                                                   kern, ip, is_surface(ip), quad, pot, refl, d, have)
                           END IF
                        END IF
                        nrows = nrows + 1
                     END DO
                  END DO
               END DO
            END DO
         END DO
      END DO
      CLOSE (unit)
      IF (direct_on) CLOSE (unit_full)
      WRITE (*, '(a,I0,a,I0,a,I0,a)') ' QTF field points: ', nrows, ' rows at ', npt, ' points, ', nkept, ' frequencies'
      DEALLOCATE (phi_store, grad_store, w_store, k_store, h_store, b_store, rho_store, oufr_store, px, is_surface, beta_deg)
      IF (direct_diff) DEALLOCATE (dir_store, dir_have)
      IF (direct_sum) DEALLOCATE (dir_store_sum, dir_have_sum)
      IF (ALLOCATED(rad_have)) DEALLOCATE (rad_have)
      IF (ALLOCATED(rad_have_sum)) DEALLOCATE (rad_have_sum)
      active = .FALSE.
      radiated_on = .FALSE.
      direct_on = .FALSE.
      direct_diff = .FALSE.
      direct_sum = .FALSE.
   END SUBROUTINE qfield_finish

   SUBROUTINE write_header(unit, complete)
      ! The header of QTF_FieldPartial.txt (complete = .FALSE.), unchanged by the direct part, or of QTF_Field.txt
      ! (.TRUE.): the columns of the partial file and of the direct part, and what the direct part holds and leaves out
      INTEGER, INTENT(IN) :: unit
      LOGICAL, INTENT(IN) :: complete
      INTEGER :: ip
      CHARACTER(LEN=:), ALLOCATABLE :: rad_text

      ! With qtf_motion = 0 the words of Phase 7c-2, so that the header stays as it was, to the byte
      IF (radiated_on) THEN
         rad_text = 'radiated wave of the second-order motion included where has_rad = 1 (Phase 6c)'
      ELSE
         rad_text = 'radiated wave of the second-order motion not included (Phase 6c)'
      END IF
      IF (complete) THEN
         WRITE (unit, '(a)') '# QTF_Field.txt: second-order pressure and elevation at the field points of ControlFile.in: '// &
            'the rows of QTF_FieldPartial.txt and the direct part'
         WRITE (unit, '(a)') '# Re Im, Re_quad Im_quad and Re_inc Im_inc: as in QTF_FieldPartial.txt, the partial value (the'// &
            ' quadratic part and the second-order incident wave).'
         IF (.NOT. direct_sum) THEN
            ! No sum column: the header lines of Phase 7a, so that QTF_Field.txt of qtf = 1 stays as it was, to the byte
            WRITE (unit, '(a)') '# Re_dir Im_dir: the direct part, i rho Omega phi_D (pressure) and (i Omega/g) phi_D (elevation),'// &
               ' Omega = omega_1 - omega_2, per unit amplitude pair, conjugated and scaled as Re Im. Add it to Re Im.'
            WRITE (unit, '(a)') '# phi_D of the direct solve with body forcing; free-surface forcing not included (it is in .12d with qtf_terms = full);'// &
               ' '//rad_text
            WRITE (unit, '(a)') '# has_direct: 1 where the direct solve covers the row (difference frequency, omega_1 /='// &
               ' omega_2), else 0 with Re_dir = Im_dir = 0 (mean rows and sum rows).'
         ELSE
            WRITE (unit, '(a)') '# Re_dir Im_dir: the direct part, i rho Omega phi_D (pressure) and (i Omega/g) phi_D (elevation),'// &
               ' Omega = omega_1 - omega_2 (K = 1) or omega_1 + omega_2 (K = 2), per unit amplitude pair, conjugated and scaled'// &
               ' as Re Im. Add it to Re Im.'
            WRITE (unit, '(a)') '# phi_D of the direct solve with body forcing, for both kernels; free-surface forcing not'// &
               ' included (it is in .12d and .12s with qtf_terms = full); '//rad_text
            WRITE (unit, '(a)') '# has_direct: 1 where the direct solve covers the row (K = 1: omega_1 /= omega_2; K = 2: the'// &
               ' sum pairs of QTF_Direct.12s), else 0 with Re_dir = Im_dir = 0 (mean rows and rows outside the direct solve).'
         END IF
         IF (radiated_on) THEN
            ! The radiated wave: what Re_dir Im_dir add, and the column that marks the rows left out of the motion
            WRITE (unit, '(a)') '# has_rad: 1 where Re_dir Im_dir also hold the radiated wave of the second-order motion,'// &
               ' i rho Omega phi_R (pressure) and (i Omega/g) phi_R (elevation), phi_R = sum_j (-i Omega xi_j) psi_j, xi of'// &
               ' QTF_Motion.12d (K = 1) or .12s (K = 2); else 0: no direct part, or the Omega of the row is left out of the'// &
               ' motion (rcond < 1e-12, see the header of QTF_Motion) and Re_dir Im_dir hold phi_D only.'
            WRITE (unit, '(a)') '# Columns: PER1 PER2 BETA1 BETA2 K IPT TYPE Mod Pha Re Im Re_quad Im_quad Re_inc Im_inc'// &
               ' Re_dir Im_dir has_direct has_rad'
         ELSE
            WRITE (unit, '(a)') '# Columns: PER1 PER2 BETA1 BETA2 K IPT TYPE Mod Pha Re Im Re_quad Im_quad Re_inc Im_inc'// &
               ' Re_dir Im_dir has_direct'
         END IF
         WRITE (unit, '(a)') '# K: 1 difference frequency (T-), 2 sum frequency (T+). TYPE: 1 pressure, 2 elevation.'
         WRITE (unit, '(a)') '# e^{+i omega t} as the other WAMIT-format files; pressure / (rho g A^2 / L), elevation / (A^2 /'// &
            ' L), L = Reference_body_length.'
         WRITE (unit, '(a)') '# POINT IPT X Y Z TYPE'
         DO ip = 1, npt
            WRITE (unit, '(a,I6,3ES14.6,I6)') '# POINT', ip, px(:, ip), MERGE(2, 1, is_surface(ip))
         END DO
         RETURN
      END IF
      WRITE (unit, '(a)') '# QTF_FieldPartial.txt: PARTIAL second-order pressure and elevation at the field points of '// &
         'ControlFile.in'
      WRITE (unit, '(a)') '# OMITTED: the second-order diffracted and radiated waves. Near the bodies the omitted part is'// &
         ' not small;'
      WRITE (unit, '(a)') '#          they are small only where those waves are small against the incident wave.'
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
