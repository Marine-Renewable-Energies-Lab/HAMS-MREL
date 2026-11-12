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
MODULE QTFFreeSurface
   ! Free-surface integral of the second-order potential in the difference-frequency QTF (Phase 6a). The forcing of the
   ! second-order potential phi_D (diffracted and radiated part) on z = 0 is
   !   -Omega^2 phi_D + g d phi_D/dz = Q_D,   Q_D = Q[tot_1, tot_2] - Q[inc_1, inc_2],
   ! with Q of SecondOrderIncident (free_surface_forcing): the incident-incident part is carried by the incident
   ! second-order potential, so only the difference of the total and the incident forcing enters here.
   !
   ! The first-order fields the integral needs (Task 6): at the inner free-surface points, the total phi, grad phi and
   ! phi_zz = -lap_h phi from phi at four horizontal offsets (fs_inner_fields; offsets shrunk next to a waterline,
   ! fs_inner_steps); outside the circle r = R, the partial waves of the scattered potential (fs_circle_modes), from samples
   ! on the cylinder by the all-centroid route of the boundary integral formula (fs_sample_phi), so that the samples are a
   ! sum of the free-surface modes, with a check that the modes give the samples back (fs_mode_miss, fs_mode_warning).

   USE SecondOrderIncident, ONLY: free_surface_forcing
   USE NearFieldDrift, ONLY: KERNEL_DIFF
   USE SecondOrderTypes, ONLY: panel_mesh, first_order_state, cs_points, cs_fields, CS_FREE
   USE PartialWaves, ONLY: pw_set, pw_radial, pw_decompose, pw_eval_point, pw_radial_at, pw_bessel_k_scaled, &
                           pw_bessel_k_scaled_c, pw_hankel_scaled
   USE FieldVelocity, ONLY: incident_at
   USE FirstOrderState, ONLY: wave_gravity
   USE WaterlineField, ONLY: field_influence_row, field_potential
   USE QTFField, ONLY: body_winding
   USE HAMS_mod, ONLY: NTHREAD
   USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: fs_point_field, fs_forcing_diff
   PUBLIC :: fs_inner_fields, fs_inner_steps, fs_sample_phi, fs_circle_modes, fs_mode_miss, fs_mode_warning, fs_auto_modes
   PUBLIC :: FS_MODE_TOL, fs_mode_scales, fs_auto_evanescent
   PUBLIC :: FS_SIGN, fs_inner_pair, fs_outer_pair, fs_outer_annulus, fs_outer_tail, fs_outer_ring, fs_outer_rfar
   PUBLIC :: fs_outer_sliver, fs_radius_message, fs_pair_order_message
   PUBLIC :: fs_pair_terms, fs_coverage_check

   ! The mode check: the partial waves must give the samples on the circle to this fraction of their largest value
   REAL(8), PARAMETER :: FS_MODE_TOL = 1D-3
   ! The smallest scale of the mode check, as a fraction of the largest sample over all states (fs_mode_scales)
   REAL(8), PARAMETER :: FS_SCALE_FLOOR = 1D-6
   ! The coverage check of a mesh from FSInnerMesh.pnl (fs_coverage_check): a warning when it leaves out (or counts twice)
   ! more than FS_COVER_WARN of the free surface inside R, a stop above FS_COVER_STOP either way
   REAL(8), PARAMETER :: FS_COVER_WARN = 1D-2, FS_COVER_STOP = 1D-1

   ! The sign s_FS of T_FS = s_FS (i Omega rho / g) int psi_k Q_D dS: +1, settled by the direct solve of the Task 2
   ! reference (fs_sign_from_bvp: 7e-6 with +1, 1.25 with -1)
   REAL(8), PARAMETER :: FS_SIGN = 1D0
   COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
   REAL(8), PARAMETER :: PI = 3.14159265358979323846D0
   ! The outer integral (Task 7): Gauss-Legendre points per radial interval, Gauss-Laguerre points per tail contour, the
   ! decay e^-EV_DROP below which an evanescent mode is dropped, ln(1e12) and the periods of the slowest phase class
   ! that set r_far (fs_outer_rfar), and the points per chord of the sliver (fs_outer_sliver)
   INTEGER, PARAMETER :: NGL = 8, NLAG = 48, NSL_TH = 8, NSL_R = 4
   REAL(8), PARAMETER :: EV_DROP = 40D0, RFAR_DECAY = 27.631021115928547D0, RFAR_WAVES = 20D0

   TYPE :: outer_data
      ! One pair of the outer integral (outer_setup): modes of 1, 2 and psi; plane-wave modes |n| <= ni1, ni2
      INTEGER :: m1 = 0, m2 = 0, mp = 0, nk = 0, n1 = 0, n2 = 0, np = 0, ni1 = 0, ni2 = 0, mx = 0
      REAL(8) :: w1 = 0D0, w2 = 0D0, g = 0D0, k1 = 0D0, k2 = 0D0, kp = 0D0, beta1 = 0D0, beta2 = 0D0
      REAL(8) :: radius = 0D0, h = 0D0, centre(2) = 0D0, dzi1 = 0D0, dzi2 = 0D0
      LOGICAL :: inc1 = .FALSE., inc2 = .FALSE.
      COMPLEX(8) :: amp1 = (0D0, 0D0), amp2c = (0D0, 0D0)
      REAL(8), ALLOCATABLE :: kev1(:), kev2(:), kevp(:), d1(:), e1(:), d2(:), e2(:), dp(:), ep(:)
      COMPLEX(8), ALLOCATABLE :: a1(:, :, :), a2(:, :, :), a2c(:, :, :), ap(:, :, :)
      REAL(8), ALLOCATABLE :: kr1(:, :), kr2(:, :), krp(:, :)
   END TYPE outer_data

   TYPE :: outer_work
      ! The work arrays of one thread of the outer integral (work_alloc)
      COMPLEX(8), ALLOCATABLE :: h1(:), h2(:), hp(:), j1(:), j2(:), kc(:)
      REAL(8), ALLOCATABLE :: kq1(:, :), dkq1(:, :), kq2(:, :), dkq2(:, :), kqp(:, :), dkqp(:, :)
      COMPLEX(8), ALLOCATABLE :: f1(:, :), c2(:, :), ps(:, :), tmp(:, :), fi(:, :), ci(:, :)
   END TYPE outer_work

   TYPE :: fs_point_field
      ! One first-order state at one point of z = 0: potential, gradient, and d^2 phi / dz^2.
      COMPLEX(8) :: phi = (0D0, 0D0)
      COMPLEX(8) :: grad(3) = (0D0, 0D0)
      COMPLEX(8) :: phizz = (0D0, 0D0)
   END TYPE fs_point_field

CONTAINS

   PURE COMPLEX(8) FUNCTION fs_forcing_diff(w1, tot1, inc1, w2, tot2, inc2, g)
      ! Q_D = Q[tot1, tot2] - Q[inc1, inc2] at the difference frequency (KERNEL_DIFF), for the states 1 (omega w1) and 2
      ! (omega w2), total and incident, at one point of z = 0.
      REAL(8), INTENT(IN) :: w1, w2, g
      TYPE(fs_point_field), INTENT(IN) :: tot1, inc1, tot2, inc2

      fs_forcing_diff = free_surface_forcing(KERNEL_DIFF, w1, tot1%phi, tot1%grad, tot1%phizz, &
                                             w2, tot2%phi, tot2%grad, tot2%phizz, g) &
                        - free_surface_forcing(KERNEL_DIFF, w1, inc1%phi, inc1%grad, inc1%phizz, &
                                               w2, inc2%phi, inc2%grad, inc2%phizz, g)
   END FUNCTION fs_forcing_diff

   SUBROUTINE fs_inner_fields(mesh, pts, states, step, is_multi, f, delta)
      ! The total first-order fields of every state at the inner points of the free surface: phi, grad phi and phi_zz.
      ! phi_P comes from the boundary integral formula (stock route, fs_sample_phi) at each point x and at its four
      ! horizontal offsets x +- delta e_x, x +- delta e_y, all five with x as decision point (near or far, and the series of
      ! FINGREEN3D), so that no switch between them enters a difference. Then, second order in delta [derived: central
      ! differences],
      !   d_x phi_P = [phi_P(x + delta e_x) - phi_P(x - delta e_x)] / (2 delta),   d_y phi_P likewise,
      !   lap_h phi_P = [phi_P(x + delta e_x) + phi_P(x - delta e_x) + phi_P(x + delta e_y) + phi_P(x - delta e_y)
      !                  - 4 phi_P(x)] / delta^2,
      ! and phi_P,zz = -lap_h phi_P (Laplace). d_z phi_P is nu phi_P at a point from z = -delta up (the linear free-surface
      ! condition, as fields_green), and the central difference of phi_P at x -+ delta e_z further down. 5 rows per point
      ! at z = 0 (7 below). The incident part is added exactly: phi_I, grad phi_I and phi_I,zz = k^2 phi_I.
      ! The rounding of the Green function enters lap_h divided by delta^2, so delta must not be too small: the driver
      ! passes step = 2 cs_fd_step, with which lap_h is that of the difference of central differences over delta/2
      ! (fields_green at cs_fd_step) to rounding, at 5 rows per point in place of 25 (fs_laplace_residual).
      ! delta(pts%n), if given, is the offset of each point; otherwise fs_inner_steps gives it: step x the point size, smaller
      ! next to a waterline. All states must be at the current global frequency.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(cs_points), INTENT(IN) :: pts
      TYPE(first_order_state), INTENT(IN) :: states(:)
      REAL(8), INTENT(IN) :: step
      LOGICAL, INTENT(IN) :: is_multi
      TYPE(fs_point_field), INTENT(OUT) :: f(pts%n, SIZE(states))
      REAL(8), INTENT(IN), OPTIONAL :: delta(:)
      ! The stencil: the centre, +x, -x, +y, -y, and below z = -delta also +z, -z
      REAL(8), PARAMETER :: OFF(3, 7) = RESHAPE([0D0, 0D0, 0D0, 1D0, 0D0, 0D0, -1D0, 0D0, 0D0, 0D0, 1D0, 0D0, 0D0, -1D0, 0D0, &
                                                 0D0, 0D0, 1D0, 0D0, 0D0, -1D0], [3, 7])
      TYPE(cs_points) :: ex
      REAL(8), ALLOCATABLE :: d(:), dec(:, :)
      INTEGER, ALLOCATABLE :: first(:), nst(:)
      COMPLEX(8), ALLOCATABLE :: ph(:, :)
      COMPLEX(8) :: phi_i, grad_i(3), p(7)
      INTEGER :: i, j, ist, k, nshrunk

      IF (step <= 0D0) THEN
         WRITE (*, '(a)') ' fs_inner_fields: the offset step must be positive.'
         ERROR STOP 1
      END IF
      ALLOCATE (d(pts%n))
      IF (PRESENT(delta)) THEN
         IF (SIZE(delta) /= pts%n) THEN
            WRITE (*, '(a)') ' fs_inner_fields: delta must have one offset per point.'
            ERROR STOP 1
         END IF
         d = delta
      ELSE
         CALL fs_inner_steps(mesh, pts, step, d, nshrunk)
      END IF
      IF (pts%n == 0) RETURN

      ! The stencil points of all inner points, in one list: point i has nst(i) of them from first(i)
      ALLOCATE (first(pts%n), nst(pts%n))
      k = 0
      DO i = 1, pts%n
         first(i) = k + 1
         nst(i) = MERGE(5, 7, pts%x(3, i) >= -d(i))
         k = k + nst(i)
      END DO
      ex%n = k
      ALLOCATE (ex%x(3, ex%n), ex%nrm(3, ex%n), ex%w(ex%n), ex%psize(ex%n), ex%kind(ex%n), ex%body(ex%n), dec(3, ex%n))
      DO i = 1, pts%n
         DO j = 1, nst(i)
            ex%x(:, first(i) + j - 1) = pts%x(:, i) + d(i)*OFF(:, j)
            dec(:, first(i) + j - 1) = pts%x(:, i)
         END DO
      END DO
      ex%nrm = 0D0
      ex%nrm(3, :) = -1D0
      ex%w = 0D0
      ex%psize = 0D0
      ex%kind = CS_FREE
      ex%body = 0
      ALLOCATE (ph(ex%n, SIZE(states)))
      CALL fs_sample_phi(mesh, ex, states, is_multi, ph, decide=dec)

      DO ist = 1, SIZE(states)
         DO i = 1, pts%n
            p(1:nst(i)) = ph(first(i):first(i) + nst(i) - 1, ist)
            CALL incident_at(states(ist), pts%x(:, i), phi_i, grad_i)
            f(i, ist)%phi = p(1) + phi_i
            f(i, ist)%grad(1) = (p(2) - p(3))/(2D0*d(i)) + grad_i(1)
            f(i, ist)%grad(2) = (p(4) - p(5))/(2D0*d(i)) + grad_i(2)
            IF (nst(i) == 7) THEN
               f(i, ist)%grad(3) = (p(6) - p(7))/(2D0*d(i)) + grad_i(3)
            ELSE
               ! Each part exactly, as fields_green: nu = omega^2/g of the solver for phi_P, the incident wave its own
               f(i, ist)%grad(3) = states(ist)%nu*p(1) + grad_i(3)
            END IF
            f(i, ist)%phizz = -(p(2) + p(3) + p(4) + p(5) - 4D0*p(1))/d(i)**2 + states(ist)%k**2*phi_i
         END DO
      END DO
   END SUBROUTINE fs_inner_fields

   SUBROUTINE fs_inner_steps(mesh, pts, step, delta, nshrunk, nfloor)
      ! The offset delta of each inner point for fs_inner_fields: step x its size, halved until two conditions hold, but not
      ! below DELTA_FLOOR x step x its size (Green-function rounding divided by delta^2 would take over):
      !   - no point of the stencil, x +- delta e_x and x +- delta e_y, lies in a waterplane (body_winding of QTFField 1/2
      !     or more for a body), where the boundary integral formula gives no fluid value;
      !   - nor does a point at CLEARANCE delta from x in any of 8 directions (axes and diagonals): next to the hull the
      !     field varies on the scale of the distance to it, so a stencil that only just clears the waterplane differences
      !     it badly (fs_waterline_offsets: phi_zz 23 % off at delta = 0.25 m, 0.3 m from the waterline, against 1.4e-3
      !     with the clearance, delta = 0.016 m).
      ! nshrunk counts the points whose delta was shrunk, nfloor those held at the floor (still clear of the waterplane,
      ! but closer than the clearance asks). A stencil still in a waterplane at the floor stops the run, naming the point.
      ! Frequency-independent: the driver calls it once.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(cs_points), INTENT(IN) :: pts
      REAL(8), INTENT(IN) :: step
      REAL(8), INTENT(OUT) :: delta(pts%n)
      INTEGER, INTENT(OUT) :: nshrunk
      INTEGER, INTENT(OUT), OPTIONAL :: nfloor
      REAL(8), PARAMETER :: CLEARANCE = 16D0, DELTA_FLOOR = 1D-3
      REAL(8) :: d, d0, dmin
      INTEGER :: i, nf, ibad
      LOGICAL :: in_wp, at_floor

      nshrunk = 0
      nf = 0
      ibad = 0
      !$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(d, d0, dmin, in_wp, at_floor) REDUCTION(+:nshrunk, nf) &
      !$OMP REDUCTION(MAX:ibad) SCHEDULE(DYNAMIC, 4)
      DO i = 1, pts%n
         d0 = step*pts%psize(i)
         dmin = DELTA_FLOOR*d0
         d = d0
         at_floor = .FALSE.
         DO
            in_wp = stencil_in_waterplane(mesh, pts%x(:, i), d)
            IF (.NOT. (in_wp .OR. clearance_in_waterplane(mesh, pts%x(:, i), CLEARANCE*d))) EXIT
            IF (at_floor) EXIT
            d = 0.5D0*d
            IF (d <= dmin) THEN
               d = dmin
               at_floor = .TRUE.
            END IF
         END DO
         IF (in_wp) ibad = MAX(ibad, i)
         IF (d < d0) nshrunk = nshrunk + 1
         IF (at_floor) nf = nf + 1
         delta(i) = d
      END DO
      !$OMP END PARALLEL DO
      IF (ibad > 0) THEN
         WRITE (*, '(a,I0,a,3F12.4,a)') ' QTF free surface: the difference stencil of inner point ', ibad, ' at (', &
            pts%x(:, ibad), ') lies in a waterplane even at the smallest offset; the point is on a waterline. Check the'// &
            ' inner free-surface mesh (FSInnerMesh.pnl or qtf_fs_rows).'
         ERROR STOP 1
      END IF
      IF (PRESENT(nfloor)) nfloor = nf
   END SUBROUTINE fs_inner_steps

   LOGICAL FUNCTION stencil_in_waterplane(mesh, x, d)
      ! A point of the stencil of fs_inner_fields at x with offset d, x +- d e_x and x +- d e_y, lies in a waterplane
      ! (body_winding >= 1/2)
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: x(3), d
      REAL(8), PARAMETER :: OFF(2, 4) = RESHAPE([1D0, 0D0, -1D0, 0D0, 0D0, 1D0, 0D0, -1D0], [2, 4])

      stencil_in_waterplane = any_in_waterplane(mesh, x, d, OFF)
   END FUNCTION stencil_in_waterplane

   LOGICAL FUNCTION clearance_in_waterplane(mesh, x, r)
      ! A point at distance r from x along the axes or the diagonals lies in a waterplane
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: x(3), r
      REAL(8), PARAMETER :: S = 0.70710678118654752D0
      REAL(8), PARAMETER :: OFF(2, 8) = RESHAPE([1D0, 0D0, -1D0, 0D0, 0D0, 1D0, 0D0, -1D0, S, S, S, -S, -S, S, -S, -S], [2, 8])

      clearance_in_waterplane = any_in_waterplane(mesh, x, r, OFF)
   END FUNCTION clearance_in_waterplane

   LOGICAL FUNCTION any_in_waterplane(mesh, x, d, off)
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: x(3), d, off(:, :)
      REAL(8), ALLOCATABLE :: wind(:)
      INTEGER :: j

      any_in_waterplane = .FALSE.
      DO j = 1, SIZE(off, 2)
         CALL body_winding(mesh, [x(1) + d*off(1, j), x(2) + d*off(2, j), x(3)], wind)
         IF (MAXVAL(wind) >= 0.5D0) THEN
            any_in_waterplane = .TRUE.
            RETURN
         END IF
      END DO
   END FUNCTION any_in_waterplane

   SUBROUTINE fs_sample_phi(mesh, pts, states, is_multi, phi, centroid_only, decide)
      ! The perturbation potential phi_P of every state at every point of pts, any kind and any z <= 0, from one row of the
      ! boundary integral formula per point (field_influence_row in field mode); no gradient. phi(pts%n, SIZE(states)).
      ! With centroid_only every panel takes its centroid value (near_panel = .FALSE.), so that phi_P is a sum of the
      ! free-surface modes outside the hull: the samples of the partial waves on the circle. Otherwise the stock route of
      ! fields_green. decide(3, pts%n), if given, is the decision point of each point (near or far, series of FINGREEN3D),
      ! else the point itself: the points of a difference stencil pass their centre (fs_inner_fields). All states must be
      ! at the current global frequency.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(cs_points), INTENT(IN) :: pts
      TYPE(first_order_state), INTENT(IN) :: states(:)
      LOGICAL, INTENT(IN) :: is_multi
      COMPLEX(8), INTENT(OUT) :: phi(pts%n, SIZE(states))
      LOGICAL, INTENT(IN), OPTIONAL :: centroid_only
      REAL(8), INTENT(IN), OPTIONAL :: decide(:, :)
      COMPLEX(8), ALLOCATABLE :: row_g(:, :), row_dg(:, :)
      REAL(8) :: xd(3)
      LOGICAL :: near
      INTEGER :: i, ist

      near = .TRUE.
      IF (PRESENT(centroid_only)) near = .NOT. centroid_only
      IF (PRESENT(decide)) THEN
         IF (SIZE(decide, 1) /= 3 .OR. SIZE(decide, 2) /= pts%n) THEN
            WRITE (*, '(a)') ' fs_sample_phi: decide must have one point per point of pts.'
            ERROR STOP 1
         END IF
      END IF
      ! As fields_green: the points are independent, the rows are per thread, the panel loop runs on one thread here
      !$OMP PARALLEL NUM_THREADS(NTHREAD) PRIVATE(row_g, row_dg, ist, xd)
      ALLOCATE (row_g(mesh%nbase, mesh%nsys), row_dg(mesh%nbase, mesh%nsys))
      !$OMP DO SCHEDULE(DYNAMIC, 4)
      DO i = 1, pts%n
         xd = pts%x(:, i)
         IF (PRESENT(decide)) xd = decide(:, i)
         CALL field_influence_row(mesh, pts%x(:, i), is_multi, row_g, row_dg, xd, field_mode=.TRUE., near_panel=near)
         DO ist = 1, SIZE(states)
            phi(i, ist) = field_potential(mesh, states(ist), row_g, row_dg)
         END DO
      END DO
      !$OMP END DO
      DEALLOCATE (row_g, row_dg)
      !$OMP END PARALLEL
   END SUBROUTINE fs_sample_phi

   SUBROUTINE fs_circle_modes(mesh, circle, states, is_multi, centre, radius, h, k0, kev, mmax, nev, pw, miss)
      ! The partial waves of the perturbation potential of every state outside the circle r = radius about centre: phi_P at
      ! the samples circle of pw_sample_points(centre, radius, h, mmax, nev), by the centroid route (fs_sample_phi), so that
      ! the samples are a sum of the free-surface modes, then pw_decompose with the propagating wave number k0 and the
      ! evanescent kev(1:nev) (pw_evanescent_from_wvn) of the current frequency. miss is the mode check of each state
      ! (fs_mode_miss), relative to fs_mode_scales: each state's own largest sample, so that every state is held to the
      ! tolerance of itself. The samples are evaluated in the fluid whatever their kind (CS_WALL from pw_sample_points).
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(cs_points), INTENT(IN) :: circle
      TYPE(first_order_state), INTENT(IN) :: states(:)
      LOGICAL, INTENT(IN) :: is_multi
      REAL(8), INTENT(IN) :: centre(2), radius, h, k0, kev(:)
      INTEGER, INTENT(IN) :: mmax, nev
      TYPE(pw_set), INTENT(OUT) :: pw(SIZE(states))
      REAL(8), INTENT(OUT) :: miss(SIZE(states))
      COMPLEX(8), ALLOCATABLE :: samples(:, :)
      REAL(8) :: scale(SIZE(states))
      INTEGER :: ist

      ALLOCATE (samples(circle%n, SIZE(states)))
      CALL fs_sample_phi(mesh, circle, states, is_multi, samples, centroid_only=.TRUE.)
      scale = fs_mode_scales(samples)
      DO ist = 1, SIZE(states)
         CALL pw_decompose(samples(:, ist), centre, radius, h, k0, kev, mmax, nev, pw(ist))
         miss(ist) = fs_mode_miss(pw(ist), circle, samples(:, ist), scale(ist))
      END DO
   END SUBROUTINE fs_circle_modes

   PURE FUNCTION fs_mode_scales(samples) RESULT(scale)
      ! The scale of the mode check of each state (column) of samples: its own largest |sample|, but at least FS_SCALE_FLOOR
      ! of the largest over all states, so that a state that is zero up to rounding (psi of yaw on an axisymmetric body) is
      ! not held to 1e-3 of its rounding, while every real state is held to 1e-3 of itself.
      COMPLEX(8), INTENT(IN) :: samples(:, :)
      REAL(8) :: scale(SIZE(samples, 2))
      REAL(8) :: big
      INTEGER :: k

      big = 0D0
      IF (SIZE(samples) > 0) big = MAXVAL(ABS(samples))
      DO k = 1, SIZE(samples, 2)
         scale(k) = 0D0
         IF (SIZE(samples, 1) > 0) scale(k) = MAXVAL(ABS(samples(:, k)))
         scale(k) = MAX(scale(k), FS_SCALE_FLOOR*big)
      END DO
   END FUNCTION fs_mode_scales

   REAL(8) FUNCTION fs_mode_miss(pw, pts, samples, scale)
      ! The mode check: max |phi_pw(x_i) - samples(i)| over the points of pts (the samples on the circle, at depth), relative
      ! to scale if given and positive, else to max |samples(i)|; phi_pw by pw_eval_point with the radial functions at r = R
      ! computed once. HUGE if a value is not finite; the absolute miss if the scale is 0.
      TYPE(pw_set), INTENT(IN) :: pw
      TYPE(cs_points), INTENT(IN) :: pts
      COMPLEX(8), INTENT(IN) :: samples(:)
      REAL(8), INTENT(IN), OPTIONAL :: scale
      TYPE(pw_radial) :: rad
      COMPLEX(8) :: phi
      REAL(8) :: worst, big, e
      INTEGER :: i

      CALL pw_radial_at(pw, pw%radius, rad)
      worst = 0D0
      big = 0D0
      DO i = 1, pts%n
         CALL pw_eval_point(pw, pts%x(:, i), phi, rad=rad)
         e = ABS(phi - samples(i))
         IF (.NOT. (IEEE_IS_FINITE(e) .AND. IEEE_IS_FINITE(ABS(samples(i))))) THEN
            fs_mode_miss = HUGE(1D0)
            RETURN
         END IF
         worst = MAX(worst, e)
         big = MAX(big, ABS(samples(i)))
      END DO
      IF (PRESENT(scale)) THEN
         IF (scale > 0D0) big = scale
      END IF
      fs_mode_miss = worst
      IF (big > 0D0) fs_mode_miss = worst/big
   END FUNCTION fs_mode_miss

   FUNCTION fs_mode_warning(miss) RESULT(msg)
      ! The log line of the mode check when miss (fs_mode_miss) is above FS_MODE_TOL, else ''
      REAL(8), INTENT(IN) :: miss
      CHARACTER(LEN=:), ALLOCATABLE :: msg
      CHARACTER(LEN=160) :: line

      msg = ''
      IF (miss <= FS_MODE_TOL) RETURN
      WRITE (line, '(a,ES9.3,a)') 'QTF free surface: the partial waves miss the field on the circle by ', miss, &
         '; raise qtf_fs_modes or qtf_fs_evanescent'
      msg = TRIM(line)
   END FUNCTION fs_mode_warning

   INTEGER FUNCTION fs_auto_modes(modes, k0_max, radius)
      ! The angular modes M of the partial waves: qtf_fs_modes (modes) if positive, else ceil(k0_max R) + 10, with k0_max the
      ! largest first-order wave number of the run (Jacobi-Anger: J_m(kR) is negligible for m well above kR)
      INTEGER, INTENT(IN) :: modes
      REAL(8), INTENT(IN) :: k0_max, radius

      IF (modes > 0) THEN
         fs_auto_modes = modes
      ELSE
         fs_auto_modes = CEILING(k0_max*radius) + 10
      END IF
   END FUNCTION fs_auto_modes

   INTEGER FUNCTION fs_auto_evanescent(nev, h, radius, r_wl)
      ! The evanescent modes N_ev of the partial waves: qtf_fs_evanescent (nev) if positive, else
      ! max(8, ceil(9 h / (pi (R - r_wl)))), r_wl the largest distance from the centre to a waterline node. The slowest
      ! evanescent mode kept, k_N > (N - 1/2) pi / h ~ N pi / h, then decays by about e^-9 (1e-4) from the waterline to R,
      ! and the first one left out decays faster still (ruling of the Task 6 fix round). -1 when R <= r_wl.
      INTEGER, INTENT(IN) :: nev
      REAL(8), INTENT(IN) :: h, radius, r_wl

      IF (nev > 0) THEN
         fs_auto_evanescent = nev
      ELSE IF (.NOT. radius > r_wl) THEN
         ! The circle cuts a body (fs_radius_message stops the run): no division by R - r_wl <= 0
         fs_auto_evanescent = -1
      ELSE
         fs_auto_evanescent = MAX(8, CEILING(9D0*h/(PI*(radius - r_wl))))
      END IF
   END FUNCTION fs_auto_evanescent


   ! ---------------------------------------------------------------------------------------- the integrals (Task 7)

   FUNCTION fs_radius_message(radius, r_wl) RESULT(msg)
      ! The stop message when the circle R does not clear the waterlines (R <= r_wl, the largest distance from the centre
      ! to a waterline node; the equality too, since R - r_wl sets the evanescent modes), else ''.
      REAL(8), INTENT(IN) :: radius, r_wl
      CHARACTER(LEN=:), ALLOCATABLE :: msg
      CHARACTER(LEN=300) :: line

      msg = ''
      IF (radius > r_wl) RETURN
      WRITE (line, '(a,F0.4,a,F0.4,a)') 'the circle R = ', radius, ' m does not exceed the waterline radius r_wl = ', r_wl, &
         ' m (the largest distance from the centre to a waterline node): it would cut a body; raise qtf_fs_radius.'
      msg = TRIM(line)
   END FUNCTION fs_radius_message

   FUNCTION fs_pair_order_message(w1, w2, omega) RESULT(msg)
      ! The stop message of the outer integral when the pair is not ordered w1 > w2 with Omega = omega > 0 (the
      ! difference frequency; psi_k enters as H^(1)(K r), so the phase class k1 - k2 + K of the tail must be positive),
      ! else ''. Omega = 0 is not refused here: the routines return 0 for it before asking.
      REAL(8), INTENT(IN) :: w1, w2, omega
      CHARACTER(LEN=:), ALLOCATABLE :: msg
      CHARACTER(LEN=300) :: line

      msg = ''
      IF (w1 > w2 .AND. omega > 0D0) RETURN
      WRITE (line, '(a,ES12.5,a,ES12.5,a,ES12.5,a)') 'the outer free-surface integral needs the pair ordered w1 > w2 '// &
         'with Omega > 0 (w1 = ', w1, ', w2 = ', w2, ', Omega = ', omega, '): the caller must order the pair.'
      msg = TRIM(line)
   END FUNCTION fs_pair_order_message

   SUBROUTINE fs_pair_terms(w1, f1, inc1, pw1, beta1, k1, w2, f2, inc2, pw2, beta2, k2, psi, pw_psi, wts, np, g, rho, &
                            t_in, t_out, r_far)
      ! The free-surface terms of T_pot of the pair of states 1 (omega w1, heading beta1 in radians, incident wave number
      ! k1 = pw1%k0) and 2, as QTFDriver adds them (qtf_terms = full):
      !   t_in   fs_inner_pair over the inner points: total fields f1, f2, incident fields inc1, inc2 (phi_zz = k^2 phi),
      !          psi_k there psi(n, 6 nbody), panel areas wts;
      !   t_out  fs_outer_pair from R to fs_outer_rfar and its tail, plus fs_outer_sliver over the gap between the polygon of
      !          np chords that a generated inner mesh ends on and the circle (np < 3, a mesh from FSInnerMesh.pnl: none).
      ! psi and pw_psi are psi_k at |Omega|, Omega = w1 - w2, as the store keeps them. r_far, if present, is the r_far used
      ! (0 if none).
      !   Omega = 0 (w1 == w2, the diagonal of the QTF): no free-surface term; zeros, and no outer routine is called
      !   (fs_outer_rfar and fs_outer_ring stop on Omega = 0).
      !   w1 < w2: the pair is evaluated as (2, 1), at -Omega > 0 with the same psi, and conjugated, T(1, 2) = conj T(2, 1)
      !   [derived]: Q_D(2, 1) = conj Q_D(1, 2) (pair_product of KERNEL_DIFF) and psi(-Omega) = conj psi(Omega), as
      !   potential_pair takes for the hull terms (fs_inner_pair_hermitian checks it on the inner sum). So the outer
      !   routines only see ordered pairs.
      REAL(8), INTENT(IN) :: w1, beta1, k1, w2, beta2, k2, wts(:), g, rho
      TYPE(fs_point_field), INTENT(IN) :: f1(:), inc1(:), f2(:), inc2(:)
      TYPE(pw_set), INTENT(IN) :: pw1, pw2, pw_psi(:)
      COMPLEX(8), INTENT(IN) :: psi(:, :)
      INTEGER, INTENT(IN) :: np
      COMPLEX(8), INTENT(OUT) :: t_in(:), t_out(:)
      REAL(8), INTENT(OUT), OPTIONAL :: r_far

      t_in = (0D0, 0D0)
      t_out = (0D0, 0D0)
      IF (PRESENT(r_far)) r_far = 0D0
      IF (w1 == w2) RETURN
      IF (w1 > w2) THEN
         CALL ordered(w1, f1, inc1, pw1, beta1, k1, w2, f2, inc2, pw2, beta2, k2)
      ELSE
         CALL ordered(w2, f2, inc2, pw2, beta2, k2, w1, f1, inc1, pw1, beta1, k1)
         t_in = CONJG(t_in)
         t_out = CONJG(t_out)
      END IF

   CONTAINS

      SUBROUTINE ordered(wa, fa, ia, pa, ba, ka, wb, fb, ib, pb, bb, kb)
         ! The terms of the pair (a, b) with wa > wb
         REAL(8), INTENT(IN) :: wa, ba, ka, wb, bb, kb
         TYPE(fs_point_field), INTENT(IN) :: fa(:), ia(:), fb(:), ib(:)
         TYPE(pw_set), INTENT(IN) :: pa, pb
         COMPLEX(8) :: sl(SIZE(t_out))
         REAL(8) :: rf

         CALL fs_inner_pair(wa, fa, ia, wb, fb, ib, psi, wts, g, rho, t_in)
         rf = fs_outer_rfar(pa, ka, pb, kb, pw_psi)
         CALL fs_outer_pair(wa, pa, ba, ka, wb, pb, bb, kb, pw_psi, wa - wb, g, rho, rf, t_out)
         IF (np >= 3) THEN
            CALL fs_outer_sliver(wa, pa, ba, ka, wb, pb, bb, kb, pw_psi, wa - wb, g, rho, np, sl)
            t_out = t_out + sl
         END IF
         IF (PRESENT(r_far)) r_far = rf
      END SUBROUTINE ordered

   END SUBROUTINE fs_pair_terms

   SUBROUTINE fs_coverage_check(w, radius, awp, miss, msg, stops)
      ! The coverage of an inner free-surface mesh from FSInnerMesh.pnl (Task 7 review, minor 4). The outer integral starts
      ! at the circle r = R, and a file mesh gets no sliver (fs_outer_sliver needs the polygon of the generated mesh), so
      ! the file mesh must reach R. With w the panel areas and awp the waterplane area of all bodies,
      !   miss = (pi R^2 - awp - sum w) / (pi R^2 - awp),
      ! the fraction of the free surface inside R that the mesh leaves out (< 0: more area than there is, from overlapping
      ! panels or panels over a waterplane). |miss| > FS_COVER_WARN: msg is a warning; |miss| > FS_COVER_STOP: msg is a
      ! stop message and stops is true (the mesh does not reach R, or it counts part of the free surface twice). A mesh that ends on a polygon inscribed in the circle leaves
      ! out the polygon gap, about (2 pi / n)^2 / 6 of pi R^2 for n chords.
      REAL(8), INTENT(IN) :: w(:), radius, awp
      REAL(8), INTENT(OUT) :: miss
      CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: msg
      LOGICAL, INTENT(OUT) :: stops
      CHARACTER(LEN=400) :: line
      REAL(8) :: region

      region = PI*radius**2 - awp
      miss = (region - SUM(w))/region
      stops = ABS(miss) > FS_COVER_STOP
      msg = ''
      IF (stops .AND. miss < 0D0) THEN
         WRITE (line, '(a,F10.3,a,F8.3,a)') 'QTF free surface: FSInnerMesh.pnl has more area than the free surface inside'// &
            ' R =', radius, ' m, by', -100D0*miss, ' % (overlapping panels, or panels over a waterplane), so part of the'// &
            ' free surface would be counted twice. Correct the mesh.'
         msg = TRIM(line)
      ELSE IF (stops) THEN
         WRITE (line, '(a,F10.3,a,F8.3,a)') 'QTF free surface: FSInnerMesh.pnl does not reach the circle R =', radius, &
            ' m: it leaves out', 100D0*miss, ' % of the free surface inside R (the circle minus the waterplanes), which no'// &
            ' part of the integral covers. Extend the mesh to the circle.'
         msg = TRIM(line)
      ELSE IF (ABS(miss) > FS_COVER_WARN) THEN
         WRITE (line, '(a,F8.3,a,F10.3,a)') 'QTF free surface warning: the inner mesh of FSInnerMesh.pnl leaves out', &
            100D0*miss, ' % of the free surface inside R =', radius, ' m (the circle minus the waterplanes; < 0: it'// &
            ' counts some twice). The outer integral starts at R: extend the mesh to the circle, without overlaps.'
         msg = TRIM(line)
      END IF
   END SUBROUTINE fs_coverage_check

   SUBROUTINE fs_inner_pair(w1, f1, inc1, w2, f2, inc2, psi, wts, g, rho, t_in)
      ! The inner part of T_FS for the pair of states 1 (omega w1) and 2 (omega w2): the sum over the inner points of
      !   t_in(k) = s_FS (i Omega rho / g) sum_i wts(i) psi(i, k) Q_D(x_i),   Omega = w1 - w2,
      ! Q_D = fs_forcing_diff of the total fields f1, f2 and the incident fields inc1, inc2 at the points (fs_inner_fields,
      ! incident_at). psi(n, 6 nbody) holds psi_k at the points (or, for Phase 7b, any modal potential); wts the panel areas.
      REAL(8), INTENT(IN) :: w1, w2, g, rho, wts(:)
      TYPE(fs_point_field), INTENT(IN) :: f1(:), inc1(:), f2(:), inc2(:)
      COMPLEX(8), INTENT(IN) :: psi(:, :)
      COMPLEX(8), INTENT(OUT) :: t_in(:)
      COMPLEX(8) :: q
      INTEGER :: i, k

      IF (SIZE(f1) /= SIZE(wts) .OR. SIZE(f2) /= SIZE(wts) .OR. SIZE(inc1) /= SIZE(wts) .OR. SIZE(inc2) /= SIZE(wts) &
          .OR. SIZE(psi, 1) /= SIZE(wts) .OR. SIZE(t_in) /= SIZE(psi, 2)) THEN
         WRITE (*, '(a)') ' fs_inner_pair: the fields, psi, the weights and t_in do not match in size.'
         ERROR STOP 1
      END IF
      t_in = (0D0, 0D0)
      DO i = 1, SIZE(wts)
         q = wts(i)*fs_forcing_diff(w1, f1(i), inc1(i), w2, f2(i), inc2(i), g)
         DO k = 1, SIZE(t_in)
            t_in(k) = t_in(k) + psi(i, k)*q
         END DO
      END DO
      t_in = FS_SIGN*CI*(w1 - w2)*rho/g*t_in
   END SUBROUTINE fs_inner_pair

   SUBROUTINE fs_outer_pair(w1, pw1, beta1, k1, w2, pw2, beta2, k2, pw_psi, omega, g, rho, r_far, t_out)
      ! The outer part of T_FS, r > R, for the pair of states 1 and 2:
      !   t_out(k) = s_FS (i Omega rho / g) int_{r > R} psi_k Q_D dS,
      ! from the partial waves of the scattered potentials pw1 (at w1) and pw2 (at w2), the incident waves (heading beta1 in
      ! radians, wave number k1 = pw1%k0; k1 <= 0: no incident wave; the same for 2) and the partial waves pw_psi(6 nbody)
      ! of psi_k at Omega = omega (for Phase 7b, any set of modal potentials at Omega). In r > R
      !   Q_D = Q[I_1, S_2] + Q[S_1, I_2] + Q[S_1, S_2]
      ! (Q sesquilinear: free_surface_forcing with KERNEL_DIFF; the incident-incident part belongs to the incident
      ! second-order potential). Every part is a modal field: S by its partial waves, I by Jacobi-Anger,
      ! e^{ikr cos(theta - beta)} = sum_n i^n J_n(kr) e^{in(theta - beta)}, so that the angular integral is the exact sum over
      ! the triples of modes with a - b + c = 0 (fs_outer_ring); for a part with a plane wave this is
      ! 2 pi sum_n g_n i^n J_n(kr) e^{+in beta}. Radially: Gauss-Legendre from R to r_far (fs_outer_annulus), then the
      ! propagating tail from r_far to infinity on rotated contours per phase class (fs_outer_tail). r_far from
      ! fs_outer_rfar, or larger. Omega = 0 gives 0.
      REAL(8), INTENT(IN) :: w1, beta1, k1, w2, beta2, k2, omega, g, rho, r_far
      TYPE(pw_set), INTENT(IN) :: pw1, pw2, pw_psi(:)
      COMPLEX(8), INTENT(OUT) :: t_out(:)
      COMPLEX(8) :: tail(SIZE(t_out))

      IF (.NOT. r_far > pw1%radius) THEN
         WRITE (*, '(a)') ' fs_outer_pair: r_far must exceed the radius R of the partial waves.'
         ERROR STOP 1
      END IF
      CALL fs_outer_annulus(w1, pw1, beta1, k1, w2, pw2, beta2, k2, pw_psi, omega, g, rho, pw1%radius, r_far, t_out)
      CALL fs_outer_tail(w1, pw1, beta1, k1, w2, pw2, beta2, k2, pw_psi, omega, g, rho, r_far, tail)
      t_out = t_out + tail
   END SUBROUTINE fs_outer_pair

   REAL(8) FUNCTION fs_outer_rfar(pw1, k1, pw2, k2, pw_psi)
      ! The radius where the outer integral turns from the real axis to the tail contours: the slowest evanescent mode of
      ! all sets has decayed below 1e-12, e^{-k_1 (r_far - R)} = 1e-12 (the Task 2 reference: 5e-4 was left at 150 m beyond
      ! R in 100 m of water), and r_far - R is at least 20 periods 2 pi / kappa_min of the slowest phase class of the tail:
      ! kappa = k_1 - k_2 + K always, -k_1 - k_2 + K with incident wave 1 (k1 > 0), k_1 + k_2 + K with incident wave 2.
      ! Only the classes that occur count (k_1 - k_2 - K, for one, never does; it tends to 0 in shallow water).
      TYPE(pw_set), INTENT(IN) :: pw1, pw2, pw_psi(:)
      REAL(8), INTENT(IN) :: k1, k2
      REAL(8) :: kevmin, kmin, ka, kb, kp

      kevmin = HUGE(1D0)
      IF (pw1%nev > 0) kevmin = MIN(kevmin, pw1%kev(1))
      IF (pw2%nev > 0) kevmin = MIN(kevmin, pw2%kev(1))
      kp = 0D0
      IF (SIZE(pw_psi) > 0) THEN
         IF (pw_psi(1)%nev > 0) kevmin = MIN(kevmin, pw_psi(1)%kev(1))
         kp = pw_psi(1)%k0
      END IF
      ka = pw1%k0
      kb = pw2%k0
      IF (.NOT. ka > kb) THEN
         WRITE (*, '(a)') ' fs_outer_rfar: the pair must be ordered w1 > w2 (k1 > k2); the caller must order the pair.'
         ERROR STOP 1
      END IF
      kmin = ABS(ka - kb + kp)
      IF (k1 > 0D0) kmin = MIN(kmin, ABS(-ka - kb + kp))
      IF (k2 > 0D0) kmin = MIN(kmin, ABS(ka + kb + kp))
      fs_outer_rfar = pw1%radius + RFAR_WAVES*2D0*PI/kmin
      IF (kevmin < HUGE(1D0)) fs_outer_rfar = MAX(fs_outer_rfar, pw1%radius + RFAR_DECAY/kevmin)
   END FUNCTION fs_outer_rfar

   SUBROUTINE fs_outer_annulus(w1, pw1, beta1, k1, w2, pw2, beta2, k2, pw_psi, omega, g, rho, r_a, r_b, t)
      ! s_FS (i Omega rho / g) int_{r_a < r < r_b} psi_k Q_D r dr dtheta (the arguments of fs_outer_pair), r_a >= R: the
      ! ring integral of fs_outer_ring by Gauss-Legendre of order NGL on intervals of at most 1/8 of the shortest
      ! wavelength of the integrand, 2 pi / (k_1 + k_2 + K), at most 2 / k_n of the fastest evanescent mode still alive (a
      ! mode is dropped where it has decayed by e^-40), and at most 2 r / (M + 1). OpenMP over the intervals; each interval
      ! keeps its own sum, and these are added in order after the parallel loop, so that the result does not depend on the
      ! threads or the schedule (the same bits run to run).
      REAL(8), INTENT(IN) :: w1, beta1, k1, w2, beta2, k2, omega, g, rho, r_a, r_b
      TYPE(pw_set), INTENT(IN) :: pw1, pw2, pw_psi(:)
      COMPLEX(8), INTENT(OUT) :: t(:)
      TYPE(outer_data) :: od
      REAL(8), ALLOCATABLE :: edge(:)
      COMPLEX(8), ALLOCATABLE :: part(:, :)
      REAL(8) :: xg(NGL), wg(NGL)
      INTEGER :: ni, i

      t = (0D0, 0D0)
      IF (omega == 0D0 .OR. .NOT. r_b > r_a .OR. SIZE(t) == 0) RETURN
      CALL stop_if_unordered(w1, w2, omega)
      CALL outer_setup(w1, pw1, beta1, k1, w2, pw2, beta2, k2, pw_psi, g, od)
      IF (SIZE(t) /= od%nk) THEN
         WRITE (*, '(a)') ' fs_outer_annulus: t must have one entry per set of pw_psi.'
         ERROR STOP 1
      END IF
      IF (r_a < od%radius*(1D0 - 1D-12)) THEN
         WRITE (*, '(a)') ' fs_outer_annulus: r_a must not be inside the circle R of the partial waves.'
         ERROR STOP 1
      END IF
      CALL radial_intervals(od, r_a, r_b, edge, ni)
      CALL gauss_legendre01(NGL, xg, wg)
      ALLOCATE (part(od%nk, ni))
      !$OMP PARALLEL NUM_THREADS(NTHREAD) DEFAULT(SHARED)
      CALL annulus_worker(od, edge, ni, xg, wg, part)
      !$OMP END PARALLEL
      DO i = 1, ni
         t = t + part(:, i)
      END DO
      t = FS_SIGN*CI*omega*rho/g*t
   END SUBROUTINE fs_outer_annulus

   SUBROUTINE annulus_worker(od, edge, ni, xg, wg, part)
      ! The share of one thread of fs_outer_annulus (an orphaned OMP DO): its own work arrays; part(:, i) the sum over the
      ! Gauss points of interval i, each interval written by one thread only
      TYPE(outer_data), INTENT(IN) :: od
      REAL(8), INTENT(IN) :: edge(:), xg(:), wg(:)
      INTEGER, INTENT(IN) :: ni
      COMPLEX(8), INTENT(INOUT) :: part(:, :)
      TYPE(outer_work) :: wk
      COMPLEX(8), ALLOCATABLE :: ring(:), acc(:)
      REAL(8) :: r, dr
      INTEGER :: i, q

      CALL work_alloc(od, wk)
      ALLOCATE (ring(od%nk), acc(od%nk))
      !$OMP DO SCHEDULE(DYNAMIC, 8) PRIVATE(i, q, r, dr)
      DO i = 1, ni
         dr = edge(i + 1) - edge(i)
         acc = (0D0, 0D0)
         DO q = 1, SIZE(xg)
            r = edge(i) + xg(q)*dr
            CALL ring_real(od, r, ring, wk)
            acc = acc + (wg(q)*dr*r)*ring
         END DO
         part(:, i) = acc
      END DO
      !$OMP END DO
   END SUBROUTINE annulus_worker

   SUBROUTINE fs_outer_tail(w1, pw1, beta1, k1, w2, pw2, beta2, k2, pw_psi, omega, g, rho, r_far, t)
      ! s_FS (i Omega rho / g) int_{r > r_far} psi_k Q_D dS (the arguments of fs_outer_pair): the propagating modes only (the
      ! evanescent ones have decayed, fs_outer_rfar). Writing J_n = (H^(1)_n + H^(2)_n) / 2 for the plane waves, and
      ! conj(H^(1)(k_2 r)) = H^(2)(k_2 r) (analytic in r) for state 2, the ring integral is a sum of three phase classes
      ! e^{i kappa r} G(r), G slowly varying:
      !   kappa_a = k_1 - k_2 + K: S_1 S_2 psi, (H^(1) part of I_1) S_2 psi, S_1 (H^(2) part of I_2) psi;
      !   kappa_b = -k_1 - k_2 + K: (H^(2) part of I_1) S_2 psi;   kappa_c = k_1 + k_2 + K: S_1 (H^(1) part of I_2) psi.
      ! Each class is analytic in r off the origin, so (Cauchy, Jordan) its integral is that along r = r_far + i s t,
      ! s = sign(kappa), t > 0, where e^{i kappa r} = e^{i kappa r_far} e^{-|kappa| t}:
      !   int_{r_far}^inf e^{i kappa r} G r dr = i s e^{i kappa r_far} int_0^inf e^{-|kappa| t} G(r) r dt,
      ! by Gauss-Laguerre of order NLAG, with the exponential factors of the Hankel functions taken out (pw_hankel_scaled),
      ! as the Task 2 reference does. Checked against the real-axis integral in fs_outer_tail (r_far against 2 r_far).
      REAL(8), INTENT(IN) :: w1, beta1, k1, w2, beta2, k2, omega, g, rho, r_far
      TYPE(pw_set), INTENT(IN) :: pw1, pw2, pw_psi(:)
      COMPLEX(8), INTENT(OUT) :: t(:)
      TYPE(outer_data) :: od
      TYPE(outer_work) :: wk
      REAL(8) :: xl(NLAG), wl(NLAG), kap(3), s, tq
      COMPLEX(8) :: rq, acc(SIZE(t)), ring(SIZE(t))
      INTEGER :: cls, q

      t = (0D0, 0D0)
      IF (omega == 0D0 .OR. SIZE(t) == 0) RETURN
      CALL stop_if_unordered(w1, w2, omega)
      CALL outer_setup(w1, pw1, beta1, k1, w2, pw2, beta2, k2, pw_psi, g, od)
      IF (SIZE(t) /= od%nk) THEN
         WRITE (*, '(a)') ' fs_outer_tail: t must have one entry per set of pw_psi.'
         ERROR STOP 1
      END IF
      CALL work_alloc(od, wk)
      CALL gauss_laguerre(NLAG, xl, wl)
      kap = [od%k1 - od%k2 + od%kp, -od%k1 - od%k2 + od%kp, od%k1 + od%k2 + od%kp]
      DO cls = 1, 3
         IF (cls == 2 .AND. .NOT. od%inc1) CYCLE
         IF (cls == 3 .AND. .NOT. od%inc2) CYCLE
         IF (kap(cls) == 0D0) THEN
            WRITE (*, '(a)') ' fs_outer_tail: a phase class of the tail does not oscillate (kappa = 0); the integral diverges.'
            ERROR STOP 1
         END IF
         s = SIGN(1D0, kap(cls))
         acc = (0D0, 0D0)
         DO q = 1, NLAG
            tq = xl(q)/ABS(kap(cls))
            rq = CMPLX(r_far, s*tq, KIND=8)
            CALL ring_tail(od, cls, rq, ring, wk)
            acc = acc + (wl(q)*rq)*ring
         END DO
         t = t + CI*s*EXP(CI*kap(cls)*r_far)/ABS(kap(cls))*acc
      END DO
      t = FS_SIGN*CI*omega*rho/g*t
   END SUBROUTINE fs_outer_tail

   SUBROUTINE fs_outer_ring(w1, pw1, beta1, k1, w2, pw2, beta2, k2, pw_psi, g, r, ring)
      ! The ring integral ring(k) = int_0^2pi psi_k Q_D dtheta at radius r >= R (the integrand of the radial integral, without
      ! r and the factor s_FS i Omega rho / g), the arguments as fs_outer_pair. With the modal amplitudes at r of field 1,
      ! X_a e^{ia theta}, of conj(field 2), C_b e^{-ib theta}, and of psi_k, P_c e^{ic theta},
      !   ring = 2 pi sum_c P_c sum_a Q(a, b = a + c),
      ! Q(a, b) = free_surface_forcing (KERNEL_DIFF) of the mode pair (gradient (d_r, i a / r, d_z) of each; pair_product
      ! conjugates state 2, so it is given conj(C_b)). Exact for the band-limited modal fields: the same as sampling theta at
      ! enough points; for the plane waves, a runs over |a| <= M_2 + M_psi, all that the triples reach (Jacobi-Anger).
      REAL(8), INTENT(IN) :: w1, beta1, k1, w2, beta2, k2, g, r
      TYPE(pw_set), INTENT(IN) :: pw1, pw2, pw_psi(:)
      COMPLEX(8), INTENT(OUT) :: ring(:)
      TYPE(outer_data) :: od
      TYPE(outer_work) :: wk

      CALL stop_if_unordered(w1, w2, w1 - w2)
      CALL outer_setup(w1, pw1, beta1, k1, w2, pw2, beta2, k2, pw_psi, g, od)
      IF (SIZE(ring) /= od%nk) THEN
         WRITE (*, '(a)') ' fs_outer_ring: ring must have one entry per set of pw_psi.'
         ERROR STOP 1
      END IF
      CALL work_alloc(od, wk)
      CALL ring_real(od, r, ring, wk)
   END SUBROUTINE fs_outer_ring

   SUBROUTINE fs_outer_sliver(w1, pw1, beta1, k1, w2, pw2, beta2, k2, pw_psi, omega, g, rho, np, t)
      ! s_FS (i Omega rho / g) int psi_k Q_D dS over the gap between the circle R and the polygon of np chords inscribed in
      ! it with vertices at theta_j = 2 pi (j - 1) / np about the centre: the outer edge of the generated inner mesh
      ! (build_fs_inner, np_circle). The partial waves are used a little inside R there (at most R (1 - cos(pi/np))), where
      ! their series still converge. Per chord: Gauss-Legendre of order NSL_TH in theta and NSL_R in r from the chord to the
      ! circle; the fields pointwise from the modal amplitudes at r, Q_D by fs_forcing_diff. np < 3 (a mesh from
      ! FSInnerMesh.pnl, np_circle = 0): t = 0.
      REAL(8), INTENT(IN) :: w1, beta1, k1, w2, beta2, k2, omega, g, rho
      TYPE(pw_set), INTENT(IN) :: pw1, pw2, pw_psi(:)
      INTEGER, INTENT(IN) :: np
      COMPLEX(8), INTENT(OUT) :: t(:)
      TYPE(outer_data) :: od
      TYPE(outer_work) :: wk
      TYPE(fs_point_field) :: s1, s2, i1, i2, t1, t2
      REAL(8) :: xt(NSL_TH), wt(NSL_TH), xr(NSL_R), wr(NSL_R), th0, dth, th, rho_c, r, dxy(2)
      COMPLEX(8) :: q, e
      INTEGER :: j, it, ir, k, c

      t = (0D0, 0D0)
      IF (omega == 0D0 .OR. np < 3 .OR. SIZE(t) == 0) RETURN
      CALL stop_if_unordered(w1, w2, omega)
      CALL outer_setup(w1, pw1, beta1, k1, w2, pw2, beta2, k2, pw_psi, g, od)
      IF (SIZE(t) /= od%nk) THEN
         WRITE (*, '(a)') ' fs_outer_sliver: t must have one entry per set of pw_psi.'
         ERROR STOP 1
      END IF
      CALL work_alloc(od, wk)
      CALL gauss_legendre01(NSL_TH, xt, wt)
      CALL gauss_legendre01(NSL_R, xr, wr)
      dth = 2D0*PI/np
      DO j = 1, np
         th0 = dth*(j - 1)
         DO it = 1, NSL_TH
            th = th0 + xt(it)*dth
            rho_c = od%radius*COS(0.5D0*dth)/COS(th - th0 - 0.5D0*dth)
            DO ir = 1, NSL_R
               r = rho_c + xr(ir)*(od%radius - rho_c)
               CALL modes_real(od, r, wk)
               s1 = modal_point(od%m1, wk%f1, r, th, .FALSE.)
               s2 = modal_point(od%m2, wk%c2, r, th, .TRUE.)
               dxy = [r*COS(th), r*SIN(th)]
               i1 = plane_point(od%inc1, od%amp1, od%k1, od%beta1, od%dzi1, dxy)
               i2 = plane_point(od%inc2, CONJG(od%amp2c), od%k2, od%beta2, od%dzi2, dxy)
               t1%phi = s1%phi + i1%phi
               t1%grad = s1%grad + i1%grad
               t1%phizz = s1%phizz + i1%phizz
               t2%phi = s2%phi + i2%phi
               t2%grad = s2%grad + i2%grad
               t2%phizz = s2%phizz + i2%phizz
               q = fs_forcing_diff(od%w1, t1, i1, od%w2, t2, i2, od%g)*(wt(it)*dth*wr(ir)*(od%radius - rho_c)*r)
               DO k = 1, od%nk
                  e = (0D0, 0D0)
                  DO c = -od%mp, od%mp
                     e = e + wk%ps(c, k)*EXP(CI*(c*th))
                  END DO
                  t(k) = t(k) + e*q
               END DO
            END DO
         END DO
      END DO
      t = FS_SIGN*CI*omega*rho/g*t
   END SUBROUTINE fs_outer_sliver

   ! ------------------------------------------------------------------------------------------ outer integral, internals

   SUBROUTINE outer_setup(w1, pw1, beta1, k1, w2, pw2, beta2, k2, pw_psi, g, od)
      ! The data of the outer integral of one pair: mode counts, wave numbers, the incident amplitudes at the centre
      ! (-i g_w / omega e^{i k (c . e_beta)}, g_w = wave_gravity as incident_at), the vertical factors per mode
      ! (phi_z = d phi, phi_zz = e phi on z = 0: d = k0 tanh k0 h and -k_n tan k_n h, e = k0^2 and -k_n^2), and the
      ! evanescent amplitudes at R, a(m, n) K_m(k_n R), with e^{k_n R} K_m(k_n R), so that the radial factor
      ! K_m(k_n r) / K_m(k_n R) = [e^x K_m](k_n r) / [e^x K_m](k_n R) e^{-k_n (r - R)} needs no K_m beyond 600 however far out
      ! r is (review M3 of Task 5). Stops if the sets do not share the circle and depth, the psi sets their modes, or an
      ! incident wave number differs from that of its partial waves.
      REAL(8), INTENT(IN) :: w1, beta1, k1, w2, beta2, k2, g
      TYPE(pw_set), INTENT(IN) :: pw1, pw2, pw_psi(:)
      TYPE(outer_data), INTENT(OUT) :: od
      REAL(8) :: gw
      INTEGER :: k

      od%nk = SIZE(pw_psi)
      IF (od%nk == 0) THEN
         WRITE (*, '(a)') ' QTF free surface: the outer integral needs at least one set of partial waves of psi.'
         ERROR STOP 1
      END IF
      IF (.NOT. (same_circle(pw1, pw2) .AND. same_circle(pw1, pw_psi(1)))) THEN
         WRITE (*, '(a)') ' QTF free surface: the partial waves of the pair and of psi are not on the same circle and depth.'
         ERROR STOP 1
      END IF
      DO k = 2, od%nk
         IF (.NOT. same_circle(pw_psi(1), pw_psi(k)) .OR. pw_psi(k)%mmax /= pw_psi(1)%mmax .OR. &
             pw_psi(k)%nev /= pw_psi(1)%nev .OR. pw_psi(k)%k0 /= pw_psi(1)%k0) THEN
            WRITE (*, '(a)') ' QTF free surface: the partial waves of the psi_k do not share their modes.'
            ERROR STOP 1
         END IF
      END DO
      od%w1 = w1
      od%w2 = w2
      od%g = g
      od%m1 = pw1%mmax
      od%m2 = pw2%mmax
      od%mp = pw_psi(1)%mmax
      od%n1 = pw1%nev
      od%n2 = pw2%nev
      od%np = pw_psi(1)%nev
      od%k1 = pw1%k0
      od%k2 = pw2%k0
      od%kp = pw_psi(1)%k0
      od%h = pw1%h
      od%radius = pw1%radius
      od%centre = pw1%centre
      od%beta1 = beta1
      od%beta2 = beta2
      od%inc1 = k1 > 0D0
      od%inc2 = k2 > 0D0
      od%ni1 = od%m2 + od%mp
      od%ni2 = od%m1 + od%mp
      od%mx = MAX(od%m1, od%m2, od%mp, od%ni1, od%ni2)
      IF (od%inc1) THEN
         IF (ABS(k1 - od%k1) > 1D-9*k1) THEN
            WRITE (*, '(a)') ' QTF free surface: the incident wave number of state 1 is not that of its partial waves.'
            ERROR STOP 1
         END IF
         gw = wave_gravity(w1, od%k1, od%h)
         od%amp1 = -CI*gw/w1*EXP(CI*od%k1*(od%centre(1)*COS(beta1) + od%centre(2)*SIN(beta1)))
         od%dzi1 = od%k1*TANH(od%k1*od%h)
      END IF
      IF (od%inc2) THEN
         IF (ABS(k2 - od%k2) > 1D-9*k2) THEN
            WRITE (*, '(a)') ' QTF free surface: the incident wave number of state 2 is not that of its partial waves.'
            ERROR STOP 1
         END IF
         gw = wave_gravity(w2, od%k2, od%h)
         od%amp2c = CONJG(-CI*gw/w2*EXP(CI*od%k2*(od%centre(1)*COS(beta2) + od%centre(2)*SIN(beta2))))
         od%dzi2 = od%k2*TANH(od%k2*od%h)
      END IF
      ALLOCATE (od%a1(-od%m1:od%m1, 0:od%n1, 1), od%a2(-od%m2:od%m2, 0:od%n2, 1), od%a2c(-od%m2:od%m2, 0:od%n2, 1))
      ALLOCATE (od%ap(-od%mp:od%mp, 0:od%np, od%nk))
      CALL set_amplitudes(pw1, od%a1(:, :, 1), od%kev1, od%kr1, od%d1, od%e1)
      CALL set_amplitudes(pw2, od%a2(:, :, 1), od%kev2, od%kr2, od%d2, od%e2)
      DO k = 1, od%nk
         CALL set_amplitudes(pw_psi(k), od%ap(:, :, k), od%kevp, od%krp, od%dp, od%ep)
      END DO
      ! State 2 enters conjugated: in the tail its propagating amplitudes are conj(a2) on H^(2) (analytic in r)
      od%a2c = (0D0, 0D0)
      od%a2c(:, 0, 1) = CONJG(od%a2(:, 0, 1))
   END SUBROUTINE outer_setup

   SUBROUTINE stop_if_unordered(w1, w2, omega)
      ! ERROR STOP with fs_pair_order_message unless w1 > w2 and omega > 0
      REAL(8), INTENT(IN) :: w1, w2, omega
      CHARACTER(LEN=:), ALLOCATABLE :: msg

      msg = fs_pair_order_message(w1, w2, omega)
      IF (LEN(msg) > 0) THEN
         WRITE (*, '(a,a)') ' QTF free surface: ', msg
         ERROR STOP 1
      END IF
   END SUBROUTINE stop_if_unordered

   LOGICAL FUNCTION same_circle(a, b)
      TYPE(pw_set), INTENT(IN) :: a, b

      same_circle = ABS(a%radius - b%radius) <= 1D-9*a%radius .AND. ABS(a%h - b%h) <= 1D-9*a%h .AND. &
                    NORM2(a%centre - b%centre) <= 1D-9*a%radius
   END FUNCTION same_circle

   SUBROUTINE set_amplitudes(pw, a, kev, kr, d, e)
      ! a(-M:M, 0:N): the propagating amplitudes of pw as they are, the evanescent ones at R (a K_m(k_n R)); kev the
      ! evanescent wave numbers; kr(0:M+1, N) = e^{k_n R} K_m(k_n R); d, e(0:N) the vertical factors of phi_z and phi_zz.
      TYPE(pw_set), INTENT(IN) :: pw
      COMPLEX(8), INTENT(OUT) :: a(-pw%mmax:pw%mmax, 0:pw%nev)
      REAL(8), ALLOCATABLE, INTENT(OUT) :: kev(:), kr(:, :), d(:), e(:)
      INTEGER :: m, n

      ALLOCATE (kev(MAX(pw%nev, 1)), kr(0:pw%mmax + 1, MAX(pw%nev, 1)), d(0:pw%nev), e(0:pw%nev))
      kev = 0D0
      kr = 1D0
      a(:, 0) = pw%a(:, 0)
      d(0) = pw%k0*TANH(pw%k0*pw%h)
      e(0) = pw%k0**2
      DO n = 1, pw%nev
         kev(n) = pw%kev(n)
         CALL pw_bessel_k_scaled(pw%mmax + 1, pw%kev(n)*pw%radius, kr(:, n))
         DO m = -pw%mmax, pw%mmax
            a(m, n) = pw%a(m, n)*(kr(ABS(m), n)*EXP(-pw%kev(n)*pw%radius))
         END DO
         d(n) = -pw%kev(n)*TAN(pw%kev(n)*pw%h)
         e(n) = -pw%kev(n)**2
      END DO
   END SUBROUTINE set_amplitudes

   SUBROUTINE work_alloc(od, wk)
      ! The work arrays of one thread for ring_real, ring_tail and modes_real
      TYPE(outer_data), INTENT(IN) :: od
      TYPE(outer_work), INTENT(OUT) :: wk

      ALLOCATE (wk%h1(0:od%m1 + 1), wk%h2(0:od%m2 + 1), wk%hp(0:od%mp + 1), wk%j1(0:od%ni1 + 1), wk%j2(0:od%ni2 + 1))
      ALLOCATE (wk%kc(0:od%mx + 1))
      ALLOCATE (wk%kq1(0:od%m1, MAX(od%n1, 1)), wk%dkq1(0:od%m1, MAX(od%n1, 1)))
      ALLOCATE (wk%kq2(0:od%m2, MAX(od%n2, 1)), wk%dkq2(0:od%m2, MAX(od%n2, 1)))
      ALLOCATE (wk%kqp(0:od%mp, MAX(od%np, 1)), wk%dkqp(0:od%mp, MAX(od%np, 1)))
      ALLOCATE (wk%f1(-od%m1:od%m1, 4), wk%c2(-od%m2:od%m2, 4), wk%ps(-od%mp:od%mp, od%nk), wk%tmp(-od%mp:od%mp, 4))
      ALLOCATE (wk%fi(-od%ni1:od%ni1, 4), wk%ci(-od%ni2:od%ni2, 4))
      wk%kq1 = 0D0
      wk%dkq1 = 0D0
      wk%kq2 = 0D0
      wk%dkq2 = 0D0
      wk%kqp = 0D0
      wk%dkqp = 0D0
   END SUBROUTINE work_alloc

   SUBROUTINE radial_intervals(od, r_a, r_b, edge, ni)
      ! The interval edges of fs_outer_annulus: from r_a to r_b in steps of the smallest of 1/8 of 2 pi / (k_1 + k_2 + K),
      ! 2 / k_n of the fastest evanescent mode alive at the start of the step, and 2 r / (M + 1), M the most modes of a set
      TYPE(outer_data), INTENT(IN) :: od
      REAL(8), INTENT(IN) :: r_a, r_b
      REAL(8), ALLOCATABLE, INTENT(OUT) :: edge(:)
      INTEGER, INTENT(OUT) :: ni
      REAL(8) :: r, dl, kal, st
      INTEGER :: pass

      dl = 2D0*PI/(od%k1 + od%k2 + od%kp)/8D0
      DO pass = 1, 2
         ni = 0
         r = r_a
         IF (pass == 2) edge(1) = r_a
         DO WHILE (r < r_b)
            kal = MAX(fastest_alive(od%kev1, od%n1, r), fastest_alive(od%kev2, od%n2, r), fastest_alive(od%kevp, od%np, r))
            st = MIN(dl, 2D0*r/(od%mx + 1))
            IF (kal > 0D0) st = MIN(st, 2D0/kal)
            r = MIN(r + st, r_b)
            ni = ni + 1
            IF (pass == 2) edge(ni + 1) = r
         END DO
         IF (pass == 1) ALLOCATE (edge(ni + 1))
      END DO
      edge(ni + 1) = r_b
   CONTAINS
      REAL(8) FUNCTION fastest_alive(kev, nev, x)
         REAL(8), INTENT(IN) :: kev(:), x
         INTEGER, INTENT(IN) :: nev
         INTEGER :: n

         fastest_alive = 0D0
         DO n = 1, nev
            IF (kev(n)*(x - od%radius) < EV_DROP) fastest_alive = MAX(fastest_alive, kev(n))
         END DO
      END FUNCTION fastest_alive
   END SUBROUTINE radial_intervals

   SUBROUTINE radial_real(k0, kev, nev, mmax, radius, kr, r, hv, kq, dkq, nal, kc)
      ! The radial functions of one set at real r: hv(0:mmax+1) = H_m(k0 r) (J by the intrinsic, Y by its upward
      ! recurrence from Y_0 and Y_1, which is stable); for the evanescent modes still alive (k_n (r - R) < 40, the first nal)
      ! kq(m, n) = K_m(k_n r) / K_m(k_n R) and dkq(m, n) = its d/dr, from the scaled K of pw_bessel_k_scaled_c.
      INTEGER, INTENT(IN) :: nev, mmax
      REAL(8), INTENT(IN) :: k0, kev(:), radius, kr(0:, :), r
      COMPLEX(8), INTENT(OUT) :: hv(0:mmax + 1)
      REAL(8), INTENT(INOUT) :: kq(0:, :), dkq(0:, :)
      INTEGER, INTENT(OUT) :: nal
      COMPLEX(8), INTENT(INOUT) :: kc(0:mmax + 1)
      REAL(8) :: x, y0, y1, y2, ex, kn
      INTEGER :: m, n

      x = k0*r
      y0 = BESSEL_Y0(x)
      y1 = BESSEL_Y1(x)
      hv(0) = CMPLX(BESSEL_J0(x), y0, KIND=8)
      hv(1) = CMPLX(BESSEL_J1(x), y1, KIND=8)
      DO m = 1, mmax
         y2 = (2D0*m/x)*y1 - y0
         y0 = y1
         y1 = y2
         hv(m + 1) = CMPLX(BESSEL_JN(m + 1, x), y2, KIND=8)
      END DO
      nal = 0
      DO n = 1, nev
         kn = kev(n)
         IF (.NOT. kn*(r - radius) < EV_DROP) EXIT
         nal = n
         CALL pw_bessel_k_scaled_c(mmax + 1, CMPLX(kn*r, 0D0, KIND=8), kc)
         ex = EXP(-kn*(r - radius))
         kq(0, n) = REAL(kc(0))/kr(0, n)*ex
         dkq(0, n) = -kn*REAL(kc(1))/kr(0, n)*ex
         DO m = 1, mmax
            kq(m, n) = REAL(kc(m))/kr(m, n)*ex
            dkq(m, n) = -0.5D0*kn*(REAL(kc(m - 1)) + REAL(kc(m + 1)))/kr(m, n)*ex
         END DO
      END DO
   END SUBROUTINE radial_real

   PURE SUBROUTINE apply_modes(mmax, nev, nset, kset, a, k0, hv, kq, dkq, nal, d, e, f)
      ! The modal amplitudes at r of set kset of a(:, :, 1..nset): f(m, 1..4) = phi, d phi / dr, phi_z, phi_zz of mode m
      ! (factor e^{im theta}), from the propagating radial function hv (H_m(k0 r), or a scaled Hankel function on a tail
      ! contour) and the nal evanescent ratios kq, dkq of radial_real. H_{-m} = (-1)^m H_m; H_m' = (H_{m-1} - H_{m+1}) / 2,
      ! H_0' = -H_1.
      INTEGER, INTENT(IN) :: mmax, nev, nset, kset, nal
      COMPLEX(8), INTENT(IN) :: a(-mmax:mmax, 0:nev, nset), hv(0:mmax + 1)
      REAL(8), INTENT(IN) :: k0, kq(0:, :), dkq(0:, :), d(0:nev), e(0:nev)
      COMPLEX(8), INTENT(OUT) :: f(-mmax:mmax, 4)
      COMPLEX(8) :: p, hd
      REAL(8) :: sg
      INTEGER :: m, am, n

      DO m = -mmax, mmax
         am = ABS(m)
         sg = 1D0
         IF (m < 0 .AND. MODULO(am, 2) == 1) sg = -1D0
         IF (am == 0) THEN
            hd = -hv(1)
         ELSE
            hd = 0.5D0*(hv(am - 1) - hv(am + 1))
         END IF
         p = a(m, 0, kset)*(sg*hv(am))
         f(m, 1) = p
         f(m, 2) = a(m, 0, kset)*(sg*k0*hd)
         f(m, 3) = d(0)*p
         f(m, 4) = e(0)*p
         DO n = 1, nal
            p = a(m, n, kset)*kq(am, n)
            f(m, 1) = f(m, 1) + p
            f(m, 2) = f(m, 2) + a(m, n, kset)*dkq(am, n)
            f(m, 3) = f(m, 3) + d(n)*p
            f(m, 4) = f(m, 4) + e(n)*p
         END DO
      END DO
   END SUBROUTINE apply_modes

   PURE SUBROUTINE plane_modes(ni, amp, k, beta, s, jv, dz, f)
      ! The modal amplitudes of a plane wave, |n| <= ni: amp (s i)^n e^{-i s n beta} times J_n(kr), k J_n'(kr), dz J_n and
      ! k^2 J_n, with jv(0:ni+1) the radial function (J_n(kr), or half a scaled Hankel function on a tail contour).
      ! s = +1: the incident wave amp e^{ikr cos(theta - beta)} (Jacobi-Anger); s = -1: its conjugate continued in r,
      ! conj(amp) sum (-i)^n J_n(kr) e^{in beta} e^{-in theta}, amp given conjugated.
      INTEGER, INTENT(IN) :: ni
      COMPLEX(8), INTENT(IN) :: amp, jv(0:ni + 1)
      REAL(8), INTENT(IN) :: k, beta, s, dz
      COMPLEX(8), INTENT(OUT) :: f(-ni:ni, 4)
      COMPLEX(8) :: c, jd
      REAL(8) :: sg
      INTEGER :: n, an

      DO n = -ni, ni
         an = ABS(n)
         sg = 1D0
         IF (n < 0 .AND. MODULO(an, 2) == 1) sg = -1D0
         IF (an == 0) THEN
            jd = -jv(1)
         ELSE
            jd = 0.5D0*(jv(an - 1) - jv(an + 1))
         END IF
         c = amp*EXP(CI*(s*n*(0.5D0*PI - beta)))
         f(n, 1) = c*(sg*jv(an))
         f(n, 2) = c*(sg*k*jd)
         f(n, 3) = dz*f(n, 1)
         f(n, 4) = k**2*f(n, 1)
      END DO
   END SUBROUTINE plane_modes

   PURE SUBROUTINE ring_q(r, w1, m1, f1, w2, m2, c2, mp, nk, ps, g, out)
      ! out(k) += int_0^2pi psi_k Q[1, 2] dtheta = 2 pi sum_c ps(c, k) sum_a Q(a, a + c), for field 1 with amplitudes
      ! f1(a, :) on e^{ia theta}, conj(field 2) with c2(b, :) on e^{-ib theta}, psi_k with ps(c, k) on e^{ic theta}; Q(a, b)
      ! is free_surface_forcing (KERNEL_DIFF) of the two modes at r (complex on a tail contour): gradients (d_r, i a / r,
      ! d_z) and (d_r, -i b / r, d_z) of the conjugate, given to it conjugated since pair_product conjugates state 2.
      COMPLEX(8), INTENT(IN) :: r
      REAL(8), INTENT(IN) :: w1, w2, g
      INTEGER, INTENT(IN) :: m1, m2, mp, nk
      COMPLEX(8), INTENT(IN) :: f1(-m1:m1, 4), c2(-m2:m2, 4), ps(-mp:mp, nk)
      COMPLEX(8), INTENT(INOUT) :: out(nk)
      COMPLEX(8) :: tc, g1(3), g2(3), rinv
      INTEGER :: c, a, b, k

      rinv = 1D0/r
      DO c = -mp, mp
         tc = (0D0, 0D0)
         DO a = MAX(-m1, -m2 - c), MIN(m1, m2 - c)
            b = a + c
            g1(1) = f1(a, 2)
            g1(2) = CI*a*rinv*f1(a, 1)
            g1(3) = f1(a, 3)
            g2(1) = CONJG(c2(b, 2))
            g2(2) = CONJG(-CI*b*rinv*c2(b, 1))
            g2(3) = CONJG(c2(b, 3))
            tc = tc + free_surface_forcing(KERNEL_DIFF, w1, f1(a, 1), g1, f1(a, 4), w2, CONJG(c2(b, 1)), g2, &
                                           CONJG(c2(b, 4)), g)
         END DO
         DO k = 1, nk
            out(k) = out(k) + (2D0*PI)*ps(c, k)*tc
         END DO
      END DO
   END SUBROUTINE ring_q

   SUBROUTINE modes_real(od, r, wk)
      ! The modal amplitudes at real r of field 1 (wk%f1), conj(field 2) (wk%c2) and the psi_k (wk%ps), evanescent modes
      ! included
      TYPE(outer_data), INTENT(IN) :: od
      REAL(8), INTENT(IN) :: r
      TYPE(outer_work), INTENT(INOUT) :: wk
      INTEGER :: nal, k, c

      CALL radial_real(od%k1, od%kev1, od%n1, od%m1, od%radius, od%kr1, r, wk%h1, wk%kq1, wk%dkq1, nal, wk%kc)
      CALL apply_modes(od%m1, od%n1, 1, 1, od%a1, od%k1, wk%h1, wk%kq1, wk%dkq1, nal, od%d1, od%e1, wk%f1)
      CALL radial_real(od%k2, od%kev2, od%n2, od%m2, od%radius, od%kr2, r, wk%h2, wk%kq2, wk%dkq2, nal, wk%kc)
      CALL apply_modes(od%m2, od%n2, 1, 1, od%a2, od%k2, wk%h2, wk%kq2, wk%dkq2, nal, od%d2, od%e2, wk%c2)
      wk%c2 = CONJG(wk%c2)
      CALL radial_real(od%kp, od%kevp, od%np, od%mp, od%radius, od%krp, r, wk%hp, wk%kqp, wk%dkqp, nal, wk%kc)
      DO k = 1, od%nk
         CALL apply_modes(od%mp, od%np, od%nk, k, od%ap, od%kp, wk%hp, wk%kqp, wk%dkqp, nal, od%dp, od%ep, wk%tmp)
         DO c = -od%mp, od%mp
            wk%ps(c, k) = wk%tmp(c, 1)
         END DO
      END DO
   END SUBROUTINE modes_real

   SUBROUTINE ring_real(od, r, ring, wk)
      ! ring(k) = int_0^2pi psi_k Q_D dtheta at real r >= R: S_1 S_2, then I_1 S_2 and S_1 I_2 by their Jacobi-Anger modes
      TYPE(outer_data), INTENT(IN) :: od
      REAL(8), INTENT(IN) :: r
      COMPLEX(8), INTENT(OUT) :: ring(:)
      TYPE(outer_work), INTENT(INOUT) :: wk
      COMPLEX(8) :: rc
      INTEGER :: n

      CALL modes_real(od, r, wk)
      rc = CMPLX(r, 0D0, KIND=8)
      ring = (0D0, 0D0)
      CALL ring_q(rc, od%w1, od%m1, wk%f1, od%w2, od%m2, wk%c2, od%mp, od%nk, wk%ps, od%g, ring)
      IF (od%inc1) THEN
         DO n = 0, od%ni1 + 1
            wk%j1(n) = BESSEL_JN(n, od%k1*r)
         END DO
         CALL plane_modes(od%ni1, od%amp1, od%k1, od%beta1, 1D0, wk%j1, od%dzi1, wk%fi)
         CALL ring_q(rc, od%w1, od%ni1, wk%fi, od%w2, od%m2, wk%c2, od%mp, od%nk, wk%ps, od%g, ring)
      END IF
      IF (od%inc2) THEN
         DO n = 0, od%ni2 + 1
            wk%j2(n) = BESSEL_JN(n, od%k2*r)
         END DO
         CALL plane_modes(od%ni2, od%amp2c, od%k2, od%beta2, -1D0, wk%j2, od%dzi2, wk%ci)
         CALL ring_q(rc, od%w1, od%m1, wk%f1, od%w2, od%ni2, wk%ci, od%mp, od%nk, wk%ps, od%g, ring)
      END IF
   END SUBROUTINE ring_real

   SUBROUTINE ring_tail(od, cls, r, ring, wk)
      ! The slowly varying part G(r) of phase class cls (fs_outer_tail) at complex r: the ring integral of the class with
      ! every Hankel function scaled (pw_hankel_scaled), the propagating modes only
      TYPE(outer_data), INTENT(IN) :: od
      INTEGER, INTENT(IN) :: cls
      COMPLEX(8), INTENT(IN) :: r
      COMPLEX(8), INTENT(OUT) :: ring(:)
      TYPE(outer_work), INTENT(INOUT) :: wk
      INTEGER :: k, c, n

      ring = (0D0, 0D0)
      ! psi_k: H^(1)(K r)
      CALL pw_hankel_scaled(1, od%mp + 1, od%kp*r, wk%hp)
      DO k = 1, od%nk
         CALL apply_modes(od%mp, od%np, od%nk, k, od%ap, od%kp, wk%hp, wk%kqp, wk%dkqp, 0, od%dp, od%ep, wk%tmp)
         DO c = -od%mp, od%mp
            wk%ps(c, k) = wk%tmp(c, 1)
         END DO
      END DO
      ! conj(S_2): conj(a_2) H^(2)(k_2 r), in classes a and b
      IF (cls == 1 .OR. cls == 2) THEN
         CALL pw_hankel_scaled(2, od%m2 + 1, od%k2*r, wk%h2)
         CALL apply_modes(od%m2, od%n2, 1, 1, od%a2c, od%k2, wk%h2, wk%kq2, wk%dkq2, 0, od%d2, od%e2, wk%c2)
      END IF
      ! S_1: a_1 H^(1)(k_1 r), in classes a and c
      IF (cls == 1 .OR. cls == 3) THEN
         CALL pw_hankel_scaled(1, od%m1 + 1, od%k1*r, wk%h1)
         CALL apply_modes(od%m1, od%n1, 1, 1, od%a1, od%k1, wk%h1, wk%kq1, wk%dkq1, 0, od%d1, od%e1, wk%f1)
      END IF
      SELECT CASE (cls)
      CASE (1)
         CALL ring_q(r, od%w1, od%m1, wk%f1, od%w2, od%m2, wk%c2, od%mp, od%nk, wk%ps, od%g, ring)
         IF (od%inc1) THEN
            CALL pw_hankel_scaled(1, od%ni1 + 1, od%k1*r, wk%j1)
            DO n = 0, od%ni1 + 1
               wk%j1(n) = 0.5D0*wk%j1(n)
            END DO
            CALL plane_modes(od%ni1, od%amp1, od%k1, od%beta1, 1D0, wk%j1, od%dzi1, wk%fi)
            CALL ring_q(r, od%w1, od%ni1, wk%fi, od%w2, od%m2, wk%c2, od%mp, od%nk, wk%ps, od%g, ring)
         END IF
         IF (od%inc2) THEN
            CALL pw_hankel_scaled(2, od%ni2 + 1, od%k2*r, wk%j2)
            DO n = 0, od%ni2 + 1
               wk%j2(n) = 0.5D0*wk%j2(n)
            END DO
            CALL plane_modes(od%ni2, od%amp2c, od%k2, od%beta2, -1D0, wk%j2, od%dzi2, wk%ci)
            CALL ring_q(r, od%w1, od%m1, wk%f1, od%w2, od%ni2, wk%ci, od%mp, od%nk, wk%ps, od%g, ring)
         END IF
      CASE (2)
         CALL pw_hankel_scaled(2, od%ni1 + 1, od%k1*r, wk%j1)
         DO n = 0, od%ni1 + 1
            wk%j1(n) = 0.5D0*wk%j1(n)
         END DO
         CALL plane_modes(od%ni1, od%amp1, od%k1, od%beta1, 1D0, wk%j1, od%dzi1, wk%fi)
         CALL ring_q(r, od%w1, od%ni1, wk%fi, od%w2, od%m2, wk%c2, od%mp, od%nk, wk%ps, od%g, ring)
      CASE (3)
         CALL pw_hankel_scaled(1, od%ni2 + 1, od%k2*r, wk%j2)
         DO n = 0, od%ni2 + 1
            wk%j2(n) = 0.5D0*wk%j2(n)
         END DO
         CALL plane_modes(od%ni2, od%amp2c, od%k2, od%beta2, -1D0, wk%j2, od%dzi2, wk%ci)
         CALL ring_q(r, od%w1, od%m1, wk%f1, od%w2, od%ni2, wk%ci, od%mp, od%nk, wk%ps, od%g, ring)
      END SELECT
   END SUBROUTINE ring_tail

   FUNCTION modal_point(mmax, f, r, th, conj_field) RESULT(p)
      ! The point field at (r, theta) of modal amplitudes f (modes_real): phi, Cartesian gradient, phi_zz. conj_field: f
      ! holds the amplitudes of the conjugate on e^{-im theta} (wk%c2), and the field itself is returned.
      INTEGER, INTENT(IN) :: mmax
      COMPLEX(8), INTENT(IN) :: f(-mmax:mmax, 4)
      REAL(8), INTENT(IN) :: r, th
      LOGICAL, INTENT(IN) :: conj_field
      TYPE(fs_point_field) :: p
      COMPLEX(8) :: e, ph, ur, ut, uz, zz, x(4)
      INTEGER :: m, j

      ph = (0D0, 0D0)
      ur = (0D0, 0D0)
      ut = (0D0, 0D0)
      uz = (0D0, 0D0)
      zz = (0D0, 0D0)
      DO m = -mmax, mmax
         DO j = 1, 4
            x(j) = f(m, j)
            IF (conj_field) x(j) = CONJG(f(m, j))
         END DO
         e = EXP(CI*(m*th))
         ph = ph + x(1)*e
         ur = ur + x(2)*e
         ut = ut + CI*m*x(1)*e/r
         uz = uz + x(3)*e
         zz = zz + x(4)*e
      END DO
      p%phi = ph
      p%grad(1) = COS(th)*ur - SIN(th)*ut
      p%grad(2) = SIN(th)*ur + COS(th)*ut
      p%grad(3) = uz
      p%phizz = zz
   END FUNCTION modal_point

   FUNCTION plane_point(on, amp, k, beta, dz, dxy) RESULT(p)
      ! The incident wave amp e^{ik (dxy . e_beta)} at the offset dxy from the centre (amp: its value at the centre); zero if
      ! not on
      LOGICAL, INTENT(IN) :: on
      COMPLEX(8), INTENT(IN) :: amp
      REAL(8), INTENT(IN) :: k, beta, dz, dxy(2)
      TYPE(fs_point_field) :: p

      p%phi = (0D0, 0D0)
      p%grad = (0D0, 0D0)
      p%phizz = (0D0, 0D0)
      IF (.NOT. on) RETURN
      p%phi = amp*EXP(CI*k*(dxy(1)*COS(beta) + dxy(2)*SIN(beta)))
      p%grad(1) = CI*k*COS(beta)*p%phi
      p%grad(2) = CI*k*SIN(beta)*p%phi
      p%grad(3) = dz*p%phi
      p%phizz = k**2*p%phi
   END FUNCTION plane_point

   SUBROUTINE gauss_legendre01(n, x, w)
      ! Gauss-Legendre nodes and weights of order n on [0, 1] (Newton on P_n)
      INTEGER, INTENT(IN) :: n
      REAL(8), INTENT(OUT) :: x(n), w(n)
      REAL(8) :: t, pn, dpn, dt
      INTEGER :: i, it

      DO i = 1, n
         t = COS(PI*(i - 0.25D0)/(n + 0.5D0))
         DO it = 1, 100
            CALL legendre_p(n, t, pn, dpn)
            dt = pn/dpn
            t = t - dt
            IF (ABS(dt) < 1D-16) EXIT
         END DO
         CALL legendre_p(n, t, pn, dpn)
         x(i) = 0.5D0*(1D0 - t)
         w(i) = 1D0/((1D0 - t*t)*dpn*dpn)
      END DO
   END SUBROUTINE gauss_legendre01

   PURE SUBROUTINE legendre_p(n, t, pn, dpn)
      ! P_n(t) and P_n'(t) by the three-term recurrence
      INTEGER, INTENT(IN) :: n
      REAL(8), INTENT(IN) :: t
      REAL(8), INTENT(OUT) :: pn, dpn
      REAL(8) :: p0, p2
      INTEGER :: j

      p0 = 1D0
      pn = t
      DO j = 2, n
         p2 = ((2*j - 1)*t*pn - (j - 1)*p0)/j
         p0 = pn
         pn = p2
      END DO
      dpn = n*(t*pn - p0)/(t*t - 1D0)
   END SUBROUTINE legendre_p

   SUBROUTINE gauss_laguerre(n, x, w)
      ! Gauss-Laguerre nodes and weights of order n (weight e^{-x} on [0, inf)) from the eigenvalues of the Jacobi matrix
      ! (Golub and Welsch 1969): diagonal 2i - 1, off-diagonal i; w_i = (first component of eigenvector i)^2.
      INTEGER, INTENT(IN) :: n
      REAL(8), INTENT(OUT) :: x(n), w(n)
      REAL(8) :: e(MAX(n - 1, 1)), z(n, n), work(MAX(2*n - 2, 1))
      INTEGER :: i, info
      EXTERNAL :: DSTEQR

      DO i = 1, n
         x(i) = 2D0*i - 1D0
      END DO
      DO i = 1, n - 1
         e(i) = REAL(i, 8)
      END DO
      CALL DSTEQR('I', n, x, e, z, n, work, info)
      IF (info /= 0) THEN
         WRITE (*, '(a)') ' gauss_laguerre: the eigenvalues of the Jacobi matrix did not converge (DSTEQR).'
         ERROR STOP 1
      END IF
      DO i = 1, n
         w(i) = z(1, i)**2
      END DO
   END SUBROUTINE gauss_laguerre

END MODULE QTFFreeSurface
