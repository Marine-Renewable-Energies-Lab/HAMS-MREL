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
   USE PartialWaves, ONLY: pw_set, pw_radial, pw_decompose, pw_eval_point, pw_radial_at
   USE FieldVelocity, ONLY: incident_at
   USE WaterlineField, ONLY: field_influence_row, field_potential
   USE QTFField, ONLY: body_winding
   USE HAMS_mod, ONLY: NTHREAD
   USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: fs_point_field, fs_forcing_diff
   PUBLIC :: fs_inner_fields, fs_inner_steps, fs_sample_phi, fs_circle_modes, fs_mode_miss, fs_mode_warning, fs_auto_modes
   PUBLIC :: FS_MODE_TOL, fs_mode_scales, fs_auto_evanescent

   ! The mode check: the partial waves must give the samples on the circle to this fraction of their largest value
   REAL(8), PARAMETER :: FS_MODE_TOL = 1D-3
   ! The smallest scale of the mode check, as a fraction of the largest sample over all states (fs_mode_scales)
   REAL(8), PARAMETER :: FS_SCALE_FLOOR = 1D-6

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
      ! and the first one left out decays faster still (ruling of the Task 6 fix round).
      INTEGER, INTENT(IN) :: nev
      REAL(8), INTENT(IN) :: h, radius, r_wl
      REAL(8), PARAMETER :: PI = 3.14159265358979323846D0

      IF (nev > 0) THEN
         fs_auto_evanescent = nev
      ELSE
         fs_auto_evanescent = MAX(8, CEILING(9D0*h/(PI*(radius - r_wl))))
      END IF
   END FUNCTION fs_auto_evanescent

END MODULE QTFFreeSurface
