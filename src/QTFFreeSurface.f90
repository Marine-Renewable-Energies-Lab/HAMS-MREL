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
   ! phi_zz = -lap_h phi by fields_green with four horizontal offsets (fs_inner_fields; offsets shrunk next to a waterline,
   ! fs_inner_steps); outside the circle r = R, the partial waves of the scattered potential (fs_circle_modes), from samples
   ! on the cylinder by the all-centroid route of the boundary integral formula (fs_sample_phi), so that the samples are a
   ! sum of the free-surface modes, with a check that the modes give the samples back (fs_mode_miss, fs_mode_warning).

   USE SecondOrderIncident, ONLY: free_surface_forcing
   USE NearFieldDrift, ONLY: KERNEL_DIFF
   USE SecondOrderTypes, ONLY: panel_mesh, first_order_state, cs_points, cs_fields, CS_FREE
   USE PartialWaves, ONLY: pw_set, pw_radial, pw_decompose, pw_eval_point, pw_radial_at
   USE FieldVelocity, ONLY: fields_green, incident_at
   USE WaterlineField, ONLY: field_influence_row, field_potential
   USE QTFField, ONLY: body_winding
   USE HAMS_mod, ONLY: NTHREAD
   USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: fs_point_field, fs_forcing_diff
   PUBLIC :: fs_inner_fields, fs_inner_steps, fs_sample_phi, fs_circle_modes, fs_mode_miss, fs_mode_warning, fs_auto_modes
   PUBLIC :: FS_MODE_TOL

   ! The mode check: the partial waves must give the samples on the circle to this fraction of their largest value
   REAL(8), PARAMETER :: FS_MODE_TOL = 1D-3

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
      ! fields_green (stock route) gives phi_P and grad phi_P at each point x and at its four horizontal offsets
      ! x +- delta e_x, x +- delta e_y. The horizontal Laplacian of phi_P is the central difference of the gradient there,
      !   lap_h phi_P = [d_x phi_P(x + delta e_x) - d_x phi_P(x - delta e_x) + d_y phi_P(x + delta e_y) - d_y phi_P(x - delta e_y)]
      !                 / (2 delta),
      ! and phi_P,zz = -lap_h phi_P (Laplace). The offsets difference over delta again (fields_green), so the stencil reaches
      ! x +- 2 delta and the diagonals x +- delta e_x +- delta e_y. All points of a stencil take the near-or-far and series
      ! decisions of x (decide of fields_green), so that no switch between them enters the differences. The incident part is
      ! added exactly: phi_I, grad phi_I and phi_I,zz = k^2 phi_I. At z = 0, d phi_P/dz = nu phi_P (fields_green).
      ! delta(pts%n), if given, is the offset of each point; otherwise fs_inner_steps gives it: step x the point size, smaller
      ! next to a waterline. All states must be at the current global frequency, as for fields_green.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(cs_points), INTENT(IN) :: pts
      TYPE(first_order_state), INTENT(IN) :: states(:)
      REAL(8), INTENT(IN) :: step
      LOGICAL, INTENT(IN) :: is_multi
      TYPE(fs_point_field), INTENT(OUT) :: f(pts%n, SIZE(states))
      REAL(8), INTENT(IN), OPTIONAL :: delta(:)
      ! The five points of a stencil: the centre, then +x, -x, +y, -y
      REAL(8), PARAMETER :: OFF(2, 5) = RESHAPE([0D0, 0D0, 1D0, 0D0, -1D0, 0D0, 0D0, 1D0, 0D0, -1D0], [2, 5])
      TYPE(cs_points) :: ex
      TYPE(cs_fields), ALLOCATABLE :: fl(:)
      TYPE(first_order_state), ALLOCATABLE :: ps(:)
      REAL(8), ALLOCATABLE :: d(:), dec(:, :)
      COMPLEX(8) :: phi_i, grad_i(3), lap
      INTEGER :: i, j, j0, ist, nshrunk

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

      ex%n = 5*pts%n
      ALLOCATE (ex%x(3, ex%n), ex%nrm(3, ex%n), ex%w(ex%n), ex%psize(ex%n), ex%kind(ex%n), ex%body(ex%n), dec(3, ex%n))
      DO i = 1, pts%n
         DO j = 1, 5
            j0 = 5*(i - 1) + j
            ex%x(:, j0) = pts%x(:, i) + [d(i)*OFF(1, j), d(i)*OFF(2, j), 0D0]
            ex%psize(j0) = d(i)/step
            dec(:, j0) = pts%x(:, i)
         END DO
      END DO
      ex%nrm = 0D0
      ex%nrm(3, :) = -1D0
      ex%w = 0D0
      ex%kind = CS_FREE
      ex%body = 0

      ! The perturbation part only: the incident wave is added exactly below
      ps = states
      DO ist = 1, SIZE(ps)
         ps(ist)%amp = 0D0
      END DO
      ALLOCATE (fl(SIZE(states)))
      CALL fields_green(mesh, ex, ps, step, is_multi, fl, decide=dec)

      DO ist = 1, SIZE(states)
         DO i = 1, pts%n
            j0 = 5*(i - 1)
            lap = (fl(ist)%grad(1, j0 + 2) - fl(ist)%grad(1, j0 + 3) + fl(ist)%grad(2, j0 + 4) - fl(ist)%grad(2, j0 + 5)) &
                  /(2D0*d(i))
            CALL incident_at(states(ist), pts%x(:, i), phi_i, grad_i)
            f(i, ist)%phi = fl(ist)%phi(j0 + 1) + phi_i
            f(i, ist)%grad = fl(ist)%grad(:, j0 + 1) + grad_i
            f(i, ist)%phizz = -lap + states(ist)%k**2*phi_i
         END DO
      END DO
   END SUBROUTINE fs_inner_fields

   SUBROUTINE fs_inner_steps(mesh, pts, step, delta, nshrunk)
      ! The offset delta of each inner point for fs_inner_fields: step x its size, halved (at most MAX_HALVINGS times) until
      ! no point of its stencil lies in a waterplane, where fields_green gives no fluid value, and no point of the stencil
      ! scaled by CLEARANCE either. The stencil points are x +- delta and x +- 2 delta along x and y, and
      ! x +- delta e_x +- delta e_y; a point lies in a waterplane when body_winding (QTFField) is 1/2 or more for a body.
      ! The scaled stencil keeps the waterline about CLEARANCE delta away (to the 22.5 degrees between its 8 directions):
      ! next to the hull the field varies on the scale of the distance to it, so a stencil that only just clears the
      ! waterplane differences it badly (fs_waterline_offsets: 23 % in phi_zz at 2 delta = 0.83 of the distance).
      ! nshrunk counts the points whose delta was shrunk. The points next to a waterline so keep central differences, with
      ! a smaller step. Frequency-independent: the driver calls it once.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(cs_points), INTENT(IN) :: pts
      REAL(8), INTENT(IN) :: step
      REAL(8), INTENT(OUT) :: delta(pts%n)
      INTEGER, INTENT(OUT) :: nshrunk
      INTEGER, PARAMETER :: MAX_HALVINGS = 30
      REAL(8), PARAMETER :: CLEARANCE = 16D0
      REAL(8) :: d
      INTEGER :: i, it

      nshrunk = 0
      !$OMP PARALLEL DO NUM_THREADS(NTHREAD) PRIVATE(d, it) REDUCTION(+:nshrunk) SCHEDULE(DYNAMIC, 4)
      DO i = 1, pts%n
         d = step*pts%psize(i)
         DO it = 1, MAX_HALVINGS
            IF (.NOT. (stencil_in_waterplane(mesh, pts%x(:, i), d) .OR. &
                       stencil_in_waterplane(mesh, pts%x(:, i), CLEARANCE*d))) EXIT
            d = 0.5D0*d
         END DO
         IF (d < step*pts%psize(i)) nshrunk = nshrunk + 1
         delta(i) = d
      END DO
      !$OMP END PARALLEL DO
   END SUBROUTINE fs_inner_steps

   LOGICAL FUNCTION stencil_in_waterplane(mesh, x, d)
      ! A point of the stencil of fs_inner_fields at x with offset d lies in a waterplane (body_winding >= 1/2)
      TYPE(panel_mesh), INTENT(IN) :: mesh
      REAL(8), INTENT(IN) :: x(3), d
      REAL(8), PARAMETER :: OFF(2, 12) = RESHAPE([1D0, 0D0, -1D0, 0D0, 0D0, 1D0, 0D0, -1D0, 2D0, 0D0, -2D0, 0D0, &
                                                  0D0, 2D0, 0D0, -2D0, 1D0, 1D0, 1D0, -1D0, -1D0, 1D0, -1D0, -1D0], [2, 12])
      REAL(8), ALLOCATABLE :: wind(:)
      INTEGER :: j

      stencil_in_waterplane = .FALSE.
      DO j = 1, 12
         CALL body_winding(mesh, [x(1) + d*OFF(1, j), x(2) + d*OFF(2, j), x(3)], wind)
         IF (MAXVAL(wind) >= 0.5D0) THEN
            stencil_in_waterplane = .TRUE.
            RETURN
         END IF
      END DO
   END FUNCTION stencil_in_waterplane

   SUBROUTINE fs_sample_phi(mesh, pts, states, is_multi, phi, centroid_only)
      ! The perturbation potential phi_P of every state at every point of pts, any kind and any z <= 0, from one row of the
      ! boundary integral formula per point (field_influence_row in field mode, with the point as its own decision point, as
      ! fields_green takes phi at a centre); no gradient. phi(pts%n, SIZE(states)). With centroid_only every panel takes
      ! its centroid value (near_panel = .FALSE.), so that phi_P is a sum of the free-surface modes outside the hull: the
      ! samples of the partial waves on the circle. Otherwise the stock route of fields_green. All states must be at the
      ! current global frequency.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(cs_points), INTENT(IN) :: pts
      TYPE(first_order_state), INTENT(IN) :: states(:)
      LOGICAL, INTENT(IN) :: is_multi
      COMPLEX(8), INTENT(OUT) :: phi(pts%n, SIZE(states))
      LOGICAL, INTENT(IN), OPTIONAL :: centroid_only
      COMPLEX(8), ALLOCATABLE :: row_g(:, :), row_dg(:, :)
      LOGICAL :: near
      INTEGER :: i, ist

      near = .TRUE.
      IF (PRESENT(centroid_only)) near = .NOT. centroid_only
      ! As fields_green: the points are independent, the rows are per thread, the panel loop runs on one thread here
      !$OMP PARALLEL NUM_THREADS(NTHREAD) PRIVATE(row_g, row_dg, ist)
      ALLOCATE (row_g(mesh%nbase, mesh%nsys), row_dg(mesh%nbase, mesh%nsys))
      !$OMP DO SCHEDULE(DYNAMIC, 4)
      DO i = 1, pts%n
         CALL field_influence_row(mesh, pts%x(:, i), is_multi, row_g, row_dg, pts%x(:, i), field_mode=.TRUE., near_panel=near)
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
      ! (fs_mode_miss): how far the partial waves are from the samples, relative to the largest sample of all the states,
      ! so that a state that is zero by symmetry (psi of yaw on an axisymmetric body) does not warn on its rounding. The
      ! samples are evaluated in the fluid whatever their kind (CS_WALL from pw_sample_points).
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(cs_points), INTENT(IN) :: circle
      TYPE(first_order_state), INTENT(IN) :: states(:)
      LOGICAL, INTENT(IN) :: is_multi
      REAL(8), INTENT(IN) :: centre(2), radius, h, k0, kev(:)
      INTEGER, INTENT(IN) :: mmax, nev
      TYPE(pw_set), INTENT(OUT) :: pw(SIZE(states))
      REAL(8), INTENT(OUT) :: miss(SIZE(states))
      COMPLEX(8), ALLOCATABLE :: samples(:, :)
      REAL(8) :: big
      INTEGER :: ist

      ALLOCATE (samples(circle%n, SIZE(states)))
      CALL fs_sample_phi(mesh, circle, states, is_multi, samples, centroid_only=.TRUE.)
      big = 0D0
      IF (SIZE(samples) > 0) big = MAXVAL(ABS(samples))
      DO ist = 1, SIZE(states)
         CALL pw_decompose(samples(:, ist), centre, radius, h, k0, kev, mmax, nev, pw(ist))
         miss(ist) = fs_mode_miss(pw(ist), circle, samples(:, ist), big)
      END DO
   END SUBROUTINE fs_circle_modes

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

END MODULE QTFFreeSurface
