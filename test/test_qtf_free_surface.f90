!> Unit tests of QTFFreeSurface and of the free-surface forcing kernel of SecondOrderIncident (Phase 6a, Tasks 1 and 6).
!> The field tests of Task 6 set module globals (the mesh arrays and the wave of WaveDyn_mod), so the suite runs serially
!> (see unit_tester.f90).
module test_qtf_free_surface
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderIncident, only : incident2, incident2_pair, incident2_phi, incident2_grad, free_surface_forcing
  use FirstOrderState, only : incident_phi, incident_grad, wave_gravity
  use NearFieldDrift, only : KERNEL_DIFF
  use QTFFreeSurface, only : fs_point_field, fs_forcing_diff, fs_inner_fields, fs_inner_steps, fs_circle_modes, &
                             fs_mode_warning, fs_auto_modes, FS_MODE_TOL, fs_sample_phi, fs_mode_scales, fs_auto_evanescent
  use HAMS_mod, only : NTHREAD
  use WaveDyn_mod, only : BETA, W1, WK
  use SecondOrderTypes, only : so_options, panel_mesh, first_order_state, cs_points, cs_fields, CS_FREE
  use FieldVelocity, only : fields_green
  use QTFField, only : body_winding
  use PartialWaves, only : pw_set, pw_sample_points, pw_evanescent_from_wvn
  use second_order_fixtures, only : make_cylinder, maccamy_fuchs_state
  use test_waterline_field, only : set_wave, load_panels
  implicit none
  private

  public :: collect_qtf_free_surface

  real(8), parameter :: g = 9.81d0, pi = 3.14159265358979323846d0, deg = pi / 180d0
  ! Depth of the cylinder fixture of the field tests (as pw_hams_field: 5 m of water under the keel)
  real(8), parameter :: fx_depth = 15d0

contains

subroutine collect_qtf_free_surface(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("fs_forcing_matches_incident2", test_matches_incident2), &
    new_unittest("fs_forcing_swap", test_swap), &
    new_unittest("fs_forcing_incident_only_zero", test_incident_only_zero), &
    new_unittest("fs_laplace_residual", test_laplace_residual), &
    new_unittest("fs_waterline_offsets", test_waterline_offsets), &
    new_unittest("fs_mode_check_warns", test_mode_check_warns), &
    new_unittest("fs_mode_check_per_state", test_mode_check_per_state), &
    new_unittest("fs_auto_evanescent", test_auto_evanescent) &
    ]
end subroutine collect_qtf_free_surface

!> Wave number of omega in depth h, from omega^2 = g k tanh(k h).
real(8) function wave_number(omega, h)
  real(8), intent(in) :: omega, h
  real(8) :: k, nu
  integer :: it
  nu = omega**2 / g
  k = max(nu, sqrt(nu / h))
  do it = 1, 100
    k = k - (k * tanh(k * h) - nu) / (tanh(k * h) + k * h / cosh(k * h)**2)
  end do
  wave_number = k
end function wave_number

!> Incident wave (omega, heading in degrees) as a point field at x.
function wave_field(omega, beta, h, x) result(f)
  real(8), intent(in) :: omega, beta, h, x(3)
  type(fs_point_field) :: f
  real(8) :: k
  k = wave_number(omega, h)
  f%phi = incident_phi(x, omega, k, h, 1d0, beta * deg, g)
  call incident_grad(x, omega, k, h, 1d0, beta * deg, g, f%grad)
  f%phizz = k**2 * f%phi
end function wave_field

!> A reproducible complex field from a seed.
function pseudo_field(seed) result(f)
  integer, intent(in) :: seed
  type(fs_point_field) :: f
  integer :: m
  f%phi = cmplx(sin(1.3d0 * seed), cos(0.7d0 * seed + 1d0), 8)
  do m = 1, 3
    f%grad(m) = cmplx(cos(0.9d0 * seed + m), sin(1.7d0 * seed * m + 0.3d0), 8)
  end do
  f%phizz = cmplx(sin(2.1d0 * seed + 0.5d0), cos(0.4d0 * seed), 8)
end function pseudo_field

subroutine test_matches_incident2(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: depth = 100d0, w1 = 0.6d0, w2 = 0.8d0, b1 = 0d0, b2 = 30d0
  real(8), parameter :: pts(3, 3) = reshape([0d0, 0d0, 0d0, 4.2d0, -1.7d0, 0d0, -9d0, 6.5d0, 0d0], [3, 3])
  type(incident2) :: inc
  type(fs_point_field) :: f1, f2
  complex(8) :: q, expected, gr(3)
  real(8) :: k1, k2
  integer :: ip, ia

  k1 = wave_number(w1, depth)
  k2 = wave_number(w2, depth)
  ! both orders of the pair
  do ia = 1, 2
    if (ia == 1) then
      call incident2_pair(w1, k1, b1 * deg, w2, k2, b2 * deg, depth, g, inc)
    else
      call incident2_pair(w2, k2, b2 * deg, w1, k1, b1 * deg, depth, g, inc)
    end if
    do ip = 1, 3
      if (ia == 1) then
        f1 = wave_field(w1, b1, depth, pts(:, ip))
        f2 = wave_field(w2, b2, depth, pts(:, ip))
        q = free_surface_forcing(KERNEL_DIFF, w1, f1%phi, f1%grad, f1%phizz, w2, f2%phi, f2%grad, f2%phizz, g)
      else
        f1 = wave_field(w2, b2, depth, pts(:, ip))
        f2 = wave_field(w1, b1, depth, pts(:, ip))
        q = free_surface_forcing(KERNEL_DIFF, w2, f1%phi, f1%grad, f1%phizz, w1, f2%phi, f2%grad, f2%phizz, g)
      end if
      call incident2_grad(inc, pts(:, ip), gr)
      expected = -inc%omega**2 * incident2_phi(inc, pts(:, ip)) + g * gr(3)
      call check(error, abs(q - expected) <= 1d-10 * abs(expected), &
                 "free_surface_forcing differs from -Omega^2 phi_I2 + g d phi_I2/dz")
      if (allocated(error)) return
    end do
  end do
end subroutine test_matches_incident2

subroutine test_swap(error)
  type(error_type), allocatable, intent(out) :: error
  type(fs_point_field) :: t1, i1, t2, i2
  complex(8) :: a, b
  integer :: is

  do is = 1, 4
    t1 = pseudo_field(10 * is + 1)
    i1 = pseudo_field(10 * is + 2)
    t2 = pseudo_field(10 * is + 3)
    i2 = pseudo_field(10 * is + 4)
    a = fs_forcing_diff(0.9d0, t1, i1, 0.55d0, t2, i2, g)
    b = fs_forcing_diff(0.55d0, t2, i2, 0.9d0, t1, i1, g)
    call check(error, abs(a - conjg(b)) <= 1d-14 * max(abs(a), 1d-300), "fs_forcing_diff(i, j) is not conj(fs_forcing_diff(j, i))")
    if (allocated(error)) return
  end do
end subroutine test_swap

subroutine test_incident_only_zero(error)
  type(error_type), allocatable, intent(out) :: error
  type(fs_point_field) :: f1, f2

  f1 = pseudo_field(5)
  f2 = pseudo_field(8)
  call check(error, fs_forcing_diff(0.9d0, f1, f1, 0.55d0, f2, f2, g) == (0d0, 0d0), "forcing with total = incident is not exactly 0")
end subroutine test_incident_only_zero

! ---------------------------------------------------------------------------------------------- fields (Task 6)

!> The 448-panel sharp cylinder of the QTF fixtures (radius 10, draft 10, sharp bottom edge) in depth 15, with the
!> MacCamy-Fuchs scattered wave on its panels and amp = 0, so that only the panel integrals enter (as in pw_hams_field).
!> Sets the module globals of the wave (set_wave): omega 1.1, heading 0.6.
subroutine cylinder_fixture(mesh, states, error)
  type(panel_mesh), intent(out) :: mesh
  type(first_order_state), intent(out) :: states(1)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: full

  call set_wave(fx_depth)
  call make_cylinder(10d0, 10d0, 32, 8, 6, full, .true.)
  call load_panels(full, 0, 0, mesh)
  call check(error, mesh%n == 448, "the fixture is not the 448-panel cylinder")
  if (allocated(error)) return
  call maccamy_fuchs_state(mesh, 10d0, fx_depth, WK, W1, 1d0, BETA, wave_gravity(W1, WK, fx_depth), 1000d0, states(1))
  states(1)%amp = 0d0
end subroutine cylinder_fixture

!> Points of kind CS_FREE at x(3, n), all of size psz.
subroutine free_points(x, psz, pts)
  real(8), intent(in) :: x(:, :), psz
  type(cs_points), intent(out) :: pts

  pts%n = size(x, 2)
  allocate (pts%x(3, pts%n), pts%nrm(3, pts%n), pts%w(pts%n), pts%psize(pts%n), pts%kind(pts%n), pts%body(pts%n))
  pts%x = x
  pts%nrm = 0d0
  pts%nrm(3, :) = -1d0
  pts%w = 1d0
  pts%psize = psz
  pts%kind = CS_FREE
  pts%body = 0
end subroutine free_points

!> At points below z = 0 (-2 m, at r = 12, 14.5 and 19 m from the cylinder axis), the horizontal Laplacian of
!> fs_inner_fields (5-point, phi_zz = -lap_h phi) plus d2 phi / dz2 from the 3-point difference of phi at z -+ delta
!> (fs_sample_phi with the same decision point) is below 1e-3 of |lap_h phi|: Laplace's equation, so the offsets
!> difference the right field. delta = step x point size = 0.1 x 1 m, the driver's 2 x cs_fd_step for a 1 m point.
!> Both differences are O(delta^2) in truncation, but the rounding of the Green function enters divided by delta^2, so
!> there is a best delta, not a limit delta -> 0. Measured in fix round 1, residual relative at r = 12, 14.5, 19 m:
!>   delta 0.02: 4.2e-3 2.0e-3 1.6e-2    delta 0.04: 2.1e-3 1.5e-3 2.4e-3    delta 0.05: 1.4e-3 6.4e-4 1.8e-3
!>   delta 0.08: 1.0e-3 7.2e-4 5.6e-4    delta 0.10: 6.6e-4 1.8e-4 7.2e-4    delta 0.20: 1.4e-3 2.4e-3 1.5e-3
!>   delta 0.30: 2.4e-3 2.2e-3 4.6e-4
!> At 0.1 these are the values of the 25-row stencil of 9fdd193 at 0.05 (the same differences, to rounding). The gate is
!> on delta = 0.1; delta = 0.05 is printed only.
subroutine test_laplace_residual(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: zp = -2d0, psz = 1d0, steps(2) = [0.1d0, 0.05d0]
  real(8), parameter :: rp(3) = [12d0, 14.5d0, 19d0], thp(3) = [0.3d0, 2.0d0, -1.9d0]
  type(panel_mesh) :: mesh
  type(first_order_state) :: states(1)
  type(cs_points) :: pts, vp
  type(fs_point_field) :: f(3, 1)
  complex(8) :: ph(9, 1)
  real(8) :: x(3, 3), xv(3, 9), dec(3, 9), d, worst(2), lap_min
  complex(8) :: lap_h, phizz_fd
  integer :: i, is, nthread_saved
  character(len=200) :: info

  call cylinder_fixture(mesh, states, error)
  if (allocated(error)) return
  nthread_saved = NTHREAD
  NTHREAD = 4
  worst = 0d0
  lap_min = huge(1d0)
  do is = 1, 2
    d = steps(is) * psz
    do i = 1, 3
      x(:, i) = [rp(i) * cos(thp(i)), rp(i) * sin(thp(i)), zp]
      xv(:, 3 * i - 2) = x(:, i)
      xv(:, 3 * i - 1) = x(:, i) + [0d0, 0d0, d]
      xv(:, 3 * i) = x(:, i) - [0d0, 0d0, d]
      dec(:, 3 * i - 2:3 * i) = spread(x(:, i), 2, 3)
    end do
    call free_points(x, psz, pts)
    call free_points(xv, psz, vp)
    call fs_inner_fields(mesh, pts, states, steps(is), .false., f)
    call fs_sample_phi(mesh, vp, states, .false., ph, decide=dec)
    do i = 1, 3
      lap_h = -f(i, 1)%phizz
      phizz_fd = (ph(3 * i - 1, 1) + ph(3 * i, 1) - 2d0 * ph(3 * i - 2, 1)) / d**2
      if (is == 1) lap_min = min(lap_min, abs(lap_h))
      worst(is) = max(worst(is), abs(lap_h + phizz_fd) / max(abs(lap_h), tiny(1d0)))
      print '(4x,a,f5.2,a,f6.2,a,2es12.4,a,es10.2)', "fs_laplace_residual delta", d, " at r =", rp(i), &
        ": |lap_h phi|, |phi_zz| ", abs(lap_h), abs(phizz_fd), ", residual relative ", &
        abs(lap_h + phizz_fd) / max(abs(lap_h), tiny(1d0))
    end do
  end do
  NTHREAD = nthread_saved
  call check(error, lap_min > 0d0, "the horizontal Laplacian of fs_inner_fields is zero")
  if (allocated(error)) return
  write (info, '(a,2es10.3)') "lap_h phi + phi_zz relative to |lap_h phi| at delta = 0.1, 0.05 m ", worst
  call check(error, worst(1) <= 1d-3, trim(info))
end subroutine test_laplace_residual

!> True if a point of the stencil of fs_inner_fields at x with offset d (x +- d e_x, x +- d e_y) lies in a waterplane.
!> An independent copy, on purpose, of the private test of QTFFreeSurface: the test checks the module's offsets with its
!> own oracle, so the two must be changed together if the stencil changes.
logical function stencil_in_waterplane(mesh, x, d)
  type(panel_mesh), intent(in) :: mesh
  real(8), intent(in) :: x(3), d
  real(8), allocatable :: wind(:)
  real(8) :: off(2, 4)
  integer :: j

  off = reshape([d, 0d0, -d, 0d0, 0d0, d, 0d0, -d], [2, 4])
  stencil_in_waterplane = .false.
  do j = 1, 4
    call body_winding(mesh, [x(1) + off(1, j), x(2) + off(2, j), 0d0], wind)
    if (maxval(wind) >= 0.5d0) stencil_in_waterplane = .true.
  end do
end function stencil_in_waterplane

!> An inner point about 0.3 m from the waterline of the cylinder (r = 10.3 m) with a point size of 10 m and step 0.05:
!> the offsets of 0.5 m reach into the waterplane. fs_inner_steps shrinks delta there, and only there (a point at
!> r = 25 m keeps 0.5 m), until no stencil point is in the waterplane and the waterline is about 16 delta away.
!> fs_inner_fields with it agrees with a reference at delta = 0.0025 m (point size 0.05 m) to 1e-2 in phi_zz and 1e-3 in
!> phi and grad phi. Two stencils that do not: 0.5 m, inside the waterplane (phi_zz off by more than 10 %), and 0.25 m,
!> the largest that only just clears it (off by more than 1e-2: the reason for the clearance of 16 offsets).
!> Measured (fix round 1): delta 0.0156 m after shrinking, phi_zz 1.4e-3 off; just clear (0.25 m) 2.3e-1; unshrunk 59.
subroutine test_waterline_offsets(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: step = 0.05d0, psz = 10d0, th0 = 0.05d0
  type(panel_mesh) :: mesh
  type(first_order_state) :: states(1)
  type(cs_points) :: pts, one, ref_pt
  type(fs_point_field) :: f(2, 1), fref(1, 1), fbad(1, 1), fclear(1, 1)
  real(8) :: x(3, 2), delta(2), e_zz, e_phi, e_grad, e_bad, e_clear
  integer :: nshrunk, nfloor, nthread_saved
  character(len=200) :: info

  call cylinder_fixture(mesh, states, error)
  if (allocated(error)) return
  x(:, 1) = [10.3d0 * cos(th0), 10.3d0 * sin(th0), 0d0]
  x(:, 2) = [25d0 * cos(th0 + 1d0), 25d0 * sin(th0 + 1d0), 0d0]
  call free_points(x, psz, pts)
  call check(error, stencil_in_waterplane(mesh, x(:, 1), step * psz), &
             "the test point's stencil at 0.5 m does not reach the waterplane: the test does not test")
  if (allocated(error)) return
  call check(error, .not. stencil_in_waterplane(mesh, x(:, 1), 0.25d0), &
             "the stencil at 0.25 m reaches the waterplane: the just-clear comparison does not test")
  if (allocated(error)) return

  call fs_inner_steps(mesh, pts, step, delta, nshrunk, nfloor)
  print '(4x,a,2f10.5,a,i0,a,i0)', "fs_inner_steps next to the waterline and away from it: delta ", delta, ", shrunk ", &
    nshrunk, ", at the floor ", nfloor
  call check(error, nshrunk == 1 .and. nfloor == 0, "points with a shrunk delta, or at the floor")
  if (allocated(error)) return
  call check(error, delta(2) == step * psz, "delta away from the waterline changed")
  if (allocated(error)) return
  call check(error, delta(1) < step * psz .and. delta(1) > 0d0, "delta next to the waterline not shrunk")
  if (allocated(error)) return
  call check(error, .not. stencil_in_waterplane(mesh, x(:, 1), delta(1)), "the shrunk stencil still reaches the waterplane")
  if (allocated(error)) return

  nthread_saved = NTHREAD
  NTHREAD = 4
  call fs_inner_fields(mesh, pts, states, step, .false., f)
  call free_points(x(:, 1:1), 0.05d0, ref_pt)
  call fs_inner_fields(mesh, ref_pt, states, step, .false., fref)
  call free_points(x(:, 1:1), psz, one)
  call fs_inner_fields(mesh, one, states, step, .false., fbad, delta=[step * psz])
  call fs_inner_fields(mesh, one, states, step, .false., fclear, delta=[0.25d0])
  NTHREAD = nthread_saved

  e_zz = abs(f(1, 1)%phizz - fref(1, 1)%phizz) / abs(fref(1, 1)%phizz)
  e_phi = abs(f(1, 1)%phi - fref(1, 1)%phi) / abs(fref(1, 1)%phi)
  e_grad = norm2(abs(f(1, 1)%grad - fref(1, 1)%grad)) / norm2(abs(fref(1, 1)%grad))
  e_bad = abs(fbad(1, 1)%phizz - fref(1, 1)%phizz) / abs(fref(1, 1)%phizz)
  e_clear = abs(fclear(1, 1)%phizz - fref(1, 1)%phizz) / abs(fref(1, 1)%phizz)
  print '(4x,a,3es10.2,a,2es10.2)', "next to the waterline against delta = 0.0025 m: phi_zz, phi, grad ", e_zz, e_phi, &
    e_grad, "; phi_zz unshrunk and just clear ", e_bad, e_clear
  write (info, '(a,3es10.3)') "shrunk stencil against the reference (phi_zz, phi, grad) ", e_zz, e_phi, e_grad
  call check(error, e_zz <= 1d-2 .and. e_phi <= 1d-3 .and. e_grad <= 1d-3, trim(info))
  if (allocated(error)) return
  write (info, '(a,es10.3)') "the stencil inside the waterplane agrees with the reference: the test does not test ", e_bad
  call check(error, e_bad > 0.1d0, trim(info))
  if (allocated(error)) return
  write (info, '(a,es10.3)') "the just-clear stencil agrees with the reference: the clearance is not needed ", e_clear
  call check(error, e_clear > 1d-2, trim(info))
end subroutine test_waterline_offsets

!> With qtf_fs_modes = 2 on the fixture (circle R = 15 m about the axis, nev = 8), the partial waves cannot give the
!> samples on the circle (k0 R = 2.3 needs |m| up to about 13), and fs_mode_warning gives the warning of the brief. With
!> the automatic M = ceil(k0 R) + 10 they do, to 1e-3, and there is no warning.
subroutine test_mode_check_warns(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: radius = 15d0
  integer, parameter :: nev = 8
  type(so_options) :: opts
  type(panel_mesh) :: mesh
  type(first_order_state) :: states(1)
  type(cs_points) :: circle
  type(pw_set) :: pw(1)
  real(8) :: kev(nev), miss(1), miss_auto
  integer :: mmax, nthread_saved
  character(len=:), allocatable :: msg
  character(len=200) :: info

  call cylinder_fixture(mesh, states, error)
  if (allocated(error)) return
  call check(error, fs_auto_modes(0, WK, radius), ceiling(WK * radius) + 10, "automatic M")
  if (allocated(error)) return
  call pw_evanescent_from_wvn(nev, kev)
  nthread_saved = NTHREAD
  NTHREAD = 4

  opts%qtf_fs_modes = 2
  mmax = fs_auto_modes(opts%qtf_fs_modes, WK, radius)
  call check(error, mmax, 2, "M with qtf_fs_modes = 2")
  if (allocated(error)) then
    NTHREAD = nthread_saved
    return
  end if
  call pw_sample_points([0d0, 0d0], radius, fx_depth, mmax, nev, circle)
  call fs_circle_modes(mesh, circle, states, .false., [0d0, 0d0], radius, fx_depth, WK, kev, mmax, nev, pw, miss)
  msg = fs_mode_warning(miss(1))
  print '(4x,a,es10.2,a,a)', "mode check with M = 2: miss ", miss(1), ", log: ", msg
  call check(error, index(msg, "QTF free surface: the partial waves miss the field on the circle by ") == 1 .and. &
             index(msg, "; raise qtf_fs_modes or qtf_fs_evanescent") > 0, "no warning with qtf_fs_modes = 2: '"//msg//"'")
  if (allocated(error)) then
    NTHREAD = nthread_saved
    return
  end if

  opts%qtf_fs_modes = 0
  mmax = fs_auto_modes(opts%qtf_fs_modes, WK, radius)
  call pw_sample_points([0d0, 0d0], radius, fx_depth, mmax, nev, circle)
  call fs_circle_modes(mesh, circle, states, .false., [0d0, 0d0], radius, fx_depth, WK, kev, mmax, nev, pw, miss)
  NTHREAD = nthread_saved
  miss_auto = miss(1)
  msg = fs_mode_warning(miss_auto)
  print '(4x,a,i0,a,es10.2)', "mode check with the automatic M = ", mmax, ": miss ", miss_auto
  write (info, '(a,es10.3,a,a)') "automatic M misses the samples by ", miss_auto, " or warns: ", msg
  call check(error, miss_auto <= FS_MODE_TOL .and. miss_auto > 0d0 .and. len(msg) == 0, trim(info))
end subroutine test_mode_check_warns

!> The mode check holds every state to 1e-3 of itself (review I1). Three states on the fixture, R = 15 m, M = 6, nev = 8:
!>   1. the MacCamy-Fuchs scattered wave, which 6 angular modes resolve (k0 R = 2.3);
!>   2. the field of one hull panel near the waterline (phi_P = 1e-4 there, 0 elsewhere): local, so 6 modes truncate it,
!>      and 1e-4 times smaller than state 1 on the circle;
!>   3. state 2 times 1e-12: zero up to rounding, as psi of yaw on an axisymmetric body.
!> State 2 must warn, though its miss is far below 1e-3 of the largest sample of all states (the scale of 9fdd193, under
!> which it was silent); states 1 and 3 must not (state 3 is held to 1e-6 of the largest, fs_mode_scales).
subroutine test_mode_check_per_state(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: radius = 15d0
  integer, parameter :: nev = 8, mmax = 6
  type(panel_mesh) :: mesh
  type(first_order_state) :: states(3), one(1)
  type(cs_points) :: circle
  type(pw_set) :: pw(3)
  complex(8), allocatable :: samples(:, :)
  real(8) :: kev(nev), miss(3), sc(3), big, old2
  integer :: j, nthread_saved
  character(len=200) :: info

  call cylinder_fixture(mesh, one, error)
  if (allocated(error)) return
  states = one(1)
  j = minloc(norm2(mesh%xc - spread([10d0, 0d0, -0.5d0], 2, mesh%n), dim=1), dim=1)
  states(2)%phi_p = (0d0, 0d0)
  states(2)%dphi_p_dn = (0d0, 0d0)
  states(2)%phi_p(j) = (1d-4, 0d0)
  states(3) = states(2)
  states(3)%phi_p(j) = (1d-16, 0d0)
  call pw_evanescent_from_wvn(nev, kev)
  call pw_sample_points([0d0, 0d0], radius, fx_depth, mmax, nev, circle)
  nthread_saved = NTHREAD
  NTHREAD = 4
  call fs_circle_modes(mesh, circle, states, .false., [0d0, 0d0], radius, fx_depth, WK, kev, mmax, nev, pw, miss)
  allocate (samples(circle%n, 3))
  call fs_sample_phi(mesh, circle, states, .false., samples, centroid_only=.true.)
  NTHREAD = nthread_saved
  sc = fs_mode_scales(samples)
  big = maxval(abs(samples))
  old2 = miss(2) * sc(2) / big
  print '(4x,a,3es10.2,a,es10.2,a,3es10.2)', "mode check per state: miss ", miss, "; state 2 against the largest of all ", &
    old2, "; scales ", sc
  write (info, '(a,3es10.3)') "mode check per state (1 resolved, 2 truncated, 3 rounding): ", miss
  call check(error, miss(1) <= FS_MODE_TOL .and. miss(2) > FS_MODE_TOL .and. miss(3) <= FS_MODE_TOL, trim(info))
  if (allocated(error)) return
  write (info, '(a,es10.3)') "state 2 would warn under the scale of all states too: the test does not test ", old2
  call check(error, old2 <= FS_MODE_TOL, trim(info))
  if (allocated(error)) return
  call check(error, len(fs_mode_warning(maxval(miss))) > 0, "no warning for the truncated state")
end subroutine test_mode_check_per_state

!> The automatic N_ev (qtf_fs_evanescent = 0): max(8, ceil(9 h / (pi (R - r_wl)))), a given value as it is.
subroutine test_auto_evanescent(error)
  type(error_type), allocatable, intent(out) :: error

  ! 9 x 100 / (5 pi) = 57.30
  call check(error, fs_auto_evanescent(0, 100d0, 15d0, 10d0), 58, "h 100, R 15, r_wl 10")
  if (allocated(error)) return
  ! 9 x 15 / (5 pi) = 8.59
  call check(error, fs_auto_evanescent(0, 15d0, 15d0, 10d0), 9, "h 15, R 15, r_wl 10")
  if (allocated(error)) return
  ! 9 x 10 / (20 pi) = 1.43: the floor of 8
  call check(error, fs_auto_evanescent(0, 10d0, 30d0, 10d0), 8, "h 10, R 30, r_wl 10")
  if (allocated(error)) return
  call check(error, fs_auto_evanescent(5, 100d0, 15d0, 10d0), 5, "a given N_ev")
end subroutine test_auto_evanescent


end module test_qtf_free_surface
