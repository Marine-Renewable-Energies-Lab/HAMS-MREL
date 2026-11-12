!> Unit tests of QTFFreeSurface and of the free-surface forcing kernel of SecondOrderIncident (Phase 6a, Tasks 1, 6 and 7).
!> The field tests of Task 6 set module globals (the mesh arrays and the wave of WaveDyn_mod), so the suite runs serially
!> (see unit_tester.f90).
module test_qtf_free_surface
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderIncident, only : incident2, incident2_pair, incident2_phi, incident2_grad, free_surface_forcing
  use NearFieldDrift, only : KERNEL_DIFF
  use QTFFreeSurface, only : fs_point_field, fs_forcing_diff, fs_inner_fields, fs_inner_steps, fs_circle_modes, &
                             fs_mode_warning, fs_auto_modes, FS_MODE_TOL, fs_sample_phi, fs_mode_scales, fs_auto_evanescent, &
                             fs_inner_pair, fs_outer_pair, fs_outer_annulus, fs_outer_tail, fs_outer_ring, fs_outer_rfar, &
                             fs_outer_sliver, fs_radius_message, FS_SIGN, fs_pair_order_message
  use HAMS_mod, only : NTHREAD, ISOL, IRSP
  use WaveDyn_mod, only : AMP, BETA, H, NK, TP, V, W1, WK, WVN
  use Body_mod, only : XW
  use SecondOrderTypes, only : so_options, panel_mesh, first_order_state, cs_points, cs_fields, CS_FREE, surface_geom
  use FieldVelocity, only : fields_green, incident_at
  use QTFField, only : body_winding
  use PartialWaves, only : pw_set, pw_sample_points, pw_evanescent_from_wvn, pw_eval_surface, pw_radial, pw_radial_at, &
                           pw_bessel_k
  use FirstOrderState, only : incident_phi, incident_grad, wave_gravity
  use SurfaceGeometry, only : build_surface_geom, DEFAULT_SHARP_DEG
  use ControlSurfaceMesh, only : build_fs_inner
  use second_order_fixtures, only : make_cylinder, maccamy_fuchs_state
  use test_waterline_field, only : set_wave, load_panels
  use Const_mod, only : g_const => G
  implicit none
  private

  public :: collect_qtf_free_surface

  real(8), parameter :: g = 9.81d0, pi = 3.14159265358979323846d0, deg = pi / 180d0
  ! Depth of the cylinder fixture of the field tests (as pw_hams_field: 5 m of water under the keel)
  real(8), parameter :: fx_depth = 15d0
  ! The synthetic partial waves of the outer-integral tests (Task 7): depth, circle, angular and evanescent modes
  real(8), parameter :: rho = 1025d0
  real(8), parameter :: sy_h = 20d0, sy_radius = 12d0, sy_centre(2) = [1.5d0, -0.7d0]
  integer, parameter :: sy_mmax = 3, sy_nev = 1

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
    new_unittest("fs_auto_evanescent", test_auto_evanescent), &
    new_unittest("fs_radius_at_waterline_stops", test_radius_at_waterline), &
    new_unittest("fs_inner_steps_floor", test_inner_steps_floor), &
    new_unittest("fs_inner_pair_sum", test_inner_pair_sum), &
    new_unittest("fs_outer_plane_wave_angular", test_outer_plane_wave_angular), &
    new_unittest("fs_outer_vs_brute_force", test_outer_vs_brute_force), &
    new_unittest("fs_outer_tail", test_outer_tail), &
    new_unittest("fs_outer_pair_order", test_pair_order), &
    new_unittest("fs_inner_quadrature_converges", test_inner_quadrature), &
    new_unittest("fs_inner_outer_radius", test_inner_outer_radius) &
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


! ---------------------------------------------------------------------------------------------- integrals (Task 7)

!> A circle R <= r_wl, the waterline radius, stops the run (strictly: the equality too), and fs_auto_evanescent does not
!> divide by R - r_wl there (it returns -1). fs_radius_message gives the stop message, empty when R clears the waterline.
subroutine test_radius_at_waterline(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=:), allocatable :: m_eq, m_in, m_ok

  m_eq = fs_radius_message(15d0, 15d0)
  m_in = fs_radius_message(15d0, 15d0 + 1d-9)
  m_ok = fs_radius_message(15d0, 15d0 - 1d-9)
  print '(4x,a,a)', "R = r_wl: ", m_eq
  call check(error, index(m_eq, "waterline") > 0 .and. index(m_in, "waterline") > 0, "R <= r_wl does not stop: "//m_eq)
  if (allocated(error)) return
  call check(error, len(m_ok) == 0, "R just outside r_wl stops: "//m_ok)
  if (allocated(error)) return
  call check(error, fs_auto_evanescent(0, 15d0, 15d0, 15d0), -1, "automatic N_ev at R = r_wl")
  if (allocated(error)) return
  call check(error, fs_auto_evanescent(0, 15d0, 15d0, 16d0), -1, "automatic N_ev at R < r_wl")
end subroutine test_radius_at_waterline

!> fs_inner_steps at its floor (review of Task 6): a point 2 mm outside a waterline node of the fixture, point size 10 m,
!> step 0.05. Its stencil at the floor, 1e-3 x 0.5 m, clears the waterplane, but the clearance of 16 offsets never does, so
!> delta stops at the floor, the point counts in nshrunk and nfloor, and the run does not stop.
subroutine test_inner_steps_floor(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: step = 0.05d0, psz = 10d0, gap = 2d-3
  type(panel_mesh) :: mesh
  type(first_order_state) :: states(1)
  type(cs_points) :: pts
  real(8) :: x(3, 1), delta(1), rn
  integer :: j, jn, nshrunk, nfloor
  character(len=200) :: info

  call cylinder_fixture(mesh, states, error)
  if (allocated(error)) return
  jn = 0
  do j = 1, mesh%nnode
    rn = hypot(mesh%xnode(1, j), mesh%xnode(2, j))
    if (abs(mesh%xnode(3, j)) < 1d-9 .and. abs(rn - 10d0) < 1d-9) then
      jn = j
      exit
    end if
  end do
  call check(error, jn > 0, "no waterline node on r = 10 m")
  if (allocated(error)) return
  x(:, 1) = [mesh%xnode(1, jn) * (1d0 + gap / 10d0), mesh%xnode(2, jn) * (1d0 + gap / 10d0), 0d0]
  call free_points(x, psz, pts)
  call fs_inner_steps(mesh, pts, step, delta, nshrunk, nfloor)
  print '(4x,a,es12.4,a,i0,a,i0)', "fs_inner_steps 2 mm off a waterline node: delta ", delta(1), ", shrunk ", nshrunk, &
    ", at the floor ", nfloor
  write (info, '(a,es12.4,2i3)') "delta, nshrunk, nfloor at the floor ", delta(1), nshrunk, nfloor
  call check(error, nshrunk == 1 .and. nfloor == 1 .and. abs(delta(1) - 1d-3 * step * psz) <= 1d-12 * step * psz, trim(info))
  if (allocated(error)) return
  call check(error, .not. stencil_in_waterplane(mesh, x(:, 1), delta(1)), "the stencil at the floor reaches the waterplane")
end subroutine test_inner_steps_floor

!> fs_inner_pair is the sum s_FS (i Omega rho / g) sum_i w_i psi_k(x_i) Q_D(x_i), Q_D = fs_forcing_diff, with s_FS = +1.
subroutine test_inner_pair_sum(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w1 = 0.9d0, w2 = 0.55d0, wts(3) = [0.7d0, 1.9d0, 0.25d0]
  type(fs_point_field) :: f1(3), i1(3), f2(3), i2(3)
  complex(8) :: psi(3, 2), t(2), ref(2)
  integer :: i, k

  do i = 1, 3
    f1(i) = pseudo_field(10 * i + 1)
    i1(i) = pseudo_field(10 * i + 2)
    f2(i) = pseudo_field(10 * i + 3)
    i2(i) = pseudo_field(10 * i + 4)
    do k = 1, 2
      psi(i, k) = cmplx(cos(0.3d0 * i + k), sin(1.1d0 * i * k), 8)
    end do
  end do
  call fs_inner_pair(w1, f1, i1, w2, f2, i2, psi, wts, g, rho, t)
  ref = (0d0, 0d0)
  do i = 1, 3
    ref = ref + wts(i) * psi(i, :) * fs_forcing_diff(w1, f1(i), i1(i), w2, f2(i), i2(i), g)
  end do
  ref = ref * cmplx(0d0, (w1 - w2) * rho / g, 8)
  call check(error, FS_SIGN == 1d0, "s_FS is not +1")
  if (allocated(error)) return
  call check(error, maxval(abs(t - ref)) <= 1d-14 * maxval(abs(ref)), "fs_inner_pair is not the weighted sum")
end subroutine test_inner_pair_sum

!> The angular integral of the parts with a plane wave, Jacobi-Anger (int e^{ikr cos(theta - beta)} g dtheta =
!> 2 pi sum g_n i^n J_n(kr) e^{+in beta}), against the trapezoid rule on 2048 angles of psi_k Q[1, 2] from the pointwise
!> fields, at kr = 5 and 50, to 1e-10 relative: incident 1 with scattered 2 (no scattered 1, no incident 2), then
!> scattered 1 with incident 2. Synthetic partial waves (M = 3, one evanescent mode, h = 20, R = 12, off-origin centre).
subroutine test_outer_plane_wave_angular(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w1 = 1d0, w2 = 0.8d0, b1 = 0.3d0, b2 = -0.5d0, krs(2) = [5d0, 50d0]
  integer, parameter :: nth = 2048
  type(pw_set) :: pw1, pw2, pz1, pz2, psi(2)
  type(fs_point_field) :: fa, fb
  complex(8) :: ring(2), direct(2), q, ps, gr(3), zz
  real(8) :: k, r, x(3), th, worst, err
  integer :: part, ik, j, l
  character(len=200) :: info

  call synthetic_set(w1, w2, pw1, pw2, psi)
  pz1 = pw1
  pz1%a = (0d0, 0d0)
  pz2 = pw2
  pz2%a = (0d0, 0d0)
  worst = 0d0
  do part = 1, 2
    do ik = 1, 2
      if (part == 1) then
        k = pw1%k0
        r = krs(ik) / k
        call fs_outer_ring(w1, pz1, b1, k, w2, pw2, b2, 0d0, psi, g, r, ring)
      else
        k = pw2%k0
        r = krs(ik) / k
        call fs_outer_ring(w1, pw1, b1, 0d0, w2, pz2, b2, k, psi, g, r, ring)
      end if
      direct = (0d0, 0d0)
      do j = 0, nth - 1
        th = 2d0 * pi * j / nth
        x = [sy_centre(1) + r * cos(th), sy_centre(2) + r * sin(th), 0d0]
        if (part == 1) then
          fa = incident_point(w1, k, b1, x)
          fb = scattered_point(pw2, x)
        else
          fa = scattered_point(pw1, x)
          fb = incident_point(w2, k, b2, x)
        end if
        q = free_surface_forcing(KERNEL_DIFF, w1, fa%phi, fa%grad, fa%phizz, w2, fb%phi, fb%grad, fb%phizz, g)
        do l = 1, 2
          call pw_eval_surface(psi(l), x, ps, gr, zz)
          direct(l) = direct(l) + ps * q * (2d0 * pi / nth)
        end do
      end do
      err = maxval(abs(ring - direct)) / maxval(abs(direct))
      print '(4x,a,i0,a,f5.1,a,2es12.4,a,es10.2)', "plane-wave part ", part, " at kr =", krs(ik), ": |ring| ", abs(ring), &
        ", relative error ", err
      worst = max(worst, err)
    end do
  end do
  write (info, '(a,es10.3)') "Jacobi-Anger angular integral against the trapezoid rule, worst relative ", worst
  call check(error, worst <= 1d-10, trim(info))
end subroutine test_outer_plane_wave_angular

!> fs_outer_annulus from R to R + 200 m (no tail) against a brute-force 2D quadrature of the same integrand on a polar grid
!> (Gauss-Legendre, 10 points per metre in r; trapezoid on 64 or 256 angles), with the fields pointwise from
!> pw_eval_surface and the exact incident wave and Q_D from fs_forcing_diff, to 1e-8 relative. Synthetic modal fields
!> (M = 3, one evanescent mode): first with no incident wave (the brief's case), then with both incident waves.
subroutine test_outer_vs_brute_force(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w1 = 1d0, w2 = 0.8d0, b1 = 0.3d0, b2 = -0.5d0
  type(pw_set) :: pw1, pw2, psi(2)
  complex(8) :: t(2), tb(2)
  real(8) :: kk1, kk2, err(2)
  integer :: ic, nth
  character(len=200) :: info

  call synthetic_set(w1, w2, pw1, pw2, psi)
  do ic = 1, 2
    if (ic == 1) then
      kk1 = 0d0
      kk2 = 0d0
      nth = 64
    else
      kk1 = pw1%k0
      kk2 = pw2%k0
      nth = 256
    end if
    call fs_outer_annulus(w1, pw1, b1, kk1, w2, pw2, b2, kk2, psi, w1 - w2, g, rho, sy_radius, sy_radius + 200d0, t)
    call brute_outer(w1, pw1, b1, kk1, w2, pw2, b2, kk2, psi, w1 - w2, sy_radius, sy_radius + 200d0, 1d0, 10, nth, tb)
    err(ic) = maxval(abs(t - tb)) / maxval(abs(tb))
    print '(4x,a,i0,a,2es14.6,a,es10.2)', "outer annulus case ", ic, ": |T| ", abs(tb), ", relative error ", err(ic)
  end do
  write (info, '(a,2es10.3)') "fs_outer_annulus against the brute-force polar grid (no incident, incident) ", err
  call check(error, maxval(err) <= 1d-8, trim(info))
end subroutine test_outer_vs_brute_force

!> The tail beyond r_far (contour rotation per phase class, complex-argument Hankel functions) against the real-axis
!> integral: the outer integral with r_far (fs_outer_rfar) and with 2 r_far agree to 1e-6 relative, i.e. tail(r_far) =
!> annulus(r_far, 2 r_far) + tail(2 r_far), where the annulus uses the real-argument intrinsics only. The annulus alone
!> is above 1e-4 of the total, so a missing tail fails. fs_outer_pair = annulus(R, r_far) + tail(r_far). Two pairs on the
!> synthetic fields with both incident waves: Omega = 0.2 and 0.05 rad/s (small K r_far for psi).
subroutine test_outer_tail(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: wa(2) = [1d0, 0.85d0], wb(2) = [0.8d0, 0.8d0], b1 = 0.3d0, b2 = -0.5d0
  type(pw_set) :: pw1, pw2, psi(2)
  complex(8) :: a_in(2), a_mid(2), tl1(2), tl2(2), tp(2), t1(2)
  real(8) :: rf, err, sig, comp, worst_err, worst_comp, least_sig, om
  integer :: ip
  character(len=200) :: info

  worst_err = 0d0
  worst_comp = 0d0
  least_sig = huge(1d0)
  do ip = 1, 2
    om = wa(ip) - wb(ip)
    call synthetic_set(wa(ip), wb(ip), pw1, pw2, psi)
    rf = fs_outer_rfar(pw1, pw1%k0, pw2, pw2%k0, psi)
    call fs_outer_annulus(wa(ip), pw1, b1, pw1%k0, wb(ip), pw2, b2, pw2%k0, psi, om, g, rho, sy_radius, rf, a_in)
    call fs_outer_annulus(wa(ip), pw1, b1, pw1%k0, wb(ip), pw2, b2, pw2%k0, psi, om, g, rho, rf, 2d0 * rf, a_mid)
    call fs_outer_tail(wa(ip), pw1, b1, pw1%k0, wb(ip), pw2, b2, pw2%k0, psi, om, g, rho, rf, tl1)
    call fs_outer_tail(wa(ip), pw1, b1, pw1%k0, wb(ip), pw2, b2, pw2%k0, psi, om, g, rho, 2d0 * rf, tl2)
    call fs_outer_pair(wa(ip), pw1, b1, pw1%k0, wb(ip), pw2, b2, pw2%k0, psi, om, g, rho, rf, tp)
    t1 = a_in + tl1
    err = maxval(abs(tl1 - a_mid - tl2)) / maxval(abs(t1))
    sig = maxval(abs(a_mid)) / maxval(abs(t1))
    comp = maxval(abs(tp - t1)) / maxval(abs(t1))
    print '(4x,a,f5.2,a,f10.1,a,2es12.4,a,es10.2,a,es10.2,a,es10.2)', "tail, Omega", om, ": r_far", rf, &
      " m, |T| ", abs(t1), ", |tail(r_far)| / |T| ", maxval(abs(tl1)) / maxval(abs(t1)), ", r_far vs 2 r_far ", err, &
      ", annulus r_far..2 r_far ", sig
    worst_err = max(worst_err, err)
    worst_comp = max(worst_comp, comp)
    least_sig = min(least_sig, sig)
  end do
  write (info, '(a,es10.3)') "outer integral with r_far against 2 r_far, worst relative ", worst_err
  call check(error, worst_err <= 1d-6, trim(info))
  if (allocated(error)) return
  write (info, '(a,es10.3)') "the annulus r_far..2 r_far is below 1e-4 of the total: the test does not test ", least_sig
  call check(error, least_sig > 1d-4, trim(info))
  if (allocated(error)) return
  write (info, '(a,es10.3)') "fs_outer_pair is not annulus + tail ", worst_comp
  call check(error, worst_comp <= 1d-12, trim(info))
end subroutine test_outer_tail

!> The mesh-level radius check: on the 448-panel cylinder in depth 15, inner(R1) + outer(R1) = inner(R2) + outer(R2) for
!> R1 = 15 m and R2 = 20 m, to 1e-2 relative (the brief's gate) on the inner mesh of panel size lambda/10 at omega 3.11.
!> Fields: the MacCamy-Fuchs wave at (1.1, heading 0.6) and (0.9, heading 0.2), amp 1 (the incident wave exactly, the
!> scattered part from the hull panels); psi: the scattered MacCamy-Fuchs wave at Omega = 0.2 (heading 0), a field that
!> satisfies the free-surface condition at Omega and radiates, from the hull panels. Inner: the generated inner mesh
!> (build_fs_inner), fs_inner_fields and fs_sample_phi as the driver; outer: the partial waves of the circle samples
!> (fs_circle_modes), fs_outer_pair to fs_outer_rfar, and the sliver between the inner polygon and the circle
!> (fs_outer_sliver). The sum over the whole free surface outside the body cannot depend on R.
!> The frequencies are om_a, om_b, not w1, w2: a local w1 would hide the global W1 of WaveDyn_mod, and the states would be
!> built at the wrong frequency (the cause of the 2 % plateau of the first version, fix round 1 of Task 7).
!> Measured (fix round 1), relative difference against the inner mesh of omega 1.1 (the default panel size), 2.2, 3.11 and
!> 4.4: 8.5e-2, 1.28e-2, 1.90e-3, 2.35e-3. It falls as the inner panel sum converges (second order,
!> fs_inner_quadrature_converges) down to about 2e-3.
subroutine test_inner_outer_radius(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: om_a = 1.1d0, om_b = 0.9d0, b1 = 0.6d0, b2 = 0.2d0, radii(2) = [15d0, 20d0], step = 0.04d0
  real(8), parameter :: meshes(3) = [1.1d0, 2.2d0, 3.1112698372208d0]
  real(8) :: errs(3)
  integer :: im
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(so_options) :: opts
  type(first_order_state) :: st(1), sp(1), s1(1), s2(1)
  type(cs_points) :: pts, circle
  type(fs_point_field), allocatable :: f1(:, :), f2(:, :), i1(:), i2(:)
  complex(8), allocatable :: psi(:, :)
  real(8), allocatable :: delta(:), kev(:)
  type(pw_set) :: pw1(1), pw2(1), pwp(1)
  complex(8) :: tin(1), tout(1), tsl(1), total(2), nosliver(2), phi, grad(3)
  real(8) :: centre(2), radius, miss(1), k1, k2, k1max, rf, err, err_ns
  integer :: ir, i, mmax, nev, nshrunk, npc, nthread_saved
  logical :: ok
  character(len=:), allocatable :: msg
  character(len=200) :: info

  call cylinder_fixture(mesh, st, error)
  if (allocated(error)) return
  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "geometry: "//msg)
  if (allocated(error)) return
  nthread_saved = NTHREAD
  NTHREAD = 8
  do im = 1, 3
  do ir = 1, 2
    opts%qtf_fs_radius = radii(ir)
    call build_fs_inner(mesh, geom, opts, fx_depth, 'fs_task7_no_input_dir', pts, centre, radius, ok, msg, &
                        omega_max=meshes(im), np_circle=npc)
    call check(error, ok .and. npc >= 8, "inner mesh: "//msg)
    if (allocated(error)) exit
    call set_frequency(om_a, b1, fx_depth)
    k1max = WK
    mmax = fs_auto_modes(0, k1max, radius)
    nev = fs_auto_evanescent(0, fx_depth, radius, 10d0)
    if (allocated(delta)) deallocate (delta, psi, f1, f2, i1, i2)
    if (allocated(kev)) deallocate (kev)
    allocate (delta(pts%n), psi(pts%n, 1), f1(pts%n, 1), f2(pts%n, 1), i1(pts%n), i2(pts%n), kev(nev))
    call fs_inner_steps(mesh, pts, step, delta, nshrunk)
    call pw_sample_points(centre, radius, fx_depth, mmax, nev, circle)

    ! psi at Omega: the scattered MacCamy-Fuchs wave, perturbation only
    call set_frequency(om_a - om_b, 0d0, fx_depth)
    call maccamy_fuchs_state(mesh, 10d0, fx_depth, WK, W1, 1d0, 0d0, wave_gravity(W1, WK, fx_depth), rho, sp(1))
    sp(1)%amp = 0d0
    call fs_sample_phi(mesh, pts, sp, .false., psi)
    call pw_evanescent_from_wvn(nev, kev)
    call fs_circle_modes(mesh, circle, sp, .false., centre, radius, fx_depth, WK, kev, mmax, nev, pwp, miss)
    print '(4x,a,f5.1,a,i0,a,i0,a,i0,a,i0,a,es9.2)', "R =", radius, ": ", pts%n, " inner points (polygon of ", npc, &
      "), M ", mmax, ", N_ev ", nev, "; psi mode check ", miss(1)

    ! The two first-order states
    call set_frequency(om_a, b1, fx_depth)
    k1 = WK
    call maccamy_fuchs_state(mesh, 10d0, fx_depth, WK, W1, 1d0, b1, wave_gravity(W1, WK, fx_depth), rho, s1(1))
    call fs_inner_fields(mesh, pts, s1, step, .false., f1, delta)
    do i = 1, pts%n
      call incident_at(s1(1), pts%x(:, i), phi, grad)
      i1(i)%phi = phi
      i1(i)%grad = grad
      i1(i)%phizz = k1**2 * phi
    end do
    call pw_evanescent_from_wvn(nev, kev)
    call fs_circle_modes(mesh, circle, s1, .false., centre, radius, fx_depth, WK, kev, mmax, nev, pw1, miss)
    call set_frequency(om_b, b2, fx_depth)
    k2 = WK
    call maccamy_fuchs_state(mesh, 10d0, fx_depth, WK, W1, 1d0, b2, wave_gravity(W1, WK, fx_depth), rho, s2(1))
    call fs_inner_fields(mesh, pts, s2, step, .false., f2, delta)
    do i = 1, pts%n
      call incident_at(s2(1), pts%x(:, i), phi, grad)
      i2(i)%phi = phi
      i2(i)%grad = grad
      i2(i)%phizz = k2**2 * phi
    end do
    call pw_evanescent_from_wvn(nev, kev)
    call fs_circle_modes(mesh, circle, s2, .false., centre, radius, fx_depth, WK, kev, mmax, nev, pw2, miss)

    call fs_inner_pair(om_a, f1(:, 1), i1, om_b, f2(:, 1), i2, psi, pts%w, g, rho, tin)
    rf = fs_outer_rfar(pw1(1), k1, pw2(1), k2, pwp)
    call fs_outer_pair(om_a, pw1(1), b1, k1, om_b, pw2(1), b2, k2, pwp, om_a - om_b, g, rho, rf, tout)
    call fs_outer_sliver(om_a, pw1(1), b1, k1, om_b, pw2(1), b2, k2, pwp, om_a - om_b, g, rho, npc, tsl)
    total(ir) = tin(1) + tout(1) + tsl(1)
    nosliver(ir) = tin(1) + tout(1)
    print '(4x,a,f5.1,a,2es13.5,a,2es13.5,a,2es13.5,a,f8.1)', "R =", radius, ": inner ", tin(1), ", outer ", tout(1), &
      ", sliver ", tsl(1), ", r_far ", rf
  end do
  if (allocated(error)) exit
  errs(im) = abs(total(1) - total(2)) / abs(total(1))
  err_ns = abs(nosliver(1) - nosliver(2)) / abs(nosliver(1))
  print '(4x,a,f4.1,a,2es13.5,a,es10.2,a,es10.2)', "inner mesh of omega ", meshes(im), ": total ", total(1), &
    ", R 15 against R 20: relative ", errs(im), "; without the sliver ", err_ns
  end do
  NTHREAD = nthread_saved
  if (allocated(error)) return
  err = errs(1)
  write (info, '(a,3es10.3)') "inner + outer at R = 15 against R = 20, relative (inner mesh of omega 1.1, 2.2, 3.11) ", errs
  call check(error, errs(3) <= 1d-2, trim(info))
end subroutine test_inner_outer_radius

! ---------------------------------------------------------------------------------------------- helpers (Task 7)

!> The wave globals of a frequency, as set_wave of test_waterline_field, at omega and heading in depth.
subroutine set_frequency(omega, heading, depth)
  real(8), intent(in) :: omega, heading, depth
  external :: DISPERSION

  ISOL = 1
  IRSP = 0
  AMP = 1d0
  BETA = heading
  XW = 0d0
  W1 = omega
  TP = 2d0 * pi / W1
  H = depth
  V = W1**2 / g_const
  call DISPERSION(WVN, NK, W1, H)
  WK = WVN(1)
end subroutine set_frequency

!> A synthetic pw_set at omega (k0 and kev from DISPERSION in depth sy_h, circle sy_radius about sy_centre): every mode
!> is a reproducible complex number of size 1/(1 + |m|) at r = R (half that for the evanescent mode).
subroutine synthetic_pw(omega, mmax, nev, seed, pw)
  real(8), intent(in) :: omega
  integer, intent(in) :: mmax, nev, seed
  type(pw_set), intent(out) :: pw
  real(8) :: wvn_l(nev + 1), kv(0:mmax)
  complex(8) :: c, hk
  integer :: m, n
  external :: DISPERSION

  call DISPERSION(wvn_l, nev + 1, omega, sy_h)
  pw%mmax = mmax
  pw%nev = nev
  pw%k0 = wvn_l(1)
  pw%h = sy_h
  pw%radius = sy_radius
  pw%centre = sy_centre
  allocate (pw%kev(nev), pw%a(-mmax:mmax, 0:nev))
  pw%kev = wvn_l(2:nev + 1)
  do n = 0, nev
    if (n > 0) call pw_bessel_k(mmax, pw%kev(n) * sy_radius, kv)
    do m = -mmax, mmax
      c = cmplx(sin(1.3d0 * seed + 0.7d0 * m + 1.1d0 * n), cos(0.9d0 * seed - 0.4d0 * m + 0.3d0 * n), 8) / (1 + abs(m))
      if (n == 0) then
        hk = cmplx(bessel_jn(abs(m), pw%k0 * sy_radius), bessel_yn(abs(m), pw%k0 * sy_radius), 8)
        if (m < 0 .and. mod(abs(m), 2) == 1) hk = -hk
        pw%a(m, 0) = c / hk
      else
        pw%a(m, n) = 0.5d0 * c / kv(abs(m))
      end if
    end do
  end do
end subroutine synthetic_pw

!> The synthetic fields of the outer tests: states 1 and 2 at w1 and w2, and two psi at Omega = w1 - w2.
subroutine synthetic_set(w1, w2, pw1, pw2, psi)
  real(8), intent(in) :: w1, w2
  type(pw_set), intent(out) :: pw1, pw2, psi(2)

  call synthetic_pw(w1, sy_mmax, sy_nev, 1, pw1)
  call synthetic_pw(w2, sy_mmax, sy_nev, 2, pw2)
  call synthetic_pw(w1 - w2, sy_mmax, sy_nev, 3, psi(1))
  call synthetic_pw(w1 - w2, sy_mmax, sy_nev, 4, psi(2))
end subroutine synthetic_set

!> The incident wave (omega, wave number k, heading beta in radians, depth sy_h, unit amplitude) at x as a point field;
!> zero for k <= 0.
function incident_point(omega, k, beta, x, depth) result(f)
  real(8), intent(in) :: omega, k, beta, x(3)
  real(8), intent(in), optional :: depth
  type(fs_point_field) :: f
  real(8) :: gw, hh

  f%phi = (0d0, 0d0)
  f%grad = (0d0, 0d0)
  f%phizz = (0d0, 0d0)
  if (k <= 0d0) return
  hh = sy_h
  if (present(depth)) hh = depth
  gw = wave_gravity(omega, k, hh)
  f%phi = incident_phi(x, omega, k, hh, 1d0, beta, gw)
  call incident_grad(x, omega, k, hh, 1d0, beta, gw, f%grad)
  f%phizz = k**2 * f%phi
end function incident_point

!> The field of the partial waves pw at x (z = 0) as a point field, from pw_eval_surface.
function scattered_point(pw, x, rad) result(f)
  type(pw_set), intent(in) :: pw
  real(8), intent(in) :: x(3)
  type(pw_radial), intent(in), optional :: rad
  type(fs_point_field) :: f

  if (present(rad)) then
    call pw_eval_surface(pw, x, f%phi, f%grad, f%phizz, rad)
  else
    call pw_eval_surface(pw, x, f%phi, f%grad, f%phizz)
  end if
end function scattered_point

!> a + b, componentwise
function add_fields(a, b) result(f)
  type(fs_point_field), intent(in) :: a, b
  type(fs_point_field) :: f

  f%phi = a%phi + b%phi
  f%grad = a%grad + b%grad
  f%phizz = a%phizz + b%phizz
end function add_fields

!> Gauss-Legendre nodes and weights of order n on [0, 1] (Newton on P_n).
subroutine gauss_legendre01(n, x, w)
  integer, intent(in) :: n
  real(8), intent(out) :: x(n), w(n)
  real(8) :: t, p0, p1, p2, dp
  integer :: i, j, it

  do i = 1, n
    t = cos(pi * (i - 0.25d0) / (n + 0.5d0))
    do it = 1, 100
      p0 = 1d0
      p1 = t
      do j = 2, n
        p2 = ((2 * j - 1) * t * p1 - (j - 1) * p0) / j
        p0 = p1
        p1 = p2
      end do
      dp = n * (t * p1 - p0) / (t * t - 1d0)
      t = t - p1 / dp
      if (abs(p1 / dp) < 1d-16) exit
    end do
    p0 = 1d0
    p1 = t
    do j = 2, n
      p2 = ((2 * j - 1) * t * p1 - (j - 1) * p0) / j
      p0 = p1
      p1 = p2
    end do
    dp = n * (t * p1 - p0) / (t * t - 1d0)
    x(i) = 0.5d0 * (1d0 - t)
    w(i) = 1d0 / ((1d0 - t * t) * dp * dp)
  end do
end subroutine gauss_legendre01

!> Brute force: the outer integral of s_FS (i Omega rho / g) psi_k Q_D over r_a < r < r_b on a polar grid about sy_centre,
!> Gauss-Legendre of order nq on panels of about dr in r, the trapezoid rule on nth angles; the fields pointwise
!> (pw_eval_surface, the exact incident wave), Q_D by fs_forcing_diff.
subroutine brute_outer(w1, pw1, b1, k1, w2, pw2, b2, k2, psi, omega, r_a, r_b, dr, nq, nth, t)
  ! The circle centre and depth are those of pw1
  real(8), intent(in) :: w1, b1, k1, w2, b2, k2, omega, r_a, r_b, dr
  type(pw_set), intent(in) :: pw1, pw2, psi(:)
  integer, intent(in) :: nq, nth
  complex(8), intent(out) :: t(size(psi))
  real(8) :: xq(nq), wq(nq), r, x(3), th, len_p
  type(pw_radial) :: rad1, rad2, radp
  type(fs_point_field) :: i1, i2, t1, t2
  complex(8) :: q, ps, gr(3), zz
  integer :: ip, iq, j, k, np

  call gauss_legendre01(nq, xq, wq)
  np = max(1, nint((r_b - r_a) / dr))
  len_p = (r_b - r_a) / np
  t = (0d0, 0d0)
  do ip = 1, np
    do iq = 1, nq
      r = r_a + (ip - 1 + xq(iq)) * len_p
      call pw_radial_at(pw1, r, rad1)
      call pw_radial_at(pw2, r, rad2)
      call pw_radial_at(psi(1), r, radp)
      do j = 0, nth - 1
        th = 2d0 * pi * j / nth
        x = [pw1%centre(1) + r * cos(th), pw1%centre(2) + r * sin(th), 0d0]
        i1 = incident_point(w1, k1, b1, x, pw1%h)
        i2 = incident_point(w2, k2, b2, x, pw1%h)
        t1 = add_fields(scattered_point(pw1, x, rad1), i1)
        t2 = add_fields(scattered_point(pw2, x, rad2), i2)
        q = fs_forcing_diff(w1, t1, i1, w2, t2, i2, g)
        do k = 1, size(psi)
          call pw_eval_surface(psi(k), x, ps, gr, zz, radp)
          t(k) = t(k) + wq(iq) * len_p * r * (2d0 * pi / nth) * ps * q
        end do
      end do
    end do
  end do
  t = t * cmplx(0d0, omega * rho / g, 8)
end subroutine brute_outer



! ---------------------------------------------------------------------------------------------- Task 7 fix round 1

!> The outer routines refuse a pair that is not ordered w1 > w2 with Omega > 0 (fs_pair_order_message gives the stop
!> message; Omega = 0 is returned as 0 before the check by the routines).
subroutine test_pair_order(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=:), allocatable :: m_ok, m_rev, m_neg, m_eq

  m_ok = fs_pair_order_message(1d0, 0.8d0, 0.2d0)
  m_rev = fs_pair_order_message(0.8d0, 1d0, 0.2d0)
  m_neg = fs_pair_order_message(1d0, 0.8d0, -0.2d0)
  m_eq = fs_pair_order_message(0.8d0, 0.8d0, 0d0)
  print '(4x,a,a)', "w1 < w2: ", m_rev
  call check(error, len(m_ok) == 0, "an ordered pair is refused: "//m_ok)
  if (allocated(error)) return
  call check(error, index(m_rev, "order the pair") > 0 .and. index(m_neg, "order the pair") > 0 .and. &
             index(m_eq, "order the pair") > 0, "an unordered pair or Omega <= 0 is not refused")
end subroutine test_pair_order

!> The partial waves at R of the fixture: psi (scattered MacCamy-Fuchs at Omega = om_a - om_b, heading 0) and the
!> scattered parts of the states at (om_a, b1) and (om_b, b2), from the circle samples (fs_circle_modes), M and N_ev
!> automatic. The frequencies are not called w1, w2 here: W1 is the global frequency of WaveDyn_mod.
subroutine fixture_modes(mesh, radius, om_a, b1, om_b, b2, pw1, pw2, pwp, k1, k2)
  type(panel_mesh), intent(in) :: mesh
  real(8), intent(in) :: radius, om_a, b1, om_b, b2
  type(pw_set), intent(out) :: pw1(1), pw2(1), pwp(1)
  real(8), intent(out) :: k1, k2
  type(first_order_state) :: s1(1), s2(1), sp(1)
  type(cs_points) :: circle
  real(8), allocatable :: kev(:)
  real(8) :: miss(1)
  integer :: mmax, nev

  call set_frequency(om_a, b1, fx_depth)
  mmax = fs_auto_modes(0, WK, radius)
  nev = fs_auto_evanescent(0, fx_depth, radius, 10d0)
  allocate (kev(nev))
  call pw_sample_points([0d0, 0d0], radius, fx_depth, mmax, nev, circle)
  call set_frequency(om_a - om_b, 0d0, fx_depth)
  call maccamy_fuchs_state(mesh, 10d0, fx_depth, WK, W1, 1d0, 0d0, wave_gravity(W1, WK, fx_depth), rho, sp(1))
  sp(1)%amp = 0d0
  call pw_evanescent_from_wvn(nev, kev)
  call fs_circle_modes(mesh, circle, sp, .false., [0d0, 0d0], radius, fx_depth, WK, kev, mmax, nev, pwp, miss)
  call set_frequency(om_a, b1, fx_depth)
  k1 = WK
  call maccamy_fuchs_state(mesh, 10d0, fx_depth, WK, W1, 1d0, b1, wave_gravity(W1, WK, fx_depth), rho, s1(1))
  call pw_evanescent_from_wvn(nev, kev)
  call fs_circle_modes(mesh, circle, s1, .false., [0d0, 0d0], radius, fx_depth, WK, kev, mmax, nev, pw1, miss)
  call set_frequency(om_b, b2, fx_depth)
  k2 = WK
  call maccamy_fuchs_state(mesh, 10d0, fx_depth, WK, W1, 1d0, b2, wave_gravity(W1, WK, fx_depth), rho, s2(1))
  call pw_evanescent_from_wvn(nev, kev)
  call fs_circle_modes(mesh, circle, s2, .false., [0d0, 0d0], radius, fx_depth, WK, kev, mmax, nev, pw2, miss)
end subroutine fixture_modes

!> psi Q_D at x (z = 0) of the partial waves pw1, pw2 (scattered), the incident waves and psi pwp, pointwise:
!> pw_eval_surface (phi_zz exact per mode), incident_point, fs_forcing_diff
complex(8) function modal_integrand(x, om_a, pw1, b1, k1, om_b, pw2, b2, k2, pwp) result(f)
  real(8), intent(in) :: x(3), om_a, b1, k1, om_b, b2, k2
  type(pw_set), intent(in) :: pw1, pw2, pwp
  type(fs_point_field) :: i1, i2, t1, t2
  complex(8) :: ps, gr(3), zz

  i1 = incident_point(om_a, k1, b1, x, fx_depth)
  i2 = incident_point(om_b, k2, b2, x, fx_depth)
  t1 = add_fields(scattered_point(pw1, x), i1)
  t2 = add_fields(scattered_point(pw2, x), i2)
  call pw_eval_surface(pwp, x, ps, gr, zz)
  f = ps * fs_forcing_diff(om_a, t1, i1, om_b, t2, i2, g)
end function modal_integrand

!> A C-infinity bump on (15.5, 19.5) m, 1 at 17.5 m
real(8) function radial_bump(r)
  real(8), intent(in) :: r
  real(8) :: u

  radial_bump = 0d0
  u = (r - 17.5d0) / 2d0
  if (abs(u) < 1d0) radial_bump = exp(1d0 - 1d0 / (1d0 - u * u))
end function radial_bump

!> The panel sum of the generated inner mesh (build_fs_inner, R = 20 m) converges as O(h^2) on a smooth integrand (review
!> of Task 7, diagnostics i and ii): the partial-wave integrand of the fixture at R = 15 m (pw_eval_surface, exact phi_zz,
!> the exact incident waves, fs_forcing_diff) under a smooth radial bump on (15.5, 19.5) m, summed over the points of the
!> meshes with panel size ps, ps/2 and ps/4 (lambda/10 at omega 2.2, 2.2 sqrt 2, 4.4), against a polar Gauss reference
!> (10 points per 0.25 m, 512 angles). Then the same with the evanescent modes n >= 4 zeroed.
!> Measured: 1.41e-2, 1.82e-3, 3.21e-4 (all modes) and 1.41e-2, 1.82e-3, 3.20e-4 (n < 4): second order, and not driven
!> by the steep evanescent modes. Gates: the finest below 1e-3, and at least 16x below the coarsest.
subroutine test_inner_quadrature(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: om_a = 1.1d0, om_b = 0.9d0, b1 = 0.6d0, b2 = 0.2d0
  real(8), parameter :: wm(3) = [2.2d0, 3.1112698372208d0, 4.4d0]
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(so_options) :: opts
  type(first_order_state) :: st(1)
  type(pw_set) :: pw1(1), pw2(1), pwp(1), z1, z2, zp
  type(cs_points) :: pts
  real(8) :: k1, k2, centre(2), radius, xq(10), wq(10), r, th, x(3), rr, e(2, 3)
  complex(8) :: ref(2), msum(2)
  integer :: im, i, j, l, npc, nthread_saved
  logical :: ok
  character(len=:), allocatable :: msg
  character(len=200) :: info

  call cylinder_fixture(mesh, st, error)
  if (allocated(error)) return
  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  nthread_saved = NTHREAD
  NTHREAD = 8
  call fixture_modes(mesh, 15d0, om_a, b1, om_b, b2, pw1, pw2, pwp, k1, k2)
  NTHREAD = nthread_saved
  z1 = pw1(1)
  z2 = pw2(1)
  zp = pwp(1)
  z1%a(:, 4:) = (0d0, 0d0)
  z2%a(:, 4:) = (0d0, 0d0)
  zp%a(:, 4:) = (0d0, 0d0)
  call gauss_legendre01(10, xq, wq)
  ref = (0d0, 0d0)
  do i = 1, 16
    do j = 1, 10
      r = 15.5d0 + (i - 1 + xq(j)) * 0.25d0
      do l = 0, 511
        th = 2d0 * pi * l / 512
        x = [r * cos(th), r * sin(th), 0d0]
        ref(1) = ref(1) + radial_bump(r) * wq(j) * 0.25d0 * r * (2d0 * pi / 512) * &
          modal_integrand(x, om_a, pw1(1), b1, k1, om_b, pw2(1), b2, k2, pwp(1))
        ref(2) = ref(2) + radial_bump(r) * wq(j) * 0.25d0 * r * (2d0 * pi / 512) * &
          modal_integrand(x, om_a, z1, b1, k1, om_b, z2, b2, k2, zp)
      end do
    end do
  end do
  opts%qtf_fs_radius = 20d0
  do im = 1, 3
    call build_fs_inner(mesh, geom, opts, fx_depth, 'fs_task7_no_input_dir', pts, centre, radius, ok, msg, &
                        omega_max=wm(im), np_circle=npc)
    call check(error, ok, "inner mesh: "//msg)
    if (allocated(error)) return
    msum = (0d0, 0d0)
    do i = 1, pts%n
      rr = norm2(pts%x(1:2, i) - centre)
      if (radial_bump(rr) == 0d0) cycle
      msum(1) = msum(1) + pts%w(i) * radial_bump(rr) * modal_integrand(pts%x(:, i), om_a, pw1(1), b1, k1, om_b, &
        pw2(1), b2, k2, pwp(1))
      msum(2) = msum(2) + pts%w(i) * radial_bump(rr) * modal_integrand(pts%x(:, i), om_a, z1, b1, k1, om_b, z2, b2, k2, zp)
    end do
    e(:, im) = abs(msum - ref) / abs(ref)
    print '(4x,a,f6.3,a,i0,a,2es11.3)', "inner mesh of omega ", wm(im), " (", pts%n, &
      " points): bump sum against polar, all modes and n < 4 ", e(:, im)
  end do
  write (info, '(a,6es10.3)') "panel sum of a smooth integrand at ps, ps/2, ps/4 (all modes, n < 4) ", e
  call check(error, maxval(e(:, 3)) <= 1d-3 .and. minval(e(:, 1) / e(:, 3)) >= 16d0, trim(info))
end subroutine test_inner_quadrature

end module test_qtf_free_surface
