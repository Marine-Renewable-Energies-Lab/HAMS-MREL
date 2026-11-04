!> Unit tests for the far-field mean drift, against the exact MacCamy-Fuchs solution (no BEM solve).
module test_far_field
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderTypes, only : panel_mesh, first_order_state, so_options
  use FarFieldDrift, only : far_amplitude, group_ratio, far_field_pair
  use second_order_fixtures, only : make_cylinder, maccamy_fuchs_a, maccamy_fuchs_state
  implicit none
  private

  public :: collect_far_field

  real(8), parameter :: pi = 3.141592653589793d0
  real(8), parameter :: g = 9.81d0, rho = 1025d0
  real(8), parameter :: radius = 5d0, depth = 10d0, omega = 1d0

contains

subroutine collect_far_field(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("kochin_matches_series", test_kochin_matches_series), &
    new_unittest("drift_matches_analytic", test_drift_matches_analytic), &
    new_unittest("energy_identity", test_energy_identity), &
    new_unittest("lateral_zero", test_lateral_zero), &
    new_unittest("deep_limit", test_deep_limit), &
    new_unittest("derivative_fd", test_derivative_fd), &
    new_unittest("conjugate_symmetry", test_conjugate_symmetry), &
    new_unittest("diagonal_consistent", test_diagonal_consistent), &
    new_unittest("reference_shift", test_reference_shift), &
    new_unittest("group_ratio_limits", test_group_ratio_limits) &
    ]
end subroutine collect_far_field

! ---------------------------------------------------------------------------------------------- helpers

!> Wave number from omega^2/g = k tanh(kh), by Newton iteration.
real(8) function wave_number(h)
  real(8), intent(in) :: h
  real(8) :: nu, k, fk
  integer :: it
  nu = omega**2 / g
  k = max(nu, sqrt(nu / h))
  do it = 1, 100
    fk = k * tanh(k * h) - nu
    k = k - fk / (tanh(k * h) + k * h / cosh(k * h)**2)
    if (abs(fk) < 1d-15) exit
  end do
  wave_number = k
end function wave_number

!> Bottom-mounted cylinder (a = 5, h = 10) and its MacCamy-Fuchs state at heading beta, per unit amplitude.
subroutine setup_cylinder(mesh, k, beta, s)
  type(panel_mesh), intent(out) :: mesh
  real(8), intent(out) :: k
  real(8), intent(in) :: beta
  type(first_order_state), intent(out) :: s
  call make_cylinder(radius, depth, 120, 20, 1, mesh, .false.)
  k = wave_number(depth)
  call maccamy_fuchs_state(mesh, radius, depth, k, omega, 1d0, beta, g, rho, s)
end subroutine setup_cylinder

!> Trapezoid integral of |f|^2 weight(theta) on [0, 2pi) of the analytic series, with n points.
real(8) function series_integral(k, wcos, n)
  real(8), intent(in) :: k
  logical, intent(in) :: wcos   ! false: int |a|^2, true: int |a|^2 cos(theta)
  integer, intent(in) :: n
  integer :: j
  real(8) :: th
  series_integral = 0d0
  do j = 1, n
    th = 2d0 * pi * real(j - 1, 8) / real(n, 8)
    series_integral = series_integral + abs(maccamy_fuchs_a(th, radius, k, 1d0, 0d0))**2 * merge(cos(th), 1d0, wcos)
  end do
  series_integral = series_integral * 2d0 * pi / real(n, 8)
end function series_integral

! ---------------------------------------------------------------------------------------------- tests

subroutine test_kochin_matches_series(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(first_order_state) :: s
  complex(8) :: a, da, ref
  real(8) :: k, th, err, amax
  integer :: j

  call setup_cylinder(mesh, k, 0d0, s)
  err = 0d0
  amax = 0d0
  do j = 0, 35
    th = 2d0 * pi * real(j, 8) / 36d0
    call far_amplitude(mesh, s, th, a, da)
    ref = maccamy_fuchs_a(th, radius, k, 1d0, 0d0)
    err = max(err, abs(a - ref))
    amax = max(amax, abs(ref))
  end do
  call check(error, err / amax < 1d-2, "Kochin amplitude differs from the MacCamy-Fuchs series")
end subroutine test_kochin_matches_series

subroutine test_drift_matches_analytic(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(first_order_state) :: s
  type(so_options) :: opts
  complex(8) :: f(3)
  real(8) :: k, fref, kappa

  call setup_cylinder(mesh, k, 0d0, s)
  call far_field_pair(mesh, s, s, opts, f)
  kappa = 0.5d0 * rho * g * group_ratio(k, depth)
  fref = kappa * (series_integral(k, .false., 720) - series_integral(k, .true., 720))
  call check(error, abs(real(f(1)) - fref) < 2d-2 * fref, "drift force differs from the analytic MacCamy-Fuchs value")
end subroutine test_drift_matches_analytic

subroutine test_energy_identity(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(first_order_state) :: s
  complex(8) :: a, da
  real(8) :: k, x_int, x_fwd
  integer, parameter :: n = 720
  integer :: j

  call setup_cylinder(mesh, k, 0d0, s)
  x_int = 0d0
  do j = 1, n
    call far_amplitude(mesh, s, 2d0 * pi * real(j - 1, 8) / real(n, 8), a, da)
    x_int = x_int + abs(a)**2
  end do
  x_int = x_int * 2d0 * pi / real(n, 8)
  call far_amplitude(mesh, s, 0d0, a, da)
  x_fwd = -2d0 * sqrt(2d0 * pi / k) * real(a * exp(cmplx(0d0, pi / 4d0, kind=8)))
  call check(error, abs(x_int - x_fwd) < 1d-2 * x_int, "forward scattering does not match the scattered energy")
end subroutine test_energy_identity

subroutine test_lateral_zero(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(first_order_state) :: s
  type(so_options) :: opts
  complex(8) :: f(3)
  real(8) :: k

  call setup_cylinder(mesh, k, 0d0, s)
  call far_field_pair(mesh, s, s, opts, f)
  call check(error, abs(f(2)) < 1d-6 * abs(f(1)), "lateral drift force is not zero for a symmetric body")
  if (allocated(error)) return
  call check(error, abs(f(3)) / radius < 1d-6 * abs(f(1)), "yaw moment is not zero for a symmetric body")
end subroutine test_lateral_zero

subroutine test_deep_limit(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(first_order_state) :: sd, sf
  complex(8) :: ad, daf, af, dad
  real(8) :: k, th, err, amax
  integer :: i, j

  call make_cylinder(radius, 3d0, 40, 6, 3, mesh, .true.)
  k = omega**2 / g
  sd%omega = omega; sd%k = k; sd%nu = k; sd%h = -1d0; sd%amp = 1d0; sd%beta = 0d0; sd%g = g; sd%rho = rho
  sf = sd
  sf%h = 30d0 / k
  allocate(sd%phi_p(mesh%n), sd%dphi_p_dn(mesh%n))
  do i = 1, mesh%n
    sd%phi_p(i) = cmplx(cos(0.3d0 * mesh%xc(1, i)), sin(0.2d0 * mesh%xc(2, i) + mesh%xc(3, i)), kind=8)
    sd%dphi_p_dn(i) = cmplx(sin(0.1d0 * mesh%xc(3, i)), cos(0.25d0 * mesh%xc(1, i) - 0.1d0 * mesh%xc(2, i)), kind=8)
  end do
  sf%phi_p = sd%phi_p
  sf%dphi_p_dn = sd%dphi_p_dn
  err = 0d0
  amax = 0d0
  do j = 0, 7
    th = 2d0 * pi * real(j, 8) / 8d0
    call far_amplitude(mesh, sd, th, ad, dad)
    call far_amplitude(mesh, sf, th, af, daf)
    err = max(err, abs(ad - af))
    amax = max(amax, abs(ad))
  end do
  call check(error, amax > 0d0 .and. err / amax < 1d-8, "large-kh amplitude differs from deep water")
end subroutine test_deep_limit

subroutine test_derivative_fd(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(first_order_state) :: s
  complex(8) :: a, da, ap, dap, am, dam
  real(8) :: k, th, dth, err, dmax
  integer :: j

  call setup_cylinder(mesh, k, 0.4d0, s)
  dth = 1d-5
  err = 0d0
  dmax = 0d0
  do j = 0, 7
    th = 2d0 * pi * real(j, 8) / 8d0 + 0.1d0
    call far_amplitude(mesh, s, th, a, da)
    call far_amplitude(mesh, s, th + dth, ap, dap)
    call far_amplitude(mesh, s, th - dth, am, dam)
    err = max(err, abs(da - (ap - am) / (2d0 * dth)))
    dmax = max(dmax, abs(da))
  end do
  call check(error, dmax > 0d0 .and. err / dmax < 1d-6, "da differs from finite differences of a")
end subroutine test_derivative_fd

subroutine test_conjugate_symmetry(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(first_order_state) :: s1, s2
  type(so_options) :: opts
  complex(8) :: f12(3), f21(3)
  real(8) :: k

  call setup_cylinder(mesh, k, 0d0, s1)
  call maccamy_fuchs_state(mesh, radius, depth, k, omega, 1d0, 0.7d0, g, rho, s2)
  call far_field_pair(mesh, s1, s2, opts, f12)
  call far_field_pair(mesh, s2, s1, opts, f21)
  call check(error, maxval(abs(f12 - conjg(f21))) <= 1d-10 * maxval(abs(f12)), "F(b1,b2) is not conj F(b2,b1)")
end subroutine test_conjugate_symmetry

subroutine test_diagonal_consistent(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(first_order_state) :: s
  type(so_options) :: opts
  complex(8) :: f(3), a, da
  real(8) :: k, kappa, beta, fx(2), ex(2), x, sq, th
  integer :: j, n

  beta = 0.4d0
  call setup_cylinder(mesh, k, beta, s)
  call far_field_pair(mesh, s, s, opts, f)
  call check(error, maxval(abs(aimag(f))) <= 1d-12 * maxval(abs(real(f))), "diagonal load is not real")
  if (allocated(error)) return

  ! Single-heading form kappa [X e_beta - int |a|^2 e_theta], X = -2 sqrt(2pi/k) Re[a(beta) e^{i pi/4}]
  kappa = 0.5d0 * rho * g * group_ratio(k, depth)
  sq = sqrt(2d0 * pi / k)
  call far_amplitude(mesh, s, beta, a, da)
  x = -2d0 * sq * real(a * exp(cmplx(0d0, pi / 4d0, kind=8)))
  n = 720
  ex = 0d0
  do j = 1, n
    th = 2d0 * pi * real(j - 1, 8) / real(n, 8)
    call far_amplitude(mesh, s, th, a, da)
    ex = ex + abs(a)**2 * [cos(th), sin(th)] * 2d0 * pi / real(n, 8)
  end do
  fx = kappa * (x * [cos(beta), sin(beta)] - ex)
  call check(error, maxval(abs(real(f(1:2)) - fx)) <= 1d-10 * maxval(abs(fx)), "diagonal load differs from the single-heading form")
end subroutine test_diagonal_consistent

subroutine test_reference_shift(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, moved
  type(first_order_state) :: s, sm
  type(so_options) :: opts
  complex(8) :: f(3), fm(3)
  real(8) :: k, d(3), beta, dxf
  complex(8) :: phase

  beta = 0d0
  call setup_cylinder(mesh, k, beta, s)
  d = [3d0, 4d0, 0d0]
  moved = mesh
  moved%xnode = moved%xnode + spread(d, 2, moved%nnode)
  moved%xc = moved%xc + spread(d, 2, moved%n)
  ! The load is unchanged and the moment about the global origin gains d x F (the body's load acts at its new place).
  ! The same physical wave reaches the moved body later by the phase e^{ik d.e_beta}
  phase = exp(cmplx(0d0, k * (d(1) * cos(beta) + d(2) * sin(beta)), kind=8))
  sm = s
  sm%phi_p = s%phi_p * phase
  sm%dphi_p_dn = s%dphi_p_dn * phase
  call far_field_pair(mesh, s, s, opts, f)
  call far_field_pair(moved, sm, sm, opts, fm)
  dxf = d(1) * real(f(2)) - d(2) * real(f(1))   ! (d x F)_z
  call check(error, maxval(abs(fm(1:2) - f(1:2))) <= 1d-6 * abs(f(1)), "force changed under a shift of the body")
  if (allocated(error)) return
  call check(error, abs(real(fm(3)) - (real(f(3)) + dxf)) <= 1d-6 * abs(f(1)) * norm2(d), &
             "yaw moment did not change by (d x F)_z under a shift of the body")
end subroutine test_reference_shift

subroutine test_group_ratio_limits(error)
  type(error_type), allocatable, intent(out) :: error
  call check(error, group_ratio(1d0, -1d0) == 0.5d0, "deep-water group ratio is not 1/2")
  if (allocated(error)) return
  call check(error, abs(group_ratio(1d-4, 1d0) - 1d0) < 1d-6, "shallow-water group ratio is not 1")
  if (allocated(error)) return
  call check(error, abs(group_ratio(1d0, 30d0) - 0.5d0) < 1d-12, "large-kh group ratio is not 1/2")
end subroutine test_group_ratio_limits

end module test_far_field
