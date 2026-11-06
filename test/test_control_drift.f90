!> Unit tests for the control-surface mean drift, Phase A (F_x, F_y, M_z), with the fields on C built directly, no BEM solve:
!> the incident wave alone, which carries no mean load, and the exact MacCamy-Fuchs field of a fixed bottom-mounted cylinder,
!> whose drift force is known from the far field.
module test_control_drift
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderTypes, only : panel_mesh, first_order_state, cs_points, cs_fields
  use FirstOrderState, only : incident_phi, incident_grad, wave_gravity
  use FarFieldDrift, only : group_ratio
  use ControlSurfaceMesh, only : cylinder_surface, box_surface, add_surface_points
  use ControlSurfaceDrift, only : control_surface_pair, CS_T_LINE, CS_T_FLUX, CS_T_PRESSURE
  use second_order_fixtures, only : maccamy_fuchs_field, maccamy_fuchs_drift
  implicit none
  private

  public :: collect_control_drift

  real(8), parameter :: pi = 3.141592653589793d0
  real(8), parameter :: g = 9.81d0, rho = 1025d0

  ! Incident wave alone: C of radius 30 and depth 30 in h = 100. The limit is on the load over rho g A^2 R (A = 1).
  real(8), parameter :: inc_radius = 30d0, inc_depth = 30d0, inc_h = 100d0
  real(8), parameter :: inc_half(2) = [30d0, 20d0]   ! the box C; R is its larger half-length
  real(8), parameter :: inc_omegas(3) = [0.8d0, 1.2d0, 2.0d0]
  ! Panel size of the coarse C; the fine C has half of it. At 1 m and 2 rad/s the pair term converges at a ratio 2.97 per
  ! halving, not yet in its asymptotic range (measured; the wall rows must resolve the e^{2kz} decay, 2k = 0.82/m)
  real(8), parameter :: inc_panel = 0.5d0
  real(8), parameter :: incident_limit = 1d-3
  real(8), parameter :: halving_ratio = 3d0           ! the error must fall by at least this when the panel size is halved

  ! MacCamy-Fuchs: cylinder a = 5 in h = 10, C to the seabed
  real(8), parameter :: mf_a = 5d0, mf_h = 10d0
  real(8), parameter :: mf_omegas(2) = [1.0d0, 1.6d0]
  real(8), parameter :: mf_panel = 0.5d0
  real(8), parameter :: mf_limit = 5d-3               ! F_x against the far-field value, relative
  real(8), parameter :: mf_lateral_limit = 1d-6       ! |F_y| and |M_z| over F_x
  real(8), parameter :: mf_c_limit = 5d-3             ! F_x of different C against each other, relative
  real(8), parameter :: symmetry_limit = 1d-12

contains

subroutine collect_control_drift(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("incident_only_zero", test_incident_only_zero), &
    new_unittest("incident_only_zero_box", test_incident_only_zero_box), &
    new_unittest("maccamy_fuchs_exact", test_maccamy_fuchs_exact), &
    new_unittest("maccamy_fuchs_independent_of_c", test_maccamy_fuchs_independent_of_c), &
    new_unittest("pair_symmetry", test_pair_symmetry) &
    ]
end subroutine collect_control_drift

! ---------------------------------------------------------------------------------------------- helpers

!> Wave number from omega^2/g = k tanh(kh), by Newton iteration.
real(8) function wave_number(omega, h)
  real(8), intent(in) :: omega, h
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

!> A state with only the wave data the control-surface loads read (omega, g, rho), and the rest for the fields.
subroutine wave_state(omega, h, beta, s)
  real(8), intent(in) :: omega, h, beta
  type(first_order_state), intent(out) :: s

  s%omega = omega
  s%k = wave_number(omega, h)
  s%nu = omega**2 / g
  s%h = h
  s%amp = 1d0
  s%beta = beta
  s%g = g
  s%rho = rho
end subroutine wave_state

!> Points of a cylinder C (shape 1) or box C (shape 2) centred on the origin, as body 1.
subroutine surface_points(shape, size_a, size_b, depth, panel, h, pts)
  integer, intent(in) :: shape
  real(8), intent(in) :: size_a, size_b, depth, panel, h
  type(cs_points), intent(out) :: pts
  type(panel_mesh) :: cm

  if (shape == 1) then
    call cylinder_surface([0d0, 0d0], size_a, depth, panel, h, cm)
  else
    call box_surface([0d0, 0d0], [size_a, size_b], 0d0, depth, panel, h, cm)
  end if
  call add_surface_points(cm, 1, pts)
end subroutine surface_points

!> The incident wave of s at every point of pts, with g_w = wave_gravity, as FieldVelocity takes it.
subroutine incident_fields(pts, s, f)
  type(cs_points), intent(in) :: pts
  type(first_order_state), intent(in) :: s
  type(cs_fields), intent(out) :: f
  real(8) :: g_w
  integer :: i

  allocate (f%phi(pts%n), f%grad(3, pts%n))
  g_w = wave_gravity(s%omega, s%k, s%h)
  do i = 1, pts%n
    f%phi(i) = incident_phi(pts%x(:, i), s%omega, s%k, s%h, s%amp, s%beta, g_w)
    call incident_grad(pts%x(:, i), s%omega, s%k, s%h, s%amp, s%beta, g_w, f%grad(:, i))
  end do
end subroutine incident_fields

!> The exact MacCamy-Fuchs field of the cylinder a = mf_a at every point of pts.
subroutine mf_fields(pts, s, f)
  type(cs_points), intent(in) :: pts
  type(first_order_state), intent(in) :: s
  type(cs_fields), intent(out) :: f
  integer :: i

  allocate (f%phi(pts%n), f%grad(3, pts%n))
  do i = 1, pts%n
    call maccamy_fuchs_field(pts%x(:, i), mf_a, s%h, s%k, s%omega, s%amp, s%beta, s%g, f%phi(i), f%grad(:, i))
  end do
end subroutine mf_fields

!> (F_x, F_y, M_z) of body 1 about the origin, the sum of the three Phase A terms.
function loads(pts, s1, s2, f1, f2) result(f)
  type(cs_points), intent(in) :: pts
  type(first_order_state), intent(in) :: s1, s2
  type(cs_fields), intent(in) :: f1, f2
  complex(8) :: f(3)
  complex(8) :: terms(6, 6, 1)
  real(8) :: xr(3, 1)

  xr = 0d0
  call control_surface_pair(pts, xr, s1, s2, f1, f2, terms)
  f(1) = terms(1, CS_T_LINE, 1) + terms(1, CS_T_FLUX, 1) + terms(1, CS_T_PRESSURE, 1)
  f(2) = terms(2, CS_T_LINE, 1) + terms(2, CS_T_FLUX, 1) + terms(2, CS_T_PRESSURE, 1)
  f(3) = terms(6, CS_T_LINE, 1) + terms(6, CS_T_FLUX, 1) + terms(6, CS_T_PRESSURE, 1)
end function loads

!> Worst load of the incident wave alone, over the frequencies: of the mean load of each heading given, |F_x|, |F_y|, |M_z| over
!> rho g R (worst_mean), and of the pair term of each heading with heading 0, |F_x|, |F_y| over rho g R and |M_z| over
!> rho g R^2 (worst_pair).
subroutine incident_worst(shape, panel, betas, worst_mean, worst_pair)
  integer, intent(in) :: shape
  real(8), intent(in) :: panel, betas(:)
  real(8), intent(out) :: worst_mean, worst_pair
  type(cs_points) :: pts
  type(first_order_state) :: s, s0
  type(cs_fields) :: f, f0
  complex(8) :: load(3)
  real(8) :: r
  integer :: io, ib

  if (shape == 1) then
    call surface_points(1, inc_radius, 0d0, inc_depth, panel, inc_h, pts)
    r = inc_radius
  else
    call surface_points(2, inc_half(1), inc_half(2), inc_depth, panel, inc_h, pts)
    r = maxval(inc_half)
  end if
  worst_mean = 0d0
  worst_pair = 0d0
  do io = 1, size(inc_omegas)
    call wave_state(inc_omegas(io), inc_h, 0d0, s0)
    call incident_fields(pts, s0, f0)
    do ib = 1, size(betas)
      call wave_state(inc_omegas(io), inc_h, betas(ib), s)
      call incident_fields(pts, s, f)
      load = loads(pts, s, s, f, f)
      worst_mean = max(worst_mean, maxval(abs(load)) / (rho * g * r))
      if (betas(ib) == 0d0) cycle
      load = loads(pts, s, s0, f, f0)
      ! The moment over rho g R^2, so that it has the units of the forces over rho g R
      worst_pair = max(worst_pair, maxval(abs(load(1:2))) / (rho * g * r), abs(load(3)) / (rho * g * r**2))
    end do
  end do
end subroutine incident_worst

!> The mean load of one incident wave has integrands that do not vary in x and y (|phi|^2 of a progressive wave depends on z
!> only), so on a closed C made of vertical rows it cancels to round-off at any panel size [measured: 4e-16 to 2e-15]. The
!> panel-size convergence is therefore checked on the pair term of two incident headings, whose integrands vary around C and
!> whose exact value is also zero [derived: the identity is bilinear in the two fields]. The limit holds on the mean loads at
!> both panel sizes and on the pair term of the fine C.
subroutine check_incident(shape, betas, name, error)
  integer, intent(in) :: shape
  real(8), intent(in) :: betas(:)
  character(len=*), intent(in) :: name
  type(error_type), allocatable, intent(out) :: error
  real(8) :: coarse, fine, coarse_pair, fine_pair

  call incident_worst(shape, inc_panel, betas, coarse, coarse_pair)
  call incident_worst(shape, 0.5d0 * inc_panel, betas, fine, fine_pair)
  print '(4x,a,": mean load / (rho g R) ",es10.3," at panel ",f5.2,", ",es10.3," at panel ",f5.2)', &
    name, coarse, inc_panel, fine, 0.5d0 * inc_panel
  print '(4x,a,": pair term with heading 0 / (rho g R) ",es10.3," at panel ",f5.2,", ",es10.3," at panel ",f5.2, &
    &", ratio ",f6.2)', name, coarse_pair, inc_panel, fine_pair, 0.5d0 * inc_panel, coarse_pair / fine_pair
  call check(error, coarse < incident_limit .and. fine < incident_limit, name//": the incident wave carries a mean load")
  if (allocated(error)) return
  ! The coarse C (0.5 m) is checked by the ratio only: its pair term is 1.9e-3 at 2 rad/s (measured)
  call check(error, fine_pair < incident_limit, name//": two incident waves carry a pair load")
  if (allocated(error)) return
  call check(error, coarse_pair / fine_pair >= halving_ratio, name//": the error does not fall with the panel size")
end subroutine check_incident

!> F_x of the MacCamy-Fuchs field at heading 0 on the given C, and the worst |F_y|, |M_z| over it.
subroutine mf_load(shape, size_a, size_b, omega, fx, lateral)
  integer, intent(in) :: shape
  real(8), intent(in) :: size_a, size_b, omega
  real(8), intent(out) :: fx, lateral
  type(cs_points) :: pts
  type(first_order_state) :: s
  type(cs_fields) :: f
  complex(8) :: load(3)

  call surface_points(shape, size_a, size_b, mf_h, mf_panel, mf_h, pts)
  call wave_state(omega, mf_h, 0d0, s)
  call mf_fields(pts, s, f)
  load = loads(pts, s, s, f, f)
  fx = real(load(1), 8)
  lateral = max(abs(load(2)), abs(load(3))) / abs(fx)
end subroutine mf_load

! ---------------------------------------------------------------------------------------------- tests

!> The incident wave alone carries no mean load through a closed C [derived: the mean momentum flux of a wave with no body is
!> divergence-free, and on z = 0 the term of F is a line integral that cancels the Gamma_c term].
subroutine test_incident_only_zero(error)
  type(error_type), allocatable, intent(out) :: error
  call check_incident(1, [0d0, 0.7d0], "incident_only_zero", error)
end subroutine test_incident_only_zero

subroutine test_incident_only_zero_box(error)
  type(error_type), allocatable, intent(out) :: error
  call check_incident(2, [0.7d0], "incident_only_zero_box", error)
end subroutine test_incident_only_zero_box

!> The exact field of a fixed bottom-mounted cylinder gives the far-field drift force [Chen 2007, section 5: the middle-field
!> and far-field forms are equal for the exact field].
subroutine test_maccamy_fuchs_exact(error)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: fx, lateral, fref, k
  integer :: io

  do io = 1, size(mf_omegas)
    call mf_load(1, 8d0, 0d0, mf_omegas(io), fx, lateral)
    k = wave_number(mf_omegas(io), mf_h)
    fref = maccamy_fuchs_drift(mf_a, k, 0.5d0 * rho * g * group_ratio(k, mf_h), 720)
    print '(4x,a,f5.2,a,es14.6,a,es14.6,a,es10.3,a,es10.3)', "maccamy_fuchs_exact omega ", mf_omegas(io), ": F_x ", fx, &
      " far ", fref, " relative ", abs(fx - fref) / fref, " lateral ", lateral
    call check(error, abs(fx - fref) < mf_limit * fref, "F_x differs from the MacCamy-Fuchs far-field value")
    if (allocated(error)) return
    call check(error, lateral < mf_lateral_limit, "F_y or M_z is not zero for the symmetric field")
    if (allocated(error)) return
  end do
end subroutine test_maccamy_fuchs_exact

!> F_x does not depend on C: cylinders of radius 6, 8 and 15, and a box of half-lengths 7 x 9.
subroutine test_maccamy_fuchs_independent_of_c(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: radii(3) = [6d0, 8d0, 15d0]
  real(8) :: fx(4), lateral, spread
  integer :: io, ic

  do io = 1, size(mf_omegas)
    do ic = 1, 3
      call mf_load(1, radii(ic), 0d0, mf_omegas(io), fx(ic), lateral)
    end do
    call mf_load(2, 7d0, 9d0, mf_omegas(io), fx(4), lateral)
    spread = (maxval(fx) - minval(fx)) / maxval(abs(fx))
    print '(4x,a,f5.2,a,4es14.6,a,es10.3)', "maccamy_fuchs_independent_of_c omega ", mf_omegas(io), ": F_x ", fx, &
      " spread ", spread
    call check(error, spread < mf_c_limit, "F_x depends on the control surface")
    if (allocated(error)) return
  end do
end subroutine test_maccamy_fuchs_independent_of_c

!> The pair term of (s1, s2) is the conjugate of that of (s2, s1), term by term [derived: pair_minus(f2, g2, f1, g1) is the
!> conjugate of pair_minus(f1, g1, f2, g2), and the weights are real].
subroutine test_pair_symmetry(error)
  type(error_type), allocatable, intent(out) :: error
  type(cs_points) :: pts
  type(first_order_state) :: s1, s2
  type(cs_fields) :: f1, f2
  complex(8) :: t12(6, 6, 1), t21(6, 6, 1)
  real(8) :: xr(3, 1), worst

  call surface_points(1, 8d0, 0d0, mf_h, 1d0, mf_h, pts)
  call wave_state(1.2d0, mf_h, 0d0, s1)
  call wave_state(1.2d0, mf_h, 0.7d0, s2)
  call mf_fields(pts, s1, f1)
  call mf_fields(pts, s2, f2)
  xr(:, 1) = [1d0, -2d0, 0d0]
  call control_surface_pair(pts, xr, s1, s2, f1, f2, t12)
  call control_surface_pair(pts, xr, s2, s1, f2, f1, t21)
  worst = maxval(abs(t12 - conjg(t21))) / maxval(abs(t12))
  print '(4x,a,es10.3)', "pair_symmetry: worst |T12 - conj(T21)| / max |T12| ", worst
  call check(error, maxval(abs(t12)) > 0d0, "The pair term is zero")
  if (allocated(error)) return
  call check(error, worst < symmetry_limit, "The pair term is not conjugate-symmetric")
end subroutine test_pair_symmetry

end module test_control_drift
