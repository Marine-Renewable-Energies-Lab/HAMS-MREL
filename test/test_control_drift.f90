!> Unit tests for the control-surface mean drift, with the fields on C and F built directly, no BEM solve: the incident wave
!> alone, which carries no mean load (Phase A, F_x, F_y, M_z; Phase B, F_z, M_x, M_y), and the exact MacCamy-Fuchs field of a
!> fixed bottom-mounted cylinder, whose drift force is known from the far field.
module test_control_drift
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderTypes, only : panel_mesh, first_order_state, cs_points, cs_fields, CS_WALL, CS_FREE
  use FirstOrderState, only : incident_phi, incident_grad, wave_gravity
  use FarFieldDrift, only : group_ratio
  use ControlSurfaceMesh, only : cylinder_surface, box_surface, add_surface_points
  use ControlSurfaceDrift, only : control_surface_pair, CS_T_LINE, CS_T_FLUX, CS_T_PRESSURE, CS_T_FREE
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
  ! F_z of the exact field over its largest term, 1.96e-3 measured at panel 0.25 (5.8e-3 / 1.5e-3 of rho g a at 0.5 / 0.25,
  ! a ratio of 3.84): 2x measured
  real(8), parameter :: mf_vertical_limit = 4d-3

contains

subroutine collect_control_drift(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("incident_only_zero", test_incident_only_zero), &
    new_unittest("incident_only_zero_box", test_incident_only_zero_box), &
    new_unittest("incident_only_vertical_zero", test_incident_only_vertical_zero), &
    new_unittest("maccamy_fuchs_exact", test_maccamy_fuchs_exact), &
    new_unittest("maccamy_fuchs_independent_of_c", test_maccamy_fuchs_independent_of_c), &
    new_unittest("maccamy_fuchs_vertical_zero", test_maccamy_fuchs_vertical_zero), &
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

!> The free surface F inside Gamma_c of a cylinder C with a bottom, from the centre out: the bottom panels of C moved up to
!> z = 0, with the F normal (0, 0, -1). The bottom of C is the polygon of Gamma_c cut in rings from the centre, so F covers
!> exactly the disc inside Gamma_c.
subroutine add_free_disc(pts)
  type(cs_points), intent(inout) :: pts
  type(cs_points) :: old
  logical :: bottom(pts%n)
  integer :: nf, n0

  bottom = pts%kind == CS_WALL .and. pts%nrm(3, :) > 0.5d0
  nf = count(bottom)
  old = pts
  n0 = old%n
  pts%n = n0 + nf
  deallocate (pts%x, pts%nrm, pts%w, pts%psize, pts%kind, pts%body)
  allocate (pts%x(3, pts%n), pts%nrm(3, pts%n), pts%w(pts%n), pts%psize(pts%n), pts%kind(pts%n), pts%body(pts%n))
  pts%x(:, 1:n0) = old%x
  pts%nrm(:, 1:n0) = old%nrm
  pts%w(1:n0) = old%w
  pts%psize(1:n0) = old%psize
  pts%kind(1:n0) = old%kind
  pts%body(1:n0) = old%body
  pts%x(1, n0 + 1:) = pack(old%x(1, :), bottom)
  pts%x(2, n0 + 1:) = pack(old%x(2, :), bottom)
  pts%x(3, n0 + 1:) = 0d0
  pts%nrm(1:2, n0 + 1:) = 0d0
  pts%nrm(3, n0 + 1:) = -1d0
  pts%w(n0 + 1:) = pack(old%w, bottom)
  pts%psize(n0 + 1:) = pack(old%psize, bottom)
  pts%kind(n0 + 1:) = CS_FREE
  pts%body(n0 + 1:) = 1
end subroutine add_free_disc

!> (F_z, M_x, M_y) of body 1 about the origin: every term, with no hull (the hull terms are absent).
function vertical_loads(pts, s1, s2, f1, f2) result(f)
  type(cs_points), intent(in) :: pts
  type(first_order_state), intent(in) :: s1, s2
  type(cs_fields), intent(in) :: f1, f2
  complex(8) :: f(3)
  complex(8) :: terms(6, 6, 1)
  real(8) :: xr(3, 1)

  xr = 0d0
  call control_surface_pair(pts, xr, s1, s2, f1, f2, terms)
  f = sum(terms(3:5, :, 1), dim=2)
end function vertical_loads

!> As incident_worst, for the cylinder C with F, and F_z over rho g R, M_x and M_y over rho g R^2. Also returns the largest
!> |F_z| of the F term alone in the pair term over rho g R (f_scale), to show that the cancellation is not trivial. (In the
!> mean load of one wave the F integrand cancels by itself in deep water: k = nu makes eta Phi_zt + |grad Phi|^2/2 zero.)
subroutine vertical_worst(panel, betas, worst_mean, worst_pair, f_scale)
  real(8), intent(in) :: panel, betas(:)
  real(8), intent(out) :: worst_mean, worst_pair, f_scale
  type(cs_points) :: pts
  type(first_order_state) :: s, s0
  type(cs_fields) :: f, f0
  complex(8) :: load(3), terms(6, 6, 1)
  real(8) :: r, xr(3, 1)
  integer :: io, ib

  call surface_points(1, inc_radius, 0d0, inc_depth, panel, inc_h, pts)
  call add_free_disc(pts)
  r = inc_radius
  xr = 0d0
  worst_mean = 0d0
  worst_pair = 0d0
  f_scale = 0d0
  do io = 1, size(inc_omegas)
    call wave_state(inc_omegas(io), inc_h, 0d0, s0)
    call incident_fields(pts, s0, f0)
    do ib = 1, size(betas)
      call wave_state(inc_omegas(io), inc_h, betas(ib), s)
      call incident_fields(pts, s, f)
      load = vertical_loads(pts, s, s, f, f)
      worst_mean = max(worst_mean, abs(load(1)) / (rho * g * r), maxval(abs(load(2:3))) / (rho * g * r**2))
      if (betas(ib) == 0d0) cycle
      load = vertical_loads(pts, s, s0, f, f0)
      call control_surface_pair(pts, xr, s, s0, f, f0, terms)
      f_scale = max(f_scale, abs(terms(3, CS_T_FREE, 1)) / (rho * g * r))
      worst_pair = max(worst_pair, abs(load(1)) / (rho * g * r), maxval(abs(load(2:3))) / (rho * g * r**2))
    end do
  end do
end subroutine vertical_worst

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

!> Phase B: with F inside Gamma_c, the incident wave alone carries no mean F_z, M_x or M_y through C and F [derived: the
!> mean vertical momentum flux of a wave with no body is divergence-free, Chen (2007) eqs. (37)-(38) with no hull]. As for
!> F_x, the mean load of one wave cancels at any panel size, and the convergence is checked on the pair term with heading 0.
subroutine test_incident_only_vertical_zero(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: betas(2) = [0d0, 0.7d0]
  real(8) :: coarse, fine, coarse_pair, fine_pair, f_coarse, f_fine

  call vertical_worst(inc_panel, betas, coarse, coarse_pair, f_coarse)
  call vertical_worst(0.5d0 * inc_panel, betas, fine, fine_pair, f_fine)
  print '(4x,a,es10.3," at panel ",f5.2,", ",es10.3," at panel ",f5.2,"; F term of the pair alone ",es10.3)', &
    "incident_only_vertical_zero: mean load / (rho g R) ", coarse, inc_panel, fine, 0.5d0 * inc_panel, f_fine
  print '(4x,a,es10.3," at panel ",f5.2,", ",es10.3," at panel ",f5.2,", ratio ",f6.2)', &
    "incident_only_vertical_zero: pair term with heading 0 / (rho g R) ", coarse_pair, inc_panel, fine_pair, &
    0.5d0 * inc_panel, coarse_pair / fine_pair
  call check(error, f_fine > 1d1 * incident_limit, "incident_only_vertical_zero: the F term is too small to test")
  if (allocated(error)) return
  call check(error, coarse < incident_limit .and. fine < incident_limit, &
             "incident_only_vertical_zero: the incident wave carries a mean vertical load")
  if (allocated(error)) return
  ! As in check_incident, the coarse C is checked by the ratio only
  call check(error, fine_pair < incident_limit, "incident_only_vertical_zero: two incident waves carry a vertical pair load")
  if (allocated(error)) return
  call check(error, coarse_pair / fine_pair >= halving_ratio, &
             "incident_only_vertical_zero: the error does not fall with the panel size")
end subroutine test_incident_only_vertical_zero

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

!> Polar annulus between the radii r_in and r_out about the origin, with np sectors (the angles of cylinder_surface) and
!> nr rings, appended to pts as panels of the given kind and normal at height z. Quadrilateral panels: the area by the shoelace
!> formula, the point at the mean of the four corners.
subroutine add_annulus(pts, r_in, r_out, np, nr, z, kind, nrm)
  type(cs_points), intent(inout) :: pts
  real(8), intent(in) :: r_in, r_out, z, nrm(3)
  integer, intent(in) :: np, nr, kind
  type(cs_points) :: old
  real(8) :: c(2, 4), r0, r1, t0, t1
  integer :: n0, i, j, m, k

  old = pts
  n0 = old%n
  pts%n = n0 + np * nr
  deallocate (pts%x, pts%nrm, pts%w, pts%psize, pts%kind, pts%body)
  allocate (pts%x(3, pts%n), pts%nrm(3, pts%n), pts%w(pts%n), pts%psize(pts%n), pts%kind(pts%n), pts%body(pts%n))
  pts%x(:, 1:n0) = old%x
  pts%nrm(:, 1:n0) = old%nrm
  pts%w(1:n0) = old%w
  pts%psize(1:n0) = old%psize
  pts%kind(1:n0) = old%kind
  pts%body(1:n0) = old%body
  m = n0
  do j = 1, nr
    r0 = r_in + (r_out - r_in) * real(j - 1, 8) / real(nr, 8)
    r1 = r_in + (r_out - r_in) * real(j, 8) / real(nr, 8)
    do i = 1, np
      t0 = 2d0 * pi * real(i - 1, 8) / real(np, 8)
      t1 = 2d0 * pi * real(i, 8) / real(np, 8)
      c(:, 1) = r0 * [cos(t0), sin(t0)]
      c(:, 2) = r1 * [cos(t0), sin(t0)]
      c(:, 3) = r1 * [cos(t1), sin(t1)]
      c(:, 4) = r0 * [cos(t1), sin(t1)]
      m = m + 1
      pts%x(1:2, m) = sum(c, dim=2) / 4d0
      pts%x(3, m) = z
      pts%nrm(:, m) = nrm
      pts%w(m) = 0d0
      do k = 1, 4
        pts%w(m) = pts%w(m) + 0.5d0 * (c(1, k) * c(2, mod(k, 4) + 1) - c(1, mod(k, 4) + 1) * c(2, k))
      end do
      pts%psize(m) = r1 - r0
      pts%kind(m) = kind
      pts%body(m) = 1
    end do
  end do
end subroutine add_annulus

!> The exact MacCamy-Fuchs field of the fixed bottom-mounted cylinder a = mf_a: C of radius 8 down to the sea bed, F the
!> annulus between the hull and Gamma_c, and the sea bed annulus below it closing D, with its normal up (into D), as C panels.
!> The hull is wall-sided and fixed, so the Gamma and hydrostatic terms are zero and the pressure on the hull has no vertical
!> part: F_z must be zero [derived]. Returns |F_z| / (rho g a) of the worst frequency, and the largest |F_z| of one term.
subroutine mf_vertical(panel, worst, scale)
  real(8), intent(in) :: panel
  real(8), intent(out) :: worst, scale
  real(8), parameter :: rc = 8d0
  type(cs_points) :: pts
  type(first_order_state) :: s
  type(cs_fields) :: f
  complex(8) :: terms(6, 6, 1)
  real(8) :: xr(3, 1)
  integer :: io, np, nr

  call surface_points(1, rc, 0d0, mf_h, panel, mf_h, pts)
  np = max(8, ceiling(2d0 * pi * rc / panel))
  nr = max(1, ceiling((rc - mf_a) / panel))
  call add_annulus(pts, mf_a, rc, np, nr, 0d0, CS_FREE, [0d0, 0d0, -1d0])
  call add_annulus(pts, mf_a, rc, np, nr, -mf_h, CS_WALL, [0d0, 0d0, 1d0])
  xr = 0d0
  worst = 0d0
  scale = 0d0
  do io = 1, size(mf_omegas)
    call wave_state(mf_omegas(io), mf_h, 0d0, s)
    call mf_fields(pts, s, f)
    call control_surface_pair(pts, xr, s, s, f, f, terms)
    print '(4x,a,f5.2,a,f5.2,a,6es11.3)', "maccamy_fuchs_vertical_zero panel ", panel, " omega ", mf_omegas(io), &
      ": F_z terms / (rho g a) ", real(terms(3, :, 1), 8) / (rho * g * mf_a)
    worst = max(worst, abs(sum(terms(3, :, 1))) / (rho * g * mf_a))
    scale = max(scale, maxval(abs(terms(3, :, 1))) / (rho * g * mf_a))
  end do
end subroutine mf_vertical

!> F_z of the exact MacCamy-Fuchs field through C, F and the sea bed patch is zero (mf_vertical), and the error falls with
!> the panel size.
subroutine test_maccamy_fuchs_vertical_zero(error)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: coarse, fine, scale

  call mf_vertical(mf_panel, coarse, scale)
  call mf_vertical(0.5d0 * mf_panel, fine, scale)
  print '(4x,a,es10.3,a,es10.3,a,es10.3)', "maccamy_fuchs_vertical_zero: |F_z| / (rho g a) ", coarse, " and ", fine, &
    " at half the panel size; largest term ", scale
  call check(error, fine < mf_vertical_limit * scale, "maccamy_fuchs_vertical_zero: F_z of the exact field is not zero")
  if (allocated(error)) return
  call check(error, coarse / fine >= halving_ratio, "maccamy_fuchs_vertical_zero: the error does not fall with the panel size")
end subroutine test_maccamy_fuchs_vertical_zero

end module test_control_drift
