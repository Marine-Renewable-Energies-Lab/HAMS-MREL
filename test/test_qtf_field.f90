!> Unit tests of QTFField: second-order pressure and elevation of a pair of waves at a point (quadratic part plus the
!> second-order incident potential).
module test_qtf_field
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use QTFField, only : field_pair, build_field_points
  use SecondOrderTypes, only : panel_mesh, cs_points
  use second_order_fixtures, only : make_cylinder, make_two_bodies
  use SecondOrderIncident, only : incident2, incident2_pair, incident2_phi
  use FirstOrderState, only : incident_phi, incident_grad
  use NearFieldDrift, only : KERNEL_DIFF, KERNEL_SUM
  implicit none
  private

  public :: collect_qtf_field

  real(8), parameter :: g = 9.81d0, rho = 1025d0, pi = 3.14159265358979323846d0, deg = pi / 180d0
  complex(8), parameter :: ci = (0d0, 1d0)

contains

subroutine collect_qtf_field(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("field_stokes_elevation", test_stokes_elevation), &
    new_unittest("field_stokes_pressure", test_stokes_pressure), &
    new_unittest("field_mean_pressure", test_mean_pressure), &
    new_unittest("field_swap_symmetry", test_swap_symmetry), &
    new_unittest("field_bichromatic_time", test_bichromatic_time), &
    new_unittest("field_sum_continuity", test_sum_continuity), &
    new_unittest("points_classified", test_points_classified), &
    new_unittest("points_near_zero_is_surface", test_points_near_zero), &
    new_unittest("points_above_rejected", test_points_above), &
    new_unittest("points_on_hull_rejected", test_points_on_hull), &
    new_unittest("points_inside_rejected", test_points_inside), &
    new_unittest("points_inside_waterplane_rejected", test_points_waterplane), &
    new_unittest("points_two_bodies", test_points_two_bodies) &
    ]
end subroutine collect_qtf_field

!> Wave number of omega, from omega^2 = g k tanh(k h) (h < 0: deep water).
real(8) function wave_number(omega, h)
  real(8), intent(in) :: omega, h
  real(8) :: k, nu
  integer :: it
  nu = omega**2 / g
  k = nu
  if (h < 0d0) then
    wave_number = k
    return
  end if
  k = max(nu, sqrt(nu / h))
  do it = 1, 100
    k = k - (k * tanh(k * h) - nu) / (tanh(k * h) + k * h / cosh(k * h)**2)
  end do
  wave_number = k
end function wave_number

!> Second-order incident potential of the pair (omega, heading in degrees) of unit waves in depth h, under the kernel.
subroutine make_pair(w1, b1, w2, b2, h, inc, kernel)
  real(8), intent(in) :: w1, b1, w2, b2, h
  type(incident2), intent(out) :: inc
  integer, intent(in) :: kernel
  call incident2_pair(w1, wave_number(w1, h), b1 * deg, w2, wave_number(w2, h), b2 * deg, h, g, inc, kernel)
end subroutine make_pair

!> First-order potential and gradient of the unit wave (omega, heading in degrees) at x.
subroutine wave_at(w, b, h, x, phi, u)
  real(8), intent(in) :: w, b, h, x(3)
  complex(8), intent(out) :: phi, u(3)
  real(8) :: k
  k = wave_number(w, h)
  phi = incident_phi(x, w, k, h, 1d0, b * deg, g)
  call incident_grad(x, w, k, h, 1d0, b * deg, g, u)
end subroutine wave_at

!> Quadratic and incident parts of the pair (w1, b1), (w2, b2) of unit waves at x.
subroutine pair_total(kernel, w1, b1, w2, b2, h, x, surface, quad, pot)
  integer, intent(in) :: kernel
  real(8), intent(in) :: w1, b1, w2, b2, h, x(3)
  logical, intent(in) :: surface
  complex(8), intent(out) :: quad, pot
  type(incident2) :: inc
  complex(8) :: phi1, u1(3), phi2, u2(3)
  call wave_at(w1, b1, h, x, phi1, u1)
  call wave_at(w2, b2, h, x, phi2, u2)
  call make_pair(w1, b1, w2, b2, h, inc, kernel)
  call field_pair(kernel, w1, phi1, u1, w2, phi2, u2, inc, x, surface, rho, g, quad, pot)
end subroutine pair_total

!> Second-order Stokes elevation (Dean and Dalrymple 1991), finite depth and deep water.
subroutine test_stokes_elevation(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w = 0.9d0, b = 20d0
  real(8), parameter :: x(3) = [0.7d0, 0.3d0, 0d0]
  complex(8) :: quad, pot, expected, ph
  real(8) :: k, hs
  integer :: ih

  do ih = 1, 2
    hs = merge(12d0, -1d0, ih == 1)
    k = wave_number(w, hs)
    call pair_total(KERNEL_SUM, w, b, w, b, hs, x, .true., quad, pot)
    ph = exp(2d0 * ci * k * (x(1) * cos(b * deg) + x(2) * sin(b * deg)))
    if (ih == 1) then
      expected = 0.25d0 * k * cosh(k * hs) * (2d0 + cosh(2d0 * k * hs)) / sinh(k * hs)**3 * ph
    else
      expected = 0.5d0 * k * ph
    end if
    call check(error, abs(quad + pot - expected) <= 1d-10 * abs(expected), "second-order Stokes elevation")
    if (allocated(error)) return
  end do
end subroutine test_stokes_elevation

!> Second-order pressure of a Stokes wave from Bernoulli with the Stokes potential (Dean and Dalrymple 1991).
subroutine test_stokes_pressure(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w = 0.9d0, b = 20d0, h = 12d0
  real(8), parameter :: x(3) = [0.7d0, 0.3d0, -3d0]
  complex(8) :: quad, pot, expected
  real(8) :: k

  k = wave_number(w, h)
  call pair_total(KERNEL_SUM, w, b, w, b, h, x, .false., quad, pot)
  expected = rho * (0.75d0 * w**2 * cosh(2d0 * k * (x(3) + h)) / sinh(k * h)**4 - g * k / (2d0 * sinh(2d0 * k * h))) &
             * exp(2d0 * ci * k * (x(1) * cos(b * deg) + x(2) * sin(b * deg)))
  call check(error, abs(quad + pot - expected) <= 1d-10 * abs(expected), "second-order Stokes pressure")
end subroutine test_stokes_pressure

subroutine test_mean_pressure(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w = 0.9d0, h = 12d0
  real(8), parameter :: x(3) = [0.7d0, 0.3d0, -3d0]
  type(incident2) :: inc
  complex(8) :: phi1, u1(3), quad, pot, expected
  real(8) :: k

  ! random complex first-order data, the same wave twice
  phi1 = (0.3d0, -1.1d0)
  u1 = [(0.2d0, 0.9d0), (-1.3d0, 0.1d0), (0.5d0, -0.6d0)]
  call make_pair(0.8d0, 0d0, 0.8d0, 0d0, h, inc, KERNEL_DIFF)
  call field_pair(KERNEL_DIFF, 0.8d0, phi1, u1, 0.8d0, phi1, u1, inc, x, .false., rho, g, quad, pot)
  expected = -0.25d0 * rho * sum(abs(u1)**2)
  call check(error, abs(quad - expected) <= 1d-14 * abs(expected), "mean pressure of one wave")
  if (allocated(error)) return
  call check(error, abs(pot) == 0d0, "potential part on the diagonal is not zero")
  if (allocated(error)) return

  k = wave_number(w, h)
  call pair_total(KERNEL_DIFF, w, 0d0, w, 0d0, h, x, .false., quad, pot)
  expected = -rho * g * k * cosh(2d0 * k * (x(3) + h)) / (2d0 * sinh(2d0 * k * h))
  call check(error, abs(quad - expected) <= 1d-12 * abs(expected), "mean pressure of the incident wave")
  if (allocated(error)) return
  call check(error, abs(pot) == 0d0, "incident mean pressure has a potential part")
end subroutine test_mean_pressure

subroutine test_swap_symmetry(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w1 = 0.8d0, w2 = 1.1d0, h = 12d0
  real(8), parameter :: x(3) = [0.4d0, -0.9d0, -2d0]
  type(incident2) :: a, b
  complex(8) :: phi1, u1(3), phi2, u2(3), q12, p12, q21, p21
  integer :: kernel, is
  logical :: surf

  phi1 = (0.3d0, -1.1d0); phi2 = (-0.7d0, 0.4d0)
  u1 = [(0.2d0, 0.9d0), (-1.3d0, 0.1d0), (0.5d0, -0.6d0)]
  u2 = [(0.8d0, -0.2d0), (0.1d0, 0.7d0), (-0.9d0, -0.4d0)]
  do kernel = KERNEL_DIFF, KERNEL_SUM
    call make_pair(w1, 0d0, w2, 30d0, h, a, kernel)
    call make_pair(w2, 30d0, w1, 0d0, h, b, kernel)
    do is = 1, 2
      surf = is == 1
      call field_pair(kernel, w1, phi1, u1, w2, phi2, u2, a, x, surf, rho, g, q12, p12)
      call field_pair(kernel, w2, phi2, u2, w1, phi1, u1, b, x, surf, rho, g, q21, p21)
      if (kernel == KERNEL_DIFF) then
        q12 = conjg(q12)
        p12 = conjg(p12)
      end if
      call check(error, abs(q21 - q12) <= 1d-13 * abs(q12), "swap symmetry of the quadratic part")
      if (allocated(error)) return
      call check(error, abs(p21 - p12) <= 1d-13 * abs(p12), "swap symmetry of the potential part")
      if (allocated(error)) return
    end do
  end do
end subroutine test_swap_symmetry

!> Time-domain Bernoulli expressions from the first-order waves and the second-order potential (all four ordered pairs,
!> both kernels), sampled over the common period. The e^{-i Omega t} amplitudes are 2 a1 a2 (quad + pot) of the kernels.
subroutine test_bichromatic_time(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: nt = 64
  real(8), parameter :: amp(2) = [1d0, 0.7d0], ws(2) = [1d0, 0.8d0], bs(2) = [0d0, 30d0], h = 12d0
  real(8), parameter :: xs(3, 2) = reshape([0.7d0, 0.3d0, 0d0, 0.7d0, 0.3d0, -2d0], [3, 2])
  type(incident2) :: inc
  complex(8) :: phi(2), u(3, 2), e(2), c(2, 2, 2), comp(2), quad, pot, expected
  real(8) :: x(3), t, uu(3), phi2t, eta1, wzt, val, om(2), period
  integer :: ip, it, i, j, m, kn, ik
  logical :: surf

  om = [ws(1) - ws(2), ws(1) + ws(2)]
  period = 2d0 * pi / (ws(1) - ws(2))
  do ip = 1, 2
    x = xs(:, ip)
    surf = ip == 1
    do i = 1, 2
      call wave_at(ws(i), bs(i), h, x, phi(i), u(:, i))
    end do
    do kn = KERNEL_DIFF, KERNEL_SUM
      do i = 1, 2
        do j = 1, 2
          call make_pair(ws(i), bs(i), ws(j), bs(j), h, inc, kn)
          c(i, j, kn) = incident2_phi(inc, x)
        end do
      end do
    end do
    comp = (0d0, 0d0)
    do it = 0, nt - 1
      t = period * real(it, 8) / real(nt, 8)
      do i = 1, 2
        e(i) = amp(i) * exp(-ci * ws(i) * t)
      end do
      do m = 1, 3
        uu(m) = real(sum(e * u(m, :)))
      end do
      phi2t = 0d0
      do i = 1, 2
        do j = 1, 2
          phi2t = phi2t + real(-ci * (ws(i) - ws(j)) * amp(i) * amp(j) * c(i, j, KERNEL_DIFF) &
                               * exp(-ci * (ws(i) - ws(j)) * t))
          phi2t = phi2t + real(-ci * (ws(i) + ws(j)) * amp(i) * amp(j) * c(i, j, KERNEL_SUM) &
                               * exp(-ci * (ws(i) + ws(j)) * t))
        end do
      end do
      if (surf) then
        eta1 = real(sum(ci * ws / g * e * phi))
        wzt = real(sum(-ci * ws * e * u(3, :)))
        val = -(phi2t + 0.5d0 * sum(uu**2) + eta1 * wzt) / g
      else
        val = -rho * (phi2t + 0.5d0 * sum(uu**2))
      end if
      do ik = 1, 2
        comp(ik) = comp(ik) + val * exp(ci * om(ik) * t)
      end do
    end do
    ! the field holds Re[C e^{-i Omega t}], so the projection on e^{+i Omega t} gives C / 2
    comp = 2d0 * comp / real(nt, 8)
    do kn = KERNEL_DIFF, KERNEL_SUM
      call pair_total(kn, ws(1), bs(1), ws(2), bs(2), h, x, surf, quad, pot)
      expected = 2d0 * amp(1) * amp(2) * (quad + pot)
      call check(error, abs(comp(kn) - expected) <= 1d-10 * abs(expected), &
                 "time-sampled second-order field differs from the pair kernel")
      if (allocated(error)) return
    end do
  end do
end subroutine test_bichromatic_time

subroutine test_sum_continuity(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w = 0.9d0, b = 20d0, h = 12d0
  real(8), parameter :: x(3) = [0.7d0, 0.3d0, 0d0]
  complex(8) :: quad, pot, expected
  real(8) :: k

  k = wave_number(w, h)
  call pair_total(KERNEL_SUM, w, b, w + 1d-5, b, h, x, .true., quad, pot)
  expected = 0.25d0 * k * cosh(k * h) * (2d0 + cosh(2d0 * k * h)) / sinh(k * h)**3 &
             * exp(2d0 * ci * k * (x(1) * cos(b * deg) + x(2) * sin(b * deg)))
  call check(error, abs(quad + pot - expected) <= 1d-4 * abs(expected), "sum-frequency pair does not tend to the Stokes value")
end subroutine test_sum_continuity

!> Cylinder of radius 5, draft 5, with bottom.
subroutine one_cylinder(mesh)
  type(panel_mesh), intent(out) :: mesh
  call make_cylinder(5d0, 5d0, 24, 10, 5, mesh, .true.)
end subroutine one_cylinder

subroutine test_points_classified(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  call one_cylinder(mesh)
  call build_field_points(mesh, reshape([15d0, 0d0, 0d0, 0d0, -20d0, -3d0], [3, 2]), pts, surface, ok, msg)
  call check(error, ok, "classified points rejected: " // msg)
  if (allocated(error)) return
  call check(error, surface(1) .and. .not. surface(2), "surface flags")
  if (allocated(error)) return
  call check(error, all(pts%psize > 0d0) .and. all(pts%body == 1), "psize or body")
end subroutine test_points_classified

subroutine test_points_near_zero(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  call one_cylinder(mesh)
  call build_field_points(mesh, reshape([15d0, 0d0, -1d-7], [3, 1]), pts, surface, ok, msg)
  call check(error, ok .and. surface(1) .and. pts%x(3, 1) == 0d0, "z = -1e-7 is not a surface point at z = 0")
end subroutine test_points_near_zero

subroutine test_points_above(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  call one_cylinder(mesh)
  call build_field_points(mesh, reshape([15d0, 0d0, 0.1d0], [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "above the free surface") > 0, "above the free surface not rejected")
end subroutine test_points_above

subroutine test_points_on_hull(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  integer :: i
  call one_cylinder(mesh)
  do i = 1, mesh%n
    if (abs(mesh%nrm(3, i)) < 0.1d0) exit
  end do
  call build_field_points(mesh, reshape(mesh%xc(:, i), [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "hull") > 0, "point on the hull not rejected")
end subroutine test_points_on_hull

subroutine test_points_inside(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  call one_cylinder(mesh)
  call build_field_points(mesh, reshape([0d0, 0d0, -2d0], [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "inside body 1") > 0, "point inside not rejected")
end subroutine test_points_inside

subroutine test_points_waterplane(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  call one_cylinder(mesh)
  call build_field_points(mesh, reshape([0d0, 0d0, 0d0], [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "inside") > 0, "centre of the waterplane not rejected")
  if (allocated(error)) return
  call build_field_points(mesh, reshape([2d0, 1d0, 0d0], [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "inside") > 0, "waterplane point not rejected")
end subroutine test_points_waterplane

subroutine test_points_two_bodies(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, m1
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  call one_cylinder(m1)
  call make_two_bodies(m1, 40d0, mesh)
  call build_field_points(mesh, reshape([40d0, 0d0, -2d0], [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "body 2") > 0, "point inside body 2 not rejected")
  if (allocated(error)) return
  call build_field_points(mesh, reshape([20d0, 0d0, 0d0], [3, 1]), pts, surface, ok, msg)
  call check(error, ok .and. pts%body(1) >= 1, "point between the bodies rejected")
end subroutine test_points_two_bodies

end module test_qtf_field
