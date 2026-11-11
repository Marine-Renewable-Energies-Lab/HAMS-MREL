!> Unit tests of QTFFreeSurface and of the free-surface forcing kernel of SecondOrderIncident (Phase 6a, Task 1).
module test_qtf_free_surface
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderIncident, only : incident2, incident2_pair, incident2_phi, incident2_grad, free_surface_forcing
  use FirstOrderState, only : incident_phi, incident_grad
  use NearFieldDrift, only : KERNEL_DIFF
  use QTFFreeSurface, only : fs_point_field, fs_forcing_diff
  implicit none
  private

  public :: collect_qtf_free_surface

  real(8), parameter :: g = 9.81d0, pi = 3.14159265358979323846d0, deg = pi / 180d0

contains

subroutine collect_qtf_free_surface(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("fs_forcing_matches_incident2", test_matches_incident2), &
    new_unittest("fs_forcing_swap", test_swap), &
    new_unittest("fs_forcing_incident_only_zero", test_incident_only_zero) &
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

end module test_qtf_free_surface
