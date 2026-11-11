!> Unit tests of SecondOrderIncident: the second-order incident potential of a pair of unit-amplitude waves.
module test_second_order_incident
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderIncident, only : incident2, incident2_pair, incident2_phi, incident2_grad
  use FirstOrderState, only : incident_phi, incident_grad
  implicit none
  private

  public :: collect_so_incident

  real(8), parameter :: g = 9.81d0, pi = 3.14159265358979323846d0, deg = pi / 180d0
  complex(8), parameter :: ci = (0d0, 1d0)

contains

subroutine collect_so_incident(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("incident2_laplace", test_laplace), &
    new_unittest("incident2_bottom", test_bottom), &
    new_unittest("incident2_free_surface", test_free_surface), &
    new_unittest("incident2_values", test_values), &
    new_unittest("incident2_hermitian", test_hermitian), &
    new_unittest("incident2_diagonal", test_diagonal) &
    ]
end subroutine collect_so_incident

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

!> The pair of waves (omega, heading in degrees) in depth h.
subroutine make_pair(w1, b1, w2, b2, h, inc)
  real(8), intent(in) :: w1, b1, w2, b2, h
  type(incident2), intent(out) :: inc
  call incident2_pair(w1, wave_number(w1, h), b1 * deg, w2, wave_number(w2, h), b2 * deg, h, g, inc)
end subroutine make_pair

subroutine test_laplace(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: depth = 100d0
  real(8), parameter :: zs(3) = [0d0, -0.25d0 * depth, -0.5d0 * depth]
  type(incident2) :: inc
  complex(8) :: lap
  real(8) :: x(3), d, hs
  integer :: ih, ip, m

  do ih = 1, 2
    hs = merge(depth, -1d0, ih == 1)
    call make_pair(1.2d0, 0d0, 0.5d0, 30d0, hs, inc)
    d = 1d-3 / inc%kmag
    do ip = 1, 5
      x = [1.3d0 * ip - 2d0, 0.7d0 - 0.4d0 * ip, zs(mod(ip - 1, 3) + 1)]
      lap = -6d0 * incident2_phi(inc, x)
      do m = 1, 3
        x(m) = x(m) + d
        lap = lap + incident2_phi(inc, x)
        x(m) = x(m) - 2d0 * d
        lap = lap + incident2_phi(inc, x)
        x(m) = x(m) + d
      end do
      lap = lap / d**2
      call check(error, abs(lap) <= 1d-6 * abs(inc%coef) * inc%kmag**2, "Laplacian of the second-order potential is not zero")
      if (allocated(error)) return
    end do
  end do
end subroutine test_laplace

subroutine test_bottom(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: depth = 100d0
  type(incident2) :: inc
  complex(8) :: gr(3)
  call make_pair(1.2d0, 0d0, 0.5d0, 30d0, depth, inc)
  call incident2_grad(inc, [2d0, -1d0, -depth], gr)
  call check(error, abs(gr(3)) <= 1d-12 * abs(inc%coef) * inc%kmag, "vertical velocity at the bottom is not zero")
end subroutine test_bottom

!> The exact nonlinear forcing -d/dt |grad phi|^2 + (1/g) phi_t d/dz (phi_tt + g phi_z) is sampled in time from the two
!> analytic first-order waves. Its e^{-i Omega t} amplitude is 2 a_i a_j (-Omega^2 phi2 + g phi2_z) of the second-order
!> potential: the sampled forcing holds both the (i, j) and (j, i) terms. The time factor is e^{-i omega t}, as in incident_phi.
subroutine test_free_surface(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: depth = 100d0, amp(2) = [1d0, 0.7d0], dw = 0.1d0
  integer, parameter :: nt = 64, ncase = 2
  real(8), parameter :: ws(2, ncase) = reshape([1d0, 0.8d0, 1.2d0, 0.5d0], [2, ncase])
  real(8), parameter :: bs(2, ncase) = reshape([0d0, 0d0, 0d0, 30d0], [2, ncase])
  real(8), parameter :: x(3) = [3d0, 2d0, 0d0]
  type(incident2) :: inc
  complex(8) :: phi(2), u(3, 2), s(2), comp, expected, gr(3), e(2), phi2
  real(8) :: t, force, omega_d, hs, ut(3), uu(3), phit, ssum, k(2)
  integer :: ic, ih, it, i, m

  do ic = 1, ncase
    do ih = 1, 2
      hs = merge(depth, -1d0, ih == 1)
      call make_pair(ws(1, ic), bs(1, ic), ws(2, ic), bs(2, ic), hs, inc)
      do i = 1, 2
        k(i) = wave_number(ws(i, ic), hs)
        phi(i) = incident_phi(x, ws(i, ic), k(i), hs, 1d0, bs(i, ic) * deg, g)
        call incident_grad(x, ws(i, ic), k(i), hs, 1d0, bs(i, ic) * deg, g, u(:, i))
        ! first-order waves have phi_zz = k^2 phi
        s(i) = -ws(i, ic)**2 * u(3, i) + g * k(i)**2 * phi(i)
      end do
      omega_d = ws(1, ic) - ws(2, ic)
      comp = (0d0, 0d0)
      do it = 0, nt - 1
        t = 2d0 * pi / dw * real(it, 8) / real(nt, 8)
        do i = 1, 2
          e(i) = amp(i) * exp(-ci * ws(i, ic) * t)
        end do
        do m = 1, 3
          uu(m) = real(sum(e * u(m, :)))
          ut(m) = real(sum(-ci * ws(:, ic) * e * u(m, :)))
        end do
        phit = real(sum(-ci * ws(:, ic) * e * phi))
        ssum = real(sum(e * s))
        force = -2d0 * sum(uu * ut) + phit * ssum / g
        comp = comp + force * exp(ci * omega_d * t)
      end do
      ! F contains Re[C e^{-i Omega t}], so the projection on e^{+i Omega t} gives C / 2
      comp = 2d0 * comp / real(nt, 8)
      phi2 = incident2_phi(inc, x)
      call incident2_grad(inc, x, gr)
      expected = 2d0 * amp(1) * amp(2) * (-omega_d**2 * phi2 + g * gr(3))
      call check(error, abs(comp - expected) <= 1d-10 * abs(expected), &
                 "time-sampled forcing differs from the free-surface condition of the second-order potential")
      if (allocated(error)) return
    end do
  end do
end subroutine test_free_surface

subroutine test_values(error)
  type(error_type), allocatable, intent(out) :: error
  type(incident2) :: inc

  call make_pair(1d0, 0d0, 0.8d0, 0d0, -1d0, inc)
  call check(error, abs(abs(inc%coef) - 0.5d0) <= 1d-10, "deep-water coefficient, 1.0 and 0.8")
  if (allocated(error)) return
  call make_pair(1d0, 0d0, 0.8d0, 0d0, 100d0, inc)
  call check(error, abs(abs(inc%coef) - 0.500733d0) <= 1d-5, "h = 100 coefficient, 1.0 and 0.8")
  if (allocated(error)) return
  call make_pair(1.2d0, 0d0, 0.5d0, 0d0, -1d0, inc)
  call check(error, abs(abs(inc%coef) - 0.6d0) <= 1d-10, "deep-water coefficient, 1.2 and 0.5")
  if (allocated(error)) return
  call make_pair(1.2d0, 0d0, 0.5d0, 0d0, 100d0, inc)
  call check(error, abs(abs(inc%coef) - 0.604943d0) <= 1d-5, "h = 100 coefficient, 1.2 and 0.5")
end subroutine test_values

subroutine test_hermitian(error)
  type(error_type), allocatable, intent(out) :: error
  type(incident2) :: a, b
  integer :: ih
  real(8) :: hs

  do ih = 1, 2
    hs = merge(100d0, -1d0, ih == 1)
    call make_pair(1.2d0, 10d0, 0.5d0, 30d0, hs, a)
    call make_pair(0.5d0, 30d0, 1.2d0, 10d0, hs, b)
    call check(error, abs(b%coef - conjg(a%coef)) <= 1d-13 * abs(a%coef), "swapping the pair does not conjugate the coefficient")
    if (allocated(error)) return
    call check(error, maxval(abs(b%kvec + a%kvec)) <= 1d-14 .and. abs(b%omega + a%omega) <= 1d-14, &
               "swapping the pair does not reverse K and Omega")
    if (allocated(error)) return
  end do
end subroutine test_hermitian

subroutine test_diagonal(error)
  type(error_type), allocatable, intent(out) :: error
  type(incident2) :: inc
  integer :: ih
  real(8) :: hs

  do ih = 1, 2
    hs = merge(10d0, -1d0, ih == 1)
    call make_pair(0.9d0, 20d0, 0.9d0, 20d0, hs, inc)
    call check(error, inc%coef == (0d0, 0d0), "coefficient on the exact diagonal is not zero")
    if (allocated(error)) return
    call make_pair(0.9d0, 0d0, 0.9d0, 40d0, hs, inc)
    call check(error, abs(inc%coef) <= huge(1d0), "NaN or Inf at equal frequency")
    if (allocated(error)) return
    call check(error, inc%kmag > 0d0, "equal frequency, other headings: K vanishes")
    if (allocated(error)) return
    ! At Omega = 0 the first term of Q has the factor Omega, and in the second f and g are in quadrature (phi_t against
    ! d/dz (phi_tt + g phi_z) of waves with the same phase at the origin), so Q and the coefficient vanish: no set-down
    ! from equal frequencies, in any depth.
    call check(error, abs(inc%coef) <= 1d-10, "equal frequency, other headings: coefficient is not zero")
    if (allocated(error)) return
  end do
end subroutine test_diagonal

end module test_second_order_incident
