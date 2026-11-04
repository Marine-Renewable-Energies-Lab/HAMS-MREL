!> Geometry fixtures shared by the second-order unit tests.
module second_order_fixtures
  use SecondOrderTypes, only : panel_mesh, first_order_state
  use FirstOrderState, only : depth_factor, incident_phi, incident_grad
  implicit none
  private

  public :: make_cylinder, set_panel, maccamy_fuchs_a, maccamy_fuchs_state, maccamy_fuchs_phi, maccamy_fuchs_drift, &
            maccamy_fuchs_waterline

contains

!> Fills panel i from its nodes (already stored in mesh%xnode and listed in mesh%ncon) and a given unit normal.
!> The centroid is the vertex average, the area comes from the diagonals, and nj uses mesh%xr(:,1) as moment origin.
subroutine set_panel(mesh, i, nrm)
  type(panel_mesh), intent(inout) :: mesh
  integer, intent(in) :: i
  real(8), intent(in) :: nrm(3)
  real(8) :: v(3,4), d1(3), d2(3), c(3), x(3)
  integer :: m

  do m = 1, mesh%ncn(i)
    v(:,m) = mesh%xnode(:, mesh%ncon(m,i))
  end do
  mesh%xc(:,i) = sum(v(:,1:mesh%ncn(i)), dim=2) / real(mesh%ncn(i), 8)
  if (mesh%ncn(i) == 4) then
    d1 = v(:,3) - v(:,1)
    d2 = v(:,4) - v(:,2)
  else
    d1 = v(:,2) - v(:,1)
    d2 = v(:,3) - v(:,1)
  end if
  c = [d1(2)*d2(3) - d1(3)*d2(2), d1(3)*d2(1) - d1(1)*d2(3), d1(1)*d2(2) - d1(2)*d2(1)]
  mesh%area(i) = 0.5d0 * sqrt(sum(c**2))
  mesh%psize(i) = sqrt(mesh%area(i))
  mesh%nrm(:,i) = nrm
  x = mesh%xc(:,i) - mesh%xr(:,1)
  mesh%nj(1:3,i) = nrm
  mesh%nj(4,i) = x(2)*nrm(3) - x(3)*nrm(2)
  mesh%nj(5,i) = x(3)*nrm(1) - x(1)*nrm(3)
  mesh%nj(6,i) = x(1)*nrm(2) - x(2)*nrm(1)
end subroutine set_panel

!> Vertical cylinder about the z axis with its free surface at z = 0 (open top), single body, no symmetry.
!> Normals point into the body. Side wall: ntheta x nz panels. Bottom (optional): nr rings of ntheta panels,
!> the innermost ring being triangles.
subroutine make_cylinder(radius, draft, ntheta, nz, nr, mesh, with_bottom)
  real(8), intent(in) :: radius, draft
  integer, intent(in) :: ntheta, nz, nr
  type(panel_mesh), intent(out) :: mesh
  logical, intent(in) :: with_bottom
  real(8), parameter :: pi = 3.141592653589793d0
  integer :: nside, nnode, npan, it, iz, ir, ip, itn, i
  real(8) :: th, thm, r

  nside = ntheta * (nz + 1)
  nnode = nside
  npan = ntheta * nz
  if (with_bottom) then
    nnode = nside + ntheta * (nr - 1) + 1
    npan = npan + ntheta * nr
  end if
  mesh%n = npan
  mesh%nbase = npan
  mesh%nnode = nnode
  mesh%nbody = 1
  mesh%nsys = 1
  allocate(mesh%xc(3,npan), mesh%nrm(3,npan), mesh%area(npan), mesh%psize(npan), mesh%nj(6,npan))
  allocate(mesh%xnode(3,nnode), mesh%ncn(npan), mesh%ncon(4,npan), mesh%body(npan), mesh%base(npan), mesh%image(npan))
  allocate(mesh%xr(3,1))
  mesh%xr = 0d0
  mesh%body = 1
  mesh%image = 1
  mesh%base = [(i, i = 1, npan)]
  mesh%ncon = 0

  ! Side nodes: index it + ntheta*iz, iz = 0..nz
  do iz = 0, nz
    do it = 1, ntheta
      th = 2d0 * pi * real(it - 1, 8) / real(ntheta, 8)
      mesh%xnode(:, it + ntheta*iz) = [radius*cos(th), radius*sin(th), -draft*real(iz, 8)/real(nz, 8)]
    end do
  end do
  if (with_bottom) then
    do ir = 1, nr - 1
      r = radius * real(ir, 8) / real(nr, 8)
      do it = 1, ntheta
        th = 2d0 * pi * real(it - 1, 8) / real(ntheta, 8)
        mesh%xnode(:, nside + (ir-1)*ntheta + it) = [r*cos(th), r*sin(th), -draft]
      end do
    end do
    mesh%xnode(:, nnode) = [0d0, 0d0, -draft]
  end if

  ip = 0
  do iz = 1, nz
    do it = 1, ntheta
      itn = mod(it, ntheta) + 1
      ip = ip + 1
      mesh%ncn(ip) = 4
      mesh%ncon(1:4, ip) = [it + ntheta*(iz-1), it + ntheta*iz, itn + ntheta*iz, itn + ntheta*(iz-1)]
      thm = 2d0 * pi * (real(it - 1, 8) + 0.5d0) / real(ntheta, 8)
      call set_panel(mesh, ip, [-cos(thm), -sin(thm), 0d0])
    end do
  end do
  if (with_bottom) then
    do ir = 1, nr
      do it = 1, ntheta
        itn = mod(it, ntheta) + 1
        ip = ip + 1
        if (ir == 1) then
          mesh%ncn(ip) = 3
          mesh%ncon(1:3, ip) = [nnode, bottom_node(1, it), bottom_node(1, itn)]
        else
          mesh%ncn(ip) = 4
          mesh%ncon(1:4, ip) = [bottom_node(ir-1, it), bottom_node(ir, it), bottom_node(ir, itn), bottom_node(ir-1, itn)]
        end if
        call set_panel(mesh, ip, [0d0, 0d0, 1d0])
      end do
    end do
  end if

contains

  integer function bottom_node(ring, j)
    integer, intent(in) :: ring, j
    if (ring == nr) then
      bottom_node = j + ntheta*nz
    else
      bottom_node = nside + (ring-1)*ntheta + j
    end if
  end function bottom_node

end subroutine make_cylinder

!> Hankel function of the first kind H_n^(1)(x).
complex(8) function hankel1(n, x)
  integer, intent(in) :: n
  real(8), intent(in) :: x
  hankel1 = cmplx(bessel_jn(n, x), bessel_yn(n, x), kind=8)
end function hankel1

!> Derivative H_n^(1)'(x) = H_{n-1}(x) - (n/x) H_n(x), using H_{-1} = -H_1.
complex(8) function hankel1_prime(n, x)
  integer, intent(in) :: n
  real(8), intent(in) :: x
  if (n == 0) then
    hankel1_prime = -hankel1(1, x)
  else
    hankel1_prime = hankel1(n - 1, x) - real(n, 8) / x * hankel1(n, x)
  end if
end function hankel1_prime

!> Derivative J_n'(x) = J_{n-1}(x) - (n/x) J_n(x), using J_{-1} = -J_1.
real(8) function bessel_j_prime(n, x)
  integer, intent(in) :: n
  real(8), intent(in) :: x
  if (n == 0) then
    bessel_j_prime = -bessel_jn(1, x)
  else
    bessel_j_prime = bessel_jn(n - 1, x) - real(n, 8) / x * bessel_jn(n, x)
  end if
end function bessel_j_prime

!> Far-field amplitude a(theta) of the scattered wave of a fixed bottom-mounted cylinder of radius a_rad
!> (MacCamy-Fuchs), for an incident elevation amp e^{ik(x cos(beta) + y sin(beta))}:
!> a = -amp sqrt(2/(pi k)) e^{-i pi/4} sum_n eps_n J_n'(ka)/H_n'(ka) cos n(theta - beta).
complex(8) function maccamy_fuchs_a(theta, a_rad, k, amp, beta)
  real(8), intent(in) :: theta, a_rad, k, amp, beta
  real(8), parameter :: pi = 3.141592653589793d0
  complex(8) :: acc
  integer :: n, nmax
  real(8) :: eps

  nmax = int(k * a_rad) + 20
  acc = (0d0, 0d0)
  do n = 0, nmax
    eps = merge(1d0, 2d0, n == 0)
    acc = acc + eps * bessel_j_prime(n, k * a_rad) / hankel1_prime(n, k * a_rad) * cos(real(n, 8) * (theta - beta))
  end do
  maccamy_fuchs_a = -amp * sqrt(2d0 / (pi * k)) * exp(cmplx(0d0, -pi / 4d0, kind=8)) * acc
end function maccamy_fuchs_a

!> Analytic drift force of the same cylinder at heading 0 per unit amplitude squared: kappa int |a|^2 (1 - cos theta)
!> dtheta, with kappa = (rho g / 2) (C_g / C_p), by the trapezoid rule with n points.
real(8) function maccamy_fuchs_drift(a_rad, k, kappa, n)
  real(8), intent(in) :: a_rad, k, kappa
  integer, intent(in) :: n
  real(8), parameter :: pi = 3.141592653589793d0
  real(8) :: th
  integer :: j

  maccamy_fuchs_drift = 0d0
  do j = 1, n
    th = 2d0 * pi * real(j - 1, 8) / real(n, 8)
    maccamy_fuchs_drift = maccamy_fuchs_drift + abs(maccamy_fuchs_a(th, a_rad, k, 1d0, 0d0))**2 * (1d0 - cos(th))
  end do
  maccamy_fuchs_drift = kappa * maccamy_fuchs_drift * 2d0 * pi / real(n, 8)
end function maccamy_fuchs_drift

!> Series coefficients of the total potential on r = a: c_n = eps_n i^n 2i / (pi ka H_n'(ka)).
subroutine maccamy_fuchs_coefficients(ka, cn)
  real(8), intent(in) :: ka
  complex(8), allocatable, intent(out) :: cn(:)
  complex(8), parameter :: ci = (0d0, 1d0)
  real(8), parameter :: pi = 3.141592653589793d0
  real(8) :: eps
  integer :: n, nmax

  nmax = int(ka) + 20
  allocate(cn(0:nmax))
  do n = 0, nmax
    eps = merge(1d0, 2d0, n == 0)
    cn(n) = eps * ci**n * 2d0 * ci / (pi * ka * hankel1_prime(n, ka))
  end do
end subroutine maccamy_fuchs_coefficients

!> Analytic x force of the waterline term of the same cylinder at heading 0, divided by rho g, per unit amplitude squared:
!> (1/2) int n'_x (1/2) |zeta|^2 a dtheta with n'_x = -cos(theta), by the trapezoid rule with n points.
real(8) function maccamy_fuchs_waterline(a_rad, h, k, omega, g, n)
  real(8), intent(in) :: a_rad, h, k, omega, g
  integer, intent(in) :: n
  real(8), parameter :: pi = 3.141592653589793d0
  real(8) :: th, zeta
  integer :: j

  maccamy_fuchs_waterline = 0d0
  do j = 1, n
    th = 2d0 * pi * real(j - 1, 8) / real(n, 8)
    zeta = omega / g * abs(maccamy_fuchs_phi([a_rad * cos(th), a_rad * sin(th), 0d0], a_rad, h, k, omega, 1d0, 0d0, g))
    maccamy_fuchs_waterline = maccamy_fuchs_waterline - cos(th) * 0.25d0 * zeta**2
  end do
  maccamy_fuchs_waterline = maccamy_fuchs_waterline * a_rad * 2d0 * pi / real(n, 8)
end function maccamy_fuchs_waterline

!> Exact total potential on r = a at the azimuth and z of the point x(3), which need not lie exactly on r = a.
complex(8) function maccamy_fuchs_phi(x, a_rad, h, k, omega, amp, beta, g)
  real(8), intent(in) :: x(3), a_rad, h, k, omega, amp, beta, g
  complex(8), allocatable :: cn(:)
  complex(8) :: ssum
  real(8) :: tr, zf, dzf
  integer :: n

  call maccamy_fuchs_coefficients(k * a_rad, cn)
  tr = atan2(x(2), x(1)) - beta
  ssum = (0d0, 0d0)
  do n = 0, ubound(cn, 1)
    ssum = ssum + cn(n) * cos(real(n, 8) * tr)
  end do
  call depth_factor(k, x(3), h, zf, dzf)
  maccamy_fuchs_phi = cmplx(0d0, -g * amp / omega, kind=8) * zf * ssum
end function maccamy_fuchs_phi

!> Analytic first-order state of a fixed bottom-mounted cylinder (axis z, radius a_rad), no BEM solve. The total
!> potential on r = a is evaluated at each panel centroid's azimuth and z; phi_i and grad_i are the incident wave
!> at the centroid, phi_p = phi_total - phi_i, dphi_p_dn = -grad_i . n (fixed body), and grad_p = grad_total - grad_i,
!> where the total gradient on r = a has only azimuthal and vertical parts. h < 0 is deep water.
subroutine maccamy_fuchs_state(mesh, a_rad, h, k, omega, amp, beta, g, rho, s)
  type(panel_mesh), intent(in) :: mesh
  real(8), intent(in) :: a_rad, h, k, omega, amp, beta, g, rho
  type(first_order_state), intent(out) :: s
  complex(8), allocatable :: cn(:)
  complex(8) :: pref, ssum, dssum, grad_t(3)
  real(8) :: th, zf, dzf, tr, phi_dir(3)
  integer :: n, nmax, i

  s%omega = omega
  s%k = k
  s%nu = omega**2 / g
  s%h = h
  s%amp = amp
  s%beta = beta
  s%g = g
  s%rho = rho
  allocate(s%phi_i(mesh%n), s%phi_p(mesh%n), s%dphi_p_dn(mesh%n), s%dphi_dn(mesh%n))
  allocate(s%grad_i(3, mesh%n), s%grad_p(3, mesh%n), s%xi(6, mesh%nbody))
  s%xi = (0d0, 0d0)
  s%dphi_dn = (0d0, 0d0)

  call maccamy_fuchs_coefficients(k * a_rad, cn)
  nmax = ubound(cn, 1)
  pref = cmplx(0d0, -g * amp / omega, kind=8)

  do i = 1, mesh%n
    th = atan2(mesh%xc(2, i), mesh%xc(1, i))
    tr = th - beta
    ssum = (0d0, 0d0)
    dssum = (0d0, 0d0)
    do n = 0, nmax
      ssum = ssum + cn(n) * cos(real(n, 8) * tr)
      dssum = dssum - cn(n) * real(n, 8) * sin(real(n, 8) * tr)
    end do
    call depth_factor(k, mesh%xc(3, i), h, zf, dzf)
    phi_dir = [-sin(th), cos(th), 0d0]
    grad_t = pref * zf * dssum / a_rad * phi_dir
    grad_t(3) = pref * dzf * ssum

    s%phi_i(i) = incident_phi(mesh%xc(:, i), omega, k, h, amp, beta, g)
    call incident_grad(mesh%xc(:, i), omega, k, h, amp, beta, g, s%grad_i(:, i))
    s%phi_p(i) = pref * zf * ssum - s%phi_i(i)
    s%dphi_p_dn(i) = -sum(s%grad_i(:, i) * mesh%nrm(:, i))
    s%grad_p(:, i) = grad_t - s%grad_i(:, i)
  end do
end subroutine maccamy_fuchs_state

end module second_order_fixtures
