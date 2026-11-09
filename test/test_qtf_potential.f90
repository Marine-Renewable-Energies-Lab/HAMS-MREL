!> Unit tests of QTFPotential: the body-forcing terms of the potential part of the difference-frequency QTF.
!>
!> The second-derivative term int psi n.[(X.grad)u] dS is computed in the form that Stokes' theorem gives, with first
!> derivatives only. These tests compare it with the direct form, which uses the analytic Hessian of an incident wave, on
!> analytic fields (no solver). psi is a harmonic function, e^{a z} cos(a (x cos g_k + y sin g_k) + c_k), a different one
!> for each mode column k, so that a mix-up of columns shows. Column 1 is e^{0.05 z} cos(0.05 x + 0.3).
module test_qtf_potential
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderTypes, only : panel_mesh, surface_geom, first_order_state
  use FirstOrderState, only : incident_phi, incident_grad
  use SurfaceGeometry, only : build_surface_geom, DEFAULT_SHARP_DEG
  use SurfaceVelocity, only : velocity_differencing
  use NearFieldDrift, only : pair_product, rigid_displacement, waterline_extrapolate, waterline_point, KERNEL_DIFF, KERNEL_SUM
  use QTFPotential, only : hessian_term_converted, psi_waterline, potential_pair, psi_surface_gradient, &
                           waterline_velocity, NPOT_TERMS, TP_FK, TP_DIFFRACTION, TP_B1, TP_B2, TP_B3
  use second_order_fixtures, only : make_cylinder, set_panel
  use SecondOrderIncident, only : incident2, incident2_pair, incident2_grad
  implicit none
  private

  public :: collect_qtf_potential

  real(8), parameter :: g = 9.81d0, depth = 100d0, pi = 3.14159265358979323846d0
  integer, parameter :: nmode = 6
  real(8), parameter :: a_psi = 0.05d0, rho = 1025d0
  ! The first-order states of the pair: frequency and heading
  real(8), parameter :: w1 = 0.8d0, b1 = 0d0, w2 = 0.7d0, b2 = 0.3d0
  complex(8), parameter :: ci = (0d0, 1d0)

  ! Cylinder radius 10, draft 10: 448, 1792 and 7168 panels
  integer, parameter :: cyl_ntheta(3) = [32, 64, 128], cyl_nz(3) = [8, 16, 32], cyl_nr(3) = [6, 12, 24]
  ! Sphere radius 5, centre depth 20: 128, 512 and 2048 panels
  integer, parameter :: sph_nlat(3) = [8, 16, 32]

contains

subroutine collect_qtf_potential(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("stokes_identity_smooth", test_identity_smooth), &
    new_unittest("stokes_identity_submerged", test_identity_submerged), &
    new_unittest("stokes_identity_differenced", test_identity_differenced), &
    new_unittest("stokes_identity_differenced_velocity", test_identity_differenced_velocity), &
    new_unittest("stokes_identity_waterline_velocity", test_identity_waterline_velocity), &
    new_unittest("stokes_kernel", test_kernel), &
    new_unittest("pot_fixed_body", test_pot_fixed_body), &
    new_unittest("pot_hermitian", test_pot_hermitian), &
    new_unittest("pot_equal_frequency", test_pot_equal_frequency), &
    new_unittest("pot_b3_direct", test_pot_b3_direct), &
    new_unittest("pot_b1_direct", test_pot_b1_direct), &
    new_unittest("pot_deep_scaling", test_pot_deep_scaling), &
    new_unittest("pot_fk_uniform", test_pot_fk_uniform), &
    new_unittest("pot_psi_missing", test_pot_psi_missing), &
    new_unittest("pot_signs_disturbance", test_pot_signs_disturbance) &
    ]
end subroutine collect_qtf_potential

! ---------------------------------------------------------------------------------------------- helpers

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

!> Rigid motions of the pair, one body. Pair 1 is X = xi + alpha x x with xi = (0.3, -0.2, 0.5), alpha = (0.01, 0.02, -0.015).
subroutine motions(x1, x2)
  complex(8), intent(out) :: x1(6, 1), x2(6, 1)
  x1(:, 1) = [(0.3d0, 0d0), (-0.2d0, 0d0), (0.5d0, 0d0), (0.01d0, 0d0), (0.02d0, 0d0), (-0.015d0, 0d0)]
  x2(:, 1) = [(0.1d0, 0.2d0), (0d0, 0.3d0), (-0.4d0, 0d0), (0d0, -0.02d0), (0.01d0, 0d0), (0.03d0, 0.01d0)]
end subroutine motions

!> psi of mode column k at x and its gradient.
subroutine psi_field(k, x, psi, grad)
  integer, intent(in) :: k
  real(8), intent(in) :: x(3)
  complex(8), intent(out) :: psi, grad(3)
  real(8) :: gam, c, arg, e
  gam = 0.5d0 * real(k - 1, 8)
  c = 0.3d0 + 0.4d0 * real(k - 1, 8)
  arg = a_psi * (x(1) * cos(gam) + x(2) * sin(gam)) + c
  e = exp(a_psi * x(3))
  psi = e * cos(arg)
  grad = a_psi * e * [-sin(arg) * cos(gam), -sin(arg) * sin(gam), cos(arg)]
end subroutine psi_field

!> psi and its analytic gradient at the panel centroids, for all mode columns.
subroutine psi_panels(mesh, psi, grad_psi)
  type(panel_mesh), intent(in) :: mesh
  complex(8), allocatable, intent(out) :: psi(:,:), grad_psi(:,:,:)
  integer :: i, k
  allocate(psi(mesh%n, nmode), grad_psi(3, mesh%n, nmode))
  do k = 1, nmode
    do i = 1, mesh%n
      call psi_field(k, mesh%xc(:, i), psi(i, k), grad_psi(:, i, k))
    end do
  end do
end subroutine psi_panels

!> grad psi from velocity_differencing: the tangential part from the panel values, the normal part exact.
subroutine psi_differenced(mesh, geom, psi, grad_psi)
  type(panel_mesh), intent(in) :: mesh
  type(surface_geom), intent(in) :: geom
  complex(8), intent(in) :: psi(:,:)
  complex(8), intent(inout) :: grad_psi(:,:,:)
  type(first_order_state) :: s
  integer :: k, i, nfb
  allocate(s%phi_p(mesh%n), s%dphi_p_dn(mesh%n), s%grad_p(3, mesh%n))
  do k = 1, nmode
    s%phi_p = psi(:, k)
    do i = 1, mesh%n
      s%dphi_p_dn(i) = sum(mesh%nrm(:, i) * grad_psi(:, i, k))
    end do
    call velocity_differencing(mesh, geom, s, nfb)
    grad_psi(:, :, k) = s%grad_p
  end do
end subroutine psi_differenced

!> Total velocity of the incident wave (omega, heading beta) at the panel centroids.
subroutine velocity(mesh, omega, beta, u)
  type(panel_mesh), intent(in) :: mesh
  real(8), intent(in) :: omega, beta
  complex(8), allocatable, intent(out) :: u(:,:)
  integer :: i
  allocate(u(3, mesh%n))
  do i = 1, mesh%n
    call incident_grad(mesh%xc(:, i), omega, wave_number(omega, depth), depth, 1d0, beta, g, u(:, i))
  end do
end subroutine velocity

!> As velocity, but as production gets it: the incident wave put in phi_p with its exact normal derivative, and the
!> gradient from velocity_differencing (one-sided on the waterline and bottom-edge panels).
subroutine velocity_differenced(mesh, geom, omega, beta, u)
  type(panel_mesh), intent(in) :: mesh
  type(surface_geom), intent(in) :: geom
  real(8), intent(in) :: omega, beta
  complex(8), allocatable, intent(out) :: u(:,:)
  type(first_order_state) :: s
  integer :: i, nfb
  real(8) :: k
  call velocity(mesh, omega, beta, u)
  k = wave_number(omega, depth)
  allocate(s%phi_p(mesh%n), s%dphi_p_dn(mesh%n), s%grad_p(3, mesh%n))
  do i = 1, mesh%n
    s%phi_p(i) = incident_phi(mesh%xc(:, i), omega, k, depth, 1d0, beta, g)
    s%dphi_p_dn(i) = sum(mesh%nrm(:, i) * u(:, i))
  end do
  call velocity_differencing(mesh, geom, s, nfb)
  u = s%grad_p
end subroutine velocity_differenced

!> As velocity_differenced, and u at the waterline segment midpoints from waterline_velocity, as potential_pair builds it.
!> The whole wave is the disturbance here (amp = 0, so nothing is taken as exact incident wave): phi at the waterline comes
!> from waterline_extrapolate of the panel values, u_z from the free-surface condition u_z = (omega^2/g) phi, the tangential
!> part from differencing phi along the waterline. The normal velocity at the waterline is the exact n.u of the wave, as
!> the body condition gives it in production from the rigid motion.
subroutine velocity_differenced_wl(mesh, geom, omega, beta, u, u_wl)
  type(panel_mesh), intent(in) :: mesh
  type(surface_geom), intent(in) :: geom
  real(8), intent(in) :: omega, beta
  complex(8), allocatable, intent(out) :: u(:,:), u_wl(:,:)
  type(first_order_state) :: s
  complex(8), allocatable :: phi_wl(:), vn_wl(:)
  complex(8) :: ue(3)
  integer :: i, nfb, iw
  real(8) :: k
  call velocity(mesh, omega, beta, u)
  k = wave_number(omega, depth)
  s%omega = omega
  s%k = k
  s%h = depth
  s%amp = 0d0
  s%beta = beta
  s%g = g
  s%nu = omega**2 / g
  allocate(s%phi_p(mesh%n), s%dphi_p_dn(mesh%n), s%grad_p(3, mesh%n))
  do i = 1, mesh%n
    s%phi_p(i) = incident_phi(mesh%xc(:, i), omega, k, depth, 1d0, beta, g)
    s%dphi_p_dn(i) = sum(mesh%nrm(:, i) * u(:, i))
  end do
  call velocity_differencing(mesh, geom, s, nfb)
  u = s%grad_p
  allocate(phi_wl(geom%nwl), vn_wl(geom%nwl), u_wl(3, geom%nwl))
  call waterline_extrapolate(mesh, geom, s, phi_wl)
  do iw = 1, geom%nwl
    call incident_grad(waterline_point(geom, iw), omega, k, depth, 1d0, beta, g, ue)
    vn_wl(iw) = sum(mesh%nrm(:, geom%wl_panel(iw)) * ue)
  end do
  call waterline_velocity(mesh, geom, s, phi_wl, vn_wl, u_wl)
end subroutine velocity_differenced_wl

!> n . (grad u) along each axis j: (n.H)_j, with H the analytic Hessian of the incident wave at x.
!> phi = C Z(z) e^{i (kx x + ky y)}, so H_xj = i kx phi_j, H_yj = i ky phi_j, H_zz = k^2 phi.
function normal_hessian(x, nrm, omega, beta) result(nh)
  real(8), intent(in) :: x(3), nrm(3), omega, beta
  complex(8) :: nh(3)
  complex(8) :: phi, gr(3), hs(3, 3)
  real(8) :: k, kx, ky
  k = wave_number(omega, depth)
  kx = k * cos(beta)
  ky = k * sin(beta)
  phi = incident_phi(x, omega, k, depth, 1d0, beta, g)
  call incident_grad(x, omega, k, depth, 1d0, beta, g, gr)
  hs(1, :) = ci * kx * gr
  hs(2, :) = ci * ky * gr
  hs(3, :) = [ci * kx * gr(3), ci * ky * gr(3), k * k * phi]
  nh = matmul(nrm, hs)
end function normal_hessian

!> The direct form: int psi_k n.[(X.grad)u] dS = int psi_k sum_j P(X_j, (n.H)_j) dS, with P the kernel pair.
subroutine direct_form(mesh, psi, x1, x2, kernel, out)
  type(panel_mesh), intent(in) :: mesh
  complex(8), intent(in) :: psi(:,:), x1(:,:), x2(:,:)
  integer, intent(in) :: kernel
  complex(8), intent(out) :: out(:)
  complex(8) :: d1(3), d2(3), h1(3), h2(3), p
  integer :: i, j, b
  out = (0d0, 0d0)
  do i = 1, mesh%n
    b = mesh%body(i)
    d1 = rigid_displacement(x1(:, b), mesh%xc(:, i) - mesh%xr(:, b))
    d2 = rigid_displacement(x2(:, b), mesh%xc(:, i) - mesh%xr(:, b))
    h1 = normal_hessian(mesh%xc(:, i), mesh%nrm(:, i), w1, b1)
    h2 = normal_hessian(mesh%xc(:, i), mesh%nrm(:, i), w2, b2)
    p = (0d0, 0d0)
    do j = 1, 3
      p = p + pair_product(kernel, d1(j), h1(j), d2(j), h2(j))
    end do
    out = out + psi(i, :) * p * mesh%area(i)
  end do
end subroutine direct_form

subroutine geometry(mesh, geom, error)
  type(panel_mesh), intent(in) :: mesh
  type(surface_geom), intent(out) :: geom
  type(error_type), allocatable, intent(out) :: error
  logical :: ok
  character(len=:), allocatable :: msg
  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "build_surface_geom stopped: " // msg)
end subroutine geometry

!> Closed sphere of radius r about c, latitude-longitude panels: nlat bands of 2 nlat panels, triangles at the poles.
!> Normals are the flat-panel normals, pointing into the body. The centroids are area centroids, not the vertex average that
!> HAMS (CalPanelCentre) and set_panel use. On these trapezoids the vertex average lies O(h) off the area centroid, which
!> makes the one-point rule O(h^2) with a large constant. The closed-surface integral nearly cancels, so this constant
!> dominates the relative difference: measured 2.4e-1, 7.3e-2, 1.9e-2 at 128, 512, 2048 panels with the vertex average
!> (4.8e-3 at 8192), against 8.1e-3, 2.2e-3, 5.5e-4 with the area centroid. The test checks the identity, so it takes the
!> area centroid.
subroutine make_sphere(r, c, nlat, mesh, vertex_average)
  real(8), intent(in) :: r, c(3)
  integer, intent(in) :: nlat
  type(panel_mesh), intent(out) :: mesh
  logical, intent(in), optional :: vertex_average   ! keep the vertex-average centroid of set_panel, as HAMS does
  integer :: nlon, nnode, npan, il, jl, jn, ip, i
  real(8) :: th, ph, v(3, 4), d1(3), d2(3), nv(3), xc(3)

  nlon = 2 * nlat
  nnode = 2 + (nlat - 1) * nlon
  npan = nlat * nlon
  mesh%n = npan; mesh%nbase = npan; mesh%nnode = nnode; mesh%nbody = 1; mesh%nsys = 1
  allocate(mesh%xc(3,npan), mesh%nrm(3,npan), mesh%area(npan), mesh%psize(npan), mesh%nj(6,npan))
  allocate(mesh%xnode(3,nnode), mesh%ncn(npan), mesh%ncon(4,npan), mesh%body(npan), mesh%base(npan), mesh%image(npan))
  allocate(mesh%xr(3,1))
  mesh%xr(:, 1) = c
  mesh%body = 1
  mesh%image = 1
  mesh%base = [(i, i = 1, npan)]
  mesh%ncon = 0

  mesh%xnode(:, 1) = c + [0d0, 0d0, r]
  mesh%xnode(:, nnode) = c - [0d0, 0d0, r]
  do il = 1, nlat - 1
    th = pi * real(il, 8) / real(nlat, 8)
    do jl = 1, nlon
      ph = 2d0 * pi * real(jl - 1, 8) / real(nlon, 8)
      mesh%xnode(:, ring_node(il, jl)) = c + r * [sin(th) * cos(ph), sin(th) * sin(ph), cos(th)]
    end do
  end do

  ip = 0
  do il = 1, nlat
    do jl = 1, nlon
      jn = mod(jl, nlon) + 1
      ip = ip + 1
      if (il == 1) then
        mesh%ncn(ip) = 3
        mesh%ncon(1:3, ip) = [1, ring_node(1, jl), ring_node(1, jn)]
      else if (il == nlat) then
        mesh%ncn(ip) = 3
        mesh%ncon(1:3, ip) = [nnode, ring_node(nlat - 1, jn), ring_node(nlat - 1, jl)]
      else
        mesh%ncn(ip) = 4
        mesh%ncon(1:4, ip) = [ring_node(il - 1, jl), ring_node(il, jl), ring_node(il, jn), ring_node(il - 1, jn)]
      end if
      do i = 1, mesh%ncn(ip)
        v(:, i) = mesh%xnode(:, mesh%ncon(i, ip))
      end do
      xc = sum(v(:, 1:mesh%ncn(ip)), dim=2) / real(mesh%ncn(ip), 8)
      if (mesh%ncn(ip) == 4) then
        d1 = v(:, 3) - v(:, 1)
        d2 = v(:, 4) - v(:, 2)
      else
        d1 = v(:, 2) - v(:, 1)
        d2 = v(:, 3) - v(:, 1)
      end if
      nv = [d1(2)*d2(3) - d1(3)*d2(2), d1(3)*d2(1) - d1(1)*d2(3), d1(1)*d2(2) - d1(2)*d2(1)]
      nv = nv / sqrt(sum(nv**2))
      if (sum(nv * (xc - c)) > 0d0) nv = -nv
      call set_panel(mesh, ip, nv)
      if (.not. present(vertex_average)) then
        mesh%xc(:, ip) = area_centroid(v, mesh%ncn(ip))
      else if (.not. vertex_average) then
        mesh%xc(:, ip) = area_centroid(v, mesh%ncn(ip))
      end if
    end do
  end do

contains

  !> Area centroid of a planar polygon of 3 or 4 nodes, from its triangles about node 1.
  function area_centroid(v, nc) result(xa)
    real(8), intent(in) :: v(3, 4)
    integer, intent(in) :: nc
    real(8) :: xa(3), cr(3), e1(3), e2(3), at, asum
    integer :: m
    xa = 0d0
    asum = 0d0
    do m = 2, nc - 1
      e1 = v(:, m) - v(:, 1)
      e2 = v(:, m + 1) - v(:, 1)
      cr = [e1(2)*e2(3) - e1(3)*e2(2), e1(3)*e2(1) - e1(1)*e2(3), e1(1)*e2(2) - e1(2)*e2(1)]
      at = 0.5d0 * sqrt(sum(cr**2))
      xa = xa + at * (v(:, 1) + v(:, m) + v(:, m + 1)) / 3d0
      asum = asum + at
    end do
    xa = xa / asum
  end function area_centroid

  integer function ring_node(ring, j)
    integer, intent(in) :: ring, j
    ring_node = 1 + (ring - 1) * nlon + j
  end function ring_node

end subroutine make_sphere

!> Relative difference of the converted and direct forms for pair (1, 2) and for pair (1, 1), the larger of the two.
!> differenced: grad psi from velocity_differencing instead of the analytic gradient.
!> differenced_u (default .false.): u from velocity_differencing too (velocity_differenced), the direct form stays exact.
!> waterline_u (default .false.): with differenced_u, also u at the waterline from waterline_velocity (free-surface
!> condition, body condition and differencing along the waterline) instead of the extrapolation inside the routine.
real(8) function identity_difference(mesh, geom, differenced, differenced_u, waterline_u) result(rel)
  type(panel_mesh), intent(in) :: mesh
  type(surface_geom), intent(in) :: geom
  logical, intent(in) :: differenced
  logical, intent(in), optional :: differenced_u, waterline_u
  complex(8), allocatable :: psi(:,:), grad_psi(:,:,:), psi_wl(:,:), u1(:,:), u2(:,:), uw1(:,:), uw2(:,:)
  complex(8) :: x1(6, 1), x2(6, 1), conv(nmode), dir(nmode)
  logical :: du, dw

  du = .false.
  if (present(differenced_u)) du = differenced_u
  dw = .false.
  if (present(waterline_u)) dw = waterline_u
  call psi_panels(mesh, psi, grad_psi)
  if (differenced) call psi_differenced(mesh, geom, psi, grad_psi)
  allocate(psi_wl(geom%nwl, nmode))
  call psi_waterline(mesh, geom, psi, grad_psi, psi_wl)
  call motions(x1, x2)
  call get_velocity(w1, b1, u1, uw1)
  call get_velocity(w2, b2, u2, uw2)

  if (dw) then
    call hessian_term_converted(mesh, geom, psi, grad_psi, psi_wl, x1, u1, x2, u2, KERNEL_DIFF, conv, uw1, uw2)
  else
    call hessian_term_converted(mesh, geom, psi, grad_psi, psi_wl, x1, u1, x2, u2, KERNEL_DIFF, conv)
  end if
  call direct_form(mesh, psi, x1, x2, KERNEL_DIFF, dir)
  rel = norm2(abs(conv - dir)) / norm2(abs(dir))

  ! Pair (1, 1): the brief's case on its own, with u2 = u1 at w1
  call get_velocity(w1, b1, u2, uw2)
  if (dw) then
    call hessian_term_converted(mesh, geom, psi, grad_psi, psi_wl, x1, u1, x1, u2, KERNEL_DIFF, conv, uw1, uw2)
  else
    call hessian_term_converted(mesh, geom, psi, grad_psi, psi_wl, x1, u1, x1, u2, KERNEL_DIFF, conv)
  end if
  call direct_form_same(dir)
  rel = max(rel, norm2(abs(conv - dir)) / norm2(abs(dir)))

contains

  subroutine get_velocity(omega, beta, u, uw)
    real(8), intent(in) :: omega, beta
    complex(8), allocatable, intent(out) :: u(:,:), uw(:,:)
    if (dw) then
      call velocity_differenced_wl(mesh, geom, omega, beta, u, uw)
    else if (du) then
      call velocity_differenced(mesh, geom, omega, beta, u)
    else
      call velocity(mesh, omega, beta, u)
    end if
  end subroutine get_velocity

  subroutine direct_form_same(out)
    complex(8), intent(out) :: out(:)
    complex(8) :: d(3), h(3), p
    integer :: i, j
    out = (0d0, 0d0)
    do i = 1, mesh%n
      d = rigid_displacement(x1(:, 1), mesh%xc(:, i) - mesh%xr(:, 1))
      h = normal_hessian(mesh%xc(:, i), mesh%nrm(:, i), w1, b1)
      p = (0d0, 0d0)
      do j = 1, 3
        p = p + pair_product(KERNEL_DIFF, d(j), h(j), d(j), h(j))
      end do
      out = out + psi(i, :) * p * mesh%area(i)
    end do
  end subroutine direct_form_same

end function identity_difference

! ---------------------------------------------------------------------------------------------- tests

!> Truncated cylinder with its waterline, analytic grad psi: the converted form converges to the direct form. Both forms are
!> integrals over the same flat panels, on which Stokes' theorem holds exactly, so their difference is the quadrature error
!> of the centroid and midpoint rules, O(h^2): a factor 4 per refinement. This also fixes the sign of the waterline term.
subroutine test_identity_smooth(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8) :: e(3)
  integer :: m
  character(len=200) :: line

  do m = 1, 3
    call make_cylinder(10d0, 10d0, cyl_ntheta(m), cyl_nz(m), cyl_nr(m), mesh, .true.)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    e(m) = identity_difference(mesh, geom, .false.)
  end do
  write(line, '(a,3es11.3)') "relative differences at 448, 1792, 7168 panels: ", e
  print '(4x,a)', trim(line)
  call check(error, e(1) / e(2) >= 3.5d0 .and. e(2) / e(3) >= 3.5d0, "falls by less than 3.5x per refinement; " // line)
  if (allocated(error)) return
  call check(error, e(3) < 1d-3, "not below 1e-3 at 7168 panels; " // line)
end subroutine test_identity_smooth

!> Closed sphere, no waterline: the surface terms alone carry the identity. With the vertex-average centroids that HAMS
!> uses the difference is larger (see make_sphere), so for them only the O(h^2) order is checked.
subroutine test_identity_submerged(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8) :: e(3), ev(3)
  integer :: m
  character(len=200) :: line

  do m = 1, 3
    call make_sphere(5d0, [0d0, 0d0, -20d0], sph_nlat(m), mesh)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    call check(error, geom%nwl == 0, "the submerged sphere has waterline segments")
    if (allocated(error)) return
    e(m) = identity_difference(mesh, geom, .false.)
    call make_sphere(5d0, [0d0, 0d0, -20d0], sph_nlat(m), mesh, vertex_average=.true.)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    ev(m) = identity_difference(mesh, geom, .false.)
  end do
  write(line, '(a,3es11.3)') "relative differences at 128, 512, 2048 panels: ", e
  print '(4x,a)', trim(line)
  call check(error, e(2) < e(1) .and. e(3) < e(2), "does not fall with refinement; " // line)
  if (allocated(error)) return
  call check(error, e(3) < 1d-3, "not below 1e-3 at the finest mesh; " // line)
  if (allocated(error)) return
  write(line, '(a,3es11.3)') "vertex-average centroids, relative differences: ", ev
  print '(4x,a)', trim(line)
  call check(error, ev(2) / ev(3) >= 3.5d0, "vertex-average centroids: falls by less than 3.5x from 512 to 2048; " // line)
end subroutine test_identity_submerged

!> As stokes_identity_smooth, with grad psi from velocity_differencing.
subroutine test_identity_differenced(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8) :: e(3)
  integer :: m
  character(len=200) :: line

  do m = 1, 3
    call make_cylinder(10d0, 10d0, cyl_ntheta(m), cyl_nz(m), cyl_nr(m), mesh, .true.)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    e(m) = identity_difference(mesh, geom, .true.)
  end do
  write(line, '(a,3es11.3)') "relative differences at 448, 1792, 7168 panels: ", e
  print '(4x,a)', trim(line)
  call check(error, e(2) < e(1) .and. e(3) < e(2), "does not fall with refinement; " // line)
  if (allocated(error)) return
  call check(error, e(3) < 1d-2, "not below 1e-2 at 7168 panels; " // line)
end subroutine test_identity_differenced

!> The production path: grad psi and u both from velocity_differencing, against the exact direct form. This converges
!> at O(h) only: measured 5.3e-2, 3.0e-2, 1.6e-2 (ratios 1.75, 1.87). The cause is the waterline row. There the
!> neighbours of a panel all lie below or beside it, so the fit of velocity_differencing is one-sided and its vertical
!> derivative has an O(h) error. That error enters the line term directly and again through the extrapolation of u to
!> the waterline. With the exact u on the waterline panels alone (everything else still differenced) the order is
!> O(h^2): 4.1e-3, 1.1e-3, 2.9e-4. The surface terms and the one-sided bottom-edge panels do not limit the order. Here
!> the whole incident wave is differenced, the worst case. In production the incident part can be evaluated exactly at
!> the waterline, leaving only the scattered part one-sided. This test keeps the extrapolation inside
!> hessian_term_converted (no u1_wl, u2_wl), so the assertion is O(h) convergence, with the limit set from the measured
!> values. potential_pair passes u at the waterline from waterline_velocity instead, which is O(h^2):
!> stokes_identity_waterline_velocity.
subroutine test_identity_differenced_velocity(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8) :: e(3)
  integer :: m
  character(len=200) :: line

  do m = 1, 3
    call make_cylinder(10d0, 10d0, cyl_ntheta(m), cyl_nz(m), cyl_nr(m), mesh, .true.)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    e(m) = identity_difference(mesh, geom, .true., differenced_u=.true.)
  end do
  write(line, '(a,3es11.3)') "relative differences at 448, 1792, 7168 panels: ", e
  print '(4x,a)', trim(line)
  call check(error, e(1) / e(2) >= 1.5d0 .and. e(2) / e(3) >= 1.5d0, "falls by less than 1.5x per refinement; " // line)
  if (allocated(error)) return
  call check(error, e(3) < 2.5d-2, "not below 2.5e-2 at 7168 panels; " // line)
end subroutine test_identity_differenced_velocity

!> The pair (1, 2) is the kernel applied to single-frequency evaluations L(a, b), and swapping 1 and 2 conjugates it.
!> L(a, b) = int psi n.[(a.grad) b] dS is bilinear without conjugation. KERNEL_SUM with a zero partner gives
!> pair_plus(a, 0, 0, b) = a b / 4, so L(a, b) = 4 T_SUM(a, 0; 0, b). Then T_DIFF(1, 2) = [L(X1, u2*) + L(X2*, u1)] / 4.
!> The swap gives the conjugate for real psi, as here: T_DIFF(2, 1) = sum psi P(2, 1) = sum psi conj(P(1, 2)).
subroutine test_kernel(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  complex(8), allocatable :: psi(:,:), grad_psi(:,:,:), psi_wl(:,:), u1(:,:), u2(:,:), zu(:,:)
  complex(8) :: x1(6, 1), x2(6, 1), zx(6, 1), t12(nmode), t21(nmode), la(nmode), lb(nmode), ref(nmode)
  real(8) :: scale

  call make_cylinder(10d0, 10d0, cyl_ntheta(1), cyl_nz(1), cyl_nr(1), mesh, .true.)
  call geometry(mesh, geom, error)
  if (allocated(error)) return
  call psi_panels(mesh, psi, grad_psi)
  allocate(psi_wl(geom%nwl, nmode))
  call psi_waterline(mesh, geom, psi, grad_psi, psi_wl)
  call motions(x1, x2)
  call velocity(mesh, w1, b1, u1)
  call velocity(mesh, w2, b2, u2)
  allocate(zu(3, mesh%n))
  zu = (0d0, 0d0)
  zx = (0d0, 0d0)

  call hessian_term_converted(mesh, geom, psi, grad_psi, psi_wl, x1, u1, x2, u2, KERNEL_DIFF, t12)
  call hessian_term_converted(mesh, geom, psi, grad_psi, psi_wl, x2, u2, x1, u1, KERNEL_DIFF, t21)
  call hessian_term_converted(mesh, geom, psi, grad_psi, psi_wl, x1, zu, zx, conjg(u2), KERNEL_SUM, la)
  call hessian_term_converted(mesh, geom, psi, grad_psi, psi_wl, conjg(x2), zu, zx, u1, KERNEL_SUM, lb)
  ref = (4d0 * la + 4d0 * lb) / 4d0
  scale = maxval(abs(t12))

  call check(error, scale > 0d0, "the pair term is zero")
  if (allocated(error)) return
  call check(error, maxval(abs(t12 - ref)) <= 1d-12 * scale, "pair (1, 2) is not the kernel of single-frequency terms")
  if (allocated(error)) return
  call check(error, maxval(abs(t21 - conjg(t12))) <= 1d-12 * scale, "swapping 1 and 2 does not conjugate the pair term")
end subroutine test_kernel

!> The production path with u at the waterline built by waterline_velocity (ruling R7): grad psi and the panel u from
!> velocity_differencing, u at the waterline from the free-surface condition (u_z = (omega^2/g) phi, phi from
!> waterline_extrapolate), the body condition (normal part) and differencing phi along the waterline (tangential part).
!> This removes the one-sided vertical derivative of the waterline row, so the order returns to O(h^2).
subroutine test_identity_waterline_velocity(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8) :: e(3)
  integer :: m, iw
  character(len=200) :: line

  do m = 1, 3
    call make_cylinder(10d0, 10d0, cyl_ntheta(m), cyl_nz(m), cyl_nr(m), mesh, .true.)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    call check(error, all(geom%wl_next > 0) .and. all(geom%wl_prev > 0), "a waterline segment has no loop neighbour")
    if (allocated(error)) return
    call check(error, all(geom%wl_prev(geom%wl_next) == [(iw, iw = 1, geom%nwl)]), "wl_prev is not the inverse of wl_next")
    if (allocated(error)) return
    e(m) = identity_difference(mesh, geom, .true., differenced_u=.true., waterline_u=.true.)
  end do
  write(line, '(a,3es11.3)') "relative differences at 448, 1792, 7168 panels: ", e
  print '(4x,a)', trim(line)
  call check(error, e(1) / e(2) >= 3d0 .and. e(2) / e(3) >= 3d0, "falls by less than 3x per refinement; " // line)
end subroutine test_identity_waterline_velocity

! ---------------------------------------------------------------------------------------------- potential_pair

!> Wave number of omega in depth h, deep water for h < 0.
real(8) function wave_number_h(omega, h)
  real(8), intent(in) :: omega, h
  if (h < 0d0) then
    wave_number_h = omega**2 / g
  else
    wave_number_h = wave_number(omega, h)
  end if
end function wave_number_h

!> First-order state of a unit incident wave (omega, beta) in depth h, with rigid motions xi and no disturbance
!> potential, and its total potential at the waterline.
subroutine make_state(mesh, geom, omega, beta, h, xi, s, phi_wl)
  type(panel_mesh), intent(in) :: mesh
  type(surface_geom), intent(in) :: geom
  real(8), intent(in) :: omega, beta, h
  complex(8), intent(in) :: xi(6, 1)
  type(first_order_state), intent(out) :: s
  complex(8), allocatable, intent(out) :: phi_wl(:)
  integer :: i
  s%omega = omega
  s%k = wave_number_h(omega, h)
  s%nu = omega**2 / g
  s%h = h
  s%amp = 1d0
  s%beta = beta
  s%g = g
  s%rho = rho
  allocate(s%phi_i(mesh%n), s%phi_p(mesh%n), s%dphi_p_dn(mesh%n), s%dphi_dn(mesh%n))
  allocate(s%grad_i(3, mesh%n), s%grad_p(3, mesh%n), s%xi(6, 1))
  s%xi = xi
  s%phi_p = (0d0, 0d0)
  s%grad_p = (0d0, 0d0)
  do i = 1, mesh%n
    s%phi_i(i) = incident_phi(mesh%xc(:, i), omega, s%k, h, 1d0, beta, g)
    call incident_grad(mesh%xc(:, i), omega, s%k, h, 1d0, beta, g, s%grad_i(:, i))
    s%dphi_dn(i) = sum(mesh%nrm(:, i) * s%grad_i(:, i))
  end do
  s%dphi_p_dn = (0d0, 0d0)
  allocate(phi_wl(geom%nwl))
  call waterline_extrapolate(mesh, geom, s, phi_wl)
end subroutine make_state

!> Synthetic psi: the harmonic columns of psi_panels, times e^{0.3 i k} when complex_psi, with the gradient that
!> potential_pair's caller builds (psi_surface_gradient) and the waterline values.
subroutine synthetic_psi(mesh, geom, complex_psi, psi, grad_psi, psi_wl)
  type(panel_mesh), intent(in) :: mesh
  type(surface_geom), intent(in) :: geom
  logical, intent(in) :: complex_psi
  complex(8), allocatable, intent(out) :: psi(:,:), grad_psi(:,:,:), psi_wl(:,:)
  integer :: k
  call psi_panels(mesh, psi, grad_psi)
  if (complex_psi) then
    do k = 1, nmode
      psi(:, k) = psi(:, k) * exp(ci * 0.3d0 * real(k, 8))
    end do
  end if
  call psi_surface_gradient(mesh, geom, psi, grad_psi)
  allocate(psi_wl(geom%nwl, nmode))
  call psi_waterline(mesh, geom, psi, grad_psi, psi_wl)
end subroutine synthetic_psi

subroutine cylinder_448(mesh, geom, error)
  type(panel_mesh), intent(out) :: mesh
  type(surface_geom), intent(out) :: geom
  type(error_type), allocatable, intent(out) :: error
  call make_cylinder(10d0, 10d0, cyl_ntheta(1), cyl_nz(1), cyl_nr(1), mesh, .true.)
  call geometry(mesh, geom, error)
end subroutine cylinder_448

!> With xi = 0 in both states, B1, B2 and B3 are zero exactly: each is a product with a motion.
subroutine test_pot_fixed_body(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s1, s2
  complex(8), allocatable :: psi(:,:), grad_psi(:,:,:), psi_wl(:,:), wl1(:), wl2(:)
  complex(8) :: zx(6, 1), terms(nmode, NPOT_TERMS)

  call cylinder_448(mesh, geom, error)
  if (allocated(error)) return
  zx = (0d0, 0d0)
  call make_state(mesh, geom, w1, b1, depth, zx, s1, wl1)
  call make_state(mesh, geom, w2, b2, depth, zx, s2, wl2)
  call synthetic_psi(mesh, geom, .true., psi, grad_psi, psi_wl)
  call potential_pair(mesh, geom, s1, s2, wl1, wl2, psi, grad_psi, psi_wl, w1 - w2, rho, g, terms)
  call check(error, all(terms(:, TP_B1) == (0d0, 0d0)) .and. all(terms(:, TP_B2) == (0d0, 0d0)) &
             .and. all(terms(:, TP_B3) == (0d0, 0d0)), "B1, B2 or B3 is not zero for a fixed body")
  if (allocated(error)) return
  call check(error, maxval(abs(terms(:, TP_FK))) > 0d0 .and. maxval(abs(terms(:, TP_DIFFRACTION))) > 0d0, &
             "the incident terms are zero")
end subroutine test_pot_fixed_body

!> Swapping the states gives the conjugate when psi at Omega < 0 is taken as conj psi(|Omega|), with complex psi so that a
!> missing conjugation shows. Each term on its own, moving body, finite depth, different headings.
subroutine test_pot_hermitian(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s1, s2
  complex(8), allocatable :: psi(:,:), grad_psi(:,:,:), psi_wl(:,:), wl1(:), wl2(:)
  complex(8) :: x1(6, 1), x2(6, 1), t12(nmode, NPOT_TERMS), t21(nmode, NPOT_TERMS)
  integer :: j
  character(len=80) :: line

  call cylinder_448(mesh, geom, error)
  if (allocated(error)) return
  call motions(x1, x2)
  call make_state(mesh, geom, w1, b1, depth, x1, s1, wl1)
  call make_state(mesh, geom, w2, b2, depth, x2, s2, wl2)
  call synthetic_psi(mesh, geom, .true., psi, grad_psi, psi_wl)
  call potential_pair(mesh, geom, s1, s2, wl1, wl2, psi, grad_psi, psi_wl, 0.1d0, rho, g, t12)
  call potential_pair(mesh, geom, s2, s1, wl2, wl1, psi, grad_psi, psi_wl, 0.1d0, rho, g, t21)
  do j = 1, NPOT_TERMS
    write(line, '(a,i0)') "term ", j
    call check(error, maxval(abs(t12(:, j))) > 0d0, trim(line) // " is zero")
    if (allocated(error)) return
    call check(error, maxval(abs(t21(:, j) - conjg(t12(:, j)))) <= 1d-12 * maxval(abs(t12(:, j))), &
               trim(line) // ": (2, 1) is not conj (1, 2)")
    if (allocated(error)) return
  end do
end subroutine test_pot_hermitian

!> Equal frequencies, different headings (Omega = 0, K /= 0): every term is zero and finite.
subroutine test_pot_equal_frequency(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s1, s2
  complex(8), allocatable :: psi(:,:), grad_psi(:,:,:), psi_wl(:,:), wl1(:), wl2(:)
  complex(8) :: x1(6, 1), x2(6, 1), terms(nmode, NPOT_TERMS)

  call cylinder_448(mesh, geom, error)
  if (allocated(error)) return
  call motions(x1, x2)
  call make_state(mesh, geom, w1, 0d0, depth, x1, s1, wl1)
  call make_state(mesh, geom, w1, 40d0 * pi / 180d0, depth, x2, s2, wl2)
  call synthetic_psi(mesh, geom, .true., psi, grad_psi, psi_wl)
  terms = (1d0, 1d0)
  call potential_pair(mesh, geom, s1, s2, wl1, wl2, psi, grad_psi, psi_wl, 0.1d0, rho, g, terms)
  call check(error, all(terms == (0d0, 0d0)), "the terms are not zero at Omega = 0")
end subroutine test_pot_equal_frequency

!> B3 = int psi_k n.d/dt(H x) dS, H x = (1/2) alpha x (alpha x x), against the Omega component of the time series of
!> (1/2) alpha(t) x (alpha(t) x x), alpha(t) = Re[alpha_1 e^{-i w1 t} + alpha_2 e^{-i w2 t}], over one beat period (64
!> samples). With w1 = 0.8, w2 = 0.7 every frequency of the product is a multiple of Omega below 32 Omega, so the
!> projection is exact. Its Omega part is Re[(T12 + conj T21) e^{-i Omega t}] = Re[2 T12 e^{-i Omega t}], so
!> T12 = (1/N) sum h(t_n) e^{i Omega t_n}; d/dt gives -i Omega.
subroutine test_pot_b3_direct(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s1, s2
  complex(8), allocatable :: psi(:,:), grad_psi(:,:,:), psi_wl(:,:), wl1(:), wl2(:)
  complex(8) :: x1(6, 1), x2(6, 1), terms(nmode, NPOT_TERMS), ref(nmode), acc(3), bi
  integer, parameter :: nt = 64
  real(8) :: om, t, a(3), x(3), hx(3), ax(3)
  integer :: i, n

  call cylinder_448(mesh, geom, error)
  if (allocated(error)) return
  call motions(x1, x2)
  call make_state(mesh, geom, w1, b1, depth, x1, s1, wl1)
  call make_state(mesh, geom, w2, b2, depth, x2, s2, wl2)
  call synthetic_psi(mesh, geom, .true., psi, grad_psi, psi_wl)
  call potential_pair(mesh, geom, s1, s2, wl1, wl2, psi, grad_psi, psi_wl, 0.1d0, rho, g, terms)

  om = w1 - w2
  ref = (0d0, 0d0)
  do i = 1, mesh%n
    x = mesh%xc(:, i) - mesh%xr(:, 1)
    acc = (0d0, 0d0)
    do n = 0, nt - 1
      t = 2d0 * pi * real(n, 8) / (real(nt, 8) * om)
      a = real(x1(4:6, 1) * exp(-ci * w1 * t) + x2(4:6, 1) * exp(-ci * w2 * t))
      ax = cross_r(a, x)
      hx = 0.5d0 * cross_r(a, ax)
      acc = acc + hx * exp(ci * om * t)
    end do
    bi = sum(mesh%nrm(:, i) * (-ci * om * acc / real(nt, 8)))
    ref = ref + psi(i, :) * bi * mesh%area(i)
  end do
  ref = ci * om * rho * ref
  call check(error, maxval(abs(ref)) > 0d0, "the reference B3 is zero")
  if (allocated(error)) return
  call check(error, maxval(abs(terms(:, TP_B3) - ref)) <= 1d-10 * maxval(abs(ref)), "B3 differs from the time series")

contains

  pure function cross_r(p, q) result(c)
    real(8), intent(in) :: p(3), q(3)
    real(8) :: c(3)
    c = [p(2)*q(3) - p(3)*q(2), p(3)*q(1) - p(1)*q(3), p(1)*q(2) - p(2)*q(1)]
  end function cross_r

end subroutine test_pot_b3_direct

!> B1 = int psi_k (Xdot - u).(alpha x n) dS against the Omega component of the time series of the product of the real
!> signals Xdot(t) - u(t) and alpha(t), projected as in pot_b3_direct. Fixes the order of the triple product and the
!> sign of u.
subroutine test_pot_b1_direct(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s1, s2
  complex(8), allocatable :: psi(:,:), grad_psi(:,:,:), psi_wl(:,:), wl1(:), wl2(:)
  complex(8) :: x1(6, 1), x2(6, 1), terms(nmode, NPOT_TERMS), ref(nmode), acc, v1(3), v2(3)
  integer, parameter :: nt = 64
  real(8) :: om, t, a(3), v(3), x(3), an(3)
  integer :: i, n

  call cylinder_448(mesh, geom, error)
  if (allocated(error)) return
  call motions(x1, x2)
  call make_state(mesh, geom, w1, b1, depth, x1, s1, wl1)
  call make_state(mesh, geom, w2, b2, depth, x2, s2, wl2)
  call synthetic_psi(mesh, geom, .true., psi, grad_psi, psi_wl)
  call potential_pair(mesh, geom, s1, s2, wl1, wl2, psi, grad_psi, psi_wl, 0.1d0, rho, g, terms)

  om = w1 - w2
  ref = (0d0, 0d0)
  do i = 1, mesh%n
    x = mesh%xc(:, i) - mesh%xr(:, 1)
    v1 = -ci * w1 * rigid_displacement(x1(:, 1), x) - s1%grad_i(:, i)
    v2 = -ci * w2 * rigid_displacement(x2(:, 1), x) - s2%grad_i(:, i)
    acc = (0d0, 0d0)
    do n = 0, nt - 1
      t = 2d0 * pi * real(n, 8) / (real(nt, 8) * om)
      a = real(x1(4:6, 1) * exp(-ci * w1 * t) + x2(4:6, 1) * exp(-ci * w2 * t))
      v = real(v1 * exp(-ci * w1 * t) + v2 * exp(-ci * w2 * t))
      an = [a(2)*mesh%nrm(3, i) - a(3)*mesh%nrm(2, i), a(3)*mesh%nrm(1, i) - a(1)*mesh%nrm(3, i), &
            a(1)*mesh%nrm(2, i) - a(2)*mesh%nrm(1, i)]
      acc = acc + sum(v * an) * exp(ci * om * t)
    end do
    ref = ref + psi(i, :) * acc / real(nt, 8) * mesh%area(i)
  end do
  ref = ci * om * rho * ref
  call check(error, maxval(abs(ref)) > 0d0, "the reference B1 is zero")
  if (allocated(error)) return
  call check(error, maxval(abs(terms(:, TP_B1) - ref)) <= 1d-10 * maxval(abs(ref)), "B1 differs from the time series")
end subroutine test_pot_b1_direct

!> Deep water, collinear waves, fixed synthetic psi and motions: every bracket of T_pot has a finite limit as Omega -> 0
!> (phi_I2 tends to i omega/2 per unit amplitudes), so T_pot = i Omega rho [...] is proportional to Omega.
subroutine test_pot_deep_scaling(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s1, s2
  complex(8), allocatable :: psi(:,:), grad_psi(:,:,:), psi_wl(:,:), wl1(:), wl2(:)
  complex(8) :: x1(6, 1), x2(6, 1), terms(nmode, NPOT_TERMS)
  real(8), parameter :: dw(3) = [0.02d0, 0.01d0, 0.005d0], wj = 0.8d0
  real(8) :: r(3)
  integer :: m
  character(len=200) :: line

  call cylinder_448(mesh, geom, error)
  if (allocated(error)) return
  call motions(x1, x2)
  call synthetic_psi(mesh, geom, .false., psi, grad_psi, psi_wl)
  call make_state(mesh, geom, wj, 0d0, -1d0, x2, s2, wl2)
  do m = 1, 3
    call make_state(mesh, geom, wj + dw(m), 0d0, -1d0, x1, s1, wl1)
    call potential_pair(mesh, geom, s1, s2, wl1, wl2, psi, grad_psi, psi_wl, dw(m), rho, g, terms)
    r(m) = norm2(abs(sum(terms, dim=2))) / dw(m)
  end do
  write(line, '(a,3es13.5)') "|T_pot| / Omega at Omega = 0.02, 0.01, 0.005: ", r
  print '(4x,a)', trim(line)
  call check(error, abs(r(2) / r(1) - 1d0) <= 0.05d0 .and. abs(r(3) / r(1) - 1d0) <= 0.05d0, &
             "|T_pot| is not proportional to Omega within 5%; " // line)
end subroutine test_pot_deep_scaling

!> phi_I2 replaced by a constant c: the FK force is i Omega rho c int n dS. On a closed sphere that is 0; on the truncated
!> cylinder only the bottom contributes, n_z = +1 (out of the fluid), so F_z = i Omega rho c A_wp. The diffraction term
!> is 0, since dphi_I2/dn = 0.
subroutine test_pot_fk_uniform(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s1, s2
  complex(8), allocatable :: psi(:,:), grad_psi(:,:,:), psi_wl(:,:), wl1(:), wl2(:)
  complex(8) :: zx(6, 1), terms(nmode, NPOT_TERMS)
  complex(8), parameter :: c = (0.7d0, -0.4d0)
  real(8) :: om, ref

  om = w1 - w2
  zx = (0d0, 0d0)
  call make_sphere(5d0, [0d0, 0d0, -20d0], sph_nlat(1), mesh)
  call geometry(mesh, geom, error)
  if (allocated(error)) return
  call make_state(mesh, geom, w1, b1, depth, zx, s1, wl1)
  call make_state(mesh, geom, w2, b2, depth, zx, s2, wl2)
  call synthetic_psi(mesh, geom, .true., psi, grad_psi, psi_wl)
  call potential_pair(mesh, geom, s1, s2, wl1, wl2, psi, grad_psi, psi_wl, 0.1d0, rho, g, terms, phi2_constant=c)
  ref = om * rho * abs(c) * sum(mesh%area)
  call check(error, norm2(abs(terms(1:3, TP_FK))) <= 1d-12 * ref, "sphere: the FK force of a constant is not zero")
  if (allocated(error)) return
  call check(error, all(terms(:, TP_DIFFRACTION) == (0d0, 0d0)), "sphere: the diffraction term of a constant is not zero")
  if (allocated(error)) return

  call cylinder_448(mesh, geom, error)
  if (allocated(error)) return
  call make_state(mesh, geom, w1, b1, depth, zx, s1, wl1)
  call make_state(mesh, geom, w2, b2, depth, zx, s2, wl2)
  call synthetic_psi(mesh, geom, .true., psi, grad_psi, psi_wl)
  call potential_pair(mesh, geom, s1, s2, wl1, wl2, psi, grad_psi, psi_wl, 0.1d0, rho, g, terms, phi2_constant=c)
  ref = om * rho * abs(c) * geom%awp(1)
  call check(error, abs(terms(1, TP_FK)) <= 1d-10 * ref .and. abs(terms(2, TP_FK)) <= 1d-10 * ref, &
             "cylinder: F_x or F_y of a constant is not zero")
  if (allocated(error)) return
  call check(error, abs(terms(3, TP_FK) - ci * om * rho * c * geom%awp(1)) <= 1d-10 * ref, &
             "cylinder: F_z is not i Omega rho c A_wp")
end subroutine test_pot_fk_uniform

!> A psi solved at another Omega, or a psi that was never stored (all zero), is refused: ok = .false., terms zero.
subroutine test_pot_psi_missing(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s1, s2
  complex(8), allocatable :: psi(:,:), grad_psi(:,:,:), psi_wl(:,:), wl1(:), wl2(:)
  complex(8) :: x1(6, 1), x2(6, 1), terms(nmode, NPOT_TERMS)
  logical :: ok

  call cylinder_448(mesh, geom, error)
  if (allocated(error)) return
  call motions(x1, x2)
  call make_state(mesh, geom, w1, b1, depth, x1, s1, wl1)
  call make_state(mesh, geom, w2, b2, depth, x2, s2, wl2)
  call synthetic_psi(mesh, geom, .true., psi, grad_psi, psi_wl)
  call potential_pair(mesh, geom, s1, s2, wl1, wl2, psi, grad_psi, psi_wl, 0.1d0, rho, g, terms, ok=ok)
  call check(error, ok, "a matching psi is refused")
  if (allocated(error)) return
  call potential_pair(mesh, geom, s1, s2, wl1, wl2, psi, grad_psi, psi_wl, 0.2d0, rho, g, terms, ok=ok)
  call check(error, .not. ok .and. all(terms == (0d0, 0d0)), "psi of another Omega is not refused")
  if (allocated(error)) return
  psi = (0d0, 0d0)
  call potential_pair(mesh, geom, s1, s2, wl1, wl2, psi, grad_psi, psi_wl, 0.1d0, rho, g, terms, ok=ok)
  call check(error, .not. ok .and. all(terms == (0d0, 0d0)), "a psi that was never stored is not refused")
end subroutine test_pot_psi_missing

!> Signs of B2 and of the incident-diffraction term, with a disturbance in both states (review I1). phi_p is a harmonic
!> synthetic field with its analytic gradient in grad_p, so u = grad_i + grad_p differs from the incident wave, and the
!> waterline potential comes from waterline_extrapolate of that state. The references are built from the public parts:
!>   B2 = -i Omega rho hessian_term_converted(u = grad_i + grad_p, u_wl from waterline_velocity with v_n = -i w X.n)
!>   diffraction = -i Omega rho sum psi n.grad phi_I2 A
!> A flipped sign of either term, or u without grad_p, fails here.
subroutine test_pot_signs_disturbance(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s1, s2
  type(incident2) :: inc
  complex(8), allocatable :: psi(:,:), grad_psi(:,:,:), psi_wl(:,:), wl1(:), wl2(:), u1(:,:), u2(:,:), uw1(:,:), uw2(:,:)
  complex(8) :: x1(6, 1), x2(6, 1), terms(nmode, NPOT_TERMS), hs(nmode), ref(nmode), gr2(3)
  real(8) :: om
  integer :: i

  call cylinder_448(mesh, geom, error)
  if (allocated(error)) return
  call motions(x1, x2)
  call make_state(mesh, geom, w1, b1, depth, x1, s1, wl1)
  call make_state(mesh, geom, w2, b2, depth, x2, s2, wl2)
  call add_disturbance(s1, 2, (0.6d0, 0.3d0), wl1)
  call add_disturbance(s2, 4, (-0.2d0, 0.5d0), wl2)
  call synthetic_psi(mesh, geom, .true., psi, grad_psi, psi_wl)
  call potential_pair(mesh, geom, s1, s2, wl1, wl2, psi, grad_psi, psi_wl, 0.1d0, rho, g, terms)
  om = w1 - w2

  ! B2
  allocate(u1(3, mesh%n), u2(3, mesh%n), uw1(3, geom%nwl), uw2(3, geom%nwl))
  u1 = s1%grad_i + s1%grad_p
  u2 = s2%grad_i + s2%grad_p
  call waterline_velocity(mesh, geom, s1, wl1, normal_velocity(s1), uw1)
  call waterline_velocity(mesh, geom, s2, wl2, normal_velocity(s2), uw2)
  call hessian_term_converted(mesh, geom, psi, grad_psi, psi_wl, x1, u1, x2, u2, KERNEL_DIFF, hs, uw1, uw2)
  ref = -ci * om * rho * hs
  call check(error, maxval(abs(ref)) > 0d0, "the reference B2 is zero")
  if (allocated(error)) return
  call check(error, maxval(abs(terms(:, TP_B2) - ref)) <= 1d-12 * maxval(abs(ref)), &
             "B2 is not -i Omega rho times the converted Hessian term of u = grad_i + grad_p")
  if (allocated(error)) return

  ! Incident-diffraction
  call incident2_pair(w1, s1%k, b1, w2, s2%k, b2, depth, g, inc)
  ref = (0d0, 0d0)
  do i = 1, mesh%n
    call incident2_grad(inc, mesh%xc(:, i), gr2)
    ref = ref + psi(i, :) * sum(mesh%nrm(:, i) * gr2) * mesh%area(i)
  end do
  ref = -ci * om * rho * ref
  call check(error, maxval(abs(ref)) > 0d0, "the reference diffraction term is zero")
  if (allocated(error)) return
  call check(error, maxval(abs(terms(:, TP_DIFFRACTION) - ref)) <= 1d-12 * maxval(abs(ref)), &
             "the diffraction term is not -i Omega rho int psi dphi_I2/dn dS")

contains

  !> phi_p = c psi_field(k), grad_p its analytic gradient, and the total waterline potential again.
  subroutine add_disturbance(s, k, c, phi_wl)
    type(first_order_state), intent(inout) :: s
    integer, intent(in) :: k
    complex(8), intent(in) :: c
    complex(8), intent(inout) :: phi_wl(:)
    complex(8) :: p, gp(3)
    integer :: ip
    do ip = 1, mesh%n
      call psi_field(k, mesh%xc(:, ip), p, gp)
      s%phi_p(ip) = c * p
      s%grad_p(:, ip) = c * gp
    end do
    call waterline_extrapolate(mesh, geom, s, phi_wl)
  end subroutine add_disturbance

  !> -i omega X.n at the waterline midpoints, with the panel normal
  function normal_velocity(s) result(vn)
    type(first_order_state), intent(in) :: s
    complex(8) :: vn(geom%nwl)
    integer :: iw
    do iw = 1, geom%nwl
      vn(iw) = -ci * s%omega * sum(mesh%nrm(:, geom%wl_panel(iw)) &
                                   * rigid_displacement(s%xi(:, 1), waterline_point(geom, iw) - mesh%xr(:, 1)))
    end do
  end function normal_velocity

end subroutine test_pot_signs_disturbance

end module test_qtf_potential
