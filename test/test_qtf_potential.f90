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
  use NearFieldDrift, only : pair_product, rigid_displacement, KERNEL_DIFF, KERNEL_SUM
  use QTFPotential, only : hessian_term_converted, psi_waterline
  use second_order_fixtures, only : make_cylinder, set_panel
  implicit none
  private

  public :: collect_qtf_potential

  real(8), parameter :: g = 9.81d0, depth = 100d0, pi = 3.14159265358979323846d0
  integer, parameter :: nmode = 6
  real(8), parameter :: a_psi = 0.05d0
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
    new_unittest("stokes_kernel", test_kernel) &
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
subroutine make_sphere(r, c, nlat, mesh)
  real(8), intent(in) :: r, c(3)
  integer, intent(in) :: nlat
  type(panel_mesh), intent(out) :: mesh
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
      mesh%xc(:, ip) = area_centroid(v, mesh%ncn(ip))
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
real(8) function identity_difference(mesh, geom, differenced) result(rel)
  type(panel_mesh), intent(in) :: mesh
  type(surface_geom), intent(in) :: geom
  logical, intent(in) :: differenced
  complex(8), allocatable :: psi(:,:), grad_psi(:,:,:), psi_wl(:,:), u1(:,:), u2(:,:)
  complex(8) :: x1(6, 1), x2(6, 1), conv(nmode), dir(nmode)

  call psi_panels(mesh, psi, grad_psi)
  if (differenced) call psi_differenced(mesh, geom, psi, grad_psi)
  allocate(psi_wl(geom%nwl, nmode))
  call psi_waterline(mesh, geom, psi, grad_psi, psi_wl)
  call motions(x1, x2)
  call velocity(mesh, w1, b1, u1)
  call velocity(mesh, w2, b2, u2)

  call hessian_term_converted(mesh, geom, psi, grad_psi, psi_wl, x1, u1, x2, u2, KERNEL_DIFF, conv)
  call direct_form(mesh, psi, x1, x2, KERNEL_DIFF, dir)
  rel = norm2(abs(conv - dir)) / norm2(abs(dir))

  ! Pair (1, 1): the brief's case on its own, with u2 = u1 at w1
  call velocity(mesh, w1, b1, u2)
  call hessian_term_converted(mesh, geom, psi, grad_psi, psi_wl, x1, u1, x1, u2, KERNEL_DIFF, conv)
  call direct_form_same(dir)
  rel = max(rel, norm2(abs(conv - dir)) / norm2(abs(dir)))

contains

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

!> Closed sphere, no waterline: the surface terms alone carry the identity.
subroutine test_identity_submerged(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8) :: e(3)
  integer :: m
  character(len=200) :: line

  do m = 1, 3
    call make_sphere(5d0, [0d0, 0d0, -20d0], sph_nlat(m), mesh)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    call check(error, geom%nwl == 0, "the submerged sphere has waterline segments")
    if (allocated(error)) return
    e(m) = identity_difference(mesh, geom, .false.)
  end do
  write(line, '(a,3es11.3)') "relative differences at 128, 512, 2048 panels: ", e
  print '(4x,a)', trim(line)
  call check(error, e(2) < e(1) .and. e(3) < e(2), "does not fall with refinement; " // line)
  if (allocated(error)) return
  call check(error, e(3) < 1d-3, "not below 1e-3 at the finest mesh; " // line)
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

end module test_qtf_potential
