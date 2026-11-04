!> Unit tests for the incident wave helpers, symmetry image expansion and first-order state conversion.
module test_first_order_state
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderTypes, only : panel_mesh, first_order_state
  use FirstOrderState, only : depth_factor, incident_phi, incident_grad, expand_images, build_state
  use second_order_fixtures, only : make_cylinder, set_panel
  implicit none
  private

  public :: collect_first_order_state, collect_first_order_state_serial

  real(8), parameter :: g = 9.80665d0
  real(8), parameter :: pi = 3.141592653589793d0
  complex(8), parameter :: ci = (0d0, 1d0)

contains

subroutine collect_first_order_state(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("incident_elevation", test_incident_elevation), &
    new_unittest("incident_grad_fd", test_incident_grad_fd), &
    new_unittest("incident_laplace_and_fs", test_incident_laplace_and_fs), &
    new_unittest("incident_large_kh", test_incident_large_kh), &
    new_unittest("build_state_conjugates", test_build_state_conjugates), &
    new_unittest("build_state_fixed", test_build_state_fixed), &
    new_unittest("build_state_all_bodies", test_build_state_all_bodies), &
    new_unittest("build_state_image_index", test_build_state_image_index), &
    new_unittest("expand_images_y", test_expand_images_y), &
    new_unittest("expand_images_x", test_expand_images_x), &
    new_unittest("expand_images_single", test_expand_images_single), &
    new_unittest("cylinder_fixture", test_cylinder_fixture) &
    ]
end subroutine collect_first_order_state

!> Sets WaveDyn_mod globals through VINP and DINP, so it must not run in parallel with other tests.
subroutine collect_first_order_state_serial(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("matches_hams_vinp_dinp", test_matches_hams_vinp_dinp) &
    ]
end subroutine collect_first_order_state_serial

! ---------------------------------------------------------------------------------------------- helpers

!> Root of omega^2 = g k tanh(k h) by Newton iteration. h < 0 means deep water.
real(8) function wavenumber(omega, h)
  real(8), intent(in) :: omega, h
  integer :: it
  real(8) :: k, f, df, t

  k = omega**2 / g
  if (h < 0d0) then
    wavenumber = k
    return
  end if
  do it = 1, 100
    t = tanh(k*h)
    f = g*k*t - omega**2
    df = g*t + g*k*h*(1d0 - t**2)
    k = k - f/df
    if (abs(f) < 1d-14) exit
  end do
  wavenumber = k
end function wavenumber

logical function close_c(a, b, tol)
  complex(8), intent(in) :: a, b
  real(8), intent(in) :: tol
  close_c = abs(a - b) <= tol * max(abs(b), 1d-300)
end function close_c

!> Incident potential and gradient at x, with the wavenumber from the dispersion relation.
subroutine incident_at(omega, h, amp, beta, x, phi, grad)
  real(8), intent(in) :: omega, h, amp, beta, x(3)
  complex(8), intent(out) :: phi, grad(3)
  real(8) :: k

  k = wavenumber(omega, h)
  phi = incident_phi(x, omega, k, h, amp, beta, g)
  call incident_grad(x, omega, k, h, amp, beta, g, grad)
end subroutine incident_at

!> Two flat panels in one or two bodies, built by hand so the expected values can be checked on paper.
subroutine two_panel_mesh(nbody, mesh)
  integer, intent(in) :: nbody
  type(panel_mesh), intent(out) :: mesh
  integer :: i

  mesh%n = 2; mesh%nbase = 2; mesh%nbody = nbody; mesh%nsys = 1; mesh%nnode = 0
  allocate(mesh%xc(3,2), mesh%nrm(3,2), mesh%area(2), mesh%psize(2), mesh%nj(6,2))
  allocate(mesh%xnode(3,0), mesh%ncn(2), mesh%ncon(4,2), mesh%body(2), mesh%base(2), mesh%image(2), mesh%xr(3,nbody))
  mesh%xr = 0d0
  mesh%ncn = 4; mesh%ncon = 0
  mesh%area = [1d0, 2d0]; mesh%psize = [1d0, 1.5d0]
  mesh%xc(:,1) = [0.3d0, -0.2d0, -0.5d0]
  mesh%xc(:,2) = [-0.4d0, 0.6d0, -1.5d0]
  mesh%nrm(:,1) = [1d0, 0d0, 0d0]
  mesh%nrm(:,2) = [0.6d0, 0d0, -0.8d0]
  mesh%nj = 0d0
  mesh%nj(1:3,:) = mesh%nrm
  mesh%nj(3,1) = 0.5d0     ! a generalized (non-unit) component, so the sum is not just n
  mesh%nj(4,2) = 0.25d0
  mesh%body = [(min(i, nbody), i = 1, 2)]
  mesh%base = [1, 2]
  mesh%image = 1
end subroutine two_panel_mesh

! ---------------------------------------------------------------------------------------------- incident wave

subroutine test_incident_elevation(error)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: hs(2), omega, k, amp
  complex(8) :: eta
  integer :: i

  omega = 1.2d0; amp = 1.7d0
  hs = [-1d0, 15d0]
  do i = 1, 2
    k = wavenumber(omega, hs(i))
    eta = ci * omega / g * incident_phi([0d0, 0d0, 0d0], omega, k, hs(i), amp, 0.7d0, g)
    call check(error, abs(eta - amp) < 1d-12 * amp, "elevation differs from amp at the origin")
    if (allocated(error)) return
  end do
end subroutine test_incident_elevation

subroutine test_incident_grad_fd(error)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: hs(2), omega, k, amp, beta, x(3), e(3), step
  complex(8) :: grad(3), fd
  integer :: i, d

  omega = 1.2d0; amp = 1.7d0; beta = 0.7d0; step = 1d-5
  x = [1.3d0, -0.8d0, -2.2d0]
  hs = [-1d0, 15d0]
  do i = 1, 2
    k = wavenumber(omega, hs(i))
    call incident_grad(x, omega, k, hs(i), amp, beta, g, grad)
    do d = 1, 3
      e = 0d0; e(d) = step
      fd = (incident_phi(x + e, omega, k, hs(i), amp, beta, g) - incident_phi(x - e, omega, k, hs(i), amp, beta, g)) &
           / (2d0 * step)
      call check(error, close_c(grad(d), fd, 1d-6), "incident_grad differs from central differences")
      if (allocated(error)) return
    end do
  end do
end subroutine test_incident_grad_fd

subroutine test_incident_laplace_and_fs(error)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: hs(2), omega, k, amp, beta, x(3), e(3), step, nu
  complex(8) :: lap, p0, grad(3)
  integer :: i, d

  omega = 1.2d0; amp = 1.7d0; beta = 0.7d0; step = 1d-3
  nu = omega**2 / g
  hs = [-1d0, 15d0]
  do i = 1, 2
    k = wavenumber(omega, hs(i))
    x = [1.3d0, -0.8d0, -2.2d0]
    p0 = incident_phi(x, omega, k, hs(i), amp, beta, g)
    lap = -6d0 * p0
    do d = 1, 3
      e = 0d0; e(d) = step
      lap = lap + incident_phi(x + e, omega, k, hs(i), amp, beta, g) + incident_phi(x - e, omega, k, hs(i), amp, beta, g)
    end do
    lap = lap / step**2
    call check(error, abs(lap) < 1d-6 * k**2 * abs(p0), "incident potential is not harmonic")
    if (allocated(error)) return

    ! free surface condition d(phi)/dz = nu phi at z = 0
    x(3) = 0d0
    p0 = incident_phi(x, omega, k, hs(i), amp, beta, g)
    call incident_grad(x, omega, k, hs(i), amp, beta, g, grad)
    call check(error, close_c(grad(3), nu * p0, 1d-6), "free surface condition violated")
    if (allocated(error)) return
  end do
end subroutine test_incident_laplace_and_fs

subroutine test_incident_large_kh(error)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: k, omega, x(3), zf, dzf
  complex(8) :: pf, pd, gf(3), gd(3)

  k = 5d0; omega = sqrt(g * k)
  x = [0.4d0, 1.1d0, -0.6d0]
  pf = incident_phi(x, omega, k, 400d0, 1d0, 0.3d0, g)
  pd = incident_phi(x, omega, k, -1d0, 1d0, 0.3d0, g)
  call incident_grad(x, omega, k, 400d0, 1d0, 0.3d0, g, gf)
  call incident_grad(x, omega, k, -1d0, 1d0, 0.3d0, g, gd)
  call check(error, abs(pf) < huge(1d0) .and. abs(pf - pd) <= 1d-10 * abs(pd), "large-kh potential differs from deep water")
  if (allocated(error)) return
  call check(error, maxval(abs(gf - gd)) <= 1d-10 * maxval(abs(gd)), "large-kh gradient differs from deep water")
  if (allocated(error)) return
  call depth_factor(k, -0.6d0, 400d0, zf, dzf)
  call check(error, abs(zf - exp(-3d0)) <= 1d-10 * exp(-3d0) .and. abs(dzf - k*zf) <= 1d-10 * abs(k*zf), &
             "depth_factor at large kh")
end subroutine test_incident_large_kh

! ---------------------------------------------------------------------------------------------- build_state

subroutine test_build_state_conjugates(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(first_order_state) :: st
  complex(8) :: mxpot(2,7,1), dspl(6), p0, gr(3)
  real(8) :: omega, k

  call two_panel_mesh(1, mesh)
  omega = 2d0; k = wavenumber(omega, -1d0)
  mxpot = (0d0, 0d0)
  mxpot(1,7,1) = (1d0, 0d0)
  mxpot(1,1,1) = (1d0, 1d0)
  mxpot(1,3,1) = (0d0, 1d0)
  dspl = (0d0, 0d0)
  dspl(1) = (1d0, 2d0)      ! e^{+i w t} amplitudes, so xi is the conjugate
  dspl(3) = (0d0, 1d0)
  call build_state(mesh, omega, k, omega**2/g, -1d0, 1d0, 0.4d0, g, 1000d0, mxpot, dspl, .false., st)

  call check(error, close_c(st%xi(1,1), (1d0, -2d0), 1d-14) .and. close_c(st%xi(3,1), (0d0, -1d0), 1d-14), &
             "xi is not conj(dspl)")
  if (allocated(error)) return
  ! hand values: phi_p = i*1 + (-2i)[(1-2i)(1+i) + (-i)(i)] = -2 - 7i ; dphi_dn = (-2i)[(1-2i)*1 + (-i)*0.5] = -5 - 2i
  call check(error, close_c(st%phi_p(1), (-2d0, -7d0), 1d-13), "phi_p differs from hand value")
  if (allocated(error)) return
  call check(error, close_c(st%dphi_dn(1), (-5d0, -2d0), 1d-13), "dphi_dn differs from hand value")
  if (allocated(error)) return

  call incident_at(omega, -1d0, 1d0, 0.4d0, mesh%xc(:,1), p0, gr)
  call check(error, close_c(st%phi_i(1), p0, 1d-14), "phi_i differs from incident_phi")
  if (allocated(error)) return
  call check(error, maxval(abs(st%grad_i(:,1) - gr)) <= 1d-14 * maxval(abs(gr)), "grad_i differs from incident_grad")
  if (allocated(error)) return
  call check(error, close_c(st%dphi_p_dn(1), st%dphi_dn(1) - sum(gr * mesh%nrm(:,1)), 1d-13), "dphi_p_dn")
  if (allocated(error)) return
  call check(error, allocated(st%grad_p), "grad_p not allocated")
  if (allocated(error)) return
  call check(error, all(shape(st%grad_p) == [3, 2]), "grad_p is not (3,n)")
end subroutine test_build_state_conjugates

subroutine test_build_state_fixed(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(first_order_state) :: st
  complex(8) :: mxpot(2,7,1), dspl(6), gr(3), p0
  real(8) :: omega, k

  call two_panel_mesh(1, mesh)
  omega = 2d0; k = wavenumber(omega, -1d0)
  mxpot = (0.3d0, -0.2d0)
  dspl = (1d0, 2d0)
  call build_state(mesh, omega, k, omega**2/g, -1d0, 1d0, 0.4d0, g, 1000d0, mxpot, dspl, .true., st)
  call check(error, all(abs(st%xi) == 0d0), "xi not zero for a fixed body")
  if (allocated(error)) return
  call check(error, all(abs(st%dphi_dn) == 0d0), "dphi_dn not zero for a fixed body")
  if (allocated(error)) return
  ! with no motion the perturbation normal derivative is minus the incident one
  call incident_at(omega, -1d0, 1d0, 0.4d0, mesh%xc(:,2), p0, gr)
  call check(error, close_c(st%dphi_p_dn(2), -sum(gr * mesh%nrm(:,2)), 1d-13), "dphi_p_dn for a fixed body")
  if (allocated(error)) return
  ! phi_p keeps only the scattered part, i * mxpot(:, 7)
  call check(error, close_c(st%phi_p(2), ci * mxpot(2,7,1), 1d-14), "phi_p for a fixed body")
end subroutine test_build_state_fixed

subroutine test_build_state_all_bodies(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(first_order_state) :: st
  complex(8) :: mxpot(2,13,1), dspl(12), expected
  real(8) :: omega, k

  call two_panel_mesh(2, mesh)       ! panel 1 on body 1, panel 2 on body 2
  omega = 1.5d0; k = wavenumber(omega, -1d0)
  mxpot = (0d0, 0d0)
  mxpot(1,1,1) = (1d0, 0d0)          ! body 1 surge, on a body 1 panel
  mxpot(1,7,1) = (0d0, 2d0)          ! body 2 surge, non-zero on a body 1 panel
  mxpot(1,13,1) = (0.5d0, 0d0)
  dspl = (0d0, 0d0)
  dspl(1) = (1d0, 1d0)
  dspl(7) = (2d0, -1d0)
  call build_state(mesh, omega, k, omega**2/g, -1d0, 1d0, 0d0, g, 1000d0, mxpot, dspl, .false., st)

  expected = ci * mxpot(1,13,1) - ci * omega * (conjg(dspl(1)) * mxpot(1,1,1) + conjg(dspl(7)) * mxpot(1,7,1))
  call check(error, close_c(st%phi_p(1), expected, 1d-13), "radiation sum does not run over all bodies")
  if (allocated(error)) return
  ! dphi_dn uses only the motion of the panel's own body. Panel 1 is on body 1, with nj(:,1) = (1, 0, 0.5, 0, 0, 0)
  expected = -ci * omega * conjg(dspl(1)) * 1d0
  call check(error, close_c(st%dphi_dn(1), expected, 1d-13), "dphi_dn picked up motion of another body")
  if (allocated(error)) return
  ! Panel 2 is on body 2 (modes 7..12), with nj(:,2) = (0.6, 0, -0.8, 0.25, 0, 0)
  expected = -ci * omega * conjg(dspl(7)) * 0.6d0
  call check(error, close_c(st%dphi_dn(2), expected, 1d-13), "dphi_dn of the second body")
end subroutine test_build_state_all_bodies

subroutine test_build_state_image_index(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(first_order_state) :: st
  complex(8) :: mxpot(1,7,2), dspl(6)
  real(8) :: xc(1,3), dxyz(1,6), area(1), psize(1), xnode(4,3), xr(3,1)
  integer :: ncn(1), ncon(1,4), body(1)
  real(8) :: omega, k

  xc(1,:) = [0.5d0, 0.7d0, -1d0]
  dxyz(1,:) = [0d0, 1d0, 0d0, 0d0, 0d0, 0.7d0]
  area = 1d0; psize = 1d0; ncn = 4; ncon(1,:) = [1, 2, 3, 4]; body = 1; xr = 0d0
  xnode = 0d0
  call expand_images(xc, dxyz, area, psize, xnode, ncn, ncon, body, xr, 0, 1, 2, mesh)
  omega = 1d0; k = wavenumber(omega, -1d0)
  mxpot = (0d0, 0d0)
  mxpot(1,7,1) = (1d0, 0d0)
  mxpot(1,7,2) = (0d0, 3d0)
  dspl = (0d0, 0d0)
  call build_state(mesh, omega, k, omega**2/g, -1d0, 1d0, 0d0, g, 1000d0, mxpot, dspl, .false., st)
  call check(error, close_c(st%phi_p(1), ci * (1d0, 0d0), 1d-14), "image 1 potential")
  if (allocated(error)) return
  call check(error, close_c(st%phi_p(2), ci * (0d0, 3d0), 1d-14), "image 2 potential")
end subroutine test_build_state_image_index

! ---------------------------------------------------------------------------------------------- expand_images

!> Half of a cube of side 2 centered at the origin (y >= 0), normals into the cube. Twelve panels.
subroutine half_cube(mesh)
  type(panel_mesh), intent(out) :: mesh
  integer :: ip

  mesh%n = 12; mesh%nbase = 12; mesh%nnode = 48; mesh%nbody = 1; mesh%nsys = 1
  allocate(mesh%xc(3,12), mesh%nrm(3,12), mesh%area(12), mesh%psize(12), mesh%nj(6,12))
  allocate(mesh%xnode(3,48), mesh%ncn(12), mesh%ncon(4,12), mesh%body(12), mesh%base(12), mesh%image(12), mesh%xr(3,1))
  mesh%xr = 0d0; mesh%body = 1; mesh%image = 1; mesh%ncn = 4
  ip = 0
  ! x = +-1 faces: y in [0,1], z split in two
  call rect([1d0, 0d0, -1d0], [0d0, 1d0, 0d0], [0d0, 0d0, 1d0], [-1d0, 0d0, 0d0])
  call rect([1d0, 0d0, 0d0], [0d0, 1d0, 0d0], [0d0, 0d0, 1d0], [-1d0, 0d0, 0d0])
  call rect([-1d0, 0d0, -1d0], [0d0, 1d0, 0d0], [0d0, 0d0, 1d0], [1d0, 0d0, 0d0])
  call rect([-1d0, 0d0, 0d0], [0d0, 1d0, 0d0], [0d0, 0d0, 1d0], [1d0, 0d0, 0d0])
  ! z = +-1 faces: y in [0,1], x split in two
  call rect([-1d0, 0d0, 1d0], [1d0, 0d0, 0d0], [0d0, 1d0, 0d0], [0d0, 0d0, -1d0])
  call rect([0d0, 0d0, 1d0], [1d0, 0d0, 0d0], [0d0, 1d0, 0d0], [0d0, 0d0, -1d0])
  call rect([-1d0, 0d0, -1d0], [1d0, 0d0, 0d0], [0d0, 1d0, 0d0], [0d0, 0d0, 1d0])
  call rect([0d0, 0d0, -1d0], [1d0, 0d0, 0d0], [0d0, 1d0, 0d0], [0d0, 0d0, 1d0])
  ! y = 1 face: 2 x 2 panels
  call rect([-1d0, 1d0, -1d0], [1d0, 0d0, 0d0], [0d0, 0d0, 1d0], [0d0, -1d0, 0d0])
  call rect([0d0, 1d0, -1d0], [1d0, 0d0, 0d0], [0d0, 0d0, 1d0], [0d0, -1d0, 0d0])
  call rect([-1d0, 1d0, 0d0], [1d0, 0d0, 0d0], [0d0, 0d0, 1d0], [0d0, -1d0, 0d0])
  call rect([0d0, 1d0, 0d0], [1d0, 0d0, 0d0], [0d0, 0d0, 1d0], [0d0, -1d0, 0d0])
  mesh%base = [(ip, ip = 1, 12)]

contains

  subroutine rect(p0, e1, e2, nv)
    real(8), intent(in) :: p0(3), e1(3), e2(3), nv(3)
    ip = ip + 1
    mesh%xnode(:, 4*ip-3) = p0
    mesh%xnode(:, 4*ip-2) = p0 + e1
    mesh%xnode(:, 4*ip-1) = p0 + e1 + e2
    mesh%xnode(:, 4*ip)   = p0 + e2
    mesh%ncon(:, ip) = [4*ip-3, 4*ip-2, 4*ip-1, 4*ip]
    call set_panel(mesh, ip, nv)
  end subroutine rect

end subroutine half_cube

!> Hands a panel_mesh to expand_images in the layout of the HAMS global arrays.
subroutine expand_from(base, isx, isy, nsys, mesh)
  type(panel_mesh), intent(in) :: base
  integer, intent(in) :: isx, isy, nsys
  type(panel_mesh), intent(out) :: mesh

  call expand_images(transpose(base%xc), transpose(base%nj), base%area, base%psize, transpose(base%xnode), &
                     base%ncn, transpose(base%ncon), base%body, base%xr, isx, isy, nsys, mesh)
end subroutine expand_from

subroutine check_mirror(error, base, mesh, fx, fy)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh), intent(in) :: base, mesh
  real(8), intent(in) :: fx, fy
  real(8) :: f(3), fj(6)
  integer :: nb, i, j

  nb = base%n
  f = [fx, fy, 1d0]
  fj = [fx, fy, 1d0, fy, fx, fx*fy]
  call check(error, mesh%n == 2*nb .and. mesh%nbase == nb .and. mesh%nsys == 2 .and. mesh%nnode == 2*base%nnode, "sizes")
  if (allocated(error)) return
  do i = 1, nb
    call check(error, all(mesh%xc(:,i) == base%xc(:,i)) .and. all(mesh%nj(:,i) == base%nj(:,i)), "image 1 is not a copy")
    if (allocated(error)) return
    call check(error, mesh%base(i) == i .and. mesh%image(i) == 1 .and. mesh%base(nb+i) == i .and. mesh%image(nb+i) == 2, &
               "base and image indices")
    if (allocated(error)) return
    call check(error, all(mesh%xc(:,nb+i) == f * base%xc(:,i)), "image centroid is not mirrored")
    if (allocated(error)) return
    call check(error, all(mesh%nrm(:,nb+i) == f * base%nrm(:,i)), "image normal is not mirrored")
    if (allocated(error)) return
    do j = 1, 6
      call check(error, mesh%nj(j,nb+i) == fj(j) * base%nj(j,i), "image nj factor")
      if (allocated(error)) return
    end do
    call check(error, mesh%area(nb+i) == base%area(i) .and. mesh%psize(nb+i) == base%psize(i) &
               .and. mesh%body(nb+i) == base%body(i), "image area, size or body")
    if (allocated(error)) return
    call check(error, all(mesh%ncon(:,nb+i) == base%ncon(:,i) + base%nnode) .and. mesh%ncn(nb+i) == base%ncn(i), &
               "image connectivity is not offset")
    if (allocated(error)) return
    do j = 1, mesh%ncn(nb+i)
      call check(error, all(mesh%xnode(:, mesh%ncon(j,nb+i)) == f * base%xnode(:, base%ncon(j,i))), &
                 "image node coordinates")
      if (allocated(error)) return
    end do
  end do
end subroutine check_mirror

subroutine test_expand_images_y(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: base, mesh
  real(8) :: closure(3)
  integer :: i

  call half_cube(base)
  call expand_from(base, 0, 1, 2, mesh)
  call check_mirror(error, base, mesh, 1d0, -1d0)
  if (allocated(error)) return
  do i = 1, 3
    closure(i) = sum(mesh%area * mesh%nrm(i,:))
  end do
  call check(error, maxval(abs(closure)) < 1d-14, "sum of area*normal over the closed cube is not zero")
  if (allocated(error)) return
  call check(error, abs(sum(mesh%area) - 24d0) < 1d-13, "cube surface area")
end subroutine test_expand_images_y

subroutine test_expand_images_x(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: base, mesh

  call half_cube(base)
  call expand_from(base, 1, 0, 2, mesh)      ! ISX = 1, ISY = 0: x is mirrored
  call check_mirror(error, base, mesh, -1d0, 1d0)
end subroutine test_expand_images_x

subroutine test_expand_images_single(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: base, mesh

  call half_cube(base)
  call expand_from(base, 0, 0, 1, mesh)
  call check(error, mesh%n == 12 .and. mesh%nsys == 1 .and. mesh%nnode == 48, "single system sizes")
  if (allocated(error)) return
  call check(error, all(mesh%image == 1) .and. all(mesh%base == base%base), "single system indices")
  if (allocated(error)) return
  call check(error, all(mesh%xc == base%xc) .and. all(mesh%nj == base%nj) .and. all(mesh%ncon == base%ncon) &
             .and. all(mesh%xnode == base%xnode), "single system is not a copy")
end subroutine test_expand_images_single

subroutine test_cylinder_fixture(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  real(8) :: r, d, side, bottom
  integer :: i

  r = 2d0; d = 1.5d0
  call make_cylinder(r, d, 48, 6, 4, mesh, .true.)
  call check(error, mesh%n == 48*6 + 48*4, "cylinder panel count")
  if (allocated(error)) return
  ! normals point into the body: radially inward on the wall, upward on the bottom
  do i = 1, mesh%n
    if (abs(mesh%xc(3,i) + d) > 1d-12) then
      call check(error, sum(mesh%nrm(1:2,i) * mesh%xc(1:2,i)) < 0d0, "wall normal does not point inward")
    else
      call check(error, mesh%nrm(3,i) == 1d0, "bottom normal does not point up")
    end if
    if (allocated(error)) return
  end do
  bottom = sum(mesh%area, mask=abs(mesh%xc(3,:) + d) < 1d-12)
  side = sum(mesh%area) - bottom
  call check(error, abs(bottom - pi*r**2) < 5d-3 * pi*r**2 .and. abs(side - 2d0*pi*r*d) < 5d-3 * 2d0*pi*r*d, &
             "cylinder areas")
  if (allocated(error)) return
  ! with the open waterplane, the vertical component of sum(area*n) is the waterplane area
  call check(error, abs(sum(mesh%area * mesh%nrm(3,:)) - pi*r**2) < 5d-3 * pi*r**2, "vertical closure")
end subroutine test_cylinder_fixture

! ---------------------------------------------------------------------------------------------- HAMS cross-check

subroutine test_matches_hams_vinp_dinp(error)
  use WaveDyn_mod, only : H, AMP, V, W1, WK
  use PatcVelct, only : VINP, DINP
  type(error_type), allocatable, intent(out) :: error
  real(8) :: hs(2), omega, k, beta, x(3)
  complex(8) :: dpx, dpy, dpz, gr(3), v0, p
  integer :: i

  omega = 1.2d0; beta = 0.6d0
  x = [2d0, -1.5d0, -3d0]
  hs = [70d0, -1d0]
  do i = 1, 2
    k = wavenumber(omega, hs(i))
    H = hs(i); AMP = 1.3d0; V = omega**2 / g; W1 = omega; WK = k
    v0 = VINP(x(1), x(2), x(3), 0d0, 0d0, beta)
    call DINP(x(1), x(2), x(3), 0d0, 0d0, beta, dpx, dpy, dpz)
    p = incident_phi(x, omega, k, hs(i), AMP, beta, g)
    call check(error, close_c(ci * v0, p, 1d-11), "incident_phi differs from i*VINP")
    if (allocated(error)) return
    call incident_grad(x, omega, k, hs(i), AMP, beta, g, gr)
    ! DINP is the gradient of VINP, so i*DINP is the gradient of the physical potential
    call check(error, close_c(ci * dpx, gr(1), 1d-11) .and. close_c(ci * dpy, gr(2), 1d-11) &
               .and. close_c(ci * dpz, gr(3), 1d-11), "incident_grad differs from i*DINP")
    if (allocated(error)) return
  end do
end subroutine test_matches_hams_vinp_dinp

end module test_first_order_state
