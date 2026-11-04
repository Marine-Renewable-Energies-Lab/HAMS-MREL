!> Unit tests for the body velocity by tangential differencing.
module test_surface_velocity
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderTypes, only : panel_mesh, surface_geom, first_order_state
  use FirstOrderState, only : incident_phi, incident_grad
  use SurfaceGeometry, only : build_surface_geom, DEFAULT_SHARP_DEG
  use SurfaceVelocity, only : velocity_differencing
  use second_order_fixtures, only : make_cylinder, set_panel
  implicit none
  private

  public :: collect_surface_velocity

  real(8), parameter :: g = 9.81d0, omega = 1.2d0, amp = 1d0, beta = 0.4d0
  real(8), parameter :: radius = 1d0, draft = 1d0

contains

subroutine collect_surface_velocity(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("linear_field_exact", test_linear_field_exact), &
    new_unittest("known_normal_preserved", test_known_normal_preserved), &
    new_unittest("incident_field_converges", test_incident_field_converges), &
    new_unittest("waterline_row_one_sided", test_waterline_row_one_sided), &
    new_unittest("fallback_counted", test_fallback_counted) &
    ]
end subroutine collect_surface_velocity

! ---------------------------------------------------------------------------------------------- helpers

!> Horizontal plate z = -1 of nx*ny unit panels with normal +z.
subroutine make_plate(nx, ny, mesh)
  integer, intent(in) :: nx, ny
  type(panel_mesh), intent(out) :: mesh
  integer :: n, ix, iy, ip

  n = nx * ny
  mesh%n = n; mesh%nbase = n; mesh%nnode = (nx + 1) * (ny + 1); mesh%nbody = 1; mesh%nsys = 1
  allocate(mesh%xc(3,n), mesh%nrm(3,n), mesh%area(n), mesh%psize(n), mesh%nj(6,n))
  allocate(mesh%xnode(3,mesh%nnode), mesh%ncn(n), mesh%ncon(4,n), mesh%body(n), mesh%base(n), mesh%image(n))
  allocate(mesh%xr(3,1))
  mesh%xr = 0d0
  mesh%body = 1
  mesh%image = 1
  mesh%base = [(ip, ip = 1, n)]
  mesh%ncon = 0
  do iy = 0, ny
    do ix = 0, nx
      mesh%xnode(:, 1 + ix + (nx + 1) * iy) = [real(ix, 8), real(iy, 8), -1d0]
    end do
  end do
  ip = 0
  do iy = 0, ny - 1
    do ix = 0, nx - 1
      ip = ip + 1
      mesh%ncn(ip) = 4
      mesh%ncon(1:4, ip) = [1 + ix + (nx + 1) * iy, 2 + ix + (nx + 1) * iy, 2 + ix + (nx + 1) * (iy + 1), &
                            1 + ix + (nx + 1) * (iy + 1)]
      call set_panel(mesh, ip, [0d0, 0d0, 1d0])
    end do
  end do
end subroutine make_plate

subroutine alloc_state(n, s)
  integer, intent(in) :: n
  type(first_order_state), intent(out) :: s
  allocate(s%phi_p(n), s%dphi_p_dn(n), s%grad_p(3, n))
  s%grad_p = (0d0, 0d0)
end subroutine alloc_state

subroutine geometry(mesh, geom, error)
  type(panel_mesh), intent(in) :: mesh
  type(surface_geom), intent(out) :: geom
  type(error_type), allocatable, intent(out) :: error
  logical :: ok
  character(len=:), allocatable :: msg

  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "build_surface_geom stopped: " // msg)
end subroutine geometry

!> Sets phi_p to the incident wave at the centroids and dphi_p_dn to its exact normal derivative.
subroutine incident_as_disturbance(mesh, k, s)
  type(panel_mesh), intent(in) :: mesh
  real(8), intent(in) :: k
  type(first_order_state), intent(inout) :: s
  complex(8) :: grad(3)
  integer :: i

  do i = 1, mesh%n
    s%phi_p(i) = incident_phi(mesh%xc(:, i), omega, k, -1d0, amp, beta, g)
    call incident_grad(mesh%xc(:, i), omega, k, -1d0, amp, beta, g, grad)
    s%dphi_p_dn(i) = sum(grad * mesh%nrm(:, i))
  end do
end subroutine incident_as_disturbance

!> Area-weighted RMS of |grad_p - incident_grad| over the panels listed in sel (all panels if sel is absent).
real(8) function rms_error(mesh, s, k, sel)
  type(panel_mesh), intent(in) :: mesh
  type(first_order_state), intent(in) :: s
  real(8), intent(in) :: k
  integer, intent(in), optional :: sel(:)
  complex(8) :: grad(3)
  real(8) :: num, den
  integer :: i, m, np

  num = 0d0; den = 0d0
  np = mesh%n
  if (present(sel)) np = size(sel)
  do m = 1, np
    i = m
    if (present(sel)) i = sel(m)
    call incident_grad(mesh%xc(:, i), omega, k, -1d0, amp, beta, g, grad)
    num = num + mesh%area(i) * sum(abs(s%grad_p(:, i) - grad)**2)
    den = den + mesh%area(i)
  end do
  rms_error = sqrt(num / den)
end function rms_error

! ---------------------------------------------------------------------------------------------- tests

subroutine test_linear_field_exact(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s
  complex(8), parameter :: c = (0.7d0, -1.3d0), b = (0.4d0, 0.9d0)
  real(8) :: err
  integer :: i, nfb

  call make_plate(5, 4, mesh)
  call geometry(mesh, geom, error)
  if (allocated(error)) return
  call alloc_state(mesh%n, s)
  ! phi = c x + b z has the exact gradient (c, 0, b)
  do i = 1, mesh%n
    s%phi_p(i) = c * mesh%xc(1, i) + b * mesh%xc(3, i)
    s%dphi_p_dn(i) = b
  end do

  call velocity_differencing(mesh, geom, s, nfb)

  err = 0d0
  do i = 1, mesh%n
    err = max(err, abs(s%grad_p(1, i) - c), abs(s%grad_p(2, i)), abs(s%grad_p(3, i) - b))
  end do
  call check(error, err < 1d-10, "gradient of a linear field is not exact on a flat plate")
  if (allocated(error)) return
  call check(error, nfb == 0, "a plate with a full neighbour set should need no fallback")
end subroutine test_linear_field_exact

subroutine test_known_normal_preserved(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s
  real(8) :: err
  integer :: i, nfb

  call make_cylinder(radius, draft, 24, 6, 3, mesh, .true.)
  call geometry(mesh, geom, error)
  if (allocated(error)) return
  call alloc_state(mesh%n, s)
  do i = 1, mesh%n
    s%phi_p(i) = cmplx(sin(0.3d0 * real(i, 8)), cos(0.7d0 * real(i, 8)), kind=8)
    s%dphi_p_dn(i) = cmplx(cos(1.1d0 * real(i, 8)), sin(0.5d0 * real(i, 8)), kind=8)
  end do

  call velocity_differencing(mesh, geom, s, nfb)

  err = 0d0
  do i = 1, mesh%n
    err = max(err, abs(sum(s%grad_p(:, i) * mesh%nrm(:, i)) - s%dphi_p_dn(i)))
  end do
  call check(error, err < 1d-12, "normal part of grad_p differs from dphi_p_dn")
end subroutine test_known_normal_preserved

subroutine test_incident_field_converges(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: ntheta(3) = [16, 32, 64]
  real(8), parameter :: k = omega**2 / g
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s
  real(8) :: e(3)
  integer :: m, nfb

  ! Sides only. Panels around and down double with each refinement.
  do m = 1, 3
    call make_cylinder(radius, draft, ntheta(m), ntheta(m) / 4, 1, mesh, .false.)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    call alloc_state(mesh%n, s)
    call incident_as_disturbance(mesh, k, s)
    call velocity_differencing(mesh, geom, s, nfb)
    e(m) = rms_error(mesh, s, k)
  end do

  call check(error, e(1) / e(2) >= 1.5d0, "error ratio from 16 to 32 panels around is below 1.5")
  if (allocated(error)) return
  call check(error, e(2) / e(3) >= 1.5d0, "error ratio from 32 to 64 panels around is below 1.5")
end subroutine test_incident_field_converges

subroutine test_waterline_row_one_sided(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: ntheta(3) = [16, 32, 64]
  real(8), parameter :: k = omega**2 / g
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s
  integer, allocatable :: top(:)
  complex(8) :: grad(3)
  real(8) :: e(3), ref
  integer :: m, i, nfb

  ! The top row is one-sided, so its error is first order while the interior is second order.
  do m = 1, 3
    call make_cylinder(radius, draft, ntheta(m), ntheta(m) / 4, 1, mesh, .false.)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    call alloc_state(mesh%n, s)
    call incident_as_disturbance(mesh, k, s)
    call velocity_differencing(mesh, geom, s, nfb)
    ! Panels are numbered ring by ring from the waterline down, so the top row is 1..ntheta
    top = [(i, i = 1, ntheta(m))]
    e(m) = rms_error(mesh, s, k, top)
  end do

  ref = 0d0
  do i = 1, size(top)
    call incident_grad(mesh%xc(:, i), omega, k, -1d0, amp, beta, g, grad)
    ref = ref + mesh%area(i) * sum(abs(grad)**2)
  end do
  ref = sqrt(ref / sum(mesh%area(top)))

  call check(error, e(1) / e(2) >= 1.6d0, "top-row error ratio from 16 to 32 panels around is below 1.6")
  if (allocated(error)) return
  call check(error, e(2) / e(3) >= 1.6d0, "top-row error ratio from 32 to 64 panels around is below 1.6")
  if (allocated(error)) return
  call check(error, e(3) < 0.05d0 * ref, "top-row error at 64 panels around exceeds 5% of the incident gradient")
  if (allocated(error)) then
    write(error%message, '(a,es10.3,a,es10.3)') "top-row rms error ", e(3), " vs incident gradient rms ", ref
  end if
end subroutine test_waterline_row_one_sided

subroutine test_fallback_counted(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s
  integer :: nfb, i

  ! One isolated panel: normal part only, counted once
  allocate(mesh%xc(3,1), mesh%nrm(3,1), mesh%psize(1))
  mesh%n = 1
  mesh%xc = 0d0
  mesh%nrm(:,1) = [0d0, 0d0, 1d0]
  mesh%psize = 1d0
  geom%nbr_start = [1, 1]
  allocate(geom%nbr(0))
  call alloc_state(1, s)
  s%phi_p = (1d0, 2d0)
  s%dphi_p_dn = (0.5d0, -0.25d0)
  call velocity_differencing(mesh, geom, s, nfb)
  call check(error, nfb == 1, "isolated panel should give n_fallback = 1")
  if (allocated(error)) return
  call check(error, all(s%grad_p(1:2,1) == (0d0, 0d0)) .and. s%grad_p(3,1) == (0.5d0, -0.25d0), &
             "isolated panel should keep only the normal part")
  if (allocated(error)) return

  ! Three panels in a row: the ends have one neighbour and the middle two collinear ones, so all three fall back
  deallocate(mesh%xc, mesh%nrm, mesh%psize, geom%nbr, geom%nbr_start, s%phi_p, s%dphi_p_dn, s%grad_p)
  allocate(mesh%xc(3,3), mesh%nrm(3,3), mesh%psize(3))
  mesh%n = 3
  do i = 1, 3
    mesh%xc(:,i) = [real(i, 8), 0d0, 0d0]
    mesh%nrm(:,i) = [0d0, 0d0, 1d0]
  end do
  mesh%psize = 1d0
  geom%nbr_start = [1, 2, 4, 5]
  geom%nbr = [2, 1, 3, 2]
  call alloc_state(3, s)
  s%phi_p = [(1d0, 0d0), (2d0, 0d0), (4d0, 0d0)]
  s%dphi_p_dn = (0d0, 0d0)
  call velocity_differencing(mesh, geom, s, nfb)
  call check(error, nfb == 3, "collinear neighbours should all fall back")
  if (allocated(error)) return
  call check(error, all(s%grad_p == (0d0, 0d0)), "fallback gradient should be zero here and free of NaN")
end subroutine test_fallback_counted

end module test_surface_velocity
