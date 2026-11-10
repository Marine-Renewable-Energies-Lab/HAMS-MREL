!> Unit tests of SharpEdge: the sharp edges of a hull and the panel weights for the velocity singularity there.
module test_sharp_edge
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderTypes, only : panel_mesh
  use FirstOrderState, only : expand_images
  use SharpEdge, only : sharp_edge_weights, edge_exponent, edge_strip_factor, DEFAULT_EDGE_DEG
  use second_order_fixtures, only : make_cylinder, set_panel
  implicit none
  private

  public :: collect_sharp_edge

  real(8), parameter :: pi = 3.14159265358979323846d0
  ! f of the first row (rho_c = h/2) and of the second row (rho_c = 3h/2) for p = -1/3, from the closed forms
  ! 2^(2/3) / (4/3) and (2^(2/3) - 1) (3/2) (3/2)^(1/3)
  real(8), parameter :: f_first = 2d0**(2d0/3d0) * 0.75d0
  real(8), parameter :: f_second = (2d0**(2d0/3d0) - 1d0) * 1.5d0 * 1.5d0**(1d0/3d0)

contains

subroutine collect_sharp_edge(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("edge_strip_factor", test_strip_factor), &
    new_unittest("edge_exponent", test_exponent), &
    new_unittest("edge_box", test_box), &
    new_unittest("edge_cylinder", test_cylinder), &
    new_unittest("edge_cylinder_half", test_cylinder_half) &
  ]
end subroutine collect_sharp_edge

!> f of the first and second rows of equal width against the closed forms, 1.1906 and 1.0086 for p = -1/3; the width does
!> not matter, and p = 0 gives 1
subroutine test_strip_factor(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: p = -1d0/3d0
  real(8) :: h
  integer :: m

  do m = 1, 3
    h = 0.1d0 * 7d0**(m - 1)
    call check(error, abs(edge_strip_factor(0.5d0*h, h, p) - f_first) <= 1d-14, "first row against the closed form")
    if (allocated(error)) return
    call check(error, abs(edge_strip_factor(1.5d0*h, h, p) - f_second) <= 1d-14, "second row against the closed form")
    if (allocated(error)) return
  end do
  call check(error, abs(f_first - 1.1906d0) < 5d-5 .and. abs(f_second - 1.0086d0) < 5d-5, "1.1906 and 1.0086")
  if (allocated(error)) return
  call check(error, abs(edge_strip_factor(0.5d0, 1d0, 0d0) - 1d0) <= 1d-15, "p = 0 gives 1")
end subroutine test_strip_factor

!> p = pi/alpha_f - 1 with alpha_f = 2 pi minus the interior angle: -1/3 at 90 deg, -1/5 at 135 deg, 0 at 180 deg
subroutine test_exponent(error)
  type(error_type), allocatable, intent(out) :: error

  call check(error, abs(edge_exponent(0.5d0*pi) + 1d0/3d0) <= 1d-15, "90 deg edge: p = -1/3")
  if (allocated(error)) return
  call check(error, abs(edge_exponent(0.75d0*pi) + 0.2d0) <= 1d-15, "135 deg edge: p = -1/5")
  if (allocated(error)) return
  call check(error, abs(edge_exponent(pi)) <= 1d-15, "flat: p = 0")
end subroutine test_exponent

!> A closed box (submerged, 4 x 3 x 2 panels a side) has 12 edges of 4 (nx + ny + nz) panel edges, all 90 deg: a panel on
!> one edge has f = 1.1906, a corner panel on two edges its square, the others 1
subroutine test_box(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: nx = 4, ny = 3, nz = 2
  type(panel_mesh) :: mesh
  real(8), allocatable :: f(:)
  integer :: nedge(1), nseg(1), nweighted(1), i, k, nexp
  real(8) :: expected

  call make_closed_box(8d0, 6d0, 4d0, nx, ny, nz, mesh)
  allocate(f(mesh%n))
  call sharp_edge_weights(mesh, DEFAULT_EDGE_DEG, f, nedge, nseg, nweighted)
  call check(error, nedge(1), 12, "edges of the box")
  if (allocated(error)) return
  call check(error, nseg(1), 4*(nx + ny + nz), "panel edges on them")
  if (allocated(error)) return
  nexp = 0
  do i = 1, mesh%n
    k = edges_touched(mesh, i)
    expected = f_first**k
    if (k > 0) nexp = nexp + 1
    call check(error, abs(f(i) - expected) <= 1d-13, "weight of a box panel")
    if (allocated(error)) return
  end do
  call check(error, nweighted(1), nexp, "panels weighted")
end subroutine test_box

!> The cylinder of the QTF fixtures (32 x 8 wall, 6 bottom rings): one edge, the circular bottom edge of 32 panel edges; the
!> wall panels meet at 11.25 deg and the bottom panels are coplanar, so no seam is an edge. The bottom row of the wall and
!> the outer ring of the bottom are weighted with 2^p / (p + 1), p from the interior angle (a little above 90 deg on a
!> polygon), the others not
subroutine test_cylinder(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  real(8), allocatable :: f(:)
  integer :: nedge(1), nseg(1), nweighted(1)

  call make_cylinder(10d0, 10d0, 32, 8, 6, mesh, .true.)
  allocate(f(mesh%n))
  call sharp_edge_weights(mesh, DEFAULT_EDGE_DEG, f, nedge, nseg, nweighted)
  call check_cylinder(error, mesh, f, nedge(1), nseg(1), nweighted(1))
end subroutine test_cylinder

!> The half of that cylinder (y > 0) expanded by its image: the same edge, panel edges and weights as the whole mesh; the
!> edges in the symmetry plane, between a panel and its image, are not sharp
subroutine test_cylinder_half(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: full, half, mesh
  real(8), allocatable :: f(:), ff(:)
  integer :: nedge(1), nseg(1), nweighted(1), i, j

  call make_cylinder(10d0, 10d0, 32, 8, 6, full, .true.)
  call half_mesh(full, mesh)
  allocate(f(mesh%n), ff(full%n))
  call sharp_edge_weights(mesh, DEFAULT_EDGE_DEG, f, nedge, nseg, nweighted)
  call check_cylinder(error, mesh, f, nedge(1), nseg(1), nweighted(1))
  if (allocated(error)) return
  call sharp_edge_weights(full, DEFAULT_EDGE_DEG, ff, nedge, nseg, nweighted)
  do i = 1, mesh%n
    j = minloc(sum((full%xc - spread(mesh%xc(:, i), 2, full%n))**2, dim=1), dim=1)
    call check(error, abs(f(i) - ff(j)) <= 1d-13, "weight of the half mesh differs from the whole mesh")
    if (allocated(error)) return
  end do
contains
  subroutine half_mesh(src, out)
    type(panel_mesh), intent(in) :: src
    type(panel_mesh), intent(out) :: out
    integer, allocatable :: idx(:)
    integer :: k
    idx = pack([(k, k = 1, src%n)], src%xc(2, :) > 0d0)
    call expand_images(transpose(src%xc(:, idx)), transpose(src%nrm(:, idx)), transpose(src%nj(:, idx)), src%area(idx), &
                       src%psize(idx), transpose(src%xnode), src%ncn(idx), transpose(src%ncon(:, idx)), src%body(idx), &
                       src%xr, 0, 1, 2, out)
  end subroutine half_mesh
end subroutine test_cylinder_half

subroutine check_cylinder(error, mesh, f, nedge, nseg, nweighted)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh), intent(in) :: mesh
  real(8), intent(in) :: f(:)
  integer, intent(in) :: nedge, nseg, nweighted
  real(8) :: p, fe, zbot
  integer :: i
  logical :: first_row

  call check(error, nedge, 1, "edges of the cylinder")
  if (allocated(error)) return
  call check(error, nseg, 32, "panel edges on the bottom edge")
  if (allocated(error)) return
  call check(error, nweighted, 64, "panels weighted: the wall row and the bottom ring at the edge")
  if (allocated(error)) return
  ! interior angle of the polygon edge: 90 deg between wall and bottom, measured between the normals
  p = edge_exponent(0.5d0*pi)
  fe = 2d0**p / (p + 1d0)
  zbot = minval(mesh%xnode(3, :))
  do i = 1, mesh%n
    first_row = any(abs(mesh%xnode(3, mesh%ncon(1:mesh%ncn(i), i)) - zbot) < 1d-9) .and. &
                (abs(mesh%nrm(3, i)) < 0.5d0 .or. maxval(norm2(mesh%xnode(1:2, mesh%ncon(1:mesh%ncn(i), i)), dim=1)) > 9.999d0)
    if (first_row) then
      call check(error, abs(f(i) - fe) <= 1d-12, "weight of a panel at the bottom edge")
    else
      call check(error, f(i) == 1d0, "weight of a panel away from the edge")
    end if
    if (allocated(error)) return
  end do
end subroutine check_cylinder

!> Number of box edges (lines where two faces meet) that panel i has a panel edge on
integer function edges_touched(mesh, i) result(k)
  type(panel_mesh), intent(in) :: mesh
  integer, intent(in) :: i
  real(8) :: lo(3), hi(3), x(3)
  integer :: c, m, ax
  logical :: on

  lo = minval(mesh%xnode, dim=2)
  hi = maxval(mesh%xnode, dim=2)
  k = 0
  do c = 1, 3   ! the two axes along the face, other than its normal axis
    if (abs(mesh%nrm(c, i)) > 0.5d0) cycle
    do m = 1, 2
      on = .true.
      ax = c
      ! the panel has a node on the face boundary lo(ax) or hi(ax): then one of its edges lies on that box edge
      x = merge(lo, hi, m == 1)
      on = any(abs(mesh%xnode(ax, mesh%ncon(1:mesh%ncn(i), i)) - x(ax)) < 1d-9)
      if (on) k = k + 1
    end do
  end do
end function edges_touched

!> Closed box [-a/2, a/2] x [-b/2, b/2] x [-d - 1, -1] with nx x ny x nz panels, normals into the body
subroutine make_closed_box(a, b, d, nx, ny, nz, mesh)
  real(8), intent(in) :: a, b, d
  integer, intent(in) :: nx, ny, nz
  type(panel_mesh), intent(out) :: mesh
  integer :: ix, iy, iz, ip, n, nn, i

  n = 2*(nx*ny + ny*nz + nx*nz)
  nn = (nx + 1)*(ny + 1)*(nz + 1)
  mesh%n = n; mesh%nbase = n; mesh%nnode = nn; mesh%nbody = 1; mesh%nsys = 1
  allocate (mesh%xc(3, n), mesh%nrm(3, n), mesh%area(n), mesh%psize(n), mesh%nj(6, n))
  allocate (mesh%xnode(3, nn), mesh%ncn(n), mesh%ncon(4, n), mesh%body(n), mesh%base(n), mesh%image(n))
  allocate (mesh%xr(3, 1))
  mesh%xr = 0d0; mesh%ncon = 0; mesh%body = 1; mesh%image = 1
  mesh%base = [(i, i = 1, n)]
  do iz = 0, nz
    do iy = 0, ny
      do ix = 0, nx
        mesh%xnode(:, id(ix, iy, iz)) = [-0.5d0*a + a*real(ix, 8)/real(nx, 8), -0.5d0*b + b*real(iy, 8)/real(ny, 8), &
                                         -1d0 - d*real(iz, 8)/real(nz, 8)]
      end do
    end do
  end do
  ip = 0
  do iy = 0, ny - 1
    do ix = 0, nx - 1
      call quad(id(ix, iy, nz), id(ix+1, iy, nz), id(ix+1, iy+1, nz), id(ix, iy+1, nz), [0d0, 0d0, 1d0])
      call quad(id(ix, iy, 0), id(ix+1, iy, 0), id(ix+1, iy+1, 0), id(ix, iy+1, 0), [0d0, 0d0, -1d0])
    end do
  end do
  do iz = 0, nz - 1
    do iy = 0, ny - 1
      call quad(id(0, iy, iz), id(0, iy+1, iz), id(0, iy+1, iz+1), id(0, iy, iz+1), [1d0, 0d0, 0d0])
      call quad(id(nx, iy, iz), id(nx, iy+1, iz), id(nx, iy+1, iz+1), id(nx, iy, iz+1), [-1d0, 0d0, 0d0])
    end do
    do ix = 0, nx - 1
      call quad(id(ix, 0, iz), id(ix+1, 0, iz), id(ix+1, 0, iz+1), id(ix, 0, iz+1), [0d0, 1d0, 0d0])
      call quad(id(ix, ny, iz), id(ix+1, ny, iz), id(ix+1, ny, iz+1), id(ix, ny, iz+1), [0d0, -1d0, 0d0])
    end do
  end do

contains

  integer function id(jx, jy, jz)
    integer, intent(in) :: jx, jy, jz
    id = 1 + jx + (nx + 1)*(jy + (ny + 1)*jz)
  end function id

  subroutine quad(n1, n2, n3, n4, nrm)
    integer, intent(in) :: n1, n2, n3, n4
    real(8), intent(in) :: nrm(3)
    ip = ip + 1
    mesh%ncn(ip) = 4
    mesh%ncon(:, ip) = [n1, n2, n3, n4]
    call set_panel(mesh, ip, nrm)
  end subroutine quad

end subroutine make_closed_box

end module test_sharp_edge
