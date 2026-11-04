!> Unit tests for panel neighbours, waterline extraction and per-body waterplane and volume properties.
module test_surface_geometry
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderTypes, only : panel_mesh, surface_geom
  use FirstOrderState, only : expand_images
  use SurfaceGeometry, only : build_surface_geom, lid_rim_side_ratio, DEFAULT_SHARP_DEG
  use second_order_fixtures, only : make_cylinder, set_panel
  implicit none
  private

  public :: collect_surface_geometry

  real(8), parameter :: pi = 3.141592653589793d0
  real(8), parameter :: wl_tol = 1d-6
  real(8), parameter :: THIN_ENOUGH = 0.5d0   ! of a uniform whole disc in test_lid_rim_symmetry_planes

contains

subroutine collect_surface_geometry(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("cylinder_waterline_length", test_cylinder_waterline_length), &
    new_unittest("waterline_normal_into_body", test_waterline_normal_into_body), &
    new_unittest("flared_normal_is_unit", test_flared_normal_is_unit), &
    new_unittest("multi_body_loops", test_multi_body_loops), &
    new_unittest("submerged_body_no_waterline", test_submerged_body_no_waterline), &
    new_unittest("symmetry_images_join", test_symmetry_images_join), &
    new_unittest("nonconforming_neighbours", test_nonconforming_neighbours), &
    new_unittest("conforming_no_proximity", test_conforming_no_proximity), &
    new_unittest("lid_rim_ratio", test_lid_rim_ratio), &
    new_unittest("lid_rim_symmetry_planes", test_lid_rim_symmetry_planes), &
    new_unittest("lid_rim_repeated_node_fans", test_lid_rim_repeated_node_fans), &
    new_unittest("lid_rim_only_across_rim", test_lid_rim_only_across_rim), &
    new_unittest("sharp_edge_excluded", test_sharp_edge_excluded), &
    new_unittest("wp_properties_box", test_wp_properties_box), &
    new_unittest("lid_panels_rejected", test_lid_panels_rejected), &
    new_unittest("horizontal_normal_rejected", test_horizontal_normal_rejected), &
    new_unittest("waterline_off_zero", test_waterline_off_zero), &
    new_unittest("waterline_gap_reported", test_waterline_gap_reported), &
    new_unittest("triangles_supported", test_triangles_supported) &
    ]
end subroutine collect_surface_geometry

! ---------------------------------------------------------------------------------------------- helpers

logical function close_r(a, b, tol)
  real(8), intent(in) :: a, b, tol
  close_r = abs(a - b) <= tol * max(abs(b), 1d-300)
end function close_r

!> Perimeter and area of the regular n-gon inscribed in a circle of radius r.
real(8) function polygon_perimeter(n, r)
  integer, intent(in) :: n
  real(8), intent(in) :: r
  polygon_perimeter = 2d0 * real(n, 8) * r * sin(pi / real(n, 8))
end function polygon_perimeter

real(8) function polygon_area(n, r)
  integer, intent(in) :: n
  real(8), intent(in) :: r
  polygon_area = 0.5d0 * real(n, 8) * r**2 * sin(2d0 * pi / real(n, 8))
end function polygon_area

!> Allocates a single-body mesh of n panels and nnode nodes with xr at the origin. Nodes and panels are set by the caller.
subroutine alloc_mesh(n, nnode, mesh)
  integer, intent(in) :: n, nnode
  type(panel_mesh), intent(out) :: mesh
  integer :: i

  mesh%n = n; mesh%nbase = n; mesh%nnode = nnode; mesh%nbody = 1; mesh%nsys = 1
  allocate(mesh%xc(3,n), mesh%nrm(3,n), mesh%area(n), mesh%psize(n), mesh%nj(6,n))
  allocate(mesh%xnode(3,nnode), mesh%ncn(n), mesh%ncon(4,n), mesh%body(n), mesh%base(n), mesh%image(n))
  allocate(mesh%xr(3,1))
  mesh%xr = 0d0
  mesh%xnode = 0d0
  mesh%ncon = 0
  mesh%body = 1
  mesh%image = 1
  mesh%base = [(i, i = 1, n)]
end subroutine alloc_mesh

subroutine set_quad(mesh, ip, n1, n2, n3, n4, nrm)
  type(panel_mesh), intent(inout) :: mesh
  integer, intent(in) :: ip, n1, n2, n3, n4
  real(8), intent(in) :: nrm(3)
  mesh%ncn(ip) = 4
  mesh%ncon(:, ip) = [n1, n2, n3, n4]
  call set_panel(mesh, ip, nrm)
end subroutine set_quad

!> Open-topped box x in [x0, x0+a], y in [y0, y0+b], z in [-d, 0], wetted faces only, with nx*ny*nz cells and the
!> reference point xr. Nodes form a full grid, so the interior and top-face grid nodes are unused except the rim.
subroutine make_box(x0, y0, a, b, d, nx, ny, nz, xr, mesh)
  real(8), intent(in) :: x0, y0, a, b, d, xr(3)
  integer, intent(in) :: nx, ny, nz
  type(panel_mesh), intent(out) :: mesh
  integer :: ix, iy, iz, ip

  call alloc_mesh(nx*ny + 2*ny*nz + 2*nx*nz, (nx+1)*(ny+1)*(nz+1), mesh)
  mesh%xr(:,1) = xr
  do iz = 0, nz
    do iy = 0, ny
      do ix = 0, nx
        mesh%xnode(:, id(ix, iy, iz)) = [x0 + a*real(ix, 8)/real(nx, 8), y0 + b*real(iy, 8)/real(ny, 8), &
                                         -d*real(iz, 8)/real(nz, 8)]
      end do
    end do
  end do
  ip = 0
  do iy = 0, ny - 1
    do ix = 0, nx - 1
      ip = ip + 1
      call set_quad(mesh, ip, id(ix,iy,nz), id(ix+1,iy,nz), id(ix+1,iy+1,nz), id(ix,iy+1,nz), [0d0, 0d0, 1d0])
    end do
  end do
  do iz = 0, nz - 1
    do iy = 0, ny - 1
      ip = ip + 1
      call set_quad(mesh, ip, id(0,iy,iz), id(0,iy+1,iz), id(0,iy+1,iz+1), id(0,iy,iz+1), [1d0, 0d0, 0d0])
      ip = ip + 1
      call set_quad(mesh, ip, id(nx,iy,iz), id(nx,iy+1,iz), id(nx,iy+1,iz+1), id(nx,iy,iz+1), [-1d0, 0d0, 0d0])
    end do
    do ix = 0, nx - 1
      ip = ip + 1
      call set_quad(mesh, ip, id(ix,0,iz), id(ix+1,0,iz), id(ix+1,0,iz+1), id(ix,0,iz+1), [0d0, 1d0, 0d0])
      ip = ip + 1
      call set_quad(mesh, ip, id(ix,ny,iz), id(ix+1,ny,iz), id(ix+1,ny,iz+1), id(ix,ny,iz+1), [0d0, -1d0, 0d0])
    end do
  end do

contains

  integer function id(ix, iy, iz)
    integer, intent(in) :: ix, iy, iz
    id = 1 + ix + (nx + 1) * (iy + (ny + 1) * iz)
  end function id

end subroutine make_box

!> Copy of src with only the panels flagged in keep. Nodes are unchanged.
subroutine subset_mesh(src, keep, mesh)
  type(panel_mesh), intent(in) :: src
  logical, intent(in) :: keep(:)
  type(panel_mesh), intent(out) :: mesh
  integer, allocatable :: idx(:)
  integer :: i

  idx = pack([(i, i = 1, src%n)], keep)
  call alloc_mesh(size(idx), src%nnode, mesh)
  mesh%xnode = src%xnode
  mesh%xr = src%xr
  mesh%xc = src%xc(:, idx)
  mesh%nrm = src%nrm(:, idx)
  mesh%area = src%area(idx)
  mesh%psize = src%psize(idx)
  mesh%nj = src%nj(:, idx)
  mesh%ncn = src%ncn(idx)
  mesh%ncon = src%ncon(:, idx)
end subroutine subset_mesh

!> Joins two single-body meshes into a two-body mesh, the second shifted by dx along x.
subroutine join_bodies(m1, m2, dx, mesh)
  type(panel_mesh), intent(in) :: m1, m2
  real(8), intent(in) :: dx
  type(panel_mesh), intent(out) :: mesh
  integer :: n1, k

  n1 = m1%n
  call alloc_mesh(m1%n + m2%n, m1%nnode + m2%nnode, mesh)
  deallocate(mesh%xr)
  allocate(mesh%xr(3,2))
  mesh%nbody = 2
  mesh%nbase = mesh%n
  mesh%xnode(:, 1:m1%nnode) = m1%xnode
  mesh%xnode(:, m1%nnode+1:) = m2%xnode
  mesh%xnode(1, m1%nnode+1:) = mesh%xnode(1, m1%nnode+1:) + dx
  mesh%xr(:,1) = m1%xr(:,1)
  mesh%xr(:,2) = m2%xr(:,1) + [dx, 0d0, 0d0]
  mesh%xc(:, 1:n1) = m1%xc
  mesh%xc(:, n1+1:) = m2%xc
  mesh%xc(1, n1+1:) = mesh%xc(1, n1+1:) + dx
  mesh%nrm(:, 1:n1) = m1%nrm
  mesh%nrm(:, n1+1:) = m2%nrm
  mesh%area(1:n1) = m1%area
  mesh%area(n1+1:) = m2%area
  mesh%psize(1:n1) = m1%psize
  mesh%psize(n1+1:) = m2%psize
  mesh%nj(:, 1:n1) = m1%nj
  mesh%nj(:, n1+1:) = m2%nj
  mesh%ncn(1:n1) = m1%ncn
  mesh%ncn(n1+1:) = m2%ncn
  mesh%ncon(:, 1:n1) = m1%ncon
  do k = 1, m2%n
    mesh%ncon(1:m2%ncn(k), n1+k) = m2%ncon(1:m2%ncn(k), k) + m1%nnode
  end do
  mesh%body(n1+1:) = 2
end subroutine join_bodies

!> Number of times a mesh node is an endpoint of a waterline segment.
integer function endpoint_count(geom, x)
  type(surface_geom), intent(in) :: geom
  real(8), intent(in) :: x(3)
  integer :: i

  endpoint_count = 0
  do i = 1, geom%nwl
    if (maxval(abs(geom%wl_x1(:,i) - x)) < 1d-10) endpoint_count = endpoint_count + 1
    if (maxval(abs(geom%wl_x2(:,i) - x)) < 1d-10) endpoint_count = endpoint_count + 1
  end do
end function endpoint_count

logical function has_nbr(geom, i, j)
  type(surface_geom), intent(in) :: geom
  integer, intent(in) :: i, j
  has_nbr = any(geom%nbr(geom%nbr_start(i):geom%nbr_start(i+1)-1) == j)
end function has_nbr

! ---------------------------------------------------------------------------------------------- tests

subroutine test_cylinder_waterline_length(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  logical :: ok
  character(len=:), allocatable :: msg
  integer, parameter :: nt = 48
  real(8) :: th
  integer :: it

  call make_cylinder(2d0, 1d0, nt, 2, 3, mesh, .true.)
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "build failed")
  if (allocated(error)) return
  call check(error, geom%nwl == nt, "segment count")
  if (allocated(error)) return
  call check(error, close_r(sum(geom%wl_len), polygon_perimeter(nt, 2d0), 1d-12), "waterline length")
  if (allocated(error)) return
  do it = 1, nt
    th = 2d0 * pi * real(it - 1, 8) / real(nt, 8)
    call check(error, endpoint_count(geom, [2d0*cos(th), 2d0*sin(th), 0d0]) == 2, "node in exactly two segments")
    if (allocated(error)) return
  end do
  call check(error, geom%n_wl_gaps(1) == 0 .and. len(msg) == 0, "closed loop has no gaps or warnings")
  if (allocated(error)) return
  ! counter-clockwise loop with positive waterplane area; closed cylinder volume
  call check(error, close_r(geom%awp(1), polygon_area(nt, 2d0), 1d-12), "waterplane area")
  if (allocated(error)) return
  call check(error, close_r(geom%vol(1), polygon_area(nt, 2d0) * 1d0, 1d-12), "volume")
  if (allocated(error)) return
  call check(error, size(geom%nbr_start) == mesh%n + 1 .and. size(geom%nbr) == geom%nbr_start(mesh%n + 1) - 1, &
             "neighbour list shape")
end subroutine test_cylinder_waterline_length

subroutine test_waterline_normal_into_body(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  logical :: ok
  character(len=:), allocatable :: msg
  integer :: i

  call make_cylinder(2d0, 1d0, 32, 2, 2, mesh, .true.)
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok .and. geom%nwl == 32, "build")
  if (allocated(error)) return
  do i = 1, geom%nwl
    call check(error, abs(norm2(geom%wl_nh(:,i)) - 1d0) < 1d-14 .and. geom%wl_nh(3,i) == 0d0, "unit horizontal")
    if (allocated(error)) return
    call check(error, dot_product([0d0, 0d0, geom%wl_xm(3,i)] - geom%wl_xm(:,i), geom%wl_nh(:,i)) > 0d0, &
               "points toward the axis")
    if (allocated(error)) return
    call check(error, geom%wl_nz(i) == 0d0 .and. geom%wl_panel(i) >= 1 .and. geom%wl_panel(i) <= mesh%n, "panel data")
    if (allocated(error)) return
  end do
end subroutine test_waterline_normal_into_body

!> A hull that is twice as wide at the waterline as at the bottom, so the wall normal has a vertical part.
subroutine test_flared_normal_is_unit(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  logical :: ok
  character(len=:), allocatable :: msg
  real(8) :: d1(3), d2(3), c(3), n(3), xh(3)
  integer :: i, k, p

  call make_cylinder(2d0, 1d0, 32, 2, 2, mesh, .false.)
  do k = 1, mesh%nnode
    mesh%xnode(1:2, k) = mesh%xnode(1:2, k) * (1d0 + 0.5d0 * mesh%xnode(3, k))
  end do
  do i = 1, mesh%n
    d1 = mesh%xnode(:, mesh%ncon(3,i)) - mesh%xnode(:, mesh%ncon(1,i))
    d2 = mesh%xnode(:, mesh%ncon(4,i)) - mesh%xnode(:, mesh%ncon(2,i))
    c = [d1(2)*d2(3) - d1(3)*d2(2), d1(3)*d2(1) - d1(1)*d2(3), d1(1)*d2(2) - d1(2)*d2(1)]
    n = c / norm2(c)
    xh = sum(mesh%xnode(:, mesh%ncon(1:4,i)), dim=2) / 4d0
    xh(3) = 0d0
    if (dot_product(n, xh) > 0d0) n = -n   ! into the body
    call set_panel(mesh, i, n)
  end do
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok .and. geom%nwl == 32, "build")
  if (allocated(error)) return
  do i = 1, geom%nwl
    p = geom%wl_panel(i)
    call check(error, abs(geom%wl_nz(i)) > 0.5d0 .and. geom%wl_nz(i) == mesh%nrm(3,p), "nz from the panel")
    if (allocated(error)) return
    call check(error, abs(norm2(geom%wl_nh(:,i)) - 1d0) < 1d-14 .and. geom%wl_nh(3,i) == 0d0, "unit horizontal")
    if (allocated(error)) return
    call check(error, abs(dot_product(geom%wl_nh(1:2,i), mesh%nrm(1:2,p)) - norm2(mesh%nrm(1:2,p))) < 1d-14, &
               "horizontal direction of the panel normal")
    if (allocated(error)) return
  end do
end subroutine test_flared_normal_is_unit

subroutine test_multi_body_loops(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: m1, m2, mesh
  type(surface_geom) :: geom
  logical :: ok
  character(len=:), allocatable :: msg

  call make_cylinder(2d0, 1d0, 48, 2, 2, m1, .true.)
  call make_cylinder(1d0, 0.5d0, 24, 1, 2, m2, .true.)
  call join_bodies(m1, m2, 10d0, mesh)
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok .and. geom%nwl == 72, "build and total segments")
  if (allocated(error)) return
  call check(error, count(geom%wl_body == 1) == 48 .and. count(geom%wl_body == 2) == 24, "segments per body")
  if (allocated(error)) return
  call check(error, close_r(sum(geom%wl_len, mask=geom%wl_body == 1), polygon_perimeter(48, 2d0), 1d-12), "body 1 length")
  if (allocated(error)) return
  call check(error, close_r(sum(geom%wl_len, mask=geom%wl_body == 2), polygon_perimeter(24, 1d0), 1d-12), "body 2 length")
  if (allocated(error)) return
  call check(error, all(geom%n_wl_gaps == 0), "no gaps")
  if (allocated(error)) return
  call check(error, close_r(geom%awp(2), polygon_area(24, 1d0), 1d-12) .and. abs(geom%xf(2)) < 1d-12, &
             "body 2 waterplane about its own xr")
  if (allocated(error)) return
  call check(error, close_r(geom%vol(2), polygon_area(24, 1d0) * 0.5d0, 1d-12), "body 2 volume")
end subroutine test_multi_body_loops

subroutine test_submerged_body_no_waterline(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  logical :: ok
  character(len=:), allocatable :: msg

  call make_cylinder(2d0, 1d0, 24, 2, 2, mesh, .true.)
  mesh%xnode(3, :) = mesh%xnode(3, :) - 3d0
  mesh%xc(3, :) = mesh%xc(3, :) - 3d0
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "build")
  if (allocated(error)) return
  call check(error, geom%nwl == 0 .and. geom%awp(1) == 0d0 .and. geom%n_wl_gaps(1) == 0, "no waterline")
end subroutine test_submerged_body_no_waterline

!> A half hull expanded by symmetry has reversed winding on the image, and a triangle there carries an offset in ncon(4).
subroutine test_symmetry_images_join(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: full, half, mesh
  type(surface_geom) :: gfull, geom
  logical :: ok
  character(len=:), allocatable :: msg
  integer, parameter :: nt = 24
  real(8), allocatable :: xc_b(:,:), dxyz_b(:,:), xnode_b(:,:)

  call make_cylinder(2d0, 1d0, nt, 2, 3, full, .true.)
  call subset_mesh(full, full%xc(2, :) > 0d0, half)
  xc_b = transpose(half%xc)
  dxyz_b = transpose(half%nj)
  xnode_b = transpose(half%xnode)
  call expand_images(xc_b, dxyz_b(:,1:3), dxyz_b, half%area, half%psize, xnode_b, half%ncn, transpose(half%ncon), half%body, &
                     half%xr, 0, 1, 2, mesh)
  call check(error, mesh%n == 2 * half%n .and. any(mesh%ncn == 3 .and. mesh%ncon(4,:) /= 0), "image triangles carry ncon(4)")
  if (allocated(error)) return
  call build_surface_geom(full, wl_tol, DEFAULT_SHARP_DEG, gfull, ok, msg)
  call check(error, ok, "build of the full hull")
  if (allocated(error)) return
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok .and. geom%nwl == nt .and. geom%n_wl_gaps(1) == 0 .and. len(msg) == 0, "joined loop, no gaps")
  if (allocated(error)) return
  call check(error, close_r(sum(geom%wl_len), polygon_perimeter(nt, 2d0), 1d-12), "perimeter of the full hull")
  if (allocated(error)) return
  call check(error, close_r(geom%awp(1), gfull%awp(1), 1d-12) .and. close_r(geom%vol(1), gfull%vol(1), 1d-12), &
             "waterplane area and volume of the full hull")
  if (allocated(error)) return
  ! the image panel across the plane is a neighbour of the base panel
  call check(error, has_nbr(geom, 1, half%n + 1), "neighbour across the symmetry plane")
end subroutine test_symmetry_images_join

!> Wall in the plane y = 0 with three panels on the left of x = 0.4 and two on the right. No panel of one side
!> shares a node with the middle panel of the other.
subroutine test_nonconforming_neighbours(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  logical :: ok
  character(len=:), allocatable :: msg
  real(8), parameter :: nrm(3) = [0d0, 1d0, 0d0]
  integer :: iz

  call alloc_mesh(5, 14, mesh)
  do iz = 0, 3   ! left nodes 1..8
    mesh%xnode(:, 1 + 2*iz) = [0d0, 0d0, -real(iz, 8)/3d0]
    mesh%xnode(:, 2 + 2*iz) = [0.4d0, 0d0, -real(iz, 8)/3d0]
  end do
  do iz = 0, 2   ! right nodes 9..14
    mesh%xnode(:, 9 + 2*iz) = [0.4d0, 0d0, -0.5d0*real(iz, 8)]
    mesh%xnode(:, 10 + 2*iz) = [0.8d0, 0d0, -0.5d0*real(iz, 8)]
  end do
  call set_quad(mesh, 1, 1, 2, 4, 3, nrm)
  call set_quad(mesh, 2, 3, 4, 6, 5, nrm)
  call set_quad(mesh, 3, 5, 6, 8, 7, nrm)
  call set_quad(mesh, 4, 9, 10, 12, 11, nrm)
  call set_quad(mesh, 5, 11, 12, 14, 13, nrm)
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "build")
  if (allocated(error)) return
  call check(error, has_nbr(geom, 2, 1) .and. has_nbr(geom, 2, 3), "conforming neighbours")
  if (allocated(error)) return
  call check(error, has_nbr(geom, 2, 4) .and. has_nbr(geom, 2, 5), "neighbours across the junction")
  if (allocated(error)) return
  call check(error, has_nbr(geom, 4, 2) .and. has_nbr(geom, 5, 2), "symmetric")
  if (allocated(error)) return
  call check(error, .not. has_nbr(geom, 2, 2), "not its own neighbour")
end subroutine test_nonconforming_neighbours

!> A 5 x 5 flat grid below the waterline, with a panel size large enough for centroids two columns apart to be in reach.
!> The panels with an open edge get them, the panels with none do not.
subroutine test_conforming_no_proximity(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: ng = 5
  real(8), parameter :: nrm(3) = [0d0, 1d0, 0d0]
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  logical :: ok
  character(len=:), allocatable :: msg
  integer :: ix, iz

  call alloc_mesh(ng*ng, (ng + 1)*(ng + 1), mesh)
  do iz = 0, ng
    do ix = 0, ng
      mesh%xnode(:, node(ix, iz)) = [real(ix, 8), 0d0, -1d0 - real(iz, 8)]
    end do
  end do
  do iz = 0, ng - 1
    do ix = 0, ng - 1
      call set_quad(mesh, panel(ix, iz), node(ix, iz), node(ix + 1, iz), node(ix + 1, iz + 1), node(ix, iz + 1), nrm)
    end do
  end do
  mesh%psize = 2d0
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "build")
  if (allocated(error)) return
  call check(error, has_nbr(geom, panel(2, 2), panel(1, 1)) .and. has_nbr(geom, panel(2, 2), panel(3, 2)), &
             "a conforming panel keeps the panels sharing a node")
  if (allocated(error)) return
  call check(error, .not. has_nbr(geom, panel(2, 2), panel(0, 2)) .and. .not. has_nbr(geom, panel(2, 2), panel(4, 2)), &
             "a conforming panel has a proximity neighbour")
  if (allocated(error)) return
  call check(error, geom%nbr_start(panel(2, 2) + 1) - geom%nbr_start(panel(2, 2)) == 8, &
             "a conforming panel of the grid has 8 neighbours")
  if (allocated(error)) return
  call check(error, has_nbr(geom, panel(0, 2), panel(2, 2)), "a panel with an open edge has lost its proximity neighbour")

contains

  integer function node(ix, iz)
    integer, intent(in) :: ix, iz
    node = 1 + ix + (ng + 1)*iz
  end function node

  integer function panel(ix, iz)
    integer, intent(in) :: ix, iz
    panel = 1 + ix + ng*iz
  end function panel

end subroutine test_conforming_no_proximity

!> A lid of 3 x 3 panels: a thin rim column gives the ratio, a repeated node is not a zero side, and no panels give 1.
subroutine test_lid_rim_ratio(error)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: xyz(16, 3)
  integer :: ncn(9), ncon(9, 4), ix, iy, ip
  real(8), parameter :: tol = 1d-12

  ncn = 4
  ip = 0
  do iy = 0, 2
    do ix = 0, 2
      ip = ip + 1
      ncon(ip, :) = [1 + ix + 4*iy, 2 + ix + 4*iy, 6 + ix + 4*iy, 5 + ix + 4*iy]
    end do
  end do

  call fill([0d0, 1d0, 2d0, 3d0], [0d0, 1d0, 2d0, 3d0])
  call check(error, abs(lid_rim_side_ratio(xyz, ncn, ncon, 0, 0) - 1d0) < tol, "uniform lid")
  if (allocated(error)) return

  call fill([0d0, 0.2d0, 1.2d0, 2.2d0], [0d0, 1d0, 2d0, 3d0])
  call check(error, abs(lid_rim_side_ratio(xyz, ncn, ncon, 0, 0) - 0.2d0) < tol, "thin rim column")
  if (allocated(error)) return

  call fill([0d0, 1d0, 2d0, 3d0], [0d0, 1d0, 2d0, 3d0])
  ncon(1, 4) = ncon(1, 3)
  call check(error, lid_rim_side_ratio(xyz, ncn, ncon, 0, 0) > 0.5d0, "a repeated node")
  if (allocated(error)) return

  call check(error, lid_rim_side_ratio(xyz, ncn(1:0), ncon(1:0, :), 0, 0) == 1d0, "no panels")

contains

  subroutine fill(x, y)
    real(8), intent(in) :: x(4), y(4)
    integer :: i, j

    do j = 0, 3
      do i = 0, 3
        xyz(1 + i + 4*j, :) = [x(i + 1), y(j + 1), 0d0]
      end do
    end do
  end subroutine fill

end subroutine test_lid_rim_ratio

!> A disc of uniform rings as the waterplane mesh arrays, over the angles first to first + span. The fans at the centre are
!> triangles, or with repeated quadrilaterals written with their last node twice.
subroutine make_disc(ntheta, nr, first, span, repeated, xyz, ncn, ncon)
  integer, intent(in) :: ntheta, nr
  real(8), intent(in) :: first, span
  logical, intent(in) :: repeated
  real(8), allocatable, intent(out) :: xyz(:, :)
  integer, allocatable, intent(out) :: ncn(:), ncon(:, :)
  integer :: i, j, ip, closed

  ! Rings 0 (the rim) to nr - 1, then the centre. A whole disc closes on its first column of nodes.
  closed = merge(ntheta, ntheta + 1, abs(span - 2d0*pi) < 1d-12)
  allocate (xyz(closed*nr + 1, 3), ncn(ntheta*nr), ncon(ntheta*nr, 4))
  do i = 0, nr - 1
    do j = 0, closed - 1
      xyz(node(i, j), :) = 10d0*real(nr - i, 8)/real(nr, 8)* &
                           [cos(first + span*real(j, 8)/real(ntheta, 8)), sin(first + span*real(j, 8)/real(ntheta, 8)), 0d0]
    end do
  end do
  xyz(closed*nr + 1, :) = 0d0
  ip = 0
  do j = 0, ntheta - 1
    do i = 0, nr - 2
      ip = ip + 1
      ncn(ip) = 4
      ncon(ip, :) = [node(i + 1, j + 1), node(i + 1, j), node(i, j), node(i, j + 1)]
    end do
    ip = ip + 1
    if (repeated) then
      ncn(ip) = 4
      ncon(ip, :) = [closed*nr + 1, node(nr - 1, j), node(nr - 1, j + 1), node(nr - 1, j + 1)]
    else
      ncn(ip) = 3
      ncon(ip, :) = [closed*nr + 1, node(nr - 1, j), node(nr - 1, j + 1), 0]
    end if
  end do

contains

  integer function node(ring, column)
    integer, intent(in) :: ring, column
    node = ring*closed + mod(column, closed) + 1
  end function node

end subroutine make_disc

!> The edges on a symmetry plane continue on the image of the mesh, so they are not rim: half and quarter discs have the
!> ratio of the whole disc, with the plane named by ISX or ISY as in HAMS. Without the flag, those edges do count.
subroutine test_lid_rim_symmetry_planes(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), allocatable :: xyz(:, :)
  integer, allocatable :: ncn(:), ncon(:, :)
  real(8) :: whole
  real(8), parameter :: tol = 1d-9

  call make_disc(32, 6, 0d0, 2d0*pi, .false., xyz, ncn, ncon)
  whole = lid_rim_side_ratio(xyz, ncn, ncon, 0, 0)
  call check(error, whole > THIN_ENOUGH, "a uniform whole disc")
  if (allocated(error)) return

  call make_disc(16, 6, 0d0, pi, .false., xyz, ncn, ncon)
  call check(error, abs(lid_rim_side_ratio(xyz, ncn, ncon, 0, 1) - whole) < tol, "half disc, plane y = 0, ISY = 1")
  if (allocated(error)) return
  call check(error, lid_rim_side_ratio(xyz, ncn, ncon, 0, 0) < whole - tol, "half disc without the symmetry flag")
  if (allocated(error)) return
  call check(error, lid_rim_side_ratio(xyz, ncn, ncon, 1, 0) < whole - tol, "half disc, ISX = 1 names the other plane")
  if (allocated(error)) return

  call make_disc(16, 6, -0.5d0*pi, pi, .false., xyz, ncn, ncon)
  call check(error, abs(lid_rim_side_ratio(xyz, ncn, ncon, 1, 0) - whole) < tol, "half disc, plane x = 0, ISX = 1")
end subroutine test_lid_rim_symmetry_planes

!> A triangle written as a quadrilateral has an edge from a node to itself, which belongs to no other panel but is no rim:
!> the thin fans at the centre must not count.
subroutine test_lid_rim_repeated_node_fans(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), allocatable :: xyz(:, :)
  integer, allocatable :: ncn(:), ncon(:, :)
  real(8) :: triangles, repeated

  call make_disc(96, 18, 0d0, 2d0*pi, .false., xyz, ncn, ncon)
  triangles = lid_rim_side_ratio(xyz, ncn, ncon, 0, 0)
  call make_disc(96, 18, 0d0, 2d0*pi, .true., xyz, ncn, ncon)
  repeated = lid_rim_side_ratio(xyz, ncn, ncon, 0, 0)
  call check(error, abs(repeated - triangles) < 1d-12, "fans of repeated nodes change the ratio")
end subroutine test_lid_rim_repeated_node_fans

!> Thin across the rim draws the ratio down. Panels long across the rim, narrow along it, do not: see the README.
subroutine test_lid_rim_only_across_rim(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), allocatable :: xyz(:, :)
  integer, allocatable :: ncn(:), ncon(:, :)
  real(8), parameter :: tol = 1d-9

  ! 96 x 3: rim arc 0.654, ring width 3.33
  call make_disc(96, 3, 0d0, 2d0*pi, .false., xyz, ncn, ncon)
  call check(error, abs(lid_rim_side_ratio(xyz, ncn, ncon, 0, 0) - 1d0) < tol, "long across the rim")
  if (allocated(error)) return
  ! 32 x 12: rim arc 1.96, ring width 0.833
  call make_disc(32, 12, 0d0, 2d0*pi, .false., xyz, ncn, ncon)
  call check(error, abs(lid_rim_side_ratio(xyz, ncn, ncon, 0, 0) - (10d0/12d0)/(20d0*sin(pi/32d0))) < tol, "thin across the rim")
end subroutine test_lid_rim_only_across_rim

subroutine test_sharp_edge_excluded(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom, geom_all
  logical :: ok
  character(len=:), allocatable :: msg
  integer :: i, j, k, nside
  real(8) :: cosmin

  call make_box(0d0, 0d0, 2d0, 2d0, 1d0, 2, 2, 2, [0d0, 0d0, 0d0], mesh)
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "build")
  if (allocated(error)) return
  cosmin = cos(DEFAULT_SHARP_DEG * pi / 180d0)
  nside = 0
  do i = 1, mesh%n
    do k = geom%nbr_start(i), geom%nbr_start(i+1) - 1
      j = geom%nbr(k)
      call check(error, dot_product(mesh%nrm(:,i), mesh%nrm(:,j)) >= cosmin - 1d-12, "neighbour beyond the sharp angle")
      if (allocated(error)) return
    end do
    if (mesh%nrm(1,i) == 1d0) nside = nside + geom%nbr_start(i+1) - geom%nbr_start(i)
  end do
  ! panels of the x = 0 face see the other panels of that face (3 each) and nothing else
  call check(error, nside == 4 * 3, "side panels see only their own face")
  if (allocated(error)) return
  ! with the angle limit lifted the bottom and the other walls come back
  call build_surface_geom(mesh, wl_tol, 180d0, geom_all, ok, msg)
  call check(error, size(geom_all%nbr) > size(geom%nbr), "larger angle admits more neighbours")
end subroutine test_sharp_edge_excluded

subroutine test_wp_properties_box(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  logical :: ok
  character(len=:), allocatable :: msg
  real(8), parameter :: x0 = 1d0, y0 = -2d0, a = 3d0, b = 2d0, d = 1.5d0
  real(8), parameter :: xr(3) = [0.3d0, -0.2d0, 0.1d0]
  real(8) :: u0, u1, v0, v1

  call make_box(x0, y0, a, b, d, 2, 3, 2, xr, mesh)
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok .and. geom%nwl == 2*2 + 2*3, "build")
  if (allocated(error)) return
  u0 = x0 - xr(1); u1 = u0 + a
  v0 = y0 - xr(2); v1 = v0 + b
  call check(error, close_r(geom%awp(1), a*b, 1d-12), "awp")
  if (allocated(error)) return
  call check(error, close_r(geom%xf(1), 0.5d0*(u0 + u1), 1d-12) .and. close_r(geom%yf(1), 0.5d0*(v0 + v1), 1d-12), "xf, yf")
  if (allocated(error)) return
  call check(error, close_r(geom%l11(1), b * (u1**3 - u0**3) / 3d0, 1d-12), "l11")
  if (allocated(error)) return
  call check(error, close_r(geom%l22(1), a * (v1**3 - v0**3) / 3d0, 1d-12), "l22")
  if (allocated(error)) return
  call check(error, close_r(geom%l12(1), 0.5d0*(u1**2 - u0**2) * 0.5d0*(v1**2 - v0**2), 1d-12), "l12")
  if (allocated(error)) return
  call check(error, close_r(geom%vol(1), a*b*d, 1d-12), "vol")
  if (allocated(error)) return
  call check(error, close_r(geom%xb(1), 0.5d0*(u0 + u1), 1d-12) .and. close_r(geom%yb(1), 0.5d0*(v0 + v1), 1d-12) .and. &
             close_r(geom%zb(1), -0.5d0*d - xr(3), 1d-12), "xb, yb, zb")
end subroutine test_wp_properties_box

subroutine test_lid_panels_rejected(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  logical :: ok
  character(len=:), allocatable :: msg

  call alloc_mesh(2, 8, mesh)
  mesh%xnode(:, 1) = [0d0, 0d0, 0d0]
  mesh%xnode(:, 2) = [0d0, 0d0, -1d0]
  mesh%xnode(:, 3) = [0d0, 1d0, -1d0]
  mesh%xnode(:, 4) = [0d0, 1d0, 0d0]
  mesh%xnode(:, 5) = [1d0, 0d0, 0d0]
  mesh%xnode(:, 6) = [2d0, 0d0, 0d0]
  mesh%xnode(:, 7) = [2d0, 1d0, 0d0]
  mesh%xnode(:, 8) = [1d0, 1d0, 0d0]
  call set_quad(mesh, 1, 1, 2, 3, 4, [1d0, 0d0, 0d0])
  call set_quad(mesh, 2, 5, 6, 7, 8, [0d0, 0d0, 1d0])
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, .not. ok, "a lid panel must stop the run")
  if (allocated(error)) return
  call check(error, geom%nwl == 0 .and. .not. allocated(geom%wl_x1) .and. .not. allocated(geom%awp) .and. &
             .not. allocated(geom%nbr) .and. .not. allocated(geom%n_wl_gaps), "geom is empty after a stop")
  if (allocated(error)) return
  call check(error, index(msg, "panel 2 ") > 0 .and. index(msg, "lid") > 0, "message names the panel")
end subroutine test_lid_panels_rejected

!> A waterline edge on a panel whose normal is vertical has no horizontal direction.
subroutine test_horizontal_normal_rejected(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  logical :: ok
  character(len=:), allocatable :: msg

  call alloc_mesh(1, 3, mesh)
  mesh%xnode(:, 1) = [0d0, 0d0, 0d0]
  mesh%xnode(:, 2) = [1d0, 0d0, 0d0]
  mesh%xnode(:, 3) = [0d0, 0.5d0, -0.3d0]
  mesh%ncn(1) = 3
  mesh%ncon(1:3, 1) = [1, 2, 3]
  call set_panel(mesh, 1, [0d0, 0d0, 1d0])
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, .not. ok, "vertical normal on the waterline must stop the run")
  if (allocated(error)) return
  call check(error, geom%nwl == 0 .and. .not. allocated(geom%wl_x1) .and. .not. allocated(geom%awp) .and. &
             .not. allocated(geom%nbr) .and. .not. allocated(geom%n_wl_gaps), "geom is empty after a stop")
  if (allocated(error)) return
  call check(error, index(msg, "panel 1 ") > 0 .and. index(msg, "horizontal") > 0, "message names the panel")
end subroutine test_horizontal_normal_rejected

subroutine test_waterline_off_zero(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  logical :: ok
  character(len=:), allocatable :: msg
  integer, parameter :: nt = 24
  real(8), parameter :: size_ = 4d0   ! largest dimension of the hull: the diameter

  ! top nodes just below z = 0 are still on the waterline
  call make_cylinder(2d0, 1d0, nt, 2, 2, mesh, .true.)
  mesh%xnode(3, 1:nt) = -1d-8 * size_
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok .and. geom%nwl == nt .and. len(msg) == 0, "1e-8 below z = 0 is on the waterline")
  if (allocated(error)) return

  ! 1 percent of the body size below z = 0 is not: an open hull that does not reach the surface
  call make_cylinder(2d0, 1d0, nt, 2, 2, mesh, .true.)
  mesh%xnode(3, 1:nt) = -1d-2 * size_
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok .and. geom%nwl == 0, "no segments, but not a stop")
  if (allocated(error)) return
  call check(error, index(msg, "no waterline") > 0, "warns that no waterline was found")
end subroutine test_waterline_off_zero

subroutine test_waterline_gap_reported(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: full, mesh
  type(surface_geom) :: geom
  logical :: ok, keep(24)
  character(len=:), allocatable :: msg

  call make_cylinder(2d0, 1d0, 24, 1, 2, full, .false.)
  keep = .true.
  keep(5) = .false.   ! removes one waterline segment
  call subset_mesh(full, keep, mesh)
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "a gap is a warning, not a stop")
  if (allocated(error)) return
  call check(error, geom%nwl == 23 .and. geom%n_wl_gaps(1) == 1, "one gap")
  if (allocated(error)) return
  call check(error, index(msg, "gap") > 0, "message warns about the gap")
end subroutine test_waterline_gap_reported

!> Top row of the wall split into triangles (ncon(4) = 0), the bottom keeps its triangles and quads.
subroutine test_triangles_supported(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: cyl, mesh
  type(surface_geom) :: geom
  logical :: ok
  character(len=:), allocatable :: msg
  integer, parameter :: nt = 32
  integer :: ip, q(4), i
  real(8) :: nrm(3)

  call make_cylinder(2d0, 1d0, nt, 2, 3, cyl, .true.)
  call alloc_mesh(cyl%n + nt, cyl%nnode, mesh)
  mesh%xnode = cyl%xnode
  mesh%xr = cyl%xr
  mesh%xc(:, 1:cyl%n) = cyl%xc
  mesh%nrm(:, 1:cyl%n) = cyl%nrm
  mesh%area(1:cyl%n) = cyl%area
  mesh%psize(1:cyl%n) = cyl%psize
  mesh%nj(:, 1:cyl%n) = cyl%nj
  mesh%ncn(1:cyl%n) = cyl%ncn
  mesh%ncon(:, 1:cyl%n) = cyl%ncon
  do ip = 1, nt   ! quad (q1 top, q2, q3, q4 top) becomes (q1, q3, q4) and (q1, q2, q3)
    q = cyl%ncon(:, ip)
    nrm = cyl%nrm(:, ip)
    mesh%ncn(ip) = 3
    mesh%ncon(:, ip) = [q(1), q(3), q(4), 0]
    call set_panel(mesh, ip, nrm)
    mesh%ncn(cyl%n + ip) = 3
    mesh%ncon(:, cyl%n + ip) = [q(1), q(2), q(3), 0]
    call set_panel(mesh, cyl%n + ip, nrm)
  end do
  call build_surface_geom(mesh, wl_tol, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok .and. geom%nwl == nt .and. geom%n_wl_gaps(1) == 0, "build")
  if (allocated(error)) return
  call check(error, close_r(sum(geom%wl_len), polygon_perimeter(nt, 2d0), 1d-12), "waterline length")
  if (allocated(error)) return
  do i = 1, geom%nwl
    call check(error, mesh%ncn(geom%wl_panel(i)) == 3, "waterline on a triangle")
    if (allocated(error)) return
  end do
  call check(error, close_r(geom%vol(1), polygon_area(nt, 2d0), 1d-12), "volume with triangles")
end subroutine test_triangles_supported

end module test_surface_geometry
