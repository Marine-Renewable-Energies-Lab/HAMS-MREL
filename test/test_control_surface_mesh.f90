!> Unit tests for the generated control surfaces C, the line Gamma_c and the checks against the other bodies.
module test_control_surface_mesh
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderTypes, only : panel_mesh, surface_geom, so_options, cs_points, CS_WALL, CS_LINE, CS_FREE
  use SurfaceGeometry, only : build_surface_geom, DEFAULT_SHARP_DEG
  use ControlSurfaceMesh, only : cylinder_surface, box_surface, add_surface_points, build_control_surface, &
                                 check_control_surface, write_control_surface, read_control_mesh, &
                                 free_surface_strips
  use second_order_fixtures, only : make_cylinder, set_panel
  implicit none
  private

  public :: collect_control_surface_mesh

  real(8), parameter :: pi = 3.141592653589793d0

contains

subroutine collect_control_surface_mesh(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("cylinder_area_and_normals", test_cylinder_area_and_normals), &
    new_unittest("cylinder_to_seabed", test_cylinder_to_seabed), &
    new_unittest("box_area_and_normals", test_box_area_and_normals), &
    new_unittest("gamma_c_line", test_gamma_c_line), &
    new_unittest("generated_c_rejects_other_body", test_generated_c_rejects_other_body), &
    new_unittest("separate_bodies_accepted", test_separate_bodies_accepted), &
    new_unittest("close_c_rejected_by_distance", test_close_c_rejected_by_distance), &
    new_unittest("box_built_for_each_body", test_box_built_for_each_body), &
    new_unittest("file_c_missing_file_stops", test_file_c_missing_file_stops), &
    new_unittest("file_c_round_trip", test_file_c_round_trip), &
    new_unittest("file_c_wrong_orientation_stops", test_file_c_wrong_orientation_stops), &
    new_unittest("written_mesh_has_hull_format", test_written_mesh_has_hull_format), &
    new_unittest("strips_area", test_strips_area), &
    new_unittest("strips_box_hull", test_strips_box_hull), &
    new_unittest("not_star_shaped_stops", test_not_star_shaped_stops), &
    new_unittest("with_free_appends_f", test_with_free_appends_f), &
    new_unittest("fs_rows_default", test_fs_rows_default), &
    new_unittest("fs_file_takes_precedence", test_fs_file_takes_precedence) &
    ]
end subroutine collect_control_surface_mesh

! ---------------------------------------------------------------------------------------------- helpers

!> Wall area (panels with |nz| < 0.5) and bottom area (nz > 0.5) of a control surface.
subroutine areas(cm, awall, abottom)
  type(panel_mesh), intent(in) :: cm
  real(8), intent(out) :: awall, abottom
  integer :: i
  awall = 0d0
  abottom = 0d0
  do i = 1, cm%n
    if (abs(cm%nrm(3, i)) < 0.5d0) then
      awall = awall + cm%area(i)
    else
      abottom = abottom + cm%area(i)
    end if
  end do
end subroutine areas

!> Two copies of the same hull, the second shifted by dx along x, as bodies 1 and 2 of one mesh.
subroutine two_bodies(m1, dx, mesh)
  type(panel_mesh), intent(in) :: m1
  real(8), intent(in) :: dx
  type(panel_mesh), intent(out) :: mesh
  integer :: n, nn, i

  n = m1%n
  nn = m1%nnode
  mesh%n = 2*n; mesh%nbase = 2*n; mesh%nnode = 2*nn; mesh%nbody = 2; mesh%nsys = 1
  allocate (mesh%xc(3, 2*n), mesh%nrm(3, 2*n), mesh%area(2*n), mesh%psize(2*n), mesh%nj(6, 2*n))
  allocate (mesh%xnode(3, 2*nn), mesh%ncn(2*n), mesh%ncon(4, 2*n), mesh%body(2*n), mesh%base(2*n), mesh%image(2*n))
  allocate (mesh%xr(3, 2))
  mesh%xc(:, 1:n) = m1%xc;  mesh%xc(:, n+1:) = m1%xc;  mesh%xc(1, n+1:) = m1%xc(1, :) + dx
  mesh%nrm(:, 1:n) = m1%nrm;  mesh%nrm(:, n+1:) = m1%nrm
  mesh%area(1:n) = m1%area;  mesh%area(n+1:) = m1%area
  mesh%psize(1:n) = m1%psize;  mesh%psize(n+1:) = m1%psize
  mesh%nj(:, 1:n) = m1%nj;  mesh%nj(:, n+1:) = m1%nj
  mesh%xnode(:, 1:nn) = m1%xnode;  mesh%xnode(:, nn+1:) = m1%xnode;  mesh%xnode(1, nn+1:) = m1%xnode(1, :) + dx
  mesh%ncn(1:n) = m1%ncn;  mesh%ncn(n+1:) = m1%ncn
  mesh%ncon(:, 1:n) = m1%ncon
  mesh%ncon(:, n+1:) = merge(m1%ncon + nn, 0, m1%ncon > 0)
  mesh%body(1:n) = 1;  mesh%body(n+1:) = 2
  mesh%base = [(i, i = 1, 2*n)]
  mesh%image = 1
  mesh%xr = 0d0
  mesh%xr(1, 2) = dx
end subroutine two_bodies

!> Builds C for every body of mesh with the given options and checks it.
subroutine build_and_check(mesh, opts, generated, pts, cms, ok_build, ok_check, msg)
  type(panel_mesh), intent(in) :: mesh
  type(so_options), intent(in) :: opts
  logical, intent(in) :: generated
  type(cs_points), intent(out) :: pts
  type(panel_mesh), allocatable, intent(out) :: cms(:)
  logical, intent(out) :: ok_build, ok_check
  character(len=:), allocatable, intent(out) :: msg
  type(surface_geom) :: geom
  character(len=:), allocatable :: msg0

  ok_check = .false.
  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok_build, msg0)
  if (.not. ok_build) then
    msg = msg0
    return
  end if
  call build_control_surface(mesh, geom, opts, -1d0, '.', pts, cms, ok_build, msg)
  if (.not. ok_build) return
  call check_control_surface(mesh, pts, cms, generated, ok_check, msg)
end subroutine build_and_check

! ---------------------------------------------------------------------------------------------- tests

subroutine test_cylinder_area_and_normals(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: cm
  real(8) :: aw, ab
  integer :: i

  call cylinder_surface([0d0, 0d0], 20d0, 15d0, 2d0, 100d0, cm)
  call check(error, cm%n > 0 .and. cm%nsys == 1 .and. cm%nbase == cm%n .and. cm%nbody == 1, "mesh bookkeeping")
  if (allocated(error)) return
  call areas(cm, aw, ab)
  call check(error, abs(aw/(2d0*pi*20d0*15d0) - 1d0) < 5d-3, "wall area")
  if (allocated(error)) return
  call check(error, abs(ab/(pi*20d0**2) - 1d0) < 5d-3, "bottom area")
  if (allocated(error)) return
  do i = 1, cm%n
    if (abs(cm%nrm(3, i)) < 0.5d0) then
      call check(error, dot_product(cm%nrm(:, i), cm%xc(:, i)) < 0d0, "wall normal towards the axis")
    else
      call check(error, cm%nrm(3, i) > 0.99d0, "bottom normal up")
    end if
    if (allocated(error)) return
    call check(error, abs(norm2(cm%nrm(:, i)) - 1d0) < 1d-12 .and. cm%area(i) > 0d0 .and. cm%psize(i) > 0d0, "panel data")
    if (allocated(error)) return
  end do
  ! the top wall row is about half a panel size thick
  call check(error, minval(abs(cm%xnode(3, :)), mask=abs(cm%xnode(3, :)) > 1d-12) < 1.5d0, "graded top row")
end subroutine test_cylinder_area_and_normals

subroutine test_cylinder_to_seabed(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: cm

  call cylinder_surface([1d0, 2d0], 20d0, 12d0, 2d0, 10d0, cm)
  call check(error, all(cm%nrm(3, :) <= 0.5d0), "no bottom panel")
  if (allocated(error)) return
  call check(error, abs(minval(cm%xnode(3, :)) + 10d0) < 1d-12, "wall reaches the seabed")
  if (allocated(error)) return
  call check(error, maxval(cm%xnode(3, :)) == 0d0, "wall reaches the free surface")
end subroutine test_cylinder_to_seabed

subroutine test_box_area_and_normals(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: cm
  real(8) :: aw, ab
  integer :: i

  call box_surface([5d0, -3d0], [12d0, 7d0], 0.3d0, 8d0, 1d0, 50d0, cm)
  call areas(cm, aw, ab)
  call check(error, abs(aw/(2d0*(24d0 + 14d0)*8d0) - 1d0) < 5d-3, "wall area")
  if (allocated(error)) return
  call check(error, abs(ab/(24d0*14d0) - 1d0) < 5d-3, "bottom area")
  if (allocated(error)) return
  do i = 1, cm%n
    if (abs(cm%nrm(3, i)) < 0.5d0) then
      call check(error, dot_product(cm%nrm(1:2, i), cm%xc(1:2, i) - [5d0, -3d0]) < 0d0, "wall normal towards the centre")
    else
      call check(error, cm%nrm(3, i) > 0.99d0, "bottom normal up")
    end if
    if (allocated(error)) return
  end do
  ! the box is rotated: every node is inside the circle through the corners
  call check(error, maxval(hypot(cm%xnode(1, :) - 5d0, cm%xnode(2, :) + 3d0)) < hypot(12d0, 7d0) + 1d-9, "extent")
end subroutine test_box_area_and_normals

subroutine test_gamma_c_line(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: cm
  type(cs_points) :: pts
  real(8) :: len_sum
  integer :: i, nl, nw

  call cylinder_surface([0d0, 0d0], 20d0, 15d0, 2d0, 100d0, cm)
  call add_surface_points(cm, 3, pts)
  nw = count(pts%kind == CS_WALL)
  call check(error, nw == cm%n .and. pts%n > nw, "wall points are the panels")
  if (allocated(error)) return
  nl = 0
  len_sum = 0d0
  do i = 1, pts%n
    call check(error, pts%body(i) == 3, "body number")
    if (allocated(error)) return
    if (pts%kind(i) /= CS_LINE) cycle
    nl = nl + 1
    len_sum = len_sum + pts%w(i)
    call check(error, abs(pts%x(3, i)) < 1d-12, "on the free surface")
    if (allocated(error)) return
    call check(error, abs(pts%nrm(3, i)) < 1d-12 .and. abs(norm2(pts%nrm(:, i)) - 1d0) < 1d-12, "horizontal unit normal")
    if (allocated(error)) return
    call check(error, dot_product(pts%nrm(:, i), pts%x(:, i)) < 0d0, "into F, towards the axis")
    if (allocated(error)) return
  end do
  call check(error, nl > 0 .and. abs(len_sum/(2d0*pi*20d0) - 1d0) < 5d-3, "Gamma_c length")
end subroutine test_gamma_c_line

subroutine test_generated_c_rejects_other_body(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: m1, mesh
  type(so_options) :: opts
  type(cs_points) :: pts
  type(panel_mesh), allocatable :: cms(:)
  logical :: okb, okc
  character(len=:), allocatable :: msg

  call make_cylinder(5d0, 3d0, 16, 2, 2, m1, .true.)
  call two_bodies(m1, 12d0, mesh)
  opts%cs_margin = 10d0
  call build_and_check(mesh, opts, .true., pts, cms, okb, okc, msg)
  call check(error, okb, "build")
  if (allocated(error)) return
  call check(error, .not. okc, "overlap must be rejected")
  if (allocated(error)) return
  call check(error, index(msg, 'body 1') > 0 .and. index(msg, 'body 2') > 0, "message names both bodies")
end subroutine test_generated_c_rejects_other_body

subroutine test_separate_bodies_accepted(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: m1, mesh
  type(so_options) :: opts
  type(cs_points) :: pts
  type(panel_mesh), allocatable :: cms(:)
  logical :: okb, okc
  character(len=:), allocatable :: msg

  call make_cylinder(5d0, 3d0, 16, 2, 2, m1, .true.)
  call two_bodies(m1, 40d0, mesh)
  opts%cs_margin = 6d0
  call build_and_check(mesh, opts, .true., pts, cms, okb, okc, msg)
  call check(error, okb .and. okc, "well separated bodies pass: "//msg)
  if (allocated(error)) return
  call check(error, size(cms) == 2 .and. any(pts%body == 1) .and. any(pts%body == 2), "one C per body")
  if (allocated(error)) return
  ! C of body 2 is centred on body 2
  call check(error, abs(sum(cms(2)%xnode(1, :))/cms(2)%nnode - 40d0) < 1d0, "C of body 2 around body 2")
end subroutine test_separate_bodies_accepted

subroutine test_close_c_rejected_by_distance(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(so_options) :: opts
  type(cs_points) :: pts
  type(panel_mesh), allocatable :: cms(:)
  logical :: okb, okc
  character(len=:), allocatable :: msg

  call make_cylinder(5d0, 3d0, 16, 2, 2, mesh, .true.)
  opts%cs_margin = 0.1d0          ! C wall almost touches the hull
  call build_and_check(mesh, opts, .true., pts, cms, okb, okc, msg)
  call check(error, okb, "build")
  if (allocated(error)) return
  call check(error, .not. okc, "C within a panel size of the hull must be rejected")
  if (allocated(error)) return
  call check(error, index(msg, 'panel size') > 0, "message explains the distance rule")
end subroutine test_close_c_rejected_by_distance

subroutine test_box_built_for_each_body(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: m1, mesh
  type(so_options) :: opts
  type(cs_points) :: pts
  type(panel_mesh), allocatable :: cms(:)
  logical :: okb, okc
  character(len=:), allocatable :: msg

  call make_cylinder(5d0, 3d0, 16, 2, 2, m1, .true.)
  call two_bodies(m1, 40d0, mesh)
  opts%cs_shape = 2
  opts%cs_margin = 6d0
  opts%cs_bottom_margin = 3d0
  call build_and_check(mesh, opts, .true., pts, cms, okb, okc, msg)
  call check(error, okb .and. okc, "box C passes: "//msg)
  if (allocated(error)) return
  call check(error, abs(minval(cms(1)%xnode(3, :)) + 6d0) < 1d-12, "depth is the deepest node plus the bottom margin")
  if (allocated(error)) return
  call check(error, maxval(cms(1)%xnode(1, :)) > 10.9d0 .and. maxval(cms(1)%xnode(1, :)) < 11.1d0, "half length plus margin")
end subroutine test_box_built_for_each_body

subroutine test_file_c_missing_file_stops(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(so_options) :: opts
  type(cs_points) :: pts
  type(panel_mesh), allocatable :: cms(:)
  logical :: okb, okc
  character(len=:), allocatable :: msg

  call make_cylinder(5d0, 3d0, 16, 2, 2, mesh, .true.)
  opts%cs_shape = 3
  call build_and_check(mesh, opts, .false., pts, cms, okb, okc, msg)
  call check(error, .not. okb, "missing file stops")
  if (allocated(error)) return
  call check(error, msg == 'cs_shape = file needs CSMesh_1.pnl in the input directory.', "message: "//msg)
end subroutine test_file_c_missing_file_stops

!> A generated cylinder C, written and read back as a file C: nodes and areas agree, and Gamma_c is found.
subroutine test_file_c_round_trip(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, cm0(1)
  type(panel_mesh), allocatable :: cms(:)
  type(so_options) :: opts
  type(cs_points) :: pts
  logical :: okb, okc
  character(len=:), allocatable :: msg
  integer :: i

  call make_cylinder(5d0, 3d0, 16, 2, 2, mesh, .true.)
  call cylinder_surface([0d0, 0d0], 11d0, 8d0, 2d0, -1d0, cm0(1))
  call write_control_surface('.', cm0)
  opts%cs_shape = 3
  call build_and_check(mesh, opts, .false., pts, cms, okb, okc, msg)
  call check(error, okb, "file C read: "//msg)
  if (allocated(error)) then
    call remove_cs_file()
    return
  end if
  call check(error, okc, "distance check passes: "//msg)
  if (.not. allocated(error)) call check(error, cms(1)%n == cm0(1)%n .and. cms(1)%nnode == cm0(1)%nnode, "counts")
  if (.not. allocated(error)) then
    ! the file has 12 digits after the point
    call check(error, maxval(abs(cms(1)%xnode - cm0(1)%xnode)) < 1d-9, "node coordinates")
  end if
  if (.not. allocated(error)) then
    do i = 1, cm0(1)%n
      if (abs(cms(1)%area(i) - cm0(1)%area(i)) > 1d-7*cm0(1)%area(i)) then
        call check(error, .false., "panel area")
        exit
      end if
    end do
  end if
  if (.not. allocated(error)) call check(error, any(pts%kind == CS_LINE) .and. count(pts%kind == CS_WALL) == cm0(1)%n, &
                                         "points of the file C")
  call remove_cs_file()
end subroutine test_file_c_round_trip

!> The same file with every panel's node order reversed has normals pointing away from the body.
subroutine test_file_c_wrong_orientation_stops(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, cm0(1)
  type(panel_mesh), allocatable :: cms(:)
  type(so_options) :: opts
  type(cs_points) :: pts
  logical :: okb, okc
  character(len=:), allocatable :: msg
  integer :: i, nv

  call make_cylinder(5d0, 3d0, 16, 2, 2, mesh, .true.)
  call cylinder_surface([0d0, 0d0], 11d0, 8d0, 2d0, -1d0, cm0(1))
  do i = 1, cm0(1)%n
    nv = cm0(1)%ncn(i)
    cm0(1)%ncon(1:nv, i) = cm0(1)%ncon(nv:1:-1, i)
  end do
  call write_control_surface('.', cm0)
  opts%cs_shape = 3
  call build_and_check(mesh, opts, .false., pts, cms, okb, okc, msg)
  call remove_cs_file()
  call check(error, .not. okb, "reversed normals stop")
  if (allocated(error)) return
  call check(error, msg == 'CSMesh_1.pnl: the panel normals must point towards the body.', "message: "//msg)
end subroutine test_file_c_wrong_orientation_stops

subroutine remove_cs_file()
  integer :: u, ios
  open (newunit=u, file='./CSMesh_1.pnl', status='old', iostat=ios)
  if (ios == 0) close (u, status='delete')
end subroutine remove_cs_file

subroutine test_written_mesh_has_hull_format(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: cms(1)
  integer :: u, ios, np, nn, isx, isy, k, idx, nv, v(4)
  real(8) :: x(3)

  call box_surface([0d0, 0d0], [6d0, 4d0], 0d0, 4d0, 2d0, 100d0, cms(1))
  call write_control_surface('.', cms)
  open (newunit=u, file='./CSMesh_1.pnl', status='old', action='read', iostat=ios)
  call check(error, ios == 0, "file written")
  if (allocated(error)) return
  ! the same reading sequence as HAMS_Prog uses for HullMesh.pnl
  do k = 1, 3
    read (u, *)
  end do
  read (u, *) np, nn, isx, isy
  call check(error, np == cms(1)%n .and. nn == cms(1)%nnode .and. isx == 0 .and. isy == 0, "header counts")
  if (.not. allocated(error)) then
    read (u, *)
    read (u, *)
    do k = 1, nn
      read (u, *) idx, x
      if (abs(x(1) - cms(1)%xnode(1, k)) > 1d-5 .or. abs(x(3) - cms(1)%xnode(3, k)) > 1d-5 .or. idx /= k) then
        call check(error, .false., "node record")
        exit
      end if
    end do
  end if
  if (.not. allocated(error)) then
    do k = 1, 3
      read (u, *)
    end do
    do k = 1, np
      v = 0
      read (u, *) idx, nv, v(1:cms(1)%ncn(k))
      if (idx /= k .or. nv /= cms(1)%ncn(k) .or. any(v /= cms(1)%ncon(:, k))) then
        call check(error, .false., "panel record")
        exit
      end if
    end do
  end if
  close (u, status='delete')
end subroutine test_written_mesh_has_hull_format

! ---------------------------------------------------------------------------------------------- free surface F

!> Open-topped box hull x in [-a/2, a/2], y in [-b/2, b/2], z in [-d, 0], centred on the origin, with normals into
!> the body (the layout of make_box in test_surface_geometry).
subroutine make_box_hull(a, b, d, nx, ny, nz, mesh)
  real(8), intent(in) :: a, b, d
  integer, intent(in) :: nx, ny, nz
  type(panel_mesh), intent(out) :: mesh
  integer :: ix, iy, iz, ip, n, nn, i

  n = nx*ny + 2*ny*nz + 2*nx*nz
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
                                         -d*real(iz, 8)/real(nz, 8)]
      end do
    end do
  end do
  ip = 0
  do iy = 0, ny - 1
    do ix = 0, nx - 1
      call quad(id(ix, iy, nz), id(ix+1, iy, nz), id(ix+1, iy+1, nz), id(ix, iy+1, nz), [0d0, 0d0, 1d0])
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

end subroutine make_box_hull

!> F between a hull and a cylinder C of the given radius about the origin (panel size 1), with nrows rows.
subroutine strips_for(mesh, radius, nrows, fm, ok, msg)
  type(panel_mesh), intent(in) :: mesh
  real(8), intent(in) :: radius
  integer, intent(in) :: nrows
  type(panel_mesh), intent(out) :: fm
  logical, intent(out) :: ok
  character(len=:), allocatable, intent(out) :: msg
  type(surface_geom) :: geom
  type(panel_mesh) :: cm
  type(cs_points) :: pts

  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  if (.not. ok) return
  call cylinder_surface([0d0, 0d0], radius, 8d0, 1d0, -1d0, cm)
  call add_surface_points(cm, 1, pts)
  call free_surface_strips(mesh, geom, 1, pts, nrows, fm, ok, msg)
end subroutine strips_for

!> Every F panel lies in z = 0, with the normal (0, 0, -1) and a positive area.
subroutine check_flat_down(fm, error)
  type(panel_mesh), intent(in) :: fm
  type(error_type), allocatable, intent(out) :: error
  call check(error, fm%n > 0, "F has panels")
  if (allocated(error)) return
  call check(error, maxval(abs(fm%xnode(3, :))) < 1d-12, "F in z = 0")
  if (allocated(error)) return
  call check(error, all(fm%nrm(3, :) < -0.999999d0) .and. all(fm%area > 0d0), "F normals down, positive areas")
end subroutine check_flat_down

subroutine test_strips_area(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, fm
  logical :: ok
  character(len=:), allocatable :: msg
  real(8) :: r
  integer :: i

  call make_cylinder(5d0, 3d0, 32, 2, 2, mesh, .true.)
  call strips_for(mesh, 15d0, 4, fm, ok, msg)
  call check(error, ok, "strips built: "//msg)
  if (allocated(error)) return
  call check_flat_down(fm, error)
  if (allocated(error)) return
  call check(error, abs(sum(fm%area)/(pi*(15d0**2 - 5d0**2)) - 1d0) < 5d-3, "area of F")
  if (allocated(error)) return
  do i = 1, fm%n
    r = hypot(fm%xc(1, i), fm%xc(2, i))
    call check(error, r > 5d0 .and. r < 15d0, "centroid between the waterline and Gamma_c")
    if (allocated(error)) return
  end do
end subroutine test_strips_area

subroutine test_strips_box_hull(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, fm
  logical :: ok
  character(len=:), allocatable :: msg

  call make_box_hull(20d0, 10d0, 5d0, 10, 5, 2, mesh)
  call strips_for(mesh, 25d0, 4, fm, ok, msg)
  call check(error, ok, "strips built: "//msg)
  if (allocated(error)) return
  call check_flat_down(fm, error)
  if (allocated(error)) return
  call check(error, abs(sum(fm%area)/(pi*25d0**2 - 200d0) - 1d0) < 5d-3, "area of F")
end subroutine test_strips_box_hull

!> An L-shaped waterline, (0,0) (30,0) (30,10) (10,10) (10,30) (0,30), is not star-shaped about the centre (15,15)
!> of its bounding box.
subroutine test_not_star_shaped_stops(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, cm, fm
  type(surface_geom) :: geom
  type(cs_points) :: pts
  logical :: ok
  character(len=:), allocatable :: msg
  real(8) :: v(2, 6)
  integer :: k

  call make_cylinder(5d0, 3d0, 16, 2, 2, mesh, .true.)
  v = reshape([0d0, 0d0, 30d0, 0d0, 30d0, 10d0, 10d0, 10d0, 10d0, 30d0, 0d0, 30d0], [2, 6])
  geom%nwl = 6
  allocate (geom%wl_x1(3, 6), geom%wl_x2(3, 6), geom%wl_xm(3, 6), geom%wl_len(6), geom%wl_body(6))
  do k = 1, 6
    geom%wl_x1(:, k) = [v(:, k), 0d0]
    geom%wl_x2(:, k) = [v(:, mod(k, 6) + 1), 0d0]
  end do
  geom%wl_xm = 0.5d0*(geom%wl_x1 + geom%wl_x2)
  geom%wl_len = norm2(geom%wl_x2 - geom%wl_x1, dim=1)
  geom%wl_body = 1
  call cylinder_surface([15d0, 15d0], 40d0, 8d0, 2d0, -1d0, cm)
  call add_surface_points(cm, 1, pts)
  call free_surface_strips(mesh, geom, 1, pts, 3, fm, ok, msg)
  call check(error, .not. ok, "L-shaped waterline stops")
  if (allocated(error)) return
  call check(error, msg == 'body 1: the waterline is not star-shaped about its centre; give FSMesh_1.pnl.', &
             "message: "//msg)
end subroutine test_not_star_shaped_stops

!> build_control_surface with with_free appends F as CS_FREE points; without it, the points are those of C only.
subroutine test_with_free_appends_f(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(so_options) :: opts
  type(cs_points) :: p0, p1
  type(panel_mesh), allocatable :: cms(:)
  logical :: ok
  character(len=:), allocatable :: msg
  integer :: n0

  call make_cylinder(5d0, 3d0, 32, 2, 2, mesh, .true.)
  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  opts%cs_margin = 6d0
  opts%cs_panel_size = 1d0
  call build_control_surface(mesh, geom, opts, -1d0, '.', p0, cms, ok, msg)
  call check(error, ok .and. count(p0%kind == CS_FREE) == 0, "no F without with_free")
  if (allocated(error)) return
  call build_control_surface(mesh, geom, opts, -1d0, '.', p1, cms, ok, msg, with_free=.true.)
  call check(error, ok, "build with F: "//msg)
  if (allocated(error)) return
  n0 = p0%n
  call check(error, p1%n > n0 .and. count(p1%kind == CS_FREE) == p1%n - n0, "F points appended after C")
  if (allocated(error)) return
  call check(error, all(p1%x(:, 1:n0) == p0%x) .and. all(p1%w(1:n0) == p0%w) .and. all(p1%kind(1:n0) == p0%kind), &
             "C points unchanged")
  if (allocated(error)) return
  call check(error, all(p1%nrm(1, n0+1:) == 0d0) .and. all(p1%nrm(2, n0+1:) == 0d0) .and. all(p1%nrm(3, n0+1:) == -1d0) &
             .and. all(p1%body(n0+1:) == 1) .and. all(abs(p1%x(3, n0+1:)) < 1d-12), "F points: normal, body, z")
  if (allocated(error)) return
  ! Gamma_c has the radius 5 + 6 = 11
  call check(error, abs(sum(p1%w(n0+1:))/(pi*(11d0**2 - 5d0**2)) - 1d0) < 5d-3, "area of F")
end subroutine test_with_free_appends_f

!> The default row count is nint((Gamma_c distance - waterline distance)/cs_panel_size) = nint((10.99 - 4.99)/1) = 6.
subroutine test_fs_rows_default(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(so_options) :: opts
  type(cs_points) :: pd, p3, p6
  type(panel_mesh), allocatable :: cms(:)
  logical :: ok
  character(len=:), allocatable :: msg

  call make_cylinder(5d0, 3d0, 32, 2, 2, mesh, .true.)
  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  opts%cs_margin = 6d0
  opts%cs_panel_size = 1d0
  call build_control_surface(mesh, geom, opts, -1d0, '.', pd, cms, ok, msg, with_free=.true.)
  opts%cs_fs_rows = 3
  call build_control_surface(mesh, geom, opts, -1d0, '.', p3, cms, ok, msg, with_free=.true.)
  opts%cs_fs_rows = 6
  call build_control_surface(mesh, geom, opts, -1d0, '.', p6, cms, ok, msg, with_free=.true.)
  call check(error, count(p3%kind == CS_FREE) > 0 .and. count(p6%kind == CS_FREE) == 2*count(p3%kind == CS_FREE), &
             "cs_fs_rows sets the rows")
  if (allocated(error)) return
  call check(error, count(pd%kind == CS_FREE) == count(p6%kind == CS_FREE), "default is 6 rows")
end subroutine test_fs_rows_default

!> FSMesh_1.pnl in the input directory replaces the strips: here a single 2 x 2 panel.
subroutine test_fs_file_takes_precedence(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(so_options) :: opts
  type(cs_points) :: pts
  type(panel_mesh), allocatable :: cms(:)
  logical :: ok
  character(len=:), allocatable :: msg
  integer :: u, ios

  open (newunit=u, file='./FSMesh_1.pnl', status='replace', action='write')
  write (u, '(a)') '    --------------Hull Mesh File---------------', ' ', &
                   '    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry'
  write (u, '(4i12)') 1, 4, 0, 0
  write (u, '(a)') ' ', '    # Start Definition of Node Coordinates'
  write (u, '(i5,3f20.12)') 1, 6d0, -1d0, 0d0
  write (u, '(i5,3f20.12)') 2, 6d0, 1d0, 0d0
  write (u, '(i5,3f20.12)') 3, 8d0, 1d0, 0d0
  write (u, '(i5,3f20.12)') 4, 8d0, -1d0, 0d0
  write (u, '(a)') '    # End Definition of Node Coordinates', ' ', '  # Start Definition of Node Relations'
  write (u, '(i5,i9,4i14)') 1, 4, 1, 2, 3, 4
  write (u, '(a)') '    # End Definition of Node Relations', ' ', '    --------------End Hull Mesh File---------------'
  close (u)

  call make_cylinder(5d0, 3d0, 32, 2, 2, mesh, .true.)
  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  opts%cs_margin = 6d0
  call build_control_surface(mesh, geom, opts, -1d0, '.', pts, cms, ok, msg, with_free=.true.)
  open (newunit=u, file='./FSMesh_1.pnl', status='old', iostat=ios)
  if (ios == 0) close (u, status='delete')
  call check(error, ok, "build with the file F: "//msg)
  if (allocated(error)) return
  call check(error, count(pts%kind == CS_FREE) == 1, "one F point from the file")
  if (allocated(error)) return
  call check(error, abs(sum(pts%w, mask=pts%kind == CS_FREE) - 4d0) < 1d-9, "its area")
end subroutine test_fs_file_takes_precedence

end module test_control_surface_mesh
