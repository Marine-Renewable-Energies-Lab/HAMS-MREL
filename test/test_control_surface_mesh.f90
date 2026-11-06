!> Unit tests for the generated control surfaces C, the line Gamma_c and the checks against the other bodies.
module test_control_surface_mesh
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderTypes, only : panel_mesh, surface_geom, so_options, cs_points, CS_WALL, CS_LINE
  use SurfaceGeometry, only : build_surface_geom, DEFAULT_SHARP_DEG
  use ControlSurfaceMesh, only : cylinder_surface, box_surface, add_surface_points, build_control_surface, &
                                 check_control_surface, write_control_surface
  use second_order_fixtures, only : make_cylinder
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
    new_unittest("file_shape_not_yet_available", test_file_shape_not_yet_available), &
    new_unittest("written_mesh_has_hull_format", test_written_mesh_has_hull_format) &
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

subroutine test_file_shape_not_yet_available(error)
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
  call check(error, .not. okb .and. len(msg) > 0, "file shape reports that it is not available")
end subroutine test_file_shape_not_yet_available

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

end module test_control_surface_mesh
