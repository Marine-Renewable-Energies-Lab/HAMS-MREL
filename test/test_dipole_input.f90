!> Input of dipole (zero-thickness) panels: DipoleMesh.pnl, read by DIPOLE_SETUP (ReadPanelMesh), and the checks of the
!> v1 limits (plan 2026-11-24-dipole-panels-first-order, Task 3). Fixture: test/data/dipole/disc_deep/input, a disc of
!> radius 1 m at z = -50 m, 12 rings graded toward the edge by 32 sectors, HAMS normal +z. Uses module arrays: serial.
module test_dipole_input

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use ReadPanelMesh, only : DIPOLE_SETUP, DIPOLE_CHECK_OPTIONS, DIPOLE_CHECK_GEOMETRY, DIPOLE_RELEASE, &
                            DIPOLE_CHECK_CLEARANCE, DIPOLE_CHECK_COPLANAR
  use DipoleMesh_mod, only : DNELEM, DNTND, DPXYZ, DPXYZ_C, DNRM, DDS, DPNSZ, DNCN
  use Body_mod, only : XR

  implicit none
  private
  public :: collect_dipole_input

  character(len=*), parameter :: disc_dir = "test/data/dipole/disc_deep/input"
  real(8), parameter :: pi = 3.141592653589793d0

contains

  subroutine collect_dipole_input(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("dipole_mesh_read", test_read), &
      new_unittest("dipole_mesh_rejects_surface_piercing", test_piercing), &
      new_unittest("dipole_mesh_rejects_coincident", test_coincident), &
      new_unittest("dipole_rejects_symmetry", test_symmetry), &
      new_unittest("dipole_rejects_second_order", test_second_order), &
      new_unittest("no_dipole_file_no_change", test_no_file), &
      new_unittest("dipole_rejects_near_free_surface", test_clearance), &
      new_unittest("dipole_rejects_coplanar_overlap", test_coplanar) &
    ]
  end subroutine collect_dipole_input

  !> Counts, total area (the 32-gon of radius 1, 16 sin(2 pi/32)), normals +z, heave normals 1 and roll/pitch normals
  !> from the centroids about XR, sizes positive.
  subroutine test_read(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: ok
    character(len=200) :: msg
    real(8) :: xr_save(3), area
    integer :: i
    xr_save = XR
    XR = [0d0, 0d0, -10d0]
    call DIPOLE_SETUP(disc_dir, 1, 1, 1, .false., ok, msg)
    call check(error, ok, "disc read: "//trim(msg))
    if (allocated(error)) goto 9
    call check(error, DNELEM == 384 .and. DNTND == 385, "disc: 384 panels, 385 nodes")
    if (allocated(error)) goto 9
    call check(error, count(DNCN == 3) == 32 .and. count(DNCN == 4) == 352, "disc: 32 triangles, 352 quadrilaterals")
    if (allocated(error)) goto 9
    area = sum(DDS)
    call check(error, abs(area - 16d0*sin(2d0*pi/32d0)) < 1d-10, "disc: total area of the 32-gon")
    if (allocated(error)) goto 9
    call check(error, all(abs(DNRM(:,3) - 1d0) < 1d-12) .and. all(abs(DNRM(:,1)) < 1d-12), "disc: normals +z")
    if (allocated(error)) goto 9
    do i = 1, DNELEM
      if (abs(DNRM(i,4) - (DPXYZ_C(i,2) - XR(2))) > 1d-12 .or. abs(DNRM(i,5) + (DPXYZ_C(i,1) - XR(1))) > 1d-12) then
        call check(error, .false., "disc: rotational normals (x - XR) x n")
        goto 9
      end if
    end do
    call check(error, all(DPNSZ > 0d0) .and. all(abs(DPXYZ_C(:,3) + 50d0) < 1d-12), "disc: sizes and centroids")
9   call DIPOLE_RELEASE
    XR = xr_save
  end subroutine test_read

  subroutine test_piercing(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: pxyz(3,3), pc(1,3), psz(1), hc(1,3)
    logical :: ok
    character(len=200) :: msg
    pxyz = reshape([0d0, 1d0, 0d0, 0d0, 0d0, 1d0, -1d0, -1d0, 0d0], [3,3])
    pc(1,:) = sum(pxyz, dim=1)/3d0
    psz = 1d0
    hc(1,:) = [5d0, 5d0, -5d0]
    call DIPOLE_CHECK_GEOMETRY(pxyz, 3, pc, psz, 1, hc, 1, ok, msg)
    call check(error, .not. ok .and. index(msg, "not below z = 0") > 0, "a node on z = 0 is rejected: "//trim(msg))
    if (allocated(error)) return
    pxyz(3,3) = -1d0
    call DIPOLE_CHECK_GEOMETRY(pxyz, 3, pc, psz, 1, hc, 1, ok, msg)
    call check(error, ok, "fully submerged and apart: accepted")
  end subroutine test_piercing

  subroutine test_coincident(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: pxyz(3,3), pc(1,3), psz(1), hc(2,3)
    logical :: ok
    character(len=200) :: msg
    pxyz = reshape([0d0, 1d0, 0d0, 0d0, 0d0, 1d0, -1d0, -1d0, -1d0], [3,3])
    pc(1,:) = sum(pxyz, dim=1)/3d0
    psz = 1d0
    hc(1,:) = [5d0, 5d0, -5d0]
    hc(2,:) = pc(1,:) + 1d-8
    call DIPOLE_CHECK_GEOMETRY(pxyz, 3, pc, psz, 1, hc, 2, ok, msg)
    call check(error, .not. ok .and. index(msg, "coincides with hull panel 2") > 0, "coincident centroid: "//trim(msg))
  end subroutine test_coincident

  subroutine test_symmetry(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: ok
    character(len=200) :: msg
    call DIPOLE_CHECK_OPTIONS(2, 1, 1, .false., ok, msg)
    call check(error, .not. ok .and. index(msg, "symmetry") > 0, "symmetry rejected: "//trim(msg))
    if (allocated(error)) return
    call DIPOLE_CHECK_OPTIONS(1, 2, 1, .false., ok, msg)
    call check(error, .not. ok .and. index(msg, "ISOLV") > 0, "GMRES rejected: "//trim(msg))
    if (allocated(error)) return
    call DIPOLE_CHECK_OPTIONS(1, 1, 2, .false., ok, msg)
    call check(error, .not. ok .and. index(msg, "one body") > 0, "multi-body rejected: "//trim(msg))
    if (allocated(error)) return
    call DIPOLE_CHECK_OPTIONS(1, 1, 1, .false., ok, msg)
    call check(error, ok, "v1 options accepted")
    if (allocated(error)) return
    ! through DIPOLE_SETUP with the fixture: stops before reading
    call DIPOLE_SETUP(disc_dir, 2, 1, 1, .false., ok, msg)
    call check(error, .not. ok .and. DNELEM == 0 .and. .not. allocated(DPXYZ), "setup with symmetry: not read")
  end subroutine test_symmetry

  subroutine test_second_order(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: ok
    character(len=200) :: msg
    call DIPOLE_CHECK_OPTIONS(1, 1, 1, .true., ok, msg)
    call check(error, .not. ok .and. index(msg, "first order only") > 0, "second order rejected: "//trim(msg))
  end subroutine test_second_order

  !> A directory with no DipoleMesh.pnl: nothing read, nothing allocated, DNELEM = 0.
  subroutine test_no_file(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: ok
    character(len=200) :: msg
    call DIPOLE_SETUP("test/data/config1/input", 2, 2, 1, .true., ok, msg)
    call check(error, ok .and. DNELEM == 0 .and. DNTND == 0, "no file: DNELEM = 0, OK")
    if (allocated(error)) return
    call check(error, .not. (allocated(DPXYZ) .or. allocated(DPXYZ_C) .or. allocated(DNRM) .or. allocated(DDS) &
                             .or. allocated(DPNSZ) .or. allocated(DNCN)), "no file: no dipole array allocated")
  end subroutine test_no_file

  !> A unit square panel at depth zc, nodes clockwise seen from +z (HAMS normal +z).
  subroutine square(x0, y0, zc, b, pxyz, ncon)
    real(8), intent(in) :: x0, y0, zc, b
    real(8), intent(out) :: pxyz(4,3)
    integer, intent(out) :: ncon(1,4)
    pxyz(:,1) = [x0, x0, x0 + b, x0 + b]
    pxyz(:,2) = [y0, y0 + b, y0 + b, y0]
    pxyz(:,3) = zc
    ncon(1,:) = [1, 2, 3, 4]
  end subroutine square

  !> A node less than one panel size below z = 0 stops the run (design note Sec. 3.2, 6).
  subroutine test_clearance(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: pxyz(4,3), psz(1)
    integer :: ncon(1,4), ncn(1)
    logical :: ok
    character(len=200) :: msg
    ncn = 4
    call square(0d0, 0d0, -0.9d0, 1d0, pxyz, ncon)
    psz = sqrt(2d0)
    call DIPOLE_CHECK_CLEARANCE(pxyz, ncn, ncon, psz, 1, ok, msg)
    call check(error, .not. ok .and. index(msg, "panel size") > 0, "0.9 m deep, size 1.41 m: rejected: "//trim(msg))
    if (allocated(error)) return
    call square(0d0, 0d0, -1.5d0, 1d0, pxyz, ncon)
    call DIPOLE_CHECK_CLEARANCE(pxyz, ncn, ncon, psz, 1, ok, msg)
    call check(error, ok, "1.5 m deep, size 1.41 m: accepted")
  end subroutine test_clearance

  !> A sheet panel laid over a coplanar hull panel (overlapping) stops the run; side by side, or in a parallel plane
  !> 0.1 m away, it does not. The hull panel normal is -z, the sheet's +z: the check ignores the normal direction.
  subroutine test_coplanar(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: dxyz(4,3), hxyz(4,3), dc(1,3), hc(1,3), dn(1,3), hn(1,3), dsz(1), hsz(1)
    integer :: dcon(1,4), hcon(1,4), dncn(1), hncn(1)
    logical :: ok
    character(len=200) :: msg
    dncn = 4
    hncn = 4
    dn(1,:) = [0d0, 0d0, 1d0]
    hn(1,:) = [0d0, 0d0, -1d0]
    dsz = sqrt(2d0)
    hsz = sqrt(2d0)
    call square(0d0, 0d0, -5d0, 1d0, hxyz, hcon)
    hc(1,:) = [0.5d0, 0.5d0, -5d0]
    ! overlapping: the sheet centroid (0.75, 0.5) is inside the hull panel
    call square(0.25d0, 0d0, -5d0, 1d0, dxyz, dcon)
    dc(1,:) = [0.75d0, 0.5d0, -5d0]
    call DIPOLE_CHECK_COPLANAR(dxyz, dncn, dcon, dc, dn, dsz, 1, hxyz, hncn, hcon, hc, hn, hsz, 1, ok, msg)
    call check(error, .not. ok .and. index(msg, "overlaps hull panel 1") > 0, "overlap rejected: "//trim(msg))
    if (allocated(error)) return
    ! side by side
    call square(1d0, 0d0, -5d0, 1d0, dxyz, dcon)
    dc(1,:) = [1.5d0, 0.5d0, -5d0]
    call DIPOLE_CHECK_COPLANAR(dxyz, dncn, dcon, dc, dn, dsz, 1, hxyz, hncn, hcon, hc, hn, hsz, 1, ok, msg)
    call check(error, ok, "side by side accepted")
    if (allocated(error)) return
    ! over it but 0.1 m below: not coplanar
    call square(0d0, 0d0, -5.1d0, 1d0, dxyz, dcon)
    dc(1,:) = [0.5d0, 0.5d0, -5.1d0]
    call DIPOLE_CHECK_COPLANAR(dxyz, dncn, dcon, dc, dn, dsz, 1, hxyz, hncn, hcon, hc, hn, hsz, 1, ok, msg)
    call check(error, ok, "parallel plane accepted")
  end subroutine test_coplanar

end module test_dipole_input
