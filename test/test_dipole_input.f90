!> Input of dipole (zero-thickness) panels: the dipole block of HullMesh.pnl (NDIPOLE, the fifth field of line 4), read
!> by HULL_HEADER and ReadBodyMesh and set up by DIPOLE_SETUP (ReadPanelMesh), and the checks of the v1 limits. Fixtures:
!> test/data/dipole/disc_deep/input, a disc of radius 1 m at z = -50 m, 12 rings graded toward the edge by 32 sectors,
!> HAMS normal +z, as a HullMesh.pnl without hull panels; test/data/dipole/column_plate/nolid/input, a column with an
!> annulus plate (hull and dipole panels on one node list). Uses module arrays: serial.
module test_dipole_input

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use ReadPanelMesh, only : DIPOLE_SETUP, DIPOLE_CHECK_OPTIONS, DIPOLE_CHECK_GEOMETRY, DIPOLE_RELEASE, &
                            DIPOLE_CHECK_CLEARANCE, DIPOLE_CHECK_COPLANAR, HULL_HEADER, ReadBodyMesh, CalNormals
  use DipoleMesh_mod, only : DNELEM, DNTND, DPXYZ, DPXYZ_C, DNRM, DDS, DPNSZ, DNCN, DNCON
  use PanelMesh_mod, only : NELEM, NTND, ISX, ISY, XYZ, DS, PNSZ, NCN, NCON, XYZ_P, DXYZ_P, GEONRM_P
  use Body_mod, only : XR

  implicit none
  private
  public :: collect_dipole_input

  character(len=*), parameter :: disc_dir = "test/data/dipole/disc_deep/input"
  character(len=*), parameter :: plate_dir = "test/data/dipole/column_plate/nolid/input"
  character(len=*), parameter :: plain_dir = "test/data/image_term/image_cos/input"
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
      new_unittest("no_dipole_block_no_change", test_no_file), &
      new_unittest("dipole_rejects_near_free_surface", test_clearance), &
      new_unittest("dipole_rejects_coplanar_overlap", test_coplanar), &
      new_unittest("hull_header_fifth_field", test_header), &
      new_unittest("dipole_block_shared_nodes", test_shared_nodes), &
      new_unittest("dipole_block_rejects_bad_node", test_bad_node) &
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
    call read_hull(disc_dir, ok)
    call check(error, ok .and. NELEM == 0 .and. DNELEM == 384, "disc: header 0 hull panels, NDIPOLE = 384")
    if (allocated(error)) goto 9
    call DIPOLE_SETUP(1, 1, 1, .false., ok, msg)
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
    call release_hull
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
    ! through DIPOLE_SETUP with the fixture: stops before setting up the sheet, and frees the panels read
    call read_hull(disc_dir, ok)
    call check(error, ok .and. DNELEM == 384, "disc read")
    if (.not. allocated(error)) then
      call DIPOLE_SETUP(2, 1, 1, .false., ok, msg)
      call check(error, .not. ok .and. DNELEM == 0 .and. .not. (allocated(DPXYZ) .or. allocated(DNCN)), &
                 "setup with symmetry: not set up")
    end if
    call DIPOLE_RELEASE
    call release_hull
  end subroutine test_symmetry

  subroutine test_second_order(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: ok
    character(len=200) :: msg
    call DIPOLE_CHECK_OPTIONS(1, 1, 1, .true., ok, msg)
    call check(error, .not. ok .and. index(msg, "first order only") > 0, "second order rejected: "//trim(msg))
  end subroutine test_second_order

  !> A HullMesh.pnl with four fields on line 4: no dipole block read, nothing allocated, DNELEM = 0, and DIPOLE_SETUP
  !> accepts it whatever the options.
  subroutine test_no_file(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: ok
    character(len=200) :: msg
    call read_hull(plain_dir, ok)
    call check(error, ok .and. NELEM == 180 .and. NTND == 181 .and. DNELEM == 0, "four fields: 180 panels, DNELEM = 0")
    if (allocated(error)) goto 9
    call DIPOLE_SETUP(2, 2, 1, .true., ok, msg)
    call check(error, ok .and. DNELEM == 0 .and. DNTND == 0, "no dipole block: DNELEM = 0, OK")
    if (allocated(error)) goto 9
    call check(error, .not. (allocated(DPXYZ) .or. allocated(DPXYZ_C) .or. allocated(DNRM) .or. allocated(DDS) &
                             .or. allocated(DPNSZ) .or. allocated(DNCN) .or. allocated(DNCON)), &
               "no dipole block: no dipole array allocated")
9   call release_hull
  end subroutine test_no_file

  !> Line 4 of HullMesh.pnl: four fields give NDIPOLE = 0, a fifth is NDIPOLE; a negative NDIPOLE or fewer than four
  !> numbers are rejected.
  subroutine test_header(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: ne, nt, kx, ky, nd
    logical :: ok
    call HULL_HEADER("         192         744           0           0", ne, nt, kx, ky, nd, ok)
    call check(error, ok .and. ne == 192 .and. nt == 744 .and. kx == 0 .and. ky == 0 .and. nd == 0, "four fields")
    if (allocated(error)) return
    call HULL_HEADER("  192  1128  0  1  96", ne, nt, kx, ky, nd, ok)
    call check(error, ok .and. ne == 192 .and. nt == 1128 .and. kx == 0 .and. ky == 1 .and. nd == 96, "five fields")
    if (allocated(error)) return
    call HULL_HEADER("  192  744  0  0  0", ne, nt, kx, ky, nd, ok)
    call check(error, ok .and. nd == 0, "fifth field 0")
    if (allocated(error)) return
    call HULL_HEADER("  192  744  0  0  -1", ne, nt, kx, ky, nd, ok)
    call check(error, .not. ok, "negative NDIPOLE rejected")
    if (allocated(error)) return
    call HULL_HEADER("  192  744", ne, nt, kx, ky, nd, ok)
    call check(error, .not. ok, "two fields rejected")
  end subroutine test_header

  !> The column with an annulus plate on one node list: 192 hull panels on nodes 1 to 744, 96 dipole panels on nodes 745
  !> to 1128. DIPOLE_SETUP takes the 384 sheet nodes out of the list in order, so DPXYZ is XYZ(745:1128) and the
  !> renumbered panel k has nodes 4k-3 to 4k (as written). The plate lies at z = -1 outside the column: the checks pass.
  subroutine test_shared_nodes(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: ok
    character(len=200) :: msg
    integer :: k
    call read_hull(plate_dir, ok)
    call check(error, ok .and. NELEM == 192 .and. NTND == 1128 .and. DNELEM == 96, "header 192 1128 0 0 96")
    if (allocated(error)) goto 9
    call check(error, all(DNCON(:,1:4) > 744), "dipole panels refer to nodes after the hull nodes")
    if (allocated(error)) goto 9
    call DIPOLE_SETUP(1, 1, 1, .false., ok, msg)
    call check(error, ok, "column plate set up: "//trim(msg))
    if (allocated(error)) goto 9
    call check(error, DNTND == 384 .and. all(DPXYZ == XYZ(745:1128,:)), "sheet nodes 745 to 1128 copied in order")
    if (allocated(error)) goto 9
    do k = 1, DNELEM
      if (DNCN(k) /= 4 .or. any(DNCON(k,1:4) /= [4*k-3, 4*k-2, 4*k-1, 4*k])) then
        call check(error, .false., "renumbered connectivity")
        goto 9
      end if
    end do
    call check(error, all(abs(DPXYZ_C(:,3) + 1d0) < 1d-12) .and. all(DDS > 0d0), "sheet at z = -1")
9   call DIPOLE_RELEASE
    call release_hull
  end subroutine test_shared_nodes

  !> A dipole panel that refers to a node outside the node list stops the run.
  subroutine test_bad_node(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: ok
    character(len=200) :: msg
    call read_hull(plate_dir, ok)
    call check(error, ok .and. DNELEM == 96, "column plate read")
    if (allocated(error)) goto 9
    DNCON(5,2) = NTND + 1
    call DIPOLE_SETUP(1, 1, 1, .false., ok, msg)
    call check(error, .not. ok .and. index(msg, "dipole panel 5 ") > 0 .and. index(msg, "outside the node list") > 0 &
               .and. DNELEM == 0 .and. .not. allocated(DNCN), "node outside the list rejected: "//trim(msg))
9   call DIPOLE_RELEASE
    call release_hull
  end subroutine test_bad_node

  !> Reads DIR/HullMesh.pnl as HAMS_Prog does for one body without the lid: line 4 by HULL_HEADER (DNELEM = NDIPOLE),
  !> the nodes, hull panels and dipole block by ReadBodyMesh (unit 2), then the hull panel properties by CalNormals.
  subroutine read_hull(dir, ok)
    character(len=*), intent(in) :: dir
    logical, intent(out) :: ok
    character(len=512) :: line
    integer :: k, ios, nd
    call release_hull
    open(2, file=dir//"/HullMesh.pnl", status="old", action="read", iostat=ios)
    ok = (ios == 0)
    if (.not. ok) return
    do k = 1, 3
      read(2, *)
    end do
    read(2, '(a)') line
    call HULL_HEADER(line, NELEM, NTND, ISX, ISY, nd, ok)
    if (.not. ok) then
      close(2)
      return
    end if
    DNELEM = nd
    do k = 1, 2
      read(2, *)
    end do
    allocate(XYZ(NTND,3), DS(NELEM), PNSZ(NELEM), NCN(NELEM), NCON(NELEM,4))
    allocate(XYZ_P(NELEM,3), DXYZ_P(NELEM,6), GEONRM_P(NELEM,3))
    call ReadBodyMesh
    close(2)
    if (NELEM > 0) call CalNormals(0, dir)
  end subroutine read_hull

  !> Frees the hull arrays read_hull allocates.
  subroutine release_hull
    if (allocated(XYZ)) deallocate(XYZ)
    if (allocated(DS)) deallocate(DS)
    if (allocated(PNSZ)) deallocate(PNSZ)
    if (allocated(NCN)) deallocate(NCN)
    if (allocated(NCON)) deallocate(NCON)
    if (allocated(XYZ_P)) deallocate(XYZ_P)
    if (allocated(DXYZ_P)) deallocate(DXYZ_P)
    if (allocated(GEONRM_P)) deallocate(GEONRM_P)
  end subroutine release_hull

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
