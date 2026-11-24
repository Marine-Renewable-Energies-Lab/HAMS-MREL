!> Test geometry for dipole panels, built straight into the solver's module arrays: a column (vertical cylinder of radius
!> 1 m and draft T) with a horizontal plate of outer radius B at its bottom, in three forms
!>   SHEET   column wall, the column's bottom disc (hull), and the annulus 1 <= r <= B as a dipole sheet at z = -T;
!>   CLOSED  column wall, and the plate as a closed body of thickness TH: top annulus at z = -T, outer rim, and a bottom
!>           disc 0 <= r <= B at z = -T - TH (the shape the closed-body heave plate has);
!>   COLUMN  column wall and bottom disc only;
!> and optionally a lid (the waterplane disc r < 1 at z = 0) for IRSP = 1. Panels are quadrilaterals (triangles at the
!> centre) with NT sectors; the wall has NZ rows graded toward the bottom, discs and the annulus NR rings graded toward
!> their outer edge (radii a sin(pi i / 2 N), design note Sec. 6). Hull normals point out of the fluid, as CalTransNormals
!> gives them; each panel's node order is set from the wanted normal.
module dipole_bodies

  use PanelMesh_mod, only : NELEM, NTND, NSYS, ISX, ISY, XYZ, NCN, NCON, XYZ_P, DS, PNSZ, DXYZ_P, GEONRM_P
  use Inerfs_mod, only : iNELEM, iNTND, iXYZ, iNCN, iNCON, iXYZ_P, iDXYZ_P, iDS, iPNSZ
  use Potentials_mod, only : RKBN, CGRN, PKBN, DGRN, MXPOT
  use DipoleMesh_mod, only : DNELEM, DNTND, DPXYZ, DNCN, DNCON, DPXYZ_C, DNRM, DDS, DPNSZ
  use Body_mod, only : XR
  use ReadPanelMesh, only : DIPOLE_RELEASE
  use NormalProcess, only : CalPanelCentre, CalPanelArea, CalPanelChartLength, CalTransNormals, CalRotNormals

  implicit none
  private
  public :: build_body, release_body, SHEET, CLOSED, COLUMN

  integer, parameter :: SHEET = 1, CLOSED = 2, COLUMN = 3
  real(8), parameter :: pi = 3.141592653589793d0

  ! growing panel lists: hull (h), sheet (d), lid (l); four corners each (a triangle repeats its last corner)
  real(8), allocatable :: hp(:,:,:), dp(:,:,:), lp(:,:,:)
  integer, allocatable :: hn(:), dn(:), ln(:)
  integer :: nh, nd, nl

contains

  !> kind: SHEET, CLOSED or COLUMN. t draft, b plate radius, th thickness (CLOSED), nt sectors, nz wall rows, nr rings
  !> of each disc and of the annulus, lid: add the waterplane lid (nr rings).
  subroutine build_body(kind, t, b, th, nt, nz, nr, lid)
    integer, intent(in) :: kind, nt, nz, nr
    real(8), intent(in) :: t, b, th
    logical, intent(in) :: lid
    integer :: i
    real(8) :: zb
    call release_body
    allocate(hp(4,3,4000), dp(4,3,4000), lp(4,3,4000), hn(4000), dn(4000), ln(4000))
    nh = 0
    nd = 0
    nl = 0
    ! column wall, z from 0 to -t, rows graded toward the bottom; normal into the column (-radial)
    do i = 1, nz
      call ring_band(1d0, -t*sin(0.5d0*pi*(i - 1)/nz), -t*sin(0.5d0*pi*i/nz), nt, 'h')
    end do
    select case (kind)
    case (SHEET)
      call disc_rings(0d0, 1d0, -t, nt, nr, +1d0, 'h')            ! column bottom, normal +z (into the column)
      call disc_rings(1d0, b, -t, nt, nr, +1d0, 'd')              ! the sheet
    case (COLUMN)
      call disc_rings(0d0, 1d0, -t, nt, nr, +1d0, 'h')
    case (CLOSED)
      zb = -t - th
      call disc_rings(1d0, b, -t, nt, nr, -1d0, 'h')              ! plate top, normal -z (into the plate)
      call ring_band(b, -t, zb, nt, 'h')                          ! outer rim, normal -radial
      call disc_rings(0d0, 1d0, zb, nt, nr, +1d0, 'h')            ! plate bottom, normal +z
      call disc_rings(1d0, b, zb, nt, nr, +1d0, 'h')
    end select
    if (lid) call disc_rings(0d0, 1d0, 0d0, nt, nr, +1d0, 'l')
    call to_modules
  end subroutine build_body

  !> A vertical band at radius r between z1 > z2, normal toward the axis.
  subroutine ring_band(r, z1, z2, nt, which)
    real(8), intent(in) :: r, z1, z2
    integer, intent(in) :: nt
    character, intent(in) :: which
    integer :: j
    real(8) :: t0, t1, c(4,3), want(3)
    do j = 0, nt - 1
      t0 = 2d0*pi*j/nt
      t1 = 2d0*pi*(j + 1)/nt
      c(1,:) = [r*cos(t0), r*sin(t0), z1]
      c(2,:) = [r*cos(t1), r*sin(t1), z1]
      c(3,:) = [r*cos(t1), r*sin(t1), z2]
      c(4,:) = [r*cos(t0), r*sin(t0), z2]
      want = -[cos(0.5d0*(t0 + t1)), sin(0.5d0*(t0 + t1)), 0d0]
      call add_panel(c, 4, want, which)
    end do
  end subroutine ring_band

  !> Rings r0 <= r <= r1 at height z with radii graded toward r1 (r0 + (r1 - r0) sin(pi i / 2 nr)); a centre disc (r0 = 0)
  !> starts with triangles. Normal sgn * z.
  subroutine disc_rings(r0, r1, z, nt, nr, sgn, which)
    real(8), intent(in) :: r0, r1, z, sgn
    integer, intent(in) :: nt, nr
    character, intent(in) :: which
    integer :: i, j
    real(8) :: ra, rb, t0, t1, c(4,3)
    do i = 1, nr
      ra = r0 + (r1 - r0)*sin(0.5d0*pi*(i - 1)/nr)
      rb = r0 + (r1 - r0)*sin(0.5d0*pi*i/nr)
      do j = 0, nt - 1
        t0 = 2d0*pi*j/nt
        t1 = 2d0*pi*(j + 1)/nt
        if (ra < 1d-12) then
          c(1,:) = [0d0, 0d0, z]
          c(2,:) = [rb*cos(t0), rb*sin(t0), z]
          c(3,:) = [rb*cos(t1), rb*sin(t1), z]
          c(4,:) = c(3,:)
          call add_panel(c, 3, [0d0, 0d0, sgn], which)
        else
          c(1,:) = [ra*cos(t0), ra*sin(t0), z]
          c(2,:) = [rb*cos(t0), rb*sin(t0), z]
          c(3,:) = [rb*cos(t1), rb*sin(t1), z]
          c(4,:) = [ra*cos(t1), ra*sin(t1), z]
          call add_panel(c, 4, [0d0, 0d0, sgn], which)
        end if
      end do
    end do
  end subroutine disc_rings

  !> Store a panel with its nodes ordered so that CalTransNormals gives the normal along WANT: that normal is opposite
  !> to the right-hand normal of the node order (design note item 1.1).
  subroutine add_panel(c, n, want, which)
    real(8), intent(in) :: c(4,3), want(3)
    integer, intent(in) :: n
    character, intent(in) :: which
    real(8) :: m(3), u(3), w(3), cc(4,3)
    if (n == 4) then
      u = c(3,:) - c(1,:)
      w = c(4,:) - c(2,:)
    else
      u = c(2,:) - c(1,:)
      w = c(3,:) - c(1,:)
    end if
    m = [u(2)*w(3) - u(3)*w(2), u(3)*w(1) - u(1)*w(3), u(1)*w(2) - u(2)*w(1)]
    cc = c
    if (dot_product(m, want) > 0d0) cc(1:n,:) = c(n:1:-1,:)
    select case (which)
    case ('h')
      nh = nh + 1
      hp(:,:,nh) = cc
      hn(nh) = n
    case ('d')
      nd = nd + 1
      dp(:,:,nd) = cc
      dn(nd) = n
    case ('l')
      nl = nl + 1
      lp(:,:,nl) = cc
      ln(nl) = n
    end select
  end subroutine add_panel

  !> Panels to the solver's arrays (each panel with its own nodes), then the geometry routines of CalNormals.
  subroutine to_modules
    integer :: k, j
    NELEM = nh
    NTND = 4*nh
    NSYS = 1
    ISX = 0
    ISY = 0
    XR = 0d0
    allocate(XYZ(NTND,3), NCN(NELEM), NCON(NELEM,4), XYZ_P(NELEM,3), DS(NELEM), PNSZ(NELEM), DXYZ_P(NELEM,6))
    allocate(GEONRM_P(NELEM,3))
    NCON = 0
    do k = 1, nh
      NCN(k) = hn(k)
      do j = 1, hn(k)
        XYZ(4*(k - 1) + j,:) = hp(j,:,k)
        NCON(k,j) = 4*(k - 1) + j
      end do
    end do
    call CalPanelCentre(XYZ, NTND, NELEM, NCN, NCON, XYZ_P)
    call CalPanelArea(XYZ, NTND, NELEM, NCN, NCON, DS)
    call CalPanelChartLength(XYZ, NTND, NELEM, NCN, NCON, PNSZ)
    call CalTransNormals(XYZ, NTND, NELEM, NCN, NCON, DXYZ_P)
    call CalRotNormals(XR, XYZ_P, DXYZ_P, NELEM)
    GEONRM_P = DXYZ_P(:,1:3)
    allocate(RKBN(NELEM,NELEM,1,4), CGRN(NELEM,NELEM,1,4))

    DNELEM = nd
    DNTND = 4*nd
    if (nd > 0) then
      allocate(DPXYZ(DNTND,3), DNCN(DNELEM), DNCON(DNELEM,4), DPXYZ_C(DNELEM,3), DNRM(DNELEM,6), DDS(DNELEM))
      allocate(DPNSZ(DNELEM))
      DNCON = 0
      do k = 1, nd
        DNCN(k) = dn(k)
        do j = 1, dn(k)
          DPXYZ(4*(k - 1) + j,:) = dp(j,:,k)
          DNCON(k,j) = 4*(k - 1) + j
        end do
      end do
      call CalPanelCentre(DPXYZ, DNTND, DNELEM, DNCN, DNCON, DPXYZ_C)
      call CalPanelArea(DPXYZ, DNTND, DNELEM, DNCN, DNCON, DDS)
      call CalPanelChartLength(DPXYZ, DNTND, DNELEM, DNCN, DNCON, DPNSZ)
      call CalTransNormals(DPXYZ, DNTND, DNELEM, DNCN, DNCON, DNRM)
      call CalRotNormals(XR, DPXYZ_C, DNRM, DNELEM)
    end if

    iNELEM = nl
    iNTND = 4*nl
    if (nl > 0) then
      allocate(iXYZ(iNTND,3), iNCN(iNELEM), iNCON(iNELEM,4), iXYZ_P(iNELEM,3), iDXYZ_P(iNELEM,6), iDS(iNELEM))
      allocate(iPNSZ(iNELEM))
      iNCON = 0
      do k = 1, nl
        iNCN(k) = ln(k)
        do j = 1, ln(k)
          iXYZ(4*(k - 1) + j,:) = lp(j,:,k)
          iNCON(k,j) = 4*(k - 1) + j
        end do
      end do
      call CalPanelCentre(iXYZ, iNTND, iNELEM, iNCN, iNCON, iXYZ_P)
      call CalPanelArea(iXYZ, iNTND, iNELEM, iNCN, iNCON, iDS)
      call CalPanelChartLength(iXYZ, iNTND, iNELEM, iNCN, iNCON, iPNSZ)
      call CalTransNormals(iXYZ, iNTND, iNELEM, iNCN, iNCON, iDXYZ_P)
      call CalRotNormals(XR, iXYZ_P, iDXYZ_P, iNELEM)
      allocate(PKBN(iNELEM,NELEM,1,4), DGRN(iNELEM,NELEM,1,4))
    end if
    deallocate(hp, dp, lp, hn, dn, ln)
  end subroutine to_modules

  !> Free everything build_body allocated (and what other suites may have left in these arrays).
  subroutine release_body
    if (allocated(XYZ)) deallocate(XYZ)
    if (allocated(NCN)) deallocate(NCN)
    if (allocated(NCON)) deallocate(NCON)
    if (allocated(XYZ_P)) deallocate(XYZ_P)
    if (allocated(DS)) deallocate(DS)
    if (allocated(PNSZ)) deallocate(PNSZ)
    if (allocated(DXYZ_P)) deallocate(DXYZ_P)
    if (allocated(GEONRM_P)) deallocate(GEONRM_P)
    if (allocated(RKBN)) deallocate(RKBN)
    if (allocated(CGRN)) deallocate(CGRN)
    if (allocated(PKBN)) deallocate(PKBN)
    if (allocated(DGRN)) deallocate(DGRN)
    if (allocated(iXYZ)) deallocate(iXYZ)
    if (allocated(iNCN)) deallocate(iNCN)
    if (allocated(iNCON)) deallocate(iNCON)
    if (allocated(iXYZ_P)) deallocate(iXYZ_P)
    if (allocated(iDXYZ_P)) deallocate(iDXYZ_P)
    if (allocated(iDS)) deallocate(iDS)
    if (allocated(iPNSZ)) deallocate(iPNSZ)
    if (allocated(MXPOT)) deallocate(MXPOT)
    if (allocated(hp)) deallocate(hp, dp, lp, hn, dn, ln)
    iNELEM = 0
    NELEM = 0
    call DIPOLE_RELEASE
  end subroutine release_body

end module dipole_bodies
