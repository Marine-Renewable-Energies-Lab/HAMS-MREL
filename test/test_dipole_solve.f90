!> Solve with dipole (zero-thickness) panels (DipoleAssb; notes/specs/2026-11-24-dipole-panels-design.md, Sec. 3-4).
!> A disc of radius a, broadside, deep under the free surface, with no hull: the added mass in an infinite fluid is
!> (8/3) rho a^3 (Lamb 1932, Sec. 108). HAMS has no infinite-fluid G; its two frequency limits (V = 0: G = 1/r + 1/r1,
!> V = -1: G = 1/r - 1/r1) bracket that value and close in on it with depth (design note Sec. 4). The suite sets the
!> mesh arrays of the solver, which are module globals, so it runs serially.
module test_dipole_solve

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use PanelMesh_mod, only : NELEM, XYZ, NCN, NCON, XYZ_P, DS, PNSZ, DXYZ_P, GEONRM_P
  use Potentials_mod, only : MXPOT
  use DipoleMesh_mod, only : DNELEM, DNTND, DPXYZ, DNCN, DNCON, DPXYZ_C, DNRM, DDS, DPNSZ
  use WaveDyn_mod, only : V, H, WVN, WK, W1, AMP, BETA
  use Body_mod, only : XR, XW
  use HAMS_mod, only : NTHREAD
  use DipoleAssb, only : DIPOLE_ASSB_LEFT, DIPOLE_ASSB_RBC, DIPOLE_ASSB_DBC, DIPOLE_SOLVE, DIPOLE_FORCE_SUM, &
                         DIPOLE_RELEASE_ASSB
  use ReadPanelMesh, only : DIPOLE_RELEASE
  use NormalProcess, only : CalPanelCentre, CalPanelArea, CalPanelChartLength, CalTransNormals, CalRotNormals

  implicit none
  private
  public :: collect_dipole_solve

  real(8), parameter :: pi = 3.141592653589793d0

contains

  subroutine collect_dipole_solve(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("disc_infinite_fluid_added_mass", test_disc), &
      new_unittest("disc_normal_flip", test_flip) &
    ]
  end subroutine collect_dipole_solve

  !> A disc of radius a at depth d, NR rings with radii a sin(pi i / 2 NR) (graded toward the rim, design note Sec. 6)
  !> and NT sectors, into DipoleMesh_mod; FLIP reverses the node order (HAMS normal -z instead of +z). No hull panels.
  subroutine make_disc(a, d, nr, nt, flip)
    real(8), intent(in) :: a, d
    integer, intent(in) :: nr, nt
    logical, intent(in) :: flip
    integer :: i, j, jn, k, e
    real(8) :: r, t
    call DIPOLE_RELEASE
    call clear_hull
    DNTND = 1 + nr*nt
    DNELEM = nt*nr
    allocate(DPXYZ(DNTND,3), DNCN(DNELEM), DNCON(DNELEM,4))
    allocate(DPXYZ_C(DNELEM,3), DNRM(DNELEM,6), DDS(DNELEM), DPNSZ(DNELEM))
    DPXYZ(1,:) = [0d0, 0d0, -d]
    do i = 1, nr
      r = a*sin(0.5d0*pi*i/nr)
      do j = 0, nt - 1
        t = 2d0*pi*j/nt
        DPXYZ(1 + (i - 1)*nt + j + 1,:) = [r*cos(t), r*sin(t), -d]
      end do
    end do
    DNCON = 0
    e = 0
    do j = 0, nt - 1
      jn = mod(j + 1, nt)
      e = e + 1
      DNCN(e) = 3
      DNCON(e,1:3) = [1, 2 + jn, 2 + j]                       ! clockwise seen from +z: HAMS normal +z
    end do
    do i = 1, nr - 1
      do j = 0, nt - 1
        jn = mod(j + 1, nt)
        e = e + 1
        DNCN(e) = 4
        DNCON(e,1:4) = [2 + (i - 1)*nt + j, 2 + (i - 1)*nt + jn, 2 + i*nt + jn, 2 + i*nt + j]
      end do
    end do
    if (flip) then
      do k = 1, DNELEM
        DNCON(k,1:DNCN(k)) = DNCON(k,DNCN(k):1:-1)
      end do
    end if
    call CalPanelCentre(DPXYZ, DNTND, DNELEM, DNCN, DNCON, DPXYZ_C)
    call CalPanelArea(DPXYZ, DNTND, DNELEM, DNCN, DNCON, DDS)
    call CalPanelChartLength(DPXYZ, DNTND, DNELEM, DNCN, DNCON, DPNSZ)
    call CalTransNormals(DPXYZ, DNTND, DNELEM, DNCN, DNCON, DNRM)
    XR = [0d0, 0d0, -d]
    call CalRotNormals(XR, DPXYZ_C, DNRM, DNELEM)
  end subroutine make_disc

  !> No hull panels: zero-size hull arrays.
  subroutine clear_hull
    NELEM = 0
    if (allocated(XYZ)) deallocate(XYZ)
    if (allocated(NCN)) deallocate(NCN)
    if (allocated(NCON)) deallocate(NCON)
    if (allocated(XYZ_P)) deallocate(XYZ_P)
    if (allocated(DS)) deallocate(DS)
    if (allocated(PNSZ)) deallocate(PNSZ)
    if (allocated(DXYZ_P)) deallocate(DXYZ_P)
    if (allocated(GEONRM_P)) deallocate(GEONRM_P)
    allocate(XYZ(0,3), NCN(0), NCON(0,4), XYZ_P(0,3), DS(0), PNSZ(0), DXYZ_P(0,6), GEONRM_P(0,3))
  end subroutine clear_hull

  !> Radiation at the limit vlim (0 or -1) in deep water: A33/(rho a^3) = Re F3, F the sheet's force sum (Sec. 4).
  function heave_added_mass(vlim) result(a33)
    real(8), intent(in) :: vlim
    real(8) :: a33
    complex(8), allocatable :: amat(:,:), brmat(:,:)
    integer, allocatable :: ipiv(:)
    complex(8) :: f(6)
    integer :: ntot
    V = vlim
    H = -1d0
    WVN = 0d0
    ntot = NELEM + DNELEM
    if (allocated(MXPOT)) deallocate(MXPOT)
    allocate(MXPOT(ntot,7,1), amat(ntot,ntot), brmat(ntot,6), ipiv(ntot))
    call DIPOLE_ASSB_LEFT(amat, ipiv, ntot)
    call DIPOLE_ASSB_RBC(brmat, ntot)
    call DIPOLE_SOLVE(amat, ipiv, brmat, ntot, 6, 1)
    call DIPOLE_FORCE_SUM(3, f)
    a33 = real(f(3))
  end function heave_added_mass

  !> At d = 10 a and 20 a, both limits, centroid collocation: constant dipole panels converge as O(h) here (the 2D
  !> analogue has error exactly +1/N; desk check 2026-11-24: +8.2 % at 384 panels, +4.4 % at 1536, +2.7 % at 3072).
  !> So the gate is the Richardson extrapolation 2 A(1536) - A(384) (12 x 32 and 24 x 64 panels, ring radii graded to
  !> the rim) within 1 % of 8/3, as agreed with the hub. Also: the two limits differ by less than 0.1 % (the image
  !> effect at these depths), and less at 20 a than at 10 a.
  subroutine test_disc(error)
    real(8), parameter :: exact = 8d0/3d0, depths(2) = [10d0, 20d0]
    type(error_type), allocatable, intent(out) :: error
    real(8) :: a0(2), ai(2), ext, spread(2)
    integer :: k, m, nthread_saved
    integer, parameter :: nrs(2) = [12, 24], nts(2) = [32, 64]
    character(len=200) :: msg
    nthread_saved = NTHREAD
    NTHREAD = 4
    do k = 1, 2
      do m = 1, 2
        call make_disc(1d0, depths(k), nrs(m), nts(m), .false.)
        a0(m) = heave_added_mass(0d0)
        ai(m) = heave_added_mass(-1d0)
      end do
      spread(k) = abs(a0(2) - ai(2))/exact
      do m = 1, 2
        ext = merge(2d0*a0(2) - a0(1), 2d0*ai(2) - ai(1), m == 1)
        write(msg, '(a,f5.1,a,i2,a,2f10.6,a,f10.6,a,es10.2)') 'disc d/a =', depths(k), ', limit', m, &
             ': A33/(rho a^3) at 384 and 1536 panels', merge(a0, ai, m == 1), ', Richardson', ext, &
             ', rel err', (ext - exact)/exact
        print '(2x,a)', trim(msg)
        call check(error, abs(ext - exact) <= 1d-2*exact, "Richardson within 1 %: "//trim(msg))
        if (allocated(error)) goto 9
      end do
      call check(error, spread(k) < 1d-3, "V = 0 and V = -1 within 0.1 % at this depth")
      if (allocated(error)) goto 9
    end do
    print '(2x,a,2es10.2)', 'disc: |A(V=0) - A(V=-1)|/(8/3) at d/a = 10, 20:', spread
    call check(error, spread(2) < spread(1), "the image effect falls with depth")
9   NTHREAD = nthread_saved
    call cleanup
  end subroutine test_disc

  !> Review Focus 3: the same disc with the node order reversed (normal -z): the jumps change sign, A33 does not;
  !> at a finite frequency also B33 and the diffraction heave force (sum of -mu n_D dS) do not.
  subroutine test_flip(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: f1(6), f2(6), x1(6), x2(6)
    real(8) :: mu1, mu2
    integer :: nthread_saved
    nthread_saved = NTHREAD
    NTHREAD = 4
    call make_disc(1d0, 3d0, 6, 16, .false.)
    call solve_finite(f1, x1, mu1)
    call make_disc(1d0, 3d0, 6, 16, .true.)
    call solve_finite(f2, x2, mu2)
    print '(2x,a,2es14.6,a,2es14.6)', 'flip: F3', f1(3), ' vs', f2(3)
    call check(error, abs(mu1 + mu2) <= 1d-10*abs(mu1), "jump changes sign with the normal")
    if (allocated(error)) goto 9
    call check(error, abs(f1(3) - f2(3)) <= 1d-10*abs(f1(3)), "radiation force unchanged by the normal flip")
    if (allocated(error)) goto 9
    call check(error, abs(x1(3) - x2(3)) <= 1d-10*abs(x1(3)), "diffraction force unchanged by the normal flip")
9   NTHREAD = nthread_saved
    call cleanup
  contains
    subroutine solve_finite(f, x, mu)
      complex(8), intent(out) :: f(6), x(6)
      real(8), intent(out) :: mu
      complex(8), allocatable :: amat(:,:), brmat(:,:), bdmat(:,:)
      integer, allocatable :: ipiv(:)
      integer :: ntot
      V = 0.5d0
      WK = V
      W1 = sqrt(V*9.80665d0)
      H = -1d0
      WVN = 0d0
      AMP = 1d0
      BETA = 0.3d0
      XW = 0d0
      ntot = NELEM + DNELEM
      if (allocated(MXPOT)) deallocate(MXPOT)
      allocate(MXPOT(ntot,7,1), amat(ntot,ntot), brmat(ntot,6), bdmat(ntot,1), ipiv(ntot))
      call DIPOLE_ASSB_LEFT(amat, ipiv, ntot)
      call DIPOLE_ASSB_RBC(brmat, ntot)
      call DIPOLE_SOLVE(amat, ipiv, brmat, ntot, 6, 1)
      call DIPOLE_ASSB_DBC(bdmat(:,1), ntot)
      call DIPOLE_SOLVE(amat, ipiv, bdmat, ntot, 1, 7)
      call DIPOLE_FORCE_SUM(3, f)
      call DIPOLE_FORCE_SUM(7, x)
      mu = real(MXPOT(1,3,1))
    end subroutine solve_finite
  end subroutine test_flip

  subroutine cleanup
    call DIPOLE_RELEASE
    call DIPOLE_RELEASE_ASSB
    if (allocated(MXPOT)) deallocate(MXPOT)
    V = 0d0
  end subroutine cleanup

end module test_dipole_solve
