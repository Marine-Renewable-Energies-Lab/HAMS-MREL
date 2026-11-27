!> Dipole panels on several bodies (notes/specs/2026-11-27-dipole-multibody-design.md, Sec. 8): the reader of the
!> dipole block of HullMesh_N.pnl, the setup of the combined sheet, the option checks, and the multi-body solve
!> (DipoleMulti) against the single-body one (DipoleAssb) for one body, then two bodies: reciprocity of the 12 x 12
!> radiation forces, two bodies far apart against one, and the Haskind relation for each body. Body: the column of
!> radius 1 m and draft 1 m with a plate of radius 1.5 m (test/dipole_bodies.f90), deep water. Uses the solver's module
!> arrays: serial.
module test_dipole_multi

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use dipole_bodies, only : build_body, release_body, SHEET, CLOSED, COLUMN
  use dipole_multi_bodies, only : build_bodies, release_bodies, hull_start
  use PanelMesh_mod, only : NELEM, IPIV, DS, DXYZ_P, XYZ_P, GEONRM_P, NELEM_MULTI, DXYZ_MULTI_COMB, DS_MULTI_COMB, &
                            XYZ_GLOBAL_MULTI_COMB_P, GEONRM_MULTI_COMB, IPIV_MULTI_COMB
  use Inerfs_mod, only : iNELEM, tNELEM
  use Potentials_mod, only : MXPOT, MXPOT_MULTI_COMB
  use DipolePanel_mod, only : DNELEM, DNRM, DDS, DPXYZ_C
  use DipolePanelMulti_mod
  use WaveDyn_mod, only : V, H, WVN, WK, W1, AMP, BETA, NBODY
  use Body_mod, only : XW, XR_MULTI
  use HAMS_mod, only : NTHREAD, ISOLV
  use PatcVelct, only : VINP, DINP
  use CalGreenFunc, only : CALGREEN, CALGREEN_IRR
  use CalGreenFuncMulti, only : CALGREEN_MULTI, CALGREEN_MULTI_IRR
  use DipoleAssb, only : DIPOLE_ASSB_LEFT, DIPOLE_ASSB_RBC, DIPOLE_ASSB_DBC, DIPOLE_SOLVE, DIPOLE_FORCE_SUM, &
                         DIPOLE_RELEASE_ASSB, DIPOLE_ASSB_LEFT_IRR, DIPOLE_ASSB_RHS_IRR
  use DipoleMulti, only : DIPOLE_SETUP_MULTI, DIPOLE_CHECK_OPTIONS_MULTI, DIPOLE_ASSB_LEFT_MULTI, DIPOLE_ASSB_RBC_MULTI, &
                          DIPOLE_ASSB_DBC_MULTI, DIPOLE_SOLVE_MULTI, DIPOLE_ASSB_LEFT_IRR_MULTI, &
                          DIPOLE_ASSB_RHS_IRR_MULTI, DIPOLE_FORCE_SUM_MULTI, DIPOLE_RELEASE_MULTI
  use ReadPanelMeshMulti, only : DIPOLE_READ_PANELS_MULTI

  implicit none
  private
  public :: collect_dipole_multi

  real(8), parameter :: GRAV = 9.80665d0, TDRAFT = 1d0, BPLATE = 1.5d0, PI = 3.141592653589793d0
  integer, parameter :: NT = 24, NZ = 4, NR = 4

contains

  subroutine collect_dipole_multi(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("multi_dipole_read_block", test_read_block), &
      new_unittest("multi_dipole_setup", test_setup), &
      new_unittest("multi_dipole_options", test_options), &
      new_unittest("multi_one_body_equals_single", test_one_body), &
      new_unittest("multi_one_body_lid_equals_single", test_one_body_lid), &
      new_unittest("multi_two_body_reciprocity", test_reciprocity), &
      new_unittest("multi_two_body_far_apart", test_far_apart), &
      new_unittest("multi_two_body_haskind", test_haskind), &
      new_unittest("multi_two_body_lid", test_two_body_lid) &
    ]
  end subroutine collect_dipole_multi

  subroutine set_frequency(nu)
    real(8), intent(in) :: nu
    V = nu
    H = -1d0
    WVN = 0d0
    WK = nu
    W1 = sqrt(nu*GRAV)
  end subroutine set_frequency

  !> The dipole block of body 2 of two (one panel on body 1, two on body 2) read from a file: panel sizes and node
  !> numbers offset by the nodes of the bodies before, in the range of body 2.
  subroutine test_read_block(error)
    type(error_type), allocatable, intent(out) :: error
    integer :: u
    call DIPOLE_MULTI_INIT(2, [1, 2])
    open(newunit=u, status='scratch', action='readwrite')
    write(u, '(a)') '    # Start Definition of Dipole Panels'
    write(u, '(a)') '         1    4       5       6       7       8'
    write(u, '(a)') '         2    3       8       7       9'
    rewind(u)
    call DIPOLE_READ_PANELS_MULTI(u, 2, 100)
    close(u)
    call check(error, all(DOWNER == [1, 2, 2]) .and. DSTART(2) == 2 .and. DEND(2) == 3, 'owners and ranges')
    if (allocated(error)) return
    call check(error, all(DNCN_COMB(2:3) == [4, 3]), 'nodes per panel')
    if (allocated(error)) return
    call check(error, all(DNCON_COMB(2,:) == [105, 106, 107, 108]) .and. all(DNCON_COMB(3,1:3) == [108, 107, 109]), &
               'node numbers offset by the nodes of body 1')
    call DIPOLE_MULTI_RELEASE
  end subroutine test_read_block

  !> Two copies of the column with the sheet at x = -20 and 20 m: the setup takes the sheet nodes out of the combined
  !> node list; each body's centroids, areas and normals are the single body's moved by its offset, and the rotational
  !> normals are about the body's own XR_MULTI.
  subroutine test_setup(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: ok
    character(len=200) :: msg
    integer :: b, k, nd
    real(8) :: dc, dn
    call build_bodies(SHEET, TDRAFT, BPLATE, 0d0, NT, NZ, NR, .false., [-20d0, 20d0])
    nd = DNELEM
    call DIPOLE_SETUP_MULTI(1, 1, .false., 0, ok, msg)
    call check(error, ok, 'setup: '//trim(msg))
    if (allocated(error)) return
    call check(error, DNELEM_TOTAL == 2*nd .and. DSTART(2) == nd + 1 .and. all(DOWNER(nd+1:) == 2), 'ranges')
    if (allocated(error)) return
    dc = 0d0
    dn = 0d0
    do b = 1, 2
      do k = 1, nd
        dc = max(dc, maxval(abs(DPXYZ_C_COMB(DSTART(b)+k-1,:) - DPXYZ_C(k,:) - [XR_MULTI(b,1), 0d0, 0d0])))
        dc = max(dc, abs(DDS_COMB(DSTART(b)+k-1) - DDS(k)))
        dn = max(dn, maxval(abs(DNRM_COMB(DSTART(b)+k-1,:) - DNRM(k,:))))
      end do
    end do
    print '(2x,a,2es10.2)', 'max deviation of centroids and areas, of normals (1:6) from the single body:', dc, dn
    call check(error, dc < 1d-12 .and. dn < 1d-12, 'sheet geometry of each body')
    call DIPOLE_RELEASE_MULTI
    call release_bodies
  end subroutine test_setup

  !> Symmetry, the iterative solver, second order and field points each stop the multi-body dipole path.
  subroutine test_options(error)
    type(error_type), allocatable, intent(out) :: error
    logical :: ok
    character(len=200) :: msg
    call DIPOLE_MULTI_INIT(2, [0, 4])
    call DIPOLE_CHECK_OPTIONS_MULTI(2, 1, .false., 0, ok, msg)
    call check(error, .not. ok .and. index(msg, 'symmetry') > 0 .and. index(msg, 'body 2') > 0, 'symmetry: '//trim(msg))
    if (allocated(error)) return
    call DIPOLE_CHECK_OPTIONS_MULTI(1, 2, .false., 0, ok, msg)
    call check(error, .not. ok .and. index(msg, 'ISOLV = 1') > 0, 'ISOLV: '//trim(msg))
    if (allocated(error)) return
    call DIPOLE_CHECK_OPTIONS_MULTI(1, 1, .true., 0, ok, msg)
    call check(error, .not. ok .and. index(msg, 'first order') > 0, 'second order: '//trim(msg))
    if (allocated(error)) return
    call DIPOLE_CHECK_OPTIONS_MULTI(1, 1, .false., 3, ok, msg)
    call check(error, .not. ok .and. index(msg, 'field-point') > 0, 'field points: '//trim(msg))
    if (allocated(error)) return
    call DIPOLE_CHECK_OPTIONS_MULTI(1, 1, .false., 0, ok, msg)
    call check(error, ok, 'allowed options rejected: '//trim(msg))
    call DIPOLE_MULTI_RELEASE
  end subroutine test_options

  !> Radiation (modes 1-6) and diffraction of the single-body path into MXPOT(1:ntot, 1:7). irsp = 1: with the lid.
  subroutine solve_single(irsp)
    integer, intent(in) :: irsp
    complex(8), allocatable :: af(:,:), cm(:,:), b(:,:), d(:,:)
    integer :: ntot
    ntot = NELEM + DNELEM
    if (allocated(MXPOT)) deallocate(MXPOT)
    if (allocated(IPIV)) deallocate(IPIV)
    allocate(MXPOT(ntot,7,1), IPIV(ntot,1))
    if (irsp == 0) then
      allocate(af(ntot,ntot), b(ntot,6))
      call CALGREEN
      call DIPOLE_ASSB_LEFT(af, IPIV(:,1), ntot)
      call DIPOLE_ASSB_RBC(b, ntot)
      call DIPOLE_SOLVE(af, IPIV(:,1), b, ntot, 6, 1)
      call DIPOLE_ASSB_DBC(b(:,1), ntot)
      call DIPOLE_SOLVE(af, IPIV(:,1), b(:,1:1), ntot, 1, 7)
    else
      allocate(af(ntot+iNELEM,ntot), cm(ntot,ntot), b(ntot+iNELEM,6), d(ntot,6))
      call CALGREEN_IRR
      call DIPOLE_ASSB_LEFT_IRR(af, cm, IPIV(:,1), ntot, iNELEM)
      call DIPOLE_ASSB_RHS_IRR(af, b, d, ntot, iNELEM, 6, 1)
      call DIPOLE_SOLVE(cm, IPIV(:,1), d, ntot, 6, 1)
      call DIPOLE_ASSB_RHS_IRR(af, b(:,1:1), d(:,1:1), ntot, iNELEM, 1, 2)
      call DIPOLE_SOLVE(cm, IPIV(:,1), d(:,1:1), ntot, 1, 7)
    end if
  end subroutine solve_single

  !> Radiation (modes 1 to 6 NBODY) and diffraction of the multi-body path into MXPOT_MULTI_COMB(1:ntot, :).
  subroutine solve_multi(irsp)
    integer, intent(in) :: irsp
    complex(8), allocatable :: af(:,:), cm(:,:), b(:,:), d(:,:)
    integer :: ntot, nh, nl, nb, ib, c
    logical :: ok
    character(len=200) :: msg
    nb = NBODY
    nh = sum(NELEM_MULTI)
    nl = 0
    if (irsp == 1) nl = nb*iNELEM
    tNELEM = nh + nl
    ok = .true.
    msg = ''
    ! once per geometry: the setup renumbers the sheet nodes (build_bodies frees DPXYZ_COMB for a new geometry)
    if (.not. allocated(DPXYZ_COMB)) call DIPOLE_SETUP_MULTI(1, 1, .false., 0, ok, msg)
    if (.not. ok) error stop trim(msg)
    ntot = nh + DNELEM_TOTAL
    if (allocated(MXPOT_MULTI_COMB)) deallocate(MXPOT_MULTI_COMB)
    if (allocated(IPIV_MULTI_COMB)) deallocate(IPIV_MULTI_COMB)
    allocate(MXPOT_MULTI_COMB(ntot,6*nb+1,1), IPIV_MULTI_COMB(ntot,1))
    if (irsp == 0) then
      allocate(af(ntot,ntot), b(ntot,6*nb))
      call CALGREEN_MULTI
      call DIPOLE_ASSB_LEFT_MULTI(af, IPIV_MULTI_COMB(:,1), ntot)
      do ib = 1, nb
        c = 6*(ib - 1)
        call DIPOLE_ASSB_RBC_MULTI(b(:,c+1:c+6), ntot, ib)
      end do
      call DIPOLE_SOLVE_MULTI(af, IPIV_MULTI_COMB(:,1), b, ntot, 6*nb, 1)
      call DIPOLE_ASSB_DBC_MULTI(b(:,1), ntot)
      call DIPOLE_SOLVE_MULTI(af, IPIV_MULTI_COMB(:,1), b(:,1:1), ntot, 1, 6*nb+1)
    else
      allocate(af(ntot+nl,ntot), cm(ntot,ntot), b(ntot+nl,6*nb), d(ntot,6*nb))
      call CALGREEN_MULTI_IRR
      call DIPOLE_ASSB_LEFT_IRR_MULTI(af, cm, IPIV_MULTI_COMB(:,1), ntot, nl)
      do ib = 1, nb
        c = 6*(ib - 1)
        call DIPOLE_ASSB_RHS_IRR_MULTI(af, b(:,c+1:c+6), d(:,c+1:c+6), ntot, nl, 6, 1, ib)
      end do
      call DIPOLE_SOLVE_MULTI(cm, IPIV_MULTI_COMB(:,1), d, ntot, 6*nb, 1)
      call DIPOLE_ASSB_RHS_IRR_MULTI(af, b(:,1:1), d(:,1:1), ntot, nl, 1, 2, 0)
      call DIPOLE_SOLVE_MULTI(cm, IPIV_MULTI_COMB(:,1), d(:,1:1), ntot, 1, 6*nb+1)
    end if
  end subroutine solve_multi

  !> f(6(m-1)+k, 6(n-1)+l): the radiation force sum on body m, mode k, of the potential of body n, mode l (hull and sheet
  !> of body m), from MXPOT_MULTI_COMB. A = rho Re f, B = rho omega Im f.
  subroutine multi_forces(f)
    complex(8), intent(out) :: f(:,:)
    complex(8) :: fd(6)
    integer :: m, n, l, k, j0, j1, col
    do m = 1, NBODY
      j0 = hull_start(m)
      j1 = j0 + NELEM_MULTI(m) - 1
      do n = 1, NBODY
        do l = 1, 6
          col = 6*(n - 1) + l
          call DIPOLE_FORCE_SUM_MULTI(m, col, fd)
          do k = 1, 6
            f(6*(m-1)+k, col) = sum(MXPOT_MULTI_COMB(j0:j1,col,1)*DXYZ_MULTI_COMB(j0:j1,k)*DS_MULTI_COMB(j0:j1)) + fd(k)
          end do
        end do
      end do
    end do
  end subroutine multi_forces

  !> The same sums for the single-body path (MXPOT, DipoleAssb).
  subroutine single_forces(f)
    complex(8), intent(out) :: f(6,6)
    complex(8) :: fd(6)
    integer :: j, k
    do j = 1, 6
      call DIPOLE_FORCE_SUM(j, fd)
      do k = 1, 6
        f(k,j) = sum(MXPOT(1:NELEM,j,1)*DXYZ_P(:,k)*DS) + fd(k)
      end do
    end do
  end subroutine single_forces

  !> One body through the multi-body path (NBODY = 1 in the multi arrays) gives the potentials, sheet jumps and force
  !> sums of the single-body path, radiation and diffraction, nu = 1 /m.
  subroutine test_one_body(error)
    type(error_type), allocatable, intent(out) :: error
    call one_body_case(0, error)
  end subroutine test_one_body

  !> As test_one_body with the lid (IRSP = 1), at the first heave irregular frequency of the column (nu = 2.445 /m).
  subroutine test_one_body_lid(error)
    type(error_type), allocatable, intent(out) :: error
    call one_body_case(1, error)
  end subroutine test_one_body_lid

  subroutine one_body_case(irsp, error)
    integer, intent(in) :: irsp
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: fs(6,6), fm(6,6)
    real(8) :: dx, df
    integer :: ntot, nthread_saved
    nthread_saved = NTHREAD
    NTHREAD = 4
    ISOLV = 1
    if (irsp == 0) then
      call set_frequency(1d0)
    else
      call set_frequency(2.445d0)
    end if
    AMP = 1d0
    BETA = 30d0*PI/180d0
    XW = 0d0
    call build_bodies(SHEET, TDRAFT, BPLATE, 0d0, NT, NZ, NR, irsp == 1, [0d0])
    ntot = NELEM + DNELEM
    call solve_single(irsp)
    call single_forces(fs)
    call solve_multi(irsp)
    call multi_forces(fm)
    dx = maxval(abs(MXPOT_MULTI_COMB(1:ntot,1:7,1) - MXPOT(1:ntot,1:7,1)))/maxval(abs(MXPOT(1:ntot,1:7,1)))
    df = maxval(abs(fm - fs))/maxval(abs(fs))
    NTHREAD = nthread_saved
    print '(2x,a,i1,a,2es10.2)', 'IRSP = ', irsp, ': multi against single, max rel diff of potentials, forces', dx, df
    call DIPOLE_RELEASE_ASSB
    call DIPOLE_RELEASE_MULTI
    call release_bodies
    call check(error, dx < 1d-9 .and. df < 1d-9, 'one body: the multi path equals the single path')
  end subroutine one_body_case

  !> Two columns with sheets 8 m apart (gap 5 m between the plates), nu = 1 /m: the 12 x 12 force sums f of the radiation
  !> potentials are symmetric (A and B reciprocal), the coupling blocks included.
  subroutine test_reciprocity(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: f(12,12)
    real(8) :: asym, acoup, scale
    integer :: nthread_saved
    nthread_saved = NTHREAD
    NTHREAD = 4
    ISOLV = 1
    call set_frequency(1d0)
    AMP = 1d0
    BETA = 0d0
    XW = 0d0
    call build_bodies(SHEET, TDRAFT, BPLATE, 0d0, NT, NZ, NR, .false., [-4d0, 4d0])
    call solve_multi(0)
    call multi_forces(f)
    NTHREAD = nthread_saved
    scale = maxval(abs(f))
    asym = maxval(abs(f - transpose(f)))/scale
    acoup = maxval(abs(f(1:6,7:12) - transpose(f(7:12,1:6))))/maxval(abs(f(1:6,7:12)))
    print '(2x,a,2es10.2,a,2es12.4)', 'max |f - f^T| / max |f|, coupling block alone:', asym, acoup, &
          '; f(3,9), f(9,3):', abs(f(3,9)), abs(f(9,3))
    print '(2x,a,4es14.6)', 'f(3,3), f(3,9) (Re, Im; the integration suite dipole_multi checks the program against them):', &
          f(3,3), f(3,9)
    call DIPOLE_RELEASE_MULTI
    call release_bodies
    call check(error, asym < 2d-2 .and. acoup < 5d-2, 'reciprocity of the 12 x 12 radiation forces')
  end subroutine test_reciprocity

  !> Two columns with sheets 200 m apart against one alone, nu = 1 /m: each diagonal 6 x 6 block of the force sums is
  !> close to the single body's, and the coupling blocks are as large as for the same pair of closed 0.02 m plates
  !> (hull panels only), whose interaction is that of the base code.
  subroutine test_far_apart(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: f(12,12), fs(6,6), fc(12,12), fcs(6,6)
    real(8) :: dd, dc, cd, cc
    integer :: nthread_saved
    nthread_saved = NTHREAD
    NTHREAD = 4
    ISOLV = 1
    call set_frequency(1d0)
    call build_bodies(SHEET, TDRAFT, BPLATE, 0d0, NT, NZ, NR, .false., [-100d0, 100d0])
    call solve_single(0)
    call single_forces(fs)
    call solve_multi(0)
    call multi_forces(f)
    call DIPOLE_RELEASE_ASSB
    call DIPOLE_RELEASE_MULTI
    call build_bodies(CLOSED, TDRAFT, BPLATE, 0.02d0, NT, NZ, NR, .false., [-100d0, 100d0])
    call solve_single(0)
    call single_forces(fcs)
    call solve_multi(0)
    call multi_forces(fc)
    NTHREAD = nthread_saved
    dd = max(maxval(abs(f(1:6,1:6) - fs)), maxval(abs(f(7:12,7:12) - fs)))/maxval(abs(fs))
    dc = max(maxval(abs(f(1:6,7:12))), maxval(abs(f(7:12,1:6))))/maxval(abs(fs))
    cd = max(maxval(abs(fc(1:6,1:6) - fcs)), maxval(abs(fc(7:12,7:12) - fcs)))/maxval(abs(fcs))
    cc = max(maxval(abs(fc(1:6,7:12))), maxval(abs(fc(7:12,1:6))))/maxval(abs(fcs))
    print '(2x,a,2es10.2,a,2es10.2)', 'far apart, sheets: diagonal change, coupling', dd, dc, '; closed plates:', cd, cc
    call DIPOLE_RELEASE_ASSB
    call DIPOLE_RELEASE_MULTI
    call release_bodies
    call check(error, dd < 1d-2 .and. dd < 2d0*cd + 1d-4, 'far apart: diagonal blocks as for one body')
    if (allocated(error)) return
    ! Wave coupling decays only like 1/sqrt(k R) (about 4 % at 200 m here), so the coupling of the sheets is held
    ! against that of the closed plates, not against zero.
    call check(error, abs(dc - cc) <= 0.2d0*cc, 'far apart: coupling blocks as for the closed plates')
  end subroutine test_far_apart

  !> Exciting force of each body, direct (pressure of the incident and diffracted potentials on the body's hull and
  !> sheet) against the Haskind relation (incident potential and the body's radiation potentials over all hulls and
  !> sheets), two columns with sheets 8 m apart, heading 30 deg, nu = 1 /m, modes 1, 3 and 5. The error of each body is
  !> held against that of the single body on the same mesh.
  subroutine test_haskind(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: xd(12), xh(12), sd(6), sh(6)
    real(8) :: em, es
    integer :: m, i, k, nthread_saved
    integer, parameter :: modes(3) = [1, 3, 5]
    character(len=200) :: msg
    nthread_saved = NTHREAD
    NTHREAD = 4
    ISOLV = 1
    call set_frequency(1d0)
    AMP = 1d0
    BETA = 30d0*PI/180d0
    XW = 0d0
    call build_bodies(SHEET, TDRAFT, BPLATE, 0d0, NT, NZ, NR, .false., [-4d0, 4d0])
    call solve_single(0)
    call haskind_single(sd, sh)
    call solve_multi(0)
    call haskind_multi(xd, xh)
    NTHREAD = nthread_saved
    call DIPOLE_RELEASE_ASSB
    call DIPOLE_RELEASE_MULTI
    call release_bodies
    do m = 1, 2
      do i = 1, 3
        k = modes(i)
        em = abs(xd(6*(m-1)+k) - xh(6*(m-1)+k))/abs(xd(6*(m-1)+k))
        es = abs(sd(k) - sh(k))/abs(sd(k))
        write(msg, '(a,i1,a,i1,a,4es12.4,a,es9.2,a,es9.2)') 'body ', m, ', mode ', k, ': direct, Haskind', &
             xd(6*(m-1)+k), xh(6*(m-1)+k), '; rel diff', em, '; single body', es
        print '(2x,a)', trim(msg)
        call check(error, em <= max(1d-2, 2d0*es), trim(msg))
        if (allocated(error)) return
      end do
    end do
  end subroutine test_haskind

  !> Radiation of the base multi-body path (AssbMatxMulti / AssbMatx_irrMulti, hull panels only) into
  !> MXPOT_MULTI_COMB(1:nh, 1:6 NBODY).
  subroutine solve_base_multi(irsp)
    use AssbMatxMulti, only : ASSB_LEFT_MULTI, ASSB_RBC_MULTI, RADIATION_SOLVER_MULTI
    use AssbMatx_irrMulti, only : ASSB_LEFT_IRR_MULTI, ASSB_RBC_IRR_MULTI, RADIATION_SOLVER_IRR_MULTI
    integer, intent(in) :: irsp
    complex(8), allocatable :: amat(:,:,:), brmat(:,:,:), cmat(:,:,:), drmat(:,:,:)
    integer :: nh, nl, nb, ib, j0, j1, tn
    nb = NBODY
    nh = sum(NELEM_MULTI)
    if (allocated(MXPOT_MULTI_COMB)) deallocate(MXPOT_MULTI_COMB)
    if (allocated(IPIV_MULTI_COMB)) deallocate(IPIV_MULTI_COMB)
    allocate(MXPOT_MULTI_COMB(nh,6*nb+1,1), IPIV_MULTI_COMB(nh,1))
    if (irsp == 0) then
      tNELEM = nh
      allocate(amat(nh,nh,1), brmat(nh,6*nb,1))
      brmat = (0d0, 0d0)
      call CALGREEN_MULTI
      call ASSB_LEFT_MULTI(amat, IPIV_MULTI_COMB, nh, 1)
      do ib = 1, nb
        j0 = hull_start(ib)
        j1 = j0 + NELEM_MULTI(ib) - 1
        call ASSB_RBC_MULTI(brmat(:,6*(ib-1)+1:6*ib,:), nh, 1, j0, j1)
      end do
      call RADIATION_SOLVER_MULTI(amat, brmat, IPIV_MULTI_COMB, MXPOT_MULTI_COMB, nh, 1)
    else
      nl = nb*iNELEM
      tn = nh + nl
      tNELEM = tn
      allocate(amat(tn,tn,1), brmat(tn,6*nb,1), cmat(nh,nh,1), drmat(nh,6*nb,1))
      brmat = (0d0, 0d0)
      drmat = (0d0, 0d0)
      call CALGREEN_MULTI_IRR
      call ASSB_LEFT_IRR_MULTI(amat, cmat, IPIV_MULTI_COMB, nh, nl, tn, 1)
      do ib = 1, nb
        j0 = hull_start(ib)
        j1 = j0 + NELEM_MULTI(ib) - 1
        call ASSB_RBC_IRR_MULTI(brmat(:,6*(ib-1)+1:6*ib,:), drmat(:,6*(ib-1)+1:6*ib,:), amat, nh, nl, tn, 1, j0, j1)
      end do
      call RADIATION_SOLVER_IRR_MULTI(cmat, drmat, IPIV_MULTI_COMB, MXPOT_MULTI_COMB, nh, 1)
    end if
  end subroutine solve_base_multi

  !> The lids (IRSP = 1) on two columns 8 m apart, nu = 1 /m (away from the columns' irregular frequencies).
  !> (a) Hull panels only (the plate as a closed 0.02 m plate): DipoleMulti, without and with the lids, gives the base
  !>     multi-body path (AssbMatxMulti / AssbMatx_irrMulti) to 1e-9. (A thin closed plate is no reference for the lid
  !>     itself: its interior makes the lid move the forces a lot, docs/input_files.md, Dipole panels.)
  !> (b) The lids move the coupling block of the sheets' forces no more than that of the plain columns in the base path
  !>     (within a factor 2, plus 1e-3 of the largest force), as lid_dipole_removes_irregular does for one body.
  subroutine test_two_body_lid(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: f0(12,12), f1(12,12), c0(12,12), c1(12,12), d0(12,12), d1(12,12)
    real(8) :: e0, e1, ds, dc
    integer :: nthread_saved
    character(len=240) :: msg
    nthread_saved = NTHREAD
    NTHREAD = 4
    ISOLV = 1
    call set_frequency(1d0)
    ! (a) closed plates: base path against DipoleMulti
    call build_bodies(CLOSED, TDRAFT, BPLATE, 0.02d0, NT, NZ, NR, .true., [-4d0, 4d0])
    call solve_base_multi(0)
    call multi_forces(c0)
    call solve_base_multi(1)
    call multi_forces(c1)
    call solve_multi(0)
    call multi_forces(d0)
    call solve_multi(1)
    call multi_forces(d1)
    call DIPOLE_RELEASE_MULTI
    e0 = maxval(abs(d0 - c0))/maxval(abs(c0))
    e1 = maxval(abs(d1 - c1))/maxval(abs(c1))
    ! (b) sheets (DipoleMulti) and plain columns (base path), without and with the lids
    call build_bodies(SHEET, TDRAFT, BPLATE, 0d0, NT, NZ, NR, .true., [-4d0, 4d0])
    call solve_multi(0)
    call multi_forces(f0)
    call solve_multi(1)
    call multi_forces(f1)
    call DIPOLE_RELEASE_MULTI
    ds = maxval(abs(f1(1:6,7:12) - f0(1:6,7:12)))/maxval(abs(f0))
    call build_bodies(COLUMN, TDRAFT, BPLATE, 0d0, NT, NZ, NR, .true., [-4d0, 4d0])
    call solve_base_multi(0)
    call multi_forces(c0)
    call solve_base_multi(1)
    call multi_forces(c1)
    dc = maxval(abs(c1(1:6,7:12) - c0(1:6,7:12)))/maxval(abs(c0))
    NTHREAD = nthread_saved
    call release_bodies
    write(msg, '(a,2es10.2)') 'closed plates, DipoleMulti against the base path, max rel diff without, with lids:', e0, e1
    print '(2x,a)', trim(msg)
    call check(error, e0 < 1d-9 .and. e1 < 1d-9, trim(msg))
    if (allocated(error)) return
    write(msg, '(a,2es10.2,a,4es12.4)') 'lid change of the coupling block / max |f|: sheets, columns', ds, dc, &
         '; f(3,9) of the sheets without, with lids:', f0(3,9), f1(3,9)
    print '(2x,a)', trim(msg)
    call check(error, ds <= 2d0*dc + 1d-3, trim(msg))
  end subroutine test_two_body_lid

  subroutine haskind_single(xd, xh)
    use DipolePanel_mod, only : DNRM
    complex(8), intent(out) :: xd(6), xh(6)
    complex(8) :: phii, dx, dy, dz, dn
    integer :: i
    xd = 0d0
    xh = 0d0
    do i = 1, NELEM
      phii = VINP(XYZ_P(i,1), XYZ_P(i,2), XYZ_P(i,3), XW(1), XW(2), BETA)
      call DINP(XYZ_P(i,1), XYZ_P(i,2), XYZ_P(i,3), XW(1), XW(2), BETA, dx, dy, dz)
      dn = dx*GEONRM_P(i,1) + dy*GEONRM_P(i,2) + dz*GEONRM_P(i,3)
      xd = xd + (phii + MXPOT(i,7,1))*DXYZ_P(i,:)*DS(i)
      xh = xh + (phii*DXYZ_P(i,:) - MXPOT(i,1:6,1)*dn)*DS(i)
    end do
    do i = 1, DNELEM
      call DINP(DPXYZ_C(i,1), DPXYZ_C(i,2), DPXYZ_C(i,3), XW(1), XW(2), BETA, dx, dy, dz)
      dn = dx*DNRM(i,1) + dy*DNRM(i,2) + dz*DNRM(i,3)
      xd = xd - MXPOT(NELEM+i,7,1)*DNRM(i,1:6)*DDS(i)
      xh = xh + MXPOT(NELEM+i,1:6,1)*dn*DDS(i)
    end do
  end subroutine haskind_single

  !> Body m: direct X over its hull and sheet; Haskind over all hulls (n_k only on body m's own) and all sheets.
  subroutine haskind_multi(xd, xh)
    complex(8), intent(out) :: xd(:), xh(:)
    complex(8) :: phii, dx, dy, dz, dn
    integer :: i, m, nh, j0, j1, c, nb
    nb = NBODY
    nh = sum(NELEM_MULTI)
    xd = 0d0
    xh = 0d0
    do m = 1, nb
      c = 6*(m - 1)
      j0 = hull_start(m)
      j1 = j0 + NELEM_MULTI(m) - 1
      do i = 1, nh
        associate (p => XYZ_GLOBAL_MULTI_COMB_P(i,:))
          phii = VINP(p(1), p(2), p(3), XW(1), XW(2), BETA)
          call DINP(p(1), p(2), p(3), XW(1), XW(2), BETA, dx, dy, dz)
        end associate
        dn = dx*GEONRM_MULTI_COMB(i,1) + dy*GEONRM_MULTI_COMB(i,2) + dz*GEONRM_MULTI_COMB(i,3)
        if (i >= j0 .and. i <= j1) then
          xd(c+1:c+6) = xd(c+1:c+6) + (phii + MXPOT_MULTI_COMB(i,6*nb+1,1))*DXYZ_MULTI_COMB(i,:)*DS_MULTI_COMB(i)
          xh(c+1:c+6) = xh(c+1:c+6) + phii*DXYZ_MULTI_COMB(i,:)*DS_MULTI_COMB(i)
        end if
        xh(c+1:c+6) = xh(c+1:c+6) - MXPOT_MULTI_COMB(i,c+1:c+6,1)*dn*DS_MULTI_COMB(i)
      end do
      do i = 1, DNELEM_TOTAL
        call DINP(DPXYZ_C_COMB(i,1), DPXYZ_C_COMB(i,2), DPXYZ_C_COMB(i,3), XW(1), XW(2), BETA, dx, dy, dz)
        dn = dx*DNRM_COMB(i,1) + dy*DNRM_COMB(i,2) + dz*DNRM_COMB(i,3)
        if (DOWNER(i) == m) xd(c+1:c+6) = xd(c+1:c+6) - MXPOT_MULTI_COMB(nh+i,6*nb+1,1)*DNRM_COMB(i,1:6)*DDS_COMB(i)
        xh(c+1:c+6) = xh(c+1:c+6) + MXPOT_MULTI_COMB(nh+i,c+1:c+6,1)*dn*DDS_COMB(i)
      end do
    end do
  end subroutine haskind_multi

end module test_dipole_multi
