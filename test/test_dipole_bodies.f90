!> Dipole panels with a hull (plan Tasks 5 and 6): a column of radius 1 m and draft 1 m with a plate of radius 1.5 m at
!> its bottom (test/dipole_bodies.f90). The plate as a dipole sheet against the plate as a closed body that gets thinner
!> (Review Focus 4), and the lid (IRSP = 1) with dipole columns. Deep water. Uses the solver's module arrays: serial.
module test_dipole_bodies

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use dipole_bodies, only : build_body, release_body, SHEET, CLOSED, COLUMN
  use PanelMesh_mod, only : NELEM, IPIV, DS, DXYZ_P
  use Inerfs_mod, only : iNELEM
  use Potentials_mod, only : MXPOT
  use DipoleMesh_mod, only : DNELEM
  use WaveDyn_mod, only : V, H, WVN, WK, W1
  use HAMS_mod, only : NTHREAD, ISOLV
  use CalGreenFunc, only : CALGREEN, CALGREEN_IRR
  use AssbMatx, only : ASSB_LEFT, ASSB_RBC, RADIATION_SOLVER
  use AssbMatx_irr, only : ASSB_LEFT_IRR, ASSB_RBC_IRR, RADIATION_SOLVER_IRR
  use DipoleAssb, only : DIPOLE_ASSB_LEFT, DIPOLE_ASSB_RBC, DIPOLE_SOLVE, DIPOLE_FORCE_SUM, DIPOLE_RELEASE_ASSB, &
                         DIPOLE_ASSB_LEFT_IRR, DIPOLE_ASSB_RHS_IRR, DIPOLE_LID_RESIDUAL

  implicit none
  private
  public :: collect_dipole_bodies

  real(8), parameter :: GRAV = 9.80665d0, TDRAFT = 1d0, BPLATE = 1.5d0
  integer, parameter :: NT = 24, NZ = 4, NR = 4
  real(8) :: lid_residual = 0d0

contains

  subroutine collect_dipole_bodies(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("thin_plate_converges_to_dipole", test_thin_plate), &
      new_unittest("lid_dipole_removes_irregular", test_lid), &
      new_unittest("haskind_dipole", test_haskind) &
    ]
  end subroutine collect_dipole_bodies

  subroutine set_frequency(nu)
    real(8), intent(in) :: nu
    V = nu
    H = -1d0
    WVN = 0d0
    WK = nu
    if (nu > 0d0) then
      W1 = sqrt(nu*GRAV)
    else
      W1 = 0d0
    end if
  end subroutine set_frequency

  !> f(k, j) = sum phi_j n_k dS over hull and sheet (the sheet as sum -mu n_D dS), k, j = 1..6: A = rho Re f,
  !> B = rho omega Im f (RFORCE). irsp: with the lid.
  subroutine radiation_sums(irsp, f)
    integer, intent(in) :: irsp
    complex(8), intent(out) :: f(6,6)
    complex(8), allocatable :: amat(:,:,:), brmat(:,:,:), cmat(:,:,:), drmat(:,:,:), af(:,:), cm(:,:), b(:,:), d(:,:)
    complex(8) :: fd(6)
    integer :: ntot, j, k, tn
    ntot = NELEM + DNELEM
    if (allocated(MXPOT)) deallocate(MXPOT)
    if (allocated(IPIV)) deallocate(IPIV)
    allocate(MXPOT(ntot,7,1), IPIV(ntot,1))
    if (DNELEM == 0) then
      if (irsp == 0) then
        allocate(amat(NELEM,NELEM,1), brmat(NELEM,6,1))
        call CALGREEN
        call ASSB_LEFT(amat, IPIV, NELEM, 1)
        call ASSB_RBC(brmat, NELEM, 1)
        call RADIATION_SOLVER(amat, brmat, IPIV, MXPOT, NELEM, 1)
      else
        tn = NELEM + iNELEM
        allocate(amat(tn,tn,1), brmat(tn,6,1), cmat(NELEM,NELEM,1), drmat(NELEM,6,1))
        call CALGREEN_IRR
        call ASSB_LEFT_IRR(amat, cmat, IPIV, NELEM, tn, 1)
        call ASSB_RBC_IRR(brmat, drmat, amat, NELEM, tn, 1)
        call RADIATION_SOLVER_IRR(cmat, drmat, IPIV, MXPOT, NELEM, 1)
      end if
    else
      if (irsp == 0) then
        allocate(af(ntot,ntot), b(ntot,6))
        call CALGREEN
        call DIPOLE_ASSB_LEFT(af, IPIV(:,1), ntot)
        call DIPOLE_ASSB_RBC(b, ntot)
        call DIPOLE_SOLVE(af, IPIV(:,1), b, ntot, 6, 1)
      else
        allocate(af(ntot+iNELEM,ntot), cm(ntot,ntot), b(ntot+iNELEM,6), d(ntot,6))
        call CALGREEN_IRR
        call DIPOLE_ASSB_LEFT_IRR(af, cm, IPIV(:,1), ntot, iNELEM)
        call DIPOLE_ASSB_RHS_IRR(af, b, d, ntot, iNELEM, 6, 1)
        call DIPOLE_SOLVE(cm, IPIV(:,1), d, ntot, 6, 1)
        lid_residual = DIPOLE_LID_RESIDUAL(af, b(:,3), MXPOT(1:ntot,3,1), ntot, iNELEM)
      end if
    end if
    do j = 1, 6
      do k = 1, 6
        f(k,j) = sum(MXPOT(1:NELEM,j,1)*DXYZ_P(:,k)*DS)
      end do
      call DIPOLE_FORCE_SUM(j, fd)
      f(:,j) = f(:,j) + fd
    end do
  end subroutine radiation_sums

  !> Review Focus 4: the plate as a closed body 0.2, 0.1, 0.05 and 0.02 m thick against the sheet, at the zero-frequency
  !> limit and at nu = 1 /m: A33, A55 and B33 of the closed body approach those of the sheet as it gets thinner.
  subroutine test_thin_plate(error)
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: ths(4) = [0.2d0, 0.1d0, 0.05d0, 0.02d0], nus(2) = [0d0, 1d0]
    complex(8) :: fs(6,6), fc(6,6)
    real(8) :: dev(4,2,3)
    integer :: i, iv, nthread_saved
    character(len=200) :: msg
    nthread_saved = NTHREAD
    NTHREAD = 4
    ISOLV = 1
    do iv = 1, 2
      call set_frequency(nus(iv))
      call build_body(SHEET, TDRAFT, BPLATE, 0d0, NT, NZ, NR, .false.)
      call radiation_sums(0, fs)
      write(msg, '(a,f4.1,a,3es13.5)') 'nu =', nus(iv), ' sheet: Re f33, Re f55, Im f33 =', real(fs(3,3)), real(fs(5,5)), &
           aimag(fs(3,3))
      print '(2x,a)', trim(msg)
      do i = 1, 4
        call build_body(CLOSED, TDRAFT, BPLATE, ths(i), NT, NZ, NR, .false.)
        call radiation_sums(0, fc)
        dev(i,iv,1) = abs(real(fc(3,3)) - real(fs(3,3)))/abs(real(fs(3,3)))
        dev(i,iv,2) = abs(real(fc(5,5)) - real(fs(5,5)))/abs(real(fs(5,5)))
        dev(i,iv,3) = abs(aimag(fc(3,3)) - aimag(fs(3,3)))/max(abs(aimag(fs(3,3))), 1d-300)
        write(msg, '(a,f4.1,a,f5.2,a,3es13.5,a,3es10.2)') 'nu =', nus(iv), ' closed t =', ths(i), ':', real(fc(3,3)), &
             real(fc(5,5)), aimag(fc(3,3)), '; rel diff', dev(i,iv,:)
        print '(2x,a)', trim(msg)
      end do
    end do
    NTHREAD = nthread_saved
    call DIPOLE_RELEASE_ASSB
    call release_body
    do iv = 1, 2
      call check(error, dev(4,iv,1) < dev(1,iv,1) .and. dev(4,iv,2) < dev(1,iv,2), "closed plate approaches the sheet")
      if (allocated(error)) return
      call check(error, dev(4,iv,1) < 1d-2 .and. dev(4,iv,2) < 1d-2, "closed plate 0.02 m thick: A33, A55 within 1 %")
      if (allocated(error)) return
    end do
  end subroutine test_thin_plate

  !> Plan Task 6: the column with the sheet, with and without the lid. Away from the irregular frequencies the two
  !> agree; at the first heave irregular frequency of the column (k = 2.405, nu = k coth(k T) = 2.445 /m) the run
  !> without the lid is off and the lid run is not; the lid rows are satisfied by the least-squares solution.
  subroutine test_lid(error)
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: nus(5) = [1d0, 2.2d0, 2.445d0, 2.7d0, 3.2d0]
    complex(8) :: f0(6,6), f1(6,6)
    real(8) :: a0(5), a1(5), b0(5), b1(5), res(5), ca, cb, smooth
    integer :: iv, nthread_saved
    nthread_saved = NTHREAD
    NTHREAD = 4
    ISOLV = 1
    do iv = 1, 5
      call set_frequency(nus(iv))
      call build_body(SHEET, TDRAFT, BPLATE, 0d0, NT, NZ, NR, .true.)
      call radiation_sums(0, f0)
      call radiation_sums(1, f1)
      a0(iv) = real(f0(3,3))
      a1(iv) = real(f1(3,3))
      b0(iv) = aimag(f0(3,3))
      b1(iv) = aimag(f1(3,3))
      res(iv) = lid_residual
      print '(2x,a,f6.3,a,2es13.5,a,2es13.5,a,es10.2)', 'nu =', nus(iv), ' Re f33 no lid, lid:', a0(iv), a1(iv), &
            '; Im f33:', b0(iv), b1(iv), '; lid residual', res(iv)
    end do
    ! the base code on the column alone at nu = 1, same mesh and lid: how far its own lid moves A33 and B33
    call set_frequency(nus(1))
    call build_body(COLUMN, TDRAFT, BPLATE, 0d0, NT, NZ, NR, .true.)
    call radiation_sums(0, f0)
    call radiation_sums(1, f1)
    ca = abs(real(f1(3,3)) - real(f0(3,3)))/abs(real(f0(3,3)))
    cb = abs(aimag(f1(3,3)) - aimag(f0(3,3)))/abs(aimag(f0(3,3)))
    print '(2x,a,2es10.2,a,2es10.2)', 'nu = 1, lid against no lid, rel diff of Re f33, Im f33: column (base code)', &
          ca, cb, '; column with sheet', abs(a1(1) - a0(1))/abs(a0(1)), abs(b1(1) - b0(1))/abs(b0(1))
    NTHREAD = nthread_saved
    call DIPOLE_RELEASE_ASSB
    call release_body
    call check(error, abs(a1(1) - a0(1)) < 1d-2*abs(a0(1)), "A33: lid and no lid within 1 % at nu = 1")
    if (allocated(error)) return
    call check(error, abs(b1(1) - b0(1))/abs(b0(1)) < 2d0*cb + 1d-2, &
               "B33: lid against no lid at nu = 1 no further apart than twice the base code's own lid on the column")
    if (allocated(error)) return
    ! at the irregular frequency: no spike with the lid (within 1 % of the line through 2.2 and 2.7), a spike without
    smooth = a1(2) + (a1(4) - a1(2))*(nus(3) - nus(2))/(nus(4) - nus(2))
    call check(error, abs(a1(3) - smooth) < 1d-2*abs(smooth), "lid: no spike in A33 at nu = 2.445")
    if (allocated(error)) return
    smooth = a0(2) + (a0(4) - a0(2))*(nus(3) - nus(2))/(nus(4) - nus(2))
    call check(error, abs(a0(3) - smooth) > 2d-2*abs(smooth), "no lid: A33 is off at nu = 2.445 (the irregular frequency)")
    if (allocated(error)) return
    call check(error, maxval(res) < 1d-2, "lid rows satisfied to 1 % (rms of residual over row terms)")
  end subroutine test_lid

  !> Plan Task 5: the exciting force from the diffraction potential against Haskind's relation (Newman 1962) from the
  !> radiation potentials, nu = 1 /m, heading 30 degrees, surge, heave and pitch. Without the factor i omega rho:
  !> direct sum_hull (phi_I + phi_S) n_k dS + sum_sheet (-mu_S) n_k dS; Haskind sum_hull (phi_I n_k - phi_k dphi_I/dn) dS
  !> + sum_sheet mu_k dphi_I/dn_D dS (design note Sec. 4: phi_I has no jump, dphi_I/dn flips with the face). The two differ
  !> by the discretisation; the same check on the closed plate 0.02 m thick (hull panels only) sets the scale: the sheet
  !> must agree within 1 % or within twice the closed plate's difference.
  subroutine test_haskind(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: xd(6), xh(6), cd(6), ch(6)
    real(8) :: es, ec
    integer :: i, k, nthread_saved
    integer, parameter :: modes(3) = [1, 3, 5]
    character(len=200) :: msg
    nthread_saved = NTHREAD
    NTHREAD = 4
    ISOLV = 1
    call haskind_pair(SHEET, 0d0, xd, xh)
    call haskind_pair(CLOSED, 0.02d0, cd, ch)
    NTHREAD = nthread_saved
    call DIPOLE_RELEASE_ASSB
    call release_body
    do i = 1, 3
      k = modes(i)
      es = abs(xd(k) - xh(k))/abs(xd(k))
      ec = abs(cd(k) - ch(k))/abs(cd(k))
      write(msg, '(a,i1,a,4es13.5,a,es10.2,a,es10.2)') 'Haskind, mode ', k, ': sheet direct, Haskind', xd(k), xh(k), &
           '; rel diff', es, '; closed 0.02 m plate rel diff', ec
      print '(2x,a)', trim(msg)
      call check(error, es <= max(1d-2, 2d0*ec), trim(msg))
      if (allocated(error)) return
    end do
  end subroutine test_haskind

  subroutine haskind_pair(kind, th, xd, xh)
    use PanelMesh_mod, only : XYZ_P, GEONRM_P
    use DipoleMesh_mod, only : DPXYZ_C, DNRM, DDS
    use WaveDyn_mod, only : AMP, BETA
    use Body_mod, only : XW
    use PatcVelct, only : VINP, DINP
    use DipoleAssb, only : DIPOLE_ASSB_DBC
    integer, intent(in) :: kind
    real(8), intent(in) :: th
    complex(8), intent(out) :: xd(6), xh(6)
    complex(8) :: phii, dx, dy, dz, dn
    complex(8), allocatable :: af(:,:), b(:,:), bd(:,:)
    integer :: ntot, i
    call set_frequency(1d0)
    AMP = 1d0
    BETA = 30d0*3.141592653589793d0/180d0
    XW = 0d0
    call build_body(kind, TDRAFT, BPLATE, th, NT, NZ, NR, .false.)
    ntot = NELEM + DNELEM
    if (allocated(MXPOT)) deallocate(MXPOT)
    if (allocated(IPIV)) deallocate(IPIV)
    allocate(MXPOT(ntot,7,1), IPIV(ntot,1), af(ntot,ntot), b(ntot,6), bd(ntot,1))
    call CALGREEN
    call DIPOLE_ASSB_LEFT(af, IPIV(:,1), ntot)
    call DIPOLE_ASSB_RBC(b, ntot)
    call DIPOLE_SOLVE(af, IPIV(:,1), b, ntot, 6, 1)
    call DIPOLE_ASSB_DBC(bd(:,1), ntot)
    call DIPOLE_SOLVE(af, IPIV(:,1), bd, ntot, 1, 7)
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
  end subroutine haskind_pair

end module test_dipole_bodies
