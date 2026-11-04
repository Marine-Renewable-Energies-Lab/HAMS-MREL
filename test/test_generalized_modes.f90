!> Unit tests for generalized modes (gen_mod_<body>_<dof>.txt) in DOF 1 to 3.
!> The file's values replace DXYZ_P(:,DOF) and are that DOF's normal, used as given. Code that needs the panel's direction
!> must not see them: dG/dn on the left side, dphi_I/dn in the diffraction right side, and the field-point kernel, incident
!> term and on-body test. So each of those must return the same value with and without a mode, and the radiation right
!> side must use the mode's value. Kernel arrays and potentials are synthetic: no Green function or solve is involved,
!> except for the panel integrals of the field-point routines.
!> These tests set module globals, so the suite runs serially (see unit_tester.f90).
module test_generalized_modes
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use testing_utilities, only : get_running_exe_path
  use HAMS_mod
  use WaveDyn_mod
  use Body_mod
  use PanelMesh_mod
  use Const_mod
  use Inerfs_mod
  use Potentials_mod
  use IO, only : CreateDirectory
  use ReadPanelMesh, only : CalNormals
  use ReadPanelMeshMulti, only : CalNormalsMulti
  use CalGreenFuncMulti, only : BuildCombinedArrays
  use BodyIntgr, only : BODINT_LEFT, RBC_RIGHT, DBC_RIGHT
  use BodyIntgr_irr, only : BODINT_LEFT_IRR, RBC_IRR, DBC_IRR
  use BodyIntgrMulti, only : BODINT_MULTI_LEFT, RBC_MULTI_RIGHT, DBC_MULTI_RIGHT
  use BodyIntgr_irrMulti, only : BODINT_LEFT_IRR_MULTI, RBC_IRR_MULTI, DBC_IRR_MULTI
  use PressureElevation, only : CalPotential
  use PressureElevationMulti, only : CalPotentialMulti
  implicit none
  private
  public :: collect_generalized_modes

  integer, parameter :: NE = 6, NN = 8      ! panels and nodes of one body
  real(8), parameter :: TOL = 1d-12         ! same arithmetic on the same inputs, so only round-off is allowed

contains

  subroutine collect_generalized_modes(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("single_override_stays_and_rotations_are_geometric", test_single_override), &
      new_unittest("single_left_side_ignores_mode", test_single_left), &
      new_unittest("single_diffraction_right_side_ignores_mode", test_single_dbc), &
      new_unittest("single_radiation_right_side_uses_mode", test_single_rbc), &
      new_unittest("single_field_point_ignores_mode", test_single_field_point), &
      new_unittest("multi_override_stays_and_rotations_are_geometric", test_multi_override), &
      new_unittest("multi_left_side_ignores_mode", test_multi_left), &
      new_unittest("multi_diffraction_right_side_ignores_mode", test_multi_dbc), &
      new_unittest("multi_radiation_right_side_uses_mode", test_multi_rbc), &
      new_unittest("multi_field_point_ignores_mode", test_multi_field_point) ]
  end subroutine collect_generalized_modes

  ! ---- Single body ------------------------------------------------------------------------------------------------

  subroutine test_single_override(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: plain(NE,6)
    integer :: dof, col
    call setup_single
    call read_single(0)
    plain = DXYZ_P(1:NE,:)
    do dof = 1, 3
      call read_single(dof)
      call check(error, all(DXYZ_P(1:NE,dof) == mode_values(1)), "the mode's values are not used as given")
      if (allocated(error)) return
      do col = 1, 6
        if (col == dof) cycle
        call check(error, all(DXYZ_P(1:NE,col) == plain(:,col)), "a mode changed another column of DXYZ_P")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_single_override

  subroutine test_single_left(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: plain(NE,NE,2,3), got(NE,NE,2,3)
    integer :: dof
    call setup_single
    call read_single(0)
    call single_left(plain)
    do dof = 1, 3
      call read_single(dof)
      call single_left(got)
      call check_close(plain(:,:,:,1), got(:,:,:,1), "BODINT_LEFT", dof, error)
      if (allocated(error)) return
      call check_close(plain(:,:,:,2), got(:,:,:,2), "BODINT_LEFT_IRR, hull rows", dof, error)
      if (allocated(error)) return
      call check_close(plain(:,:,:,3), got(:,:,:,3), "BODINT_LEFT_IRR, lid rows", dof, error)
      if (allocated(error)) return
    end do
  end subroutine test_single_left

  subroutine single_left(out)
    complex(8), intent(out) :: out(NE,NE,2,3)
    complex(8) :: tindp(4)
    integer :: iel, jel, flag
    do flag = 0, 1
      do jel = 1, NE
        do iel = 1, NE
          call BODINT_LEFT(1, iel, jel, tindp, flag)
          out(iel,jel,flag+1,1) = tindp(1)
          call BODINT_LEFT_IRR(1, iel, jel, tindp, 1, flag)
          out(iel,jel,flag+1,2) = tindp(1)
          call BODINT_LEFT_IRR(1, iel, jel, tindp, 3, flag)
          out(iel,jel,flag+1,3) = tindp(1)
        end do
      end do
    end do
  end subroutine single_left

  subroutine test_single_dbc(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: plain(NE,NE,2,3), got(NE,NE,2,3)
    integer :: dof
    call setup_single
    call read_single(0)
    call single_dbc(plain)
    do dof = 1, 3
      call read_single(dof)
      call single_dbc(got)
      call check_close(plain(:,:,:,1), got(:,:,:,1), "DBC_RIGHT", dof, error)
      if (allocated(error)) return
      call check_close(plain(:,:,:,2), got(:,:,:,2), "DBC_IRR, hull rows", dof, error)
      if (allocated(error)) return
      call check_close(plain(:,:,:,3), got(:,:,:,3), "DBC_IRR, lid rows", dof, error)
      if (allocated(error)) return
    end do
  end subroutine test_single_dbc

  subroutine single_dbc(out)
    complex(8), intent(out) :: out(NE,NE,2,3)
    complex(8) :: tinrd(4,4)
    integer :: iel, jel, flag
    do flag = 0, 1
      do jel = 1, NE
        do iel = 1, NE
          tinrd = (0d0, 0d0)
          call DBC_RIGHT(1, iel, jel, tinrd, flag)
          out(iel,jel,flag+1,1) = tinrd(1,1)
          tinrd = (0d0, 0d0)
          call DBC_IRR(1, iel, jel, tinrd, 1, flag)
          out(iel,jel,flag+1,2) = tinrd(1,1)
          tinrd = (0d0, 0d0)
          call DBC_IRR(1, iel, jel, tinrd, 3, flag)
          out(iel,jel,flag+1,3) = tinrd(1,1)
        end do
      end do
    end do
  end subroutine single_dbc

  subroutine test_single_rbc(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: plain(NE,NE,6,3), got(NE,NE,6,3), expect(NE,NE)
    real(8) :: values(NE)
    integer :: dof, col, iv, iel, jel
    call setup_single
    values = mode_values(1)
    call read_single(0)
    call single_rbc(plain)
    do dof = 1, 3
      call read_single(dof)
      call single_rbc(got)
      do iv = 1, 3
        do col = 1, 6
          if (col == dof) then
            ! The right side is the Green function times the normal, so scale the plain result by the panel's values.
            do jel = 1, NE
              do iel = 1, NE
                expect(iel,jel) = plain(iel,jel,col,iv)*values(jel)/plain_normal(jel, col)
              end do
            end do
          else
            expect = plain(:,:,col,iv)
          end if
          call check_close2(expect, got(:,:,col,iv), "RBC_RIGHT / RBC_IRR, column of the mode and the others", dof, error)
          if (allocated(error)) return
        end do
      end do
    end do
  end subroutine test_single_rbc

  subroutine single_rbc(out)
    complex(8), intent(out) :: out(NE,NE,6,3)
    complex(8) :: tinrd(4,6,4)
    integer :: iel, jel
    do jel = 1, NE
      do iel = 1, NE
        tinrd = (0d0, 0d0)
        call RBC_RIGHT(1, iel, jel, tinrd, 1)
        out(iel,jel,:,1) = tinrd(1,:,1)
        tinrd = (0d0, 0d0)
        call RBC_IRR(1, iel, jel, tinrd, 1, 1)
        out(iel,jel,:,2) = tinrd(1,:,1)
        tinrd = (0d0, 0d0)
        call RBC_IRR(1, iel, jel, tinrd, 3, 1)
        out(iel,jel,:,3) = tinrd(1,:,1)
      end do
    end do
  end subroutine single_rbc

  subroutine test_single_field_point(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: plain(2,7), got(2,7)
    integer :: dof, md
    call setup_single
    call read_single(0)
    call single_field(plain)
    do dof = 1, 3
      call read_single(dof)
      call single_field(got)
      do md = 1, 7
        if (md == dof) then
          ! Its own term carries the mode's value, so this column has to change.
          call check(error, abs(got(1,md) - plain(1,md)) > 1d-6*abs(plain(1,md)), &
                     "CalPotential: the mode's own column did not use the mode's values")
        else
          call check_close(reshape(plain(:,md), [2,1,1]), reshape(got(:,md), [2,1,1]), &
                           "CalPotential (kernel, incident term, on-body test)", dof, error)
        end if
        if (allocated(error)) return
      end do
    end do
  end subroutine test_single_field_point

  subroutine single_field(out)
    complex(8), intent(out) :: out(2,7)
    real(8) :: on_body(3), in_fluid(3)
    integer :: md
    on_body = XYZ_P(6,:)                    ! a centroid on the planar panel: the on-body test must give the free term 2 pi
    in_fluid = [6d0, 5d0, -1d0]
    do md = 1, 7
      call CalPotential(on_body, merge('Diffraction', 'Radiation  ', md == 7), md, out(1,md))
      call CalPotential(in_fluid, merge('Diffraction', 'Radiation  ', md == 7), md, out(2,md))
    end do
  end subroutine single_field

  ! ---- Multi body -------------------------------------------------------------------------------------------------

  subroutine test_multi_override(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: plain(2*NE,6)
    integer :: dof, col, body, first, other
    call setup_multi
    call read_multi(0, 0)
    plain = DXYZ_MULTI_COMB(1:2*NE,:)
    do body = 1, 2
      first = (body-1)*NE + 1
      other = (2-body)*NE + 1
      do dof = 1, 3
        call read_multi(body, dof)
        call check(error, all(DXYZ_MULTI_COMB(first:first+NE-1,dof) == mode_values(body)), &
                   "the mode's values are not used as given")
        if (allocated(error)) return
        do col = 1, 6
          if (col == dof) then
            ! The other body's panels keep their normal.
            call check(error, all(DXYZ_MULTI_COMB(other:other+NE-1,col) == plain(other:other+NE-1,col)), &
                       "a mode of one body changed the other body")
          else
            call check(error, all(DXYZ_MULTI_COMB(1:2*NE,col) == plain(:,col)), "a mode changed another column")
          end if
          if (allocated(error)) return
        end do
      end do
    end do
  end subroutine test_multi_override

  subroutine test_multi_left(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: plain(2*NE,2*NE,2,3), got(2*NE,2*NE,2,3)
    integer :: dof, body
    call setup_multi
    call read_multi(0, 0)
    call multi_left(plain)
    do body = 1, 2
      do dof = 1, 3
        call read_multi(body, dof)
        call multi_left(got)
        call check_close(plain(:,:,:,1), got(:,:,:,1), "BODINT_MULTI_LEFT", dof, error)
        if (allocated(error)) return
        call check_close(plain(:,:,:,2), got(:,:,:,2), "BODINT_LEFT_IRR_MULTI, hull rows", dof, error)
        if (allocated(error)) return
        call check_close(plain(:,:,:,3), got(:,:,:,3), "BODINT_LEFT_IRR_MULTI, lid rows", dof, error)
        if (allocated(error)) return
      end do
    end do
  end subroutine test_multi_left

  subroutine multi_left(out)
    complex(8), intent(out) :: out(2*NE,2*NE,2,3)
    complex(8) :: tindp(4)
    integer :: iel, jel, flag
    do flag = 0, 1
      do jel = 1, 2*NE
        do iel = 1, 2*NE
          call BODINT_MULTI_LEFT(1, iel, jel, tindp, flag)
          out(iel,jel,flag+1,1) = tindp(1)
          call BODINT_LEFT_IRR_MULTI(1, iel, jel, tindp, 1, flag)
          out(iel,jel,flag+1,2) = tindp(1)
          call BODINT_LEFT_IRR_MULTI(1, iel, jel, tindp, 3, flag)
          out(iel,jel,flag+1,3) = tindp(1)
        end do
      end do
    end do
  end subroutine multi_left

  subroutine test_multi_dbc(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: plain(2*NE,2*NE,2,3), got(2*NE,2*NE,2,3)
    integer :: dof, body
    call setup_multi
    call read_multi(0, 0)
    call multi_dbc(plain)
    do body = 1, 2
      do dof = 1, 3
        call read_multi(body, dof)
        call multi_dbc(got)
        call check_close(plain(:,:,:,1), got(:,:,:,1), "DBC_MULTI_RIGHT", dof, error)
        if (allocated(error)) return
        call check_close(plain(:,:,:,2), got(:,:,:,2), "DBC_IRR_MULTI, hull rows", dof, error)
        if (allocated(error)) return
        call check_close(plain(:,:,:,3), got(:,:,:,3), "DBC_IRR_MULTI, lid rows", dof, error)
        if (allocated(error)) return
      end do
    end do
  end subroutine test_multi_dbc

  subroutine multi_dbc(out)
    complex(8), intent(out) :: out(2*NE,2*NE,2,3)
    complex(8) :: tinrd(4,4)
    integer :: iel, jel, flag
    do flag = 0, 1
      do jel = 1, 2*NE
        do iel = 1, 2*NE
          tinrd = (0d0, 0d0)
          call DBC_MULTI_RIGHT(1, iel, jel, tinrd, flag)
          out(iel,jel,flag+1,1) = tinrd(1,1)
          tinrd = (0d0, 0d0)
          call DBC_IRR_MULTI(1, iel, jel, tinrd, 1, flag)
          out(iel,jel,flag+1,2) = tinrd(1,1)
          tinrd = (0d0, 0d0)
          call DBC_IRR_MULTI(1, iel, jel, tinrd, 3, flag)
          out(iel,jel,flag+1,3) = tinrd(1,1)
        end do
      end do
    end do
  end subroutine multi_dbc

  subroutine test_multi_rbc(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: plain(2*NE,2*NE,6,3), got(2*NE,2*NE,6,3), expect(2*NE,2*NE)
    real(8) :: values(NE)
    integer :: dof, col, iv, iel, jel, jloc, body
    call setup_multi
    call read_multi(0, 0)
    call multi_rbc(plain)
    do body = 1, 2
      values = mode_values(body)
      do dof = 1, 3
        call read_multi(body, dof)
        call multi_rbc(got)
        do iv = 1, 3
          do col = 1, 6
            expect = plain(:,:,col,iv)
            if (col == dof) then
              do jel = (body-1)*NE + 1, body*NE
                jloc = jel - (body-1)*NE
                do iel = 1, 2*NE
                  expect(iel,jel) = plain(iel,jel,col,iv)*values(jloc)/plain_normal(jloc, col)
                end do
              end do
            end if
            call check_close2(expect, got(:,:,col,iv), "RBC_MULTI_RIGHT / RBC_IRR_MULTI", dof, error)
            if (allocated(error)) return
          end do
        end do
      end do
    end do
  end subroutine test_multi_rbc

  subroutine multi_rbc(out)
    complex(8), intent(out) :: out(2*NE,2*NE,6,3)
    complex(8) :: tinrd(4,6,4)
    integer :: iel, jel
    do jel = 1, 2*NE
      do iel = 1, 2*NE
        tinrd = (0d0, 0d0)
        call RBC_MULTI_RIGHT(1, iel, jel, tinrd, 1)
        out(iel,jel,:,1) = tinrd(1,:,1)
        tinrd = (0d0, 0d0)
        call RBC_IRR_MULTI(1, iel, jel, tinrd, 1, 1)
        out(iel,jel,:,2) = tinrd(1,:,1)
        tinrd = (0d0, 0d0)
        call RBC_IRR_MULTI(1, iel, jel, tinrd, 3, 1)
        out(iel,jel,:,3) = tinrd(1,:,1)
      end do
    end do
  end subroutine multi_rbc

  subroutine test_multi_field_point(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: plain(2,13), got(2,13)
    integer :: dof, md, body, own
    call setup_multi
    call read_multi(0, 0)
    call multi_field(plain)
    do body = 1, 2
      do dof = 1, 3
        call read_multi(body, dof)
        call multi_field(got)
        own = 6*(body-1) + dof
        do md = 1, 13
          if (md == own) then
            call check(error, abs(got(1,md) - plain(1,md)) > 1d-6*abs(plain(1,md)), &
                       "CalPotentialMulti: the mode's own column did not use the mode's values")
          else
            call check_close(reshape(plain(:,md), [2,1,1]), reshape(got(:,md), [2,1,1]), &
                             "CalPotentialMulti (kernel, incident term, on-body test)", dof, error)
          end if
          if (allocated(error)) return
        end do
      end do
    end do
  end subroutine test_multi_field_point

  subroutine multi_field(out)
    complex(8), intent(out) :: out(2,13)
    real(8) :: on_body(3), in_fluid(3)
    integer :: factor(2*NE), md, body, dof, first
    on_body = XYZ_GLOBAL_MULTI_COMB_P(2*NE,:)      ! the centroid of the planar panel of body 2
    in_fluid = [6d0, 5d0, -1d0]
    factor = 0
    do body = 1, 2
      first = (body-1)*NE + 1
      factor = 0
      factor(first:first+NE-1) = 1
      do dof = 1, 6
        md = 6*(body-1) + dof
        call CalPotentialMulti(on_body, 'Radiation', dof, out(1,md), 2, body, factor)
        call CalPotentialMulti(in_fluid, 'Radiation', dof, out(2,md), 2, body, factor)
      end do
    end do
    factor = 1
    call CalPotentialMulti(on_body, 'Diffraction', 13, out(1,13), 2, 1, factor)
    call CalPotentialMulti(in_fluid, 'Diffraction', 13, out(2,13), 2, 1, factor)
  end subroutine multi_field

  ! ---- Set-up -------------------------------------------------------------------------------------------------------

  !> A skew box, so the normal's components are all non-zero on most panels and rotational columns are non-trivial.
  subroutine box_mesh(nodes, conn)
    real(8), intent(out) :: nodes(NN,3)
    integer, intent(out) :: conn(NE,4)
    nodes(1,:) = [-1.0d0, -1.1d0, -2.1d0]
    nodes(2,:) = [ 1.2d0, -0.9d0, -2.0d0]
    nodes(3,:) = [ 1.1d0,  1.3d0, -1.8d0]
    nodes(4,:) = [-0.8d0,  1.0d0, -2.2d0]
    nodes(5,:) = [-1.3d0, -0.7d0, 0d0]
    nodes(6,:) = [ 0.9d0, -1.2d0, 0d0]
    nodes(7,:) = [ 1.4d0,  0.8d0, 0d0]
    nodes(8,:) = [-1.1d0,  1.2d0, 0d0]
    ! Panel 6 lies in one tilted plane, so that a point on it passes the on-body test of the field-point code.
    nodes(5:8,3) = -0.5d0 + 0.05d0*nodes(5:8,1) - 0.03d0*nodes(5:8,2)
    conn(1,:) = [1, 2, 3, 4]
    conn(2,:) = [1, 5, 6, 2]
    conn(3,:) = [2, 6, 7, 3]
    conn(4,:) = [3, 7, 8, 4]
    conn(5,:) = [4, 8, 5, 1]
    conn(6,:) = [5, 8, 7, 6]
  end subroutine box_mesh

  !> The values of a mode: unrelated to any normal, up to a size like example 3's.
  function mode_values(body) result(values)
    integer, intent(in) :: body
    real(8) :: values(NE)
    integer :: i
    do i = 1, NE
      values(i) = (-1d0)**i*(3d0 + 2.5d0*i + body)
    end do
  end function mode_values

  !> The folder holding the mode of (body, dof), created once. Body 0 is the folder with no mode file.
  function mode_dir(body, dof) result(dir)
    integer, intent(in) :: body, dof
    character(:), allocatable :: dir
    character(len=16) :: tag
    logical :: exists, success
    integer :: unit, i
    real(8) :: values(NE)
    write(tag, '(a,i0,a,i0)') 'b', body, '_d', dof
    dir = get_running_exe_path()//'gen_modes_unit_'//trim(tag)
    inquire(file=dir//'/.ready', exist=exists)
    if (exists) return
    call CreateDirectory(dir, success)
    if (.not. success) error stop 'test_generalized_modes: cannot create the folder for mode files'
    if (body > 0) then
      values = mode_values(body)
      write(tag, '(i0,a,i0)') body, '_', dof
      open(newunit=unit, file=dir//'/gen_mod_'//trim(tag)//'.txt', status='replace', action='write')
      do i = 1, NE
        write(unit, '(i6,1x,es24.16)') i, values(i)
      end do
      close(unit)
    end if
    open(newunit=unit, file=dir//'/.ready', status='replace', action='write')
    close(unit)
  end function mode_dir

  subroutine common_state
    NSYS = 1
    ISX = 0
    ISY = 0
    ISOL = 1
    NTHREAD = 1
    H = -1d0                 ! deep water
    W1 = 1.1d0
    V = W1**2/G
    WK = V
    TP = 2d0*PI/W1
    AMP = 1d0
    BETA = 0.7d0             ! a heading with both x and y parts of the incident gradient
    XW = 0d0
  end subroutine common_state

  subroutine setup_single
    integer :: i, j, k
    real(8) :: nodes(NN,3)
    integer :: conn(NE,4)
    call common_state
    IRSP = 0
    NBODY = 1
    NELEM = NE
    NTND = NN
    INELEM = NE
    XR = [0.1d0, -0.2d0, -1d0]
    if (.not. allocated(XYZ)) then
      allocate(XYZ(NN,3), DS(NE), PNSZ(NE), NCN(NE), NCON(NE,4), XYZ_P(NE,3), DXYZ_P(NE,6), iDXYZ_P(NE,6))
      allocate(CGRN(NE,NE,1,4), RKBN(NE,NE,1,4), DGRN(NE,NE,1,4), PKBN(NE,NE,1,4), MXPOT(NE,7,1))
    end if
    call box_mesh(nodes, conn)
    XYZ = nodes
    NCN = 4
    NCON = conn
    iDXYZ_P = 0d0
    do k = 1, 4
      do j = 1, NE
        do i = 1, NE
          RKBN(i,j,1,k) = real(sin(1d0*i + 2d0*j + 3d0*k), kind(RKBN))
          PKBN(i,j,1,k) = real(cos(2d0*i + 1d0*j + 5d0*k), kind(PKBN))
          CGRN(i,j,1,k) = cmplx(cos(1d0*i - 1d0*j + k), sin(3d0*i + j - k), kind(CGRN))
          DGRN(i,j,1,k) = cmplx(sin(2d0*i - 3d0*j + k), cos(i + 2d0*j + k), kind(DGRN))
        end do
      end do
    end do
    do k = 1, 7
      do i = 1, NE
        MXPOT(i,k,1) = dcmplx(sin(1d0*i*k), cos(1d0*i + k))
      end do
    end do
  end subroutine setup_single

  subroutine read_single(dof)
    integer, intent(in) :: dof
    call CalNormals(0, mode_dir(merge(1, 0, dof > 0), dof))
  end subroutine read_single

  subroutine setup_multi
    integer :: i, j, k, b
    real(8) :: nodes(NN,3)
    integer :: conn(NE,4)
    call common_state
    IRSP = 0
    NBODY = 2
    TNELEM = 2*NE
    TNTND = 2*NN
    NELEM_PE = 2*NE
    if (.not. allocated(XYZ_GLOBAL_MULTI)) then
      allocate(NELEM_MULTI(2), NTND_MULTI(2), XR_MULTI(2,3))
      allocate(XYZ_GLOBAL_MULTI(2,NN,3), DS_MULTI(2,NE), PNSZ_MULTI(2,NE), NCN_MULTI(2,NE), NCON_MULTI(2,NE,4))
      allocate(XYZ_MULTI_P(2,NE,3), DXYZ_MULTI_P(2,NE,6))
      allocate(XYZ_GLOBAL_MULTI_COMB(2*NN,3), XYZ_GLOBAL_MULTI_COMB_P(2*NE,3), PNSZ_MULTI_COMB(2*NE), &
               NCN_MULTI_COMB(2*NE), NCON_MULTI_COMB(2*NE,4), DS_MULTI_COMB(2*NE), DXYZ_MULTI_COMB(2*NE,6), iDXYZ_MULTI_COMB(2*NE,6))
      allocate(CGRN_MULTI_COMB(2*NE,2*NE,1,4), RKBN_MULTI_COMB(2*NE,2*NE,1,4), DGRN_MULTI_COMB(2*NE,2*NE,1,4), &
               PKBN_MULTI_COMB(2*NE,2*NE,1,4), MXPOT_MULTI_COMB(2*NE,13,1))
    end if
    NELEM_MULTI = NE
    NTND_MULTI = NN
    call box_mesh(nodes, conn)
    NCN_MULTI = 4
    do b = 1, 2
      XYZ_GLOBAL_MULTI(b,:,:) = nodes
      XYZ_GLOBAL_MULTI(b,:,2) = nodes(:,2) + 10d0*(b-1)
      NCON_MULTI(b,:,:) = conn
      XR_MULTI(b,:) = [0.1d0, 10d0*(b-1), -1d0]
    end do
    iDXYZ_MULTI_COMB = 0d0
    do k = 1, 4
      do j = 1, 2*NE
        do i = 1, 2*NE
          RKBN_MULTI_COMB(i,j,1,k) = real(sin(1d0*i + 2d0*j + 3d0*k), kind(RKBN_MULTI_COMB))
          PKBN_MULTI_COMB(i,j,1,k) = real(cos(2d0*i + 1d0*j + 5d0*k), kind(PKBN_MULTI_COMB))
          CGRN_MULTI_COMB(i,j,1,k) = cmplx(cos(1d0*i - 1d0*j + k), sin(3d0*i + j - k), kind(CGRN_MULTI_COMB))
          DGRN_MULTI_COMB(i,j,1,k) = cmplx(sin(2d0*i - 3d0*j + k), cos(i + 2d0*j + k), kind(DGRN_MULTI_COMB))
        end do
      end do
    end do
    do k = 1, 13
      do i = 1, 2*NE
        MXPOT_MULTI_COMB(i,k,1) = dcmplx(sin(1d0*i*k), cos(1d0*i + k))
      end do
    end do
  end subroutine setup_multi

  subroutine read_multi(body, dof)
    integer, intent(in) :: body, dof
    integer :: b
    do b = 1, 2
      call CalNormalsMulti(b, 0, mode_dir(body, dof))
    end do
    call BuildCombinedArrays
  end subroutine read_multi

  !> The normal a panel has with no mode: the same formula as CalTransNormals, computed here from the nodes.
  function plain_normal(iel, col) result(value)
    integer, intent(in) :: iel, col
    real(8) :: value, nodes(NN,3), v13(3), v24(3), un(3)
    integer :: conn(NE,4)
    call box_mesh(nodes, conn)
    v13 = nodes(conn(iel,3),:) - nodes(conn(iel,1),:)
    v24 = nodes(conn(iel,4),:) - nodes(conn(iel,2),:)
    un = [v24(2)*v13(3) - v24(3)*v13(2), v24(3)*v13(1) - v24(1)*v13(3), v24(1)*v13(2) - v24(2)*v13(1)]
    un = un/sqrt(sum(un**2))
    value = un(col)
  end function plain_normal

  ! ---- Checks -------------------------------------------------------------------------------------------------------

  subroutine check_close(ref, got, what, dof, error)
    complex(8), intent(in) :: ref(:,:,:), got(:,:,:)
    character(len=*), intent(in) :: what
    integer, intent(in) :: dof
    type(error_type), allocatable, intent(out) :: error
    real(8) :: rel
    character(len=64) :: info
    rel = maxval(abs(got - ref))/max(maxval(abs(ref)), 1d-30)
    write(info, '(a,i0,a,es10.3)') ' with a mode in DOF ', dof, ': largest relative change ', rel
    call check(error, rel <= TOL, what//" reads the mode's values as the panel's direction,"//trim(info))
  end subroutine check_close

  subroutine check_close2(ref, got, what, dof, error)
    complex(8), intent(in) :: ref(:,:), got(:,:)
    character(len=*), intent(in) :: what
    integer, intent(in) :: dof
    type(error_type), allocatable, intent(out) :: error
    real(8) :: rel
    character(len=64) :: info
    rel = maxval(abs(got - ref))/max(maxval(abs(ref)), 1d-30)
    write(info, '(a,i0,a,es10.3)') ' with a mode in DOF ', dof, ': largest relative difference ', rel
    ! The columns are built from single-precision kernels and divided by a normal, so round-off is larger here.
    call check(error, rel <= 1d-9, what//" does not use the mode's values as given,"//trim(info))
  end subroutine check_close2

end module test_generalized_modes
