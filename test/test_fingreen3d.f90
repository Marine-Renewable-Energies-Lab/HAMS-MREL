!> Unit tests for the finite-depth Green function (FinGreen3D.f90).
!> Both tests check continuity, so they need no reference data:
!>  1. G at the self point (R=0, z=z') must equal G at R=1e-4 m (LINTON SR->0 limit).
!>  2. G must be continuous across R/h = 0.0005, where LINTON hands over to PIDCOCK.
module test_fingreen3d
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use FinGreen3D_Open, only : FINGREEN3D_R
  implicit none
  private
  public :: collect_fingreen3d

  integer, parameter :: NK = 200
  real(8), parameter :: GRAV = 9.807d0
  real(8), parameter :: OMEGAS(3) = [0.3d0, 0.8d0, 1.5d0]
  real(8), parameter :: DEPTHS(2) = [70d0, 150d0]

contains

  subroutine collect_fingreen3d(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("selfterm_continuity", test_selfterm_continuity), &
      new_unittest("linton_pidcock_continuity", test_linton_pidcock_continuity) ]
  end subroutine collect_fingreen3d

  subroutine test_selfterm_continuity(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: WVN(NK), V, W, H
    complex(8) :: g0(3), g1(3)
    integer :: i, j
    external :: DISPERSION
    do j = 1, size(DEPTHS)
      do i = 1, size(OMEGAS)
        H = DEPTHS(j); W = OMEGAS(i); V = W*W/GRAV
        call DISPERSION(WVN, NK, W, H)
        call FINGREEN3D_R(0d0,  -20d0, -20d0, V, WVN, NK, H, g0, 1)
        call FINGREEN3D_R(1d-4, -20d0, -20d0, V, WVN, NK, H, g1, 1)
        call check(error, abs(g0(1) - g1(1)) < 1d-6*abs(g1(1)), &
                   "LINTON self-term (R=0) not continuous with R=1e-4")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_selfterm_continuity

  subroutine test_linton_pidcock_continuity(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: WVN(NK), V, W, H
    complex(8) :: gl(3), gp(3)
    integer :: i, j
    external :: DISPERSION
    do j = 1, size(DEPTHS)
      do i = 1, size(OMEGAS)
        H = DEPTHS(j); W = OMEGAS(i); V = W*W/GRAV
        call DISPERSION(WVN, NK, W, H)
        call FINGREEN3D_R(0.99999d0*5d-4*H, -20d0, -24d0, V, WVN, NK, H, gl, 1)
        call FINGREEN3D_R(1.00001d0*5d-4*H, -20d0, -24d0, V, WVN, NK, H, gp, 1)
        call check(error, abs(gl(1) - gp(1)) < 1d-2*abs(gp(1)), &
                   "G jumps across R/h = 0.0005 (LINTON -> PIDCOCK)")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_linton_pidcock_continuity

end module test_fingreen3d
