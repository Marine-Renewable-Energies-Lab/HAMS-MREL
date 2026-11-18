!> Powers of the reference length L in the WAMIT-format outputs (WAMIT User Manual, section 4, "Definitions of
!> nondimensional quantities"): added mass and damping L^k with k = 3, 4, 5 for translation-translation, mixed and
!> rotation-rotation pairs; exciting force L^m with m = 2 (translation) or 3 (rotation); hydrostatic restoring L^k with
!> k = 2, 3, 4 in the same way; motions A/L^n with n = 0 (translation) or 1 (rotation). Multi-body modes 7, 8, ... follow the
!> rigid-body mode they stand for (mode 7 is surge of body 2).
module test_wamit_scaling

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use Const_mod, only : LEXP_AB, LEXP_X, LEXP_C, LEXP_MOTION

  implicit none
  private
  public :: collect_wamit_scaling

contains

  subroutine collect_wamit_scaling(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("added_mass_and_damping", test_added_mass_and_damping), &
      new_unittest("exciting_force", test_exciting_force), &
      new_unittest("hydrostatics", test_hydrostatics), &
      new_unittest("motions", test_motions), &
      new_unittest("multi_body_modes", test_multi_body_modes) &
      ]
  end subroutine collect_wamit_scaling

  subroutine test_added_mass_and_damping(error)
    type(error_type), allocatable, intent(out) :: error
    call check(error, LEXP_AB(1, 1) == 3 .and. LEXP_AB(3, 2) == 3, "translation-translation: L^3")
    if (allocated(error)) return
    call check(error, LEXP_AB(1, 5) == 4 .and. LEXP_AB(6, 2) == 4, "mixed: L^4")
    if (allocated(error)) return
    call check(error, LEXP_AB(5, 5) == 5 .and. LEXP_AB(4, 6) == 5, "rotation-rotation: L^5")
  end subroutine test_added_mass_and_damping

  subroutine test_exciting_force(error)
    type(error_type), allocatable, intent(out) :: error
    call check(error, LEXP_X(1) == 2 .and. LEXP_X(3) == 2, "translation: L^2")
    if (allocated(error)) return
    call check(error, LEXP_X(4) == 3 .and. LEXP_X(6) == 3, "rotation: L^3")
  end subroutine test_exciting_force

  subroutine test_hydrostatics(error)
    type(error_type), allocatable, intent(out) :: error
    call check(error, LEXP_C(3, 3) == 2, "C33: L^2")
    if (allocated(error)) return
    call check(error, LEXP_C(3, 4) == 3 .and. LEXP_C(5, 3) == 3, "C34, C53: L^3")
    if (allocated(error)) return
    call check(error, LEXP_C(4, 4) == 4 .and. LEXP_C(5, 6) == 4, "C44, C56: L^4")
  end subroutine test_hydrostatics

  subroutine test_motions(error)
    type(error_type), allocatable, intent(out) :: error
    call check(error, LEXP_MOTION(1) == 0 .and. LEXP_MOTION(3) == 0, "translation: A/L^0")
    if (allocated(error)) return
    call check(error, LEXP_MOTION(4) == 1 .and. LEXP_MOTION(6) == 1, "rotation: A/L")
  end subroutine test_motions

  subroutine test_multi_body_modes(error)
    type(error_type), allocatable, intent(out) :: error
    call check(error, LEXP_AB(7, 11) == 4 .and. LEXP_AB(10, 12) == 5 .and. LEXP_AB(1, 9) == 3, "A and B of modes 7-12")
    if (allocated(error)) return
    call check(error, LEXP_X(9) == 2 .and. LEXP_X(11) == 3 .and. LEXP_MOTION(10) == 1 .and. LEXP_MOTION(7) == 0, &
               "X and motions of modes 7-12")
    if (allocated(error)) return
    call check(error, LEXP_C(9, 9) == 2 .and. LEXP_C(9, 11) == 3 .and. LEXP_C(11, 11) == 4, "C of modes 7-12")
  end subroutine test_multi_body_modes

end module test_wamit_scaling
