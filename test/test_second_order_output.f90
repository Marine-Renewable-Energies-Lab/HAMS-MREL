!> Unit tests for the second-order output writers.
module test_second_order_output
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use Const_mod, only : G, RHO
  use SecondOrderOutput, only : write_drift_row, write_terms_row
  implicit none
  private

  public :: collect_second_order_output

contains

subroutine collect_second_order_output(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)
  testsuite = [ &
    new_unittest("drift_row_force", drift_row_force), &
    new_unittest("drift_row_moment", drift_row_moment), &
    new_unittest("terms_row", terms_row) &
  ]
end subroutine collect_second_order_output

!> A force is conjugated once and divided by rho g refl, and the row has the documented layout.
subroutine drift_row_force(error)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: row(7)
  integer :: index

  call write_and_read("second_order_output_force.tmp", (3d0, 4d0)*RHO*G*2d0, .false., row, index)
  call check(error, index == 2)
  if (allocated(error)) return
  call check(error, all(abs(row(1:3) - [0.5d0, 10d0, 20d0]) < 1d-6), "PER, BETA1 and BETA2 are wrong")
  if (allocated(error)) return
  call check(error, abs(row(4) - 5d0) < 1d-5, "Mod is wrong")
  if (allocated(error)) return
  call check(error, abs(row(5) - atan2(-4d0, 3d0)*180d0/acos(-1d0)) < 1d-3, "Pha is wrong")
  if (allocated(error)) return
  call check(error, abs(row(6) - 3d0) < 1d-5 .and. abs(row(7) + 4d0) < 1d-5, "Re and Im are not the conjugate")
end subroutine drift_row_force

!> A moment is divided by rho g refl^2.
subroutine drift_row_moment(error)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: row(7)
  integer :: index

  call write_and_read("second_order_output_moment.tmp", (3d0, 4d0)*RHO*G*4d0, .true., row, index)
  call check(error, abs(row(6) - 3d0) < 1d-5 .and. abs(row(7) + 4d0) < 1d-5, "The moment scaling is wrong")
end subroutine drift_row_moment

!> A term is conjugated and scaled as a load, and the row has the documented layout.
subroutine terms_row(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: path = "second_order_output_terms.tmp"
  real(8) :: a(3), b(2)
  integer :: unit, i(3)

  open (newunit=unit, file=path, status='replace', action='readwrite')
  call write_terms_row(unit, 0.5d0, 10d0, 20d0, 2, 6, 5, (3d0, 4d0)*RHO*G*4d0, 2d0, .true.)
  call write_terms_row(unit, 0.5d0, 10d0, 20d0, 2, 1, 5, (3d0, 4d0)*RHO*G*2d0, 2d0, .false.)
  rewind (unit)
  read (unit, '(3ES14.6,3I6,2ES14.6)') a, i, b
  call check(error, all(abs(a - [0.5d0, 10d0, 20d0]) < 1d-6) .and. all(i == [2, 6, 5]), "The layout is wrong")
  if (allocated(error)) return
  call check(error, abs(b(1) - 3d0) < 1d-5 .and. abs(b(2) + 4d0) < 1d-5, "The moment is not conjugated and scaled")
  if (allocated(error)) return
  read (unit, '(3ES14.6,3I6,2ES14.6)') a, i, b
  close (unit, status='delete')
  call check(error, abs(b(1) - 3d0) < 1d-5 .and. abs(b(2) + 4d0) < 1d-5, "The force is not conjugated and scaled")
end subroutine terms_row

subroutine write_and_read(path, value, moment, row, index)
  character(len=*), intent(in) :: path
  complex(8), intent(in) :: value
  logical, intent(in) :: moment
  real(8), intent(out) :: row(7)
  integer, intent(out) :: index
  integer :: unit

  open (newunit=unit, file=path, status='replace', action='readwrite')
  call write_drift_row(unit, 0.5d0, 10d0, 20d0, 2, value, 2d0, moment)
  rewind (unit)
  read (unit, '(3ES14.6,I6,4ES14.6)') row(1:3), index, row(4:7)
  close (unit, status='delete')
end subroutine write_and_read

end module test_second_order_output
