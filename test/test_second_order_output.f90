!> Unit tests for the second-order output writers.
module test_second_order_output
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use Const_mod, only : G, RHO
  use SecondOrderOutput, only : write_drift_row, write_terms_row, write_field_row, write_direct_sum_terms_header, &
                                write_qtf_terms_row
  implicit none
  private

  public :: collect_second_order_output

contains

subroutine collect_second_order_output(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)
  testsuite = [ &
    new_unittest("drift_row_force", drift_row_force), &
    new_unittest("drift_row_moment", drift_row_moment), &
    new_unittest("terms_row", terms_row), &
    new_unittest("field_row", field_row), &
    new_unittest("field_direct_columns", field_direct_columns), &
    new_unittest("direct_sum_terms_header", direct_sum_terms_header) &
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

!> A field row is conjugated, scaled by L/(rho g) for a pressure and by L for an elevation, with Mod, Pha, Re and Im of the
!> sum and of its two parts, and K, IPT and TYPE as given.
subroutine field_row(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: path = "second_order_output_field.tmp"
  complex(8), parameter :: quad = (3d0, 4d0), pot = (-1d0, 2d0)
  real(8) :: a(4), b(8)
  integer :: unit, i(3), pass

  open (newunit=unit, file=path, status='replace', action='readwrite')
  call write_field_row(unit, 0.5d0, 0.7d0, 10d0, 20d0, 2, 3, .false., quad*RHO*G, pot*RHO*G, 2d0)
  call write_field_row(unit, 0.5d0, 0.7d0, 10d0, 20d0, 1, 4, .true., quad, pot, 2d0)
  rewind (unit)
  do pass = 1, 2
    read (unit, '(4ES14.6,3I6,8ES14.6)') a, i, b
    if (pass == 1) then
      call check(error, all(i == [2, 3, 1]), "K, IPT and TYPE of the pressure row")
    else
      call check(error, all(i == [1, 4, 2]), "K, IPT and TYPE of the elevation row")
    end if
    if (allocated(error)) exit
    call check(error, all(abs(a - [0.5d0, 0.7d0, 10d0, 20d0]) < 1d-6), "PER and BETA")
    if (allocated(error)) exit
    ! L = 2: the sum is conj(quad+pot) * 2 = (4, -12), the parts (6, -8) and (-2, -4)
    call check(error, abs(b(1) - sqrt(160d0)) < 1d-5, "Mod")
    if (allocated(error)) exit
    call check(error, abs(b(2) - atan2(-12d0, 4d0)*180d0/acos(-1d0)) < 1d-3, "Pha")
    if (allocated(error)) exit
    call check(error, abs(b(3) - 4d0) < 1d-5 .and. abs(b(4) + 12d0) < 1d-5, "Re and Im are not the conjugate of the sum")
    if (allocated(error)) exit
    call check(error, abs(b(5) - 6d0) < 1d-5 .and. abs(b(6) + 8d0) < 1d-5, "Re_quad and Im_quad")
    if (allocated(error)) exit
    call check(error, abs(b(7) + 2d0) < 1d-5 .and. abs(b(8) + 4d0) < 1d-5, "Re_inc and Im_inc")
    if (allocated(error)) exit
  end do
  close (unit, status='delete')
end subroutine field_row

!> With direct, a field row gets Re and Im of the direct part and has_direct after the 15 columns of the row without it,
!> conjugated and scaled as the other parts (pressure by L/(rho g), elevation by L), and the first 15 columns are those of the
!> row without direct, to the byte. has_direct = .false. writes a zero direct part and 0.
subroutine field_direct_columns(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: path = "second_order_output_field_direct.tmp"
  complex(8), parameter :: quad = (3d0, 4d0), pot = (-1d0, 2d0), dir = (0.5d0, -1.5d0)
  character(len=400) :: plain, with, without
  real(8) :: b(2)
  integer :: unit, ih

  open (newunit=unit, file=path, status='replace', action='readwrite')
  call write_field_row(unit, 0.5d0, 0.7d0, 10d0, 20d0, 1, 4, .true., quad, pot, 2d0)
  call write_field_row(unit, 0.5d0, 0.7d0, 10d0, 20d0, 1, 4, .true., quad, pot, 2d0, dir)
  call write_field_row(unit, 0.5d0, 0.7d0, 10d0, 20d0, 1, 4, .true., quad, pot, 2d0, dir, .false.)
  call write_field_row(unit, 0.5d0, 0.7d0, 10d0, 20d0, 2, 3, .false., quad*RHO*G, pot*RHO*G, 2d0, dir*RHO*G)
  rewind (unit)
  read (unit, '(a)') plain
  read (unit, '(a)') with
  read (unit, '(a)') without
  call check(error, with(1:len_trim(plain)) == trim(plain), "The row with direct does not start with the row without it")
  if (allocated(error)) return
  call check(error, len_trim(with) == len_trim(plain) + 2*14 + 6, "The row with direct has not two reals and an integer more")
  if (allocated(error)) return
  ! L = 2: conj(dir) * 2 = (1, 3)
  read (with(len_trim(plain) + 1:), '(2ES14.6,I6)') b, ih
  call check(error, abs(b(1) - 1d0) < 1d-5 .and. abs(b(2) - 3d0) < 1d-5 .and. ih == 1, "Re_dir, Im_dir, has_direct")
  if (allocated(error)) return
  read (without(len_trim(plain) + 1:), '(2ES14.6,I6)') b, ih
  call check(error, abs(b(1)) < 1d-30 .and. abs(b(2)) < 1d-30 .and. ih == 0, "A row without a direct part is not zero and 0")
  if (allocated(error)) return
  ! The pressure row: dir given as rho g times the unit value, scaled as quad and pot by L/(rho g)
  read (unit, '(a)') with
  read (with(len_trim(plain) + 1:), '(2ES14.6,I6)') b, ih
  call check(error, abs(b(1) - 1d0) < 1d-5 .and. abs(b(2) - 3d0) < 1d-5 .and. ih == 1, "Pressure row: Re_dir, Im_dir")
  close (unit, status='delete')
end subroutine field_direct_columns

!> QTF_DirectTermsSum.txt (Phase 7c-2): four comment lines; with full the column line names D_I2, D_D, D_FS, P_BODY and P_FS
!> and a row (keys and 5 values, as QTFDriver writes it) has 15 columns; with body it names D_I2, D_D and P_BODY, not D_FS
!> or P_FS, and a row has 11 columns. (The brief of Task 4 said 13 and 9, which does not fit its own list of 5 keys and
!> 3 or 5 complex values; ruling of Task 4, see its report.)
subroutine direct_sum_terms_header(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: path = "second_order_output_sum_terms.tmp"
  character(len=600) :: line(5)
  integer :: unit, pass, i, n
  logical :: full

  do pass = 1, 2
    full = pass == 1
    open (newunit=unit, file=path, status='replace', action='readwrite')
    call write_direct_sum_terms_header(unit, full)
    if (full) then
      call write_qtf_terms_row(unit, 0.5d0, 0.7d0, 0d0, 30d0, 3, [(1d0, 2d0), (3d0, 4d0), (5d0, 6d0), (7d0, 8d0), &
                               (9d0, 1d0)], 1d0, .false.)
    else
      call write_qtf_terms_row(unit, 0.5d0, 0.7d0, 0d0, 30d0, 3, [(1d0, 2d0), (3d0, 4d0), (5d0, 6d0)], 1d0, .false.)
    end if
    rewind (unit)
    do i = 1, 5
      read (unit, '(a)') line(i)
    end do
    close (unit, status='delete')
    call check(error, all(line(1:4)(1:1) == "#") .and. line(5)(1:1) /= "#", "four comment lines, then the row")
    if (allocated(error)) return
    if (full) then
      call check(error, index(line(2), "D_I2") > 0 .and. index(line(2), "D_D") > 0 .and. index(line(2), "D_FS") > 0 .and. &
                 index(line(2), "P_BODY") > 0 .and. index(line(2), "P_FS") > 0, "full: the columns are not all named")
    else
      call check(error, index(line(2), "D_I2") > 0 .and. index(line(2), "D_D") > 0 .and. index(line(2), "P_BODY") > 0 .and. &
                 index(line(2), "D_FS") == 0 .and. index(line(2), "P_FS") == 0, "body: the columns are not named as they are")
    end if
    if (allocated(error)) return
    n = count_fields(line(5))
    call check(error, n, merge(15, 11, full), "columns of a row")
    if (allocated(error)) return
  end do
end subroutine direct_sum_terms_header

integer function count_fields(s)
  character(len=*), intent(in) :: s
  integer :: i
  logical :: blank

  count_fields = 0
  blank = .true.
  do i = 1, len_trim(s)
    if (s(i:i) /= " " .and. blank) count_fields = count_fields + 1
    blank = s(i:i) == " "
  end do
end function count_fields

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
