!> Unit tests for CreateDirectory, which runs the platform's mkdir through execute_command_line.
module test_create_directory
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use IO, only : CreateDirectory
  implicit none
  private

  public :: collect_create_directory

  ! test-drive runs tests in parallel under OpenMP, so each test needs its own directory

contains

subroutine collect_create_directory(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("creates_nested_directories", test_creates_nested), &
    new_unittest("handles_spaces_and_double_separators", test_spaces_and_double_separators), &
    new_unittest("keeps_existing_directory_contents", test_keeps_existing), &
    new_unittest("reports_failed_mkdir", test_reports_failure) &
    ]

end subroutine collect_create_directory

subroutine test_creates_nested(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_dir_nested"
  logical :: success, leftover

  call remove_directory(root)
  call CreateDirectory(root//"/Output/Hams_format", success)
  call check(error, success, "CreateDirectory reported a failure")
  if (allocated(error)) return
  call check(error, accepts_files(root//"/Output/Hams_format"), "nested directory was not created")
  if (allocated(error)) return
  inquire(file=root//"/Output/Hams_format/.hams_dir_probe", exist=leftover)
  call check(error, .not. leftover, "probe file was left behind")
  if (allocated(error)) return
  call check(error, .not. accepts_files("-p"), "a stray -p directory was created")
  call remove_directory(root)
end subroutine test_creates_nested

subroutine test_spaces_and_double_separators(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_dir_spaces"
  logical :: success

  call remove_directory(root)
  ! CreateOutputFiles makes this path from an output directory given with a trailing slash
  call CreateDirectory(root//"/with space//Hams_format", success)
  call check(error, success, "CreateDirectory reported a failure")
  if (allocated(error)) return
  call check(error, accepts_files(root//"/with space/Hams_format"), "directory was not created")
  call remove_directory(root)
end subroutine test_spaces_and_double_separators

subroutine test_keeps_existing(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_dir_existing"
  logical :: success, kept
  integer :: unit

  call remove_directory(root)
  call CreateDirectory(root, success)
  call check(error, success, "CreateDirectory reported a failure")
  if (allocated(error)) return
  open(newunit=unit, file=root//"/existing.txt", status="replace")
  close(unit)

  call CreateDirectory(root, success)
  call check(error, success, "CreateDirectory failed on an existing directory")
  if (allocated(error)) return
  inquire(file=root//"/existing.txt", exist=kept)
  call check(error, kept, "existing file was removed")
  call remove_directory(root)
end subroutine test_keeps_existing

subroutine test_reports_failure(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_dir_failure"
  logical :: success
  integer :: unit

  call remove_directory(root)
  call CreateDirectory(root, success)
  call check(error, success, "CreateDirectory reported a failure")
  if (allocated(error)) return
  open(newunit=unit, file=root//"/plain_file", status="replace")
  close(unit)

  ! mkdir fails because a regular file sits where the parent directory should be
  call CreateDirectory(root//"/plain_file/sub", success)
  call check(error, .not. success, "failed mkdir was not reported")
  call remove_directory(root)
end subroutine test_reports_failure

logical function accepts_files(dir)
  character(len=*), intent(in) :: dir
  integer :: unit, ios

  open(newunit=unit, file=dir//"/.test_probe", status="replace", iostat=ios)
  accepts_files = (ios == 0)
  if (accepts_files) close(unit, status="delete")
end function accepts_files

subroutine remove_directory(root)
  character(len=*), intent(in) :: root
  character(len=32) :: os_name
  integer :: estat

  call get_environment_variable("OS", os_name)
  if (trim(os_name) == "Windows_NT") then
    call execute_command_line("if exist "//root//" rmdir /s /q "//root, exitstat=estat)
  else
    call execute_command_line("rm -rf "//root, exitstat=estat)
  end if
end subroutine remove_directory

end module test_create_directory
