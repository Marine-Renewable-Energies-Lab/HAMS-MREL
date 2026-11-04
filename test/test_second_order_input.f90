!> Unit tests for the optional SecondOrder.in reader and the second-order option validation.
module test_second_order_input
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use IO, only : CreateDirectory
  use SecondOrderTypes, only : so_options
  use SecondOrderInput, only : read_second_order_options, validate_second_order
  implicit none
  private

  public :: collect_second_order_input

  ! test-drive runs tests in parallel under OpenMP, so each test uses its own directory

contains

subroutine collect_second_order_input(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("absent_file_inactive", test_absent_file), &
    new_unittest("defaults_applied", test_defaults), &
    new_unittest("all_keys_parsed", test_all_keys), &
    new_unittest("unknown_key_rejected", test_unknown_key), &
    new_unittest("bad_value_rejected", test_bad_value), &
    new_unittest("validate_rejects", test_validate_rejects), &
    new_unittest("validate_accepts", test_validate_accepts) &
    ]

end subroutine collect_second_order_input

subroutine test_absent_file(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_absent"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg

  call remove_directory(root)
  call make_directory(root)
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "absent file reported as an error: "//trim(msg))
  if (allocated(error)) return
  call check(error, .not. opts%active, "options active without a file")
  call remove_directory(root)
end subroutine test_absent_file

subroutine test_defaults(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_defaults"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: "mean_drift = 1"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "read failed: "//trim(msg))
  if (allocated(error)) return
  call check(error, opts%active, "options not active")
  if (allocated(error)) return
  call check(error, opts%mean_drift, 1, "mean_drift")
  if (allocated(error)) return
  call check(error, opts%n_azimuth, 360, "n_azimuth default")
  if (allocated(error)) return
  call check(error, opts%velocity_method, 1, "velocity_method default")
  if (allocated(error)) return
  call check(error, opts%do_far .and. .not. opts%do_near, "drift_methods default is far field")
  call remove_directory(root)
end subroutine test_defaults

subroutine test_all_keys(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_allkeys"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=60) :: &
    "! a comment line", &
    "wl_z_tolerance = 2.5d-6", &
    "write_terms = 2   ! trailing comment", &
    "far_tolerance = 1.0d-5", &
    "n_azimuth = 180", &
    "motion_state = fixed", &
    "waterline_offset = 0.02d0", &
    "waterline_method = field", &
    "velocity_method = source", &
    "drift_methods = far near", &
    "mean_drift = 2"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "read failed: "//trim(msg))
  if (allocated(error)) return
  call check(error, opts%active, "active")
  call check(error, opts%mean_drift, 2, "mean_drift")
  call check(error, opts%do_far .and. opts%do_near, "drift_methods")
  call check(error, opts%velocity_method, 2, "velocity_method")
  call check(error, opts%waterline_method, 2, "waterline_method")
  call check(error, abs(opts%waterline_offset - 0.02d0) < 1d-15, "waterline_offset")
  call check(error, opts%motion_fixed, "motion_state")
  call check(error, opts%n_azimuth, 180, "n_azimuth")
  call check(error, abs(opts%far_tolerance - 1d-5) < 1d-17, "far_tolerance")
  call check(error, opts%write_terms, 2, "write_terms")
  call check(error, abs(opts%wl_z_tolerance - 2.5d-6) < 1d-18, "wl_z_tolerance")
  call remove_directory(root)
end subroutine test_all_keys

subroutine test_unknown_key(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_unknown"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_mehods = far"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, .not. ok, "unknown key accepted")
  if (allocated(error)) return
  call check(error, index(msg, "drift_mehods") > 0, "message does not name the key: "//trim(msg))
  call remove_directory(root)
end subroutine test_unknown_key

subroutine test_bad_value(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_badvalue"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: "mean_drift = 1", "velocity_method = fast"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, .not. ok, "bad value accepted")
  if (allocated(error)) return
  call check(error, index(msg, "velocity_method") > 0, "message does not name the key: "//trim(msg))
  call remove_directory(root)
end subroutine test_bad_value

subroutine test_validate_rejects(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: far, near, src
  logical :: ok
  character(len=200) :: msg

  far%active = .true.; far%mean_drift = 1; far%do_far = .true.
  near = far; near%do_near = .true.
  src = far; src%velocity_method = 2

  call validate_second_order(far, 2, -1d0, 1, 0, .false., ok, msg)
  call check(error, .not. ok, "isol=2 accepted")
  if (allocated(error)) return
  call check(error, index(msg, "ISOL") > 0, "isol message: "//trim(msg))
  if (allocated(error)) return

  call validate_second_order(far, 1, 600d0, 1, 0, .false., ok, msg)
  call check(error, .not. ok, "h=600 accepted")
  if (allocated(error)) return
  call check(error, index(msg, "H<0") > 0, "depth message: "//trim(msg))
  if (allocated(error)) return

  call validate_second_order(near, 1, -1d0, 1, 0, .true., ok, msg)
  call check(error, .not. ok, "near field with generalized modes accepted")
  if (allocated(error)) return
  call check(error, index(msg, "gen_mod") > 0, "generalized-mode message: "//trim(msg))
  if (allocated(error)) return

  call validate_second_order(src, 1, -1d0, 2, 0, .false., ok, msg)
  call check(error, .not. ok, "source velocity with GMRES accepted")
  if (allocated(error)) return
  call check(error, index(msg, "ISOLV") > 0, "isolv message: "//trim(msg))
  if (allocated(error)) return

  call validate_second_order(src, 1, -1d0, 1, 1, .false., ok, msg)
  call check(error, .not. ok, "source velocity with irregular frequency removal accepted")
  if (allocated(error)) return
  call check(error, index(msg, "IRSP") > 0, "irsp message: "//trim(msg))
end subroutine test_validate_rejects

subroutine test_validate_accepts(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: far, near
  logical :: ok
  character(len=200) :: msg

  far%active = .true.; far%mean_drift = 1; far%do_far = .true.
  near = far; near%do_near = .true.

  call validate_second_order(near, 1, -1d0, 1, 0, .false., ok, msg)
  call check(error, ok, "infinite depth rejected: "//trim(msg))
  if (allocated(error)) return
  call validate_second_order(near, 1, 70d0, 1, 0, .false., ok, msg)
  call check(error, ok, "h=70 rejected: "//trim(msg))
  if (allocated(error)) return
  call validate_second_order(far, 1, 70d0, 1, 0, .true., ok, msg)
  call check(error, ok, "far field with generalized modes rejected: "//trim(msg))
end subroutine test_validate_accepts

subroutine make_directory(root)
  character(len=*), intent(in) :: root
  logical :: success

  call CreateDirectory(root, success)
  if (.not. success) error stop "could not create test directory"
end subroutine make_directory

subroutine write_input(root, lines)
  character(len=*), intent(in) :: root, lines(:)
  integer :: unit, i

  open(newunit=unit, file=root//"/SecondOrder.in", status="replace")
  do i = 1, size(lines)
    write(unit, '(a)') trim(lines(i))
  end do
  close(unit)
end subroutine write_input

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

end module test_second_order_input
