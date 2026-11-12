!> Unit tests for the optional SecondOrder.in reader and the second-order option validation.
module test_second_order_input
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use IO, only : CreateDirectory
  use SecondOrderTypes, only : so_options, QTF_TERMS_QUADRATIC, QTF_TERMS_BODY, QTF_TERMS_FULL, QTF_B2_BODY, QTF_B2_CONTROL
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
    new_unittest("waterline_offset_positive", test_offset_positive), &
    new_unittest("far_tolerance_positive", test_far_tolerance_positive), &
    new_unittest("wl_z_tolerance_positive", test_wl_z_tolerance_positive), &
    new_unittest("validate_rejects", test_validate_rejects), &
    new_unittest("validate_accepts", test_validate_accepts), &
    new_unittest("control_keys_parse", test_control_keys), &
    new_unittest("control_defaults", test_control_defaults), &
    new_unittest("control_rejections", test_control_rejections), &
    new_unittest("qtf_keys", test_qtf_keys), &
    new_unittest("qtf_terms_body", test_qtf_terms_body), &
    new_unittest("qtf_terms_full_accepted", test_qtf_terms_full_accepted), &
    new_unittest("qtf_full_deep_water_stops", test_qtf_full_deep_water_stops), &
    new_unittest("qtf_full_sum_only_stops", test_qtf_full_sum_only_stops), &
    new_unittest("qtf_fs_keys_ranges", test_qtf_fs_keys_ranges), &
    new_unittest("qtf_fs_keys_ignored_note", test_qtf_fs_keys_ignored_note), &
    new_unittest("qtf_b2_key", test_qtf_b2_key), &
    new_unittest("qtf_edge_keys", test_qtf_edge_keys), &
    new_unittest("qtf_field_key", test_qtf_field_key) &
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

subroutine test_offset_positive(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_offset"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: "mean_drift = 1", "waterline_offset = 0"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, .not. ok, "an offset of zero accepted")
  if (allocated(error)) return
  call check(error, index(msg, "waterline_offset") > 0, "message does not name the key: "//trim(msg))
  call remove_directory(root)
end subroutine test_offset_positive

subroutine test_far_tolerance_positive(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_ftol"
  character(len=*), parameter :: values(2) = [character(len=6) :: "0", "-1d-6"]
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg
  integer :: k

  do k = 1, 2
    call remove_directory(root)
    call make_directory(root)
    call write_input(root, [character(len=40) :: "mean_drift = 1", "far_tolerance = "//trim(values(k))])
    call read_second_order_options(root, opts, ok, msg)
    call check(error, .not. ok, "far_tolerance = "//trim(values(k))//" accepted")
    if (allocated(error)) return
    call check(error, index(msg, "far_tolerance") > 0, "message does not name the key: "//trim(msg))
    if (allocated(error)) return
    call remove_directory(root)
  end do
end subroutine test_far_tolerance_positive

subroutine test_wl_z_tolerance_positive(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_wltol"
  character(len=*), parameter :: values(2) = [character(len=6) :: "0", "-1d-6"]
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg
  integer :: k

  do k = 1, 2
    call remove_directory(root)
    call make_directory(root)
    call write_input(root, [character(len=40) :: "mean_drift = 1", "wl_z_tolerance = "//trim(values(k))])
    call read_second_order_options(root, opts, ok, msg)
    call check(error, .not. ok, "wl_z_tolerance = "//trim(values(k))//" accepted")
    if (allocated(error)) return
    call check(error, index(msg, "wl_z_tolerance") > 0, "message does not name the key: "//trim(msg))
    if (allocated(error)) return
    call remove_directory(root)
  end do
end subroutine test_wl_z_tolerance_positive

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

subroutine test_control_keys(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_ctrlkeys"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: &
    "mean_drift = 1", &
    "drift_methods = far control", &
    "cs_velocity = source", &
    "cs_shape = box", &
    "cs_margin = 12.5", &
    "cs_fd_step = 0.05"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "read failed: "//trim(msg))
  if (allocated(error)) return
  call check(error, opts%do_control .and. opts%do_far .and. .not. opts%do_near, "drift_methods")
  if (allocated(error)) return
  call check(error, opts%cs_velocity, 2, "cs_velocity")
  if (allocated(error)) return
  call check(error, opts%cs_shape, 2, "cs_shape")
  if (allocated(error)) return
  call check(error, abs(opts%cs_margin - 12.5d0) < 1d-12, "cs_margin")
  if (allocated(error)) return
  call check(error, abs(opts%cs_fd_step - 0.05d0) < 1d-12, "cs_fd_step")
  call remove_directory(root)
end subroutine test_control_keys

subroutine test_control_defaults(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_ctrldef"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = control"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "read failed: "//trim(msg))
  if (allocated(error)) return
  call check(error, opts%do_control .and. .not. opts%do_far .and. .not. opts%do_near, "control alone")
  if (allocated(error)) return
  call check(error, opts%cs_velocity, 1, "cs_velocity default")
  if (allocated(error)) return
  call check(error, opts%cs_shape, 1, "cs_shape default")
  if (allocated(error)) return
  call check(error, opts%cs_margin == 0d0, "cs_margin default")
  if (allocated(error)) return
  call check(error, opts%cs_bottom_margin == 0d0, "cs_bottom_margin default")
  if (allocated(error)) return
  call check(error, opts%cs_panel_size == 0d0, "cs_panel_size default")
  if (allocated(error)) return
  call check(error, opts%cs_fs_rows, 0, "cs_fs_rows default")
  if (allocated(error)) return
  call check(error, opts%cs_fd_step == 0.02d0, "cs_fd_step default")
  call remove_directory(root)
end subroutine test_control_defaults

subroutine test_control_rejections(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: ctl, src, bad
  logical :: ok
  character(len=200) :: msg

  ctl%active = .true.; ctl%mean_drift = 1; ctl%do_control = .true.
  src = ctl; src%cs_velocity = 2

  call validate_second_order(ctl, 1, -1d0, 1, 0, .true., ok, msg)
  call check(error, .not. ok, "control with generalized modes accepted")
  if (allocated(error)) return
  call check(error, trim(msg) == "drift_methods = control cannot be used with generalized modes.", &
    "generalized-mode message: "//trim(msg))
  if (allocated(error)) return

  call validate_second_order(src, 1, -1d0, 2, 0, .false., ok, msg)
  call check(error, .not. ok, "cs_velocity source with GMRES accepted")
  if (allocated(error)) return
  call check(error, trim(msg) == "cs_velocity = source needs the direct solver (ISOLV=1) in ControlFile.in.", &
    "isolv message: "//trim(msg))
  if (allocated(error)) return

  call validate_second_order(src, 1, -1d0, 1, 1, .false., ok, msg)
  call check(error, .not. ok, "cs_velocity source with IRSP accepted")
  if (allocated(error)) return
  call check(error, trim(msg) == "cs_velocity = source cannot be used with irregular frequency removal (IRSP=1) in ControlFile.in.", &
    "irsp message: "//trim(msg))
  if (allocated(error)) return

  bad = ctl; bad%cs_margin = -1d0
  call expect_reject(error, bad, "cs_margin"); if (allocated(error)) return
  bad = ctl; bad%cs_bottom_margin = -1d0
  call expect_reject(error, bad, "cs_bottom_margin"); if (allocated(error)) return
  bad = ctl; bad%cs_panel_size = -1d0
  call expect_reject(error, bad, "cs_panel_size"); if (allocated(error)) return
  bad = ctl; bad%cs_fd_step = 0d0
  call expect_reject(error, bad, "cs_fd_step"); if (allocated(error)) return
  bad = ctl; bad%cs_fd_step = 0.6d0
  call expect_reject(error, bad, "cs_fd_step"); if (allocated(error)) return
  bad = ctl; bad%cs_fs_rows = -1
  call expect_reject(error, bad, "cs_fs_rows"); if (allocated(error)) return

  bad = ctl; bad%cs_fd_step = 0.5d0
  call validate_second_order(bad, 1, -1d0, 1, 0, .false., ok, msg)
  call check(error, ok, "cs_fd_step = 0.5 rejected: "//trim(msg))
end subroutine test_control_rejections

subroutine expect_reject(error, opts, key)
  type(error_type), allocatable, intent(out) :: error
  type(so_options), intent(in) :: opts
  character(len=*), intent(in) :: key
  logical :: ok
  character(len=200) :: msg

  call validate_second_order(opts, 1, -1d0, 1, 0, .false., ok, msg)
  call check(error, .not. ok, key//" out of range accepted")
  if (allocated(error)) return
  call check(error, index(msg, key) > 0, "message does not name "//key//": "//trim(msg))
end subroutine expect_reject

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

subroutine test_qtf_keys(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_qtf"
  type(so_options) :: opts, far_only
  logical :: ok
  character(len=200) :: msg

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", "qtf = 3", &
                          "qtf_terms = quadratic", "qtf_max_dw = 0.3"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "read failed: "//trim(msg))
  if (allocated(error)) return
  call check(error, opts%qtf, 3, "qtf")
  if (allocated(error)) return
  call check(error, abs(opts%qtf_max_dw - 0.3d0) < 1d-12, "qtf_max_dw")
  if (allocated(error)) return
  call validate_second_order(opts, 1, -1d0, 1, 0, .false., ok, msg)
  call check(error, ok, "qtf with near rejected: "//trim(msg))
  if (allocated(error)) return

  call expect_read_message(error, root, [character(len=40) :: "mean_drift = 1", "qtf = 4"], &
                           "qtf must be 0, 1, 2 or 3 in SecondOrder.in.")
  if (allocated(error)) return
  call expect_read_message(error, root, [character(len=40) :: "mean_drift = 1", "qtf = 1", "qtf_max_dw = -1"], &
                           "qtf_max_dw must not be negative in SecondOrder.in.")
  if (allocated(error)) return
  call expect_read_message(error, root, [character(len=40) :: "mean_drift = 0", "qtf = 1"], "qtf needs mean_drift = 1 or 2.")
  if (allocated(error)) return
  call expect_read_message(error, root, [character(len=40) :: "mean_drift = 1", "qtf_terms = all"], &
                           'Invalid value "all" for key "qtf_terms" in SecondOrder.in.')
  if (allocated(error)) return

  far_only%active = .true.; far_only%mean_drift = 1; far_only%do_far = .true.; far_only%qtf = 1
  call validate_second_order(far_only, 1, -1d0, 1, 0, .false., ok, msg)
  call check(error, .not. ok, "qtf without near accepted")
  if (allocated(error)) return
  call check(error, trim(msg) == "qtf needs drift_methods to include near.", "near message: "//trim(msg))
  call remove_directory(root)
end subroutine test_qtf_keys

!> qtf_field is 0 (default) or 1, and 1 needs a QTF.
subroutine test_qtf_field_key(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_qtf_field"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", "qtf = 3"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "read failed: "//trim(msg))
  if (allocated(error)) return
  call check(error, opts%qtf_field, 0, "qtf_field default")
  if (allocated(error)) return
  call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", "qtf = 3", "qtf_field = 1"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "read failed: "//trim(msg))
  if (allocated(error)) return
  call check(error, opts%qtf_field, 1, "qtf_field")
  if (allocated(error)) return
  call expect_read_message(error, root, [character(len=40) :: "mean_drift = 1", "qtf = 3", "qtf_field = 2"], &
                           "qtf_field must be 0 or 1 in SecondOrder.in.")
  if (allocated(error)) return
  call expect_read_message(error, root, [character(len=40) :: "mean_drift = 1", "qtf = 0", "qtf_field = 1"], &
                           "qtf_field needs qtf = 1, 2 or 3 in SecondOrder.in.")
  call remove_directory(root)
end subroutine test_qtf_field_key

!> qtf_terms = body is read with the difference QTF (qtf = 1 or 3) and rejected with the sum frequency only; quadratic is
!> the default.
subroutine test_qtf_terms_body(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_qtf_body"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg
  integer :: q

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", "qtf = 1"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "read failed: "//trim(msg))
  if (allocated(error)) return
  call check(error, opts%qtf_terms, QTF_TERMS_QUADRATIC, "qtf_terms default")
  if (allocated(error)) return
  do q = 1, 3, 2
    write (msg, '(a,I0)') "qtf = ", q
    call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", trim(msg), &
                            "qtf_terms = body"])
    call read_second_order_options(root, opts, ok, msg)
    call check(error, ok, "qtf_terms = body rejected: "//trim(msg))
    if (allocated(error)) return
    call check(error, opts%qtf_terms, QTF_TERMS_BODY, "qtf_terms = body")
    if (allocated(error)) return
    call validate_second_order(opts, 1, -1d0, 1, 0, .false., ok, msg)
    call check(error, ok, "qtf_terms = body not validated: "//trim(msg))
    if (allocated(error)) return
  end do
  call expect_read_message(error, root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", "qtf = 2", &
                                         "qtf_terms = body"], &
                           "qtf_terms = body needs the difference QTF (qtf = 1 or 3); sum-frequency body forcing needs "// &
                           "the free-surface integrals (Phase 6).")
  call remove_directory(root)
end subroutine test_qtf_terms_body

!> qtf_terms = full is read with qtf = 1 and 3 and validated in finite depth; the free-surface keys keep their defaults
subroutine test_qtf_terms_full_accepted(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_qtf_full"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg
  integer :: q

  call remove_directory(root)
  call make_directory(root)
  do q = 1, 3, 2
    write (msg, '(a,I0)') "qtf = ", q
    call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", trim(msg), &
                            "qtf_terms = full"])
    call read_second_order_options(root, opts, ok, msg)
    call check(error, ok, "qtf_terms = full rejected: "//trim(msg))
    if (allocated(error)) return
    call check(error, opts%qtf_terms, QTF_TERMS_FULL, "qtf_terms = full")
    if (allocated(error)) return
    call check(error, opts%qtf_fs_radius == 0d0 .and. opts%qtf_fs_rows == 0 .and. opts%qtf_fs_modes == 0 .and. &
               opts%qtf_fs_evanescent == 0 .and. .not. opts%qtf_fs_given, "free-surface key defaults")
    if (allocated(error)) return
    call validate_second_order(opts, 1, 100d0, 1, 0, .false., ok, msg)
    call check(error, ok, "qtf_terms = full not validated: "//trim(msg))
    if (allocated(error)) return
  end do
  call remove_directory(root)
end subroutine test_qtf_terms_full_accepted

!> qtf_terms = full with H < 0 stops in the validation, with the message of the plan
subroutine test_qtf_full_deep_water_stops(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_qtf_full_deep"
  type(so_options) :: opts
  logical :: ok
  character(len=300) :: msg

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", "qtf = 1", "qtf_terms = full"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "read failed: "//trim(msg))
  if (allocated(error)) return
  call validate_second_order(opts, 1, -1d0, 1, 0, .false., ok, msg)
  call check(error, .not. ok, "deep water accepted with qtf_terms = full")
  if (allocated(error)) return
  call check(error, trim(msg) == "qtf_terms = full needs a finite water depth (the outer free-surface integral uses the "// &
             "finite-depth modes); set a large finite depth instead.", "message: "//trim(msg))
  call remove_directory(root)
end subroutine test_qtf_full_deep_water_stops

!> qtf_terms = full with the sum frequency only stops at reading
subroutine test_qtf_full_sum_only_stops(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_qtf_full_sum"

  call remove_directory(root)
  call make_directory(root)
  call expect_read_message(error, root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", "qtf = 2", &
                                         "qtf_terms = full"], &
                           "qtf_terms = full needs the difference QTF (qtf = 1 or 3); the sum frequency is Phase 6b.")
  call remove_directory(root)
end subroutine test_qtf_full_sum_only_stops

!> The free-surface keys are read with their values and marked as given; negative values stop
subroutine test_qtf_fs_keys_ranges(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_qtf_fs"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg
  character(len=40), parameter :: base(3) = [character(len=40) :: "mean_drift = 1", "qtf = 1", "qtf_terms = full"]

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: base, "qtf_fs_radius = 12.5", "qtf_fs_rows = 20", "qtf_fs_modes = 16", &
                          "qtf_fs_evanescent = 4"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "keys rejected: "//trim(msg))
  if (allocated(error)) return
  call check(error, opts%qtf_fs_radius == 12.5d0 .and. opts%qtf_fs_rows == 20 .and. opts%qtf_fs_modes == 16 .and. &
             opts%qtf_fs_evanescent == 4 .and. opts%qtf_fs_given, "key values")
  if (allocated(error)) return
  call expect_read_message(error, root, [character(len=40) :: base, "qtf_fs_radius = -1"], &
                           "qtf_fs_radius must not be negative in SecondOrder.in.")
  if (allocated(error)) return
  call expect_read_message(error, root, [character(len=40) :: base, "qtf_fs_rows = -1"], &
                           "qtf_fs_rows must not be negative in SecondOrder.in.")
  if (allocated(error)) return
  call expect_read_message(error, root, [character(len=40) :: base, "qtf_fs_modes = -2"], &
                           "qtf_fs_modes must not be negative in SecondOrder.in.")
  if (allocated(error)) return
  call expect_read_message(error, root, [character(len=40) :: base, "qtf_fs_evanescent = -1"], &
                           "qtf_fs_evanescent must not be negative in SecondOrder.in.")
  if (allocated(error)) return
  call expect_read_message(error, root, [character(len=40) :: base, "qtf_fs_modes = many"], &
                           'Invalid value "many" for key "qtf_fs_modes" in SecondOrder.in.')
  call remove_directory(root)
end subroutine test_qtf_fs_keys_ranges

!> The keys are marked as given also without qtf_terms = full (here with body); the driver then notes in the log that they
!> are ignored, as for qtf_b2. Not given: no note.
subroutine test_qtf_fs_keys_ignored_note(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_qtf_fs_note"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: "mean_drift = 1", "qtf = 1", "qtf_terms = body", "qtf_fs_modes = 12"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "read failed: "//trim(msg))
  if (allocated(error)) return
  call check(error, opts%qtf_fs_given .and. opts%qtf_terms == QTF_TERMS_BODY .and. opts%qtf_fs_modes == 12, &
             "given with body: ignored, marked for the note")
  if (allocated(error)) return
  call write_input(root, [character(len=40) :: "mean_drift = 1", "qtf = 1", "qtf_terms = body"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok .and. .not. opts%qtf_fs_given, "not given: no note")
  call remove_directory(root)
end subroutine test_qtf_fs_keys_ignored_note

!> qtf_b2 is control by default, reads body and control, and rejects other values
subroutine test_qtf_b2_key(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_qtf_b2"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", "qtf = 1", &
                          "qtf_terms = body"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "read failed: "//trim(msg))
  if (allocated(error)) return
  call check(error, opts%qtf_b2, QTF_B2_CONTROL, "qtf_b2 default")
  if (allocated(error)) return
  call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", "qtf = 1", &
                          "qtf_terms = body", "qtf_b2 = body"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "qtf_b2 = body rejected: "//trim(msg))
  if (allocated(error)) return
  call check(error, opts%qtf_b2, QTF_B2_BODY, "qtf_b2 = body")
  if (allocated(error)) return
  call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", "qtf = 1", &
                          "qtf_terms = body", "QTF_B2 = Control"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "qtf_b2 = control rejected: "//trim(msg))
  if (allocated(error)) return
  call check(error, opts%qtf_b2, QTF_B2_CONTROL, "qtf_b2 = control")
  if (allocated(error)) return
  call expect_read_message(error, root, [character(len=40) :: "mean_drift = 1", "qtf = 1", "qtf_terms = body", &
                                         "qtf_b2 = hull"], 'Invalid value "hull" for key "qtf_b2" in SecondOrder.in.')
  call remove_directory(root)
end subroutine test_qtf_b2_key

!> qtf_edge_weights is off by default and reads on and off; qtf_edge_angle is 45 by default and must lie in (0, 180). Each
!> key, and qtf_b2, is marked as given, for the note in the log when it has no effect
subroutine test_qtf_edge_keys(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "hams_test_so_qtf_edge"
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg

  call remove_directory(root)
  call make_directory(root)
  call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", "qtf = 1", &
                          "qtf_terms = body"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "read failed: "//trim(msg))
  if (allocated(error)) return
  call check(error, .not. opts%qtf_edge_weights .and. opts%qtf_edge_angle == 45d0 .and. .not. opts%qtf_edge_given &
             .and. .not. opts%qtf_b2_given, "defaults: off, 45 degrees, not given")
  if (allocated(error)) return
  call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", "qtf = 1", &
                          "qtf_terms = body", "QTF_Edge_Weights = Off", "qtf_edge_angle = 30", "qtf_b2 = body"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok, "qtf_edge_weights = off rejected: "//trim(msg))
  if (allocated(error)) return
  call check(error, .not. opts%qtf_edge_weights .and. opts%qtf_edge_angle == 30d0 .and. opts%qtf_edge_given &
             .and. opts%qtf_b2_given, "off, 30 degrees, given")
  if (allocated(error)) return
  call write_input(root, [character(len=40) :: "mean_drift = 1", "drift_methods = far near", "qtf = 1", &
                          "qtf_terms = body", "qtf_edge_weights = on"])
  call read_second_order_options(root, opts, ok, msg)
  call check(error, ok .and. opts%qtf_edge_weights, "qtf_edge_weights = on")
  if (allocated(error)) return
  call expect_read_message(error, root, [character(len=40) :: "mean_drift = 1", "qtf = 1", "qtf_terms = body", &
                                         "qtf_edge_weights = yes"], &
                           'Invalid value "yes" for key "qtf_edge_weights" in SecondOrder.in.')
  if (allocated(error)) return
  call expect_read_message(error, root, [character(len=40) :: "mean_drift = 1", "qtf = 1", "qtf_terms = body", &
                                         "qtf_edge_angle = 0"], &
                           'qtf_edge_angle must be between 0 and 180 degrees in SecondOrder.in.')
  call remove_directory(root)
end subroutine test_qtf_edge_keys

subroutine expect_read_message(error, root, lines, expected)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), intent(in) :: root, lines(:), expected
  type(so_options) :: opts
  logical :: ok
  character(len=200) :: msg

  call write_input(root, lines)
  call read_second_order_options(root, opts, ok, msg)
  call check(error, .not. ok, "accepted: "//expected)
  if (allocated(error)) return
  call check(error, trim(msg) == expected, "message: "//trim(msg))
end subroutine expect_read_message

end module test_second_order_input
