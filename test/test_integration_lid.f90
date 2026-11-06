!> Whole-solver checks of the damping lid input: the stops, and that a lid switched off changes nothing. The cases are described
!> in test/data/damping_lid/README.md.
module test_integration_lid

    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
    use testing_utilities

    implicit none
    private
    public :: setup_lid, collect_tests_lid

contains

    subroutine setup_lid(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("run_lid_flag_off", run_lid_flag_off), &
            new_unittest("run_lid_reference", run_lid_reference) &
            ]
    end subroutine setup_lid

    subroutine collect_tests_lid(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("stops_lid_missing", stops_lid_missing), &
            new_unittest("stops_lid_offplane", stops_lid_offplane), &
            new_unittest("stops_lid_symmetry", stops_lid_symmetry), &
            new_unittest("stops_lid_gmres", stops_lid_gmres), &
            new_unittest("stops_lid_inside", stops_lid_inside), &
            new_unittest("flag_off_matches_no_lid", flag_off_matches_no_lid) &
            ]
    end subroutine collect_tests_lid

    ! -------------------------------------------------------------------------------------------- runs

    function input_dir(group, name) result(dir)
        character(len=*), intent(in) :: group, name
        character(:), allocatable :: dir
        dir = get_repo_root_path()//"test/data/"//group//"/"//name//"/input/"
    end function input_dir

    !> The output goes to folders of this suite, so the results of the image suite are not overwritten.
    function output_dir(name) result(dir)
        character(len=*), intent(in) :: name
        character(:), allocatable :: dir
        dir = get_running_exe_path()//"output_lid_"//name//"/"
    end function output_dir

    function log_file(name) result(file)
        character(len=*), intent(in) :: name
        character(:), allocatable :: file
        file = get_running_exe_path()//"lid_"//name//".log"
    end function log_file

    subroutine run_hams(group, name, exitcode)
        character(len=*), intent(in) :: group, name
        integer, intent(out) :: exitcode
        character(:), allocatable :: command

        command = native_path(hamsexe_path//" "//input_dir(group, name)//" "//output_dir(name)//" > "//log_file(name))//" 2>&1"
        call execute_command_line(command, exitstat=exitcode, wait=.true.)
    end subroutine run_hams

    subroutine run_case(group, name, error)
        character(len=*), intent(in) :: group, name
        type(error_type), allocatable, intent(out) :: error
        character(len=16) :: code
        integer :: exitcode

        call run_hams(group, name, exitcode)
        if (exitcode /= 0) then
            write (code, '(I0)') exitcode
            call test_failed(error, "HAMS failed for "//name//" with exit code "//trim(code))
        end if
    end subroutine run_case

    !> The run must stop with a non-zero exit code, and say why in its log.
    subroutine run_case_expect_stop(name, message, error)
        character(len=*), intent(in) :: name, message
        type(error_type), allocatable, intent(out) :: error
        character(len=512) :: line
        integer :: exitcode, unit, ios
        logical :: said

        call run_hams("damping_lid", name, exitcode)
        call check(error, exitcode /= 0, name//": the run did not stop")
        if (allocated(error)) return
        said = .false.
        open (newunit=unit, file=log_file(name), status='old', action='read')
        do
            read (unit, '(a)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, message) > 0) said = .true.
        end do
        close (unit)
        call check(error, said, name//": the message '"//message//"' is not in the log")
    end subroutine run_case_expect_stop

    subroutine run_lid_flag_off(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("damping_lid", "lid_flag_off", error)
    end subroutine run_lid_flag_off

    !> The same cylinder with no lid line at all, from the folder of the image term cases.
    subroutine run_lid_reference(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("image_term", "image_cos", error)
    end subroutine run_lid_reference

    ! -------------------------------------------------------------------------------------------- tests

    subroutine stops_lid_missing(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_missing", "DampingLid.pnl input file missing", error)
    end subroutine stops_lid_missing

    subroutine stops_lid_offplane(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_offplane", "lid node is not on z = 0", error)
    end subroutine stops_lid_offplane

    subroutine stops_lid_symmetry(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_symmetry", "lid symmetry differs from the hull", error)
    end subroutine stops_lid_symmetry

    subroutine stops_lid_gmres(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_gmres", "damping lid needs Solver_type 1", error)
    end subroutine stops_lid_gmres

    subroutine stops_lid_inside(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_inside", "lid panel inside a hull", error)
    end subroutine stops_lid_inside

    !> Damping_lid 0 with a lid file present gives the added mass, damping and excitation of the run without the line, to the
    !> last printed digit.
    subroutine flag_off_matches_no_lid(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=10), parameter :: files(2) = [character(len=10) :: "AmssDamp.1", "ExcForce.3"]
        integer :: k

        do k = 1, size(files)
            call same_text(output_dir("lid_flag_off")//"Wamit_format/"//trim(files(k)), &
                           output_dir("image_cos")//"Wamit_format/"//trim(files(k)), error)
            if (allocated(error)) return
        end do
    end subroutine flag_off_matches_no_lid

    subroutine same_text(file_a, file_b, error)
        character(len=*), intent(in) :: file_a, file_b
        type(error_type), allocatable, intent(out) :: error
        character(len=512) :: line_a, line_b
        integer :: unit_a, unit_b, stat_a, stat_b, n

        open (newunit=unit_a, file=file_a, status="old", action="read", iostat=stat_a)
        open (newunit=unit_b, file=file_b, status="old", action="read", iostat=stat_b)
        if (stat_a /= 0 .or. stat_b /= 0) then
            call test_failed(error, "Cannot read "//file_a//" and "//file_b)
            return
        end if
        n = 0
        do
            read (unit_a, '(a)', iostat=stat_a) line_a
            read (unit_b, '(a)', iostat=stat_b) line_b
            if (stat_a /= 0 .or. stat_b /= 0) exit
            n = n + 1
            if (line_a /= line_b) then
                call test_failed(error, file_a//" differs from the no-lid run")
                exit
            end if
        end do
        if (.not. allocated(error)) call check(error, (stat_a /= 0) .eqv. (stat_b /= 0), file_a//": different number of rows")
        if (.not. allocated(error)) call check(error, n > 0, file_a//" is empty")
        close (unit_a)
        close (unit_b)
    end subroutine same_text

end module test_integration_lid
