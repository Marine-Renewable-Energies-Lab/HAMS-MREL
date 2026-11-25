module test_integration_config1
  
    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed, skip_test
    use testing_utilities

    implicit none
    public :: setup_config1, collect_tests_config1
    character(len=*), parameter :: test_data_dir = "test/data/config1/"
    logical :: simulation_execution_failed = .false. ! Flag to control further tests

    ! Tolerances on the scaled difference of test_hams_matrices_are_close (each entry
    ! divided by the geometric mean of its two self-block diagonals) against the
    ! stored direct-LU outputs. Outputs differ between compilers and flags, so an
    ! exact check cannot pass on every build. Measured on config1, largest value:
    !                                          OAMASS     ODAMPING
    !   ifort 2022.1 vs ifx 2022.1             1.2e-6     2.2e-4
    !   Change of g from 9.807 to 9.80665      6.5e-6     2.2e-4
    ! The damping noise sits in couplings that the array's symmetry makes zero
    ! (surge-heave, pitch-heave at 0.05 rad/s); physical entries agree to ~1e-5.
    ! Limits are about 8x (OAMASS) and 5x (ODAMPING) the largest value measured.
    real(kind=8), parameter :: oamass_tolerance = 5.0d-5
    real(kind=8), parameter :: odamping_tolerance = 1.0d-3

    contains

    !> Collect all exported tests for config1
    subroutine setup_config1(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("run_simulation_config1", run_simulation_config1) &
        ]
    end subroutine setup_config1

    !> Collect all exported tests for config1
    subroutine collect_tests_config1(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("test_oamass_config1", test_oamass_config1), &
            new_unittest("test_odamping_config1", test_odamping_config1) &
            !new_unittest("test_oexfor_config1", test_oexfor) &
        ]
    end subroutine collect_tests_config1

    subroutine run_simulation_config1(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: inputdir, outputdir, command
        character(len=16) :: errcode
        integer :: exitcode

        inputdir = get_repo_root_path()//test_data_dir//"input/"
        outputdir = get_running_exe_path()//"output_config1/"
        command = native_path(hamsexe_path//" "//inputdir//" "//outputdir//" > "//get_running_exe_path()//"config1.log")//" 2>&1"

        ! Run the simulation
        print *, "Running simulation config1"
        call execute_command_line(command, exitstat=exitcode, wait=.true.)
        if (exitcode /= 0) then
            simulation_execution_failed = .true.
            write(errcode, '(I0)') exitcode
            call test_failed(error, "Error: could not execute HAMS for integration test config1. Failed with exit code: "//trim(errcode))
        end if
    end subroutine run_simulation_config1

    subroutine test_oamass_config1(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error

        ! Skip test is simulation did not succeed
        if (simulation_execution_failed) then
            call skip_test(error, "Skipping OAMASS comparison for config1 because simulation failed")
            return
        end if

        call test_hams_matrices_are_close(get_repo_root_path()//test_data_dir//'expected_output/', &
            get_running_exe_path()//'output_config1/Hams_format/', 'OAMASS', oamass_tolerance, error)

    end subroutine test_oamass_config1

    subroutine test_odamping_config1(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error

        ! Skip test is simulation did not succeed
        if (simulation_execution_failed) then
            call skip_test(error, "Skipping ODAMPING comparison for config1 because simulation failed")
            return
        end if

        call test_hams_matrices_are_close(get_repo_root_path()//test_data_dir//'expected_output/', &
            get_running_exe_path()//'output_config1/Hams_format/', 'ODAMPING', odamping_tolerance, error)

    end subroutine test_odamping_config1



end module test_integration_config1