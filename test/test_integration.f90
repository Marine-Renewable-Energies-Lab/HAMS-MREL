module test_integration_config1
  
    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed, skip_test
    use testing_utilities

    implicit none
    public :: setup_config1, collect_tests_config1
    character(len=*), parameter :: test_data_dir = "test/data/config1/"
    logical :: simulation_execution_failed = .false. ! Flag to control further tests

    ! Relative tolerances for comparing against the stored direct-LU outputs.
    ! Outputs differ between compilers and flags, so an absolute check cannot pass
    ! on every build. Measured on config1, largest relative difference per quantity:
    !                                     OAMASS     ODAMPING
    !   ifort 2022.1 vs stored outputs    3.4e-6     2.0e-4
    !   ifx 2022.1 vs stored outputs      6.8e-7     2.7e-4
    !   ifort 2022.1 vs ifx 2022.1        3.5e-6     3.1e-4
    !   Effect of the LINTON fix (eb6144a)    2.1e-4 to 8.7e-4 per file    up to 8.0e-4
    !   Effect of the image term integral     8.6e-5 to 5.7e-3 per file    up to 9.4e-1 (near-zero couplings)
    ! OAMASS is the sharp check: its limit is about 14x the build noise and 4x below
    ! the smallest change that fix made. Damping (~0.3) is small beside the added-mass
    ! scale (~1e5), so it carries more noise and its limit only catches larger errors.
    real(kind=8), parameter :: oamass_tolerance = 5.0d-5
    real(kind=8), parameter :: odamping_tolerance = 2.0d-3

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
        character(:), allocatable :: filename_expected, filename_actual
        character(len=1) :: istr
        integer :: i

        ! Skip test is simulation did not succeed
        if (simulation_execution_failed) then
            call skip_test(error, "Skipping OAMASS comparison for config1 because simulation failed")
            return
        end if

        do i = 1, 6
            ! Convert i to str
            write(istr, '(I1)') i
            ! Get full paths to expected and actual oamass datafiles
            filename_expected = get_repo_root_path()//test_data_dir//'expected_output/OAMASS'//istr//'.txt'
            filename_actual = get_running_exe_path()//'output_config1/Hams_format/OAMASS'//istr//'.txt'
            ! Check data in files is equal
            call test_oamass_files_are_close(filename_expected, filename_actual, oamass_tolerance, error)
            if (allocated(error)) return   ! stop further comparisons on first failure
        end do

    end subroutine test_oamass_config1

    subroutine test_odamping_config1(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: filename_expected, filename_actual
        character(len=1) :: istr
        integer :: i

        ! Skip test is simulation did not succeed
        if (simulation_execution_failed) then
            call skip_test(error, "Skipping ODAMPING comparison for config1 because simulation failed")
            return
        end if

        do i = 1, 6
            ! Convert i to str
            write(istr, '(I1)') i
            ! Get full paths to expected and actual oamass datafiles
            filename_expected = get_repo_root_path()//test_data_dir//'expected_output/ODAMPING'//istr//'.txt'
            filename_actual = get_running_exe_path()//'output_config1/Hams_format/ODAMPING'//istr//'.txt'
            ! Check data in files is equal
            ! Use OAMASS because files have same format
            call test_oamass_files_are_close(filename_expected, filename_actual, odamping_tolerance, error)
            if (allocated(error)) return   ! stop further comparisons on first failure
        end do

    end subroutine test_odamping_config1



end module test_integration_config1