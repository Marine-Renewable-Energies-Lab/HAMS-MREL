module test_integration_gmres

    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed, skip_test
    use testing_utilities

    implicit none
    public :: setup_gmres, collect_tests_gmres
    character(len=*), parameter :: test_data_dir = "test/data/config1_gmres/"
    character(len=*), parameter :: baseline_dir  = "test/data/config1/"
    logical :: simulation_execution_failed = .false.

    contains

    subroutine setup_gmres(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("run_simulation_gmres", run_simulation_gmres) &
        ]
    end subroutine setup_gmres

    subroutine collect_tests_gmres(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("test_oamass_gmres_vs_lu", test_oamass_gmres), &
            new_unittest("test_odamping_gmres_vs_lu", test_odamping_gmres) &
        ]
    end subroutine collect_tests_gmres

    subroutine run_simulation_gmres(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: inputdir, outputdir, command
        character(len=16) :: errcode
        integer :: exitcode

        inputdir  = get_repo_root_path()//test_data_dir//"input/"
        outputdir = get_running_exe_path()//"output_gmres/"
        command   = native_path(hamsexe_path//" "//inputdir//" "//outputdir//" > "//get_running_exe_path()//"gmres.log")//" 2>&1"

        print *, "Running simulation config1_gmres (ISOLV=2, H-matrix accelerated GMRES)"
        call execute_command_line(command, exitstat=exitcode, wait=.true.)
        if (exitcode /= 0) then
            simulation_execution_failed = .true.
            write(errcode, '(I0)') exitcode
            call test_failed(error, "Error: could not execute HAMS for GMRES integration test. Exit code: "//trim(errcode))
        end if
    end subroutine run_simulation_gmres

    ! Cross-check GMRES output vs the direct-LU baseline shipped with config1,
    ! using the scaled check of test_hams_matrices_are_close. GMRES uses
    ! ITER_TOL = 1e-4 (relative), so the limits are looser than the LU
    ! regression test because GMRES stops at its own tolerance.
    ! Measured largest scaled difference (ifort 2022.1): OAMASS 6.5e-6.
    subroutine test_oamass_gmres(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error

        if (simulation_execution_failed) then
            call skip_test(error, "Skipping OAMASS comparison (GMRES) because simulation failed")
            return
        end if

        call test_hams_matrices_are_close(get_repo_root_path()//baseline_dir//'expected_output/', &
            get_running_exe_path()//'output_gmres/Hams_format/', 'OAMASS', 1.0d-3, error)
    end subroutine test_oamass_gmres

    ! Damping tolerance is looser than added mass. At the low frequencies of
    ! this test case the damping peak (~0.3) sits six orders of magnitude below
    ! the added-mass scale (~1e5), so a GMRES residual of 1e-4 relative to the
    ! full solution scale legitimately produces a few-percent deviation in
    ! damping. Measured largest scaled difference is 3.1e-2 (surge-roll coupling
    ! between neighbouring bodies at 0.05 rad/s, zero by symmetry); 5e-2 still
    ! catches real solver failures (sign flips, garbage, order-of-magnitude errors).
    subroutine test_odamping_gmres(error)
        implicit none
        type(error_type), allocatable, intent(out) :: error

        if (simulation_execution_failed) then
            call skip_test(error, "Skipping ODAMPING comparison (GMRES) because simulation failed")
            return
        end if

        call test_hams_matrices_are_close(get_repo_root_path()//baseline_dir//'expected_output/', &
            get_running_exe_path()//'output_gmres/Hams_format/', 'ODAMPING', 5.0d-2, error)
    end subroutine test_odamping_gmres

end module test_integration_gmres
