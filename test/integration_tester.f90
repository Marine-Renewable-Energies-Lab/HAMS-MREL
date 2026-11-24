program integration_tester
  
    ! Testing framework imports
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use testdrive_version, only : get_testdrive_version

    ! Tests imports (add only collect functions)
    use test_integration_config1, only : setup_config1, collect_tests_config1
    use test_integration_gmres,   only : setup_gmres,   collect_tests_gmres
    use test_integration_image,   only : setup_image,   collect_tests_image
    use test_integration_gen_modes, only : setup_gen_modes, collect_tests_gen_modes
    use test_integration_drift,   only : setup_drift,   collect_tests_drift
    use test_integration_drift,   only : setup_drift_motion, collect_tests_drift_motion
    use test_integration_array,   only : setup_array,   collect_tests_array
    use test_integration_dipole,  only : setup_dipole,  collect_tests_dipole
    use testing_utilities

    implicit none
    integer :: stat, is
    character(len=64) :: only_suite
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'
    integer :: major, minor, patch
    character(len=:), allocatable :: version_string

    call get_testdrive_version(major, minor, patch, version_string)
    print *, "Test-Drive version: ", version_string

    ! Collect path to HAMS executable from command line arg. Needed for integration tests
    call get_hams_path_from_command_line_args()

    ! Specify the integration tests that need to be run
    testsuites = [ &
        new_testsuite("setup_config1", setup_config1), &
        new_testsuite("test_config1", collect_tests_config1), &
        new_testsuite("setup_gmres",   setup_gmres),   &
        new_testsuite("test_gmres",    collect_tests_gmres),  &
        new_testsuite("setup_image",   setup_image),   &
        new_testsuite("test_image",    collect_tests_image),  &
        new_testsuite("setup_gen_modes", setup_gen_modes), &
        new_testsuite("test_gen_modes", collect_tests_gen_modes), &
        new_testsuite("setup_drift",   setup_drift),   &
        new_testsuite("test_drift",    collect_tests_drift),  &
        new_testsuite("setup_drift_motion", setup_drift_motion), &
        new_testsuite("test_drift_motion", collect_tests_drift_motion), &
        new_testsuite("setup_array",   setup_array),   &
        new_testsuite("test_array",    collect_tests_array),  &
        new_testsuite("setup_dipole",  setup_dipole),  &
        new_testsuite("test_dipole",   collect_tests_dipole)  &
    ]

    ! An optional second argument (after the executable) runs only the suites whose name ends with it, for example
    ! "dipole" for setup_dipole and test_dipole (while iterating); without it every suite runs.
    only_suite = ''
    if (command_argument_count() >= 2) call get_command_argument(2, only_suite)
    stat = 0
    do is = 1, size(testsuites)
        if (len_trim(only_suite) > 0) then
            if (index(testsuites(is)%name, trim(only_suite), back=.true.) == 0 .or. &
                index(testsuites(is)%name, trim(only_suite), back=.true.) + len_trim(only_suite) - 1 &
                /= len_trim(testsuites(is)%name)) cycle
        end if
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        call run_testsuite(testsuites(is)%collect, error_unit, stat)
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

end program integration_tester