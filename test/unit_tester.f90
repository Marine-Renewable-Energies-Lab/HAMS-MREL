program unit_tester
  
    ! Testing framework imports
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite, new_testsuite, testsuite_type
    use testdrive_version, only : get_testdrive_version

    ! Tests imports (add only collect functions)
    use test_suite2, only : collect_suite2
    use test_fingreen3d, only : collect_fingreen3d
    use test_infgreen3d, only : collect_infgreen3d
    use test_iterative_solvers, only : collect_iterative_solvers
    use test_create_directory, only : collect_create_directory
    use test_second_order_input, only : collect_second_order_input
    use test_first_order_state, only : collect_first_order_state, collect_first_order_state_serial
    use test_surface_geometry, only : collect_surface_geometry
    use test_surface_velocity, only : collect_surface_velocity
    use test_far_field, only : collect_far_field
    use test_second_order_output, only : collect_second_order_output

    implicit none
    integer :: stat, is
    type(testsuite_type), allocatable :: testsuites(:)
    character(len=*), parameter :: fmt = '("#", *(1x, a))'
    integer :: major, minor, patch
    character(len=:), allocatable :: version_string

    call get_testdrive_version(major, minor, patch, version_string)
    print *, "Test-Drive version: ", version_string

    ! Specify the tests that need to be run
    testsuites = [ &
        new_testsuite("suite2", collect_suite2), &
        new_testsuite("fingreen3d", collect_fingreen3d), &
        new_testsuite("infgreen3d", collect_infgreen3d), &
        new_testsuite("iterative_solvers", collect_iterative_solvers), &
        new_testsuite("create_directory", collect_create_directory), &
        new_testsuite("second_order_input", collect_second_order_input), &
        new_testsuite("first_order_state", collect_first_order_state), &
        new_testsuite("first_order_state_serial", collect_first_order_state_serial), &
        new_testsuite("surface_geometry", collect_surface_geometry), &
        new_testsuite("surface_velocity", collect_surface_velocity), &
        new_testsuite("far_field", collect_far_field), &
        new_testsuite("second_order_output", collect_second_order_output) &
    ]

    stat = 0

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        ! The solver tests share module state (HMAT_READY, the preconditioner), and the VINP/DINP cross-check sets
        ! WaveDyn_mod globals, so those suites can't run in parallel.
        call run_testsuite(testsuites(is)%collect, error_unit, stat, &
                           parallel=(testsuites(is)%name /= "iterative_solvers" .and. &
                                  testsuites(is)%name /= "first_order_state_serial"))
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

end program unit_tester