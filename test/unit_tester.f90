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
    use test_generalized_modes, only : collect_generalized_modes

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
        new_testsuite("generalized_modes", collect_generalized_modes) &
    ]

    stat = 0

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        ! These suites share module state (HMAT_READY, the preconditioner, the mesh arrays), so they can't run in parallel.
        call run_testsuite(testsuites(is)%collect, error_unit, stat, &
                           parallel=(testsuites(is)%name /= "iterative_solvers" .and. testsuites(is)%name /= "generalized_modes"))
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

end program unit_tester