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
    use test_image_term, only : collect_image_term
    use test_generalized_modes, only : collect_generalized_modes
    use test_panel_integrals, only : collect_panel_integrals
    use test_second_order_input, only : collect_second_order_input
    use test_first_order_state, only : collect_first_order_state, collect_first_order_state_serial
    use test_surface_geometry, only : collect_surface_geometry
    use test_surface_velocity, only : collect_surface_velocity
    use test_far_field, only : collect_far_field
    use test_near_field, only : collect_near_field
    use test_second_order_output, only : collect_second_order_output
    use test_waterline_field, only : collect_waterline_field
    use test_source_velocity, only : collect_source_velocity
    use test_control_surface_mesh, only : collect_control_surface_mesh
    use test_field_velocity, only : collect_field_velocity
    use test_control_drift, only : collect_control_drift
    use test_qtf, only : collect_qtf
    use test_second_order_incident, only : collect_so_incident
    use test_qtf_potential, only : collect_qtf_potential
    use test_sharp_edge, only : collect_sharp_edge
    use test_qtf_field, only : collect_qtf_field
    use test_qtf_direct, only : collect_qtf_direct

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
        new_testsuite("image_term", collect_image_term), &
        new_testsuite("generalized_modes", collect_generalized_modes), &
        new_testsuite("panel_integrals", collect_panel_integrals), &
        new_testsuite("second_order_input", collect_second_order_input), &
        new_testsuite("first_order_state", collect_first_order_state), &
        new_testsuite("first_order_state_serial", collect_first_order_state_serial), &
        new_testsuite("surface_geometry", collect_surface_geometry), &
        new_testsuite("surface_velocity", collect_surface_velocity), &
        new_testsuite("far_field", collect_far_field), &
        new_testsuite("near_field", collect_near_field), &
        new_testsuite("second_order_output", collect_second_order_output), &
        new_testsuite("waterline_field", collect_waterline_field), &
        new_testsuite("source_velocity", collect_source_velocity), &
        new_testsuite("control_surface_mesh", collect_control_surface_mesh), &
        new_testsuite("field_velocity", collect_field_velocity), &
        new_testsuite("control_drift", collect_control_drift), &
        new_testsuite("qtf", collect_qtf), &
        new_testsuite("so_incident", collect_so_incident), &
        new_testsuite("qtf_potential", collect_qtf_potential), &
        new_testsuite("sharp_edge", collect_sharp_edge), &
        new_testsuite("qtf_field", collect_qtf_field), &
        new_testsuite("qtf_direct", collect_qtf_direct) &
    ]

    stat = 0

    do is = 1, size(testsuites)
        write(error_unit, fmt) "Testing:", testsuites(is)%name
        ! These suites share module state (HMAT_READY, the preconditioner, the mesh arrays, WaveDyn_mod globals, the QTF store),
        ! so they can't run in parallel.
        call run_testsuite(testsuites(is)%collect, error_unit, stat, &
                           parallel=(testsuites(is)%name /= "iterative_solvers" .and. testsuites(is)%name /= "image_term" &
                                     .and. testsuites(is)%name /= "generalized_modes" &
                                     .and. testsuites(is)%name /= "panel_integrals" &
                                     .and. testsuites(is)%name /= "first_order_state_serial" &
                                     .and. testsuites(is)%name /= "waterline_field" &
                                     .and. testsuites(is)%name /= "source_velocity" &
                                     .and. testsuites(is)%name /= "control_surface_mesh" &
                                     .and. testsuites(is)%name /= "field_velocity" &
                                     .and. testsuites(is)%name /= "qtf" &
                                     .and. testsuites(is)%name /= "qtf_potential" &
                                     .and. testsuites(is)%name /= "qtf_direct"))
    end do

    if (stat > 0) then
        write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
        error stop
    end if

end program unit_tester