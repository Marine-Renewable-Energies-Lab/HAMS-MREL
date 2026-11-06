!> Whole-solver checks of the damping lid: the stops, that a lid switched off changes nothing, and the single-body solve with a
!> lid. The cases, the measured values and the limits set from them are in test/data/damping_lid/README.md.
module test_integration_lid

    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
    use testing_utilities

    implicit none
    private
    public :: setup_lid, collect_tests_lid

    real(8), parameter :: g = 9.80665d0, omega = 1.35d0, depth = 50d0

    ! Limits set from the measurements in test/data/damping_lid/README.md (epsilon lid, 0.2 m/s unless the case says otherwise;
    ! they replace the limits of the phi = 0 lid)
    real(8), parameter :: damping_floor = 1d-6        ! B_ii >= -damping_floor * max(abs(B)); measured: every B_ii > 0
    real(8), parameter :: symmetry_limit = 1d-2       ! abs(A_ij - A_ji) / max(abs(A)), the same for B; measured 6.1e-3
    real(8), parameter :: exact_limit = 1d-12         ! cases equal by construction, printed values, relative to the largest entry
    real(8), parameter :: haskind_factor = 2d0        ! Haskind residual of lid_ann_160_eps0 below 2 times that of no lid
    real(8), parameter :: split_floor = 1d-3          ! lid_ann_160_split differs from lid_ann_160 in B11 by more than this; measured 7.9e-2
    real(8), parameter :: order_limit = 1d-11         ! relative to max(abs(A)), max(abs(B)), max(abs(X)); measured 4.7e-12
    real(8), parameter :: tri_limit = 3d-4            ! measured 1.4e-4
    real(8), parameter :: limits_limit = 1d-10        ! measured 0 at both limits
    real(8), parameter :: rigid_limit = 2d-5          ! diffraction elevation at the lid centroids of lid_ann_160_rigid over the incident
                                                      ! amplitude; measured 7.3e-6
    real(8), parameter :: off_centroid_limit = 0.25d0 ! relative difference of the elevation at node, edge, inside from that at c1;
                                                      ! measured 0.116
    real(8), parameter :: half_limit = 3d-6           ! A, B, X of half against full, relative to the largest entry; measured 1.4e-6
    real(8), parameter :: half_field_limit = 1.5d-6   ! diffraction and radiation values at the same 18 points; measured 6.7e-7
    real(8), parameter :: irsp_limit = 7d-3           ! A11, A33, B11, B33, abs(X1) with IRSP against without; measured 3.2e-3
    real(8), parameter :: collocation_limit = 1d-6    ! force from the field pressures against ExcForce.3; measured 1.6e-7, printed to 7 digits
    real(8), parameter :: maccamy_fuchs_limit = 2d-2  ! worst of X1 and two elevations of mm_cyl_nolid, complex; measured 8.5e-3
    real(8), parameter :: mm_reference_limit = 1.2d-2 ! worst of X1 and two elevations of mm_cyl_lid_c, _c_eps1, _c_rigid against
                                                      ! the reference, complex; measured 8.49e-3, 8.42e-3, 7.07e-3. HAMS-MREL is
                                                      ! 5.38e-2 from the conjugate reference (Im(alpha) flipped) at 1 m/s, and a
                                                      ! flipped binary gives 5.69e-2
    real(8), parameter :: mm_effect_limit = 3d-2      ! worst lid-effect error of the same, relative to the reference lid effect;
                                                      ! measured 2.21e-2, 2.20e-2, 1.97e-2 (flipped binary: 2.01 and 2.04 at 0.2
                                                      ! and 1 m/s)
    real(8), parameter :: mm_power_limit = 8d-3       ! lid power of mm_cyl_lid_c and _c_eps1 from the centroid elevations against
                                                      ! the reference; measured -4.14e-3, -4.12e-3 (conjugate reference +2.08e-2,
                                                      ! +1.06e-1; flipped binary +1.66e-2, +1.01e-1)
    real(8), parameter :: two_far_limit = 1d-3        ! body 1 of lid_two_far against lid_ann_160, the limit of multi_matches_single;
                                                      ! measured 4.5e-4
    real(8), parameter :: gap_symmetry_limit = 1d-2   ! abs(A_ij - A_ji) / max(abs(A)), the same for B, 12 x 12, lid_gap; measured
                                                      ! 8.0e-3, 8.6e-3 without the lid
    real(8), parameter :: gap_damping_floor = 1d-5    ! eigenvalues of B of lid_gap >= -gap_damping_floor * max(abs(B)); measured
                                                      ! 4.6e-6, 2.3e-5 without the lid (phi = 0 lid measurement), falling with
                                                      ! the hull mesh (README)

    !> The bottom-mounted cylinder of the reference cases (mm_cyl_*): radius 2, depth 10, k = 0.15 given as a wave number, so
    !> omega = sqrt(g k tanh(kh)), heading 0, field points (0, 8, 0) and (8, 0, 0).
    real(8), parameter :: mm_omega = 1.1538937638d0
    !> Reference values in the convention of ExcForce.3 and PressureElevation.6p (over rho g A and over A), from
    !> pre_postprocessing_python/DampingLid/annular_lid_reference.py: MacCamy-Fuchs for the bare cylinder, then the lid r 2 to 4
    !> at epsilon 0.2, 1 and 1e6 m/s (lid projection, extrapolated in 1/N from N = 160, 320; spread 1.4e-7 or less at 0.2 and
    !> 1 m/s, 7.7e-6 at 1e6 m/s). Each is X1, the elevation at (0, 8, 0), the elevation at (8, 0, 0).
    complex(8), parameter :: maccamy_fuchs(3) = [(1.65550518d0, 23.42095485d0), &    ! X1
                                                 (1.01387576d0, 0.04392707d0), &     ! elevation at (0, 8, 0)
                                                 (0.29600397d0, -0.97000684d0)]      ! elevation at (8, 0, 0)
    character(len=20), parameter :: mm_lid_cases(3) = [character(len=20) :: "mm_cyl_lid_c", "mm_cyl_lid_c_eps1", &
                                                       "mm_cyl_lid_c_rigid"]
    real(8), parameter :: mm_eps(3) = [0.2d0, 1d0, 1d6]
    complex(8), parameter :: mm_reference(3,3) = reshape([ &
        (1.73472971d0, 23.39194867d0), (1.00952678d0, 0.04447564d0), (0.29053829d0, -0.96832198d0), &    ! 0.2 m/s
        (2.02520747d0, 23.23109427d0), (0.99339102d0, 0.04910955d0), (0.27074161d0, -0.95859350d0), &    ! 1 m/s
        (0.87716703d0, 20.84072211d0), (0.99990858d0, 0.17835860d0), (0.31664285d0, -0.79760753d0)], &  ! 1e6 m/s
        [3, 3])
    !> Lid power of the reference, 1/2 rho E int abs(dphi/dz)**2 dA in W per m**2 of incident amplitude, rho = 1000, at 0.2 and
    !> 1 m/s (extrapolated from N = 160, 320; spread 5.9e-8 and 1.4e-6). It equals the net power into the cylinder r = 1000 to
    !> 1.1e-11 and 2.6e-10 (N = 320).
    real(8), parameter :: mm_lid_power(2) = [5.268414d3, 2.472964d4]
    !> Lid panels of the 48 x 8 cases; their centroids are field points 3 to 386, after (0, 8, 0) and (8, 0, 0).
    integer, parameter :: n384 = 384

    !> The single-body solve cases, all run by setup_lid: those with a lid and their no-lid references lid_limits_off and
    !> lid_nolid_160.
    character(len=24), parameter :: solve_cases(22) = [character(len=24) :: "lid_ann_40", "lid_ann_160", "lid_ann_640", &
        "lid_ann_160_rev", "lid_ann_tri", "lid_touch", "lid_patches", "lid_shrink_1", "lid_shrink_2", "lid_deep", "lid_limits", &
        "lid_limits_off", "lid_hull_160", "lid_field_off_centroid", &
        "lid_ann_160_pts", "lid_ann_160_half", "lid_ann_160_irsp", &
        "lid_ann_160_eps0", "lid_ann_160_rigid", "lid_ann_160_two_dampers", "lid_ann_160_split", "lid_nolid_160"]

    !> The two-body cases, run by setup_lid: with a lid, and lid_gap_nolid without.
    character(len=16), parameter :: two_cases(5) = [character(len=16) :: "lid_two_far", "lid_gap", "lid_gap_eps0", &
        "lid_gap_eps1", "lid_gap_nolid"]

    !> Field points of lid_ann_160 and its epsilon variants: the 160 lid centroids, then (12, 0, 0), (20, 0, 0) and (12, 0, -3).
    integer, parameter :: n160 = 160

    !> A, B and X at one frequency: the full 6 x 6 matrices and the six forces of heading 0.
    type :: coefficients
        real(8) :: a(6,6), b(6,6)
        complex(8) :: x(6)
    end type coefficients

contains

    subroutine setup_lid(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("run_lid_flag_off", run_lid_flag_off), &
            new_unittest("run_lid_reference", run_lid_reference), &
            new_unittest("run_lid_solves", run_lid_solves), &
            new_unittest("run_lid_mm", run_lid_mm), &
            new_unittest("run_lid_two", run_lid_two) &
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
            new_unittest("stops_lid_value", stops_lid_value), &
            new_unittest("stops_lid_repeated", stops_lid_repeated), &
            new_unittest("stops_lid_damp_missing", stops_lid_damp_missing), &
            new_unittest("stops_lid_damp_count", stops_lid_damp_count), &
            new_unittest("stops_lid_damp_range", stops_lid_damp_range), &
            new_unittest("stops_lid_damp_order", stops_lid_damp_order), &
            new_unittest("stops_lid_damp_overlap", stops_lid_damp_overlap), &
            new_unittest("stops_lid_damp_gap", stops_lid_damp_gap), &
            new_unittest("stops_lid_damp_negative", stops_lid_damp_negative), &
            new_unittest("stops_lid_damp_value", stops_lid_damp_value), &
            new_unittest("stops_lid_zero_frequency", stops_lid_zero_frequency), &
            new_unittest("flag_off_matches_no_lid", flag_off_matches_no_lid), &
            new_unittest("lid_runs", lid_runs), &
            new_unittest("lid_damping_nonnegative", lid_damping_nonnegative), &
            new_unittest("lid_converges", lid_converges), &
            new_unittest("lid_symmetric_ab", lid_symmetric_ab), &
            new_unittest("lid_haskind_eps0", lid_haskind_eps0), &
            new_unittest("lid_order_free", lid_order_free), &
            new_unittest("lid_tri_close", lid_tri_close), &
            new_unittest("lid_shrink_to_none", lid_shrink_to_none), &
            new_unittest("lid_limits_finite", lid_limits_finite), &
            new_unittest("lid_eps0_matches_no_lid", lid_eps0_matches_no_lid), &
            new_unittest("lid_rigid_elevation_small", lid_rigid_elevation_small), &
            new_unittest("lid_two_dampers_split", lid_two_dampers_split), &
            new_unittest("lid_absorbs", lid_absorbs), &
            new_unittest("lid_field_off_centroid", lid_field_off_centroid), &
            new_unittest("lid_field_matches_collocation", lid_field_matches_collocation), &
            new_unittest("lid_half_matches_full", lid_half_matches_full), &
            new_unittest("lid_irsp_matches", lid_irsp_matches), &
            new_unittest("mm_nolid_matches_maccamy_fuchs", mm_nolid_matches_maccamy_fuchs), &
            new_unittest("mm_lid_matches_reference", mm_lid_matches_reference), &
            new_unittest("mm_lid_converges", mm_lid_converges), &
            new_unittest("mm_lid_power_matches", mm_lid_power_matches), &
            new_unittest("stops_lid_inside_two", stops_lid_inside_two), &
            new_unittest("stops_lid_irsp_two", stops_lid_irsp_two), &
            new_unittest("lid_multi_matches_single", lid_multi_matches_single), &
            new_unittest("lid_gap_symmetric_ab", lid_gap_symmetric_ab), &
            new_unittest("lid_gap_damping_nonnegative", lid_gap_damping_nonnegative), &
            new_unittest("lid_gap_eps0_matches_no_lid", lid_gap_eps0_matches_no_lid), &
            new_unittest("lid_gap_damps", lid_gap_damps) &
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

    !> Every solve case, whether or not an earlier one failed, so that every test below has its output.
    subroutine run_lid_solves(error)
        type(error_type), allocatable, intent(out) :: error
        type(error_type), allocatable :: one
        character(len=512) :: failed
        integer :: k

        failed = ""
        do k = 1, size(solve_cases)
            call run_case("damping_lid", trim(solve_cases(k)), one)
            if (allocated(one)) then
                failed = trim(failed)//" "//trim(solve_cases(k))
                deallocate (one)
            end if
        end do
        call check(error, len_trim(failed) == 0, "HAMS failed for"//trim(failed))
    end subroutine run_lid_solves

    !> The bottom-mounted cylinder of the independent reference: without lid, with the lid at three meshes (0.2 m/s), and at 1
    !> and 1e6 m/s on the finest.
    subroutine run_lid_mm(error)
        type(error_type), allocatable, intent(out) :: error
        type(error_type), allocatable :: one
        character(len=20), parameter :: names(6) = [character(len=20) :: "mm_cyl_nolid", "mm_cyl_lid_a", "mm_cyl_lid_b", &
                                                    "mm_cyl_lid_c", "mm_cyl_lid_c_eps1", "mm_cyl_lid_c_rigid"]
        character(len=128) :: failed
        integer :: k

        failed = ""
        do k = 1, size(names)
            call run_case("damping_lid", trim(names(k)), one)
            if (allocated(one)) then
                failed = trim(failed)//" "//trim(names(k))
                deallocate (one)
            end if
        end do
        call check(error, len_trim(failed) == 0, "HAMS failed for"//trim(failed))
    end subroutine run_lid_mm

    !> The two-body cases with a lid, whether or not the first failed.
    subroutine run_lid_two(error)
        type(error_type), allocatable, intent(out) :: error
        type(error_type), allocatable :: one
        character(len=128) :: failed
        integer :: k

        failed = ""
        do k = 1, size(two_cases)
            call run_case("damping_lid", trim(two_cases(k)), one)
            if (allocated(one)) then
                failed = trim(failed)//" "//trim(two_cases(k))
                deallocate (one)
            end if
        end do
        call check(error, len_trim(failed) == 0, "HAMS failed for"//trim(failed))
    end subroutine run_lid_two

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

    !> Damping_lid 2: only 0 and 1 are allowed.
    subroutine stops_lid_value(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_value", "Damping_lid must be 0 or 1", error)
    end subroutine stops_lid_value

    !> Panel 1 of the lid has a repeated vertex.
    subroutine stops_lid_repeated(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_repeated", "lid panel in DampingLid.pnl has a repeated vertex", error)
    end subroutine stops_lid_repeated

    subroutine stops_lid_damp_missing(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_damp_missing", "no damping block", error)
    end subroutine stops_lid_damp_missing

    subroutine stops_lid_damp_count(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_damp_count", "number of dampers", error)
    end subroutine stops_lid_damp_count

    subroutine stops_lid_damp_range(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_damp_range", "damper panel range", error)
    end subroutine stops_lid_damp_range

    subroutine stops_lid_damp_order(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_damp_order", "damper panel range", error)
    end subroutine stops_lid_damp_order

    subroutine stops_lid_damp_overlap(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_damp_overlap", "lid panel in two dampers", error)
    end subroutine stops_lid_damp_overlap

    subroutine stops_lid_damp_gap(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_damp_gap", "lid panel in no damper", error)
    end subroutine stops_lid_damp_gap

    subroutine stops_lid_damp_negative(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_damp_negative", "negative damping", error)
    end subroutine stops_lid_damp_negative

    subroutine stops_lid_damp_value(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_damp_value", "cannot read the damping", error)
    end subroutine stops_lid_damp_value

    !> The frequency list holds 0 and 1.35 rad/s; the lid has no equation at zero frequency.
    subroutine stops_lid_zero_frequency(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_zero_frequency", "zero frequency with a damping lid", error)
    end subroutine stops_lid_zero_frequency

    !> A lid disc inside the second of two bodies: the multi-body branch of HULL_SOLID_ANGLE.
    subroutine stops_lid_inside_two(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_inside_two", "lid panel inside a hull", error)
    end subroutine stops_lid_inside_two

    !> The multi-body lid has no interior rows (README).
    subroutine stops_lid_irsp_two(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case_expect_stop("lid_stop_irsp_two", "damping lid with several bodies needs If_remove_irr_freq 0", error)
    end subroutine stops_lid_irsp_two

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

    ! -------------------------------------------------------------------------------------------- solve: readers

    !> A and B of AmssDamp.1 at frequency freq, and X of heading 0 of ExcForce.3 at the same frequency. The zero and infinite
    !> frequency rows have no B, so they have four columns. X is left zero when freq has no excitation rows.
    subroutine read_coefficients(name, freq, c, error)
        character(len=*), intent(in) :: name
        real(8), intent(in) :: freq
        type(coefficients), intent(out) :: c
        type(error_type), allocatable, intent(out) :: error
        character(len=512) :: line
        real(8) :: f, beta, row(7)
        integer :: unit, stat, i, j, n

        c%a = 0d0
        c%b = 0d0
        c%x = (0d0, 0d0)
        n = 0
        open (newunit=unit, file=output_dir(name)//"Wamit_format/AmssDamp.1", status="old", action="read", iostat=stat)
        if (stat /= 0) then
            call test_failed(error, "No AmssDamp.1 for "//name)
            return
        end if
        do
            read (unit, '(a)', iostat=stat) line
            if (stat /= 0) exit
            if (columns(line) == 5) then
                read (line, *) f, i, j, row(1:2)
            else if (columns(line) == 4) then
                read (line, *) f, i, j, row(1)
                row(2) = 0d0
            else
                cycle
            end if
            if (abs(f - freq) > 1d-6 .or. i > 6 .or. j > 6) cycle
            c%a(i, j) = row(1)
            c%b(i, j) = row(2)
            n = n + 1
        end do
        close (unit)
        call check(error, n == 36, name//": AmssDamp.1 has not 36 entries at the frequency")
        if (allocated(error)) return
        open (newunit=unit, file=output_dir(name)//"Wamit_format/ExcForce.3", status="old", action="read", iostat=stat)
        if (stat /= 0) then
            call test_failed(error, "No ExcForce.3 for "//name)
            return
        end if
        do
            read (unit, '(a)', iostat=stat) line
            if (stat /= 0) exit
            if (columns(line) /= 7) cycle
            read (line, *) f, beta, i, row(1:4)
            if (abs(f - freq) < 1d-6 .and. abs(beta) < 1d-6 .and. i <= 6) c%x(i) = cmplx(row(3), row(4), kind=8)
        end do
        close (unit)
    end subroutine read_coefficients

    integer function columns(line)
        character(len=*), intent(in) :: line
        logical :: in_word
        integer :: i

        columns = 0
        in_word = .false.
        do i = 1, len_trim(line)
            if (line(i:i) /= " " .and. .not. in_word) columns = columns + 1
            in_word = line(i:i) /= " "
        end do
    end function columns

    !> The five values the convergence and closeness checks follow: A11, A33, B11, B33, abs(X1).
    function key_values(c) result(v)
        type(coefficients), intent(in) :: c
        real(8) :: v(5)
        v = [c%a(1, 1), c%a(3, 3), c%b(1, 1), c%b(3, 3), abs(c%x(1))]
    end function key_values

    subroutine check_below(error, name, value, limit)
        type(error_type), allocatable, intent(out) :: error
        character(len=*), intent(in) :: name
        real(8), intent(in) :: value, limit
        character(len=64) :: text

        write (*, '(4x, a, ": ", es10.3)') name, value
        write (text, '(es10.3, " is not below ", es10.3)') value, limit
        call check(error, value < limit, name//": "//trim(text))
    end subroutine check_below

    ! -------------------------------------------------------------------------------------------- solve: tests

    !> Every solve case wrote AmssDamp.1 and ExcForce.3 with no NaN or Infinity in them. The exit codes are checked by
    !> run_lid_solves.
    subroutine lid_runs(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=10), parameter :: files(2) = [character(len=10) :: "AmssDamp.1", "ExcForce.3"]
        character(len=512) :: line
        integer :: k, m, unit, stat, n

        do k = 1, size(solve_cases)
            do m = 1, size(files)
                open (newunit=unit, file=output_dir(trim(solve_cases(k)))//"Wamit_format/"//trim(files(m)), status="old", &
                      action="read", iostat=stat)
                if (stat /= 0) then
                    call test_failed(error, "No "//trim(files(m))//" for "//trim(solve_cases(k)))
                    return
                end if
                n = 0
                do
                    read (unit, '(a)', iostat=stat) line
                    if (stat /= 0) exit
                    n = n + 1
                    if (index(line, "NaN") > 0 .or. index(line, "Infinity") > 0) then
                        call test_failed(error, trim(files(m))//" of "//trim(solve_cases(k))//" holds "//trim(line))
                        close (unit)
                        return
                    end if
                end do
                close (unit)
                call check(error, n > 0, trim(files(m))//" of "//trim(solve_cases(k))//" is empty")
                if (allocated(error)) return
            end do
        end do
    end subroutine lid_runs

    !> B_ii >= 0: the body radiates and the lid absorbs, both non-negative, to the round-off of the printed values.
    subroutine lid_damping_nonnegative(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: c
        real(8) :: worst, bmax
        integer :: k, m

        worst = 0d0
        do k = 1, size(solve_cases)
            call read_coefficients(trim(solve_cases(k)), omega, c, error)
            if (allocated(error)) return
            bmax = maxval(abs(c%b))
            do m = 1, 3
                worst = max(worst, -c%b(m, m)/bmax)
            end do
        end do
        call check_below(error, "lid_damping_nonnegative: worst -B_ii/max(abs(B))", worst, damping_floor)
    end subroutine lid_damping_nonnegative

    !> The change from 160 to 640 lid panels is smaller than that from 40 to 160, in each of A11, A33, B11, B33 and abs(X1)
    !> (ratios 3.1, 4.0, 8.5, 4.3, 8.6 measured with epsilon 0.2; with the phi = 0 lid B33 did not converge at 40 panels and was
    !> not checked, README).
    subroutine lid_converges(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=6), parameter :: labels(5) = [character(len=6) :: "A11", "A33", "B11", "B33", "abs X1"]
        type(coefficients) :: c40, c160, c640
        real(8) :: v40(5), v160(5), v640(5), d1, d2
        integer :: m

        call read_coefficients("lid_ann_40", omega, c40, error)
        if (allocated(error)) return
        call read_coefficients("lid_ann_160", omega, c160, error)
        if (allocated(error)) return
        call read_coefficients("lid_ann_640", omega, c640, error)
        if (allocated(error)) return
        v40 = key_values(c40)
        v160 = key_values(c160)
        v640 = key_values(c640)
        do m = 1, 5
            d1 = abs(v160(m) - v40(m))
            d2 = abs(v640(m) - v160(m))
            write (*, '(4x, a6, " at 40, 160, 640: ", 3es14.6, "  changes ", 2es10.3, "  ratio ", es10.3)') &
                labels(m), v40(m), v160(m), v640(m), d1, d2, d1/max(d2, tiny(1d0))
            call check(error, d2 < d1, "lid_converges: "//trim(labels(m))//" changes more from 160 to 640 than from 40 to 160")
            if (allocated(error)) return
        end do
    end subroutine lid_converges

    !> A and B are symmetric, by reciprocity: the lid condition is the same homogeneous condition in every problem (spec).
    subroutine lid_symmetric_ab(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=16), parameter :: names(3) = [character(len=16) :: "lid_ann_640", "lid_touch", "lid_patches"]
        type(coefficients) :: c
        real(8) :: worst_a, worst_b, worst
        integer :: k

        worst = 0d0
        do k = 1, size(names)
            call read_coefficients(trim(names(k)), omega, c, error)
            if (allocated(error)) return
            worst_a = maxval(abs(c%a - transpose(c%a)))/maxval(abs(c%a))
            worst_b = maxval(abs(c%b - transpose(c%b)))/maxval(abs(c%b))
            write (*, '(4x, a, ": A ", es10.3, ", B ", es10.3)') trim(names(k)), worst_a, worst_b
            worst = max(worst, worst_a, worst_b)
        end do
        call check_below(error, "lid_symmetric_ab", worst, symmetry_limit)
    end subroutine lid_symmetric_ab

    !> The relative residual of the Haskind relation for heave of an axisymmetric body, as haskind_heave of the image term suite.
    real(8) function haskind_residual(c)
        type(coefficients), intent(in) :: c
        real(8) :: k, vg, predicted
        integer :: i

        k = omega*omega/g
        do i = 1, 100
            k = omega*omega/(g*tanh(k*depth))
        end do
        vg = omega/(2d0*k)*(1d0 + 2d0*k*depth/sinh(2d0*k*depth))
        predicted = k*g*abs(c%x(3))**2/(4d0*omega*vg)
        haskind_residual = abs(predicted - c%b(3, 3))/c%b(3, 3)
    end function haskind_residual

    !> The Haskind relation over the hull does not hold with epsilon > 0, since the lid absorbs (spec). With epsilon = 0 the lid
    !> is no lid, so the heave residual of lid_ann_160_eps0 stays below twice that of lid_nolid_160. lid_ann_160 is printed.
    subroutine lid_haskind_eps0(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: c0, bare, c

        call read_coefficients("lid_ann_160_eps0", omega, c0, error)
        if (allocated(error)) return
        call read_coefficients("lid_nolid_160", omega, bare, error)
        if (allocated(error)) return
        call read_coefficients("lid_ann_160", omega, c, error)
        if (allocated(error)) return
        write (*, '(4x, "Haskind heave residual: no lid ", es10.3, ", epsilon 0 ", es10.3, ", epsilon 0.2 ", es10.3)') &
            haskind_residual(bare), haskind_residual(c0), haskind_residual(c)
        call check_below(error, "lid_haskind_eps0", haskind_residual(c0), haskind_factor*haskind_residual(bare))
    end subroutine lid_haskind_eps0

    !> The largest relative difference of A, B and X of two cases, each relative to the largest entry of that matrix in b.
    real(8) function coefficient_difference(a, b)
        type(coefficients), intent(in) :: a, b
        coefficient_difference = max(maxval(abs(a%a - b%a))/maxval(abs(b%a)), maxval(abs(a%b - b%b))/maxval(abs(b%b)), &
                                     maxval(abs(a%x - b%x))/maxval(abs(b%x)))
    end function coefficient_difference

    !> epsilon = 0 gives alpha = 0: the lid rows only evaluate the representation and do not feed the hull rows, so A, B and X
    !> of lid_ann_160_eps0 are those of lid_nolid_160 (the same input with Damping_lid 0), and so are the diffraction values at
    !> every field point, on the lid (whose elevation is i omega phi / g at epsilon = 0) and off it, and the radiation pressure in
    !> the fluid at (12, 0, -3). The no-lid run writes no radiation rows for points on z = 0.
    subroutine lid_eps0_matches_no_lid(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: c0, bare
        complex(8), allocatable :: dif0(:), rad0(:,:), difb(:), radb(:,:)
        real(8) :: worst_c, worst_f, worst_r

        call read_coefficients("lid_ann_160_eps0", omega, c0, error)
        if (allocated(error)) return
        call read_coefficients("lid_nolid_160", omega, bare, error)
        if (allocated(error)) return
        worst_c = coefficient_difference(c0, bare)
        call read_field("lid_ann_160_eps0", dif0, rad0, error)
        if (allocated(error)) return
        call read_field("lid_nolid_160", difb, radb, error)
        if (allocated(error)) return
        call check(error, size(dif0) == n160 + 3 .and. size(difb) == n160 + 3 .and. size(radb, 1) == n160 + 3, &
                   "lid_eps0_matches_no_lid: not 163 points")
        if (allocated(error)) return
        worst_f = maxval(abs(dif0 - difb))/maxval(abs(difb))
        worst_r = maxval(abs(rad0(n160 + 3, :) - radb(n160 + 3, :)))/maxval(abs(radb(n160 + 3, :)))
        write (*, '(4x, "lid_eps0_matches_no_lid: A, B, X ", es10.3, ", diffraction ", es10.3, ", radiation in the fluid ", &
               &es10.3)') worst_c, worst_f, worst_r
        call check_below(error, "lid_eps0_matches_no_lid", max(worst_c, worst_f, worst_r), exact_limit)
    end subroutine lid_eps0_matches_no_lid

    !> Two dampers 1 to 80 and 81 to 160, both 0.2 m/s, are one damper of 0.2 m/s: A, B, X and every field value equal. With
    !> 1.0 m/s on 81 to 160 instead (lid_ann_160_split) B11 differs from both by more than split_floor: each panel takes its own
    !> epsilon.
    subroutine lid_two_dampers_split(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: one, two, split
        complex(8), allocatable :: dif1(:), rad1(:,:), dif2(:), rad2(:,:)
        real(8) :: worst, worst_f, d1, d2

        call read_coefficients("lid_ann_160", omega, one, error)
        if (allocated(error)) return
        call read_coefficients("lid_ann_160_two_dampers", omega, two, error)
        if (allocated(error)) return
        call read_coefficients("lid_ann_160_split", omega, split, error)
        if (allocated(error)) return
        worst = coefficient_difference(two, one)
        call read_field("lid_ann_160", dif1, rad1, error)
        if (allocated(error)) return
        call read_field("lid_ann_160_two_dampers", dif2, rad2, error)
        if (allocated(error)) return
        worst_f = max(maxval(abs(dif2 - dif1))/maxval(abs(dif1)), maxval(abs(rad2 - rad1))/maxval(abs(rad1)))
        d1 = abs(split%b(1, 1) - one%b(1, 1))/abs(one%b(1, 1))
        d2 = abs(split%b(1, 1) - two%b(1, 1))/abs(two%b(1, 1))
        write (*, '(4x, "lid_two_dampers_split: two against one ", es10.3, ", field ", es10.3, "; B11 0.2: ", es14.6, &
               &", split: ", es14.6, ", relative ", es10.3)') worst, worst_f, one%b(1, 1), split%b(1, 1), d1
        call check_below(error, "lid_two_dampers_split: two equal dampers against one", max(worst, worst_f), exact_limit)
        if (allocated(error)) return
        call check(error, min(d1, d2) > split_floor, "lid_two_dampers_split: the split lid gives the B11 of the uniform one")
    end subroutine lid_two_dampers_split

    !> Areas of the first n panels of DampingLid.pnl of a case, by the cross product of the diagonals (planar quads) or of two
    !> edges (triangles), as CalPanelArea of the solver for flat panels.
    subroutine read_lid_areas(name, n, area, error)
        character(len=*), intent(in) :: name
        integer, intent(in) :: n
        real(8), intent(out) :: area(n)
        type(error_type), allocatable, intent(out) :: error
        character(len=512) :: line
        real(8), allocatable :: node(:,:)
        real(8) :: p(3,4), v(3)
        integer :: unit, stat, np, nn, isx, isy, i, k, id, nv, ids(4)

        open (newunit=unit, file=input_dir("damping_lid", name)//"DampingLid.pnl", status="old", action="read", iostat=stat)
        if (stat /= 0) then
            call test_failed(error, "No DampingLid.pnl for "//name)
            return
        end if
        read (unit, '(a)') line
        read (unit, '(a)') line
        read (unit, '(a)') line
        read (unit, *) np, nn, isx, isy
        allocate (node(3, nn))
        read (unit, '(a)') line
        read (unit, '(a)') line
        do i = 1, nn
            read (unit, *) id, node(:, i)
        end do
        read (unit, '(a)') line
        read (unit, '(a)') line
        read (unit, '(a)') line
        do i = 1, n
            read (unit, *) id, nv, ids(1:nv)
            do k = 1, nv
                p(:, k) = node(:, ids(k))
            end do
            if (nv == 3) then
                v = 0.5d0*cross(p(:, 2) - p(:, 1), p(:, 3) - p(:, 1))
            else
                v = 0.5d0*cross(p(:, 3) - p(:, 1), p(:, 4) - p(:, 2))
            end if
            area(i) = norm2(v)
        end do
        close (unit)
    end subroutine read_lid_areas

    !> The lid absorbs. lid_damping_nonnegative checks B_ii >= 0. Here, from the elevations at the 160 lid centroids of
    !> lid_ann_160 (one-point rule per panel), the lid power P = 1/2 rho omega**2 sum_l epsilon_l abs(eta_l)**2 A_l (spec):
    !> - diffraction, per unit incident amplitude: P >= 0, and P > 0 for epsilon 0.2;
    !> - radiation of surge and heave, per unit motion: the lid power is part of the power the body radiates,
    !>   1/2 omega**2 B_jj, so rho sum_l epsilon abs(eta_l)**2 A_l <= B_jj. In the nondimensional output (B over rho omega,
    !>   eta per unit motion, REFL = 1): sum_l epsilon abs(eta_l)**2 A_l <= omega B_jj. Derived; a flipped sign of alpha would make
    !>   the lid a source and could break it. The ratio is printed and checked below 1.
    subroutine lid_absorbs(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), parameter :: eps = 0.2d0, rho = 1000d0
        type(coefficients) :: c
        complex(8), allocatable :: dif(:), rad(:,:)
        real(8) :: area(n160), power, lid_rad(2), ratio(2)
        integer :: m

        call read_lid_areas("lid_ann_160", n160, area, error)
        if (allocated(error)) return
        call read_field("lid_ann_160", dif, rad, error)
        if (allocated(error)) return
        call check(error, size(dif) == n160 + 3 .and. size(rad, 1) == n160 + 3, "lid_absorbs: not 163 points")
        if (allocated(error)) return
        call read_coefficients("lid_ann_160", omega, c, error)
        if (allocated(error)) return
        power = 0.5d0*rho*omega**2*sum(eps*abs(dif(1:n160))**2*area)
        do m = 1, 2
            lid_rad(m) = sum(eps*abs(rad(1:n160, 2*m - 1))**2*area)
            ratio(m) = lid_rad(m)/(omega*c%b(2*m - 1, 2*m - 1))
        end do
        write (*, '(4x, "lid_absorbs: diffraction lid power ", es12.4, " W per m**2 of incident amplitude; lid area ", f9.4)') &
            power, sum(area)
        write (*, '(4x, "lid_absorbs: lid power over radiated power, surge ", es10.3, ", heave ", es10.3)') ratio
        call check(error, power > 0d0, "lid_absorbs: the diffraction lid power is not positive")
        if (allocated(error)) return
        call check_below(error, "lid_absorbs: largest lid power over radiated power", maxval(ratio), 1d0)
    end subroutine lid_absorbs

    !> The vertex order of the lid panels does not matter.
    subroutine lid_order_free(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: c, r
        real(8) :: worst

        call read_coefficients("lid_ann_160", omega, c, error)
        if (allocated(error)) return
        call read_coefficients("lid_ann_160_rev", omega, r, error)
        if (allocated(error)) return
        worst = max(maxval(abs(r%a - c%a))/maxval(abs(c%a)), maxval(abs(r%b - c%b))/maxval(abs(c%b)), &
                    maxval(abs(r%x - c%x))/maxval(abs(c%x)))
        call check_below(error, "lid_order_free", worst, order_limit)
    end subroutine lid_order_free

    !> Each quad of lid_ann_160 split in two triangles gives nearly the same A33, B33 and abs(X1).
    subroutine lid_tri_close(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: c, t
        real(8) :: vc(5), vt(5), worst

        call read_coefficients("lid_ann_160", omega, c, error)
        if (allocated(error)) return
        call read_coefficients("lid_ann_tri", omega, t, error)
        if (allocated(error)) return
        vc = key_values(c)
        vt = key_values(t)
        worst = maxval(abs(vt([2, 4, 5]) - vc([2, 4, 5]))/abs(vc([2, 4, 5])))
        call check_below(error, "lid_tri_close", worst, tri_limit)
    end subroutine lid_tri_close

    !> A lid ring that narrows towards nothing moves the results towards those of the run with no lid.
    subroutine lid_shrink_to_none(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: bare, s1, s2
        real(8) :: d1, d2, vb(5)

        call read_coefficients("image_cos", omega, bare, error)
        if (allocated(error)) return
        call read_coefficients("lid_shrink_1", omega, s1, error)
        if (allocated(error)) return
        call read_coefficients("lid_shrink_2", omega, s2, error)
        if (allocated(error)) return
        vb = key_values(bare)
        d1 = maxval(abs(key_values(s1) - vb)/abs(vb))
        d2 = maxval(abs(key_values(s2) - vb)/abs(vb))
        write (*, '(4x, "lid_shrink_to_none: distance from no lid at width 0.5, 0.25: ", 2es10.3)') d1, d2
        call check(error, d2 < d1, "lid_shrink_to_none: the narrower ring is not closer to no lid")
    end subroutine lid_shrink_to_none

    !> With the 0 and infinite frequency limits on: no NaN (lid_runs), and the lid is skipped at both limits (alpha = 0 at zero
    !> frequency, phi = 0 on the free surface at infinite frequency), so A there is that of the same run with Damping_lid 0.
    subroutine lid_limits_finite(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: lid, bare
        real(8) :: d(2)
        real(8), parameter :: freqs(2) = [0d0, -1d0]
        integer :: k

        do k = 1, 2
            call read_coefficients("lid_limits", freqs(k), lid, error)
            if (allocated(error)) return
            call read_coefficients("lid_limits_off", freqs(k), bare, error)
            if (allocated(error)) return
            d(k) = maxval(abs(lid%a - bare%a))/maxval(abs(bare%a))
        end do
        write (*, '(4x, "lid_limits_finite: A against no lid at zero ", es10.3, ", at infinite frequency ", es10.3)') d
        call check_below(error, "lid_limits_finite: A at both limits against no lid", maxval(d), limits_limit)
    end subroutine lid_limits_finite

    ! -------------------------------------------------------------------------------------------- field output

    !> The rows of PressureElevation.6p of a lid case: dif(ip) the diffraction value of point ip (7 columns, value in the last
    !> two), rad(ip, m) the radiation value of mode m (14 columns, mode m in columns 2m + 1 and 2m + 2). All non-dimensional.
    !> The radiation rows of points on z = 0 are written when the lid is on.
    subroutine read_field(name, dif, rad, error)
        character(len=*), intent(in) :: name
        complex(8), allocatable, intent(out) :: dif(:), rad(:,:)
        type(error_type), allocatable, intent(out) :: error
        character(len=512) :: line
        real(8) :: row(14)
        complex(8), allocatable :: bigger(:,:)
        integer :: unit, stat, n, ip, m

        allocate (dif(0), rad(0, 6))
        open (newunit=unit, file=output_dir(name)//"Wamit_format/PressureElevation.6p", status="old", action="read", iostat=stat)
        if (stat /= 0) then
            call test_failed(error, "No PressureElevation.6p for "//name)
            return
        end if
        do
            read (unit, '(a)', iostat=stat) line
            if (stat /= 0) exit
            n = columns(line)
            if (n == 7) then
                read (line, *) row(1:7)
                ip = nint(row(3))
                if (ip > size(dif)) dif = [dif, (cmplx(0d0, 0d0, kind=8), m = size(dif) + 1, ip)]
                dif(ip) = cmplx(row(6), row(7), kind=8)
            else if (n == 14) then
                read (line, *) row(1:14)
                ip = nint(row(2))
                if (ip > size(rad, 1)) then
                    allocate (bigger(ip, 6))
                    bigger = (0d0, 0d0)
                    bigger(1:size(rad, 1), :) = rad
                    call move_alloc(bigger, rad)
                end if
                do m = 1, 6
                    rad(ip, m) = cmplx(row(2*m + 1), row(2*m + 2), kind=8)
                end do
            end if
        end do
        close (unit)
    end subroutine read_field

    !> epsilon = 1e6 m/s is close to a rigid lid: alpha tends to -nu and eta = i omega phi / (g - i omega epsilon) to 0. The total
    !> diffraction elevation at the 160 lid centroids of lid_ann_160_rigid, over the incident amplitude (1 in the non-dimensional
    !> output), is small, while the wave off the lid, at 12 and 20 m, is not.
    subroutine lid_rigid_elevation_small(error)
        type(error_type), allocatable, intent(out) :: error
        complex(8), allocatable :: dif(:), rad(:,:)

        call read_field("lid_ann_160_rigid", dif, rad, error)
        if (allocated(error)) return
        call check(error, size(dif) == n160 + 3, "lid_rigid_elevation_small: not 163 points")
        if (allocated(error)) return
        write (*, '(4x, "lid_rigid_elevation_small: |elevation| at 12 and 20 m: ", 2es10.3, ", largest radiation on the lid ", &
               &es10.3)') abs(dif(n160 + 1:n160 + 2)), maxval(abs(rad(1:n160, :)))
        call check(error, maxval(abs(dif(n160 + 1:n160 + 2))) > 0.1d0, "lid_rigid_elevation_small: no wave outside the lid")
        if (allocated(error)) return
        call check_below(error, "lid_rigid_elevation_small: diffraction", maxval(abs(dif(1:n160))), rigid_limit)
    end subroutine lid_rigid_elevation_small

    !> At a lid node, in the middle of a lid edge, inside a panel off its centroid, at a node of the radius 8 ring (0, 8, 0) and
    !> far away, the field is finite. The first three lie in or on the edge of panel 1 or next to it, so their diffraction
    !> elevation is close to that at the centroid c1 of panel 1: the relative difference is below off_centroid_limit.
    subroutine lid_field_off_centroid(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=10), parameter :: labels(15) = [character(len=10) :: "node", "edge", "inside", "c1", "c2", "c3", "c4", &
            "c5", "c6", "c7", "c8", "c9", "c10", "(0, 8, 0)", "(20, 0, 0)"]
        complex(8), allocatable :: dif(:), rad(:,:)
        real(8) :: worst
        integer :: ip

        call read_field("lid_field_off_centroid", dif, rad, error)
        if (allocated(error)) return
        call check(error, size(dif) == 15 .and. size(rad, 1) == 15, "lid_field_off_centroid: not 15 points in each")
        if (allocated(error)) return
        worst = 0d0
        do ip = 1, 15
            write (*, '(4x, a10, " diffraction ", 2es13.5, "  |radiation| ", 6es10.3)') labels(ip), dif(ip), abs(rad(ip, :))
            if (dif(ip) /= dif(ip) .or. any(rad(ip, :) /= rad(ip, :))) then
                call test_failed(error, "lid_field_off_centroid: NaN at "//trim(labels(ip)))
                return
            end if
            if (ip <= 3) worst = max(worst, abs(dif(ip) - dif(4))/abs(dif(4)))
        end do
        call check_below(error, "lid_field_off_centroid: elevation at node, edge, inside against c1", worst, off_centroid_limit)
    end subroutine lid_field_off_centroid

    !> As field_points_match_collocation of the image term suite: the pressure at every hull centroid, integrated over the hull,
    !> is the excitation force, here with a lid.
    subroutine lid_field_matches_collocation(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: xc(:,:), area(:), normal(:,:)
        complex(8), allocatable :: dif(:), rad(:,:)
        type(coefficients) :: c
        complex(8) :: force(3)
        real(8) :: worst
        integer :: i, k

        call read_panels(xc, area, normal, error)
        if (allocated(error)) return
        call read_field("lid_hull_160", dif, rad, error)
        if (allocated(error)) return
        call check(error, size(dif) == size(area), "lid_field_matches_collocation: not every panel has a pressure")
        if (allocated(error)) return
        call read_coefficients("lid_hull_160", omega, c, error)
        if (allocated(error)) return
        force = (0d0, 0d0)
        do i = 1, size(area)
            do k = 1, 3
                force(k) = force(k) + dif(i)*area(i)*normal(k, i)
            end do
        end do
        worst = 0d0
        do k = 1, 3, 2
            ! The sign of the force follows the normal of the solver, which points into the body
            worst = max(worst, abs(abs(force(k) + c%x(k)))/abs(c%x(k)))
        end do
        call check_below(error, "lid_field_matches_collocation: force from field points", worst, collocation_limit)
    end subroutine lid_field_matches_collocation

    !> The half hull with the half lid (ISY = 1) gives the A, B and X of the full hull with the full lid, and the field values at
    !> the same 18 points, which lie on both sides of y = 0, on the lid plane and in the fluid. Differences are relative to the
    !> largest entry of A, of B and of X, and of the diffraction and of the radiation values.
    subroutine lid_half_matches_full(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: full, half
        complex(8), allocatable :: dif_f(:), rad_f(:,:), dif_h(:), rad_h(:,:)
        real(8) :: worst

        call read_coefficients("lid_ann_160_pts", omega, full, error)
        if (allocated(error)) return
        call read_coefficients("lid_ann_160_half", omega, half, error)
        if (allocated(error)) return
        worst = max(maxval(abs(half%a - full%a))/maxval(abs(full%a)), maxval(abs(half%b - full%b))/maxval(abs(full%b)), &
                    maxval(abs(half%x - full%x))/maxval(abs(full%x)))
        call check_below(error, "lid_half_matches_full: A, B, X", worst, half_limit)
        if (allocated(error)) return

        call read_field("lid_ann_160_pts", dif_f, rad_f, error)
        if (allocated(error)) return
        call read_field("lid_ann_160_half", dif_h, rad_h, error)
        if (allocated(error)) return
        call check(error, size(dif_f) == 18 .and. size(dif_h) == 18 .and. size(rad_f, 1) == 18 .and. size(rad_h, 1) == 18, &
                   "lid_half_matches_full: not 18 points in each")
        if (allocated(error)) return
        worst = max(maxval(abs(dif_h - dif_f))/maxval(abs(dif_f)), maxval(abs(rad_h - rad_f))/maxval(abs(rad_f)))
        call check_below(error, "lid_half_matches_full: field values", worst, half_field_limit)
    end subroutine lid_half_matches_full

    !> IRSP = 1 with the uniform waterplane mesh of image_cos_lid, solved by least squares with the interior rows, gives the
    !> five key values of lid_ann_160 within 7e-3; measured 3.2e-3 (B11) with epsilon 0.2. The phi = 0 lid gave 1.1e-2, whose
    !> cause was left open (README).
    subroutine lid_irsp_matches(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: plain, irsp

        call read_coefficients("lid_ann_160", omega, plain, error)
        if (allocated(error)) return
        call read_coefficients("lid_ann_160_irsp", omega, irsp, error)
        if (allocated(error)) return
        call check_below(error, "lid_irsp_matches", maxval(abs(key_values(irsp) - key_values(plain))/abs(key_values(plain))), &
                         irsp_limit)
    end subroutine lid_irsp_matches

    ! -------------------------------------------------------------------------------------------- independent reference

    !> The bare bottom-mounted cylinder (24 x 10 panels, open at the sea bed) gives the MacCamy-Fuchs surge force and the
    !> elevations at (0, 8, 0) and (8, 0, 0) within 2%, as complex numbers. This checks the hull on the sea bed and the output
    !> convention of ExcForce.3 and PressureElevation.6p before the lid is compared.
    subroutine mm_nolid_matches_maccamy_fuchs(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=*), parameter :: labels(3) = [character(len=13) :: "X1", "eta(0, 8, 0)", "eta(8, 0, 0)"]
        type(coefficients) :: c
        complex(8), allocatable :: dif(:), rad(:,:)
        real(8) :: diff(3)
        integer :: i

        call read_coefficients("mm_cyl_nolid", mm_omega, c, error)
        if (allocated(error)) return
        call read_field("mm_cyl_nolid", dif, rad, error)
        if (allocated(error)) return
        call check(error, size(dif) == 2, "mm_cyl_nolid: not 2 field points")
        if (allocated(error)) return
        diff = abs([c%x(1), dif(1), dif(2)] - maccamy_fuchs)/abs(maccamy_fuchs)
        write (*, '(4x, a, 3(2x, a, " ", es10.3))') "mm_cyl_nolid", (trim(labels(i)), diff(i), i = 1, 3)
        call check_below(error, "mm_nolid_matches_maccamy_fuchs: worst", maxval(diff), maccamy_fuchs_limit)
    end subroutine mm_nolid_matches_maccamy_fuchs

    !> X1 and the two elevations of a reference case, as complex numbers.
    subroutine mm_values(name, v, error)
        character(len=*), intent(in) :: name
        complex(8), intent(out) :: v(3)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: c
        complex(8), allocatable :: dif(:), rad(:,:)

        v = (0d0, 0d0)
        call read_coefficients(name, mm_omega, c, error)
        if (allocated(error)) return
        call read_field(name, dif, rad, error)
        if (allocated(error)) return
        call check(error, size(dif) >= 2, name//": fewer than 2 field points")
        if (allocated(error)) return
        v = [c%x(1), dif(1), dif(2)]
    end subroutine mm_values

    !> The lid-effect error: the change the lid makes in HAMS-MREL (the lid case minus mm_cyl_nolid) against the change it
    !> makes in the reference (the lid reference minus MacCamy-Fuchs), relative to the reference change, worst of the three
    !> values. The hull error, common to both runs (0.85% in X1), drops out to first order.
    real(8) function lid_effect_error(lid, bare, ref)
        complex(8), intent(in) :: lid(3), bare(3), ref(3)
        lid_effect_error = maxval(abs((lid - bare) - (ref - maccamy_fuchs))/abs(ref - maccamy_fuchs))
    end function lid_effect_error

    !> The lid r 2 to 4, 48 x 8 panels, at epsilon 0.2, 1 and 1e6 m/s: X1 and the elevations at (0, 8, 0) and (8, 0, 0) within
    !> 1.2e-2 of the eigenfunction reference as complex numbers, and the lid effect within 3e-2 of the reference lid effect.
    !> The first is set by the hull error; at 1 m/s HAMS-MREL is 5.38e-2 from the conjugate reference (Im(alpha) flipped), so
    !> a flipped sign fails it. The second removes the hull error; a flipped binary gives 2.01, so it fails at 0.2 m/s too.
    subroutine mm_lid_matches_reference(error)
        type(error_type), allocatable, intent(out) :: error
        complex(8) :: v(3), bare(3)
        real(8) :: diff(3), worst, effect, worst_effect
        integer :: k

        call mm_values("mm_cyl_nolid", bare, error)
        if (allocated(error)) return
        worst = 0d0
        worst_effect = 0d0
        do k = 1, size(mm_lid_cases)
            call mm_values(trim(mm_lid_cases(k)), v, error)
            if (allocated(error)) return
            diff = abs(v - mm_reference(:, k))/abs(mm_reference(:, k))
            effect = lid_effect_error(v, bare, mm_reference(:, k))
            write (*, '(4x, a20, " eps ", es8.1, ":  X1 ", es10.3, "  eta(0, 8, 0) ", es10.3, "  eta(8, 0, 0) ", es10.3, &
                   &"  lid effect ", es10.3)') mm_lid_cases(k), mm_eps(k), diff, effect
            worst = max(worst, maxval(diff))
            worst_effect = max(worst_effect, effect)
        end do
        call check_below(error, "mm_lid_matches_reference: worst", worst, mm_reference_limit)
        if (allocated(error)) return
        call check_below(error, "mm_lid_matches_reference: worst lid effect", worst_effect, mm_effect_limit)
    end subroutine mm_lid_matches_reference

    !> The lid-effect error at 0.2 m/s falls from mm_cyl_lid_a (2 rings) to _b (4) to _c (8). The error against the reference
    !> itself does not fall: the hull error sets it (README).
    subroutine mm_lid_converges(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=12), parameter :: names(3) = [character(len=12) :: "mm_cyl_lid_a", "mm_cyl_lid_b", "mm_cyl_lid_c"]
        complex(8) :: v(3), bare(3)
        real(8) :: e(3)
        integer :: k

        call mm_values("mm_cyl_nolid", bare, error)
        if (allocated(error)) return
        do k = 1, 3
            call mm_values(names(k), v, error)
            if (allocated(error)) return
            e(k) = lid_effect_error(v, bare, mm_reference(:, 1))
        end do
        write (*, '(4x, "mm_lid_converges: lid-effect error at 2, 4, 8 rings ", 3es10.3)') e
        call check(error, e(1) > e(2) .and. e(2) > e(3), "mm_lid_converges: the lid-effect error does not fall with the rings")
    end subroutine mm_lid_converges

    !> The lid power from the elevations at the 384 lid centroids, P = 1/2 rho omega**2 sum_l epsilon abs(eta_l)**2 A_l per
    !> unit incident amplitude, is positive and within 8e-3 of the reference lid power at 0.2 and 1 m/s.
    subroutine mm_lid_power_matches(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), parameter :: rho = 1000d0
        complex(8), allocatable :: dif(:), rad(:,:)
        real(8) :: area(n384), power, rel(2)
        integer :: k

        do k = 1, 2
            call read_lid_areas(trim(mm_lid_cases(k)), n384, area, error)
            if (allocated(error)) return
            call read_field(trim(mm_lid_cases(k)), dif, rad, error)
            if (allocated(error)) return
            call check(error, size(dif) == n384 + 2, trim(mm_lid_cases(k))//": not 386 field points")
            if (allocated(error)) return
            power = 0.5d0*rho*mm_omega**2*sum(mm_eps(k)*abs(dif(3:n384 + 2))**2*area)
            call check(error, power > 0d0, "mm_lid_power_matches: the lid power is not positive")
            if (allocated(error)) return
            rel(k) = (power - mm_lid_power(k))/mm_lid_power(k)
            write (*, '(4x, a20, " lid power ", es12.5, " W, reference ", es12.5, ", relative ", es10.3)') &
                mm_lid_cases(k), power, mm_lid_power(k), rel(k)
        end do
        call check_below(error, "mm_lid_power_matches: worst", maxval(abs(rel)), mm_power_limit)
    end subroutine mm_lid_power_matches

    ! -------------------------------------------------------------------------------------------- two bodies

    !> A and B of AmssDamp.1 of a two-body case at frequency freq, 12 x 12.
    subroutine read_ab_two(name, freq, a, b, error)
        character(len=*), intent(in) :: name
        real(8), intent(in) :: freq
        real(8), intent(out) :: a(12,12), b(12,12)
        type(error_type), allocatable, intent(out) :: error
        character(len=512) :: line
        real(8) :: f, v(2)
        integer :: unit, stat, i, j, n

        a = 0d0
        b = 0d0
        n = 0
        open (newunit=unit, file=output_dir(name)//"Wamit_format/AmssDamp.1", status="old", action="read", iostat=stat)
        if (stat /= 0) then
            call test_failed(error, "No AmssDamp.1 for "//name)
            return
        end if
        do
            read (unit, '(a)', iostat=stat) line
            if (stat /= 0) exit
            if (columns(line) /= 5) cycle
            read (line, *) f, i, j, v
            if (abs(f - freq) > 1d-6 .or. i > 12 .or. j > 12) cycle
            a(i, j) = v(1)
            b(i, j) = v(2)
            n = n + 1
        end do
        close (unit)
        call check(error, n == 144, name//": AmssDamp.1 has not 144 entries at the frequency")
    end subroutine read_ab_two

    !> The rows of PressureElevation.6p of a two-body case: dif(ip) the diffraction value of point ip (5 columns, value in the
    !> last two), rad(ip, m) the radiation value of mode m = 1 to 12, mode 6(k - 1) + j of body k: a row of 14 columns, point
    !> number second and the 6 values of body 1, then a row of 12 columns with the 6 values of body 2.
    subroutine read_field_two(name, dif, rad, error)
        character(len=*), intent(in) :: name
        complex(8), allocatable, intent(out) :: dif(:), rad(:,:)
        type(error_type), allocatable, intent(out) :: error
        character(len=512) :: line
        real(8) :: row(14)
        complex(8), allocatable :: bigger(:,:)
        integer :: unit, stat, n, ip, m, last

        allocate (dif(0), rad(0, 12))
        last = 0
        open (newunit=unit, file=output_dir(name)//"Wamit_format/PressureElevation.6p", status="old", action="read", iostat=stat)
        if (stat /= 0) then
            call test_failed(error, "No PressureElevation.6p for "//name)
            return
        end if
        do
            read (unit, '(a)', iostat=stat) line
            if (stat /= 0) exit
            n = columns(line)
            if (n == 5) then
                read (line, *) row(1:5)
                ip = nint(row(3))
                if (ip > size(dif)) dif = [dif, (cmplx(0d0, 0d0, kind=8), m = size(dif) + 1, ip)]
                dif(ip) = cmplx(row(4), row(5), kind=8)
            else if (n == 14) then
                read (line, *) row(1:14)
                ip = nint(row(2))
                if (ip > size(rad, 1)) then
                    allocate (bigger(ip, 12))
                    bigger = (0d0, 0d0)
                    bigger(1:size(rad, 1), :) = rad
                    call move_alloc(bigger, rad)
                end if
                do m = 1, 6
                    rad(ip, m) = cmplx(row(2*m + 1), row(2*m + 2), kind=8)
                end do
                last = ip
            else if (n == 12 .and. last > 0) then
                read (line, *) row(1:12)
                do m = 1, 6
                    rad(last, 6 + m) = cmplx(row(2*m - 1), row(2*m), kind=8)
                end do
                last = 0
            end if
        end do
        close (unit)
    end subroutine read_field_two

    !> Eigenvalues of the symmetric matrix s, by cyclic Jacobi rotations.
    function symmetric_eigenvalues(s) result(w)
        real(8), intent(in) :: s(:,:)
        real(8) :: w(size(s, 1))
        real(8) :: m(size(s, 1), size(s, 1)), off, theta, t, c, sn, tau, apq
        integer :: n, p, q, r, sweep

        n = size(s, 1)
        m = s
        do sweep = 1, 100
            off = 0d0
            do p = 1, n - 1
                do q = p + 1, n
                    off = off + m(p, q)**2
                end do
            end do
            if (off <= 1d-30*sum(m**2)) exit
            do p = 1, n - 1
                do q = p + 1, n
                    apq = m(p, q)
                    if (apq == 0d0) cycle
                    theta = (m(q, q) - m(p, p))/(2d0*apq)
                    t = sign(1d0, theta)/(abs(theta) + sqrt(theta*theta + 1d0))
                    c = 1d0/sqrt(t*t + 1d0)
                    sn = t*c
                    tau = sn/(1d0 + c)
                    do r = 1, n
                        if (r == p .or. r == q) cycle
                        t = m(r, p)
                        m(r, p) = t - sn*(m(r, q) + tau*t)
                        m(p, r) = m(r, p)
                        m(r, q) = m(r, q) + sn*(t - tau*m(r, q))
                        m(q, r) = m(r, q)
                    end do
                    t = sn/c*apq
                    m(p, p) = m(p, p) - t
                    m(q, q) = m(q, q) + t
                    m(p, q) = 0d0
                    m(q, p) = 0d0
                end do
            end do
        end do
        do p = 1, n
            w(p) = m(p, p)
        end do
    end function symmetric_eigenvalues

    !> A lid with a second body 5 km away gives the A, B and X of body 1 of the single-body run with the same lid, in the five
    !> key values, within the limit of multi_matches_single. The far body still radiates into this one: B22 and B44, which
    !> radiate towards it, differ the most (README). The five key values differ by 4.5e-4 at most with epsilon 0.2.
    subroutine lid_multi_matches_single(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: single, multi

        call read_coefficients("lid_ann_160", omega, single, error)
        if (allocated(error)) return
        call read_coefficients("lid_two_far", omega, multi, error)
        if (allocated(error)) return
        call check_below(error, "lid_multi_matches_single", &
                         maxval(abs(key_values(multi) - key_values(single))/abs(key_values(single))), two_far_limit)
    end subroutine lid_multi_matches_single

    !> The 12 x 12 A and B of two boxes with a lid in the gap are symmetric, by reciprocity.
    subroutine lid_gap_symmetric_ab(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: a(12,12), b(12,12), worst_a, worst_b

        call read_ab_two("lid_gap", omega, a, b, error)
        if (allocated(error)) return
        worst_a = maxval(abs(a - transpose(a)))/maxval(abs(a))
        worst_b = maxval(abs(b - transpose(b)))/maxval(abs(b))
        write (*, '(4x, "lid_gap: A ", es10.3, ", B ", es10.3)') worst_a, worst_b
        call check_below(error, "lid_gap_symmetric_ab", max(worst_a, worst_b), gap_symmetry_limit)
    end subroutine lid_gap_symmetric_ab

    !> The lid absorbs, so the 12 x 12 B of the two boxes, symmetric part, has no negative eigenvalue beyond the hull
    !> discretisation: the two boxes without the lid have -2.3e-5 at this mesh, and both fall with the hull mesh (README).
    subroutine lid_gap_damping_nonnegative(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: a(12,12), b(12,12), w(12)

        call read_ab_two("lid_gap", omega, a, b, error)
        if (allocated(error)) return
        w = symmetric_eigenvalues(0.5d0*(b + transpose(b)))
        write (*, '(4x, "lid_gap: eigenvalues of B over max(abs(B)): ", 12es10.2)') w/maxval(abs(b))
        call check_below(error, "lid_gap_damping_nonnegative: -min eigenvalue/max(abs(B))", -minval(w)/maxval(abs(b)), &
                         gap_damping_floor)
    end subroutine lid_gap_damping_nonnegative

    !> X of ExcForce.3 of a two-body case at frequency freq and heading beta, the 12 modes.
    subroutine read_x_two(name, freq, beta, x, error)
        character(len=*), intent(in) :: name
        real(8), intent(in) :: freq, beta
        complex(8), intent(out) :: x(12)
        type(error_type), allocatable, intent(out) :: error
        character(len=512) :: line
        real(8) :: f, b, row(4)
        integer :: unit, stat, i, n

        x = (0d0, 0d0)
        n = 0
        open (newunit=unit, file=output_dir(name)//"Wamit_format/ExcForce.3", status="old", action="read", iostat=stat)
        if (stat /= 0) then
            call test_failed(error, "No ExcForce.3 for "//name)
            return
        end if
        do
            read (unit, '(a)', iostat=stat) line
            if (stat /= 0) exit
            if (columns(line) /= 7) cycle
            read (line, *) f, b, i, row
            if (abs(f - freq) < 1d-6 .and. abs(b - beta) < 1d-6 .and. i <= 12) then
                x(i) = cmplx(row(3), row(4), kind=8)
                n = n + 1
            end if
        end do
        close (unit)
        call check(error, n == 12, name//": ExcForce.3 has not 12 entries at the frequency and heading")
    end subroutine read_x_two

    !> epsilon = 0 on the gap lid is no lid: the 12 x 12 A and B and the 12 forces at heading 90 of lid_gap_eps0 are those of
    !> lid_gap_nolid, relative to the largest entry of each, and so is the diffraction elevation at the 83 field points.
    subroutine lid_gap_eps0_matches_no_lid(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: a0(12,12), b0(12,12), ab(12,12), bb(12,12), worst, worst_f
        complex(8) :: x0(12), xb(12)
        complex(8), allocatable :: dif0(:), rad0(:,:), difb(:), radb(:,:)

        call read_ab_two("lid_gap_eps0", omega, a0, b0, error)
        if (allocated(error)) return
        call read_ab_two("lid_gap_nolid", omega, ab, bb, error)
        if (allocated(error)) return
        call read_x_two("lid_gap_eps0", omega, 90d0, x0, error)
        if (allocated(error)) return
        call read_x_two("lid_gap_nolid", omega, 90d0, xb, error)
        if (allocated(error)) return
        worst = max(maxval(abs(a0 - ab))/maxval(abs(ab)), maxval(abs(b0 - bb))/maxval(abs(bb)), maxval(abs(x0 - xb))/maxval(abs(xb)))
        call read_field_two("lid_gap_eps0", dif0, rad0, error)
        if (allocated(error)) return
        call read_field_two("lid_gap_nolid", difb, radb, error)
        if (allocated(error)) return
        call check(error, size(dif0) == 83 .and. size(difb) == 83, "lid_gap_eps0_matches_no_lid: not 83 points")
        if (allocated(error)) return
        worst_f = maxval(abs(dif0 - difb))/maxval(abs(difb))
        write (*, '(4x, "lid_gap_eps0_matches_no_lid: A, B, X ", es10.3, ", diffraction ", es10.3)') worst, worst_f
        call check_below(error, "lid_gap_eps0_matches_no_lid", max(worst, worst_f), exact_limit)
    end subroutine lid_gap_eps0_matches_no_lid

    !> The lid damps the gap: the largest diffraction elevation at the 80 gap centroids, over the incident amplitude, falls from
    !> no lid (lid_gap_nolid) to epsilon 0.2 (lid_gap) to epsilon 1 (lid_gap_eps1), while the wave off the lid stays.
    subroutine lid_gap_damps(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=13), parameter :: names(3) = [character(len=13) :: "lid_gap_nolid", "lid_gap", "lid_gap_eps1"]
        complex(8), allocatable :: dif(:), rad(:,:)
        real(8) :: peak(3)
        integer :: k

        do k = 1, 3
            call read_field_two(trim(names(k)), dif, rad, error)
            if (allocated(error)) return
            call check(error, size(dif) == 83, trim(names(k))//": not 83 points")
            if (allocated(error)) return
            peak(k) = maxval(abs(dif(1:80)))
            write (*, '(4x, a13, ": largest |elevation| on the gap lid ", es12.5, ", off the lid ", 3es10.3)') names(k), &
                peak(k), abs(dif(81:83))
            call check(error, maxval(abs(dif(81:83))) > 0.1d0, trim(names(k))//": no wave off the lid")
            if (allocated(error)) return
        end do
        call check(error, peak(2) < peak(1) .and. peak(3) < peak(2), "lid_gap_damps: the gap elevation does not fall with epsilon")
    end subroutine lid_gap_damps

    !> Panel centroids, areas and unit normals out of the body, from the hull of the cylinder of the image term cases (copied
    !> from test_integration_image, whose own copy is private).
    subroutine read_panels(xc, area, normal, error)
        real(8), allocatable, intent(out) :: xc(:,:), area(:), normal(:,:)
        type(error_type), allocatable, intent(out) :: error
        character(len=512) :: line
        real(8), allocatable :: node(:,:)
        real(8) :: p(3,4), v(3), scale_v
        integer :: unit, stat, np, nn, isx, isy, i, k, id, nv, ids(4)

        open (newunit=unit, file=input_dir("image_term", "image_cos")//"HullMesh.pnl", status="old", action="read", iostat=stat)
        if (stat /= 0) then
            call test_failed(error, "No hull mesh of image_cos")
            return
        end if
        read (unit, '(a)') line
        read (unit, '(a)') line
        read (unit, '(a)') line
        read (unit, *) np, nn, isx, isy
        allocate (node(3, nn), xc(3, np), area(np), normal(3, np))
        read (unit, '(a)') line
        read (unit, '(a)') line
        do i = 1, nn
            read (unit, *) id, node(:, i)
        end do
        read (unit, '(a)') line
        read (unit, '(a)') line
        read (unit, '(a)') line
        do i = 1, np
            read (unit, *) id, nv, ids(1:nv)
            do k = 1, nv
                p(:, k) = node(:, ids(k))
            end do
            xc(:, i) = sum(p(:, 1:nv), dim=2)/real(nv, 8)
            if (nv == 3) then
                v = 0.5d0*cross(p(:, 2) - p(:, 1), p(:, 3) - p(:, 1))
            else
                v = 0.5d0*cross(p(:, 3) - p(:, 1), p(:, 4) - p(:, 2))
            end if
            area(i) = norm2(v)
            scale_v = 1d0/area(i)
            ! A cylinder on the z axis: out of the body is away from the axis, and down on the bottom
            if (abs(v(3)) > 0.5d0*area(i)) then
                normal(:, i) = merge(-1d0, 1d0, v(3) > 0d0)*scale_v*v
            else
                normal(:, i) = merge(1d0, -1d0, v(1)*xc(1, i) + v(2)*xc(2, i) > 0d0)*scale_v*v
            end if
        end do
        close (unit)
    end subroutine read_panels

    pure function cross(a, b) result(c)
        real(8), intent(in) :: a(3), b(3)
        real(8) :: c(3)
        c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
    end function cross

end module test_integration_lid
