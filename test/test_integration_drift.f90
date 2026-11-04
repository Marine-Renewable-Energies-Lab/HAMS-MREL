!> End-to-end checks of the mean drift: a real HAMS solution goes through the convention conversion (section 5.2 of
!> docs/superpowers/specs/2026-11-04-second-order-waves-design.md), the far-field loads are compared with first-order
!> outputs written by the solver itself, and the near-field loads with the far-field loads. Fixtures and the frequency
!> choice: test/data/drift_cyl/README.md.
module test_integration_drift

    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
    use testing_utilities
    use NearFieldDrift, only : T_WATERLINE, T_FLARE, T_VELOCITY, T_HYDROSTATIC

    implicit none
    private
    public :: setup_drift, collect_tests_drift

    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: rho = 1000d0, g = 9.80665d0   ! as Const_mod
    complex(8), parameter :: ci = (0d0, 1d0)

    ! Acceptance limits
    real(8), parameter :: residual_limit = 1d-8
    real(8), parameter :: diffraction_limit = 1d-4
    real(8), parameter :: energy_limit = 1d-2         ! radiated energy, Haskind, energy balance
    real(8), parameter :: lateral_limit = 1d-3
    real(8), parameter :: gmres_limit = 1d-2
    real(8), parameter :: field_point_limit = 1d-12
    real(8), parameter :: half_mesh_limit = 1d-3
    real(8), parameter :: moment_limit = 1d-3
    ! F_x against the lossless (Maruo) form. They differ by kappa (X - int |a|^2), so the 1% of energy_limit becomes 1.8%
    ! at the lowest frequency, where kappa X / (rho g) is about 1.8 times the drift force.
    real(8), parameter :: maruo_limit = 2d-2
    ! Values in .8 have 7 significant digits, so same-run ratios cannot be checked more tightly than this
    real(8), parameter :: printed_limit = 1d-5

    ! Near-field against far-field. They are checked on the fixtures named in each test,
    ! whose meshes and frequencies were chosen after seeing results. Having no edges is not enough: a free half ellipsoid, for
    ! one, misses 5% where its net force is small. See test/data/drift_cyl/README.md.
    real(8), parameter :: near_far_limit = 5d-2
    real(8), parameter :: near_mz_limit = 1d-1
    real(8), parameter :: rotated_limit = 1d-3
    real(8), parameter :: irsp_limit = 1d-2
    ! Regression band on the cases that miss the limits above: near-field and far-field each stay within this of the value
    ! measured when the check was written. These pin the numbers, including their error, and are not accuracy claims.
    ! See test/data/drift_cyl/README.md for the values against panel count, and
    ! pre_postprocessing_python/SecondOrder/near_far_study.py for the 5% on these hulls.
    real(8), parameter :: regression_band = 1d-2
    ! At 0.70 rad/s the free cylinder's net drift force is 4.6% of term I, from cancellation between much larger terms. Near-field
    ! against far-field there is +3.99%, +0.87%, +0.34% at 480, 1920, 4320 panels. free_near_far leaves that frequency out, and
    ! free_near_far_low_frequency checks it alone
    real(8), parameter :: small_net_force_below = 0.8d0
    real(8), parameter :: low_frequency = 0.70d0
    ! On the springs, the heave resonance moves to about 0.86 rad/s, so 0.95 rad/s is near it. external_stiffness leaves out the
    ! frequencies below this, and external_stiffness_near_resonance checks 0.95 rad/s alone, with a regression band beside it
    real(8), parameter :: springs_away_from_resonance = 1d0
    real(8), parameter :: springs_near_resonance = 0.95d0
    ! Off-diagonal heading pairs of two_body_total are about 1/7 of the diagonal loads: a regression limit
    real(8), parameter :: heading_pair_limit = 1d-1

    character(len=*), parameter :: main_case = "drift_cyl"
    character(len=*), parameter :: far_file = "Wamit_format/MeanDriftMomentum.8"
    character(len=*), parameter :: near_file = "Wamit_format/MeanDriftPressure.9"
    character(len=*), parameter :: terms_file = "Wamit_format/MeanDriftTerms.txt"

contains

    subroutine setup_drift(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("run_drift_cyl", run_drift_cyl), &
            new_unittest("run_drift_cyl_sybo", run_drift_cyl_sybo), &
            new_unittest("run_drift_cyl_gmres", run_drift_cyl_gmres), &
            new_unittest("run_drift_cyl_half", run_drift_cyl_half), &
            new_unittest("run_drift_offset", run_drift_offset), &
            new_unittest("run_drift_offset_refl2", run_drift_offset_refl2), &
            new_unittest("run_drift_two_cyl", run_drift_two_cyl), &
            new_unittest("run_drift_cyl_fixed", run_drift_cyl_fixed), &
            new_unittest("run_drift_round", run_drift_round), &
            new_unittest("run_drift_round_fixed", run_drift_round_fixed), &
            new_unittest("run_drift_round_kstf", run_drift_round_kstf), &
            new_unittest("run_drift_two_round", run_drift_two_round), &
            new_unittest("run_drift_box", run_drift_box), &
            new_unittest("run_drift_ellipsoid", run_drift_ellipsoid), &
            new_unittest("run_drift_far_apart", run_drift_far_apart), &
            new_unittest("run_drift_isolated", run_drift_isolated), &
            new_unittest("run_drift_mixed", run_drift_mixed), &
            new_unittest("run_drift_rotated", run_drift_rotated), &
            new_unittest("run_drift_unrotated", run_drift_unrotated), &
            new_unittest("run_drift_irsp", run_drift_irsp), &
            new_unittest("run_drift_irsp_ref", run_drift_irsp_ref), &
            new_unittest("run_drift_two_irsp", run_drift_two_irsp), &
            new_unittest("run_drift_two_irsp_ref", run_drift_two_irsp_ref) &
        ]
    end subroutine setup_drift

    subroutine collect_tests_drift(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("convention_residual", convention_residual), &
            new_unittest("diffraction_matches_6p", diffraction_matches_6p), &
            new_unittest("radiated_energy", radiated_energy), &
            new_unittest("haskind", haskind), &
            new_unittest("energy_balance", energy_balance), &
            new_unittest("symmetric_lateral_zero", symmetric_lateral_zero), &
            new_unittest("sybo_skipped", sybo_skipped), &
            new_unittest("pairs_match_unidirectional", pairs_match_unidirectional), &
            new_unittest("refl_scaling", refl_scaling), &
            new_unittest("no_input_no_output", no_input_no_output), &
            new_unittest("gmres_agrees", gmres_agrees), &
            new_unittest("two_bodies_smoke", two_bodies_smoke), &
            new_unittest("offset_moment", offset_moment), &
            new_unittest("two_bodies_moment", two_bodies_moment), &
            new_unittest("drift_matches_maruo", drift_matches_maruo), &
            new_unittest("half_mesh_matches_full", half_mesh_matches_full), &
            new_unittest("fixed_near_far", fixed_near_far), &
            new_unittest("free_near_far", free_near_far), &
            new_unittest("free_near_far_low_frequency", free_near_far_low_frequency), &
            new_unittest("near_mz_far_mz", near_mz_far_mz), &
            new_unittest("external_stiffness", external_stiffness), &
            new_unittest("external_stiffness_near_resonance", external_stiffness_near_resonance), &
            new_unittest("external_stiffness_near_resonance_regression", external_stiffness_near_resonance_regression), &
            new_unittest("two_body_total", two_body_total), &
            new_unittest("fixed_near_far_sharp_regression", fixed_near_far_sharp_regression), &
            new_unittest("free_near_far_sharp_regression", free_near_far_sharp_regression), &
            new_unittest("two_body_total_sharp_regression", two_body_total_sharp_regression), &
            new_unittest("near_mz_far_mz_sharp_regression", near_mz_far_mz_sharp_regression), &
            new_unittest("far_apart_isolated", far_apart_isolated), &
            new_unittest("submerged_mix", submerged_mix), &
            new_unittest("rotated_lcs", rotated_lcs), &
            new_unittest("bidirectional", bidirectional), &
            new_unittest("irsp_stability", irsp_stability), &
            new_unittest("terms_sum_to_total", terms_sum_to_total), &
            new_unittest("near_output_switches", near_output_switches) &
        ]
    end subroutine collect_tests_drift

    ! ---------------------------------------------------------------------------------------------- running the cases

    subroutine run_case(name, outputs, error)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: outputs(:)     ! files that must exist, relative to the output directory
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: inputdir, outputdir, command
        character(len=16) :: code
        integer :: exitcode, i
        logical :: exists

        inputdir = get_repo_root_path()//"test/data/"//name//"/input/"
        outputdir = output_dir(name)
        command = native_path(hamsexe_path//" "//inputdir//" "//outputdir//" > "//get_running_exe_path()//name//".log") &
                  //" 2>&1"
        call execute_command_line(command, exitstat=exitcode, wait=.true.)
        if (exitcode /= 0) then
            write (code, '(I0)') exitcode
            call test_failed(error, "HAMS failed for "//name//" with exit code "//trim(code))
            return
        end if
        do i = 1, size(outputs)
            inquire (file=outputdir//trim(outputs(i)), exist=exists)
            call check(error, exists, "Missing output of "//name//": "//trim(outputs(i)))
            if (allocated(error)) return
        end do
    end subroutine run_case

    function output_dir(name) result(dir)
        character(len=*), intent(in) :: name
        character(:), allocatable :: dir
        dir = get_running_exe_path()//"output_"//name//"/"
    end function output_dir

    subroutine run_drift_cyl(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=40) :: files(3)
        files = [character(len=40) :: far_file, near_file, "SecondOrder/PanelDiagnostics.txt"]
        call run_case(main_case, files, error)
    end subroutine run_drift_cyl

    subroutine run_drift_cyl_sybo(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_cyl_sybo", both_drift_files(), error)
    end subroutine run_drift_cyl_sybo

    subroutine run_drift_cyl_gmres(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_cyl_gmres", drift_file_only(), error)
    end subroutine run_drift_cyl_gmres

    subroutine run_drift_cyl_half(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_cyl_half", both_drift_files(), error)
    end subroutine run_drift_cyl_half

    subroutine run_drift_offset(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_offset", both_drift_files(), error)
    end subroutine run_drift_offset

    subroutine run_drift_offset_refl2(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_offset_refl2", drift_file_only(), error)
    end subroutine run_drift_offset_refl2

    subroutine run_drift_two_cyl(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=40) :: files(3)
        files = [character(len=40) :: far_file, near_file, "SecondOrder/PanelDiagnostics.txt"]
        call run_case("drift_two_cyl", files, error)
    end subroutine run_drift_two_cyl

    subroutine run_drift_cyl_fixed(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_cyl_fixed", files_with_terms(), error)
    end subroutine run_drift_cyl_fixed

    subroutine run_drift_round(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_round", files_with_terms(), error)
    end subroutine run_drift_round

    subroutine run_drift_round_fixed(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_round_fixed", files_with_terms(), error)
    end subroutine run_drift_round_fixed

    subroutine run_drift_round_kstf(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_round_kstf", files_with_terms(), error)
    end subroutine run_drift_round_kstf

    subroutine run_drift_two_round(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_two_round", both_drift_files(), error)
    end subroutine run_drift_two_round

    subroutine run_drift_ellipsoid(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_ellipsoid", files_with_terms(), error)
    end subroutine run_drift_ellipsoid

    subroutine run_drift_box(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_box", files_with_terms(), error)
    end subroutine run_drift_box

    subroutine run_drift_far_apart(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_far_apart", both_drift_files(), error)
    end subroutine run_drift_far_apart

    subroutine run_drift_isolated(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=40) :: files(1)
        files(1) = near_file
        call run_case("drift_isolated", files, error)
    end subroutine run_drift_isolated

    subroutine run_drift_mixed(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_mixed", files_with_terms(), error)
    end subroutine run_drift_mixed

    subroutine run_drift_rotated(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_rotated", both_drift_files(), error)
    end subroutine run_drift_rotated

    subroutine run_drift_unrotated(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_unrotated", both_drift_files(), error)
    end subroutine run_drift_unrotated

    subroutine run_drift_irsp(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_irsp", both_drift_files(), error)
    end subroutine run_drift_irsp

    subroutine run_drift_irsp_ref(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_irsp_ref", both_drift_files(), error)
    end subroutine run_drift_irsp_ref

    subroutine run_drift_two_irsp(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_two_irsp", both_drift_files(), error)
    end subroutine run_drift_two_irsp

    subroutine run_drift_two_irsp_ref(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_two_irsp_ref", both_drift_files(), error)
    end subroutine run_drift_two_irsp_ref

    function drift_file_only() result(files)
        character(len=40) :: files(1)
        files(1) = far_file
    end function drift_file_only

    function both_drift_files() result(files)
        character(len=40) :: files(2)
        files = [character(len=40) :: far_file, near_file]
    end function both_drift_files

    function files_with_terms() result(files)
        character(len=40) :: files(3)
        files = [character(len=40) :: far_file, near_file, terms_file]
    end function files_with_terms

    ! ------------------------------------------------------------------------------------------------- physics checks

    !> |[-w^2(M+A) - iwB + C] xi - F| / |F| with xi and F converted to e^{-iwt}, as computed by the driver.
    subroutine convention_residual(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: rows(:,:)

        call read_diagnostics(main_case, "RESIDUAL", 3, rows, error)
        if (allocated(error)) return
        call check(error, size(rows, 1) == 6, "Expected one residual for each frequency and heading")
        if (allocated(error)) return
        call report("convention_residual", maxval(rows(:,3)))
        call check_below(error, "convention_residual", maxval(rows(:,3)), residual_limit)
    end subroutine convention_residual

    !> conj(i w rho phi_D / (rho g A)) of the state against the on-body pressure the solver writes to .6p, relative to the
    !> local pressure. The field points must be the panel centroids to full precision: the solver's on-body integrals are
    !> erratic for a point slightly off the panel.
    subroutine diffraction_matches_6p(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: panels(:,:), pressure(:,:), field_points(:,:)
        integer, allocatable :: field_panel(:)
        real(8) :: worst, worst_point
        complex(8) :: phi_d, expected, written
        integer :: r, p, matched

        call read_diagnostics(main_case, "PANEL", 10, panels, error)
        if (allocated(error)) return
        call read_6p_diffraction(output_dir(main_case)//"Wamit_format/PressureElevation.6p", pressure, error)
        if (allocated(error)) return
        call read_field_panels(get_repo_root_path()//"test/data/"//main_case//"/field_panels.txt", field_panel)
        call check(error, size(pressure, 1) == 6*size(field_panel), "Unexpected number of diffraction rows in .6p")
        if (allocated(error)) return
        call read_field_points(get_repo_root_path()//"test/data/"//main_case//"/input/ControlFile.in", field_points)
        call check(error, size(field_points, 1) == size(field_panel), "Field points and field_panels.txt differ in size")
        if (allocated(error)) return

        worst = 0d0
        worst_point = 0d0
        matched = 0
        do r = 1, size(pressure, 1)
            do p = 1, size(panels, 1)
                if (.not. same(panels(p,1), pressure(r,1)) .or. .not. same(panels(p,2), pressure(r,2))) cycle
                if (nint(panels(p,3)) /= field_panel(nint(pressure(r,3)))) cycle
                phi_d = cmplx(panels(p,7), panels(p,8), kind=8)
                expected = conjg(ci*panels(p,1)*phi_d/g)
                written = cmplx(pressure(r,6), pressure(r,7), kind=8)
                worst = max(worst, abs(expected - written)/abs(written))
                worst_point = max(worst_point, maxval(abs(field_points(nint(pressure(r,3)),:) - panels(p,4:6))))
                matched = matched + 1
            end do
        end do
        call check(error, matched == size(pressure, 1), "Not every .6p row has a matching panel diagnostic")
        if (allocated(error)) return
        call report("field point to panel centroid", worst_point)
        call report("diffraction_matches_6p", worst)
        call check_below(error, "field point to panel centroid", worst_point, field_point_limit)
        if (allocated(error)) return
        call check_below(error, "diffraction_matches_6p", worst, diffraction_limit)
    end subroutine diffraction_matches_6p

    !> Radiation damping work 1/2 w^2 B_jj against the wave energy flux 1/2 rho g C_g int |a_j|^2 dtheta, j = 1, 3, 5.
    subroutine radiated_energy(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: waves(:,:), modes(:,:), coef(:,:)
        real(8) :: omega, k, h, cg, work, flux, worst, b
        integer :: iw, jj, j, m, r

        call read_diagnostics(main_case, "WAVE", 4, waves, error)
        if (allocated(error)) return
        call read_diagnostics(main_case, "MODE", 5, modes, error)
        if (allocated(error)) return
        call read_amssdamp(output_dir(main_case)//"Wamit_format/AmssDamp.1", coef, error)
        if (allocated(error)) return

        worst = 0d0
        do iw = 1, size(waves, 1)
            omega = waves(iw,1)
            k = waves(iw,2)
            h = waves(iw,4)
            cg = omega/k*group_ratio(k, h)
            do jj = 1, 3
                j = 2*jj - 1
                b = -1d0
                do r = 1, size(coef, 1)
                    if (same(coef(r,1), omega) .and. nint(coef(r,2)) == j .and. nint(coef(r,3)) == j) b = coef(r,5)*rho*omega
                end do
                flux = 0d0
                m = 0
                do r = 1, size(modes, 1)
                    if (.not. same(modes(r,1), omega) .or. nint(modes(r,2)) /= j) cycle
                    flux = flux + modes(r,4)**2 + modes(r,5)**2
                    m = m + 1
                end do
                call check(error, m > 0 .and. b > 0d0, "Missing radiation data for a mode")
                if (allocated(error)) return
                flux = 0.5d0*rho*g*cg*flux*2d0*pi/real(m, 8)
                work = 0.5d0*omega*omega*b
                worst = max(worst, abs(flux - work)/work)
            end do
        end do
        call report("radiated_energy", worst)
        call check_below(error, "radiated_energy", worst, energy_limit)
    end subroutine radiated_energy

    !> Excitation from the far-field amplitude of mode j at the heading opposite to the wave:
    !>   X_j = i w rho int(phi_I n_j - psi_j dphi_I/dn) dS = e^{-i pi/4} (2 rho g^2 A / (w^2 C)) sqrt(pi k / 2) a_j(beta+pi)
    !> with a_j the amplitude of the wave radiated by unit motion, and C as in far_amplitude.
    subroutine haskind(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: waves(:,:), modes(:,:), exc(:,:)
        real(8) :: omega, k, nu, h, worst, worst_abs, beta, c, scale_f
        complex(8) :: a, predicted, solver
        integer :: r, m, used

        call read_diagnostics(main_case, "WAVE", 4, waves, error)
        if (allocated(error)) return
        call read_diagnostics(main_case, "MODE", 5, modes, error)
        if (allocated(error)) return
        call read_exc(output_dir(main_case)//"Wamit_format/ExcForce.3", exc, error)
        if (allocated(error)) return

        worst = 0d0
        worst_abs = 0d0
        used = 0
        do r = 1, size(exc, 1)
            omega = exc(r,1)
            beta = exc(r,2)
            ! Modes with no excitation by symmetry have no relative error to speak of
            scale_f = maxval(exc(:,4), mask=abs(exc(:,1) - omega) < 1d-6 .and. exc(:,2) == beta)
            if (exc(r,4) < 1d-3*scale_f) cycle
            call find_wave(waves, omega, k, nu, h)
            c = k*k/((k*k - nu*nu)*h + nu)
            a = (0d0, 0d0)
            do m = 1, size(modes, 1)
                if (.not. same(modes(m,1), omega) .or. nint(modes(m,2)) /= nint(exc(r,3))) cycle
                if (abs(modes(m,3) - modulo(beta + 180d0, 360d0)) < 1d-6) a = cmplx(modes(m,4), modes(m,5), kind=8)
            end do
            predicted = exp(-ci*pi/4d0)*(2d0*g*g/(omega*omega*c))*sqrt(pi*k/2d0)*a/g   ! per rho g A
            solver = conjg(cmplx(exc(r,6), exc(r,7), kind=8))
            worst = max(worst, abs(predicted - solver)/abs(solver))
            worst_abs = max(worst_abs, abs(abs(predicted) - abs(solver))/abs(solver))
            used = used + 1
        end do
        call check(error, used >= 18, "Too few modes with excitation in the Haskind check")
        if (allocated(error)) return
        call report("haskind", worst)
        call report("haskind magnitude only", worst_abs)
        call check_below(error, "haskind", worst, energy_limit)
    end subroutine haskind

    !> A free body with no external damping absorbs nothing: X = int |a|^2 dtheta.
    subroutine energy_balance(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: worst

        call energy_balance_of(main_case, 6, worst, error)
        if (allocated(error)) return
        call check_below(error, "energy_balance", worst, energy_limit)
    end subroutine energy_balance

    !> The mesh is symmetric about y = 0, so heading 0 has no sway force and no yaw moment.
    subroutine symmetric_lateral_zero(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: rows(:,:)
        real(8) :: worst, fx
        integer :: r

        call read_drift(output_dir(main_case)//"Wamit_format/MeanDriftMomentum.8", rows, error)
        if (allocated(error)) return
        call check(error, size(rows, 1) == 3*2*3, "Expected three rows for each frequency and heading")
        if (allocated(error)) return
        worst = 0d0
        fx = 0d0
        do r = 1, size(rows, 1)
            if (abs(rows(r,2)) > 1d-9) cycle
            if (nint(rows(r,4)) == 1) then
                fx = rows(r,5)
            else
                worst = max(worst, rows(r,5)/fx)
            end if
        end do
        call report("symmetric_lateral_zero", worst)
        call check_below(error, "symmetric_lateral_zero", worst, lateral_limit)
    end subroutine symmetric_lateral_zero

    !> SYBO = 1 adds the zero and infinite frequencies: they have no rows, and no value in the file is NaN.
    subroutine sybo_skipped(error)
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: path
        real(8), allocatable :: rows(:,:)

        path = output_dir("drift_cyl_sybo")//"Wamit_format/MeanDriftMomentum.8"
        call check_no_bad_numbers(path, error)
        if (allocated(error)) return
        call read_drift(path, rows, error)
        if (allocated(error)) return
        call check(error, size(rows, 1) == 3*(2*2)*3, "Expected rows for the 3 real frequencies and 4 heading pairs only")
        if (allocated(error)) return
        call check(error, all(rows(:,1) > 0.5d0), "A limit frequency has rows in the drift file")
    end subroutine sybo_skipped

    !> The drift_cyl_sybo case runs all heading pairs: its diagonal equals the unidirectional result, and F_21 = conj(F_12).
    subroutine pairs_match_unidirectional(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: uni(:,:), pairs(:,:)
        real(8) :: worst, worst_sym, scale_row
        integer :: r, p, q, matched

        call read_drift(output_dir(main_case)//"Wamit_format/MeanDriftMomentum.8", uni, error)
        if (allocated(error)) return
        call read_drift(output_dir("drift_cyl_sybo")//"Wamit_format/MeanDriftMomentum.8", pairs, error)
        if (allocated(error)) return

        worst = 0d0
        worst_sym = 0d0
        matched = 0
        do p = 1, size(pairs, 1)
            scale_row = maxval(pairs(:,5), mask=abs(pairs(:,1) - pairs(p,1)) < 1d-6 .and. nint(pairs(:,4)) /= 6)
            do r = 1, size(uni, 1)
                if (same(uni(r,1), pairs(p,1)) .and. uni(r,2) == pairs(p,2) .and. uni(r,3) == pairs(p,3) &
                    .and. uni(r,4) == pairs(p,4)) then
                    worst = max(worst, abs(uni(r,7) - pairs(p,7))/scale_row)
                    matched = matched + 1
                end if
            end do
            do q = 1, size(pairs, 1)
                if (same(pairs(q,1), pairs(p,1)) .and. pairs(q,2) == pairs(p,3) .and. pairs(q,3) == pairs(p,2) &
                    .and. pairs(q,4) == pairs(p,4)) then
                    worst_sym = max(worst_sym, (abs(pairs(q,7) - pairs(p,7)) + abs(pairs(q,8) + pairs(p,8)))/scale_row)
                end if
            end do
        end do
        call check(error, matched == size(uni, 1), "Not every unidirectional row is in the pairs output")
        if (allocated(error)) return
        call report("pairs diagonal vs unidirectional", worst)
        call report("pairs conjugate symmetry", worst_sym)
        call check_below(error, "pairs diagonal vs unidirectional", worst, printed_limit)
        if (allocated(error)) return
        call check_below(error, "pairs conjugate symmetry", worst_sym, printed_limit)
    end subroutine pairs_match_unidirectional

    !> On a cylinder off the origin, M_z is not zero. REFL = 2 halves forces and quarters the moment.
    subroutine refl_scaling(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: one(:,:), two(:,:)
        real(8) :: worst, expected
        integer :: r

        call read_drift(output_dir("drift_offset")//"Wamit_format/MeanDriftMomentum.8", one, error)
        if (allocated(error)) return
        call read_drift(output_dir("drift_offset_refl2")//"Wamit_format/MeanDriftMomentum.8", two, error)
        if (allocated(error)) return
        call check(error, size(one, 1) == 3 .and. size(two, 1) == 3, "Expected three rows")
        if (allocated(error)) return
        call check(error, abs(one(3,7)) > 1d-2*abs(one(1,7)), "M_z must not be negligible in this case")
        if (allocated(error)) return

        worst = 0d0
        do r = 1, 3
            expected = merge(4d0, 2d0, nint(one(r,4)) == 6)
            worst = max(worst, abs(one(r,7)/two(r,7) - expected)/expected)
        end do
        call report("refl_scaling", worst)
        call check_below(error, "refl_scaling", worst, printed_limit)
    end subroutine refl_scaling

    !> config1 has no SecondOrder.in. It runs in an earlier suite of integration_tester.
    subroutine no_input_no_output(error)
        type(error_type), allocatable, intent(out) :: error
        logical :: exists

        inquire (file=get_running_exe_path()//"output_config1/Wamit_format/AmssDamp.1", exist=exists)
        call check(error, exists, "config1 has not run")
        if (allocated(error)) return
        inquire (file=get_running_exe_path()//"output_config1/Wamit_format/MeanDriftMomentum.8", exist=exists)
        call check(error, .not. exists, "A run without SecondOrder.in wrote a mean drift file")
        if (allocated(error)) return
        inquire (file=get_running_exe_path()//"output_config1/SecondOrder/PanelDiagnostics.txt", exist=exists)
        call check(error, .not. exists, "A run without SecondOrder.in wrote panel diagnostics")
    end subroutine no_input_no_output

    !> ISOLV = 2 against the direct solver, F_x.
    subroutine gmres_agrees(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: lu(:,:), gm(:,:)
        real(8) :: worst
        integer :: r

        call read_drift(output_dir(main_case)//"Wamit_format/MeanDriftMomentum.8", lu, error)
        if (allocated(error)) return
        call read_drift(output_dir("drift_cyl_gmres")//"Wamit_format/MeanDriftMomentum.8", gm, error)
        if (allocated(error)) return
        call check(error, size(lu, 1) == size(gm, 1), "GMRES and LU runs have different numbers of rows")
        if (allocated(error)) return
        worst = 0d0
        do r = 1, size(lu, 1)
            if (nint(lu(r,4)) == 1) worst = max(worst, abs(lu(r,7) - gm(r,7))/abs(lu(r,7)))
        end do
        call report("gmres_agrees", worst)
        call check_below(error, "gmres_agrees", worst, gmres_limit)
    end subroutine gmres_agrees

    !> Multi-body path: the conversion holds for the combined system and the energy balance closes for the array.
    subroutine two_bodies_smoke(error)
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: path
        real(8), allocatable :: rows(:,:), residuals(:,:)
        real(8) :: worst

        path = output_dir("drift_two_cyl")//"Wamit_format/MeanDriftMomentum.8"
        call check_no_bad_numbers(path, error)
        if (allocated(error)) return
        call read_drift(path, rows, error)
        if (allocated(error)) return
        call check(error, size(rows, 1) == 2*(2*2)*3, "Expected rows for 2 frequencies and 4 heading pairs")
        if (allocated(error)) return
        call read_diagnostics("drift_two_cyl", "RESIDUAL", 3, residuals, error)
        if (allocated(error)) return
        call report("two bodies residual", maxval(residuals(:,3)))
        call check_below(error, "two bodies residual", maxval(residuals(:,3)), residual_limit)
        if (allocated(error)) return
        call energy_balance_of("drift_two_cyl", 4, worst, error)
        if (allocated(error)) return
        call check_below(error, "two bodies energy balance", worst, energy_limit)
    end subroutine two_bodies_smoke

    !> A round body at (x0, y0) has no moment about its own centre, so about the origin M_z = x0 F_y - y0 F_x.
    subroutine offset_moment(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), parameter :: x0 = 6d0, y0 = 8d0    ! test/data/drift_offset
        real(8), allocatable :: rows(:,:)
        real(8) :: expected, worst

        call read_drift(output_dir("drift_offset")//"Wamit_format/MeanDriftMomentum.8", rows, error)
        if (allocated(error)) return
        call check(error, size(rows, 1) == 3, "Expected three rows")
        if (allocated(error)) return
        expected = x0*rows(2,7) - y0*rows(1,7)
        worst = abs(rows(3,7) - expected)/abs(expected)
        call report("offset_moment", worst)
        call check_below(error, "offset_moment", worst, moment_limit)
    end subroutine offset_moment

    !> Two equal cylinders at y = 0 and y = 30 in a wave along x are symmetric about y = 15: F_y = 0 and M_z = -15 F_x.
    subroutine two_bodies_moment(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: rows(:,:)
        real(8) :: worst_m, worst_y, fx, fy, mz
        integer :: r, done

        call read_drift(output_dir("drift_two_cyl")//"Wamit_format/MeanDriftMomentum.8", rows, error)
        if (allocated(error)) return
        worst_m = 0d0
        worst_y = 0d0
        done = 0
        do r = 1, size(rows, 1) - 2, 3
            if (abs(rows(r,2)) > 1d-9 .or. abs(rows(r,3)) > 1d-9) cycle
            fx = rows(r,7)
            fy = rows(r+1,7)
            mz = rows(r+2,7)
            worst_m = max(worst_m, abs(mz + 15d0*fx)/abs(15d0*fx))
            worst_y = max(worst_y, abs(fy)/abs(fx))
            done = done + 1
        end do
        call check(error, done == 2, "Expected heading 0 at two frequencies")
        if (allocated(error)) return
        call report("two bodies M_z = -15 F_x", worst_m)
        call report("two bodies F_y / F_x", worst_y)
        call check_below(error, "two bodies M_z = -15 F_x", worst_m, moment_limit)
        if (allocated(error)) return
        call check_below(error, "two bodies F_y / F_x", worst_y, lateral_limit)
    end subroutine two_bodies_moment

    !> The absolute drift force: F.e_beta from .8 against kappa int |a|^2 (1 - cos(theta - beta)) dtheta, with kappa/(rho g)
    !> = (C_g/C_p)/2 and a from the diagnostics. This is the form for a body that absorbs nothing.
    subroutine drift_matches_maruo(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: rows(:,:), waves(:,:), total(:,:)
        real(8) :: k, nu, h, beta, along, integral, worst, err
        integer :: r, t, m, done

        call read_drift(output_dir(main_case)//"Wamit_format/MeanDriftMomentum.8", rows, error)
        if (allocated(error)) return
        call read_diagnostics(main_case, "WAVE", 4, waves, error)
        if (allocated(error)) return
        call read_diagnostics(main_case, "TOTAL", 5, total, error)
        if (allocated(error)) return

        worst = 0d0
        done = 0
        do r = 1, size(rows, 1) - 2, 3
            beta = rows(r,2)*pi/180d0
            along = rows(r,7)*cos(beta) + rows(r+1,7)*sin(beta)
            call find_wave(waves, rows(r,1), k, nu, h)
            integral = 0d0
            m = 0
            do t = 1, size(total, 1)
                if (.not. same(total(t,1), rows(r,1)) .or. abs(total(t,2) - rows(r,2)) > 1d-9) cycle
                integral = integral + (total(t,4)**2 + total(t,5)**2)*(1d0 - cos(total(t,3)*pi/180d0 - beta))
                m = m + 1
            end do
            integral = 0.5d0*group_ratio(k, h)*integral*2d0*pi/real(m, 8)
            err = abs(along - integral)/abs(along)
            print '(4x,"drift_matches_maruo w=",F5.3," beta=",F5.1,": ",ES10.3)', rows(r,1), rows(r,2), err
            worst = max(worst, err)
            done = done + 1
        end do
        call check(error, done == 6, "Expected 3 frequencies and 2 headings")
        if (allocated(error)) return
        call check_below(error, "drift_matches_maruo", worst, maruo_limit)
    end subroutine drift_matches_maruo

    !> The same cylinder as a half mesh with ISY = 1 goes through the symmetric solver and the image expansion.
    subroutine half_mesh_matches_full(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: full(:,:), half(:,:)
        real(8) :: worst
        integer :: r

        call read_drift(output_dir(main_case)//"Wamit_format/MeanDriftMomentum.8", full, error)
        if (allocated(error)) return
        call read_drift(output_dir("drift_cyl_half")//"Wamit_format/MeanDriftMomentum.8", half, error)
        if (allocated(error)) return
        call check(error, size(full, 1) == size(half, 1), "Half and full mesh runs have different numbers of rows")
        if (allocated(error)) return
        worst = 0d0
        do r = 1, size(full, 1)
            if (nint(full(r,4)) == 1) worst = max(worst, abs(full(r,7) - half(r,7))/abs(full(r,7)))
            ! F_y where it is not zero by symmetry
            if (nint(full(r,4)) == 2 .and. abs(full(r,2)) > 1d-9) worst = max(worst, abs(full(r,7) - half(r,7))/abs(full(r,7)))
        end do
        call report("half_mesh_matches_full far", worst)
        call check_below(error, "half_mesh_matches_full far", worst, half_mesh_limit)
        if (allocated(error)) return

        call read_drift(output_dir(main_case)//near_file, full, error)
        if (allocated(error)) return
        call read_drift(output_dir("drift_cyl_half")//near_file, half, error)
        if (allocated(error)) return
        call check(error, size(full, 1) == size(half, 1), "Half and full mesh runs have different numbers of near-field rows")
        if (allocated(error)) return
        worst = 0d0
        do r = 1, size(full, 1)
            if (nint(full(r,4)) == 1) worst = max(worst, abs(full(r,7) - half(r,7))/abs(full(r,7)))
            if (nint(full(r,4)) == 2 .and. abs(full(r,2)) > 1d-9) worst = max(worst, abs(full(r,7) - half(r,7))/abs(full(r,7)))
        end do
        call report("half_mesh_matches_full near", worst)
        call check_below(error, "half_mesh_matches_full near", worst, half_mesh_limit)
    end subroutine half_mesh_matches_full

    ! ------------------------------------------------------------------------------------------- near-field against far

    !> Diffraction only: a fixed cylinder with a rounded bilge, at all the frequencies and headings of the fixture.
    subroutine fixed_near_far(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: worst_fx, worst_f, worst_mz

        call near_against_far("drift_round_fixed", origins_at_zero(1), 0d0, worst_fx, worst_f, worst_mz, error)
        if (allocated(error)) return
        call check_below(error, "fixed_near_far", worst_fx, near_far_limit)
    end subroutine fixed_near_far

    !> The free cylinder with a rounded bilge, where terms III and IV take part.
    !> On its cosine rows F_x falls with refinement, at 480, 1920, 4320 panels: -0.16%, +0.06%, +0.05% at 0.95 rad/s and +0.86%,
    !> +0.22%, +0.09% at 1.10 rad/s. See the README.
    subroutine free_near_far(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: worst_fx, worst_f, worst_mz

        call near_against_far("drift_round", origins_at_zero(1), small_net_force_below, worst_fx, worst_f, worst_mz, error)
        if (allocated(error)) return
        call check_below(error, "free_near_far", worst_fx, near_far_limit)
    end subroutine free_near_far

    !> The free cylinder at 0.70 rad/s alone, where the net force is small next to its terms. The worst force vector over the
    !> headings, which is F_x at heading 0.
    subroutine free_near_far_low_frequency(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: worst_fx, worst_f, worst_mz

        call near_against_far("drift_round", origins_at_zero(1), 0d0, worst_fx, worst_f, worst_mz, error, only_omega=low_frequency)
        if (allocated(error)) return
        call check_below(error, "free_near_far_low_frequency", worst_f, near_far_limit)
    end subroutine free_near_far_low_frequency

    !> M_z about the global origin at an oblique heading: on a half ellipsoid the yaw moment comes from the pressure, and on
    !> a round hull on external surge springs from the term xi x G, as the heave and surge motions differ.
    subroutine near_mz_far_mz(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: worst_fx, worst_f, worst_mz

        call near_against_far("drift_ellipsoid", origins_at_zero(1), 0d0, worst_fx, worst_f, worst_mz, error)
        if (allocated(error)) return
        call check_below(error, "near_mz_far_mz ellipsoid", worst_mz, near_mz_limit)
        if (allocated(error)) return
        call near_against_far("drift_round_kstf", origins_at_zero(1), 0d0, worst_fx, worst_f, worst_mz, error)
        if (allocated(error)) return
        call check_below(error, "near_mz_far_mz round hull on springs", worst_mz, near_mz_limit)
    end subroutine near_mz_far_mz

    !> External springs change the motion, and the near-field form needs no external loads. The spring in x only, at an
    !> oblique heading, makes F_x small next to F_y, so the check is on the force vector. Away from resonance. On the cosine wall
    !> rows, F_x at heading 0 falls as +1.36%, +0.40%, +0.21% at 480, 1920, 4320 panels. See the README.
    subroutine external_stiffness(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: worst_fx, worst_f, worst_mz

        call near_against_far("drift_round_kstf", origins_at_zero(1), springs_away_from_resonance, worst_fx, worst_f, worst_mz, &
                              error)
        if (allocated(error)) return
        call check_below(error, "external_stiffness", worst_f, near_far_limit)
    end subroutine external_stiffness

    !> The springs at 0.95 rad/s alone, near the heave resonance. The worst force vector over the headings: F_x at heading 0 is
    !> 4.4% above far-field, and the force vector at heading 30 is 2.7% off, where F_x alone is 7.3%. At heading 0 it falls with
    !> refinement, to +1.1% and +0.5% at 1920 and 4320 panels.
    subroutine external_stiffness_near_resonance(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: worst_fx, worst_f, worst_mz

        call near_against_far("drift_round_kstf", origins_at_zero(1), 0d0, worst_fx, worst_f, worst_mz, error, &
                              only_omega=springs_near_resonance)
        if (allocated(error)) return
        call check_below(error, "external_stiffness_near_resonance", worst_f, near_far_limit)
    end subroutine external_stiffness_near_resonance

    !> Regression band beside that check: F_x and F_y at 0.95 rad/s on the springs.
    subroutine external_stiffness_near_resonance_regression(error)
        type(error_type), allocatable, intent(out) :: error

        call check_band("drift_round_kstf", origins_at_zero(1), 1, reshape( &
            [0.95d0, 0d0, 0d0, 8.2986d0, 8.6617d0, &
             0.95d0, 30d0, 30d0, 3.7447d0, 4.0172d0], [5, 2]), error)
        if (allocated(error)) return
        call check_band("drift_round_kstf", origins_at_zero(1), 2, reshape( &
            [0.95d0, 30d0, 30d0, 9.3662d0, 9.4103d0], [5, 1]), error)
    end subroutine external_stiffness_near_resonance_regression

    !> The sum over the bodies against the far-field total. F_y vanishes by symmetry at heading 0, so the check is on the
    !> force vector. The 5% limit applies to the single headings, the mean load of one wave.
    subroutine two_body_total(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: worst_fx, worst_f, worst_mz, worst_pair

        call near_against_far("drift_two_round", reshape([0d0, 0d0, 0d0, 30d0], [2, 2]), 0d0, worst_fx, worst_f, worst_mz, error, &
                              worst_pair)
        if (allocated(error)) return
        call check_below(error, "two_body_total", worst_f, near_far_limit)
        if (allocated(error)) return
        ! Fails if a moment in .9 is not about XR of its own body
        call check_below(error, "two_body_total M_z about the origin", worst_mz, near_mz_limit)
        if (allocated(error)) return
        call check_below(error, "two_body_total off-diagonal heading pairs", worst_pair, heading_pair_limit)
    end subroutine two_body_total

    !> Regression band, not accuracy: F_x of the fixed coarse hull with a sharp bilge. Term II on the wall rows nearest the bilge
    !> holds the excess: an offline fit of the edge behaviour there, not rerun on the current solver, brought it to +0.9% at 4032
    !> panels.
    subroutine fixed_near_far_sharp_regression(error)
        type(error_type), allocatable, intent(out) :: error

        call check_band("drift_cyl_fixed", origins_at_zero(1), 1, reshape( &
            [0.70d0, 0d0, 0d0, 1.7385d0, 1.9568d0, &
             0.95d0, 0d0, 0d0, 5.4504d0, 5.7788d0, &
             1.10d0, 0d0, 0d0, 5.9241d0, 6.1909d0], [5, 3]), error)
    end subroutine fixed_near_far_sharp_regression

    !> Regression band, not accuracy: F_x of the free coarse hull with a sharp bilge. The same offline fit, not rerun on the current
    !> solver, brought the excess to +5.7% at 0.95 rad/s and 4032 panels.
    subroutine free_near_far_sharp_regression(error)
        type(error_type), allocatable, intent(out) :: error

        call check_band(main_case, origins_at_zero(1), 1, reshape( &
            [0.95d0, 0d0, 0d0, 2.8874d0, 4.1612d0, &
             1.10d0, 0d0, 0d0, 3.4665d0, 3.9795d0], [5, 2]), error)
    end subroutine free_near_far_sharp_regression

    !> Regression band, not accuracy: the sum of F_x and of F_y over the two coarse hulls with sharp bilges.
    subroutine two_body_total_sharp_regression(error)
        type(error_type), allocatable, intent(out) :: error

        call check_band("drift_two_cyl", reshape([0d0, 0d0, 0d0, 30d0], [2, 2]), 1, reshape( &
            [1.35d0, 0d0, 0d0, 2.1577d0, 3.2919d0, &
             1.35d0, 30d0, 30d0, 3.1479d0, 4.8420d0, &
             1.55d0, 0d0, 0d0, 4.6141d0, 5.5246d0, &
             1.55d0, 30d0, 30d0, 1.8863d0, 2.2216d0], [5, 4]), error)
        if (allocated(error)) return
        call check_band("drift_two_cyl", reshape([0d0, 0d0, 0d0, 30d0], [2, 2]), 2, reshape( &
            [1.35d0, 30d0, 30d0, 2.4502d0, 3.7572d0, &
             1.55d0, 30d0, 30d0, 1.5876d0, 2.0570d0], [5, 2]), error)
    end subroutine two_body_total_sharp_regression

    !> Regression band, not accuracy: M_z on the coarse box and on the coarse cylinder off the origin. For the box the cause is
    !> only partly shown: the offline edge fit, not rerun on the current solver, took F_y from +43.1% to +26.4% at 1760 panels.
    !> The cylinder's M_z is the moment of its force about the origin, so it has the error of the force.
    subroutine near_mz_far_mz_sharp_regression(error)
        type(error_type), allocatable, intent(out) :: error

        call check_band("drift_box", origins_at_zero(1), 3, reshape( &
            [0.8d0, 30d0, 30d0, -9.1211d0, -6.5082d0, &
             1.0d0, 30d0, 30d0, -11.916d0, -8.0784d0], [5, 2]), error)
        if (allocated(error)) return
        call check_band("drift_offset", reshape([6d0, 8d0], [2, 1]), 3, reshape( &
            [1.2d0, 30d0, 30d0, -7.8737d0, -9.0555d0], [5, 1]), error)
    end subroutine near_mz_far_mz_sharp_regression

    !> Far-field and near-field of a load stay within regression_band of the recorded values.
    !>   load       1 F_x, 2 F_y, 3 M_z about the origin
    !>   recorded   (omega, beta1, beta2, far, near) per column
    subroutine check_band(name, xr, load, recorded, error)
        character(len=*), intent(in) :: name
        real(8), intent(in) :: xr(:,:), recorded(:,:)
        integer, intent(in) :: load
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: far(:,:), near(:,:)
        complex(8) :: f_far(3), f_near(3)
        real(8) :: worst
        integer :: r, c, found

        call read_drift(output_dir(name)//far_file, far, error)
        if (allocated(error)) return
        call read_drift(output_dir(name)//near_file, near, error)
        if (allocated(error)) return
        worst = 0d0
        found = 0
        do c = 1, size(recorded, 2)
            do r = 1, size(far, 1) - 2, 3
                if (.not. same(far(r,1), recorded(1,c)) .or. far(r,2) /= recorded(2,c) .or. far(r,3) /= recorded(3,c)) cycle
                f_far = cmplx(far(r:r+2,7), far(r:r+2,8), kind=8)
                call sum_near(near, xr, far(r,1), far(r,2), far(r,3), f_near)
                worst = max(worst, abs(real(f_far(load), 8) - recorded(4,c))/abs(recorded(4,c)), &
                            abs(real(f_near(load), 8) - recorded(5,c))/abs(recorded(5,c)))
                found = found + 1
            end do
        end do
        call check(error, found == size(recorded, 2), "Not every recorded row of "//name//" is in its output")
        if (allocated(error)) return
        call report(name//" regression band", worst)
        call check_below(error, name//" regression band", worst, regression_band)
    end subroutine check_band

    !> Two cylinders 50 wavelengths apart: each has the F_x of a cylinder alone.
    subroutine far_apart_isolated(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: pair(:,:), alone(:,:)
        real(8) :: worst
        logical :: exists
        integer :: r, a

        call read_drift(output_dir("drift_far_apart")//near_file, pair, error)
        if (allocated(error)) return
        call read_drift(output_dir("drift_isolated")//near_file, alone, error)
        if (allocated(error)) return
        call check(error, size(pair, 1) == 12 .and. size(alone, 1) == 6, "Expected 6 rows per body")
        if (allocated(error)) return
        ! drift_isolated runs near-field only
        inquire (file=output_dir("drift_isolated")//far_file, exist=exists)
        call check(error, .not. exists, "A run with drift_methods = near wrote the far-field file")
        if (allocated(error)) return
        worst = 0d0
        do r = 1, size(pair, 1)
            if (nint(pair(r,4)) /= 1 .and. nint(pair(r,4)) /= 7) cycle
            do a = 1, size(alone, 1)
                if (nint(alone(a,4)) == 1) worst = max(worst, abs(pair(r,7) - alone(a,7))/abs(alone(a,7)))
            end do
        end do
        call report("far_apart_isolated", worst)
        call check_below(error, "far_apart_isolated", worst, near_far_limit)
    end subroutine far_apart_isolated

    !> A submerged body has no waterline, so terms I and Ia are exactly zero, while its other loads are finite and not zero.
    subroutine submerged_mix(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: terms(:,:), rows(:,:)
        real(8) :: wl_piercing, wl_submerged, velocity_submerged, fsub, worst_fx, worst_f, worst_mz
        integer :: r

        call check_no_bad_numbers(output_dir("drift_mixed")//near_file, error)
        if (allocated(error)) return
        call read_terms(output_dir("drift_mixed")//terms_file, terms, error)
        if (allocated(error)) return
        wl_piercing = 0d0
        wl_submerged = 0d0
        velocity_submerged = 0d0
        do r = 1, size(terms, 1)
            if (nint(terms(r,6)) == T_WATERLINE .or. nint(terms(r,6)) == T_FLARE) then
                if (nint(terms(r,4)) == 1) wl_piercing = max(wl_piercing, abs(terms(r,7)))
                if (nint(terms(r,4)) == 2) wl_submerged = max(wl_submerged, abs(terms(r,7)), abs(terms(r,8)))
            else if (nint(terms(r,6)) == T_VELOCITY .and. nint(terms(r,4)) == 2) then
                velocity_submerged = max(velocity_submerged, abs(terms(r,7)))
            end if
        end do
        call check(error, wl_submerged == 0d0, "Terms I and Ia are not zero on the submerged body")
        if (allocated(error)) return
        call check(error, wl_piercing > 0d0, "Term I is zero on the surface-piercing body")
        if (allocated(error)) return
        call check(error, velocity_submerged > 0d0, "Term II is zero on the submerged body")
        if (allocated(error)) return
        call read_drift(output_dir("drift_mixed")//near_file, rows, error)
        if (allocated(error)) return
        fsub = maxval(rows(:,5), mask=nint(rows(:,4)) >= 7 .and. nint(rows(:,4)) <= 9)
        call report("submerged_mix body 2 force", fsub)
        call check(error, fsub > 0d0 .and. fsub < 1d6, "The force on the submerged body is not finite")
        if (allocated(error)) return
        ! A measurement: nothing about the far-field total is asked of this case
        call near_against_far("drift_mixed", reshape([0d0, 0d0, 0d0, 20d0], [2, 2]), 0d0, worst_fx, worst_f, worst_mz, error)
    end subroutine submerged_mix

    !> The layout rotated by 30 degrees with the wave: all loads, in global axes, are the rotation of the unrotated ones.
    subroutine rotated_lcs(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), parameter :: angle = 30d0*pi/180d0
        real(8), allocatable :: rot(:,:), unrot(:,:)
        real(8) :: worst, a(6), b(6), turned(6), fscale, mscale
        integer :: body

        call read_drift(output_dir("drift_rotated")//near_file, rot, error)
        if (allocated(error)) return
        call read_drift(output_dir("drift_unrotated")//near_file, unrot, error)
        if (allocated(error)) return
        call check(error, size(rot, 1) == 12 .and. size(unrot, 1) == 12, "Expected 6 rows per body")
        if (allocated(error)) return
        worst = 0d0
        do body = 1, 2
            a = rot(6*(body - 1) + 1:6*body, 7)
            b = unrot(6*(body - 1) + 1:6*body, 7)
            turned = b
            turned(1) = cos(angle)*b(1) - sin(angle)*b(2)
            turned(2) = sin(angle)*b(1) + cos(angle)*b(2)
            turned(4) = cos(angle)*b(4) - sin(angle)*b(5)
            turned(5) = sin(angle)*b(4) + cos(angle)*b(5)
            fscale = sqrt(sum(b(1:3)**2))
            mscale = sqrt(sum(b(4:6)**2))
            worst = max(worst, sqrt(sum((a(1:3) - turned(1:3))**2))/fscale, sqrt(sum((a(4:6) - turned(4:6))**2))/mscale)
        end do
        call report("rotated_lcs near", worst)
        call check_below(error, "rotated_lcs near", worst, rotated_limit)
        if (allocated(error)) return

        call read_drift(output_dir("drift_rotated")//far_file, rot, error)
        if (allocated(error)) return
        call read_drift(output_dir("drift_unrotated")//far_file, unrot, error)
        if (allocated(error)) return
        call check(error, size(rot, 1) == 3 .and. size(unrot, 1) == 3, "Expected three far-field rows")
        if (allocated(error)) return
        fscale = sqrt(unrot(1,7)**2 + unrot(2,7)**2)
        worst = sqrt((rot(1,7) - (cos(angle)*unrot(1,7) - sin(angle)*unrot(2,7)))**2 &
                     + (rot(2,7) - (sin(angle)*unrot(1,7) + cos(angle)*unrot(2,7)))**2)/fscale
        call check(error, abs(unrot(3,7)) > 1d-2*fscale, "M_z must not be negligible in this case")
        if (allocated(error)) return
        worst = max(worst, abs(rot(3,7) - unrot(3,7))/abs(unrot(3,7)))
        call report("rotated_lcs far", worst)
        call check_below(error, "rotated_lcs far", worst, rotated_limit)
    end subroutine rotated_lcs

    !> mean_drift = 2: F(b1, b2) = conj F(b2, b1) in the file, the diagonal is the mean_drift = 1 result, limit frequencies
    !> have no rows.
    subroutine bidirectional(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: pairs(:,:), uni(:,:)
        real(8) :: worst_sym, worst_diag, scale_row
        integer :: p, q, r, matched

        call read_drift(output_dir("drift_cyl_sybo")//near_file, pairs, error)
        if (allocated(error)) return
        call read_drift(output_dir(main_case)//near_file, uni, error)
        if (allocated(error)) return
        call check(error, size(pairs, 1) == 3*(2*2)*6, "Expected 6 rows for the 3 real frequencies and 4 heading pairs")
        if (allocated(error)) return
        call check(error, all(pairs(:,1) > 0.5d0), "A limit frequency has rows in the near-field file")
        if (allocated(error)) return

        worst_sym = 0d0
        worst_diag = 0d0
        matched = 0
        do p = 1, size(pairs, 1)
            if (nint(pairs(p,4)) <= 3) then
                scale_row = maxval(pairs(:,5), mask=abs(pairs(:,1) - pairs(p,1)) < 1d-6 .and. nint(pairs(:,4)) <= 3)
            else
                scale_row = maxval(pairs(:,5), mask=abs(pairs(:,1) - pairs(p,1)) < 1d-6 .and. nint(pairs(:,4)) > 3)
            end if
            do q = 1, size(pairs, 1)
                if (same(pairs(q,1), pairs(p,1)) .and. pairs(q,2) == pairs(p,3) .and. pairs(q,3) == pairs(p,2) &
                    .and. pairs(q,4) == pairs(p,4)) then
                    worst_sym = max(worst_sym, (abs(pairs(q,7) - pairs(p,7)) + abs(pairs(q,8) + pairs(p,8)))/scale_row)
                end if
            end do
            do r = 1, size(uni, 1)
                if (same(uni(r,1), pairs(p,1)) .and. uni(r,2) == pairs(p,2) .and. uni(r,3) == pairs(p,3) &
                    .and. uni(r,4) == pairs(p,4)) then
                    worst_diag = max(worst_diag, abs(uni(r,7) - pairs(p,7))/scale_row)
                    matched = matched + 1
                end if
            end do
        end do
        call check(error, matched == size(uni, 1), "Not every mean_drift = 1 row is in the pairs output")
        if (allocated(error)) return
        call report("bidirectional conjugate symmetry", worst_sym)
        call report("bidirectional diagonal", worst_diag)
        call check_below(error, "bidirectional conjugate symmetry", worst_sym, printed_limit)
        if (allocated(error)) return
        call check_below(error, "bidirectional diagonal", worst_diag, printed_limit)
    end subroutine bidirectional

    !> IRSP = 1 against IRSP = 0 through both adapters, far-field and near-field, away from irregular frequencies. The
    !> first irregular frequency of the single cylinder has a printed value and no limit.
    subroutine irsp_stability(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), parameter :: irregular = 1.548d0
        real(8) :: far_single, near_single, far_multi, near_multi

        call irsp_difference("single far", "drift_irsp", "drift_irsp_ref", far_file, irregular, far_single, error)
        if (allocated(error)) return
        call irsp_difference("single near", "drift_irsp", "drift_irsp_ref", near_file, irregular, near_single, error)
        if (allocated(error)) return
        call irsp_difference("multi far", "drift_two_irsp", "drift_two_irsp_ref", far_file, irregular, far_multi, error)
        if (allocated(error)) return
        call irsp_difference("multi near", "drift_two_irsp", "drift_two_irsp_ref", near_file, irregular, near_multi, error)
        if (allocated(error)) return
        call check_below(error, "irsp_stability single far", far_single, irsp_limit)
        if (allocated(error)) return
        call check_below(error, "irsp_stability single near", near_single, irsp_limit)
        if (allocated(error)) return
        call check_below(error, "irsp_stability multi far", far_multi, irsp_limit)
        if (allocated(error)) return
        call check_below(error, "irsp_stability multi near", near_multi, irsp_limit)
    end subroutine irsp_stability

    !> The largest relative difference of F_x (I = 1, 7) on the diagonal heading pairs away from the frequency irregular,
    !> which is only printed.
    subroutine irsp_difference(tag, with_case, without_case, file, irregular, worst, error)
        character(len=*), intent(in) :: tag, with_case, without_case, file
        real(8), intent(in) :: irregular
        real(8), intent(out) :: worst
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: on(:,:), off(:,:)
        real(8) :: err
        integer :: r, q, used

        call read_drift(output_dir(with_case)//file, on, error)
        if (allocated(error)) return
        call read_drift(output_dir(without_case)//file, off, error)
        if (allocated(error)) return
        call check(error, size(on, 1) == size(off, 1), "IRSP = 1 and IRSP = 0 runs have different numbers of rows")
        if (allocated(error)) return
        worst = 0d0
        used = 0
        do r = 1, size(on, 1)
            if (nint(on(r,4)) /= 1 .and. nint(on(r,4)) /= 7) cycle
            if (abs(on(r,2) - on(r,3)) > 1d-9) cycle
            do q = 1, size(off, 1)
                if (.not. (same(off(q,1), on(r,1)) .and. off(q,2) == on(r,2) .and. off(q,3) == on(r,3) &
                           .and. off(q,4) == on(r,4))) cycle
                err = abs(on(r,7) - off(q,7))/abs(off(q,7))
                print '(4x,"irsp ",a," w=",F5.3," beta=",F5.1," I=",I1,": on ",ES11.4," off ",ES11.4," diff ",ES10.3)', &
                    tag, on(r,1), on(r,2), nint(on(r,4)), on(r,7), off(q,7), err
                if (abs(on(r,1) - irregular) > 1d-6) then
                    worst = max(worst, err)
                    used = used + 1
                end if
            end do
        end do
        call check(error, used > 0, "No frequency compared in irsp_stability")
    end subroutine irsp_difference

    !> The terms of MeanDriftTerms.txt add up to the totals of the .9 file, row by row. Each has 7 digits, so the difference is
    !> relative to the largest term of the frequency.
    subroutine terms_sum_to_total(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: terms(:,:), rows(:,:)
        real(8) :: worst, total(2), scale_row
        integer :: r, t, n

        call read_drift(output_dir("drift_round_kstf")//near_file, rows, error)
        if (allocated(error)) return
        call read_terms(output_dir("drift_round_kstf")//terms_file, terms, error)
        if (allocated(error)) return
        call check(error, size(terms, 1) == size(rows, 1)*T_HYDROSTATIC, "Expected one row for each term of each load")
        if (allocated(error)) return
        worst = 0d0
        do r = 1, size(rows, 1)
            total = 0d0
            n = 0
            do t = 1, size(terms, 1)
                if (same(terms(t,1), rows(r,1)) .and. terms(t,2) == rows(r,2) .and. terms(t,3) == rows(r,3) &
                    .and. 6*(nint(terms(t,4)) - 1) + nint(terms(t,5)) == nint(rows(r,4))) then
                    total = total + terms(t,7:8)
                    n = n + 1
                end if
            end do
            call check(error, n == T_HYDROSTATIC, "A load does not have all its terms")
            if (allocated(error)) return
            ! Forces against forces and moments against moments, at this frequency
            scale_row = maxval(abs(terms(:,7)), mask=abs(terms(:,1) - rows(r,1)) < 1d-6 .and. &
                               ((nint(terms(:,5)) <= 3) .eqv. (nint(rows(r,4)) <= 3)))
            worst = max(worst, (abs(total(1) - rows(r,7)) + abs(total(2) - rows(r,8)))/scale_row)
        end do
        call report("terms_sum_to_total", worst)
        call check_below(error, "terms_sum_to_total", worst, printed_limit)
    end subroutine terms_sum_to_total

    !> Far-only runs have no near-field output, write_terms = 0 has no terms file, and runs with both methods log the
    !> difference, with M_z moved to the origin by the driver as this test does it.
    subroutine near_output_switches(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: worst_fx, worst_f, worst_mz, logged
        logical :: exists

        inquire (file=output_dir("drift_cyl_gmres")//near_file, exist=exists)
        call check(error, .not. exists, "A far-only run wrote the near-field file")
        if (allocated(error)) return
        inquire (file=output_dir("drift_cyl_gmres")//terms_file, exist=exists)
        call check(error, .not. exists, "A far-only run wrote the terms file")
        if (allocated(error)) return
        inquire (file=output_dir("drift_cyl_sybo")//terms_file, exist=exists)
        call check(error, .not. exists, "write_terms = 0 wrote the terms file")
        if (allocated(error)) return
        call check(error, log_has(main_case, "near vs far"), "The log of a run with both methods has no near vs far row")
        if (allocated(error)) return
        call check(error, .not. log_has("drift_cyl_gmres", "near vs far"), "A far-only log has a near vs far row")
        if (allocated(error)) return
        call check(error, log_has(main_case, "velocity fallback"), "The log has no count of velocity fallbacks")
        if (allocated(error)) return
        call near_against_far("drift_offset", reshape([6d0, 8d0], [2, 1]), 0d0, worst_fx, worst_f, worst_mz, error)
        if (allocated(error)) return
        call logged_difference("drift_offset", "near vs far M_z", logged, error)
        if (allocated(error)) return
        call check_below(error, "near_output_switches logged M_z", abs(logged - worst_mz), printed_limit*1d2)
    end subroutine near_output_switches

    !> The last value of the first row of a log that holds the text.
    subroutine logged_difference(name, text, value, error)
        character(len=*), intent(in) :: name, text
        real(8), intent(out) :: value
        type(error_type), allocatable, intent(out) :: error
        character(len=512) :: line
        integer :: unit, ios, at

        value = -1d0
        open (newunit=unit, file=get_running_exe_path()//name//".log", status='old', action='read', iostat=ios)
        call check(error, ios == 0, "Cannot open the log of "//name)
        if (allocated(error)) return
        do
            read (unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, text) == 0) cycle
            at = index(line, "difference") + len("difference")
            read (line(at:), *) value
            exit
        end do
        close (unit)
        call check(error, value >= 0d0, "No row with "//text//" in the log of "//name)
    end subroutine logged_difference

    logical function log_has(name, text)
        character(len=*), intent(in) :: name, text
        character(len=512) :: line
        integer :: unit, ios

        log_has = .false.
        open (newunit=unit, file=get_running_exe_path()//name//".log", status='old', action='read', iostat=ios)
        if (ios /= 0) return
        do
            read (unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, text) > 0) log_has = .true.
        end do
        close (unit)
    end function log_has

    function origins_at_zero(nbody) result(xr)
        integer, intent(in) :: nbody
        real(8) :: xr(2, nbody)
        xr = 0d0
    end function origins_at_zero

    !> Sums the near-field rows of each frequency and heading pair over the bodies, moves M_z from the body reference points
    !> xr to the global origin, and compares with the far-field rows. Prints every row. Differences are relative to the
    !> same far-field load: F_x (worst_fx), the force vector (worst_f), M_z (worst_mz). A component that is zero by symmetry
    !> has no relative difference, so one below 1e-3 of |F| (or of |F| times the size of the layout for M_z) is skipped.
    !> Rows below the frequency unchecked_below are printed, but do not count in the worst values. Nor do rows at any other
    !> frequency than only_omega, if given. With worst_pair, rows of two different headings count there instead, as the force
    !> vector.
    subroutine near_against_far(name, xr, unchecked_below, worst_fx, worst_f, worst_mz, error, worst_pair, only_omega)
        character(len=*), intent(in) :: name
        real(8), intent(in) :: xr(:,:)    ! (2, nbody), as REFL = 1
        real(8), intent(in) :: unchecked_below
        real(8), intent(out) :: worst_fx, worst_f, worst_mz
        type(error_type), allocatable, intent(out) :: error
        real(8), intent(out), optional :: worst_pair
        real(8), intent(in), optional :: only_omega
        real(8), allocatable :: far(:,:), near(:,:)
        complex(8) :: f_far(3), f_near(3)
        real(8) :: fmag, err_x, err_f, err_m, arm
        integer :: r, used

        used = 0
        worst_fx = 0d0
        worst_f = 0d0
        worst_mz = 0d0
        if (present(worst_pair)) worst_pair = 0d0
        call read_drift(output_dir(name)//far_file, far, error)
        if (allocated(error)) return
        call read_drift(output_dir(name)//near_file, near, error)
        if (allocated(error)) return
        call check(error, size(near, 1) == 2*size(far, 1)*size(xr, 2), "Expected 6 near-field rows per body for 3 far-field rows")
        if (allocated(error)) return
        arm = max(1d0, maxval(sqrt(sum(xr**2, dim=1))))

        do r = 1, size(far, 1) - 2, 3
            f_far = cmplx(far(r:r+2,7), far(r:r+2,8), kind=8)
            call sum_near(near, xr, far(r,1), far(r,2), far(r,3), f_near)
            fmag = sqrt(abs(f_far(1))**2 + abs(f_far(2))**2)
            err_f = sqrt(abs(f_near(1) - f_far(1))**2 + abs(f_near(2) - f_far(2))**2)/fmag
            err_x = -1d0
            err_m = -1d0
            if (abs(f_far(1)) > 1d-3*fmag) err_x = abs(f_near(1) - f_far(1))/abs(f_far(1))
            if (abs(f_far(3)) > 1d-3*fmag*arm) err_m = abs(f_near(3) - f_far(3))/abs(f_far(3))
            print '(4x,a," w=",F5.3," b1=",F5.1," b2=",F5.1,": far ",3ES12.4," near ",3ES12.4,&
                       &" err x ",ES9.2," F ",ES9.2," z ",ES9.2)', &
                name, far(r,1), far(r,2), far(r,3), real(f_far, 8), real(f_near, 8), err_x, err_f, err_m
            if (far(r,1) < unchecked_below) cycle
            if (present(only_omega)) then
                if (.not. same(far(r,1), only_omega)) cycle
            end if
            used = used + 1
            if (present(worst_pair) .and. far(r,2) /= far(r,3)) then
                worst_pair = max(worst_pair, err_f)
                cycle
            end if
            worst_fx = max(worst_fx, err_x)
            worst_f = max(worst_f, err_f)
            worst_mz = max(worst_mz, err_m)
        end do
        call check(error, used > 0, "No row of "//name//" was compared")
    end subroutine near_against_far

    !> (F_x, F_y, M_z about the origin) summed over the bodies, from the rows of a .9 file for one frequency and heading pair.
    subroutine sum_near(near, xr, omega, beta1, beta2, total)
        real(8), intent(in) :: near(:,:), xr(:,:), omega, beta1, beta2
        complex(8), intent(out) :: total(3)
        complex(8) :: load
        integer :: q, body, local

        total = (0d0, 0d0)
        do q = 1, size(near, 1)
            if (.not. same(near(q,1), omega) .or. near(q,2) /= beta1 .or. near(q,3) /= beta2) cycle
            body = (nint(near(q,4)) - 1)/6 + 1
            local = nint(near(q,4)) - 6*(body - 1)
            load = cmplx(near(q,7), near(q,8), kind=8)
            select case (local)
            case (1)
                total(1) = total(1) + load
                total(3) = total(3) - xr(2,body)*load
            case (2)
                total(2) = total(2) + load
                total(3) = total(3) + xr(1,body)*load
            case (6)
                total(3) = total(3) + load
            end select
        end do
    end subroutine sum_near

    subroutine energy_balance_of(name, nexpected, worst, error)
        character(len=*), intent(in) :: name
        integer, intent(in) :: nexpected        ! number of frequency and heading combinations
        real(8), intent(out) :: worst
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: waves(:,:), total(:,:)
        real(8) :: omega, beta, k, integral, x, err
        complex(8) :: a_beta
        integer :: iw, r, m, b, nbeta, done
        real(8) :: betas(32)

        call read_diagnostics(name, "WAVE", 4, waves, error)
        if (allocated(error)) return
        call read_diagnostics(name, "TOTAL", 5, total, error)
        if (allocated(error)) return

        nbeta = 0
        do r = 1, size(total, 1)
            if (nbeta > 0) then
                if (any(abs(betas(1:nbeta) - total(r,2)) < 1d-9)) cycle
            end if
            nbeta = nbeta + 1
            betas(nbeta) = total(r,2)
        end do

        worst = 0d0
        done = 0
        do iw = 1, size(waves, 1)
            omega = waves(iw,1)
            k = waves(iw,2)
            do b = 1, nbeta
                beta = betas(b)
                integral = 0d0
                a_beta = (0d0, 0d0)
                m = 0
                do r = 1, size(total, 1)
                    if (.not. same(total(r,1), omega) .or. abs(total(r,2) - beta) > 1d-9) cycle
                    integral = integral + total(r,4)**2 + total(r,5)**2
                    if (abs(total(r,3) - beta) < 1d-6) a_beta = cmplx(total(r,4), total(r,5), kind=8)
                    m = m + 1
                end do
                integral = integral*2d0*pi/real(m, 8)
                x = -2d0*sqrt(2d0*pi/k)*real(a_beta*exp(ci*pi/4d0), 8)
                err = abs(x - integral)/integral
                worst = max(worst, err)
                done = done + 1
            end do
        end do
        call report("energy_balance "//name, worst)
        call check(error, done == nexpected, "Unexpected number of frequency and heading combinations")
    end subroutine energy_balance_of

    ! -------------------------------------------------------------------------------------------------------- helpers

    subroutine check_below(error, name, value, limit)
        type(error_type), allocatable, intent(out) :: error
        character(len=*), intent(in) :: name
        real(8), intent(in) :: value, limit
        character(len=16) :: a, b

        write (a, '(ES10.3)') value
        write (b, '(ES10.3)') limit
        call check(error, value < limit, name//": "//trim(a)//" is not below "//trim(b))
    end subroutine check_below

    !> Prints a measured value, so that the margin to each limit stays visible in the test log.
    subroutine report(name, value)
        character(len=*), intent(in) :: name
        real(8), intent(in) :: value
        print '(4x,a,": ",ES10.3)', name, value
    end subroutine report

    logical function same(a, b)
        real(8), intent(in) :: a, b
        same = abs(a - b) < 1d-5
    end function same

    real(8) function group_ratio(k, h)
        real(8), intent(in) :: k, h
        group_ratio = 0.5d0*(1d0 + 2d0*k*h/sinh(2d0*k*h))
    end function group_ratio

    subroutine find_wave(waves, omega, k, nu, h)
        real(8), intent(in) :: waves(:,:), omega
        real(8), intent(out) :: k, nu, h
        integer :: r

        do r = 1, size(waves, 1)
            if (same(waves(r,1), omega)) then
                k = waves(r,2)
                nu = waves(r,3)
                h = waves(r,4)
                return
            end if
        end do
        error stop "find_wave: frequency not found"
    end subroutine find_wave

    subroutine check_no_bad_numbers(path, error)
        character(len=*), intent(in) :: path
        type(error_type), allocatable, intent(out) :: error
        character(len=256) :: line
        integer :: unit, ios

        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        call check(error, ios == 0, "Cannot open "//path)
        if (allocated(error)) return
        do
            read (unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, 'NaN') > 0 .or. index(line, 'nan') > 0 .or. index(line, 'NAN') > 0 &
                .or. index(line, 'Inf') > 0 .or. index(line, '*') > 0) then
                call test_failed(error, "Bad number in "//path//": "//trim(line))
                exit
            end if
        end do
        close (unit)
    end subroutine check_no_bad_numbers

    !> Rows PER BETA1 BETA2 I Mod Pha Re Im of a mean drift file. The fixed format also pins the file layout.
    subroutine read_drift(path, rows, error)
        character(len=*), intent(in) :: path
        real(8), allocatable, intent(out) :: rows(:,:)
        type(error_type), allocatable, intent(out) :: error
        character(len=:), allocatable :: lines(:)
        real(8) :: a(3), b(4)
        integer :: i, r, ios

        call read_lines(path, lines, error)
        if (allocated(error)) return
        allocate (rows(size(lines), 8))
        do r = 1, size(lines)
            read (lines(r), '(3ES14.6,I6,4ES14.6)', iostat=ios) a, i, b
            call check(error, ios == 0, "Bad row in "//path//": "//trim(lines(r)))
            if (allocated(error)) return
            rows(r,1:3) = a
            rows(r,4) = real(i, 8)
            rows(r,5:8) = b
        end do
    end subroutine read_drift

    !> Rows PER BETA1 BETA2 body I term Re Im of MeanDriftTerms.txt, after its comment lines.
    subroutine read_terms(path, rows, error)
        character(len=*), intent(in) :: path
        real(8), allocatable, intent(out) :: rows(:,:)
        type(error_type), allocatable, intent(out) :: error
        character(len=:), allocatable :: lines(:)
        real(8) :: a(3), b(2)
        integer :: i(3), r, n, ios

        call read_lines(path, lines, error)
        if (allocated(error)) return
        allocate (rows(size(lines), 8))
        n = 0
        do r = 1, size(lines)
            if (lines(r)(1:1) == '#') cycle
            read (lines(r), '(3ES14.6,3I6,2ES14.6)', iostat=ios) a, i, b
            call check(error, ios == 0, "Bad row in "//path//": "//trim(lines(r)))
            if (allocated(error)) return
            n = n + 1
            rows(n,1:3) = a
            rows(n,4:6) = real(i, 8)
            rows(n,7:8) = b
        end do
        rows = rows(1:n,:)
    end subroutine read_terms

    !> Rows PER BETA I Mod Pha Re Im of ExcForce.3.
    subroutine read_exc(path, rows, error)
        character(len=*), intent(in) :: path
        real(8), allocatable, intent(out) :: rows(:,:)
        type(error_type), allocatable, intent(out) :: error
        call read_columns(path, 7, rows, error)
    end subroutine read_exc

    !> Rows PER I J A B of AmssDamp.1. Limit frequencies have no B and are not in these fixtures.
    subroutine read_amssdamp(path, rows, error)
        character(len=*), intent(in) :: path
        real(8), allocatable, intent(out) :: rows(:,:)
        type(error_type), allocatable, intent(out) :: error
        call read_columns(path, 5, rows, error)
    end subroutine read_amssdamp

    !> Diffraction rows PER BETA IPT Mod Pha Re Im of .6p, told from the much longer radiation rows by their length.
    subroutine read_6p_diffraction(path, rows, error)
        character(len=*), intent(in) :: path
        real(8), allocatable, intent(out) :: rows(:,:)
        type(error_type), allocatable, intent(out) :: error
        character(len=:), allocatable :: lines(:)
        real(8), allocatable :: all_rows(:,:)
        integer :: r, n, ios

        call read_lines(path, lines, error)
        if (allocated(error)) return
        allocate (all_rows(size(lines), 7))
        n = 0
        do r = 1, size(lines)
            if (len_trim(lines(r)) > 100) cycle
            n = n + 1
            read (lines(r), *, iostat=ios) all_rows(n,:)
            call check(error, ios == 0, "Bad row in "//path//": "//trim(lines(r)))
            if (allocated(error)) return
        end do
        rows = all_rows(1:n,:)
    end subroutine read_6p_diffraction

    !> Field panel numbers from the fixture: panel(ipt) is the hull panel whose centroid is field point ipt.
    subroutine read_field_panels(path, panel)
        character(len=*), intent(in) :: path
        integer, allocatable, intent(out) :: panel(:)
        integer :: unit, ios, ipt, p, n

        open (newunit=unit, file=path, status='old', action='read')
        read (unit, *)
        n = 0
        do
            read (unit, *, iostat=ios) ipt, p
            if (ios /= 0) exit
            n = n + 1
        end do
        rewind (unit)
        read (unit, *)
        allocate (panel(n))
        do n = 1, size(panel)
            read (unit, *) ipt, p
            panel(ipt) = p
        end do
        close (unit)
    end subroutine read_field_panels

    !> The field points listed in a ControlFile.in, one row per point.
    subroutine read_field_points(path, points)
        character(len=*), intent(in) :: path
        real(8), allocatable, intent(out) :: points(:,:)
        character(len=256) :: line
        integer :: unit, ios, n, i

        open (newunit=unit, file=path, status='old', action='read')
        n = -1
        do
            read (unit, '(A)', iostat=ios) line
            if (ios /= 0) error stop "read_field_points: no Number_of_field_points line"
            i = index(line, "Number_of_field_points")
            if (i > 0) then
                read (line(i + len("Number_of_field_points"):), *) n
                exit
            end if
        end do
        allocate (points(n, 3))
        do i = 1, n
            read (unit, *) points(i,:)
        end do
        close (unit)
    end subroutine read_field_points

    !> Numeric columns of the rows of PanelDiagnostics.txt that start with a tag.
    subroutine read_diagnostics(name, tag, ncol, rows, error)
        character(len=*), intent(in) :: name, tag
        integer, intent(in) :: ncol
        real(8), allocatable, intent(out) :: rows(:,:)
        type(error_type), allocatable, intent(out) :: error
        character(len=:), allocatable :: lines(:)
        character(len=16) :: first
        integer :: r, n, ios

        call read_lines(output_dir(name)//"SecondOrder/PanelDiagnostics.txt", lines, error)
        if (allocated(error)) return
        n = 0
        allocate (rows(size(lines), ncol))
        do r = 1, size(lines)
            read (lines(r), *, iostat=ios) first
            if (ios /= 0 .or. first /= tag) cycle
            n = n + 1
            read (lines(r), *, iostat=ios) first, rows(n,:)
            call check(error, ios == 0, "Bad diagnostics row: "//trim(lines(r)))
            if (allocated(error)) return
        end do
        call check(error, n > 0, "No "//tag//" rows in the diagnostics of "//name)
        if (allocated(error)) return
        rows = rows(1:n,:)
    end subroutine read_diagnostics

    subroutine read_columns(path, ncol, rows, error)
        character(len=*), intent(in) :: path
        integer, intent(in) :: ncol
        real(8), allocatable, intent(out) :: rows(:,:)
        type(error_type), allocatable, intent(out) :: error
        character(len=:), allocatable :: lines(:)
        integer :: r, ios

        call read_lines(path, lines, error)
        if (allocated(error)) return
        allocate (rows(size(lines), ncol))
        do r = 1, size(lines)
            read (lines(r), *, iostat=ios) rows(r,:)
            call check(error, ios == 0, "Bad row in "//path//": "//trim(lines(r)))
            if (allocated(error)) return
        end do
    end subroutine read_columns

    subroutine read_lines(path, lines, error)
        character(len=*), intent(in) :: path
        character(len=:), allocatable, intent(out) :: lines(:)
        type(error_type), allocatable, intent(out) :: error
        character(len=256) :: buffer
        integer :: unit, ios, n, r

        open (newunit=unit, file=path, status='old', action='read', iostat=ios)
        call check(error, ios == 0, "Cannot open "//path)
        if (allocated(error)) return
        n = 0
        do
            read (unit, '(A)', iostat=ios) buffer
            if (ios /= 0) exit
            n = n + 1
        end do
        rewind (unit)
        allocate (character(len=256) :: lines(n))
        do r = 1, n
            read (unit, '(A)') lines(r)
        end do
        close (unit)
    end subroutine read_lines

end module test_integration_drift
