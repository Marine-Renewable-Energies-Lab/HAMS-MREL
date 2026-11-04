!> End-to-end checks of the far-field mean drift: a real HAMS solution goes through the convention conversion
!> (section 5.2 of docs/superpowers/specs/2026-11-04-second-order-waves-design.md) and is compared with first-order
!> outputs written by the solver itself. Fixtures and the frequency choice: test/data/drift_cyl/README.md.
module test_integration_drift

    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
    use testing_utilities

    implicit none
    private
    public :: setup_drift, collect_tests_drift

    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: rho = 1000d0, g = 9.80665d0   ! as Const_mod
    complex(8), parameter :: ci = (0d0, 1d0)

    ! Acceptance limits, from the implementation plan (Task 5)
    real(8), parameter :: residual_limit = 1d-8
    real(8), parameter :: diffraction_limit = 1d-4
    real(8), parameter :: energy_limit = 1d-2         ! radiated energy, Haskind, energy balance
    real(8), parameter :: lateral_limit = 1d-3
    real(8), parameter :: gmres_limit = 1d-2
    ! Values in .8 have 7 significant digits, so same-run ratios cannot be checked more tightly than this
    real(8), parameter :: printed_limit = 1d-5

    character(len=*), parameter :: main_case = "drift_cyl"

contains

    subroutine setup_drift(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("run_drift_cyl", run_drift_cyl), &
            new_unittest("run_drift_cyl_sybo", run_drift_cyl_sybo), &
            new_unittest("run_drift_cyl_gmres", run_drift_cyl_gmres), &
            new_unittest("run_drift_offset", run_drift_offset), &
            new_unittest("run_drift_offset_refl2", run_drift_offset_refl2), &
            new_unittest("run_drift_two_cyl", run_drift_two_cyl) &
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
            new_unittest("two_bodies_smoke", two_bodies_smoke) &
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
        character(len=40) :: files(2)
        files = [character(len=40) :: "Wamit_format/MeanDriftMomentum.8", "SecondOrder/PanelDiagnostics.txt"]
        call run_case(main_case, files, error)
    end subroutine run_drift_cyl

    subroutine run_drift_cyl_sybo(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_cyl_sybo", drift_file_only(), error)
    end subroutine run_drift_cyl_sybo

    subroutine run_drift_cyl_gmres(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_cyl_gmres", drift_file_only(), error)
    end subroutine run_drift_cyl_gmres

    subroutine run_drift_offset(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_offset", drift_file_only(), error)
    end subroutine run_drift_offset

    subroutine run_drift_offset_refl2(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("drift_offset_refl2", drift_file_only(), error)
    end subroutine run_drift_offset_refl2

    subroutine run_drift_two_cyl(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=40) :: files(2)
        files = [character(len=40) :: "Wamit_format/MeanDriftMomentum.8", "SecondOrder/PanelDiagnostics.txt"]
        call run_case("drift_two_cyl", files, error)
    end subroutine run_drift_two_cyl

    function drift_file_only() result(files)
        character(len=40) :: files(1)
        files(1) = "Wamit_format/MeanDriftMomentum.8"
    end function drift_file_only

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

    !> conj(i w rho phi_D / (rho g A)) of the state against the on-body pressure the solver writes to .6p, in the units
    !> of that file (rho g A). The state's incident wave uses Const_mod's g, while the k behind VINP comes from the
    !> g = 9.807 in DISPERSION, so the two differ by 3.6e-5 of the incident pressure (spec 5.5, issue 3). Divided by the
    !> local pressure, this grows on the panels in the lee at depth, where the diffraction pressure is small.
    subroutine diffraction_matches_6p(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: panels(:,:), pressure(:,:)
        integer, allocatable :: field_panel(:)
        real(8) :: worst, worst_relative
        complex(8) :: phi_d, expected, written
        integer :: r, p, matched

        call read_diagnostics(main_case, "PANEL", 10, panels, error)
        if (allocated(error)) return
        call read_6p_diffraction(output_dir(main_case)//"Wamit_format/PressureElevation.6p", pressure, error)
        if (allocated(error)) return
        call read_field_panels(get_repo_root_path()//"test/data/"//main_case//"/field_panels.txt", field_panel)
        call check(error, size(pressure, 1) == 6*size(field_panel), "Unexpected number of diffraction rows in .6p")
        if (allocated(error)) return

        worst = 0d0
        worst_relative = 0d0
        matched = 0
        do r = 1, size(pressure, 1)
            do p = 1, size(panels, 1)
                if (.not. same(panels(p,1), pressure(r,1)) .or. .not. same(panels(p,2), pressure(r,2))) cycle
                if (nint(panels(p,3)) /= field_panel(nint(pressure(r,3)))) cycle
                phi_d = cmplx(panels(p,7), panels(p,8), kind=8)
                expected = conjg(ci*panels(p,1)*phi_d/g)
                written = cmplx(pressure(r,6), pressure(r,7), kind=8)
                worst = max(worst, abs(expected - written))
                worst_relative = max(worst_relative, abs(expected - written)/abs(written))
                matched = matched + 1
            end do
        end do
        call check(error, matched == size(pressure, 1), "Not every .6p row has a matching panel diagnostic")
        if (allocated(error)) return
        call report("diffraction_matches_6p", worst)
        call report("diffraction_matches_6p relative to the local pressure", worst_relative)
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
