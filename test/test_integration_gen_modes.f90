!> Integration test for generalized modes in DOF 1 to 3, on a small box written at run time.
!> A mode replaces one DOF's normal. So every entry of A, B and X that does not involve that DOF, and the field-point output
!> of the other modes and of diffraction, must equal the run with no mode. There is no reference for the mode's own numbers,
!> so those are checked for consistency: A and B symmetry, Haskind, and agreement with the same data placed in DOF 5, where
!> the columns read as the panel's direction are not touched.
module test_integration_gen_modes

    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed, skip_test
    use testing_utilities
    use IO, only : CreateDirectory

    implicit none
    private
    public :: setup_gen_modes, collect_tests_gen_modes

    integer, parameter :: NFREQ = 2, NHEAD = 12
    real(8), parameter :: W_MIN = 0.8d0, W_STEP = 0.3d0, HEAD_STEP = 30d0
    real(8), parameter :: GRAV = 9.80665d0, PI = 3.141592653589793d0
    real(8), parameter :: DEPTH_MULTI = 30d0
    ! The box: panels per side and half lengths. Well below its first irregular frequency, about 2.2 rad/s.
    integer, parameter :: PX = 4, PY = 3, PZ = 2
    real(8), parameter :: HX = 4d0, HY = 3d0, DRAFT = 3d0
    real(8), parameter :: BODY2_Y = 20d0

    ! Same inputs and the same factorisation, so the entries must agree to the printed digits (7 significant figures).
    real(8), parameter :: EQUAL_TOL = 1d-5
    ! Entries below this share of the largest in the matrix are compared to that floor, as they are zero by the box's symmetry.
    real(8), parameter :: FLOOR_SHARE = 1d-5
    ! GMRES stops at ITER_TOL = 1e-4, but with the same matrix each right side converges the same way.
    real(8), parameter :: GMRES_TOL = 1d-4
    ! Measured on this mesh with the fix, ifort 2022.1 and ifx 2024.2: see the values printed by the tests.
    real(8), parameter :: SYMMETRY_TOL = 5d-3
    real(8), parameter :: HASKIND_TOL = 5d-2

    type :: run_case
        character(len=4) :: name
        integer :: nbody, irsp, solver, body, dof
    end type run_case

    type :: run_result
        real(8), allocatable :: A(:,:,:), B(:,:,:)        ! (frequency, i, j)
        complex(8), allocatable :: X(:,:,:)               ! (frequency, heading, i)
        real(8), allocatable :: rad(:,:), dif(:,:)        ! field-point rows of PressureElevation.6p; rad(:,1) gives the row's body
    end type run_result

    integer, parameter :: NCASE = 15
    type(run_case), parameter :: cases(NCASE) = [ &
        run_case('s0  ', 1, 0, 1, 0, 0), run_case('s1  ', 1, 0, 1, 1, 1), run_case('s3  ', 1, 0, 1, 1, 3), &
        run_case('s5  ', 1, 0, 1, 1, 5), &
        run_case('sl0 ', 1, 1, 1, 0, 0), run_case('sl1 ', 1, 1, 1, 1, 1), run_case('sl2 ', 1, 1, 1, 1, 2), &
        run_case('m0  ', 2, 0, 1, 0, 0), run_case('m1  ', 2, 0, 1, 1, 1), run_case('m2  ', 2, 0, 1, 2, 2), &
        run_case('ml0 ', 2, 1, 1, 0, 0), run_case('ml1 ', 2, 1, 1, 1, 1), run_case('ml3 ', 2, 1, 1, 2, 3), &
        run_case('g0  ', 2, 0, 2, 0, 0), run_case('g1  ', 2, 0, 2, 1, 1) ]

    type(run_result) :: results(NCASE)
    logical :: simulation_execution_failed = .false.

contains

    subroutine setup_gen_modes(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ new_unittest("run_simulations_gen_modes", run_simulations) ]
    end subroutine setup_gen_modes

    subroutine collect_tests_gen_modes(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("single_other_modes_equal_no_mode", test_single), &
            new_unittest("single_lid_other_modes_equal_no_mode", test_single_lid), &
            new_unittest("multi_other_modes_equal_no_mode", test_multi), &
            new_unittest("multi_lid_other_modes_equal_no_mode", test_multi_lid), &
            new_unittest("multi_gmres_other_modes_equal_no_mode", test_gmres), &
            new_unittest("dof1_gives_the_same_mode_as_dof5", test_dof1_equals_dof5), &
            new_unittest("mode_a_and_b_are_symmetric", test_symmetry), &
            new_unittest("mode_satisfies_haskind", test_haskind) &
        ]
    end subroutine collect_tests_gen_modes

    subroutine run_simulations(error)
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: root, inputdir, outputdir, command
        character(len=16) :: errcode
        integer :: i, exitcode
        logical :: success

        root = get_running_exe_path()//"gen_modes/"
        do i = 1, NCASE
            inputdir = root//trim(cases(i)%name)//"_in"
            outputdir = root//trim(cases(i)%name)//"_out"
            call CreateDirectory(inputdir, success)
            if (.not. success) then
                simulation_execution_failed = .true.
                call test_failed(error, "Could not create "//inputdir)
                return
            end if
            call write_inputs(inputdir, cases(i))
            command = native_path(hamsexe_path//" "//inputdir//" "//outputdir//" > "//outputdir//".log")//" 2>&1"
            call execute_command_line(command, exitstat=exitcode, wait=.true.)
            if (exitcode /= 0) then
                simulation_execution_failed = .true.
                write(errcode, '(I0)') exitcode
                call test_failed(error, "HAMS failed on generalized mode case "//trim(cases(i)%name)//" with exit code "//trim(errcode))
                return
            end if
            call read_results(outputdir, cases(i), results(i))
        end do
    end subroutine run_simulations

    ! ---- The checks --------------------------------------------------------------------------------------------------

    subroutine test_single(error)
        type(error_type), allocatable, intent(out) :: error
        call compare_with_no_mode("s0  ", "s1  ", EQUAL_TOL, error)
        if (allocated(error)) return
        call compare_with_no_mode("s0  ", "s3  ", EQUAL_TOL, error)
    end subroutine test_single

    subroutine test_single_lid(error)
        type(error_type), allocatable, intent(out) :: error
        call compare_with_no_mode("sl0 ", "sl1 ", EQUAL_TOL, error)
        if (allocated(error)) return
        call compare_with_no_mode("sl0 ", "sl2 ", EQUAL_TOL, error)
    end subroutine test_single_lid

    subroutine test_multi(error)
        type(error_type), allocatable, intent(out) :: error
        call compare_with_no_mode("m0  ", "m1  ", EQUAL_TOL, error)
        if (allocated(error)) return
        call compare_with_no_mode("m0  ", "m2  ", EQUAL_TOL, error)
    end subroutine test_multi

    subroutine test_multi_lid(error)
        type(error_type), allocatable, intent(out) :: error
        call compare_with_no_mode("ml0 ", "ml1 ", EQUAL_TOL, error)
        if (allocated(error)) return
        call compare_with_no_mode("ml0 ", "ml3 ", EQUAL_TOL, error)
    end subroutine test_multi_lid

    subroutine test_gmres(error)
        type(error_type), allocatable, intent(out) :: error
        call compare_with_no_mode("g0  ", "g1  ", GMRES_TOL, error)
    end subroutine test_gmres

    !> Entries of A, B, X not involving the mode, and field-point values not of the mode, against the run with no mode.
    subroutine compare_with_no_mode(ref_name, name, tol, error)
        character(len=*), intent(in) :: ref_name, name
        real(8), intent(in) :: tol
        type(error_type), allocatable, intent(out) :: error
        integer :: ir, ic, m, n, i, j, k, ih
        real(8) :: worst_a, worst_b, worst_x, worst_fp
        character(len=120) :: info

        if (simulation_execution_failed) then
            call skip_test(error, "Skipping because a simulation failed")
            return
        end if
        ir = find_case(ref_name)
        ic = find_case(name)
        m = 6*(cases(ic)%body - 1) + cases(ic)%dof
        n = 6*cases(ic)%nbody
        worst_a = 0d0
        worst_b = 0d0
        worst_x = 0d0
        do k = 1, NFREQ
            do j = 1, n
                do i = 1, n
                    if (i == m .or. j == m) cycle
                    worst_a = max(worst_a, rel_diff(results(ic)%A(k,i,j), results(ir)%A(k,i,j), maxval(abs(results(ir)%A(k,:,:)))))
                    worst_b = max(worst_b, rel_diff(results(ic)%B(k,i,j), results(ir)%B(k,i,j), maxval(abs(results(ir)%B(k,:,:)))))
                end do
            end do
            do ih = 1, NHEAD
                do i = 1, n
                    if (i == m) cycle
                    worst_x = max(worst_x, rel_diff(abs(results(ic)%X(k,ih,i)), abs(results(ir)%X(k,ih,i)), &
                                                    maxval(abs(results(ir)%X(k,ih,:)))))
                end do
            end do
        end do
        worst_fp = field_point_difference(results(ic), results(ir), cases(ic))

        write(info, '(a,4(a,es9.2))') trim(name), ': largest change against no mode, A', worst_a, ', B', worst_b, &
              ', |X|', worst_x, ', field points', worst_fp
        print '(a)', ' gen_modes '//trim(info)
        call check(error, max(worst_a, worst_b, worst_x, worst_fp) <= tol, "Modes other than the generalized one changed. "//trim(info))
    end subroutine compare_with_no_mode

    !> Each row's values against the run with no mode, leaving out the mode's own pair in the radiation rows.
    function field_point_difference(run, ref, c) result(worst)
        type(run_result), intent(in) :: run, ref
        type(run_case), intent(in) :: c
        real(8) :: worst, scale
        integer :: r, col, first

        worst = 0d0
        scale = maxval(abs(ref%rad(:,2:)))
        first = 2*(c%dof - 1) + 2
        do r = 1, size(ref%rad, 1)
            do col = 2, 13
                if (nint(ref%rad(r,1)) == c%body .and. (col == first .or. col == first + 1)) cycle
                worst = max(worst, rel_diff(run%rad(r,col), ref%rad(r,col), scale))
            end do
        end do
        scale = maxval(abs(ref%dif))
        do r = 1, size(ref%dif, 1)
            do col = 1, 3
                worst = max(worst, rel_diff(run%dif(r,col), ref%dif(r,col), scale))
            end do
        end do
    end function field_point_difference

    !> DOF 5 is not read as the panel's direction, so the same data there must give the same mode as in DOF 1: its own
    !> A, B, X, and its coupling with the modes both runs keep (2, 3, 4, 6). A bad kernel moves all of these.
    subroutine test_dof1_equals_dof5(error)
        type(error_type), allocatable, intent(out) :: error
        integer :: a1, a5, k, ih, j
        integer, parameter :: kept(4) = [2, 3, 4, 6]
        real(8) :: worst
        character(len=80) :: info

        if (simulation_execution_failed) then
            call skip_test(error, "Skipping because a simulation failed")
            return
        end if
        a1 = find_case("s1  ")
        a5 = find_case("s5  ")
        worst = 0d0
        do k = 1, NFREQ
            worst = max(worst, rel_diff(results(a1)%A(k,1,1), results(a5)%A(k,5,5), abs(results(a5)%A(k,5,5))))
            worst = max(worst, rel_diff(results(a1)%B(k,1,1), results(a5)%B(k,5,5), abs(results(a5)%B(k,5,5))))
            do j = 1, size(kept)
                worst = max(worst, rel_diff(results(a1)%A(k,1,kept(j)), results(a5)%A(k,5,kept(j)), abs(results(a5)%A(k,5,5))))
                worst = max(worst, rel_diff(results(a1)%A(k,kept(j),1), results(a5)%A(k,kept(j),5), abs(results(a5)%A(k,5,5))))
                worst = max(worst, rel_diff(results(a1)%B(k,1,kept(j)), results(a5)%B(k,5,kept(j)), abs(results(a5)%B(k,5,5))))
                worst = max(worst, rel_diff(results(a1)%B(k,kept(j),1), results(a5)%B(k,kept(j),5), abs(results(a5)%B(k,5,5))))
            end do
            do ih = 1, NHEAD
                worst = max(worst, rel_diff(abs(results(a1)%X(k,ih,1)), abs(results(a5)%X(k,ih,5)), &
                                            maxval(abs(results(a5)%X(k,:,5)))))
            end do
        end do
        write(info, '(a,es9.2)') 'mode in DOF 1 against DOF 5, largest relative difference', worst
        print '(a)', ' gen_modes '//trim(info)
        call check(error, worst <= EQUAL_TOL, trim(info))
    end subroutine test_dof1_equals_dof5

    !> A(m,j) = A(j,m) and B(m,j) = B(j,m) for the generalized mode m.
    subroutine test_symmetry(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: nchk = 4
        character(len=4), parameter :: names(nchk) = ["s1  ", "sl2 ", "m1  ", "ml3 "]
        integer :: ic, ii, m, j, k
        real(8) :: worst_a, worst_b
        character(len=100) :: info

        if (simulation_execution_failed) then
            call skip_test(error, "Skipping because a simulation failed")
            return
        end if
        do ii = 1, nchk
            ic = find_case(names(ii))
            m = 6*(cases(ic)%body - 1) + cases(ic)%dof
            worst_a = 0d0
            worst_b = 0d0
            do k = 1, NFREQ
                do j = 1, 6*cases(ic)%nbody
                    ! Relative to the mode's self term and the other mode's: a coupling is only a meaningful size against these.
                    worst_a = max(worst_a, abs(results(ic)%A(k,m,j) - results(ic)%A(k,j,m)) &
                                           /sqrt(abs(results(ic)%A(k,m,m)*results(ic)%A(k,j,j))))
                    worst_b = max(worst_b, abs(results(ic)%B(k,m,j) - results(ic)%B(k,j,m)) &
                                           /sqrt(abs(results(ic)%B(k,m,m)*max(abs(results(ic)%B(k,j,j)), 1d-3*abs(results(ic)%B(k,m,m))))))
                end do
            end do
            write(info, '(a,a,es9.2,a,es9.2)') trim(names(ii)), ': mode row against column, A', worst_a, ', B', worst_b
            print '(a)', ' gen_modes '//trim(info)
            call check(error, max(worst_a, worst_b) <= SYMMETRY_TOL, "Not symmetric. "//trim(info))
            if (allocated(error)) return
        end do
    end subroutine test_symmetry

    !> B_mm = k g /(8 pi w Vg) * integral over headings of |X_m|^2, in the output units A/rho, B/(rho w), X = F/(rho g).
    subroutine test_haskind(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: nchk = 4
        character(len=4), parameter :: names(nchk) = ["s1  ", "sl2 ", "m1  ", "ml3 "]
        integer :: ic, ii, m, k
        real(8) :: w, wk, vg, predicted, worst, depth
        character(len=100) :: info

        if (simulation_execution_failed) then
            call skip_test(error, "Skipping because a simulation failed")
            return
        end if
        do ii = 1, nchk
            ic = find_case(names(ii))
            m = 6*(cases(ic)%body - 1) + cases(ic)%dof
            depth = merge(-1d0, DEPTH_MULTI, cases(ic)%nbody == 1)
            worst = 0d0
            do k = 1, NFREQ
                w = W_MIN + (k-1)*W_STEP
                call wave_number(w, depth, wk, vg)
                predicted = wk*GRAV/(8d0*PI*w*vg)*sum(abs(results(ic)%X(k,:,m))**2)*2d0*PI/NHEAD
                worst = max(worst, abs(predicted/results(ic)%B(k,m,m) - 1d0))
            end do
            write(info, '(a,a,es9.2)') trim(names(ii)), ': Haskind against B of the mode, largest relative difference', worst
            print '(a)', ' gen_modes '//trim(info)
            call check(error, worst <= HASKIND_TOL, "Haskind fails. "//trim(info))
            if (allocated(error)) return
        end do
    end subroutine test_haskind

    ! ---- Helpers -------------------------------------------------------------------------------------------------------

    !> The difference of a from the reference b, relative to b, with a floor for entries near zero.
    pure function rel_diff(a, b, scale) result(d)
        real(8), intent(in) :: a, b, scale
        real(8) :: d
        d = abs(a - b)/max(abs(b), FLOOR_SHARE*scale, 1d-30)
    end function rel_diff

    function find_case(name) result(ic)
        character(len=*), intent(in) :: name
        integer :: ic
        do ic = 1, NCASE
            if (trim(cases(ic)%name) == trim(name)) return
        end do
        error stop 'test_integration_gen_modes: unknown case'
    end function find_case

    subroutine wave_number(w, depth, wk, vg)
        real(8), intent(in) :: w, depth
        real(8), intent(out) :: wk, vg
        integer :: it
        wk = w**2/GRAV
        vg = GRAV/(2d0*w)
        if (depth < 0d0) return
        do it = 1, 100
            wk = w**2/(GRAV*tanh(wk*depth))
        end do
        vg = w/(2d0*wk)*(1d0 + 2d0*wk*depth/sinh(2d0*wk*depth))
    end subroutine wave_number

    subroutine read_results(outputdir, c, res)
        character(len=*), intent(in) :: outputdir
        type(run_case), intent(in) :: c
        type(run_result), intent(out) :: res
        integer :: unit, n, ifreq, ih, i, j, k, stat, nrow_r, nrow_d, count, ipos
        real(8) :: w, a, b, beta, amod, phs, re, im, row(14)
        character(len=512) :: line
        real(8), allocatable :: rad(:,:), dif(:,:)

        n = 6*c%nbody
        allocate(res%A(NFREQ,n,n), res%B(NFREQ,n,n), res%X(NFREQ,NHEAD,n))

        open(newunit=unit, file=outputdir//'/Wamit_format/AmssDamp.1', status='old', action='read')
        do ifreq = 1, NFREQ
            do k = 1, n*n
                read(unit, *) w, i, j, a, b
                res%A(ifreq,i,j) = a
                res%B(ifreq,i,j) = b
            end do
        end do
        close(unit)

        open(newunit=unit, file=outputdir//'/Wamit_format/ExcForce.3', status='old', action='read')
        do ifreq = 1, NFREQ
            do ih = 1, NHEAD
                do k = 1, n
                    read(unit, *) w, beta, i, amod, phs, re, im
                    res%X(ifreq,ih,i) = dcmplx(re, im)
                end do
            end do
        end do
        close(unit)

        ! Rows of 14 numbers open a radiation block (w, point, body 1), rows of 12 are the other bodies' values, and rows of 7
        ! are diffraction (w, heading, point, modulus, phase, real, imaginary). The phase is left out, as it is noisy near zero.
        allocate(rad(4000,13), dif(4000,3))
        nrow_r = 0
        nrow_d = 0
        count = 0
        open(newunit=unit, file=outputdir//'/Wamit_format/PressureElevation.6p', status='old', action='read')
        do
            read(unit, '(a)', iostat=stat) line
            if (stat /= 0) exit
            ipos = count_numbers(line)
            read(line, *) row(1:ipos)
            if (ipos == 14) then
                count = 1
                nrow_r = nrow_r + 1
                rad(nrow_r,1) = 1d0
                rad(nrow_r,2:13) = row(3:14)
            else if (ipos == 12) then
                count = count + 1
                nrow_r = nrow_r + 1
                rad(nrow_r,1) = count
                rad(nrow_r,2:13) = row(1:12)
            else if (ipos == 7) then
                nrow_d = nrow_d + 1
                dif(nrow_d,:) = [row(4), row(6), row(7)]
            end if
        end do
        close(unit)
        res%rad = rad(1:nrow_r,:)
        res%dif = dif(1:nrow_d,:)
    end subroutine read_results

    pure function count_numbers(line) result(n)
        character(len=*), intent(in) :: line
        integer :: n, i
        logical :: in_word
        n = 0
        in_word = .false.
        do i = 1, len_trim(line)
            if (line(i:i) == ' ') then
                in_word = .false.
            else if (.not. in_word) then
                in_word = .true.
                n = n + 1
            end if
        end do
    end function count_numbers

    ! ---- Inputs --------------------------------------------------------------------------------------------------------

    subroutine write_inputs(dir, c)
        character(len=*), intent(in) :: dir
        type(run_case), intent(in) :: c
        integer :: b
        character(len=8) :: tag

        call write_control(dir, c)
        do b = 1, c%nbody
            tag = ''
            if (c%nbody > 1) write(tag, '(a,i0)') '_', b
            call write_hull(dir//'/HullMesh'//trim(tag)//'.pnl')
            call write_hydrostatic(dir//'/Hydrostatic'//trim(tag)//'.in', b)
            if (c%irsp == 1) call write_lid(dir//'/WaterPlaneMesh'//trim(tag)//'.pnl')
        end do
        if (c%dof > 0) call write_mode(dir, c%body, c%dof)
    end subroutine write_inputs

    !> The reader uses fixed columns, so each label is written in a field as wide as the columns the reader skips.
    subroutine write_control(dir, c)
        character(len=*), intent(in) :: dir
        type(run_case), intent(in) :: c
        integer :: u, b
        character(len=40) :: label

        open(newunit=u, file=dir//'/ControlFile.in', status='replace', action='write')
        write(u, '(a)') '   --------------HAMS Control file---------------'
        write(u, '(a)') ''
        write(u, '(a14,f30.15)') '   Waterdepth', merge(-1d0, DEPTH_MULTI, c%nbody == 1)
        write(u, '(a)') ''
        write(u, '(a)') '   #Start Definition of Wave Frequencies'
        write(u, '(a27,i16)') '    0_inf_frequency_limits', 0
        write(u, '(a25,i16)') '    Input_frequency_type', 3
        write(u, '(a25,i17)') '    Output_frequency_type', 3
        write(u, '(a26,i16)') '    Number_of_frequencies', -NFREQ
        write(u, '(a27,f30.15)') '    Minimum_frequency_Wmin', W_MIN
        write(u, '(a19,f30.15)') '    Frequency_step', W_STEP
        write(u, '(a)') '   #End Definition of Wave Frequencies'
        write(u, '(a)') ''
        if (c%nbody > 1) then
            write(u, '(a)') '   # Start of definition of bodies for multi-body interaction'
            write(u, '(a24,i16)') '    Number_of_bodies', c%nbody
            do b = 1, c%nbody
                write(label, '(a,i0)') '    LCS_', b
                write(u, '(a26,4f12.3)') trim(label), 0d0, (b-1)*BODY2_Y, 0d0, 0d0
            end do
            write(u, '(a)') '   # End of definition of bodies for multi-body interaction'
            write(u, '(a)') ''
        end if
        write(u, '(a)') '   #Start Definition of Wave Headings'
        write(u, '(a23,i16)') '    Number_of_headings', -NHEAD
        write(u, '(a20,f30.15)') '    Minimum_heading', 0d0
        write(u, '(a17,f30.15)') '    Heading_step', HEAD_STEP
        write(u, '(a)') '   #End Definition of Wave Headings'
        write(u, '(a)') ''
        if (c%nbody == 1) then
            write(u, '(a28,3f12.3)') '    Reference_body_center', 0d0, 0d0, 0d0
        else
            do b = 1, c%nbody
                write(label, '(a,i0)') '    Reference_body_center_', b
                write(u, '(a28,3f12.3)') trim(label), 0d0, (b-1)*BODY2_Y, 0d0
            end do
        end if
        write(u, '(a26,f30.15)') '    Reference_body_length', 1d0
        write(u, '(a26,i16)') '    Wave_diffrac_solution', 1
        write(u, '(a23,i16)') '    If_remove_irr_freq', c%irsp
        write(u, '(a23,i16)') '    Number of threads', 1
        if (c%solver == 2) write(u, '(a,i0)') '    Solver_type              ', c%solver
        write(u, '(a)') ''
        write(u, '(a)') '   #Start Definition of Pressure and/or Elevation (PE)'
        write(u, '(a27,i16)') '    Number_of_field_points', 2
        ! A point in the fluid, and the centroid of a side panel of body 1 for the on-body test.
        write(u, '(3f12.5)') HX + 6d0, 1d0, -1d0
        write(u, '(3f12.5)') HX, HY/PY*(1d0 - PY), -DRAFT/PZ/2d0
        write(u, '(a)') '   #End Definition of Pressure and/or Elevation'
        write(u, '(a)') ''
        write(u, '(a)') '   ----------End HAMS Control file---------------'
        close(u)
    end subroutine write_control

    subroutine write_hydrostatic(path, b)
        character(len=*), intent(in) :: path
        integer, intent(in) :: b
        integer :: u, i, block
        real(8) :: row(6)
        character(len=32), parameter :: titles(5) = [character(len=32) :: ' Body Mass Matrix:', &
            ' External Linear Damping Matrix:', ' External Quadratic Damping:', ' Hydrostatic Restoring Matrix:', &
            ' External Restoring Matrix:']

        open(newunit=u, file=path, status='replace', action='write')
        write(u, '(a)') ' Center of Gravity:'
        write(u, '(3es24.15)') 0d0, (b-1)*BODY2_Y, 0d0
        do block = 1, 5
            write(u, '(a)') titles(block)
            do i = 1, 6
                row = 0d0
                row(i) = 1d0
                write(u, '(6(2x,e12.5))') row
            end do
        end do
        close(u)
    end subroutine write_hydrostatic

    !> Panels in a fixed order, with each one's four corners: sides, then the bottom. Panel numbers match the mode file.
    subroutine box_panels(corners)
        real(8), allocatable, intent(out) :: corners(:,:,:)        ! (panel, corner, xyz)
        integer :: i, j, k, n
        real(8) :: dx, dy, dz, x, y, z

        dx = 2d0*HX/PX
        dy = 2d0*HY/PY
        dz = DRAFT/PZ
        allocate(corners(2*PX*PZ + 2*PY*PZ + PX*PY, 4, 3))
        n = 0
        do k = 1, PZ
            z = -(k-1)*dz
            do i = 1, PX
                x = -HX + (i-1)*dx
                call add(n, [x, -HY, z], [x+dx, -HY, z], [x+dx, -HY, z-dz], [x, -HY, z-dz], [0d0, 1d0, 0d0])
                call add(n, [x, HY, z], [x+dx, HY, z], [x+dx, HY, z-dz], [x, HY, z-dz], [0d0, -1d0, 0d0])
            end do
            do j = 1, PY
                y = -HY + (j-1)*dy
                call add(n, [-HX, y, z], [-HX, y+dy, z], [-HX, y+dy, z-dz], [-HX, y, z-dz], [1d0, 0d0, 0d0])
                call add(n, [HX, y, z], [HX, y+dy, z], [HX, y+dy, z-dz], [HX, y, z-dz], [-1d0, 0d0, 0d0])
            end do
        end do
        do j = 1, PY
            y = -HY + (j-1)*dy
            do i = 1, PX
                x = -HX + (i-1)*dx
                call add(n, [x, y, -DRAFT], [x+dx, y, -DRAFT], [x+dx, y+dy, -DRAFT], [x, y+dy, -DRAFT], [0d0, 0d0, 1d0])
            end do
        end do

    contains

        !> Store a panel so that CalTransNormals gives the wanted normal, which points out of the fluid, into the body.
        subroutine add(np, p1, p2, p3, p4, normal)
            integer, intent(inout) :: np
            real(8), intent(in) :: p1(3), p2(3), p3(3), p4(3), normal(3)
            real(8) :: v13(3), v24(3), un(3)
            v13 = p3 - p1
            v24 = p4 - p2
            un = [v24(2)*v13(3) - v24(3)*v13(2), v24(3)*v13(1) - v24(1)*v13(3), v24(1)*v13(2) - v24(2)*v13(1)]
            np = np + 1
            corners(np,1,:) = p1
            corners(np,3,:) = p3
            if (dot_product(un, normal) > 0d0) then
                corners(np,2,:) = p2
                corners(np,4,:) = p4
            else
                corners(np,2,:) = p4
                corners(np,4,:) = p2
            end if
        end subroutine add

    end subroutine box_panels

    subroutine write_hull(path)
        character(len=*), intent(in) :: path
        real(8), allocatable :: corners(:,:,:)
        call box_panels(corners)
        call write_mesh(path, corners, '    --------------Hull Mesh File---------------')
    end subroutine write_hull

    subroutine write_lid(path)
        character(len=*), intent(in) :: path
        real(8), allocatable :: corners(:,:,:)
        integer :: i, j, n
        real(8) :: dx, dy, x, y

        dx = 2d0*HX/PX
        dy = 2d0*HY/PY
        allocate(corners(PX*PY,4,3))
        n = 0
        do j = 1, PY
            y = -HY + (j-1)*dy
            do i = 1, PX
                x = -HX + (i-1)*dx
                n = n + 1
                corners(n,1,:) = [x, y, 0d0]
                corners(n,2,:) = [x+dx, y, 0d0]
                corners(n,3,:) = [x+dx, y+dy, 0d0]
                corners(n,4,:) = [x, y+dy, 0d0]
            end do
        end do
        call write_mesh(path, corners, '    --------------Waterplane Mesh File---------------')
    end subroutine write_lid

    subroutine write_mesh(path, corners, title)
        character(len=*), intent(in) :: path, title
        real(8), intent(in) :: corners(:,:,:)
        integer :: u, p, c, np
        np = size(corners, 1)
        open(newunit=u, file=path, status='replace', action='write')
        write(u, '(a)') title
        write(u, '(a)') ' '
        write(u, '(a)') '    # Number of Panels, Nodes, X-Symmetry and Y-Symmetry'
        write(u, '(4i12)') np, 4*np, 0, 0
        write(u, '(a)') ' '
        write(u, '(a)') '    # Start Definition of Node Coordinates     ! node_number   x   y   z'
        do p = 1, np
            do c = 1, 4
                write(u, '(i8,3f16.8)') 4*(p-1) + c, corners(p,c,:)
            end do
        end do
        write(u, '(a)') '    # End Definition of Node Coordinates'
        write(u, '(a)') ' '
        write(u, '(a)') '    # Start Definition of Node Relations   ! panel_number  number_of_vertices   Vertex1_ID   Vertex2_ID   Vertex3_ID   (Vertex4_ID)'
        do p = 1, np
            write(u, '(6i8)') p, 4, 4*(p-1) + 1, 4*(p-1) + 2, 4*(p-1) + 3, 4*(p-1) + 4
        end do
        write(u, '(a)') '    # End Definition of Node Relations'
        write(u, '(a)') ' '
        write(u, '(a)') '    --------------End Mesh File---------------'
        close(u)
    end subroutine write_mesh

    !> The mode's values: a smooth function of the panel centroid, with the size and sign change of a real mode shape.
    subroutine write_mode(dir, body, dof)
        character(len=*), intent(in) :: dir
        integer, intent(in) :: body, dof
        real(8), allocatable :: corners(:,:,:)
        real(8) :: xc(3)
        integer :: u, p
        character(len=8) :: tag

        call box_panels(corners)
        write(tag, '(i0,a,i0)') body, '_', dof
        open(newunit=u, file=dir//'/gen_mod_'//trim(tag)//'.txt', status='replace', action='write')
        do p = 1, size(corners, 1)
            xc = sum(corners(p,:,:), dim=1)/4d0
            write(u, '(i8,es24.15)') p, 3d0 + 2d0*xc(1) - 1.5d0*xc(2) + 4d0*xc(3)
        end do
        close(u)
    end subroutine write_mode

end module test_integration_gen_modes
