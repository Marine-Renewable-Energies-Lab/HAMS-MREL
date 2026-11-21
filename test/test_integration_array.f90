!> The multi-body solve against an exact solution: four bottom-mounted cylinders at the corners of a square (drift_le_square,
!> the Linton & Evans array at test size, 64 panels per cylinder), against the multiple-scattering solution of
!> pre_postprocessing_python/SecondOrder/linton_evans_reference.py, whose rows are test/data/drift_le_square/array_reference.txt.
!> The checks are the exciting forces and moments of each body (complex, so phases and the body order count), the near-field
!> mean drift of each body, the far-field mean drift of the array and the total elevation at the field points. At 64 panels per
!> cylinder these are mesh errors, not accuracy claims: the limits are about twice the errors measured when the test was
!> written, and a wrong sign, conjugate, body index or moment centre gives errors of order 1. The converged comparison is in
!> docs/second_order_validation.md, "Linton & Evans array".
module test_integration_array

    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
    use testing_utilities

    implicit none
    private
    public :: setup_array, collect_tests_array

    character(len=*), parameter :: case_name = "drift_le_square"
    integer, parameter :: n_body = 4, n_freq = 2, n_field = 4
    ! Limits, about twice the errors measured on this mesh when the test was written (ifort 2022.1, Windows, 2026-11-20; the
    ! measured value beside each, the largest over both frequencies, 1.60 rad/s in every case but the far-field). At 1.60 rad/s
    ! (k a = 1.32) the 16 panels around a cylinder are 12 per wavelength, which sets these errors. A sign or conjugation error
    ! in X gives 2 |Im X|/|X| or 2, a body or moment-centre error the size of the loads themselves
    real(8), parameter :: exciting_force_limit = 2d-1    ! measured 1.08e-1
    real(8), parameter :: exciting_moment_limit = 3d-1   ! measured 1.42e-1
    real(8), parameter :: near_force_limit = 4.5d-1      ! measured 2.31e-1 (cylinder 2 at 1.60 rad/s)
    real(8), parameter :: near_moment_limit = 1.6d-1     ! measured 7.9e-2 (0.95 rad/s)
    real(8), parameter :: far_force_limit = 4d-2         ! measured 1.8e-2
    real(8), parameter :: elevation_limit = 2d-1         ! measured 9.2e-2
    ! F_z and M_z of a bottom-mounted circle are zero (vertical wall, radial normal), over the largest |F_x|, |F_y| of the
    ! array: measured 1.7e-11
    real(8), parameter :: zero_load_limit = 1d-8

    type :: reference_type
        real(8) :: omega(n_freq), heading
        complex(8) :: x(4, n_body, n_freq)        ! X_1, X_2, X_4, X_5 as NREL + i NIMG
        real(8) :: near(4, n_body, n_freq)        ! F_x, F_y, M_x, M_y
        real(8) :: far(3, n_freq)                 ! F_x, F_y, M_z of the array
        complex(8) :: eta(n_field, n_freq)
    end type reference_type

contains

    subroutine setup_array(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [new_unittest("run_drift_le_square", run_drift_le_square)]
    end subroutine setup_array

    subroutine collect_tests_array(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("array_exciting_exact", array_exciting_exact), &
            new_unittest("array_near_exact", array_near_exact), &
            new_unittest("array_far_exact", array_far_exact), &
            new_unittest("array_elevation_exact", array_elevation_exact) &
        ]
    end subroutine collect_tests_array

    function output_dir() result(dir)
        character(:), allocatable :: dir
        dir = get_running_exe_path()//"output_"//case_name//"/"
    end function output_dir

    subroutine run_drift_le_square(error)
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: inputdir, command
        character(len=16) :: code
        integer :: exitcode

        inputdir = get_repo_root_path()//"test/data/"//case_name//"/input/"
        command = native_path(hamsexe_path//" "//inputdir//" "//output_dir()//" > "//get_running_exe_path()//case_name &
                              //".log")//" 2>&1"
        call execute_command_line(command, exitstat=exitcode, wait=.true.)
        if (exitcode /= 0) then
            write (code, '(I0)') exitcode
            call test_failed(error, "HAMS failed for "//case_name//" with exit code "//trim(code))
        end if
    end subroutine run_drift_le_square

    ! ---------------------------------------------------------------------------------------------- the checks

    !> X_1, X_2 (forces) and X_4, X_5 (moments about the body centre) of each body against the exact values, complex:
    !> max |X - X_exact| over the largest |X| of that kind of that body.
    subroutine array_exciting_exact(error)
        type(error_type), allocatable, intent(out) :: error
        type(reference_type) :: ref
        complex(8) :: x(6*n_body, n_freq)
        real(8) :: worst_f, worst_m, scale
        integer :: f, b, j

        call read_reference(ref, error)
        if (allocated(error)) return
        call read_exciting(ref, x, error)
        if (allocated(error)) return
        worst_f = 0d0
        worst_m = 0d0
        do f = 1, n_freq
            do b = 1, n_body
                scale = max(abs(ref%x(1, b, f)), abs(ref%x(2, b, f)))
                do j = 1, 2
                    worst_f = max(worst_f, abs(x(6*(b-1) + j, f) - ref%x(j, b, f))/scale)
                end do
                scale = max(abs(ref%x(3, b, f)), abs(ref%x(4, b, f)))
                do j = 3, 4
                    worst_m = max(worst_m, abs(x(6*(b-1) + j + 1, f) - ref%x(j, b, f))/scale)
                end do
            end do
        end do
        call report("array_exciting_exact forces", worst_f)
        call report("array_exciting_exact moments", worst_m)
        call check_below(error, "array_exciting_exact forces", worst_f, exciting_force_limit)
        if (allocated(error)) return
        call check_below(error, "array_exciting_exact moments", worst_m, exciting_moment_limit)
    end subroutine array_exciting_exact

    !> Near-field mean F_x, F_y and M_x, M_y of each body (MeanDriftPressure.9) against the exact values, over the largest
    !> |F_x|, |F_y| (or |M_x|, |M_y|) of the four bodies at that frequency, since the drift of one cylinder can be small (0.05
    !> against 2.8 at 1.60 rad/s); F_z and M_z, zero for this hull, over the same force scale.
    subroutine array_near_exact(error)
        type(error_type), allocatable, intent(out) :: error
        type(reference_type) :: ref
        real(8) :: near(6*n_body, n_freq), worst_f, worst_m, worst_z, scale_f, scale_m
        integer :: f, b

        call read_reference(ref, error)
        if (allocated(error)) return
        call read_drift_rows(output_dir()//"Wamit_format/MeanDriftPressure.9", ref, 6*n_body, near, error)
        if (allocated(error)) return
        worst_f = 0d0
        worst_m = 0d0
        worst_z = 0d0
        do f = 1, n_freq
            scale_f = maxval(abs(ref%near(1:2, :, f)))
            scale_m = maxval(abs(ref%near(3:4, :, f)))
            do b = 1, n_body
                worst_f = max(worst_f, abs(near(6*(b-1) + 1, f) - ref%near(1, b, f))/scale_f, &
                              abs(near(6*(b-1) + 2, f) - ref%near(2, b, f))/scale_f)
                worst_z = max(worst_z, abs(near(6*(b-1) + 3, f))/scale_f, abs(near(6*(b-1) + 6, f))/scale_f)
                worst_m = max(worst_m, abs(near(6*(b-1) + 4, f) - ref%near(3, b, f))/scale_m, &
                              abs(near(6*(b-1) + 5, f) - ref%near(4, b, f))/scale_m)
            end do
        end do
        call report("array_near_exact F_x, F_y", worst_f)
        call report("array_near_exact M_x, M_y", worst_m)
        call report("array_near_exact F_z, M_z", worst_z)
        call check_below(error, "array_near_exact F_x, F_y", worst_f, near_force_limit)
        if (allocated(error)) return
        call check_below(error, "array_near_exact M_x, M_y", worst_m, near_moment_limit)
        if (allocated(error)) return
        call check_below(error, "array_near_exact F_z, M_z", worst_z, zero_load_limit)
    end subroutine array_near_exact

    !> Far-field mean F_x, F_y of the array (MeanDriftMomentum.8) against the exact values, over |F|.
    subroutine array_far_exact(error)
        type(error_type), allocatable, intent(out) :: error
        type(reference_type) :: ref
        real(8) :: far(6, n_freq), worst, scale
        integer :: f

        call read_reference(ref, error)
        if (allocated(error)) return
        call read_drift_rows(output_dir()//"Wamit_format/MeanDriftMomentum.8", ref, 6, far, error)
        if (allocated(error)) return
        worst = 0d0
        do f = 1, n_freq
            scale = hypot(ref%far(1, f), ref%far(2, f))
            worst = max(worst, abs(far(1, f) - ref%far(1, f))/scale, abs(far(2, f) - ref%far(2, f))/scale)
        end do
        call report("array_far_exact", worst)
        call check_below(error, "array_far_exact", worst, far_force_limit)
    end subroutine array_far_exact

    !> The total elevation at the field points (diffraction rows of PressureElevation.6p, Wave_diffrac_solution 1) against
    !> the exact one, per unit amplitude.
    subroutine array_elevation_exact(error)
        type(error_type), allocatable, intent(out) :: error
        type(reference_type) :: ref
        character(len=:), allocatable :: lines(:)
        real(8) :: v(5), worst
        integer :: r, f, p, ios, seen

        call read_reference(ref, error)
        if (allocated(error)) return
        call read_lines(output_dir()//"Wamit_format/PressureElevation.6p", lines, error)
        if (allocated(error)) return
        worst = 0d0
        seen = 0
        do r = 1, size(lines)
            if (count_fields(lines(r)) /= 5) cycle   ! radiation rows have 14 fields, their continuation lines 12
            read (lines(r), *, iostat=ios) v
            if (ios /= 0) cycle
            f = frequency_index(ref, v(1))
            p = nint(v(3))
            if (f == 0 .or. p < 1 .or. p > n_field) cycle
            worst = max(worst, abs(cmplx(v(4), v(5), kind=8) - ref%eta(p, f)))
            seen = seen + 1
        end do
        call check(error, seen == n_freq*n_field, "array_elevation_exact: expected 8 diffraction rows")
        if (allocated(error)) return
        call report("array_elevation_exact", worst)
        call check_below(error, "array_elevation_exact", worst, elevation_limit)
    end subroutine array_elevation_exact

    ! ---------------------------------------------------------------------------------------------- reading

    subroutine read_reference(ref, error)
        type(reference_type), intent(out) :: ref
        type(error_type), allocatable, intent(out) :: error
        character(len=:), allocatable :: lines(:)
        character(len=8) :: kind
        real(8) :: omega, heading, v(12)
        integer :: r, b, f, p, nf, ios

        call read_lines(get_repo_root_path()//"test/data/"//case_name//"/array_reference.txt", lines, error)
        if (allocated(error)) return
        nf = 0
        ref%omega = -1d0
        do r = 1, size(lines)
            if (lines(r)(1:1) == '!' .or. len_trim(lines(r)) == 0) cycle
            read (lines(r), *, iostat=ios) kind, omega, heading
            call check(error, ios == 0, "Bad reference row: "//trim(lines(r)))
            if (allocated(error)) return
            f = frequency_index(ref, omega)
            if (f == 0) then
                nf = nf + 1
                call check(error, nf <= n_freq, "More frequencies in the reference than expected")
                if (allocated(error)) return
                ref%omega(nf) = omega
                f = nf
            end if
            ref%heading = heading
            select case (trim(kind))
            case ("body")
                read (lines(r), *) kind, omega, heading, b, v
                ref%x(:, b, f) = [cmplx(v(1), v(2), kind=8), cmplx(v(3), v(4), kind=8), cmplx(v(5), v(6), kind=8), &
                                  cmplx(v(7), v(8), kind=8)]
                ref%near(:, b, f) = v(9:12)
            case ("array")
                read (lines(r), *) kind, omega, heading, v(1:3)
                ref%far(:, f) = v(1:3)
            case ("field")
                read (lines(r), *) kind, omega, heading, p, v(1:2)
                ref%eta(p, f) = cmplx(v(1), v(2), kind=8)
            end select
        end do
        call check(error, nf == n_freq, "Expected 2 frequencies in the reference")
    end subroutine read_reference

    !> NREL + i NIMG of ExcForce.3, by index and frequency, at the reference's heading.
    subroutine read_exciting(ref, x, error)
        type(reference_type), intent(in) :: ref
        complex(8), intent(out) :: x(:, :)
        type(error_type), allocatable, intent(out) :: error
        character(len=:), allocatable :: lines(:)
        real(8) :: v(7)
        integer :: r, f, i, ios, seen

        call read_lines(output_dir()//"Wamit_format/ExcForce.3", lines, error)
        if (allocated(error)) return
        seen = 0
        do r = 1, size(lines)
            read (lines(r), *, iostat=ios) v
            if (ios /= 0) cycle
            f = frequency_index(ref, v(1))
            i = nint(v(3))
            if (f == 0 .or. abs(v(2) - ref%heading) > 1d-6 .or. i < 1 .or. i > size(x, 1)) cycle
            x(i, f) = cmplx(v(6), v(7), kind=8)
            seen = seen + 1
        end do
        call check(error, seen == size(x), "ExcForce.3: expected 6 rows per body and frequency")
    end subroutine read_exciting

    !> Re of a mean drift file (MeanDriftPressure.9, .8), by index and frequency, beta1 = beta2 = the reference heading.
    subroutine read_drift_rows(path, ref, n_index, values, error)
        character(len=*), intent(in) :: path
        type(reference_type), intent(in) :: ref
        integer, intent(in) :: n_index
        real(8), intent(out) :: values(:, :)
        type(error_type), allocatable, intent(out) :: error
        character(len=:), allocatable :: lines(:)
        real(8) :: v(8)
        integer :: r, f, i, ios

        call read_lines(path, lines, error)
        if (allocated(error)) return
        values = 0d0
        do r = 1, size(lines)
            if (lines(r)(1:1) == '#') cycle
            read (lines(r), *, iostat=ios) v
            if (ios /= 0) cycle
            f = frequency_index(ref, v(1))
            i = nint(v(4))
            if (f == 0 .or. abs(v(2) - ref%heading) > 1d-6 .or. abs(v(3) - ref%heading) > 1d-6) cycle
            if (i < 1 .or. i > n_index) cycle
            values(i, f) = v(7)
        end do
    end subroutine read_drift_rows

    integer function frequency_index(ref, omega)
        type(reference_type), intent(in) :: ref
        real(8), intent(in) :: omega
        integer :: f
        frequency_index = 0
        do f = 1, n_freq
            if (abs(ref%omega(f) - omega) < 1d-4) frequency_index = f
        end do
    end function frequency_index

    integer function count_fields(line)
        character(len=*), intent(in) :: line
        integer :: i
        logical :: inside
        count_fields = 0
        inside = .false.
        do i = 1, len_trim(line)
            if (line(i:i) /= ' ') then
                if (.not. inside) count_fields = count_fields + 1
                inside = .true.
            else
                inside = .false.
            end if
        end do
    end function count_fields

    subroutine read_lines(path, lines, error)
        character(len=*), intent(in) :: path
        character(len=:), allocatable, intent(out) :: lines(:)
        type(error_type), allocatable, intent(out) :: error
        character(len=512) :: buffer
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
        allocate (character(len=512) :: lines(n))
        do r = 1, n
            read (unit, '(A)') lines(r)
        end do
        close (unit)
    end subroutine read_lines

    subroutine check_below(error, name, value, limit)
        type(error_type), allocatable, intent(out) :: error
        character(len=*), intent(in) :: name
        real(8), intent(in) :: value, limit
        character(len=160) :: message
        write (message, '(A,": ",ES10.3," above the limit ",ES10.3)') name, value, limit
        call check(error, value <= limit, trim(message))
    end subroutine check_below

    subroutine report(name, value)
        character(len=*), intent(in) :: name
        real(8), intent(in) :: value
        print '(4x,A,": ",ES10.3)', name, value
    end subroutine report

end module test_integration_array
