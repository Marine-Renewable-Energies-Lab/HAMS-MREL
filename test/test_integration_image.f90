!> Whole-solver checks of the panel integral of the free-surface image term, on one cylinder through every path. The cases are
!> described in test/data/image_term/README.md. The limits come from runs before and after the change, listed there.
module test_integration_image

    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
    use testing_utilities

    implicit none
    private
    public :: setup_image, collect_tests_image

    real(8), parameter :: g = 9.80665d0, omega = 1.35d0, depth = 50d0

    real(8), parameter :: field_limit = 1d-5          ! pressures and forces are printed to 7 digits
    real(8), parameter :: half_limit = 1d-4
    real(8), parameter :: gmres_limit = 1d-3          ! as test_integration_gmres
    real(8), parameter :: two_body_limit = 1d-3       ! the far body still radiates into this one: B33 differs by 5e-4
    real(8), parameter :: lid_limit = 5d-3
    real(8), parameter :: top_row_limit = 0.2d0
    real(8), parameter :: haskind_limit = 1d-2

    type :: coefficients
        real(8) :: a(6), b(6)
        complex(8) :: x(6)
    end type coefficients

contains

    subroutine setup_image(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("run_image_cos", run_image_cos), &
            new_unittest("run_image_cos_half", run_image_cos_half), &
            new_unittest("run_image_cos_gmres", run_image_cos_gmres), &
            new_unittest("run_image_cos_lid", run_image_cos_lid), &
            new_unittest("run_image_cos_two", run_image_cos_two), &
            new_unittest("run_image_uni", run_image_uni), &
            new_unittest("run_image_cos_fine", run_image_cos_fine), &
            new_unittest("run_image_uni_fine", run_image_uni_fine) &
            ]
    end subroutine setup_image

    subroutine collect_tests_image(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [ &
            new_unittest("field_points_match_collocation", field_points_match_collocation), &
            new_unittest("field_points_match_collocation_multi", field_points_match_collocation_multi), &
            new_unittest("top_rows_smooth", top_rows_smooth), &
            new_unittest("half_matches_full", half_matches_full), &
            new_unittest("gmres_matches_direct", gmres_matches_direct), &
            new_unittest("multi_matches_single", multi_matches_single), &
            new_unittest("lid_matches_no_lid", lid_matches_no_lid), &
            new_unittest("cosine_and_uniform_approach", cosine_and_uniform_approach), &
            new_unittest("haskind_heave", haskind_heave) &
            ]
    end subroutine collect_tests_image

    ! -------------------------------------------------------------------------------------------- runs

    subroutine run_case(name, error)
        character(len=*), intent(in) :: name
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: command
        character(len=16) :: code
        integer :: exitcode

        command = native_path(hamsexe_path//" "//get_repo_root_path()//"test/data/image_term/"//name//"/input/ " &
                              //output_dir(name)//" > "//get_running_exe_path()//name//".log")//" 2>&1"
        call execute_command_line(command, exitstat=exitcode, wait=.true.)
        if (exitcode /= 0) then
            write (code, '(I0)') exitcode
            call test_failed(error, "HAMS failed for "//name//" with exit code "//trim(code))
        end if
    end subroutine run_case

    function output_dir(name) result(dir)
        character(len=*), intent(in) :: name
        character(:), allocatable :: dir
        dir = get_running_exe_path()//"output_"//name//"/"
    end function output_dir

    subroutine run_image_cos(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("image_cos", error)
    end subroutine run_image_cos

    subroutine run_image_cos_half(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("image_cos_half", error)
    end subroutine run_image_cos_half

    subroutine run_image_cos_gmres(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("image_cos_gmres", error)
    end subroutine run_image_cos_gmres

    subroutine run_image_cos_lid(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("image_cos_lid", error)
    end subroutine run_image_cos_lid

    subroutine run_image_cos_two(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("image_cos_two", error)
    end subroutine run_image_cos_two

    subroutine run_image_uni(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("image_uni", error)
    end subroutine run_image_uni

    subroutine run_image_cos_fine(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("image_cos_fine", error)
    end subroutine run_image_cos_fine

    subroutine run_image_uni_fine(error)
        type(error_type), allocatable, intent(out) :: error
        call run_case("image_uni_fine", error)
    end subroutine run_image_uni_fine

    ! -------------------------------------------------------------------------------------------- readers

    !> Diagonal A and B, and X, of the first body. The first six rows of a body pair are read by their modes.
    subroutine read_coefficients(name, c, error)
        character(len=*), intent(in) :: name
        type(coefficients), intent(out) :: c
        type(error_type), allocatable, intent(out) :: error
        real(8) :: freq, beta, re, im, v1, v2
        integer :: unit, i, j, stat

        c%a = 0d0
        c%b = 0d0
        c%x = (0d0, 0d0)
        open (newunit=unit, file=output_dir(name)//"Wamit_format/AmssDamp.1", status="old", action="read", iostat=stat)
        if (stat /= 0) then
            call test_failed(error, "No AmssDamp.1 for "//name)
            return
        end if
        do
            read (unit, *, iostat=stat) freq, i, j, v1, v2
            if (stat /= 0) exit
            if (i == j .and. i <= 6) then
                c%a(i) = v1
                c%b(i) = v2
            end if
        end do
        close (unit)
        open (newunit=unit, file=output_dir(name)//"Wamit_format/ExcForce.3", status="old", action="read", iostat=stat)
        if (stat /= 0) then
            call test_failed(error, "No ExcForce.3 for "//name)
            return
        end if
        do
            read (unit, *, iostat=stat) freq, beta, i, v1, v2, re, im
            if (stat /= 0) exit
            if (i <= 6) c%x(i) = cmplx(re, im, kind=8)
        end do
        close (unit)
    end subroutine read_coefficients

    !> Panel centroids, areas and unit normals out of the body, from a mesh file.
    subroutine read_panels(name, file, xc, area, normal, error)
        character(len=*), intent(in) :: name, file
        real(8), allocatable, intent(out) :: xc(:,:), area(:), normal(:,:)
        type(error_type), allocatable, intent(out) :: error
        character(len=512) :: line
        real(8), allocatable :: node(:,:)
        real(8) :: p(3,4), v(3), scale_v
        integer :: unit, stat, np, nn, isx, isy, i, k, id, nv, ids(4)

        open (newunit=unit, file=get_repo_root_path()//"test/data/image_term/"//name//"/input/"//file, status="old", &
              action="read", iostat=stat)
        if (stat /= 0) then
            call test_failed(error, "No mesh for "//name)
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
            ! The body is a cylinder on the z axis of its own frame: out of the body is away from the axis, and down on the bottom
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

    !> Rows of a .6p file with ncol columns: the point in column icol_point, the value in the last two columns.
    subroutine read_6p(name, ncol, icol_point, values, error)
        character(len=*), intent(in) :: name
        integer, intent(in) :: ncol, icol_point
        complex(8), allocatable, intent(out) :: values(:)
        type(error_type), allocatable, intent(out) :: error
        character(len=512) :: line
        real(8) :: row(14)
        integer :: unit, stat, n, count_cols, ipoint

        allocate (values(0))
        open (newunit=unit, file=output_dir(name)//"Wamit_format/PressureElevation.6p", status="old", action="read", iostat=stat)
        if (stat /= 0) then
            call test_failed(error, "No PressureElevation.6p for "//name)
            return
        end if
        do
            read (unit, '(a)', iostat=stat) line
            if (stat /= 0) exit
            count_cols = columns(line)
            if (count_cols /= ncol) cycle
            read (line, *) row(1:ncol)
            ipoint = nint(row(icol_point))
            n = size(values)
            if (ipoint > n) values = [values, (cmplx(0d0, 0d0, kind=8), stat = n + 1, ipoint)]
            values(ipoint) = cmplx(row(ncol - 1), row(ncol), kind=8)
        end do
        close (unit)
    end subroutine read_6p

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

    !> Relative difference of two sets of diagonal A and B, over the modes that are not zero by symmetry.
    real(8) function worst_ab(c1, c2)
        type(coefficients), intent(in) :: c1, c2
        integer, parameter :: modes(5) = [1, 2, 3, 4, 5]
        integer :: m

        worst_ab = 0d0
        do m = 1, size(modes)
            worst_ab = max(worst_ab, abs(c1%a(modes(m)) - c2%a(modes(m)))/abs(c2%a(modes(m))), &
                           abs(c1%b(modes(m)) - c2%b(modes(m)))/abs(c2%b(modes(m))))
        end do
    end function worst_ab

    subroutine check_below(error, name, value, limit)
        type(error_type), allocatable, intent(out) :: error
        character(len=*), intent(in) :: name
        real(8), intent(in) :: value, limit
        character(len=64) :: text

        write (*, '(4x, a, ": ", es10.3)') name, value
        write (text, '(es10.3, " is not below ", es10.3)') value, limit
        call check(error, value < limit, name//": "//trim(text))
    end subroutine check_below

    ! -------------------------------------------------------------------------------------------- tests

    !> The potential a field point gets from the boundary integral must be the one the solver found on the panel: the pressure
    !> at every centroid, integrated over the hull, is the excitation force. Fails if collocation and field evaluation treat the
    !> image term differently.
    subroutine field_points_match_collocation(error)
        type(error_type), allocatable, intent(out) :: error

        call check_field_force("image_cos", "HullMesh.pnl", 7, 3, 1, field_limit, error)
    end subroutine field_points_match_collocation

    !> The same through the multi-body path, on the first of two bodies. Its rows have one column fewer.
    subroutine field_points_match_collocation_multi(error)
        type(error_type), allocatable, intent(out) :: error

        call check_field_force("image_cos_two", "HullMesh_1.pnl", 5, 3, 1, field_limit, error)
    end subroutine field_points_match_collocation_multi

    subroutine check_field_force(name, mesh_file, ncol, icol_point, mode, limit, error)
        character(len=*), intent(in) :: name, mesh_file
        integer, intent(in) :: ncol, icol_point, mode
        real(8), intent(in) :: limit
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: xc(:,:), area(:), normal(:,:)
        complex(8), allocatable :: pressure(:)
        type(coefficients) :: c
        complex(8) :: force(3)
        integer :: i, k
        real(8) :: worst

        call read_panels(name, mesh_file, xc, area, normal, error)
        if (allocated(error)) return
        call read_6p(name, ncol, icol_point, pressure, error)
        if (allocated(error)) return
        call check(error, size(pressure) == size(area), "Not every panel of "//name//" has a pressure")
        if (allocated(error)) return
        call read_coefficients(name, c, error)
        if (allocated(error)) return
        force = (0d0, 0d0)
        do i = 1, size(area)
            do k = 1, 3
                force(k) = force(k) + pressure(i)*area(i)*normal(k, i)
            end do
        end do
        worst = 0d0
        do k = mode, 3, 2
            ! The sign of the force follows the normal of the solver, which points into the body
            worst = max(worst, abs(abs(force(k) + c%x(k)) - 0d0)/abs(c%x(k)))
        end do
        call check_below(error, name//" force from field points", worst, limit)
    end subroutine check_field_force

    !> On rows thin towards the waterline, the slope of the radiated potential down the wall between the top two rows must follow
    !> the slope of the rows below. The one-point image gave a kink: 48% against 11% with the integral, at 720 panels.
    subroutine top_rows_smooth(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), allocatable :: xc(:,:), area(:), normal(:,:)
        complex(8), allocatable :: p(:)
        complex(8) :: slope_top, slope_next
        character(len=512) :: line
        real(8) :: row(14)
        integer :: unit, stat, ip

        call read_panels("image_cos_fine", "HullMesh.pnl", xc, area, normal, error)
        if (allocated(error)) return
        allocate (p(3))
        open (newunit=unit, file=output_dir("image_cos_fine")//"Wamit_format/PressureElevation.6p", status="old", action="read", &
              iostat=stat)
        if (stat /= 0) then
            call test_failed(error, "No PressureElevation.6p")
            return
        end if
        do
            read (unit, '(a)', iostat=stat) line
            if (stat /= 0) exit
            if (columns(line) /= 14) cycle
            read (line, *) row
            ip = nint(row(2))
            if (ip <= 3) p(ip) = cmplx(row(3), row(4), kind=8)     ! surge
        end do
        close (unit)
        slope_top = (p(1) - p(2))/(xc(3, 1) - xc(3, 2))
        slope_next = (p(2) - p(3))/(xc(3, 2) - xc(3, 3))
        call check_below(error, "top_rows_smooth", abs(slope_top - slope_next)/abs(slope_next), top_row_limit)
    end subroutine top_rows_smooth

    subroutine half_matches_full(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: full, half

        call read_coefficients("image_cos", full, error)
        if (allocated(error)) return
        call read_coefficients("image_cos_half", half, error)
        if (allocated(error)) return
        call check_below(error, "half_matches_full", worst_ab(half, full), half_limit)
    end subroutine half_matches_full

    subroutine gmres_matches_direct(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: direct, gmres

        call read_coefficients("image_cos", direct, error)
        if (allocated(error)) return
        call read_coefficients("image_cos_gmres", gmres, error)
        if (allocated(error)) return
        call check_below(error, "gmres_matches_direct", worst_ab(gmres, direct), gmres_limit)
    end subroutine gmres_matches_direct

    !> One body through the single-body path, and with a second body 5 km away through the multi-body path.
    subroutine multi_matches_single(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: single, multi

        call read_coefficients("image_cos", single, error)
        if (allocated(error)) return
        call read_coefficients("image_cos_two", multi, error)
        if (allocated(error)) return
        call check_below(error, "multi_matches_single", worst_ab(multi, single), two_body_limit)
    end subroutine multi_matches_single

    !> Below the first irregular frequency, a lid must not change much. The one-point image on the hull rows next to it gave 0.88%.
    subroutine lid_matches_no_lid(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: bare, lid

        call read_coefficients("image_cos", bare, error)
        if (allocated(error)) return
        call read_coefficients("image_cos_lid", lid, error)
        if (allocated(error)) return
        call check_below(error, "lid_matches_no_lid", worst_ab(lid, bare), lid_limit)
    end subroutine lid_matches_no_lid

    !> Cosine and uniform rows get closer as the mesh is refined, in every coefficient.
    subroutine cosine_and_uniform_approach(error)
        type(error_type), allocatable, intent(out) :: error
        type(coefficients) :: cos1, uni1, cos2, uni2
        integer :: m
        real(8) :: coarse, fine

        call read_coefficients("image_cos", cos1, error)
        if (allocated(error)) return
        call read_coefficients("image_uni", uni1, error)
        if (allocated(error)) return
        call read_coefficients("image_cos_fine", cos2, error)
        if (allocated(error)) return
        call read_coefficients("image_uni_fine", uni2, error)
        if (allocated(error)) return
        do m = 1, 5, 2
            coarse = abs(cos1%a(m) - uni1%a(m))/abs(uni1%a(m))
            fine = abs(cos2%a(m) - uni2%a(m))/abs(uni2%a(m))
            write (*, '(4x, "A mode ", i0, " cosine against uniform: ", 2es10.3)') m, coarse, fine
            call check(error, fine < coarse, "A: cosine and uniform rows do not get closer")
            if (allocated(error)) return
            coarse = abs(cos1%b(m) - uni1%b(m))/abs(uni1%b(m))
            fine = abs(cos2%b(m) - uni2%b(m))/abs(uni2%b(m))
            write (*, '(4x, "B mode ", i0, " cosine against uniform: ", 2es10.3)') m, coarse, fine
            call check(error, fine < coarse, "B: cosine and uniform rows do not get closer")
            if (allocated(error)) return
            coarse = abs(abs(cos1%x(m)) - abs(uni1%x(m)))/abs(uni1%x(m))
            fine = abs(abs(cos2%x(m)) - abs(uni2%x(m)))/abs(uni2%x(m))
            write (*, '(4x, "X mode ", i0, " cosine against uniform: ", 2es10.3)') m, coarse, fine
            call check(error, fine < coarse, "X: cosine and uniform rows do not get closer")
            if (allocated(error)) return
        end do
    end subroutine cosine_and_uniform_approach

    !> The Haskind relation for heave of an axisymmetric body: B33 = k |X3|^2 g / (4 omega Vg), non-dimensional with REFL = 1.
    subroutine haskind_heave(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=16), parameter :: names(3) = [character(len=16) :: "image_cos", "image_cos_fine", "image_uni_fine"]
        type(coefficients) :: c
        real(8) :: k, vg, worst, predicted
        integer :: i

        k = omega*omega/g
        do i = 1, 100
            k = omega*omega/(g*tanh(k*depth))
        end do
        vg = omega/(2d0*k)*(1d0 + 2d0*k*depth/sinh(2d0*k*depth))
        worst = 0d0
        do i = 1, size(names)
            call read_coefficients(trim(names(i)), c, error)
            if (allocated(error)) return
            predicted = k*g*abs(c%x(3))**2/(4d0*omega*vg)
            worst = max(worst, abs(predicted - c%b(3))/c%b(3))
        end do
        call check_below(error, "haskind_heave", worst, haskind_limit)
    end subroutine haskind_heave

end module test_integration_image
