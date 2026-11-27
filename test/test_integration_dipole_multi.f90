!> Dipole panels on several bodies through the whole program (notes/specs/2026-11-27-dipole-multibody-design.md):
!> two columns (r = 1 m, draft 1 m) with dipole annulus plates (r = 1.5 m) at x = -4 and 4 m, placed by the LCS, deep
!> water, nu = 1 /m, without and with the lids (test/data/dipole_multi/two_bodies/nolid and lid; the body is that of
!> test/data/dipole/column_plate). Checked: the 12 x 12 added mass and damping are reciprocal; body 1's A33 and the
!> coupling A39 equal those of the same solve in the unit suite dipole_multi (multi_two_body_reciprocity); the lid run
!> agrees with the run without lids at this frequency. Also three bodies in finite depth with a plate under body 2
!> (test/data/dipole_multi/three_bodies_finite) run; and dipole panels with symmetry stop.
module test_integration_dipole_multi

    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
    use testing_utilities

    implicit none
    private
    public :: setup_dipole_multi, collect_tests_dipole_multi

contains

    subroutine setup_dipole_multi(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [new_unittest("run_two_body_cases", run_cases)]
    end subroutine setup_dipole_multi

    subroutine collect_tests_dipole_multi(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [new_unittest("two_bodies_reciprocity", test_reciprocity), &
                     new_unittest("two_bodies_against_unit_suite", test_against_unit), &
                     new_unittest("two_bodies_lid", test_lid), &
                     new_unittest("three_bodies_finite_depth_runs", test_three_bodies), &
                     new_unittest("dipole_multi_symmetry_stops", test_symmetry)]
    end subroutine collect_tests_dipole_multi

    function out_dir(name) result(d)
        character(len=*), intent(in) :: name
        character(:), allocatable :: d
        d = get_running_exe_path()//"output_dipole_multi_"//trim(name)//"/"
    end function out_dir

    integer function run_case(dir, name)
        character(len=*), intent(in) :: dir, name
        character(:), allocatable :: command
        command = native_path(hamsexe_path//" "//get_repo_root_path()//dir//"/input/ "//out_dir(name)//" > " &
                              //get_running_exe_path()//"dipole_multi_"//trim(name)//".log")//" 2>&1"
        call execute_command_line(command, exitstat=run_case, wait=.true.)
    end function run_case

    subroutine run_cases(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=*), parameter :: names(3) = [character(len=14) :: "nolid", "lid", "multibody"]
        character(len=48) :: dirs(3)
        character(len=16) :: code
        integer :: k, rc
        dirs = [character(len=48) :: "test/data/dipole_multi/two_bodies/nolid", "test/data/dipole_multi/two_bodies/lid", &
                "test/data/dipole_multi/three_bodies_finite"]
        do k = 1, 3
            rc = run_case(trim(dirs(k)), trim(names(k)))
            if (rc /= 0) then
                write (code, '(I0)') rc
                call test_failed(error, "HAMS failed for "//trim(dirs(k))//" with exit code "//trim(code))
                return
            end if
        end do
    end subroutine run_cases

    !> A(i,j) and B(i,j) (WAMIT form, A/rho and B/(rho omega)) of the first frequency in AmssDamp.1 of case NAME.
    subroutine read_ab(name, n, a, b, ok)
        character(len=*), intent(in) :: name
        integer, intent(in) :: n
        real(8), intent(out) :: a(n,n), b(n,n)
        logical, intent(out) :: ok
        real(8) :: per, x, y, per1
        integer :: u, ios, i, j, cnt
        a = 0d0
        b = 0d0
        cnt = 0
        ok = .false.
        open(newunit=u, file=native_path(out_dir(name)//"Wamit_format/AmssDamp.1"), status="old", action="read", &
             iostat=ios)
        if (ios /= 0) return
        per1 = -1d0
        do
            read(u, *, iostat=ios) per, i, j, x, y
            if (ios /= 0) exit
            if (per1 < 0d0) per1 = per
            if (abs(per - per1) > 1d-9) cycle
            if (i < 1 .or. i > n .or. j < 1 .or. j > n) cycle
            a(i,j) = x
            b(i,j) = y
            cnt = cnt + 1
        end do
        close(u)
        ok = cnt == n*n
    end subroutine read_ab

    !> The 12 x 12 radiation forces f = A + i B (WAMIT form) of the two bodies, without and with the lids, are symmetric
    !> to 2 % of their largest term (the asymmetry of the panel method, as in the unit suite), the coupling block to 5 %
    !> of its largest term.
    subroutine test_reciprocity(error)
        type(error_type), allocatable, intent(out) :: error
        call reciprocity_case("nolid", error)
        if (allocated(error)) return
        call reciprocity_case("lid", error)
    end subroutine test_reciprocity

    subroutine reciprocity_case(name, error)
        character(len=*), intent(in) :: name
        type(error_type), allocatable, intent(out) :: error
        real(8) :: a(12,12), b(12,12), ea, ca
        complex(8) :: f(12,12)
        logical :: ok
        character(len=200) :: msg
        call read_ab(name, 12, a, b, ok)
        call check(error, ok, name//": no complete 12 x 12 AmssDamp.1")
        if (allocated(error)) return
        f = cmplx(a, b, kind=8)
        ea = maxval(abs(f - transpose(f)))/maxval(abs(f))
        ca = maxval(abs(f(1:6,7:12) - transpose(f(7:12,1:6))))/maxval(abs(f(1:6,7:12)))
        write(msg, '(a,a,2es10.2)') name, ': asymmetry of f = A + i B, of its coupling block', ea, ca
        print '(4x,a)', trim(msg)
        call check(error, ea < 2d-2 .and. ca < 5d-2, trim(msg))
    end subroutine reciprocity_case

    !> Body 1's A33 and B33 and the coupling A39, B39 written by the program equal those of the same solve in the unit
    !> suite dipole_multi (multi_two_body_reciprocity prints f(3,3) and f(3,9): A = Re f, B = Im f in WAMIT form), to the
    !> 6 digits printed there: the program wires DipoleMulti, the force sums and the outputs as the unit tests do.
    subroutine test_against_unit(error)
        type(error_type), allocatable, intent(out) :: error
        ! f(3,3) and f(3,9) of multi_two_body_reciprocity (Re, Im)
        real(8), parameter :: ref(4) = [6.595108d0, 1.487337d-1, 1.452779d-2, 4.159393d-2]
        real(8) :: a(12,12), b(12,12), got(4)
        logical :: ok
        character(len=200) :: msg
        call read_ab("nolid", 12, a, b, ok)
        call check(error, ok, "nolid: no complete 12 x 12 AmssDamp.1")
        if (allocated(error)) return
        got = [a(3,3), b(3,3), a(3,9), b(3,9)]
        write(msg, '(a,4es14.6,a,4es14.6)') 'A33, B33, A39, B39', got, '; unit suite', ref
        print '(4x,a)', trim(msg)
        call check(error, all(abs(got - ref) <= 2d-5*abs(ref) + 1d-7), trim(msg))
    end subroutine test_against_unit

    !> With the lids, away from the irregular frequencies of the columns (nu = 1 /m), the 12 x 12 radiation forces
    !> f = A + i B agree with those of the run without lids to 2 % of the largest term: the lids themselves move the
    !> forces by that much in the base code (plain columns: about 1.6 % of the largest term in the coupling block, unit
    !> test multi_two_body_lid; about 3 % of Im f33 for one column, lid_dipole_removes_irregular).
    subroutine test_lid(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: a0(12,12), b0(12,12), a1(12,12), b1(12,12), d
        logical :: ok0, ok1
        character(len=200) :: msg
        call read_ab("nolid", 12, a0, b0, ok0)
        call read_ab("lid", 12, a1, b1, ok1)
        call check(error, ok0 .and. ok1, "two_bodies: AmssDamp.1 missing")
        if (allocated(error)) return
        d = maxval(abs(cmplx(a1 - a0, b1 - b0, kind=8)))/maxval(abs(cmplx(a0, b0, kind=8)))
        write(msg, '(a,es10.2,a,4es14.6)') 'lids against no lids, max |f1 - f0| / max |f0|', d, '; A33, A39 without, with:', &
             a0(3,3), a1(3,3), a0(3,9), a1(3,9)
        print '(4x,a)', trim(msg)
        call check(error, d < 2d-2, trim(msg))
    end subroutine test_lid

    !> test/data/dipole_multi/three_bodies_finite: three cylinders 54 m apart (test/data/config1), 70 m water depth, a
    !> 3 m x 3 m plate of 4 dipole panels 2 m below body 2: the run completes (setup) and writes a complete, finite
    !> 18 x 18 added mass and damping.
    subroutine test_three_bodies(error)
        type(error_type), allocatable, intent(out) :: error
        real(8) :: a(18,18), b(18,18)
        logical :: ok
        call read_ab("multibody", 18, a, b, ok)
        call check(error, ok, "three_bodies_finite: no complete 18 x 18 AmssDamp.1")
        if (allocated(error)) return
        call check(error, all(abs(a) < huge(1d0)) .and. all(abs(b) < huge(1d0)) .and. all(a == a) .and. all(b == b), &
                   "three_bodies_finite: added mass or damping not finite")
    end subroutine test_three_bodies

    !> Dipole panels with the symmetry flag set (test/data/dipole_multi/two_bodies/symmetry): the run stops with the
    !> message of DIPOLE_CHECK_OPTIONS_MULTI.
    subroutine test_symmetry(error)
        type(error_type), allocatable, intent(out) :: error
        character(len=1024) :: line
        integer :: rc, u, ios
        logical :: found
        rc = run_case("test/data/dipole_multi/two_bodies/symmetry", "symmetry")
        call check(error, rc /= 0, "symmetry: the run did not stop")
        if (allocated(error)) return
        found = .false.
        open(newunit=u, file=native_path(get_running_exe_path()//"dipole_multi_symmetry.log"), status="old", &
             action="read", iostat=ios)
        if (ios == 0) then
            do
                read(u, '(a)', iostat=ios) line
                if (ios /= 0) exit
                if (index(line, "dipole panels") > 0 .and. index(line, "symmetry") > 0) found = .true.
            end do
            close(u)
        end if
        call check(error, found, "symmetry: no line with the symmetry message in the log")
    end subroutine test_symmetry

end module test_integration_dipole_multi
