!> Dipole panels (DipoleMesh.pnl) change nothing when absent: two single-body cases of test/data/image_term, without and
!> with the lid (IRSP = 0 and 1), give outputs identical line by line, to the last printed digit, to those of the base
!> commit 8f27716 of development_2, stored in test/data/dipole/no_dipole_ref (see its README). Lines are compared as text,
!> with a trailing carriage return dropped, so that git's line-ending conversion of the stored files does not matter.
module test_integration_dipole

    use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
    use testing_utilities

    implicit none
    private
    public :: setup_dipole, collect_tests_dipole

    character(len=*), parameter :: cases(2) = [character(len=13) :: "image_cos", "image_cos_lid"]
    character(len=*), parameter :: files(5) = [character(len=20) :: "AmssDamp.1", "ExcForce.3", "Hydrostat.hst", &
                                                                     "Motion.4", "PressureElevation.6p"]

contains

    subroutine setup_dipole(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [new_unittest("run_no_dipole_cases", run_cases)]
    end subroutine setup_dipole

    subroutine collect_tests_dipole(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [new_unittest("no_dipole_bit_identical", test_identical)]
    end subroutine collect_tests_dipole

    function out_dir(name) result(d)
        character(len=*), intent(in) :: name
        character(:), allocatable :: d
        d = get_running_exe_path()//"output_nodipole_"//trim(name)//"/"
    end function out_dir

    subroutine run_cases(error)
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: command
        character(len=16) :: code
        integer :: exitcode, k
        do k = 1, size(cases)
            command = native_path(hamsexe_path//" "//get_repo_root_path()//"test/data/image_term/"//trim(cases(k)) &
                                  //"/input/ "//out_dir(cases(k))//" > "//get_running_exe_path()//"nodipole_" &
                                  //trim(cases(k))//".log")//" 2>&1"
            call execute_command_line(command, exitstat=exitcode, wait=.true.)
            if (exitcode /= 0) then
                write (code, '(I0)') exitcode
                call test_failed(error, "HAMS failed for "//trim(cases(k))//" with exit code "//trim(code))
                return
            end if
        end do
    end subroutine run_cases

    !> Every line of the five WAMIT-format outputs equal to the stored base output.
    subroutine test_identical(error)
        type(error_type), allocatable, intent(out) :: error
        integer :: k, f, ua, ub, ia, ib, line
        character(len=4096) :: la, lb
        character(len=:), allocatable :: pa, pb
        do k = 1, size(cases)
            do f = 1, size(files)
                pa = get_repo_root_path()//"test/data/dipole/no_dipole_ref/"//trim(cases(k))//"/"//trim(files(f))
                pb = out_dir(cases(k))//"Wamit_format/"//trim(files(f))
                open(newunit=ua, file=native_path(pa), status="old", action="read", iostat=ia)
                open(newunit=ub, file=native_path(pb), status="old", action="read", iostat=ib)
                if (ia /= 0 .or. ib /= 0) then
                    call test_failed(error, "cannot open "//pa//" or "//pb)
                    return
                end if
                line = 0
                do
                    read(ua, '(a)', iostat=ia) la
                    read(ub, '(a)', iostat=ib) lb
                    line = line + 1
                    if (ia /= 0 .or. ib /= 0) exit
                    if (trim(strip_cr(la)) /= trim(strip_cr(lb))) exit
                end do
                close(ua)
                close(ub)
                call check(error, ia /= 0 .and. ib /= 0 .and. ia == ib, trim(cases(k))//"/"//trim(files(f)) &
                           //": differs from the base output (or has a different length)")
                if (allocated(error)) then
                    write(*, '(a,i0)') "    first difference at line ", line
                    return
                end if
            end do
        end do
    contains
        function strip_cr(s) result(t)
            character(len=*), intent(in) :: s
            character(len=len(s)) :: t
            integer :: n
            t = s
            n = len_trim(t)
            if (n > 0) then
                if (t(n:n) == achar(13)) t(n:n) = " "
            end if
        end function strip_cr
    end subroutine test_identical

end module test_integration_dipole
