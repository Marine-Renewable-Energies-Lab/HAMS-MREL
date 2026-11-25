!> Dipole panels (DipoleMesh.pnl) change nothing when absent: two single-body cases of test/data/image_term, without and
!> with the lid (IRSP = 0 and 1), give the outputs of the base commit 8f27716 of development_2, stored in
!> test/data/dipole/no_dipole_ref (see its README). Built with the same compiler (ifort on Windows) the files are equal
!> line by line; other compilers (ifx in CI) may differ in the last printed digit, so each line is compared as text and,
!> where the text differs, number by number to 1e-5 of the largest number on the line (outputs print 6 or 7 digits).
!> A trailing carriage return is dropped, so git's line-ending conversion of the stored files does not matter.
!> Also the column with a dipole plate (test/data/dipole/column_plate) through the whole program.
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
        testsuite = [new_unittest("run_no_dipole_cases", run_cases), &
                     new_unittest("run_column_plate_cases", run_plate_cases)]
    end subroutine setup_dipole

    subroutine collect_tests_dipole(testsuite)
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
        testsuite = [new_unittest("no_dipole_bit_identical", test_identical), &
                     new_unittest("column_plate_end_to_end", test_plate)]
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

    !> The column with a dipole annulus plate (test/data/dipole/column_plate), without and with the lid, deep water,
    !> omega = sqrt(g) (nu = 1 /m), through the whole program.
    subroutine run_plate_cases(error)
        type(error_type), allocatable, intent(out) :: error
        character(:), allocatable :: command
        character(len=16) :: code
        integer :: exitcode, k
        character(len=5), parameter :: pc(2) = ["nolid", "lid  "]
        do k = 1, 2
            command = native_path(hamsexe_path//" "//get_repo_root_path()//"test/data/dipole/column_plate/"//trim(pc(k)) &
                                  //"/input/ "//get_running_exe_path()//"output_column_plate_"//trim(pc(k))//"/ > " &
                                  //get_running_exe_path()//"column_plate_"//trim(pc(k))//".log")//" 2>&1"
            call execute_command_line(command, exitstat=exitcode, wait=.true.)
            if (exitcode /= 0) then
                write (code, '(I0)') exitcode
                call test_failed(error, "HAMS failed for column_plate/"//trim(pc(k))//" with exit code "//trim(code))
                return
            end if
        end do
    end subroutine run_plate_cases

    !> A33 and B33 (WAMIT form, L = 1 m) written by the program equal those of the same solve in the unit suite
    !> dipole_bodies (lid_dipole_removes_irregular, nu = 1: 6.59035 and 0.138370 without the lid, 6.58782 and 0.145919
    !> with it), to the 6 digits printed there: the program wires DipoleAssb, the force sums and the outputs as the unit
    !> tests do.
    subroutine test_plate(error)
        type(error_type), allocatable, intent(out) :: error
        real(8), parameter :: ref(2,2) = reshape([6.59035d0, 0.138370d0, 6.58782d0, 0.145919d0], [2,2])
        character(len=5), parameter :: pc(2) = ["nolid", "lid  "]
        real(8) :: per, a, b
        integer :: u, ios, i, j, k
        character(len=200) :: msg
        logical :: found
        do k = 1, 2
            open(newunit=u, file=native_path(get_running_exe_path()//"output_column_plate_"//trim(pc(k)) &
                 //"/Wamit_format/AmssDamp.1"), status="old", action="read", iostat=ios)
            if (ios /= 0) then
                call test_failed(error, "no AmssDamp.1 for column_plate/"//trim(pc(k)))
                return
            end if
            found = .false.
            do
                read(u, *, iostat=ios) per, i, j, a, b
                if (ios /= 0) exit
                if (i == 3 .and. j == 3) then
                    found = .true.
                    exit
                end if
            end do
            close(u)
            write(msg, '(a,a,a,2es14.6,a,2es14.6)') 'column_plate/', trim(pc(k)), ': A33, B33', a, b, '; unit suite', &
                 ref(:,k)
            print '(4x,a)', trim(msg)
            call check(error, found .and. abs(a - ref(1,k)) <= 2d-5*ref(1,k) .and. abs(b - ref(2,k)) <= 2d-5*ref(2,k), &
                       trim(msg))
            if (allocated(error)) return
        end do
    end subroutine test_plate

    !> Every line of the five WAMIT-format outputs equal to the stored base output (or to its printed precision).
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
                    if (trim(strip_cr(la)) /= trim(strip_cr(lb))) then
                        if (.not. lines_close(strip_cr(la), strip_cr(lb))) exit
                    end if
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
        !> The same tokens, numbers within 1e-5 of the line's largest number, words equal.
        logical function lines_close(a, b)
            character(len=*), intent(in) :: a, b
            character(len=64) :: ta(64), tb(64)
            real(8) :: va(64), vb(64), big
            logical :: na(64), nb(64)
            integer :: ka, kb, i
            call split(a, ta, va, na, ka)
            call split(b, tb, vb, nb, kb)
            lines_close = (ka == kb)
            if (.not. lines_close) return
            big = 0d0
            do i = 1, ka
                if (na(i)) big = max(big, abs(va(i)))
                if (nb(i)) big = max(big, abs(vb(i)))
            end do
            do i = 1, ka
                if (na(i) .neqv. nb(i)) then
                    lines_close = .false.
                else if (na(i)) then
                    if (abs(va(i) - vb(i)) > 1d-5*big) lines_close = .false.
                else if (ta(i) /= tb(i)) then
                    lines_close = .false.
                end if
                if (.not. lines_close) return
            end do
        end function lines_close

        subroutine split(line, tok, val, isnum, n)
            character(len=*), intent(in) :: line
            character(len=64), intent(out) :: tok(64)
            real(8), intent(out) :: val(64)
            logical, intent(out) :: isnum(64)
            integer, intent(out) :: n
            integer :: i, j, ios
            n = 0
            i = 1
            do while (i <= len_trim(line) .and. n < 64)
                if (line(i:i) == ' ') then
                    i = i + 1
                    cycle
                end if
                j = i
                do while (j <= len_trim(line))
                    if (line(j:j) == ' ') exit
                    j = j + 1
                end do
                n = n + 1
                tok(n) = line(i:j-1)
                read(tok(n), *, iostat=ios) val(n)
                isnum(n) = (ios == 0)
                i = j
            end do
        end subroutine split

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
