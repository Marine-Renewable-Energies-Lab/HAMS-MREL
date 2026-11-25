module testing_utilities

    use testdrive, only : error_type, check

    implicit none

    public :: get_running_exe_path
    public :: get_repo_root_path
    public :: native_path
    public :: test_hams_matrices_are_close

    ! Path to HAMS suppied by the user
    public :: hamsexe_path
    character(len=:), allocatable :: hamsexe_path
    public :: get_hams_path_from_command_line_args

    contains

    subroutine get_hams_path_from_command_line_args()
        integer :: n
        logical :: exists
        character(len=1024) :: arg

        n = command_argument_count()
        if (n >= 1) then 
            if (n > 1) then
                print *, "Warning: more than one argument supplied to test app."
            end if
            ! Get HAMS exe
            call get_command_argument(1, arg)
            hamsexe_path = trim(arg)
            ! Check file exists in path provided by ser
            inquire(file=hamsexe_path, exist=exists)
            if (.not. exists) then
                ! Check file exists in repo-root//path-by-user
                inquire(file=get_repo_root_path()//hamsexe_path, exist=exists)
                if (.not. exists) then
                    print *, "Error: HAMS executable needed for integration tests not found in: ", hamsexe_path
                    print *, "Expected usage: fpm test -- <path to hams relative to root of repo>"
                    stop 1
                end if
            else
                print*, "Integration tests: using HAMS executable in ", hamsexe_path
            end if
        else
            print*, "Error: path to HAMS executable needed for integration tests not found. Expected usage: fpm test -- <path to hams relative to root of repo>"
            stop 1
        end if
    end subroutine get_hams_path_from_command_line_args


    ! Commands run through cmd.exe on Windows, which reads '/' as the start of a switch:
    ! "install/bin/hamsmrel.exe" fails with "'install' is not recognized". Use '\' there.
    function native_path(path) result(native)
        implicit none
        character(len=*), intent(in) :: path
        character(:), allocatable :: native
        character(len=32) :: os
        integer :: i, stat

        native = path
        call get_environment_variable("OS", os, status=stat)
        if (stat == 0 .and. trim(os) == "Windows_NT") then
            do i = 1, len(native)
                if (native(i:i) == '/') native(i:i) = '\'
            end do
        end if
    end function native_path


    function get_running_exe_path() result(exepath)
        implicit none
        character(:), allocatable :: exepath
        character(len=1024) :: fullpath
        integer :: i

        ! Get path of running executable
        call get_command_argument(0, fullpath)

        ! Strip executable name (Find the position of the last '/' or '\')
        i = len_trim(fullpath)
        do while (i > 0)
            if (fullpath(i:i) == '/' .or. fullpath(i:i) == '\') exit
            i = i - 1
        end do
        exepath = trim(fullpath(1:i-1))//'/'
    end function get_running_exe_path


    function get_repo_root_path() result(reporootpath)
        implicit none
        character(:), allocatable :: reporootpath

        ! Build full path to top level of the repository,
        ! assuming exe is in repo/build/comp_id/test/
        ! which is fpm's default location for tests
        reporootpath = get_running_exe_path()//'../../../'
    end function get_repo_root_path


    ! Reads numeric data from a text file into a matrix using a user-specified Fortran format string.
    ! Arguments:
    !   fname (in)  - character string, path to input file
    !   fmt   (in)  - character string, Fortran format string for reading
    !   ncol  (in)  - integer, number of columns in the file per row
    !   A     (out) - allocatable real(kind=8) matrix
    subroutine read_matrix_with_format(filename, formatstr, ncol, A)
        implicit none
        character(len=*), intent(in) :: filename
        character(len=*), intent(in) :: formatstr
        integer, intent(in) :: ncol
        real(kind=8), allocatable, intent(out) :: A(:,:)

        integer, parameter :: maxrows = 10000
        real(kind=8), allocatable :: temp(:,:)
        real(kind=8) :: row(ncol)
        integer :: ios, nrows, unit

        ! Allocate temporary array
        allocate(temp(maxrows, ncol))
        nrows = 0

        ! Open file
        open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, "Error opening file:", trim(filename), "iostat=", ios
            stop 1
        end if

        ! Read file line by line into temporary array
        do
            read(unit, formatstr, iostat=ios) row
            if (ios /= 0) exit  ! Reached end-of-file
            nrows = nrows + 1
            if (nrows > maxrows) then
                print *, "Error: exceeded maximum number of rows =", maxrows
                stop 1
            end if
            temp(nrows,1:ncol) = row
        end do

        close(unit)

        ! Allocate exact-size array and copy data
        allocate(A(nrows,ncol))
        A(:,:) = temp(1:nrows,1:ncol)
        deallocate(temp)
    end subroutine read_matrix_with_format


    ! Compares the six Hams_format files of one matrix (kind = 'OAMASS' or 'ODAMPING',
    ! files kind1.txt to kind6.txt) in dirB against the reference in dirA.
    ! File i, row (frequency f, body pair a-b), column j holds M_ij^ab(f).
    ! Each entry is scaled by the self-block diagonals of the reference:
    !     max |B_ij^ab - A_ij^ab| / sqrt(D_i^a D_j^b) < tolerance,  D_i^a = |A_ii^aa|,
    ! with D floored at diag_floor times the largest diagonal at that frequency.
    ! A norm over a whole file is not used: couplings that symmetry makes zero
    ! (e.g. surge-heave in a mirror-symmetric array) carry solver noise of ~1e-4
    ! of the heave term, can exceed the file's physical entries and dominate its
    ! norm. Scaled this way they count at their true size, ~1e-6.
    ! The floor keeps near-zero diagonals (yaw of an axisymmetric hull) from
    ! blowing up the ratio.
    subroutine test_hams_matrices_are_close(dirA, dirB, kind, tolerance, error)
        implicit none
        character(len=*), intent(in) :: dirA, dirB, kind
        real(kind=8), intent(in) :: tolerance
        type(error_type), allocatable, intent(out) :: error
        real(kind=8), parameter :: diag_floor = 1.0d-6
        character(len=*), parameter :: fmt = '(F7.3,1X,F7.3,1X,6E14.5)'
        real(kind=8), allocatable :: F(:,:), MA(:,:,:), MB(:,:,:), dg(:,:)
        real(kind=8) :: d, dmax, worst, worst_freq
        integer :: i, j, nrow, npair, nbody, nfreq, ifreq, ipair, a, b, r0, worst_loc(4)
        character(len=1) :: istr
        character(len=16) :: worst_str, tolerance_str, freq_str
        character(len=64) :: where_str

        do i = 1, 6
            write(istr, '(I1)') i
            call read_matrix_with_format(dirA//kind//istr//'.txt', fmt, 8, F)
            if (i == 1) then
                nrow = size(F, 1)
                allocate(MA(nrow, 6, 6), MB(nrow, 6, 6))
            end if
            call check(error, size(F, 1) == nrow, "Row count differs in "//dirA//kind//istr//'.txt')
            if (allocated(error)) return
            MA(:, i, :) = F(:, 3:8)
            call read_matrix_with_format(dirB//kind//istr//'.txt', fmt, 8, F)
            call check(error, size(F, 1) == nrow, "Row count differs in "//dirB//kind//istr//'.txt')
            if (allocated(error)) return
            MB(:, i, :) = F(:, 3:8)
        end do

        ! Rows per frequency (column 2) give the number of body pairs
        call read_matrix_with_format(dirA//kind//'1.txt', fmt, 8, F)
        npair = 1
        do while (npair < nrow)
            if (F(npair + 1, 2) /= F(1, 2)) exit
            npair = npair + 1
        end do
        nbody = nint(sqrt(real(npair, 8)))
        call check(error, nbody * nbody == npair .and. mod(nrow, npair) == 0, &
                   "Cannot read the body pairs of "//dirA//kind//'1.txt')
        if (allocated(error)) return
        nfreq = nrow / npair

        allocate(dg(6, nbody))
        worst = -1.0d0
        worst_loc = 0
        worst_freq = 0.0d0
        do ifreq = 1, nfreq
            r0 = (ifreq - 1) * npair
            do a = 1, nbody
                do i = 1, 6
                    dg(i, a) = abs(MA(r0 + (a - 1) * nbody + a, i, i))
                end do
            end do
            dmax = maxval(dg)
            if (dmax <= 0.0d0) dmax = 1.0d0
            dg = max(dg, diag_floor * dmax)
            do ipair = 1, npair
                a = (ipair - 1) / nbody + 1
                b = mod(ipair - 1, nbody) + 1
                do i = 1, 6
                    do j = 1, 6
                        d = abs(MB(r0 + ipair, i, j) - MA(r0 + ipair, i, j)) / sqrt(dg(i, a) * dg(j, b))
                        if (d > worst) then
                            worst = d
                            worst_loc = [a, b, i, j]
                            worst_freq = F(r0 + ipair, 2)
                        end if
                    end do
                end do
            end do
        end do

        write(worst_str, '(ES10.3)') worst
        write(tolerance_str, '(ES10.3)') tolerance
        write(freq_str, '(F7.3)') worst_freq
        write(where_str, '(A,I0,A,I0,A,I0,A,I0)') "bodies ", worst_loc(1), "-", worst_loc(2), &
            ", modes ", worst_loc(3), "-", worst_loc(4)
        call check(error, worst < tolerance, "Scaled difference "//trim(worst_str)// &
                   " is not below tolerance "//trim(tolerance_str)//" ("//trim(where_str)// &
                   ", column 2 = "//trim(adjustl(freq_str))//") for "//dirB//kind//"1-6.txt")
    end subroutine test_hams_matrices_are_close

end module testing_utilities