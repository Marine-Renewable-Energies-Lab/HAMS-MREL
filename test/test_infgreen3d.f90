!> Unit test for the infinite-depth Green function (InfGreen_Appr.f90).
!> INFGREEN3D is called from OpenMP regions, so results must not depend on the number of
!> threads. The module scratch variables of HavelockGF must be THREADPRIVATE for this to hold.
module test_infgreen3d
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use INFG3D_Open, only : INFGREEN3D
  implicit none
  private
  public :: collect_infgreen3d

  integer, parameter :: NPTS = 20000
  integer, parameter :: NREPEAT = 5
  integer, parameter :: NTHREADS = 4
  real(8), parameter :: TOL = 1d-10

contains

  subroutine collect_infgreen3d(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ new_unittest("openmp_matches_serial", test_openmp_matches_serial) ]
  end subroutine collect_infgreen3d

  subroutine test_openmp_matches_serial(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8), allocatable :: serial(:,:), threaded(:,:)
    real(8), parameter :: V = 0.2d0
    integer :: i, irep

    allocate(serial(4,NPTS), threaded(4,NPTS))
    do i = 1, NPTS
      call evaluate(i, V, serial(:,i))
    end do

    do irep = 1, NREPEAT
      !$OMP PARALLEL DO NUM_THREADS(NTHREADS) PRIVATE(i) SCHEDULE(DYNAMIC,16)
      do i = 1, NPTS
        call evaluate(i, V, threaded(:,i))
      end do
      !$OMP END PARALLEL DO
      ! Not exact equality: worker threads can differ by round-off (seen: 3e-14 with gfortran on Windows).
      ! A race gives errors many orders of magnitude larger.
      call check(error, all(abs(threaded - serial) <= TOL*(abs(serial) + 1d-6)), &
                 "INFGREEN3D gives different results when called from several threads")
      if (allocated(error)) return
    end do
  end subroutine test_openmp_matches_serial

  !> Green function at a field point that varies with i, covering a wide range of R and z.
  subroutine evaluate(i, V, grn)
    integer, intent(in) :: i
    real(8), intent(in) :: V
    complex(8), intent(out) :: grn(4)
    real(8) :: x, zf
    x  = 0.1d0 + 40d0*mod(i, 97)/97d0
    zf = -0.5d0 - 30d0*mod(i, 89)/89d0
    call INFGREEN3D(x, 0d0, 0d0, 0d0, zf, -1d0, V, grn, 1)
  end subroutine evaluate

end module test_infgreen3d
