!> The lid (IRSP = 1) system is overdetermined: hull and lid rows for the hull potentials only. AssbMatx_irr and
!> AssbMatx_irrMulti solve it through the normal equations, built by NORMAL_MATRIX and NORMAL_RHS. For a complex matrix the
!> least-squares solution, the minimum of ||A x - b||, solves A^H A x = A^H b (Golub and Van Loan, Matrix Computations,
!> section 5.3). LAPACK's ZGELS solves the same least-squares problem by QR, without the normal equations.
module test_normal_equations

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use AssbMatx_irr, only : NORMAL_MATRIX, NORMAL_RHS
  use HAMS_mod, only : NTHREAD

  implicit none
  private
  public :: collect_normal_equations

  integer, parameter :: m = 9, n = 4, nrhs = 2

contains

  subroutine collect_normal_equations(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("solve_equals_least_squares", test_solve_equals_least_squares), &
      new_unittest("matrix_is_hermitian", test_matrix_is_hermitian), &
      new_unittest("real_matrix_unchanged", test_real_matrix_unchanged) &
      ]
  end subroutine collect_normal_equations

  !> A complex, full-rank, overdetermined system whose rows do not agree (b is not in the range of A).
  subroutine system(a, b)
    complex(8), intent(out) :: a(m, n), b(m, nrhs)
    integer :: i, j

    do j = 1, n
      do i = 1, m
        a(i, j) = cmplx(cos(0.7d0*i + 1.3d0*j) + merge(2d0, 0d0, i == j), sin(0.4d0*i*j + 0.2d0), kind=8)
      end do
    end do
    do j = 1, nrhs
      do i = 1, m
        b(i, j) = cmplx(sin(0.9d0*i + j), cos(0.3d0*i*j), kind=8)
      end do
    end do
  end subroutine system

  subroutine test_solve_equals_least_squares(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: a(m, n), b(m, nrhs), c(n, n), d(n, nrhs), qa(m, n), qb(m, nrhs), work(256)
    integer :: ipiv(n), info
    external :: ZGESV, ZGELS

    NTHREAD = 1
    call system(a, b)
    call NORMAL_MATRIX(a, m, n, c)
    call NORMAL_RHS(a, m, n, b, nrhs, d)
    call ZGESV(n, nrhs, c, n, ipiv, d, n, info)
    call check(error, info == 0, "ZGESV failed on the normal equations")
    if (allocated(error)) return
    qa = a
    qb = b
    call ZGELS('N', m, n, nrhs, qa, m, qb, m, work, size(work), info)
    call check(error, info == 0, "ZGELS failed")
    if (allocated(error)) return
    write (*, '(4x, a, es10.3)') "normal equations against ZGELS, max |dx|/max |x|: ", &
      maxval(abs(d - qb(1:n, :)))/maxval(abs(qb(1:n, :)))
    call check(error, maxval(abs(d - qb(1:n, :))) < 1d-10*maxval(abs(qb(1:n, :))), &
               "the normal-equation solution is not the least-squares solution")
  end subroutine test_solve_equals_least_squares

  !> A^H A is Hermitian (and positive definite for a full-rank A); A^T A is only complex symmetric.
  subroutine test_matrix_is_hermitian(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: a(m, n), b(m, nrhs), c(n, n)

    NTHREAD = 1
    call system(a, b)
    call NORMAL_MATRIX(a, m, n, c)
    call check(error, maxval(abs(c - conjg(transpose(c)))) < 1d-12*maxval(abs(c)), "the normal matrix is not Hermitian")
  end subroutine test_matrix_is_hermitian

  !> For a real matrix (the zero- and infinite-frequency limits) nothing changes: A^H A = A^T A. Compared to round-off,
  !> as the compiler may order the sums differently; the integration suite pins the limits of whole runs.
  subroutine test_real_matrix_unchanged(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8) :: a(m, n), b(m, nrhs), c(n, n), d(n, nrhs), s
    integer :: i, j, k

    NTHREAD = 1
    call system(a, b)
    a = cmplx(real(a), 0d0, kind=8)
    call NORMAL_MATRIX(a, m, n, c)
    call NORMAL_RHS(a, m, n, b, nrhs, d)
    do i = 1, n
      do j = 1, n
        s = (0d0, 0d0)
        do k = 1, m
          s = s + a(k, i)*a(k, j)
        end do
        call check(error, abs(c(i, j) - s) <= 1d-14*abs(s), "real matrix: the normal matrix changed")
        if (allocated(error)) return
      end do
      do j = 1, nrhs
        s = (0d0, 0d0)
        do k = 1, m
          s = s + b(k, j)*a(k, i)
        end do
        call check(error, abs(d(i, j) - s) <= 1d-14*abs(s), "real matrix: the right side changed")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_real_matrix_unchanged

end module test_normal_equations
