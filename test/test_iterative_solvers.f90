!> Unit tests for GMRES (IterativeSolvers.f90) and H-matrix compression (HMatrix.f90).
!>  1. A solve reported as converged must have a true residual at the solver tolerance, for
!>     general complex matrices and for the complex normal equations of the IRR path.
!>  2. H-matrix products must match the dense product, including rectangular blocks.
module test_iterative_solvers
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use IterativeSolvers, only : ZITER_SOLVE_MULTI, ITER_TOL
  use HMatrix_mod, only : HMATRIX_BUILD, HMATRIX_DESTROY, HMATVEC
  implicit none
  private
  public :: collect_iterative_solvers

  integer, parameter :: N = 300
  ! The residual estimate is taken from the Krylov recurrence, so allow a small margin over ITER_TOL.
  real(8), parameter :: RESIDUAL_MARGIN = 2d0

contains

  subroutine collect_iterative_solvers(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("gmres_true_residual_general", test_gmres_general), &
      new_unittest("gmres_true_residual_normal_equations", test_gmres_normal_equations), &
      new_unittest("hmatrix_matches_dense", test_hmatrix_matches_dense) ]
  end subroutine collect_iterative_solvers

  subroutine test_gmres_general(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8), allocatable :: A(:,:)
    integer :: iphase
    allocate(A(N,N))
    do iphase = 0, 3
      call make_matrix(0.5d0*iphase, 0.6d0, A)
      call check_solve(A, error)
      if (allocated(error)) return
    end do
  end subroutine test_gmres_general

  subroutine test_gmres_normal_equations(error)
    type(error_type), allocatable, intent(out) :: error
    complex(8), allocatable :: A(:,:), C(:,:)
    integer :: iphase
    allocate(A(N,N), C(N,N))
    do iphase = 0, 3
      call make_matrix(0.5d0*iphase, 0.3d0, A)
      C = matmul(transpose(A), A)
      call check_solve(C, error)
      if (allocated(error)) return
    end do
  end subroutine test_gmres_normal_equations

  subroutine check_solve(A, error)
    complex(8), intent(in) :: A(:,:)
    type(error_type), allocatable, intent(out) :: error
    complex(8), allocatable :: B(:,:), X(:,:), R(:,:)
    integer :: i, info
    allocate(B(N,1), X(N,1), R(N,1))
    do i = 1, N
      B(i,1) = dcmplx(sin(1.3d0*i), cos(0.7d0*i))
    end do
    call ZITER_SOLVE_MULTI(2, N, A, B, X, 1, info)
    call check(error, info == 0, "GMRES did not converge on a well-conditioned system")
    if (allocated(error)) return
    R = B - matmul(A, X)
    call check(error, sqrt(sum(abs(R)**2)) < RESIDUAL_MARGIN*ITER_TOL*sqrt(sum(abs(B)**2)), &
               "GMRES reported convergence but the true residual is above the tolerance")
  end subroutine check_solve

  !> Well-conditioned complex matrix: rotated unit diagonal plus a deterministic perturbation.
  subroutine make_matrix(phase, scale, A)
    real(8), intent(in) :: phase, scale
    complex(8), intent(out) :: A(:,:)
    integer :: i, j
    do j = 1, N
      do i = 1, N
        A(i,j) = scale*dcmplx(sin(12.9898d0*i + 78.233d0*j), cos(39.346d0*i + 11.135d0*j))*sqrt(6d0/N)
      end do
    end do
    do i = 1, N
      A(i,i) = A(i,i) + exp(dcmplx(0d0, phase))
    end do
  end subroutine make_matrix

  !> Panels in two groups of different sizes, so admissible blocks are rectangular both ways.
  subroutine test_hmatrix_matches_dense(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: N1 = 500, N2 = 200, M = N1 + N2
    real(8), allocatable :: XYZ(:,:)
    complex(8), allocatable :: A(:,:), X(:), Y(:), Yref(:)
    real(8) :: t, d
    integer :: i, j

    allocate(XYZ(M,3), A(M,M), X(M), Y(M), Yref(M))
    do i = 1, M
      t = 0.37d0*i
      XYZ(i,:) = [cos(t)*(1d0 + mod(i,17)/17d0), sin(t)*(1d0 + mod(i,13)/13d0), -mod(i,29)/29d0]
      if (i > N1) XYZ(i,1) = XYZ(i,1) + 50d0
      X(i) = dcmplx(sin(1.1d0*i), cos(2.3d0*i))
    end do
    do j = 1, M
      do i = 1, M
        d = sqrt(sum((XYZ(i,:) - XYZ(j,:))**2))
        A(i,j) = exp(dcmplx(0d0, 0.2d0*d))/(1d0 + d)
      end do
    end do
    Yref = matmul(A, X)

    call HMATRIX_BUILD(M, A, XYZ)
    call HMATVEC(M, X, Y)
    call HMATRIX_DESTROY()

    call check(error, sqrt(sum(abs(Y - Yref)**2)) < 1d-2*sqrt(sum(abs(Yref)**2)), &
               "H-matrix product differs from the dense product")
  end subroutine test_hmatrix_matches_dense

end module test_iterative_solvers
