!> The 1/r integrals over damping lid panels, as sources, against the composite Gauss quadrature of test_panel_integrals. The
!> lid arrays are module data, so this suite runs serially.
module test_damping_lid

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use Lid_mod, only : LNCN, LNCON, LXYZ
  use SingularIntgr, only : SGLINTBD_TRI2, SGLINTBD_QUAD2
  use SingularIntgrMulti, only : SGLINTBD_MULTI_TRI2, SGLINTBD_MULTI_QUAD2
  use test_panel_integrals, only : reference

  implicit none
  private
  public :: collect_damping_lid

  real(8), parameter :: tol = 1d-8

contains

  subroutine collect_damping_lid(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("lid_rankine_quad", test_lid_rankine_quad), &
      new_unittest("lid_rankine_triangle", test_lid_rankine_triangle), &
      new_unittest("lid_node_order", test_lid_node_order), &
      new_unittest("lid_multi_matches_single", test_lid_multi_matches_single) &
      ]
  end subroutine collect_damping_lid

  !> One lid panel, in the lid arrays.
  subroutine set_lid(nodes)
    real(8), intent(in) :: nodes(:, :)
    integer :: k, n

    n = size(nodes, 1)
    if (allocated(LXYZ)) deallocate (LXYZ)
    if (allocated(LNCN)) deallocate (LNCN)
    if (allocated(LNCON)) deallocate (LNCON)
    allocate (LXYZ(n, 3), LNCN(1), LNCON(1, 4))
    LXYZ = nodes
    LNCN(1) = n
    LNCON = 0
    do k = 1, n
      LNCON(1, k) = k
    end do
  end subroutine set_lid

  subroutine lid_integral(nodes, p, rk, multi)
    real(8), intent(in) :: nodes(:, :), p(3)
    real(8), intent(out) :: rk(4)
    logical, intent(in), optional :: multi
    logical :: m

    m = .false.
    if (present(multi)) m = multi
    call set_lid(nodes)
    if (size(nodes, 1) == 3) then
      if (m) then
        call SGLINTBD_MULTI_TRI2(p, 1, rk, 5)
      else
        call SGLINTBD_TRI2(p, 1, rk, 5)
      end if
    else
      if (m) then
        call SGLINTBD_MULTI_QUAD2(p, 1, rk, 5)
      else
        call SGLINTBD_QUAD2(p, 1, rk, 5)
      end if
    end if
  end subroutine lid_integral

  subroutine test_lid_rankine_quad(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: nodes(4, 3), points(4, 3)

    nodes = reshape([0d0, 2d0, 2d0, 0d0, 0d0, 0d0, 1d0, 1d0, 0d0, 0d0, 0d0, 0d0], [4, 3])
    points(1, :) = [1d0, 0.5d0, 0d0]
    points(2, :) = [2d0, 0d0, 0d0]
    points(3, :) = [0.3d0, 0.2d0, 0d0]
    points(4, :) = [5d0, 0d0, -1d0]
    ! The field points on the panel are singular for the reference, which is only for the exterior one
    call check_points_mixed(error, nodes, points)
  end subroutine test_lid_rankine_quad

  subroutine test_lid_rankine_triangle(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: nodes(3, 3), points(4, 3)

    ! Obtuse at node 2
    nodes = reshape([0d0, 1d0, 1.8d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0], [3, 3])
    points(1, :) = [(0d0 + 1d0 + 1.8d0)/3d0, 1d0/3d0, 0d0]
    points(2, :) = [1d0, 0d0, 0d0]
    points(3, :) = [1.2d0, 0.3d0, 0d0]
    points(4, :) = [5d0, 0d0, -1d0]
    call check_points_mixed(error, nodes, points)
  end subroutine test_lid_rankine_triangle

  !> The points in the panel itself are checked by the fan of triangles through the point, the others by the plain quadrature.
  subroutine check_points_mixed(error, nodes, points)
    type(error_type), allocatable, intent(out) :: error
    real(8), intent(in) :: nodes(:, :), points(:, :)
    real(8) :: rk(4), s, d(3), fan(3, 3), si
    integer :: i, k, n

    n = size(nodes, 1)
    do i = 1, size(points, 1)
      call lid_integral(nodes, points(i, :), rk)
      if (abs(points(i, 3)) > 0d0) then
        call reference(nodes, points(i, :), 48, s, d)
      else
        s = 0d0
        do k = 1, n
          fan(1, :) = nodes(k, :)
          fan(2, :) = nodes(mod(k, n) + 1, :)
          fan(3, :) = points(i, :)
          call reference(fan, points(i, :), 64, si, d)
          s = s + si
        end do
      end if
      call check(error, rk(1), s, "point "//char(48+i), rel=.true., thr=tol)
      if (allocated(error)) return
    end do
  end subroutine check_points_mixed

  subroutine test_lid_node_order(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: nodes(4, 3), perm(4, 3), rk(4), rk0(4), p(3)
    integer :: shift, k, rev

    nodes = reshape([0d0, 2d0, 2d0, 0d0, 0d0, 0d0, 1d0, 1d0, 0d0, 0d0, 0d0, 0d0], [4, 3])
    p = [5d0, 0d0, -1d0]
    call lid_integral(nodes, p, rk0)
    do rev = 0, 1
      do shift = 0, 3
        do k = 1, 4
          if (rev == 0) then
            perm(k, :) = nodes(mod(k - 1 + shift, 4) + 1, :)
          else
            perm(k, :) = nodes(mod(-(k - 1) + shift + 4, 4) + 1, :)
          end if
        end do
        call lid_integral(perm, p, rk)
        call check(error, rk(1), rk0(1), rel=.true., thr=1d-12)
        if (allocated(error)) return
      end do
    end do
  end subroutine test_lid_node_order

  subroutine test_lid_multi_matches_single(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: quad(4, 3), tri(3, 3), p(3), a(4), b(4)

    quad = reshape([0d0, 2d0, 2d0, 0d0, 0d0, 0d0, 1d0, 1d0, 0d0, 0d0, 0d0, 0d0], [4, 3])
    tri = reshape([0d0, 1d0, 1.8d0, 0d0, 0d0, 1d0, 0d0, 0d0, 0d0], [3, 3])
    p = [0.3d0, 0.2d0, -0.7d0]
    call lid_integral(quad, p, a)
    call lid_integral(quad, p, b, multi=.true.)
    call check(error, all(a == b), "quad: multi differs from single")
    if (allocated(error)) return
    call lid_integral(tri, p, a)
    call lid_integral(tri, p, b, multi=.true.)
    call check(error, all(a == b), "triangle: multi differs from single")
  end subroutine test_lid_multi_matches_single

end module test_damping_lid
