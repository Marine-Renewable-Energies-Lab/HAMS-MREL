!> The panel integral of the free-surface image term 1/r1 and of its gradient, against an independent composite Gauss
!> quadrature in this file. It uses the mesh arrays of the solver, which are module globals, so this suite runs serially.
module test_image_term

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use PanelMesh_mod, only : XYZ, NCN, NCON, PNSZ, XYZ_GLOBAL_MULTI_COMB, NCN_MULTI_COMB, NCON_MULTI_COMB, PNSZ_MULTI_COMB
  use SingularIntgr, only : IMAGE_TERM_EXACT
  use SingularIntgrMulti, only : IMAGE_TERM_EXACT_MULTI
  use FinGreen3D_Open, only : FINGREEN3D
  use INFG3D_Open, only : INFGREEN3D

  implicit none
  private
  public :: collect_image_term

  real(8), parameter :: tol = 1d-8
  ! Nodes and weights of the 8-point Gauss-Legendre rule on (-1, 1), half of them, the others by symmetry
  real(8), parameter :: gx(4) = [0.1834346424956498d0, 0.5255324099163290d0, 0.7966664774136267d0, 0.9602898564975363d0]
  real(8), parameter :: gw(4) = [0.3626837833783620d0, 0.3137066458778873d0, 0.2223810344533745d0, 0.1012285362903763d0]
  integer, parameter :: cells = 64

contains

  subroutine collect_image_term(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("thin_panel_self", test_thin_panel_self), &
      new_unittest("thin_panel_neighbours", test_thin_panel_neighbours), &
      new_unittest("square_panel", test_square_panel), &
      new_unittest("horizontal_panel_gradient", test_horizontal_panel_gradient), &
      new_unittest("triangle", test_triangle), &
      new_unittest("repeated_node_quad", test_repeated_node_quad), &
      new_unittest("infinite_frequency_sign", test_infinite_frequency_sign), &
      new_unittest("wave_part_kept", test_wave_part_kept), &
      new_unittest("multi_matches_single", test_multi_matches_single) &
      ]
  end subroutine collect_image_term

  ! ---------------------------------------------------------------------------------------------- helpers

  !> One panel in the arrays of the single-body path. nodes(k, :) is node k.
  subroutine set_mesh(nodes)
    real(8), intent(in) :: nodes(:, :)
    integer :: k

    if (allocated(XYZ)) deallocate (XYZ, NCN, NCON, PNSZ)
    allocate (XYZ(size(nodes, 1), 3), NCN(1), NCON(1, 4), PNSZ(1))
    XYZ = nodes
    NCN(1) = size(nodes, 1)
    NCON = 0
    do k = 1, size(nodes, 1)
      NCON(1, k) = k
    end do
    PNSZ(1) = longest_side(nodes)
  end subroutine set_mesh

  real(8) function longest_side(nodes)
    real(8), intent(in) :: nodes(:, :)
    integer :: k, n

    n = size(nodes, 1)
    longest_side = 0d0
    do k = 1, n
      longest_side = max(longest_side, norm2(nodes(k, :) - nodes(mod(k, n) + 1, :)))
    end do
  end function longest_side

  function centroid(nodes) result(c)
    real(8), intent(in) :: nodes(:, :)
    real(8) :: c(3)
    c = sum(nodes, dim=1)/real(size(nodes, 1), 8)
  end function centroid

  !> Integrals of 1/r1 and of d(1/r1)/d(source) over the panel, source q, r1 the distance to the mirror image of field point p
  !> in z = 0: composite Gauss on the bilinear map of the square, or of the triangle collapsed to a square.
  subroutine reference(nodes, p, s, d)
    real(8), intent(in) :: nodes(:, :), p(3)
    real(8), intent(out) :: s, d(3)
    real(8) :: a(4, 3), xa, xb, wa, wb, q(3), du(3), dv(3), jac, r, m(3), n(3)
    integer :: i, j, ia, ib

    a(1:size(nodes, 1), :) = nodes
    if (size(nodes, 1) == 3) a(4, :) = nodes(3, :)
    m = [p(1), p(2), -p(3)]
    s = 0d0
    d = 0d0
    do i = 1, cells
      do j = 1, cells
        do ia = -4, 4
          if (ia == 0) cycle
          xa = (real(i - 1, 8) + 0.5d0*(1d0 + sign(gx(abs(ia)), real(ia, 8))))/real(cells, 8)
          wa = 0.5d0*gw(abs(ia))/real(cells, 8)
          do ib = -4, 4
            if (ib == 0) cycle
            xb = (real(j - 1, 8) + 0.5d0*(1d0 + sign(gx(abs(ib)), real(ib, 8))))/real(cells, 8)
            wb = 0.5d0*gw(abs(ib))/real(cells, 8)
            q = (1d0 - xa)*(1d0 - xb)*a(1, :) + xa*(1d0 - xb)*a(2, :) + xa*xb*a(3, :) + (1d0 - xa)*xb*a(4, :)
            du = (1d0 - xb)*(a(2, :) - a(1, :)) + xb*(a(3, :) - a(4, :))
            dv = (1d0 - xa)*(a(4, :) - a(1, :)) + xa*(a(3, :) - a(2, :))
            n = [du(2)*dv(3) - du(3)*dv(2), du(3)*dv(1) - du(1)*dv(3), du(1)*dv(2) - du(2)*dv(1)]
            jac = norm2(n)
            r = norm2(q - m)
            s = s + wa*wb*jac/r
            d = d - wa*wb*jac*(q - m)/r**3
          end do
        end do
      end do
    end do
  end subroutine reference

  !> With V = 0 the Green function holds the one-point image and nothing else, so after the exchange all of it has moved to the
  !> panel integrals, which must match the quadrature, and nothing remains.
  subroutine check_panel(nodes, p, name, error)
    real(8), intent(in) :: nodes(:, :), p(3)
    character(len=*), intent(in) :: name
    type(error_type), allocatable, intent(out) :: error
    real(8) :: s, d(3), sij, dij(3), q(3), wvn(1)
    complex(8) :: grn(4)

    call set_mesh(nodes)
    call reference(nodes, p, s, d)
    q = centroid(nodes)
    wvn = 1d0
    call FINGREEN3D(q(1), p(1), q(2), p(2), q(3), p(3), 0d0, wvn, 1, 100d0, grn, 1)
    sij = 0d0
    dij = 0d0
    call IMAGE_TERM_EXACT(p, q, 0d0, 1, 1, grn, sij, dij)
    call check(error, abs(sij - s) < tol*max(abs(s), 1d-3), name//": integral of 1/r1")
    if (allocated(error)) return
    call check(error, maxval(abs(dij - d)) < tol*max(maxval(abs(d)), 1d-3), name//": integral of the gradient of 1/r1")
    if (allocated(error)) return
    call check(error, maxval(abs(grn)) < 1d-12, name//": the one-point image is left in the Green function")
  end subroutine check_panel

  !> A wall panel at x = 0, with its top edge at z = top, width w along y, height h.
  function wall_panel(w, h, top) result(nodes)
    real(8), intent(in) :: w, h, top
    real(8) :: nodes(4, 3)
    nodes(1, :) = [0d0, -0.5d0*w, top]
    nodes(2, :) = [0d0, 0.5d0*w, top]
    nodes(3, :) = [0d0, 0.5d0*w, top - h]
    nodes(4, :) = [0d0, -0.5d0*w, top - h]
  end function wall_panel

  ! ---------------------------------------------------------------------------------------------- tests

  !> The top row of a wall graded to the waterline. The one-point value here is 4.7 times the integral.
  subroutine test_thin_panel_self(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: nodes(4, 3)

    nodes = wall_panel(0.654d0, 0.0196d0, 0d0)
    call check_panel(nodes, centroid(nodes), "thin panel, its own centroid", error)
  end subroutine test_thin_panel_self

  !> Next panel along the row, the next row below, and a point off the wall.
  subroutine test_thin_panel_neighbours(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: nodes(4, 3)

    nodes = wall_panel(0.654d0, 0.0196d0, 0d0)
    call check_panel(nodes, [0d0, 0.654d0, -0.0098d0], "next in the row", error)
    if (allocated(error)) return
    call check_panel(nodes, [0d0, 0d0, -0.0489d0], "next row", error)
    if (allocated(error)) return
    call check_panel(nodes, [0.05d0, 0.02d0, -0.03d0], "off the wall", error)
    if (allocated(error)) return
    ! Field point at z = 0: the image coincides with the point itself
    call check_panel(wall_panel(0.654d0, 0.0196d0, -0.01d0), [0.2d0, 0.1d0, 0d0], "field point on the free surface", error)
  end subroutine test_thin_panel_neighbours

  subroutine test_square_panel(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: nodes(4, 3)

    nodes = wall_panel(1d0, 1d0, -2d0)
    call check_panel(nodes, centroid(nodes), "square, deep", error)
    if (allocated(error)) return
    nodes = wall_panel(1d0, 1d0, 0d0)
    call check_panel(nodes, centroid(nodes), "square, at the waterline", error)
  end subroutine test_square_panel

  !> A panel of the bottom of a shallow hull: the image is off its plane, so d/dz is not zero.
  subroutine test_horizontal_panel_gradient(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: nodes(4, 3)

    nodes(1, :) = [0d0, 0d0, -0.1d0]
    nodes(2, :) = [2d0, 0d0, -0.1d0]
    nodes(3, :) = [2d0, 1d0, -0.1d0]
    nodes(4, :) = [0d0, 1d0, -0.1d0]
    call check_panel(nodes, centroid(nodes), "horizontal, its own centroid", error)
    if (allocated(error)) return
    call check_panel(nodes, [2.5d0, 0.5d0, -0.1d0], "horizontal, next panel", error)
  end subroutine test_horizontal_panel_gradient

  subroutine test_triangle(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: nodes(3, 3)

    nodes(1, :) = [0d0, -0.3d0, 0d0]
    nodes(2, :) = [0d0, 0.3d0, 0d0]
    nodes(3, :) = [0d0, 0d0, -0.03d0]
    call check_panel(nodes, centroid(nodes), "thin triangle, its own centroid", error)
    if (allocated(error)) return
    call check_panel(nodes, [0.1d0, 0.5d0, -0.02d0], "thin triangle, a point nearby", error)
  end subroutine test_triangle

  !> A triangle written as a quadrilateral, with its last node twice, as meshes may do.
  subroutine test_repeated_node_quad(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: tri(3, 3), nodes(4, 3), p(3), q(3), s, d(3), sij, dij(3), wvn(1)
    complex(8) :: grn(4)

    tri(1, :) = [0d0, -0.3d0, 0d0]
    tri(2, :) = [0d0, 0.3d0, 0d0]
    tri(3, :) = [0d0, 0d0, -0.03d0]
    nodes(1:3, :) = tri
    nodes(4, :) = tri(3, :)
    p = [0.1d0, 0.5d0, -0.02d0]
    call reference(tri, p, s, d)
    call set_mesh(nodes)
    q = centroid(tri)
    wvn = 1d0
    call FINGREEN3D(q(1), p(1), q(2), p(2), q(3), p(3), 0d0, wvn, 1, 100d0, grn, 1)
    sij = 0d0
    dij = 0d0
    call IMAGE_TERM_EXACT(p, q, 0d0, 1, 1, grn, sij, dij)
    call check(error, abs(sij - s) < tol*abs(s), "integral of 1/r1")
    if (allocated(error)) return
    call check(error, maxval(abs(dij - d)) < tol*maxval(abs(d)), "integral of the gradient of 1/r1")
  end subroutine test_repeated_node_quad

  !> The limit of infinite frequency has the image with the opposite sign, on both depth branches.
  subroutine test_infinite_frequency_sign(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: nodes(4, 3), p(3), q(3), s, d(3), sij, dij(3), wvn(1)
    complex(8) :: grn(4)
    integer :: branch

    nodes = wall_panel(0.654d0, 0.0196d0, 0d0)
    call set_mesh(nodes)
    q = centroid(nodes)
    p = q
    call reference(nodes, p, s, d)
    wvn = 1d0
    do branch = 1, 2
      if (branch == 1) then
        call FINGREEN3D(q(1), p(1), q(2), p(2), q(3), p(3), -1d0, wvn, 1, 100d0, grn, 1)
      else
        call INFGREEN3D(q(1), p(1), q(2), p(2), q(3), p(3), -1d0, grn, 1)
      end if
      sij = 0d0
      dij = 0d0
      call IMAGE_TERM_EXACT(p, q, -1d0, 1, 1, grn, sij, dij)
      call check(error, abs(sij + s) < tol*abs(s), "integral of 1/r1 with the sign of the limit")
      if (allocated(error)) return
      call check(error, maxval(abs(dij + d)) < tol*max(maxval(abs(d)), 1d-3), "gradient with the sign of the limit")
      if (allocated(error)) return
      call check(error, maxval(abs(grn)) < 1d-12, "the one-point image is left in the Green function")
      if (allocated(error)) return
    end do
  end subroutine test_infinite_frequency_sign

  !> At a real frequency only 1/r1 moves: the wave part is as the Green function gave it.
  subroutine test_wave_part_kept(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: nodes(4, 3), p(3), q(3), sij, dij(3), r1, pq(3)
    complex(8) :: grn(4), before(4)

    nodes = wall_panel(0.654d0, 0.0196d0, 0d0)
    call set_mesh(nodes)
    q = centroid(nodes)
    p = [0d0, 0.654d0, q(3)]
    pq = [q(1) - p(1), q(2) - p(2), q(3) + p(3)]
    r1 = norm2(pq)
    call INFGREEN3D(q(1), p(1), q(2), p(2), q(3), p(3), 0.12d0, before, 1)
    grn = before
    sij = 0d0
    dij = 0d0
    call IMAGE_TERM_EXACT(p, q, 0.12d0, 1, 1, grn, sij, dij)
    call check(error, abs(grn(1) - (before(1) - 1d0/r1)) < 1d-12, "wave part of G")
    if (allocated(error)) return
    call check(error, maxval(abs(grn(2:4) - (before(2:4) + pq/r1**3))) < 1d-9, "wave part of the gradient")
  end subroutine test_wave_part_kept

  !> The two paths have their own copies of the panel integrals and of the mesh arrays.
  subroutine test_multi_matches_single(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: nodes(4, 3), tri(3, 3), p(3), q(3), s1, d1(3), s2, d2(3)
    complex(8) :: g1(4), g2(4)
    integer :: k, case

    nodes = wall_panel(0.654d0, 0.0196d0, 0d0)
    tri(1, :) = [0d0, -0.3d0, 0d0]
    tri(2, :) = [0d0, 0.3d0, 0d0]
    tri(3, :) = [0d0, 0d0, -0.03d0]
    do case = 1, 2
      if (case == 1) then
        call set_mesh(nodes)
        q = centroid(nodes)
      else
        call set_mesh(tri)
        q = centroid(tri)
      end if
      if (allocated(XYZ_GLOBAL_MULTI_COMB)) deallocate (XYZ_GLOBAL_MULTI_COMB, NCN_MULTI_COMB, NCON_MULTI_COMB, PNSZ_MULTI_COMB)
      allocate (XYZ_GLOBAL_MULTI_COMB(size(XYZ, 1), 3), NCN_MULTI_COMB(1), NCON_MULTI_COMB(1, 4), PNSZ_MULTI_COMB(1))
      XYZ_GLOBAL_MULTI_COMB = XYZ
      NCN_MULTI_COMB = NCN
      NCON_MULTI_COMB = NCON
      PNSZ_MULTI_COMB = PNSZ
      do k = 1, 2
        p = q + real(k - 1, 8)*[0.2d0, 0.3d0, -0.05d0]
        g1 = (1d0, 0d0)
        g2 = (1d0, 0d0)
        s1 = 0d0; d1 = 0d0; s2 = 0d0; d2 = 0d0
        call IMAGE_TERM_EXACT(p, q, 0.12d0, 1, 1, g1, s1, d1)
        call IMAGE_TERM_EXACT_MULTI(p, q, 0.12d0, 1, 1, g2, s2, d2)
        call check(error, s1 /= 0d0, "the single-body path did not integrate")
        if (allocated(error)) return
        call check(error, s1 == s2 .and. all(d1 == d2) .and. all(g1 == g2), "the two paths differ")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_multi_matches_single

end module test_image_term
