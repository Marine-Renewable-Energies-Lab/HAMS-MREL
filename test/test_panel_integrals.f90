!> The analytical 1/r panel integrals, against an independent composite Gauss quadrature in this file, in every node order, on
!> the single-body and multi-body copies. The
!> triangle routines place node 2 in local coordinates from ALFA, the projection of side 1-3 on side 1-2 in units of side 1-2,
!> so the triangles here cover every range of ALFA. They use the mesh arrays of the solver, so this suite runs serially.
module test_panel_integrals

  use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
  use PanelMesh_mod, only : XYZ, XYZ_P, NCN, NCON, PNSZ, ISX, ISY, XYZ_GLOBAL_MULTI_COMB, XYZ_GLOBAL_MULTI_COMB_P, &
                            iXYZ_GLOBAL_MULTI_COMB_P, NCN_MULTI_COMB, NCON_MULTI_COMB, PNSZ_MULTI_COMB
  use Inerfs_mod, only : iXYZ_P
  use SingularIntgr, only : SGLINTBD_TRI, SGLINTBD_QUAD, SGLINTBD_TRI2, SGLINTBD_QUAD2
  use SingularIntgrMulti, only : SGLINTBD_MULTI_TRI, SGLINTBD_MULTI_QUAD, SGLINTBD_MULTI_TRI2, SGLINTBD_MULTI_QUAD2

  implicit none
  private
  public :: collect_panel_integrals

  real(8), parameter :: tol = 1d-8
  ! Nodes and weights of the 8-point Gauss-Legendre rule on (-1, 1), half of them, the others by symmetry
  real(8), parameter :: gx(4) = [0.1834346424956498d0, 0.5255324099163290d0, 0.7966664774136267d0, 0.9602898564975363d0]
  real(8), parameter :: gw(4) = [0.3626837833783620d0, 0.3137066458778873d0, 0.2223810344533745d0, 0.1012285362903763d0]
  integer, parameter :: default_cells = 48
  integer, parameter :: npaths = 6
  character(len=*), parameter :: path_name(npaths) = [character(len=28) :: &
    "single, hull field", "single, lid field", "single, given point", &
    "multi, hull field", "multi, lid field", "multi, given point"]

contains

  subroutine collect_panel_integrals(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("reference_quadrature", test_reference_quadrature), &
      new_unittest("triangle_obtuse_at_node_2", test_triangle_obtuse_at_node_2), &
      new_unittest("triangle_obtuse_at_node_2_extreme", test_triangle_obtuse_at_node_2_extreme), &
      new_unittest("triangle_obtuse_at_node_1", test_triangle_obtuse_at_node_1), &
      new_unittest("triangle_obtuse_at_node_3", test_triangle_obtuse_at_node_3), &
      new_unittest("triangle_right_angled", test_triangle_right_angled), &
      new_unittest("triangle_acute", test_triangle_acute), &
      new_unittest("triangle_thin_obtuse_at_waterline", test_triangle_thin_obtuse), &
      new_unittest("triangle_self_integral", test_triangle_self), &
      new_unittest("triangle_node_order", test_triangle_node_order), &
      new_unittest("quadrilaterals", test_quadrilaterals), &
      new_unittest("quadrilateral_self_integral", test_quadrilateral_self), &
      new_unittest("repeated_node_quadrilateral", test_repeated_node) &
      ]
  end subroutine collect_panel_integrals

  ! ---------------------------------------------------------------------------------------------- reference

  !> Integrals of 1/r and of d(1/r)/d(source) over the panel, r the distance to the point t: composite Gauss on the bilinear map
  !> of the square, or of the triangle collapsed to a square.
  subroutine reference(nodes, t, cells, s, d)
    real(8), intent(in) :: nodes(:, :), t(3)
    integer, intent(in) :: cells
    real(8), intent(out) :: s, d(3)
    real(8) :: a(4, 3), xa, xb, wa, wb, q(3), du(3), dv(3), jac, r, n(3)
    integer :: i, j, ia, ib

    a(1:size(nodes, 1), :) = nodes
    if (size(nodes, 1) == 3) a(4, :) = nodes(3, :)
    s = 0d0
    d = 0d0
    do i = 1, cells
      do ia = -4, 4
        if (ia == 0) cycle
        xa = (real(i - 1, 8) + 0.5d0*(1d0 + sign(gx(abs(ia)), real(ia, 8))))/real(cells, 8)
        wa = 0.5d0*gw(abs(ia))/real(cells, 8)
        do j = 1, cells
          do ib = -4, 4
            if (ib == 0) cycle
            xb = (real(j - 1, 8) + 0.5d0*(1d0 + sign(gx(abs(ib)), real(ib, 8))))/real(cells, 8)
            wb = 0.5d0*gw(abs(ib))/real(cells, 8)
            q = (1d0 - xa)*(1d0 - xb)*a(1, :) + xa*(1d0 - xb)*a(2, :) + xa*xb*a(3, :) + (1d0 - xa)*xb*a(4, :)
            du = (1d0 - xb)*(a(2, :) - a(1, :)) + xb*(a(3, :) - a(4, :))
            dv = (1d0 - xa)*(a(4, :) - a(1, :)) + xa*(a(3, :) - a(2, :))
            n = cross(du, dv)
            jac = norm2(n)
            r = norm2(q - t)
            s = s + wa*wb*jac/r
            d = d - wa*wb*jac*(q - t)/r**3
          end do
        end do
      end do
    end do
  end subroutine reference

  !> The integral of 1/r over the panel for a point t inside it, in its plane: a fan of triangles collapsed at t, where the
  !> Jacobian cancels the singularity.
  real(8) function self_reference(nodes, t) result(s)
    real(8), intent(in) :: nodes(:, :), t(3)
    real(8) :: fan(3, 3), si, di(3)
    integer :: k, n

    n = size(nodes, 1)
    s = 0d0
    do k = 1, n
      fan(1, :) = nodes(k, :)
      fan(2, :) = nodes(mod(k, n) + 1, :)
      fan(3, :) = t
      call reference(fan, t, 16, si, di)
      s = s + si
    end do
  end function self_reference

  ! ---------------------------------------------------------------------------------------------- routines under test

  !> One panel, and the field point p, in the arrays of both paths.
  subroutine set_panel(nodes, p)
    real(8), intent(in) :: nodes(:, :), p(3)
    integer :: k, n

    n = size(nodes, 1)
    ! One by one: the other serial suite allocates only some of these
    if (allocated(XYZ)) deallocate (XYZ)
    if (allocated(XYZ_P)) deallocate (XYZ_P)
    if (allocated(NCN)) deallocate (NCN)
    if (allocated(NCON)) deallocate (NCON)
    if (allocated(PNSZ)) deallocate (PNSZ)
    if (allocated(iXYZ_P)) deallocate (iXYZ_P)
    if (allocated(XYZ_GLOBAL_MULTI_COMB)) deallocate (XYZ_GLOBAL_MULTI_COMB)
    if (allocated(XYZ_GLOBAL_MULTI_COMB_P)) deallocate (XYZ_GLOBAL_MULTI_COMB_P)
    if (allocated(iXYZ_GLOBAL_MULTI_COMB_P)) deallocate (iXYZ_GLOBAL_MULTI_COMB_P)
    if (allocated(NCN_MULTI_COMB)) deallocate (NCN_MULTI_COMB)
    if (allocated(NCON_MULTI_COMB)) deallocate (NCON_MULTI_COMB)
    if (allocated(PNSZ_MULTI_COMB)) deallocate (PNSZ_MULTI_COMB)
    allocate (XYZ(n, 3), XYZ_P(1, 3), NCN(1), NCON(1, 4), PNSZ(1), iXYZ_P(1, 3))
    allocate (XYZ_GLOBAL_MULTI_COMB(n, 3), XYZ_GLOBAL_MULTI_COMB_P(1, 3), iXYZ_GLOBAL_MULTI_COMB_P(1, 3), &
              NCN_MULTI_COMB(1), NCON_MULTI_COMB(1, 4), PNSZ_MULTI_COMB(1))
    XYZ = nodes
    NCN(1) = n
    NCON = 0
    do k = 1, n
      NCON(1, k) = k
    end do
    PNSZ(1) = longest_side(nodes)
    XYZ_P(1, :) = p
    iXYZ_P(1, :) = p
    XYZ_GLOBAL_MULTI_COMB = XYZ
    NCN_MULTI_COMB = NCN
    NCON_MULTI_COMB = NCON
    PNSZ_MULTI_COMB = PNSZ
    XYZ_GLOBAL_MULTI_COMB_P(1, :) = p
    iXYZ_GLOBAL_MULTI_COMB_P(1, :) = p
    ISX = 0
    ISY = 0
  end subroutine set_panel

  !> Run one path on the panel of set_panel.
  subroutine evaluate(path, p, s, d)
    integer, intent(in) :: path
    real(8), intent(in) :: p(3)
    real(8), intent(out) :: s, d(3)
    real(8) :: rk(4)
    logical :: tri

    tri = NCN(1) == 3
    select case (path)
    case (1, 2)
      if (tri) then
        call SGLINTBD_TRI(1, 1, 1, s, d, 2*path - 1)
      else
        call SGLINTBD_QUAD(1, 1, 1, s, d, 2*path - 1)
      end if
    case (3)
      if (tri) then
        call SGLINTBD_TRI2(p, 1, rk, 1)
      else
        call SGLINTBD_QUAD2(p, 1, rk, 1)
      end if
      s = rk(1)
      d = rk(2:4)
    case (4, 5)
      if (tri) then
        call SGLINTBD_MULTI_TRI(1, 1, 1, s, d, 2*(path - 3) - 1)
      else
        call SGLINTBD_MULTI_QUAD(1, 1, 1, s, d, 2*(path - 3) - 1)
      end if
    case (6)
      if (tri) then
        call SGLINTBD_MULTI_TRI2(p, 1, rk, 1)
      else
        call SGLINTBD_MULTI_QUAD2(p, 1, rk, 1)
      end if
      s = rk(1)
      d = rk(2:4)
    end select
  end subroutine evaluate

  ! ---------------------------------------------------------------------------------------------- geometry

  function cross(a, b) result(c)
    real(8), intent(in) :: a(3), b(3)
    real(8) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross

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

  !> ALFA as the triangle routines compute it, from the nodes as written.
  real(8) function alfa(uv)
    real(8), intent(in) :: uv(3, 2)
    alfa = dot_product(uv(2, :) - uv(1, :), uv(3, :) - uv(1, :))/sum((uv(2, :) - uv(1, :))**2)
  end function alfa

  function tri_uv(u1, v1, u2, v2, u3, v3) result(uv)
    real(8), intent(in) :: u1, v1, u2, v2, u3, v3
    real(8) :: uv(3, 2)
    uv(1, :) = [u1, v1]
    uv(2, :) = [u2, v2]
    uv(3, :) = [u3, v3]
  end function tri_uv

  !> The panel in space, from coordinates (u, v) with v >= 0. kind 1: a wall at x = 0, with v down from z = -0.1; kind 2: tilted;
  !> kind 3: horizontal; kind 4: a wall with v down from z = 0.
  subroutine embed(uv, kind, nodes, e1, en)
    real(8), intent(in) :: uv(:, :)
    integer, intent(in) :: kind
    real(8), intent(out) :: nodes(size(uv, 1), 3), e1(3), en(3)
    real(8) :: o(3), e2(3)
    integer :: k

    select case (kind)
    case (1)
      o = [0d0, 0d0, -0.1d0]; e1 = [0d0, 1d0, 0d0]; e2 = [0d0, 0d0, -1d0]
    case (2)
      o = [1d0, 2d0, -0.3d0]; e1 = [0.8d0, 0.6d0, 0d0]; e2 = [0.3d0, -0.4d0, -1d0]/sqrt(1.25d0)
    case (3)
      o = [0d0, 0d0, -0.5d0]; e1 = [1d0, 0d0, 0d0]; e2 = [0d0, 1d0, 0d0]
    case default
      o = 0d0; e1 = [0d0, 1d0, 0d0]; e2 = [0d0, 0d0, -1d0]
    end select
    do k = 1, size(uv, 1)
      nodes(k, :) = o + uv(k, 1)*e1 + uv(k, 2)*e2
    end do
    en = cross(e1, e2)
  end subroutine embed

  !> Node order number ord: the cyclic orders, then the same reversed.
  function ordered(nodes, ord) result(out)
    real(8), intent(in) :: nodes(:, :)
    integer, intent(in) :: ord
    real(8) :: out(size(nodes, 1), size(nodes, 2))
    integer :: k, n, m

    n = size(nodes, 1)
    do k = 1, n
      m = mod(k - 1 + mod(ord - 1, n), n) + 1
      if (ord > n) m = n + 1 - m
      out(k, :) = nodes(m, :)
    end do
  end function ordered

  ! ---------------------------------------------------------------------------------------------- checks

  subroutine compare(s, d, sref, dref, what, error)
    real(8), intent(in) :: s, d(3), sref, dref(3)
    character(len=*), intent(in) :: what
    type(error_type), allocatable, intent(out) :: error
    character(len=300) :: msg

    if (abs(s - sref) > tol*max(abs(sref), 1d-3) .or. maxval(abs(d - dref)) > tol*max(maxval(abs(dref)), 1d-3) &
        .or. .not. (s == s) .or. any(d /= d)) then
      write (msg, '(a, ": G ", es12.5, " against ", es12.5, ", gradient ", 3es12.4, " against ", 3es12.4)') &
        what, s, sref, d, dref
      call test_failed(error, trim(msg))
    end if
  end subroutine compare

  !> Every path and every node order, at field point p.
  subroutine check_paths(nodes, p, cells, name, error)
    real(8), intent(in) :: nodes(:, :), p(3)
    integer, intent(in) :: cells
    character(len=*), intent(in) :: name
    type(error_type), allocatable, intent(out) :: error
    real(8) :: sref, dref(3), s, d(3)
    integer :: ord, path
    character(len=2) :: ostr

    call reference(nodes, p, cells, sref, dref)
    do ord = 1, 2*size(nodes, 1)
      call set_panel(ordered(nodes, ord), p)
      write (ostr, '(i2)') ord
      do path = 1, npaths
        call evaluate(path, p, s, d)
        call compare(s, d, sref, dref, name//", order"//ostr//", "//trim(path_name(path)), error)
        if (allocated(error)) return
      end do
    end do
  end subroutine check_paths

  !> On a wall and on a tilted and a horizontal panel: off the panel and at the next panel along. check_self covers the centroid.
  subroutine check_shape(uv, name, error)
    real(8), intent(in) :: uv(:, :)
    character(len=*), intent(in) :: name
    type(error_type), allocatable, intent(out) :: error
    real(8) :: nodes(size(uv, 1), 3), e1(3), en(3), c(3), side
    integer :: kind
    character(len=1) :: kstr

    do kind = 1, 3
      call embed(uv, kind, nodes, e1, en)
      c = centroid(nodes)
      side = longest_side(nodes)
      write (kstr, '(i1)') kind
      call check_paths(nodes, c + 0.15d0*side*en, default_cells, name//", kind "//kstr//", off the panel", error)
      if (allocated(error)) return
      call check_paths(nodes, c + 1.4d0*side*e1, default_cells, name//", kind "//kstr//", next along", error)
      if (allocated(error)) return
    end do
  end subroutine check_shape

  !> The direct term of a panel at its own centroid, the case of the diagonal of the matrix: G against the fan, and the two
  !> gradient components in the plane against central differences of G. The component normal to the panel is left to the caller
  !> of these routines, which adds the solid angle, so it is zero here.
  subroutine check_self(uv, name, error)
    real(8), intent(in) :: uv(:, :)
    character(len=*), intent(in) :: name
    type(error_type), allocatable, intent(out) :: error
    real(8) :: nodes(size(uv, 1), 3), e1(3), en(3), c(3), s, d(3), sref, t(3, 2), dt(2), h, gap, a(3), b(3), en_u(3)
    integer :: kind, ord, path, m, k, n
    character(len=40) :: tag

    n = size(uv, 1)
    do kind = 1, 3
      call embed(uv, kind, nodes, e1, en)
      en_u = en/norm2(en)
      c = centroid(nodes)
      gap = huge(1d0)
      do k = 1, n
        a = nodes(k, :) - c
        b = nodes(mod(k, n) + 1, :) - c
        if (norm2(b - a) > 1d-12) gap = min(gap, norm2(cross(a, b))/norm2(b - a))
      end do
      h = 0.02d0*gap
      t(:, 1) = e1
      t(:, 2) = cross(en_u, e1)
      sref = self_reference(nodes, c)
      do m = 1, 2
        dt(m) = -(-self_reference(nodes, c + 2d0*h*t(:, m)) + 8d0*self_reference(nodes, c + h*t(:, m)) &
                  - 8d0*self_reference(nodes, c - h*t(:, m)) + self_reference(nodes, c - 2d0*h*t(:, m)))/(12d0*h)
      end do
      do ord = 1, 2*n
        call set_panel(ordered(nodes, ord), c)
        do path = 1, npaths
          call evaluate(path, c, s, d)
          write (tag, '(", kind ", i1, ", order ", i1, ", path ", i1)') kind, ord, path
          call check(error, abs(s - sref) < tol*abs(sref), name//trim(tag)//": G at its own centroid")
          if (allocated(error)) return
          call check(error, abs(dot_product(d, t(:, 1)) - dt(1)) < 1d-6*max(maxval(abs(dt)), abs(sref)) .and. &
                     abs(dot_product(d, t(:, 2)) - dt(2)) < 1d-6*max(maxval(abs(dt)), abs(sref)), &
                     name//trim(tag)//": gradient in the plane at its own centroid")
          if (allocated(error)) return
          call check(error, abs(dot_product(d, en_u)) < 1d-10*abs(sref), name//trim(tag)//": normal gradient")
          if (allocated(error)) return
        end do
      end do
    end do
  end subroutine check_self

  ! ---------------------------------------------------------------------------------------------- tests

  !> The reference itself: refining it changes nothing, it does not depend on the node order, and far off it tends to A/r.
  subroutine test_reference_quadrature(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: uv(3, 2), nodes(3, 3), e1(3), en(3), s1, d1(3), s2, d2(3), p(3)
    integer :: ord

    uv = tri_uv(0d0, 0d0, 1d0, 0d0, 1.5d0, 0.6d0)
    call embed(uv, 2, nodes, e1, en)
    p = centroid(nodes) + 0.15d0*en
    call reference(nodes, p, default_cells, s1, d1)
    call reference(nodes, p, 2*default_cells, s2, d2)
    call check(error, abs(s1 - s2) < 1d-11 .and. maxval(abs(d1 - d2)) < 1d-10, "refining the reference")
    if (allocated(error)) return
    do ord = 2, 6
      call reference(ordered(nodes, ord), p, default_cells, s2, d2)
      call check(error, abs(s1 - s2) < 1d-11 .and. maxval(abs(d1 - d2)) < 1d-10, "node order in the reference")
      if (allocated(error)) return
    end do
    ! Area 0.3 of this triangle, from its base 1 and height 0.6
    p = centroid(nodes) + 1000d0*en
    call reference(nodes, p, default_cells, s1, d1)
    call check(error, abs(s1 - 0.3d0/1000d0) < 1d-8, "far off, A/r")
  end subroutine test_reference_quadrature

  !> ALFA > 1, the case the finding is about.
  subroutine test_triangle_obtuse_at_node_2(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: uv(3, 2)

    uv = tri_uv(0d0, 0d0, 1d0, 0d0, 1.5d0, 0.6d0)
    call check(error, alfa(uv) > 1d0, "the shape has ALFA > 1")
    if (allocated(error)) return
    call check_shape(uv, "ALFA 1.5", error)
  end subroutine test_triangle_obtuse_at_node_2

  !> The largest ALFA of the repo meshes, 5.37, in example 1.
  subroutine test_triangle_obtuse_at_node_2_extreme(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: uv(3, 2)

    uv = tri_uv(0d0, 0d0, 0.2d0, 0d0, 1.074d0, 0.3d0)
    call check(error, alfa(uv) > 5d0, "the shape has ALFA > 5")
    if (allocated(error)) return
    call check_shape(uv, "ALFA 5.37", error)
  end subroutine test_triangle_obtuse_at_node_2_extreme

  !> ALFA < 0.
  subroutine test_triangle_obtuse_at_node_1(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: uv(3, 2)

    uv = tri_uv(0d0, 0d0, 1d0, 0d0, -0.5d0, 0.6d0)
    call check(error, alfa(uv) < 0d0, "the shape has ALFA < 0")
    if (allocated(error)) return
    call check_shape(uv, "ALFA -0.5", error)
  end subroutine test_triangle_obtuse_at_node_1

  !> 0 < ALFA < 1 although the triangle is obtuse. Its other cyclic orders have ALFA > 1 and ALFA < 0.
  subroutine test_triangle_obtuse_at_node_3(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: uv(3, 2)

    uv = tri_uv(0d0, 0d0, 1d0, 0d0, 0.5d0, 0.15d0)
    call check_shape(uv, "obtuse at node 3", error)
  end subroutine test_triangle_obtuse_at_node_3

  !> The right angle at node 2 is ALFA = 1, at node 1 ALFA = 0.
  subroutine test_triangle_right_angled(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: uv(3, 2)

    uv = tri_uv(0d0, 0d0, 1d0, 0d0, 1d0, 0.7d0)
    call check(error, abs(alfa(uv) - 1d0) < 1d-14, "the shape has ALFA = 1")
    if (allocated(error)) return
    call check_shape(uv, "ALFA 1", error)
    if (allocated(error)) return
    uv = tri_uv(0d0, 0d0, 1d0, 0d0, 0d0, 0.7d0)
    call check(error, abs(alfa(uv)) < 1d-14, "the shape has ALFA = 0")
    if (allocated(error)) return
    call check_shape(uv, "ALFA 0", error)
  end subroutine test_triangle_right_angled

  subroutine test_triangle_acute(error)
    type(error_type), allocatable, intent(out) :: error
    call check_shape(tri_uv(0d0, 0d0, 1d0, 0d0, 0.4d0, 0.8d0), "acute", error)
  end subroutine test_triangle_acute

  !> A thin triangle at the waterline with ALFA > 1, on a wall. It needs more cells.
  subroutine test_triangle_thin_obtuse(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: uv(3, 2), nodes(3, 3), e1(3), en(3), c(3), side

    uv = tri_uv(0d0, 0d0, 0.4d0, 0d0, 0.9d0, 0.03d0)
    call check(error, alfa(uv) > 1d0, "the shape has ALFA > 1")
    if (allocated(error)) return
    call embed(uv, 4, nodes, e1, en)
    c = centroid(nodes)
    side = longest_side(nodes)
    call check_paths(nodes, c + 0.15d0*side*en, 160, "thin, off the panel", error)
    if (allocated(error)) return
    call check_paths(nodes, c + 1.4d0*side*e1, 160, "thin, next along", error)
  end subroutine test_triangle_thin_obtuse

  subroutine test_triangle_self(error)
    type(error_type), allocatable, intent(out) :: error

    call check_self(tri_uv(0d0, 0d0, 1d0, 0d0, 1.5d0, 0.6d0), "ALFA 1.5", error)
    if (allocated(error)) return
    call check_self(tri_uv(0d0, 0d0, 1d0, 0d0, -0.5d0, 0.6d0), "ALFA -0.5", error)
    if (allocated(error)) return
    call check_self(tri_uv(0d0, 0d0, 1d0, 0d0, 1d0, 0.7d0), "ALFA 1", error)
    if (allocated(error)) return
    call check_self(tri_uv(0d0, 0d0, 1d0, 0d0, 0.4d0, 0.8d0), "acute", error)
  end subroutine test_triangle_self

  !> The same panel gives the same numbers in any node order, to rounding, without the reference.
  subroutine test_triangle_node_order(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: uvs(3, 2, 4), nodes(3, 3), e1(3), en(3), p(3), s0, d0(3), s, d(3)
    integer :: shape, kind, ord, path

    uvs(:, :, 1) = tri_uv(0d0, 0d0, 1d0, 0d0, 1.5d0, 0.6d0)
    uvs(:, :, 2) = tri_uv(0d0, 0d0, 0.2d0, 0d0, 1.074d0, 0.3d0)
    uvs(:, :, 3) = tri_uv(0d0, 0d0, 1d0, 0d0, -0.5d0, 0.6d0)
    uvs(:, :, 4) = tri_uv(0d0, 0d0, 1d0, 0d0, 0.4d0, 0.8d0)
    do shape = 1, 4
      do kind = 1, 3
        call embed(uvs(:, :, shape), kind, nodes, e1, en)
        p = centroid(nodes) + 0.15d0*en + 0.1d0*e1
        do path = 1, npaths
          do ord = 1, 6
            call set_panel(ordered(nodes, ord), p)
            call evaluate(path, p, s, d)
            if (ord == 1) then
              s0 = s
              d0 = d
            end if
            call check(error, abs(s - s0) < 1d-11*abs(s0) .and. maxval(abs(d - d0)) < 1d-10*maxval(abs(d0)), &
                       "the result depends on the node order")
            if (allocated(error)) return
          end do
        end do
      end do
    end do
  end subroutine test_triangle_node_order

  function quad_uv(u1, v1, u2, v2, u3, v3, u4, v4) result(uv)
    real(8), intent(in) :: u1, v1, u2, v2, u3, v3, u4, v4
    real(8) :: uv(4, 2)
    uv(1, :) = [u1, v1]
    uv(2, :) = [u2, v2]
    uv(3, :) = [u3, v3]
    uv(4, :) = [u4, v4]
  end function quad_uv

  !> The quadrilateral routines have no ALFA: they take local coordinates from dot products. These check that this holds.
  subroutine test_quadrilaterals(error)
    type(error_type), allocatable, intent(out) :: error

    call check_shape(quad_uv(0d0, 0d0, 1d0, 0d0, 1d0, 0.7d0, 0d0, 0.7d0), "rectangle", error)
    if (allocated(error)) return
    call check_shape(quad_uv(0d0, 0d0, 1d0, 0d0, 0.7d0, 0.6d0, 0.3d0, 0.6d0), "trapezoid", error)
    if (allocated(error)) return
    ! A parallelogram, with obtuse angles like the triangles above, and a quadrilateral with no parallel sides
    call check_shape(quad_uv(0d0, 0d0, 1d0, 0d0, 1.8d0, 0.5d0, 0.8d0, 0.5d0), "parallelogram", error)
    if (allocated(error)) return
    call check_shape(quad_uv(0d0, 0d0, 1d0, 0.1d0, 1.3d0, 0.8d0, 0.2d0, 0.5d0), "general", error)
  end subroutine test_quadrilaterals

  subroutine test_quadrilateral_self(error)
    type(error_type), allocatable, intent(out) :: error

    call check_self(quad_uv(0d0, 0d0, 1d0, 0d0, 1d0, 0.7d0, 0d0, 0.7d0), "rectangle", error)
    if (allocated(error)) return
    call check_self(quad_uv(0d0, 0d0, 1d0, 0.1d0, 1.3d0, 0.8d0, 0.2d0, 0.5d0), "general", error)
  end subroutine test_quadrilateral_self

  !> A triangle with ALFA > 1 written as a quadrilateral, with a repeated node: the quadrilateral routine is right for it. Nodes 1
  !> and 2 at the same point are left out, as the quadrilateral routine gives NaN for that form.
  subroutine test_repeated_node(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: tri(3, 3), e1(3), en(3), nodes(4, 3), p(3), sd, dd(3), s, d(3)
    integer :: kind, form, path

    do kind = 1, 3
      call embed(tri_uv(0d0, 0d0, 1d0, 0d0, 1.5d0, 0.6d0), kind, tri, e1, en)
      p = centroid(tri) + 0.15d0*en + 0.1d0*e1
      call reference(tri, p, default_cells, sd, dd)
      do form = 1, 3
        select case (form)
        case (1)  ! a b c c
          nodes(1:3, :) = tri
          nodes(4, :) = tri(3, :)
        case (2)  ! a b b c
          nodes(1:2, :) = tri(1:2, :)
          nodes(3:4, :) = tri(2:3, :)
        case (3)  ! a a b c is the one left out, so c a b c
          nodes(1, :) = tri(3, :)
          nodes(2:4, :) = tri
        end select
        call set_panel(nodes, p)
        do path = 1, npaths
          call evaluate(path, p, s, d)
          call compare(s, d, sd, dd, "repeated node", error)
          if (allocated(error)) return
        end do
      end do
    end do
  end subroutine test_repeated_node

end module test_panel_integrals
