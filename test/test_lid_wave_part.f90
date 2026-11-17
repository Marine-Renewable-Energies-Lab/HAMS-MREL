!> The wave part of the Green function integrated over a hull panel, seen from a lid point on z = 0 (IRR = 3 rows).
!> A lid point close to a waterline panel sees a strongly varying wave part, which a single value at the panel centroid times the
!> area gets wrong. WAVE_PART_PANEL integrates it; this suite checks it against an independent composite Gauss quadrature written
!> here. It uses the mesh arrays of the solver, which are module globals, so this suite runs serially.
module test_lid_wave_part

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use PanelMesh_mod, only : XYZ, NCN, NCON, PNSZ, XYZ_GLOBAL_MULTI_COMB, NCN_MULTI_COMB, NCON_MULTI_COMB, PNSZ_MULTI_COMB
  use PanelMesh_mod, only : XYZ_P, DS, NELEM, NSYS, ISX, ISY, XYZ_GLOBAL_MULTI_COMB_P, DS_MULTI_COMB, NELEM_TOTAL, iNELEM_TOTAL
  use Inerfs_mod, only : iXYZ_P, iNELEM
  use PanelMesh_mod, only : iXYZ_GLOBAL_MULTI_COMB_P
  use Potentials_mod, only : RKBN, CGRN, PKBN, DGRN, RKBN_MULTI_COMB, CGRN_MULTI_COMB, PKBN_MULTI_COMB, DGRN_MULTI_COMB
  use WaveDyn_mod, only : V_SOLVER => V, H_SOLVER => H, WVN_SOLVER => WVN
  use HAMS_mod, only : NTHREAD, GRN_KIND
  use SingularIntgr, only : WAVE_PART_PANEL, IMAGE_TERM_EXACT
  use SingularIntgrMulti, only : WAVE_PART_PANEL_MULTI
  use CalGreenFunc, only : CALGREEN_IRR
  use CalGreenFuncMulti, only : CALGREEN_MULTI_IRR
  use FinGreen3D_Open, only : FINGREEN3D
  use INFG3D_Open, only : INFGREEN3D

  implicit none
  private
  public :: collect_lid_wave_part

  integer, parameter :: NK = 200
  real(8), parameter :: GRAV = 9.80665d0
  real(8), parameter :: tol = 1d-3
  ! Nodes and weights of the 8-point Gauss-Legendre rule on (-1, 1), half of them, the others by symmetry
  real(8), parameter :: gx(4) = [0.1834346424956498d0, 0.5255324099163290d0, 0.7966664774136267d0, 0.9602898564975363d0]
  real(8), parameter :: gw(4) = [0.3626837833783620d0, 0.3137066458778873d0, 0.2223810344533745d0, 0.1012285362903763d0]
  integer, parameter :: cells = 48
  ! A waterline panel of a column, 0.44 m wide and 0.68 m high, as on the buoy that showed the problem; T = 2.5 s.
  real(8), parameter :: pw = 0.44d0, ph = 0.68d0, period = 2.5d0

contains

  subroutine collect_lid_wave_part(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("near_wall_deep", test_near_wall_deep), &
      new_unittest("near_wall_finite_depth", test_near_wall_finite_depth), &
      new_unittest("triangle_finite_depth", test_triangle_finite_depth), &
      new_unittest("limits_have_no_wave_part", test_limits_have_no_wave_part), &
      new_unittest("multi_matches_single", test_multi_matches_single), &
      new_unittest("lid_row_holds_the_integral", test_lid_row_holds_the_integral), &
      new_unittest("lid_row_holds_the_integral_multi", test_lid_row_holds_the_integral_multi), &
      new_unittest("hull_row_and_far_lid_point_keep_one_point", test_hull_row_and_far_lid_point_keep_one_point) &
      ]
  end subroutine collect_lid_wave_part

  ! ---------------------------------------------------------------------------------------------- helpers

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
    PNSZ(1) = maxval([(norm2(nodes(k, :) - nodes(mod(k, size(nodes, 1)) + 1, :)), k=1, size(nodes, 1))])
  end subroutine set_mesh

  function wall_panel() result(nodes)
    real(8) :: nodes(4, 3)
    nodes(1, :) = [0d0, -0.5d0*pw, 0d0]
    nodes(2, :) = [0d0, 0.5d0*pw, 0d0]
    nodes(3, :) = [0d0, 0.5d0*pw, -ph]
    nodes(4, :) = [0d0, -0.5d0*pw, -ph]
  end function wall_panel

  !> The wave part of G and of its source gradient at source q, field point p: the Green function with FLAG = 1 (no 1/r) less
  !> the free-surface image 1/r1, which IMAGE_TERM_EXACT integrates separately.
  subroutine wave_part(q, p, v, wvn, h, g)
    real(8), intent(in) :: q(3), p(3), v, wvn(NK), h
    complex(8), intent(out) :: g(4)
    real(8) :: pq(3), r1

    if (h < 0d0) then
      call INFGREEN3D(q(1), p(1), q(2), p(2), q(3), p(3), v, g, 1)
    else
      call FINGREEN3D(q(1), p(1), q(2), p(2), q(3), p(3), v, wvn, NK, h, g, 1)
    end if
    pq = [q(1) - p(1), q(2) - p(2), q(3) + p(3)]
    r1 = norm2(pq)
    g(1) = g(1) - 1d0/r1
    g(2:4) = g(2:4) + pq/r1**3
  end subroutine wave_part

  !> Reference: composite 8-point Gauss on the bilinear map of the panel (a triangle as a quadrilateral with a repeated node).
  subroutine reference(nodes, p, v, wvn, h, s)
    real(8), intent(in) :: nodes(:, :), p(3), v, wvn(NK), h
    complex(8), intent(out) :: s(4)
    real(8) :: a(4, 3), xa, xb, wa, wb, q(3), du(3), dv(3), n(3)
    complex(8) :: g(4)
    integer :: i, j, ia, ib

    a(1:size(nodes, 1), :) = nodes
    if (size(nodes, 1) == 3) a(4, :) = nodes(3, :)
    s = (0d0, 0d0)
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
            call wave_part(q, p, v, wvn, h, g)
            s = s + wa*wb*norm2(n)*g
          end do
        end do
      end do
    end do
  end subroutine reference

  real(8) function area(nodes)
    real(8), intent(in) :: nodes(:, :)
    real(8) :: c(3), n(3)
    integer :: k, m

    c = sum(nodes, dim=1)/real(size(nodes, 1), 8)
    n = 0d0
    m = size(nodes, 1)
    do k = 1, m
      n = n + cross(nodes(k, :) - c, nodes(mod(k, m) + 1, :) - c)
    end do
    area = 0.5d0*norm2(n)
  end function area

  function cross(a, b) result(c)
    real(8), intent(in) :: a(3), b(3)
    real(8) :: c(3)
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross

  !> Integral against the reference, and the one-point value (centroid value times area) must miss it: otherwise the case does
  !> not exercise the integration.
  subroutine check_case(nodes, p, h, name, error)
    real(8), intent(in) :: nodes(:, :), p(3), h
    character(len=*), intent(in) :: name
    type(error_type), allocatable, intent(out) :: error
    real(8) :: v, wvn(NK), c(3)
    complex(8) :: ref(4), got(4), one(4)
    real(8) :: e0, eg, o0, og
    external :: DISPERSION

    v = (2d0*acos(-1d0)/period)**2/GRAV
    wvn = 0d0
    if (h > 0d0) call DISPERSION(wvn, NK, 2d0*acos(-1d0)/period, h)
    call set_mesh(nodes)
    call reference(nodes, p, v, wvn, h, ref)
    call WAVE_PART_PANEL(p, 1, v, wvn, NK, h, got)
    c = sum(nodes, dim=1)/real(size(nodes, 1), 8)
    call wave_part(c, p, v, wvn, h, one)
    one = one*area(nodes)
    e0 = abs(got(1) - ref(1))/abs(ref(1))
    eg = norm2(abs(got(2:4) - ref(2:4)))/norm2(abs(ref(2:4)))
    o0 = abs(one(1) - ref(1))/abs(ref(1))
    og = norm2(abs(one(2:4) - ref(2:4)))/norm2(abs(ref(2:4)))
    write (*, '(4x, a, ": integral err G ", es9.2, " grad ", es9.2, "; one point err G ", es9.2, " grad ", es9.2)') &
      name, e0, eg, o0, og
    call check(error, og > 10d0*tol, name//": the one-point value is already close, the case tests nothing")
    if (allocated(error)) return
    call check(error, e0 < tol, name//": integral of the wave part of G")
    if (allocated(error)) return
    call check(error, eg < tol, name//": integral of the wave part of the gradient")
  end subroutine check_case

  ! ---------------------------------------------------------------------------------------------- tests

  !> Lid point 0.034 m from the wall (d/h = 0.05), level with the top edge, deep water.
  subroutine test_near_wall_deep(error)
    type(error_type), allocatable, intent(out) :: error
    call check_case(wall_panel(), [0.05d0*ph, 0.1d0, 0d0], -1d0, "near wall, deep water", error)
  end subroutine test_near_wall_deep

  !> The same in 30 m of water (nu*h = 19). Not deeper: above nu*h = 50 FINGREEN3D_R uses PIDCOCK down to R/h = 5e-4, and
  !> its gradient there has a kink of a few percent at that switch, which no quadrature of this order follows.
  subroutine test_near_wall_finite_depth(error)
    type(error_type), allocatable, intent(out) :: error
    call check_case(wall_panel(), [0.05d0*ph, 0.1d0, 0d0], 30d0, "near wall, h = 30 m", error)
  end subroutine test_near_wall_finite_depth

  !> A triangle with its top edge on the waterline, as in a triangulated hull.
  subroutine test_triangle_finite_depth(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: tri(3, 3)
    tri(1, :) = [0d0, -0.5d0*pw, 0d0]
    tri(2, :) = [0d0, 0.5d0*pw, 0d0]
    tri(3, :) = [0d0, 0.5d0*pw, -ph]
    call check_case(tri, [0.03d0, 0.15d0, 0d0], 30d0, "triangle, h = 30 m", error)
  end subroutine test_triangle_finite_depth

  !> At the zero- and infinite-frequency limits G has no wave part: the integral is zero, so those limits stay as they were.
  subroutine test_limits_have_no_wave_part(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: wvn(NK), v(2), h(2)
    complex(8) :: got(4)
    integer :: i, j

    call set_mesh(wall_panel())
    wvn = 1d0
    v = [0d0, -1d0]
    h = [-1d0, 100d0]
    do i = 1, 2
      do j = 1, 2
        got = (1d0, 1d0)
        call WAVE_PART_PANEL([0.03d0, 0.1d0, 0d0], 1, v(i), wvn, NK, h(j), got)
        call check(error, all(got == (0d0, 0d0)), "the wave part at a limit is not zero")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_limits_have_no_wave_part

  !> The multi-body path has its own copy of the routine and of the mesh arrays; it must give the same numbers.
  subroutine test_multi_matches_single(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: wvn(NK), v, p(3)
    complex(8) :: g1(4), g2(4)
    external :: DISPERSION

    call set_mesh(wall_panel())
    if (allocated(XYZ_GLOBAL_MULTI_COMB)) deallocate (XYZ_GLOBAL_MULTI_COMB, NCN_MULTI_COMB, NCON_MULTI_COMB, PNSZ_MULTI_COMB)
    allocate (XYZ_GLOBAL_MULTI_COMB(4, 3), NCN_MULTI_COMB(1), NCON_MULTI_COMB(1, 4), PNSZ_MULTI_COMB(1))
    XYZ_GLOBAL_MULTI_COMB = XYZ
    NCN_MULTI_COMB = NCN
    NCON_MULTI_COMB = NCON
    PNSZ_MULTI_COMB = PNSZ
    v = (2d0*acos(-1d0)/period)**2/GRAV
    call DISPERSION(wvn, NK, 2d0*acos(-1d0)/period, 100d0)
    p = [0.03d0, 0.1d0, 0d0]
    call WAVE_PART_PANEL(p, 1, v, wvn, NK, 100d0, g1)
    call WAVE_PART_PANEL_MULTI(p, 1, v, wvn, NK, 100d0, g2)
    call check(error, any(g1 /= (0d0, 0d0)), "the single-body path gave zero")
    if (allocated(error)) return
    call check(error, all(g1 == g2), "the two paths differ")
  end subroutine test_multi_matches_single

  ! ---------------------------------------------------------------------------------------------- the solver's Green tables

  !> One wall panel as the hull and one lid point, in the arrays CALGREEN_IRR reads, deep water at T = 2.5 s.
  subroutine set_solver(p)
    real(8), intent(in) :: p(3)

    call set_mesh(wall_panel())
    if (allocated(XYZ_P)) deallocate (XYZ_P, DS, iXYZ_P)
    allocate (XYZ_P(1, 3), DS(1), iXYZ_P(1, 3))
    XYZ_P(1, :) = sum(XYZ(1:4, :), dim=1)/4d0
    DS(1) = pw*ph
    iXYZ_P(1, :) = p
    NELEM = 1
    iNELEM = 1
    NSYS = 1
    ISX = 0
    ISY = 0
    NTHREAD = 1
    V_SOLVER = (2d0*acos(-1d0)/period)**2/GRAV
    H_SOLVER = -1d0
    WVN_SOLVER = 0d0
    if (allocated(RKBN)) deallocate (RKBN, CGRN, PKBN, DGRN)
    allocate (RKBN(1, 1, 1, 4), CGRN(1, 1, 1, 4), PKBN(1, 1, 1, 4), DGRN(1, 1, 1, 4))
  end subroutine set_solver

  !> The one-point value the tables held before: the Green function at the centroid, its image moved out by IMAGE_TERM_EXACT.
  function one_point(p, irr) result(g)
    real(8), intent(in) :: p(3)
    integer, intent(in) :: irr
    complex(GRN_KIND) :: g(4)
    complex(8) :: grn(4)
    real(8) :: q(3), sij, dij(3)

    q = XYZ_P(1, :)
    call INFGREEN3D(q(1), p(1), q(2), p(2), q(3), p(3), V_SOLVER, grn, 1)
    sij = 0d0
    dij = 0d0
    call IMAGE_TERM_EXACT(p, q, V_SOLVER, 1, irr, grn, sij, dij)
    g = grn
  end function one_point

  !> A lid point 0.034 m from the wall: its table entry is the panel integral over the area, not the one-point value.
  subroutine test_lid_row_holds_the_integral(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: p(3), wvn(NK)
    complex(8) :: gint(4)

    p = [0.05d0*ph, 0.1d0, 0d0]
    call set_solver(p)
    call CALGREEN_IRR
    wvn = 0d0
    call WAVE_PART_PANEL(p, 1, V_SOLVER, wvn, NK, -1d0, gint)
    call check(error, maxval(abs(DGRN(1, 1, 1, :) - gint/DS(1))) < 1d-6*maxval(abs(gint/DS(1))), &
               "the lid row does not hold the panel integral")
  end subroutine test_lid_row_holds_the_integral

  !> The same through the multi-body tables.
  subroutine test_lid_row_holds_the_integral_multi(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: p(3), wvn(NK)
    complex(8) :: gint(4)

    p = [0.05d0*ph, 0.1d0, 0d0]
    call set_solver(p)
    if (allocated(XYZ_GLOBAL_MULTI_COMB)) deallocate (XYZ_GLOBAL_MULTI_COMB, NCN_MULTI_COMB, NCON_MULTI_COMB, PNSZ_MULTI_COMB)
    if (allocated(XYZ_GLOBAL_MULTI_COMB_P)) deallocate (XYZ_GLOBAL_MULTI_COMB_P, iXYZ_GLOBAL_MULTI_COMB_P, DS_MULTI_COMB)
    if (allocated(NELEM_TOTAL)) deallocate (NELEM_TOTAL, iNELEM_TOTAL)
    if (allocated(RKBN_MULTI_COMB)) deallocate (RKBN_MULTI_COMB, CGRN_MULTI_COMB, PKBN_MULTI_COMB, DGRN_MULTI_COMB)
    allocate (XYZ_GLOBAL_MULTI_COMB(4, 3), NCN_MULTI_COMB(1), NCON_MULTI_COMB(1, 4), PNSZ_MULTI_COMB(1))
    allocate (XYZ_GLOBAL_MULTI_COMB_P(1, 3), iXYZ_GLOBAL_MULTI_COMB_P(1, 3), DS_MULTI_COMB(1), NELEM_TOTAL, iNELEM_TOTAL)
    allocate (RKBN_MULTI_COMB(1, 1, 1, 4), CGRN_MULTI_COMB(1, 1, 1, 4), PKBN_MULTI_COMB(1, 1, 1, 4), DGRN_MULTI_COMB(1, 1, 1, 4))
    XYZ_GLOBAL_MULTI_COMB = XYZ
    NCN_MULTI_COMB = NCN
    NCON_MULTI_COMB = NCON
    PNSZ_MULTI_COMB = PNSZ
    XYZ_GLOBAL_MULTI_COMB_P = XYZ_P
    iXYZ_GLOBAL_MULTI_COMB_P = iXYZ_P
    DS_MULTI_COMB = DS
    NELEM_TOTAL = 1
    iNELEM_TOTAL = 1
    call CALGREEN_MULTI_IRR
    wvn = 0d0
    call WAVE_PART_PANEL_MULTI(p, 1, V_SOLVER, wvn, NK, -1d0, gint)
    call check(error, maxval(abs(DGRN_MULTI_COMB(1, 1, 1, :) - gint/DS(1))) < 1d-6*maxval(abs(gint/DS(1))), &
               "the multi-body lid row does not hold the panel integral")
  end subroutine test_lid_row_holds_the_integral_multi

  !> Hull rows are not changed, and a lid point far from the panel compared with its size keeps the one-point value.
  subroutine test_hull_row_and_far_lid_point_keep_one_point(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: p(3)

    p = [3d0*ph, 0.1d0, 0d0]
    call set_solver(p)
    call CALGREEN_IRR
    call check(error, all(CGRN(1, 1, 1, :) == one_point(XYZ_P(1, :), 1)), "the hull row changed")
    if (allocated(error)) return
    call check(error, all(DGRN(1, 1, 1, :) == one_point(p, 3)), "the far lid point is not the one-point value")
  end subroutine test_hull_row_and_far_lid_point_keep_one_point

end module test_lid_wave_part
