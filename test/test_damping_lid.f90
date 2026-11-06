!> The 1/r integrals over damping lid panels, as sources, against the composite Gauss quadrature of test_panel_integrals, and
!> the single-layer integral of the wave Green function over a lid panel. The lid arrays and the wave quantities are module
!> data, so this suite runs serially.
module test_damping_lid

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use Lid_mod, only : LNCN, LNCON, LXYZ, LXYZ_P, LDS, LPNSZ, LNELEM, LNTND
  use WaveDyn_mod, only : V, H, WVN, NK
  use DampingLid, only : LID_GEOMETRY, LID_SINGLE_LAYER, LID_DEEP_WAVE
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
      new_unittest("lid_multi_matches_single", test_lid_multi_matches_single), &
      new_unittest("duffy_converges", test_duffy_converges), &
      new_unittest("far_rule_consistent", test_far_rule_consistent), &
      new_unittest("electrified_disc", test_electrified_disc), &
      new_unittest("finite_tends_to_deep", test_finite_tends_to_deep) &
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

  !> A lid mesh in the lid arrays, with its centroids, areas and sizes.
  subroutine set_lid_mesh(nodes, ncn, ncon)
    real(8), intent(in) :: nodes(:, :)
    integer, intent(in) :: ncn(:), ncon(:, :)

    if (allocated(LXYZ)) deallocate (LXYZ)
    if (allocated(LNCN)) deallocate (LNCN)
    if (allocated(LNCON)) deallocate (LNCON)
    if (allocated(LXYZ_P)) deallocate (LXYZ_P)
    if (allocated(LDS)) deallocate (LDS)
    if (allocated(LPNSZ)) deallocate (LPNSZ)
    LNTND = size(nodes, 1)
    LNELEM = size(ncn)
    allocate (LXYZ(LNTND, 3), LNCN(LNELEM), LNCON(LNELEM, 4))
    LXYZ = nodes
    LNCN = ncn
    LNCON = ncon
    call LID_GEOMETRY()
  end subroutine set_lid_mesh

  subroutine set_one_panel(nodes)
    real(8), intent(in) :: nodes(:, :)
    integer :: ncon(1, 4), k

    ncon = 0
    do k = 1, size(nodes, 1)
      ncon(1, k) = k
    end do
    call set_lid_mesh(nodes, [size(nodes, 1)], ncon)
  end subroutine set_one_panel

  !> Deep water (depth < 0) or finite depth at wave number nu, with the roots of the dispersion relation set as
  !> CalWaveProperts sets them.
  subroutine set_wave(nu, depth)
    real(8), intent(in) :: nu, depth
    external :: DISPERSION

    V = nu
    H = depth
    if (depth > 0d0 .and. abs(nu) > 1d-8 .and. abs(nu + 1d0) > 1d-8) call DISPERSION(WVN, NK, sqrt(9.807d0*nu), depth)
  end subroutine set_wave

  !> The 6 x 6 rule against 24 x 24, with the field point in the panel: at the centroid, at a node and off the centroid.
  subroutine test_duffy_converges(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: quad(4, 3), tri(3, 3), p(3, 3)
    complex(8) :: a, b
    integer :: shape, i

    quad = reshape([0d0, 1d0, 1d0, 0d0, 0d0, 0d0, 1d0, 1d0, 0d0, 0d0, 0d0, 0d0], [4, 3])
    tri = reshape([0d0, 1d0, 0.3d0, 0d0, 0d0, 0.9d0, 0d0, 0d0, 0d0], [3, 3])
    call set_wave(0.19d0, -1d0)
    do shape = 1, 2
      if (shape == 1) then
        call set_one_panel(quad)
        p(1, :) = LXYZ_P(1, :)
        p(2, :) = quad(3, :)
        p(3, :) = [0.7d0, 0.2d0, 0d0]
      else
        call set_one_panel(tri)
        p(1, :) = LXYZ_P(1, :)
        p(2, :) = tri(2, :)
        p(3, :) = [0.4d0, 0.2d0, 0d0]
      end if
      do i = 1, 3
        a = LID_DEEP_WAVE(p(i, :), 1, 6)
        b = LID_DEEP_WAVE(p(i, :), 1, 24)
        call check(error, abs(a - b) < 1d-6*abs(b), "shape "//char(48 + shape)//", point "//char(48 + i))
        if (allocated(error)) return
      end do
    end do
  end subroutine test_duffy_converges

  !> Across dist = 3 LPNSZ the one-point rule takes over from the Gauss rule. The one-point rule has an error of about
  !> (nu*LPNSZ)**2/24 at any distance, from the curvature of the wave, as everywhere else in the code: 1.5e-3 for the unit
  !> square at nu = 0.19 (test/data/damping_lid/README.md). So the square has side 0.5, nu*LPNSZ = 0.095.
  subroutine test_far_rule_consistent(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: quad(4, 3), p(3)
    complex(8) :: gauss, one

    quad = 0.5d0*reshape([0d0, 1d0, 1d0, 0d0, 0d0, 0d0, 1d0, 1d0, 0d0, 0d0, 0d0, 0d0], [4, 3])
    call set_wave(0.19d0, -1d0)
    call set_one_panel(quad)
    p = LXYZ_P(1, :) + [3d0*LPNSZ(1), 0d0, 0d0]
    gauss = LID_DEEP_WAVE(p, 1, 6)
    p = LXYZ_P(1, :) + [3d0*LPNSZ(1)*(1d0 + 1d-9), 0d0, 0d0]
    one = LID_DEEP_WAVE(p, 1, 6)
    call check(error, abs(gauss - one) < 1d-3*abs(gauss), "one point and Gauss differ at dist = 3 LPNSZ")
  end subroutine test_far_rule_consistent

  !> A disc of radius 1 at nu = 0, where the kernel on z = 0 is 2/r. The density of unit potential has total 1/pi
  !> (capacitance 2a/pi of the kernel 1/r). The density is singular at the rim, so the radial division sets the error:
  !> nt = nr = 8, 16 and 32, that is 64, 256 and 1024 panels.
  subroutine test_electrified_disc(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: nts(3) = [8, 16, 32], nrs(3) = [8, 16, 32]
    real(8), allocatable :: nodes(:, :), a(:, :), q(:)
    integer, allocatable :: ncn(:), ncon(:, :), ipiv(:)
    real(8) :: err(3)
    complex(8) :: sl
    integer :: m, i, j, n, info

    call set_wave(0d0, -1d0)
    do m = 1, 3
      call disc(1d0, nts(m), nrs(m), nodes, ncn, ncon)
      call set_lid_mesh(nodes, ncn, ncon)
      n = LNELEM
      allocate (a(n, n), q(n), ipiv(n))
      do j = 1, n
        do i = 1, n
          call LID_SINGLE_LAYER(LXYZ_P(i, :), j, sl)
          a(i, j) = dble(sl)
        end do
      end do
      q = 1d0
      call DGESV(n, 1, a, n, ipiv, q, n, info)
      call check(error, info == 0, "DGESV failed")
      if (allocated(error)) return
      err(m) = abs(sum(q*LDS)*acos(-1d0) - 1d0)
      deallocate (a, q, ipiv)
    end do
    call check(error, err(2) < err(1) .and. err(3) < err(2), "the charge error does not fall with refinement")
    if (allocated(error)) return
    call check(error, err(3) < 1d-2, "the charge error is above 1% on the finest mesh")
  end subroutine test_electrified_disc

  !> The disc r <= r1 on z = 0 as annulus() of make_lid_fixtures.py makes it: nt triangles about the centre node, then
  !> nt x (nr - 1) quads.
  subroutine disc(r1, nt, nr, nodes, ncn, ncon)
    real(8), intent(in) :: r1
    integer, intent(in) :: nt, nr
    real(8), allocatable, intent(out) :: nodes(:, :)
    integer, allocatable, intent(out) :: ncn(:), ncon(:, :)
    integer :: i, j, k, e
    real(8) :: r, t

    allocate (nodes(1 + nt*nr, 3), ncn(nt*nr), ncon(nt*nr, 4))
    nodes = 0d0
    do i = 1, nr
      r = r1*i/nr
      do j = 0, nt - 1
        t = 2d0*acos(-1d0)*j/nt
        nodes(2 + (i - 1)*nt + j, :) = [r*cos(t), r*sin(t), 0d0]
      end do
    end do
    ncon = 0
    e = 0
    do i = 0, nr - 1
      do j = 0, nt - 1
        k = mod(j + 1, nt)
        e = e + 1
        if (i == 0) then
          ncn(e) = 3
          ncon(e, 1:3) = [1, 2 + k, 2 + j]
        else
          ncn(e) = 4
          ncon(e, :) = [2 + (i - 1)*nt + j, 2 + (i - 1)*nt + k, 2 + i*nt + k, 2 + i*nt + j]
        end if
      end do
    end do
  end subroutine disc

  !> At large kh the finite-depth single layer of a self panel tends to the deep-water one. The test uses kh = 300, since
  !> FINGREEN3D returns NaN from kh = 360, where cosh(kh) overflows (test/data/damping_lid/README.md). V and H are restored.
  subroutine test_finite_tends_to_deep(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: quad(4, 3), p(3), v0, h0
    complex(8) :: deep, fin

    v0 = V
    h0 = H

    quad = reshape([0d0, 1d0, 1d0, 0d0, 0d0, 0d0, 1d0, 1d0, 0d0, 0d0, 0d0, 0d0], [4, 3])
    call set_one_panel(quad)
    p = LXYZ_P(1, :)
    call set_wave(1d0, -1d0)
    call LID_SINGLE_LAYER(p, 1, deep)
    call set_wave(1d0, 300d0)
    call LID_SINGLE_LAYER(p, 1, fin)
    V = v0
    H = h0
    call check(error, abs(fin - deep) < 1d-4*abs(deep), "finite depth at H = 300, nu = 1 differs from deep water")
  end subroutine test_finite_tends_to_deep

end module test_damping_lid
