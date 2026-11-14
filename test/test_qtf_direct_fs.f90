!> Unit tests of QTFDirectFS (Phase 7b, Task 3): the free-surface part of the right side of the direct second-order solve,
!> r_FS = s_D / g [sum_j gfs(i, j) w_j Q_D(xi_j) + sum_mn a(m, n, i) mom(m, n)]. The inner Green matrix against Gauss
!> quadrature of FINGREEN3D, the Graf amplitudes against FINGREEN3D outside the circle, the right side against a brute-force
!> sum, the half mesh against the whole one, and the order cut of long waves. The tests set module globals (the mesh arrays,
!> the wave, the lid), so the suite runs serially (see unit_tester.f90).
module test_qtf_direct_fs
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use HAMS_mod, only : NTHREAD, IRSP
  use WaveDyn_mod, only : H, NK, V, W1, WK, WVN
  use PanelMesh_mod, only : ISX, ISY, NSYS
  use Inerfs_mod, only : iNELEM, iXYZ_P, iDXYZ_P
  use Const_mod, only : RXY, g_const => G
  use SecondOrderTypes, only : so_options, panel_mesh, surface_geom, cs_points
  use SurfaceGeometry, only : build_surface_geom, DEFAULT_SHARP_DEG
  use ControlSurfaceMesh, only : build_fs_inner
  use FinGreen3D_Open, only : FINGREEN3D
  use PartialWaves, only : pw_bessel_k, pw_green_mmax
  use QTFDirectFS, only : direct_fs_data, direct_fs_setup, direct_fs_rhs, direct_fs_rfar, DIRECT_FS_SIGN, DIRECT_FS_TOL
  use second_order_fixtures, only : make_cylinder
  use test_waterline_field, only : set_wave, load_panels
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_positive_inf, ieee_quiet_nan
  implicit none
  private

  public :: collect_qtf_direct_fs

  real(8), parameter :: pi = 3.14159265358979323846d0
  ! The 448-panel cylinder of the QTF fixtures (radius 10, draft 10) in depth 15, the inner mesh of build_fs_inner for the
  ! largest first-order frequency fx_omax (R = 15 by the default rule), at the difference frequency fx_omega
  real(8), parameter :: fx_depth = 15d0, fx_omax = 1.1d0, fx_omega = 0.3d0
  ! Hull rows above this z are the waterline rows (the top row of the wall, 1.25 m high)
  real(8), parameter :: wl_row_z = -1.25d0

contains

subroutine collect_qtf_direct_fs(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("fs_green_matrix_matches", test_green_matrix), &
    new_unittest("fs_green_near_waterline", test_green_near_waterline), &
    new_unittest("fs_rhs_inner_only", test_rhs_inner_only), &
    new_unittest("fs_rhs_symmetry", test_rhs_symmetry), &
    new_unittest("fs_amplitudes_match_green", test_amplitudes_match_green), &
    new_unittest("fs_rhs_order_cut_finite", test_order_cut_finite), &
    new_unittest("fs_rfar_evanescent", test_rfar) &
    ]
end subroutine collect_qtf_direct_fs

! ------------------------------------------------------------------------------------------------------------- fixtures

!> The solver globals at the difference frequency omega in depth (finite), after set_wave (one thread, no lid)
subroutine set_omega(omega, depth)
  real(8), intent(in) :: omega, depth
  external :: DISPERSION

  call set_wave(depth)
  W1 = omega
  H = depth
  V = omega**2 / g_const
  call DISPERSION(WVN, NK, W1, H)
  WK = WVN(1)
end subroutine set_omega

!> The cylinder (radius 10, draft 10, ntheta x nz wall, nr bottom rings) whole (isy_in = 0) or half (1) in depth fx_depth,
!> its geometry, and the inner free surface of build_fs_inner of the whole mesh (radius: qtf_fs_radius, 0 for the default)
subroutine fixture(ntheta, nz, nr, isy_in, fs_radius, mesh, pts, fm, centre, radius, error)
  integer, intent(in) :: ntheta, nz, nr, isy_in
  real(8), intent(in) :: fs_radius
  type(panel_mesh), intent(out) :: mesh, fm
  type(cs_points), intent(out) :: pts
  real(8), intent(out) :: centre(2), radius
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: full, whole
  type(surface_geom) :: geom
  type(so_options) :: opts
  logical :: ok
  character(len=:), allocatable :: msg

  call make_cylinder(10d0, 10d0, ntheta, nz, nr, full, .true.)
  ! The inner mesh from the whole mesh, for both
  call load_panels(full, 0, 0, whole)
  call build_surface_geom(whole, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "geometry: "//msg)
  if (allocated(error)) return
  opts%qtf_fs_radius = fs_radius
  call build_fs_inner(whole, geom, opts, fx_depth, 'fs_direct_no_input_dir', pts, centre, radius, ok, msg, &
                      omega_max=fx_omax, fs_panels=fm)
  call check(error, ok, "inner mesh: "//msg)
  if (allocated(error)) return
  call check(error, fm%n == pts%n, "fs_panels does not have one panel per inner point")
  if (allocated(error)) return
  if (isy_in == 0) then
    mesh = whole
  else
    call load_panels(full, 0, isy_in, mesh)
  end if
end subroutine fixture

!> IRSP = 1 with the lid points pts (3, nl) in the single-body arrays; clear_lid undoes it
subroutine set_lid(pts)
  real(8), intent(in) :: pts(:, :)
  IRSP = 1
  iNELEM = size(pts, 2)
  if (allocated(iXYZ_P)) deallocate(iXYZ_P)
  if (allocated(iDXYZ_P)) deallocate(iDXYZ_P)
  allocate(iXYZ_P(iNELEM, 3), iDXYZ_P(iNELEM, 6))
  iXYZ_P = transpose(pts)
  iDXYZ_P = 0d0
end subroutine set_lid

subroutine clear_lid()
  IRSP = 0
  iNELEM = 0
  ISX = 0
  ISY = 0
  NSYS = 1
end subroutine clear_lid

!> Lid points on the waterplane of the cylinder: radii 3, 6, 8.5 at angles pi/12 + k pi/6 (none on y = 0); with
!> half = .true. only those with y > 0
subroutine lid_points(half, pts)
  logical, intent(in) :: half
  real(8), allocatable, intent(out) :: pts(:, :)
  real(8), parameter :: rl(3) = [3d0, 6d0, 8.5d0]
  real(8) :: th
  integer :: i, k, n

  allocate(pts(3, 36))
  n = 0
  do i = 1, 3
    do k = 0, 11
      th = pi / 12d0 + k * pi / 6d0
      if (half .and. sin(th) <= 0d0) cycle
      n = n + 1
      pts(:, n) = [rl(i) * cos(th), rl(i) * sin(th), 0d0]
    end do
  end do
  pts = pts(:, 1:n)
end subroutine lid_points

!> n-point Gauss-Legendre on [-1, 1] (Newton from the Chebyshev guesses, A&S 25.4.29)
subroutine gauss(n, x, w)
  integer, intent(in) :: n
  real(8), intent(out) :: x(n), w(n)
  real(8) :: z, p0, p1, p2, dp
  integer :: i, k, it
  do i = 1, n
    z = cos(pi * (i - 0.25d0) / (n + 0.5d0))
    do it = 1, 100
      p0 = 1d0
      p1 = z
      do k = 2, n
        p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k
        p0 = p1
        p1 = p2
      end do
      dp = n * (z * p1 - p0) / (z * z - 1d0)
      z = z - p1 / dp
      if (abs(p1 / dp) < 1d-15) exit
    end do
    x(i) = z
    w(i) = 2d0 / ((1d0 - z * z) * dp * dp)
  end do
end subroutine gauss

!> int G(x, xi) dS over panel p of fm by n x n Gauss-Legendre on its bilinear map (a triangle as a quadrilateral with
!> the last vertex twice), G from FINGREEN3D with flag 0 (the full Green function)
complex(8) function green_quadrature(fm, p, x, n) result(s)
  type(panel_mesh), intent(in) :: fm
  integer, intent(in) :: p, n
  real(8), intent(in) :: x(3)
  real(8) :: gx(n), gw(n), xv(3, 4), xq(3), du(3), ds(3), c(3), u, t
  complex(8) :: grn(4)
  integer :: k, iu, iv

  call gauss(n, gx, gw)
  do k = 1, 4
    xv(:, k) = fm%xnode(:, fm%ncon(min(k, fm%ncn(p)), p))
  end do
  s = (0d0, 0d0)
  do iu = 1, n
    u = gx(iu)
    do iv = 1, n
      t = gx(iv)
      xq = 0.25d0 * ((1 - u) * (1 - t) * xv(:, 1) + (1 + u) * (1 - t) * xv(:, 2) + (1 + u) * (1 + t) * xv(:, 3) &
                     + (1 - u) * (1 + t) * xv(:, 4))
      du = 0.25d0 * ((1 - t) * (xv(:, 2) - xv(:, 1)) + (1 + t) * (xv(:, 3) - xv(:, 4)))
      ds = 0.25d0 * ((1 - u) * (xv(:, 4) - xv(:, 1)) + (1 + u) * (xv(:, 3) - xv(:, 2)))
      c = [du(2) * ds(3) - du(3) * ds(2), du(3) * ds(1) - du(1) * ds(3), du(1) * ds(2) - du(2) * ds(1)]
      call FINGREEN3D(xq(1), x(1), xq(2), x(2), xq(3), x(3), V, WVN, NK, H, grn, 0)
      s = s + grn(1) * gw(iu) * gw(iv) * norm2(c)
    end do
  end do
end function green_quadrature

complex(8) function green_point(x, xi)
  real(8), intent(in) :: x(3), xi(3)
  complex(8) :: grn(4)
  call FINGREEN3D(xi(1), x(1), xi(2), x(2), xi(3), x(3), V, WVN, NK, H, grn, 0)
  green_point = grn(1)
end function green_point

!> A reproducible complex number from a seed
complex(8) function pseudo(seed)
  integer, intent(in) :: seed
  pseudo = cmplx(sin(1.3d0 * seed + 0.2d0), cos(0.7d0 * seed + 1d0), 8)
end function pseudo

! ------------------------------------------------------------------------------------------------------------- tests

!> The inner Green matrix at 5 (row, panel) pairs of the 448-panel cylinder away from the waterline (rows below the top
!> row of the wall), all near pairs (within 50 panel sizes, as every pair of this small free surface): gfs(i, j) against
!> (1 / w_j) int_j G dS by 20 x 20 Gauss-Legendre of FINGREEN3D (flag 0), 1e-6 relative: the Rankine part is in closed
!> form, the rest by the order of fs_remainder_order for an error of 1e-8 of the rest; the gate is that of FINGREEN3D
!> itself, which is about 4e-7 off the eigenfunction series and switches series with the distance (ledger ruling F5), so
!> two quadratures of it agree no better (measured: 2e-14 to 1.0e-7 at the five pairs). And the far branch
!> (beyond 50 panel sizes: G at the centroid, as CALGREEN): a small panel (0.05 m) appended to the inner mesh at r = 14 m,
!> against FINGREEN3D at its centroid for 5 rows on the far side of the hull, 1e-12.
subroutine test_green_matrix(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, fm, fm2
  type(cs_points) :: pts, pts2
  type(direct_fs_data) :: fs
  real(8) :: centre(2), radius, err, worst_near, worst_far, rel
  complex(8) :: ref
  integer :: rows(5), cols(5), k, i, j, nn
  character(len=200) :: info

  call fixture(32, 8, 6, 0, 0d0, mesh, pts, fm, centre, radius, error)
  if (allocated(error)) return
  call check(error, mesh%n == 448, "the fixture is not the 448-panel cylinder")
  if (allocated(error)) return
  call set_omega(fx_omega, fx_depth)
  ! Append a small square panel at (14, 0, 0): far (50 x 0.05 = 2.5 m) from the rows of the far side of the hull
  nn = fm%nnode
  fm2%n = fm%n + 1
  fm2%nnode = nn + 4
  allocate(fm2%xnode(3, nn + 4), fm2%ncn(fm2%n), fm2%ncon(4, fm2%n), fm2%nrm(3, fm2%n), fm2%xc(3, fm2%n), &
           fm2%area(fm2%n), fm2%psize(fm2%n))
  fm2%xnode(:, 1:nn) = fm%xnode
  fm2%xnode(:, nn + 1) = [13.975d0, -0.025d0, 0d0]
  fm2%xnode(:, nn + 2) = [14.025d0, -0.025d0, 0d0]
  fm2%xnode(:, nn + 3) = [14.025d0, 0.025d0, 0d0]
  fm2%xnode(:, nn + 4) = [13.975d0, 0.025d0, 0d0]
  fm2%ncn(1:fm%n) = fm%ncn
  fm2%ncon(:, 1:fm%n) = fm%ncon
  fm2%nrm(:, 1:fm%n) = fm%nrm
  fm2%ncn(fm2%n) = 4
  fm2%ncon(:, fm2%n) = [nn + 1, nn + 2, nn + 3, nn + 4]
  fm2%nrm(:, fm2%n) = [0d0, 0d0, -1d0]
  pts2%n = pts%n + 1
  allocate(pts2%x(3, pts2%n), pts2%nrm(3, pts2%n), pts2%w(pts2%n), pts2%psize(pts2%n), pts2%kind(pts2%n), &
           pts2%body(pts2%n))
  pts2%x(:, 1:pts%n) = pts%x
  pts2%nrm(:, 1:pts%n) = pts%nrm
  pts2%w(1:pts%n) = pts%w
  pts2%psize(1:pts%n) = pts%psize
  pts2%kind(1:pts%n) = pts%kind
  pts2%body(1:pts%n) = pts%body
  pts2%x(:, pts2%n) = [14d0, 0d0, 0d0]
  pts2%nrm(:, pts2%n) = [0d0, 0d0, -1d0]
  pts2%w(pts2%n) = 0.05d0**2
  pts2%psize(pts2%n) = 0.05d0
  pts2%kind(pts2%n) = pts%kind(1)
  pts2%body(pts2%n) = 1
  call direct_fs_setup(mesh, pts2, fm2, centre, radius, fx_depth, fx_omega, .false., fs)
  call check(error, fs%nrow == mesh%n .and. fs%nfs == pts2%n .and. size(fs%gfs, 1) == mesh%n .and. &
             size(fs%gfs, 2) == pts2%n, "gfs does not have the rows of the hull and the columns of the inner mesh")
  if (allocated(error)) return

  ! Near pairs away from the waterline: a bottom panel, three wall panels below the top row, a bottom panel near the rim
  rows = 0
  k = 0
  do i = 1, mesh%n
    if (mesh%xc(3, i) > wl_row_z) cycle
    if (mod(i, 37) /= 5) cycle
    k = k + 1
    rows(k) = i
    if (k == 5) exit
  end do
  call check(error, k == 5, "fewer than 5 rows away from the waterline")
  if (allocated(error)) return
  cols = [1, pts%n / 4, pts%n / 2, (3 * pts%n) / 4, pts%n]
  worst_near = 0d0
  do k = 1, 5
    i = rows(k)
    j = cols(k)
    ref = green_quadrature(fm, j, mesh%xc(:, i), 20) / pts%w(j)
    err = abs(fs%gfs(i, j) - ref) / abs(ref)
    print '(4x,a,i0,a,i0,a,3f8.3,a,es10.2)', "row ", i, ", inner panel ", j, ", row at", mesh%xc(:, i), &
      ": near rule against 20 x 20 Gauss, relative ", err
    worst_near = max(worst_near, err)
  end do
  ! The far branch
  worst_far = 0d0
  k = 0
  do i = 1, mesh%n
    if (mesh%xc(1, i) > -5d0 .or. norm2(mesh%xc(:, i) - pts2%x(:, pts2%n)) <= 50d0 * 0.05d0) cycle
    if (mod(i, 13) /= 0) cycle
    k = k + 1
    ref = green_point(mesh%xc(:, i), pts2%x(:, pts2%n))
    rel = abs(fs%gfs(i, pts2%n) - ref) / abs(ref)
    worst_far = max(worst_far, rel)
    if (k == 5) exit
  end do
  write (info, '(a,es10.3,a,es10.3,a,i0,a)') "gfs: near pairs against 20 x 20 Gauss ", worst_near, &
    ", far pairs against FINGREEN3D at the centroid ", worst_far, " (", k, " rows)"
  print '(4x,a)', trim(info)
  call check(error, k == 5 .and. worst_near <= 1d-6 .and. worst_far <= 1d-12, trim(info))
end subroutine test_green_matrix

!> The near rule at the waterline: the hull rows of the top row of the wall (centroids 0.625 m below z = 0) against the
!> inner panels next to them (the inner panel nearest to each row), gfs(i, j) w_j (Rankine part in closed form, the rest
!> by Gauss-Legendre) against 20 x 20 Gauss-Legendre of FINGREEN3D over the panel: 1e-3; the one-point value
!> G(x_i, xi_j) w_j is off by more than 1e-2 (the near rule matters there).
subroutine test_green_near_waterline(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, fm
  type(cs_points) :: pts
  type(direct_fs_data) :: fs
  real(8) :: centre(2), radius, worst, best_one, e_rule, e_one
  complex(8) :: ref, one
  integer :: i, j, nrow
  character(len=200) :: info

  call fixture(32, 8, 6, 0, 0d0, mesh, pts, fm, centre, radius, error)
  if (allocated(error)) return
  call set_omega(fx_omega, fx_depth)
  call direct_fs_setup(mesh, pts, fm, centre, radius, fx_depth, fx_omega, .false., fs)
  worst = 0d0
  best_one = huge(1d0)
  nrow = 0
  do i = 1, mesh%n
    if (mesh%xc(3, i) <= wl_row_z .or. abs(mesh%nrm(3, i)) > 0.5d0) cycle
    if (mod(i, 4) /= 1) cycle
    j = minloc(norm2(pts%x - spread(mesh%xc(:, i), 2, pts%n), dim=1), 1)
    ref = green_quadrature(fm, j, mesh%xc(:, i), 20)
    one = green_point(mesh%xc(:, i), pts%x(:, j)) * pts%w(j)
    e_rule = abs(fs%gfs(i, j) * pts%w(j) - ref) / abs(ref)
    e_one = abs(one - ref) / abs(ref)
    worst = max(worst, e_rule)
    best_one = min(best_one, e_one)
    nrow = nrow + 1
  end do
  write (info, '(i0,a,es10.3,a,es10.3)') nrow, " waterline rows: near rule against 20 x 20 Gauss, worst ", worst, &
    "; one point, smallest error ", best_one
  print '(4x,a)', trim(info)
  call check(error, nrow >= 4 .and. worst <= 1d-3 .and. best_one > 1d-2, trim(info))
end subroutine test_green_near_waterline

!> Ledger ruling F2: mom = 0 and Q_D = 1 at the inner points. r_FS (whole mesh, no lid) against the brute-force sum
!> s_D / g sum_j int_j G(x_i, .) dS by 10 x 10 Gauss-Legendre of FINGREEN3D on each inner panel: 1e-3 at the waterline
!> rows, 1e-6 at the others (every third row).
subroutine test_rhs_inner_only(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, fm
  type(cs_points) :: pts
  type(direct_fs_data) :: fs
  complex(8), allocatable :: q(:, :), mom(:, :, :), rhs(:, :, :)
  complex(8) :: ref
  real(8) :: centre(2), radius, e, worst_wl, worst_in
  integer :: i, j, nwl, nin
  character(len=200) :: info

  call fixture(32, 8, 6, 0, 0d0, mesh, pts, fm, centre, radius, error)
  if (allocated(error)) return
  call set_omega(fx_omega, fx_depth)
  call direct_fs_setup(mesh, pts, fm, centre, radius, fx_depth, fx_omega, .false., fs)
  allocate(q(pts%n, 1), mom(-fs%mg:fs%mg, 0:fs%nev, 1), rhs(mesh%nbase, 1, 1))
  q = (1d0, 0d0)
  mom = (0d0, 0d0)
  rhs = (0d0, 0d0)
  call direct_fs_rhs(fs, q, mom, rhs)
  worst_wl = 0d0
  worst_in = 0d0
  nwl = 0
  nin = 0
  do i = 1, mesh%n
    if (mod(i, 3) /= 0 .and. mesh%xc(3, i) <= wl_row_z) cycle
    ref = (0d0, 0d0)
    do j = 1, pts%n
      ref = ref + green_quadrature(fm, j, mesh%xc(:, i), 10)
    end do
    ref = DIRECT_FS_SIGN / g_const * ref
    e = abs(rhs(i, 1, 1) - ref) / abs(ref)
    if (mesh%xc(3, i) > wl_row_z) then
      worst_wl = max(worst_wl, e)
      nwl = nwl + 1
    else
      worst_in = max(worst_in, e)
      nin = nin + 1
    end if
  end do
  write (info, '(a,i0,a,es10.3,a,i0,a,es10.3)') "Q_D = 1 inside: against 10 x 10 Gauss, ", nwl, " waterline rows ", &
    worst_wl, ", ", nin, " other rows ", worst_in
  print '(4x,a)', trim(info)
  call check(error, nwl > 0 .and. nin > 0 .and. worst_wl <= 1d-3 .and. worst_in <= 1d-6, trim(info))
end subroutine test_rhs_inner_only

!> The half mesh (ISY = 1, its lid points y > 0) against the whole mesh, both with a lid, the same inner mesh and the same
!> Q_D and moments (arbitrary, two columns): the full rows of the half mesh, sum_ip RXY(im, ip) rhs(i, :, ip) / nsys
!> (as RADIATION_SOLVER expands the blocks), equal the rows of the whole mesh at the same points, 1e-7 relative to
!> the largest. The mirrored centroids of the half mesh and the centroids of the whole one differ by rounding (the nodes
!> at -theta come from their own sin and cos), and the series of FINGREEN3D, with the epsilon algorithm of Limes, jump by
!> far more than rounding when a point moves by rounding (FinGreen3D.f90 at LINTON; half_limit of test_qtf_direct is
!> 1e-8 for the same reason). Measured 1.1e-8; a wrong block layout gives order 1.
subroutine test_rhs_symmetry(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh_w, mesh_h, fm
  type(cs_points) :: pts
  type(direct_fs_data) :: fs_w, fs_h
  complex(8), allocatable :: q(:, :), mom(:, :, :), rw(:, :, :), rh(:, :, :)
  real(8), allocatable :: lid_w(:, :), lid_h(:, :)
  real(8) :: centre(2), radius, worst, big, xl(3), moved
  complex(8) :: vv(2)
  integer :: i, j, im, ip, nb, nl, nlw, m, n, col, p
  character(len=200) :: info

  call fixture(32, 8, 6, 0, 0d0, mesh_w, pts, fm, centre, radius, error)
  if (allocated(error)) return
  call set_omega(fx_omega, fx_depth)
  call lid_points(.false., lid_w)
  call set_lid(lid_w)
  call direct_fs_setup(mesh_w, pts, fm, centre, radius, fx_depth, fx_omega, .false., fs_w)
  nlw = size(lid_w, 2)
  allocate(q(pts%n, 2), mom(-fs_w%mg:fs_w%mg, 0:fs_w%nev, 2), rw(mesh_w%nbase + nlw, 2, 1))
  do col = 1, 2
    do j = 1, pts%n
      q(j, col) = pseudo(j + 1000 * col)
    end do
    do n = 0, fs_w%nev
      do m = -fs_w%mg, fs_w%mg
        mom(m, n, col) = pseudo(7 * m + 31 * n + 500 * col) / (1d0 + abs(m) + n)**2
      end do
    end do
  end do
  rw = (0d0, 0d0)
  call direct_fs_rhs(fs_w, q, mom, rw)

  call fixture(32, 8, 6, 1, 0d0, mesh_h, pts, fm, centre, radius, error)
  if (allocated(error)) return
  call set_omega(fx_omega, fx_depth)
  call lid_points(.true., lid_h)
  call set_lid(lid_h)
  call direct_fs_setup(mesh_h, pts, fm, centre, radius, fx_depth, fx_omega, .false., fs_h)
  call clear_lid()
  call check(error, fs_h%mg == fs_w%mg .and. fs_h%nev == fs_w%nev .and. fs_h%nsys == 2, &
             "the half mesh does not have the orders of the whole one, or not two blocks")
  if (allocated(error)) return
  nb = mesh_h%nbase
  nl = size(lid_h, 2)
  allocate(rh(nb + nl, 2, 2))
  rh = (0d0, 0d0)
  call direct_fs_rhs(fs_h, q, mom, rh)

  big = maxval(abs(rw))
  worst = 0d0
  moved = 0d0
  do im = 1, 2
    do i = 1, nb + nl
      vv = (0d0, 0d0)
      do ip = 1, 2
        vv = vv + RXY(im, ip) * rh(i, :, ip)
      end do
      vv = vv / 2d0
      if (i <= nb) then
        xl = mesh_h%xc(:, (im - 1) * nb + i)
        p = minloc(norm2(mesh_w%xc - spread(xl, 2, mesh_w%n), dim=1), 1)
        if (norm2(mesh_w%xc(:, p) - xl) > 1d-9) then
          call check(error, .false., "a panel of the half mesh is not a panel of the whole mesh")
          return
        end if
        moved = max(moved, norm2(mesh_w%xc(:, p) - xl))
      else
        xl = lid_h(:, i - nb)
        if (im == 2) xl(2) = -xl(2)
        p = minloc(norm2(lid_w - spread(xl, 2, nlw), dim=1), 1)
        if (norm2(lid_w(:, p) - xl) > 1d-9) then
          call check(error, .false., "a lid point of the half mesh is not a lid point of the whole mesh")
          return
        end if
        p = mesh_w%nbase + p
      end if
      worst = max(worst, maxval(abs(vv - rw(p, :, 1))) / big)
    end do
  end do
  write (info, '(a,i0,a,i0,a,es10.3,a,es10.3,a)') "half mesh (", nb, " base panels, ", nl, &
    " lid points) against the whole mesh, relative to the largest ", worst, " (centroids apart by up to ", moved, " m)"
  print '(4x,a)', trim(info)
  call check(error, worst <= 1d-7, trim(info))
end subroutine test_rhs_symmetry

!> The Graf amplitudes in the contraction form (ruling F1): sum_mn a(m, n, i) B_mn(xi), B_m0 = H_m(k0 r) e^{im theta},
!> B_mn = K_m(kn r) / K_m(kn R) e^{im theta}, |m| <= mcut, is G(x_i, xi) of FINGREEN3D at points xi on z = 0 at r = R
!> and 1.4 R, for the row nearest the circle (a waterline row), a bottom row and a lid point next to the waterline:
!> 1e-5 relative (FINGREEN3D is about 4e-7 off the eigenfunction series, ledger ruling F5). The orders follow the rules
!> of the header: mg = pw_green_mmax(r_max, R, DIRECT_FS_TOL) and the evanescent modes past nev below DIRECT_FS_TOL.
subroutine test_amplitudes_match_green(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, fm
  type(cs_points) :: pts
  type(direct_fs_data) :: fs
  real(8), allocatable :: lid(:, :), kv(:), kr(:)
  real(8) :: centre(2), radius, xi(3), r, th, worst, d, x(3)
  complex(8) :: s, ref, hm
  integer :: rows(3), ir, it, ifac, m, n, i, row
  character(len=200) :: info

  call fixture(32, 8, 6, 0, 0d0, mesh, pts, fm, centre, radius, error)
  if (allocated(error)) return
  call set_omega(fx_omega, fx_depth)
  allocate(lid(3, 1))
  lid(:, 1) = [9.5d0 * cos(0.3d0), 9.5d0 * sin(0.3d0), 0d0]
  call set_lid(lid)
  call direct_fs_setup(mesh, pts, fm, centre, radius, fx_depth, fx_omega, .false., fs)
  call clear_lid()
  d = radius - fs%r_max
  call check(error, fs%mg == pw_green_mmax(fs%r_max, radius, DIRECT_FS_TOL) .and. fs%nev >= 1, "mg is not the rule")
  if (allocated(error)) return
  call check(error, exp(-fs%kev(fs%nev) * d) <= DIRECT_FS_TOL * (1d0 - exp(-pi * d / fx_depth)), &
             "the evanescent modes past nev are not below DIRECT_FS_TOL")
  if (allocated(error)) return
  print '(4x,a,f8.4,a,f8.4,a,i0,a,i0,a,i0)', "R ", radius, ", r_max ", fs%r_max, ": mg ", fs%mg, ", nev ", fs%nev, &
    ", mcut ", fs%mcut
  rows(1) = maxloc(norm2(mesh%xc(1:2, :) - spread(centre, 2, mesh%n), dim=1), 1)
  rows(2) = minloc(mesh%xc(3, :) + 0.01d0 * norm2(mesh%xc(1:2, :), dim=1), 1)
  rows(3) = mesh%n + 1
  allocate(kv(0:fs%mg), kr(0:fs%mg))
  worst = 0d0
  do ir = 1, 3
    row = rows(ir)
    if (row <= mesh%n) then
      x = mesh%xc(:, row)
    else
      x = lid(:, 1)
    end if
    do ifac = 1, 2
      r = radius * merge(1d0, 1.4d0, ifac == 1)
      do it = 0, 2
        th = 0.4d0 + 2.1d0 * it
        xi = [centre(1) + r * cos(th), centre(2) + r * sin(th), 0d0]
        s = (0d0, 0d0)
        do m = -fs%mcut, fs%mcut
          i = abs(m)
          hm = cmplx(bessel_jn(i, fs%k0 * r), bessel_yn(i, fs%k0 * r), 8)
          if (m < 0 .and. mod(i, 2) == 1) hm = -hm
          s = s + fs%a(m, 0, row) * hm * exp(cmplx(0d0, m * th, 8))
        end do
        do n = 1, fs%nev
          call pw_bessel_k(fs%mg, fs%kev(n) * r, kv)
          call pw_bessel_k(fs%mg, fs%kev(n) * radius, kr)
          do m = -fs%mcut, fs%mcut
            s = s + fs%a(m, n, row) * kv(abs(m)) / kr(abs(m)) * exp(cmplx(0d0, m * th, 8))
          end do
        end do
        ref = green_point(x, xi)
        worst = max(worst, abs(s - ref) / abs(ref))
      end do
    end do
    print '(4x,a,3f9.4,a,es10.3)', "row at", x, ": worst so far ", worst
  end do
  write (info, '(a,es10.3)') "Graf amplitudes against FINGREEN3D outside R, worst relative ", worst
  print '(4x,a)', trim(info)
  call check(error, worst <= 1d-5, trim(info))
end subroutine test_amplitudes_match_green

!> A long wave on a circle close to the hull: Omega = 0.05 rad/s, R = 11 m about the cylinder of radius 10 (coarse,
!> 16 x 4 wall panels), so that H_m(k0 R) overflows below the order the waterline rows need (mcut < mg). The amplitudes
!> above mcut are zero, the warning measure cut_error is set, and with moments that are Inf and NaN above mcut (as
!> fs_outer_moments gives there) and finite below, the right side stays finite (no Inf x 0).
subroutine test_order_cut_finite(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, fm
  type(cs_points) :: pts
  type(direct_fs_data) :: fs
  complex(8), allocatable :: q(:, :), mom(:, :, :), rhs(:, :, :)
  real(8) :: centre(2), radius, inf, nan
  integer :: m, n, j
  character(len=200) :: info
  real(8), parameter :: om = 0.05d0

  call fixture(16, 4, 3, 0, 11d0, mesh, pts, fm, centre, radius, error)
  if (allocated(error)) return
  call set_omega(om, fx_depth)
  call direct_fs_setup(mesh, pts, fm, centre, radius, fx_depth, om, .false., fs)
  write (info, '(a,f6.2,a,f8.4,a,i0,a,i0,a,i0,a,es10.3)') "R ", radius, ", r_max ", fs%r_max, ": mg ", fs%mg, ", mcut ", &
    fs%mcut, ", nev ", fs%nev, ", cut_error ", fs%cut_error
  print '(4x,a)', trim(info)
  call check(error, fs%mcut < fs%mg .and. fs%cut_error > 0d0, "the order is not cut: "//trim(info))
  if (allocated(error)) return
  call check(error, all(fs%a(fs%mcut + 1:, :, :) == (0d0, 0d0)) .and. all(fs%a(:-fs%mcut - 1, :, :) == (0d0, 0d0)), &
             "amplitudes above mcut are not zero")
  if (allocated(error)) return
  inf = ieee_value(1d0, ieee_positive_inf)
  nan = ieee_value(1d0, ieee_quiet_nan)
  allocate(q(pts%n, 1), mom(-fs%mg:fs%mg, 0:fs%nev, 1), rhs(mesh%nbase, 1, 1))
  do j = 1, pts%n
    q(j, 1) = pseudo(j)
  end do
  do n = 0, fs%nev
    do m = -fs%mg, fs%mg
      if (abs(m) <= fs%mcut) then
        mom(m, n, 1) = pseudo(m + 17 * n) / (1d0 + abs(m) + n)**2
      else if (mod(m, 2) == 0) then
        mom(m, n, 1) = cmplx(inf, 0d0, 8)
      else
        mom(m, n, 1) = cmplx(nan, nan, 8)
      end if
    end do
  end do
  rhs = (0d0, 0d0)
  call direct_fs_rhs(fs, q, mom, rhs)
  call check(error, all(ieee_is_finite(real(rhs))) .and. all(ieee_is_finite(aimag(rhs))) .and. maxval(abs(rhs)) > 0d0, &
             "the right side is not finite with moments that are not finite above mcut")
end subroutine test_order_cut_finite

!> direct_fs_rfar raises r_far to R + 40 / kev(1), where the slowest evanescent basis function has decayed by e^-40,
!> and keeps a larger r_far
subroutine test_rfar(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, fm
  type(cs_points) :: pts
  type(direct_fs_data) :: fs
  real(8) :: centre(2), radius, r1, r2
  character(len=200) :: info

  call fixture(16, 4, 3, 0, 0d0, mesh, pts, fm, centre, radius, error)
  if (allocated(error)) return
  call set_omega(fx_omega, fx_depth)
  call direct_fs_setup(mesh, pts, fm, centre, radius, fx_depth, fx_omega, .false., fs)
  r1 = direct_fs_rfar(fs, radius + 1d0)
  r2 = direct_fs_rfar(fs, radius + 1d3)
  write (info, '(a,f10.4,a,f10.4,a,es10.3)') "r_far from R + 1: ", r1, ", from R + 1000: ", r2, &
    "; kev(1) (r_far - R) ", fs%kev(1) * (r1 - radius)
  print '(4x,a)', trim(info)
  call check(error, abs(fs%kev(1) * (r1 - radius) - 40d0) <= 1d-12 .and. r2 == radius + 1d3, trim(info))
end subroutine test_rfar

end module test_qtf_direct_fs
