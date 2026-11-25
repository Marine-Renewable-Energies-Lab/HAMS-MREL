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
  use INFG3D_Open, only : INFGREEN3D
  use direct_sum_green_mod, only : direct_sum_green
  use PartialWaves, only : pw_set, pw_bessel_k, pw_green_mmax, pw_eval_surface
  use QTFFreeSurface, only : fs_point_field, fs_forcing_diff, fs_outer_moments
  use QTFDirectFS, only : fs_remainder_order
  use QTFDirectFS, only : direct_fs_data, direct_fs_setup, direct_fs_rhs, direct_fs_rfar, DIRECT_FS_SIGN, DIRECT_FS_TOL, &
    direct_fs_estimate, direct_fs_array_mb
  use second_order_fixtures, only : make_cylinder, set_panel, maccamy_fuchs_state, read_pnl_hull
  use test_waterline_field, only : set_wave, load_panels
  ! The direct solve with the free-surface forcing (Task 4)
  use HAMS_mod, only : ISOLV
  use PanelMesh_mod, only : IPIV
  use Inerfs_mod, only : TNELEM
  use Potentials_mod, only : CGRN, RKBN, MXPOT
  use LinearMatrix_mod, only : AMAT, BRMAT
  use AssbMatx, only : ASSB_LEFT, ASSB_RBC, RADIATION_SOLVER
  use CalGreenFunc, only : CALGREEN
  use ImplementSubs, only : SetWaveFromOmega
  use SecondOrderTypes, only : first_order_state
  use FirstOrderState, only : radiation_psi, wave_gravity
  use FieldVelocity, only : incident_at
  use NearFieldDrift, only : waterline_extrapolate, KERNEL_DIFF, KERNEL_SUM
  use QTFPotential, only : potential_pair, psi_surface_gradient, psi_waterline, NPOT_TERMS
  use QTFFreeSurface, only : fs_inner_fields, fs_inner_steps, fs_sample_phi, fs_circle_modes, fs_auto_modes, &
                             fs_auto_evanescent, fs_pair_terms
  use PartialWaves, only : pw_sample_points, pw_evanescent_from_wvn
  use ControlSurfaceMesh, only : hull_radius
  use QTFDirect, only : direct_setup, direct_at_omega
  use SharpEdge, only : DEFAULT_EDGE_DEG
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
  ! Task 4. The 480-panel rounded cylinder of the 7a fixture qtf_direct_round (radius 10, draft 10, depth 100), with its
  ! first-order frequencies, the heading of the states and the water density
  character(len=*), parameter :: round_mesh = "test/data/qtf_direct_round/input/HullMesh.pnl"
  real(8), parameter :: rd_depth = 100d0, rd_omegas(3) = [0.70d0, 0.95d0, 1.10d0], rd_beta = 0.5235987755982988d0
  real(8), parameter :: rho_w = 1025d0
  ! direct_fs_equals_indirect_full: the free-surface part of T_dir (parts slice 3) against FS_inner + FS_outer of
  ! fs_pair_terms, over max |T_pot| (body and free-surface terms) of the pairs of each Omega. Measured 2.34e-3 at
  ! Omega = 0.25 and 2.83e-3 at 0.15 (max |T_FS| is 0.34 and 0.39 of max |T_pot|), well below the 1.9e-2 of the body
  ! forcing of this hull in 7a (qtf_direct_round, M_y at (0.70, 0.95)); 1.3x the larger, rounded up. With the opposite
  ! s_D: 0.67 and 0.78
  real(8), parameter :: fs_indirect_limit = 3.7d-3
  ! direct_fs_pair_order: a pair with Omega < 0 against the conjugate of the swapped pair, and the load without parts
  ! (the free-surface right side added to the body one) against the sum of the parts
  real(8), parameter :: pair_order_limit = 1d-12
  ! direct_fs_sum_equals_indirect: as fs_indirect_limit, for the sum pair (0.70, 0.70) at Omega = 1.40 (headings 0 and
  ! 30 deg) over max |T+_pot|, with the inner mesh built for Omega+ = 1.40. Measured 6.0e-3 (max |T_FS| is 0.77 of
  ! max |T+_pot|); 1.3x that, rounded up. On the inner mesh of 1.1 rad/s (100 points, too coarse for psi at Omega+) the
  ! gap is 1.54e-2, on that of 2.0 rad/s 4.5e-3: it falls with the inner mesh. With the opposite s_D: 1.55
  real(8), parameter :: fs_sum_indirect_limit = 7.9d-3

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
    new_unittest("fs_rhs_outer_contraction", test_outer_contraction), &
    new_unittest("fs_rfar_evanescent", test_rfar), &
    new_unittest("fs_memory_estimate", test_memory_estimate), &
    new_unittest("direct_fs_pair_order", test_direct_fs_pair_order), &
    new_unittest("direct_fs_equals_indirect_full", test_direct_fs_equals_indirect), &
    new_unittest("direct_fs_field_body_only", test_direct_fs_field_body_only), &
    new_unittest("direct_fs_kernel_default", test_direct_fs_kernel_default), &
    new_unittest("fs_green_matrix_large_kh", test_green_matrix_large_kh), &
    new_unittest("direct_fs_mixed_block", test_direct_fs_mixed_block), &
    new_unittest("direct_fs_sum_equals_indirect", test_direct_fs_sum_equals_indirect) &
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
!> above mcut are zero; cut_error is the bound q^(mcut+1) / ((mcut+1)(1 - q)), q = r_max / R, of the dropped tail
!> sum_{m > mcut} q^m / m, and lies above that sum (here q = 10 cos(pi/16) / 11 = 0.8916, mcut 93: 2.037e-6, the sum
!> 1.88e-6); with moments that are Inf and NaN above mcut (as fs_outer_moments gives there) and finite below, the right
!> side stays finite (no Inf x 0).
subroutine test_order_cut_finite(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, fm
  type(cs_points) :: pts
  type(direct_fs_data) :: fs
  complex(8), allocatable :: q(:, :), mom(:, :, :), rhs(:, :, :)
  real(8) :: centre(2), radius, inf, nan, qr, tail
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
  qr = fs%r_max / radius
  tail = 0d0
  do m = fs%mcut + 1, 100000
    tail = tail + qr**m / m
  end do
  write (info, '(a,es12.5,a,es12.5,a,i0)') "cut_error ", fs%cut_error, ", dropped tail sum ", tail, ", mcut ", fs%mcut
  print '(4x,a)', trim(info)
  call check(error, fs%mcut == 93 .and. abs(fs%cut_error / 2.0369d-6 - 1d0) <= 1d-3 .and. tail <= fs%cut_error, &
             "cut_error is not the bound of the dropped tail: "//trim(info))
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

!> The outer term of direct_fs_rhs against an independent reference (Task 3 review, Important 2). Two first-order states of
!> purely evanescent partial waves on the circle of the inner mesh (omega 1.0 and 0.7, so Omega = 0.3; one evanescent
!> mode, |m| <= 3 with amplitudes that differ between +m and -m), no incident waves: their Q_D decays like
!> e^{-(k1 + k2)(r - R)} outside R and has angular orders of both signs. mom from fs_outer_moments (no sliver), q = 0
!> inside; the right side at three hull rows (waterline, mid wall, bottom centre) against
!>    DIRECT_FS_SIGN / g int_R^{R + 60} int_0^{2 pi} G(x_i, xi) Q_D(xi) r dtheta dr
!> by Gauss-Legendre in r (8 points per metre) and the trapezoidal rule in theta (128 points, exact for the orders of
!> Q_D times those of G that are above 1e-10), with G from FINGREEN3D and Q_D from fs_forcing_diff of the fields of
!> pw_eval_surface. Gate 1e-5, the level of FINGREEN3D (ledger ruling F5). The contraction with a(-m) in place of a(m)
!> (the pairing of ruling F11 swapped) and with conj(a) (ZGEMM 'C') are printed and must fail by far more. Also the
!> cheap check: random finite mom, rhs against DIRECT_FS_SIGN / g sum_mn a(m, n, i) mom(m, n) by explicit loops, 1e-12.
subroutine test_outer_contraction(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: nth = 128, nint = 60, ngl = 8
  real(8), parameter :: w1s = 1d0, w2s = 0.7d0
  type(panel_mesh) :: mesh, fm
  type(cs_points) :: pts
  type(direct_fs_data) :: fs
  type(pw_set) :: pw(2)
  type(fs_point_field) :: t(2), zero
  complex(8), allocatable :: q(:, :), mom(:, :, :), rhs(:, :, :), qd(:, :)
  real(8), allocatable :: rr(:), wr(:)
  real(8) :: centre(2), radius, wv(2), x(3), xi(3), gx(ngl), gw(ngl), th, worst, worst_swap, worst_conj, rfar
  real(8) :: kk(2), kvr(0:3)
  complex(8) :: ref, sw, sc, e
  integer :: rows(3), k, m, n, i, ir, it, nr, col
  character(len=250) :: info
  external :: DISPERSION

  call fixture(32, 8, 6, 0, 0d0, mesh, pts, fm, centre, radius, error)
  if (allocated(error)) return
  call set_omega(fx_omega, fx_depth)
  call direct_fs_setup(mesh, pts, fm, centre, radius, fx_depth, fx_omega, .false., fs)
  ! The two states: one evanescent mode each, a(m, 1) K_m(k1 R) = c_m with c_{+m} /= c_{-m}
  wv = [w1s, w2s]
  do k = 1, 2
    call DISPERSION(kk, 2, wv(k), fx_depth)
    pw(k)%mmax = 3
    pw(k)%nev = 1
    pw(k)%k0 = kk(1)
    pw(k)%h = fx_depth
    pw(k)%radius = radius
    pw(k)%centre = centre
    allocate(pw(k)%kev(1), pw(k)%a(-3:3, 0:1))
    pw(k)%kev = kk(2)
    pw(k)%a = (0d0, 0d0)
    call pw_bessel_k(3, kk(2) * radius, kvr)
    do m = -3, 3
      pw(k)%a(m, 1) = cmplx(1d0 + 0.3d0 * m + 0.1d0 * k, 0.5d0 - 0.2d0 * m * k, 8) / kvr(abs(m))
    end do
  end do
  rfar = direct_fs_rfar(fs, radius + nint)
  allocate(mom(-fs%mg:fs%mg, 0:fs%nev, 1), q(pts%n, 1), rhs(mesh%nbase, 1, 1))
  call fs_outer_moments(w1s, pw(1), 0d0, 0d0, w2s, pw(2), 0d0, 0d0, fs%k0, fs%kev, fx_depth, centre, radius, fs%mg, &
                        fs%nev, fx_omega, fs%g, 0, rfar, mom)
  q = (0d0, 0d0)
  rhs = (0d0, 0d0)
  call direct_fs_rhs(fs, q, mom, rhs)

  ! Q_D on the polar grid of the reference, once for all rows
  call gauss(ngl, gx, gw)
  nr = nint * ngl
  allocate(rr(nr), wr(nr), qd(nr, nth))
  do i = 1, nint
    do k = 1, ngl
      rr((i - 1) * ngl + k) = radius + (i - 1) + 0.5d0 * (gx(k) + 1d0)
      wr((i - 1) * ngl + k) = 0.5d0 * gw(k)
    end do
  end do
  do ir = 1, nr
    do it = 1, nth
      th = 2d0 * pi * (it - 1) / nth
      xi = [centre(1) + rr(ir) * cos(th), centre(2) + rr(ir) * sin(th), 0d0]
      do k = 1, 2
        call pw_eval_surface(pw(k), xi, t(k)%phi, t(k)%grad, t(k)%phizz)
      end do
      qd(ir, it) = fs_forcing_diff(w1s, t(1), zero, w2s, t(2), zero, fs%g) * wr(ir) * rr(ir) * (2d0 * pi / nth)
    end do
  end do

  rows(1) = maxloc(mesh%xc(3, :), 1)
  rows(2) = minloc(abs(mesh%xc(3, :) + 5.625d0) + 0.01d0 * mesh%xc(1, :), 1)
  rows(3) = minloc(mesh%xc(3, :) + 0.01d0 * norm2(mesh%xc(1:2, :), dim=1), 1)
  worst = 0d0
  worst_swap = 0d0
  worst_conj = 0d0
  do k = 1, 3
    i = rows(k)
    x = mesh%xc(:, i)
    ref = (0d0, 0d0)
    do ir = 1, nr
      do it = 1, nth
        th = 2d0 * pi * (it - 1) / nth
        xi = [centre(1) + rr(ir) * cos(th), centre(2) + rr(ir) * sin(th), 0d0]
        ref = ref + green_point(x, xi) * qd(ir, it)
      end do
    end do
    ref = DIRECT_FS_SIGN / fs%g * ref
    sw = (0d0, 0d0)
    sc = (0d0, 0d0)
    do n = 0, fs%nev
      do m = -fs%mcut, fs%mcut
        sw = sw + fs%a(-m, n, i) * mom(m, n, 1)
        sc = sc + conjg(fs%a(m, n, i)) * mom(m, n, 1)
      end do
    end do
    sw = DIRECT_FS_SIGN / fs%g * sw
    sc = DIRECT_FS_SIGN / fs%g * sc
    worst = max(worst, abs(rhs(i, 1, 1) - ref) / abs(ref))
    worst_swap = max(worst_swap, abs(sw - ref) / abs(ref))
    worst_conj = max(worst_conj, abs(sc - ref) / abs(ref))
    print '(4x,a,3f8.3,a,2es12.4,a,es10.3)', "row at", x, ": reference ", ref, ", relative ", abs(rhs(i, 1, 1) - ref) / abs(ref)
  end do
  write (info, '(a,es10.3,a,es10.3,a,es10.3)') "outer contraction against quadrature of G Q_D: ", worst, &
    "; with a(-m): ", worst_swap, "; with conj(a): ", worst_conj
  print '(4x,a)', trim(info)
  call check(error, worst <= 1d-5 .and. worst_swap > 1d-2 .and. worst_conj > 1d-2, trim(info))
  if (allocated(error)) return

  ! The cheap check: random finite moments, explicit loops
  do n = 0, fs%nev
    do m = -fs%mg, fs%mg
      mom(m, n, 1) = pseudo(11 * m + 3 * n) / (1d0 + abs(m) + n)
    end do
  end do
  rhs = (0d0, 0d0)
  call direct_fs_rhs(fs, q, mom, rhs)
  worst = 0d0
  do col = 1, 1
    do i = 1, mesh%nbase
      e = (0d0, 0d0)
      do n = 0, fs%nev
        do m = -fs%mcut, fs%mcut
          e = e + fs%a(m, n, i) * mom(m, n, col)
        end do
      end do
      e = DIRECT_FS_SIGN / fs%g * e
      worst = max(worst, abs(rhs(i, col, 1) - e) / abs(e))
    end do
  end do
  write (info, '(a,es10.3)') "outer contraction against explicit loops, random moments: ", worst
  print '(4x,a)', trim(info)
  call check(error, worst <= 1d-12, trim(info))
end subroutine test_outer_contraction

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

!> The estimate of direct_fs_estimate, made before the frequency loop, against what direct_fs_setup allocates: the same
!> orders (the evanescent count of the setup is the same or lower, by the kn R cut) and the size of gfs and a to the byte
subroutine test_memory_estimate(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, fm
  type(cs_points) :: pts
  type(direct_fs_data) :: fs
  real(8) :: centre(2), radius, mb, mb_set
  integer :: mg, nev
  character(len=200) :: info

  call fixture(16, 4, 3, 0, 0d0, mesh, pts, fm, centre, radius, error)
  if (allocated(error)) return
  call set_omega(fx_omega, fx_depth)
  call direct_fs_setup(mesh, pts, fm, centre, radius, fx_depth, fx_omega, .false., fs)
  call direct_fs_estimate(mesh, pts%n, centre, radius, fx_depth, .false., mg, nev, mb)
  mb_set = 16d0 * (real(size(fs%gfs), 8) + real(size(fs%a), 8)) / 1d6
  write (info, '(a,2i5,a,2i5,a,es12.5,a,es12.5)') "estimate M, N_ev ", mg, nev, "; setup ", fs%mg, fs%nev, &
    "; MB ", mb, " against ", mb_set
  print '(4x,a)', trim(info)
  call check(error, mg == fs%mg .and. nev >= fs%nev, trim(info))
  if (allocated(error)) return
  ! with nev equal the sizes agree exactly; otherwise the estimate is the upper bound
  if (nev == fs%nev) then
    call check(error, abs(mb - mb_set) <= 1d-12 * mb_set, trim(info))
  else
    call check(error, mb >= mb_set, trim(info))
  end if
  if (allocated(error)) return
  ! the formula by hand: 16 B x nrow x (nfs + (2 mg + 1)(nev + 1))
  call check(error, abs(direct_fs_array_mb(10, 7, 3, 4) - 16d0 * 10d0 * (7d0 + 7d0 * 5d0) / 1d6) <= 1d-15, &
    "direct_fs_array_mb formula")
end subroutine test_memory_estimate

! ------------------------------------------------------------------------------------ the direct solve (Task 4)

!> Reads a hull mesh of HAMS (HullMesh.pnl, no symmetry plane) into a whole panel_mesh of one body, as make_cylinder
!> builds one: vertex-average centroid and area from the diagonals (set_panel, as HAMS), the flat-panel normal from the
!> diagonals, turned to point into the body (towards c_in, a point inside; the hull is convex)
subroutine read_hull(path, c_in, mesh, error)
  character(len=*), intent(in) :: path
  real(8), intent(in) :: c_in(3)
  type(panel_mesh), intent(out) :: mesh
  type(error_type), allocatable, intent(out) :: error
  logical :: ok
  call read_pnl_hull(path, c_in, mesh, ok)
  call check(error, ok, "cannot open "//path)
end subroutine read_hull

!> The geometry, the inner free surface of build_fs_inner (largest first-order frequency omega_max, default radius) with
!> its panels, and the partial-wave set-up of the driver (fs_setup of SecondOrderDriver): M and N_ev by the automatic
!> rules, the offsets of fs_inner_fields and the samples on the circle
subroutine fs_common(mesh, depth, omega_max, geom, pts, fm, centre, radius, npc, delta, circle, mmax, nev, error)
  type(panel_mesh), intent(in) :: mesh
  real(8), intent(in) :: depth, omega_max
  type(surface_geom), intent(out) :: geom
  type(cs_points), intent(out) :: pts, circle
  type(panel_mesh), intent(out) :: fm
  real(8), intent(out) :: centre(2), radius
  integer, intent(out) :: npc, mmax, nev
  real(8), allocatable, intent(out) :: delta(:)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: opts
  real(8) :: r_wl
  integer :: iw, nshrunk
  logical :: ok
  character(len=:), allocatable :: msg

  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "geometry: "//msg)
  if (allocated(error)) return
  call build_fs_inner(mesh, geom, opts, depth, 'fs_direct_no_input_dir', pts, centre, radius, ok, msg, &
                      omega_max=omega_max, np_circle=npc, fs_panels=fm)
  call check(error, ok, "inner mesh: "//msg)
  if (allocated(error)) return
  call set_omega(omega_max, depth)
  mmax = fs_auto_modes(0, WK, radius)
  r_wl = 0d0
  do iw = 1, geom%nwl
    r_wl = max(r_wl, norm2(geom%wl_x1(1:2, iw) - centre), norm2(geom%wl_x2(1:2, iw) - centre))
  end do
  nev = fs_auto_evanescent(0, depth, radius, max(r_wl, hull_radius(mesh, centre)))
  allocate (delta(pts%n))
  call fs_inner_steps(mesh, pts, 2d0 * opts%cs_fd_step, delta, nshrunk)
  call pw_sample_points(centre, radius, depth, mmax, nev, circle)
end subroutine fs_common

!> The first-order states of the QTF store at omegas (one heading each, betas in radians): the MacCamy-Fuchs wave of a
!> cylinder of radius 10 on the panels (a field with a scattered part; its exactness on the hull does not matter here, as
!> both methods take the same states), with rigid motions xis (6, n) for the body forcing, the waterline potential, and the
!> 6a store as the driver builds it: the fields at the inner points (fs_inner_fields), the incident fields there and the
!> partial waves of the scattered potential outside R (fs_circle_modes)
subroutine fs_states(mesh, geom, pts, delta, circle, centre, radius, depth, mmax, nev, omegas, betas, xis, nthread_fs, &
                     store, store_wl, fa, ia, pwa)
  type(panel_mesh), intent(in) :: mesh
  type(surface_geom), intent(in) :: geom
  type(cs_points), intent(in) :: pts, circle
  real(8), intent(in) :: delta(:), centre(2), radius, depth, omegas(:), betas(:)
  integer, intent(in) :: mmax, nev, nthread_fs
  complex(8), intent(in) :: xis(:, :)
  type(first_order_state), allocatable, intent(out) :: store(:, :)
  complex(8), allocatable, intent(out) :: store_wl(:, :, :)
  type(fs_point_field), allocatable, intent(out) :: fa(:, :, :), ia(:, :, :)
  type(pw_set), allocatable, intent(out) :: pwa(:, :)
  type(first_order_state) :: st(1)
  type(fs_point_field), allocatable :: f1(:, :)
  type(pw_set) :: pw1(1)
  type(so_options) :: opts
  complex(8), allocatable :: wl(:)
  complex(8) :: phi, grad(3)
  real(8) :: kev(nev), miss(1)
  integer :: i, j, ns

  ns = size(omegas)
  allocate (store(ns, 1), store_wl(geom%nwl, ns, 1), fa(pts%n, ns, 1), ia(pts%n, ns, 1), pwa(ns, 1), f1(pts%n, 1), &
            wl(geom%nwl))
  do i = 1, ns
    call set_omega(omegas(i), depth)
    NTHREAD = nthread_fs
    call maccamy_fuchs_state(mesh, 10d0, depth, WK, W1, 1d0, betas(i), wave_gravity(W1, WK, depth), rho_w, st(1))
    st(1)%xi(:, 1) = xis(:, i)
    call waterline_extrapolate(mesh, geom, st(1), wl)
    store(i, 1) = st(1)
    store_wl(:, i, 1) = wl
    call fs_inner_fields(mesh, pts, st, 2d0 * opts%cs_fd_step, .false., f1, delta)
    fa(:, i, 1) = f1(:, 1)
    do j = 1, pts%n
      call incident_at(st(1), pts%x(:, j), phi, grad)
      ia(j, i, 1)%phi = phi
      ia(j, i, 1)%grad = grad
      ia(j, i, 1)%phizz = WK**2 * phi
    end do
    call pw_evanescent_from_wvn(nev, kev)
    call fs_circle_modes(mesh, circle, st, .false., centre, radius, depth, WK, kev, mmax, nev, pw1, miss)
    pwa(i, 1) = pw1(1)
  end do
end subroutine fs_states

!> The radiation problem at omega on the ne panels of load_panels, as SolveRadiation runs it with IRSP = 0 and ISOLV = 1
!> (as radiation_at of test_qtf_direct): the stored Green function, AMAT factored, MXPOT; the wave of the solver globals
!> set to omega by SetWaveFromOmega, the depth that of the last set_omega
subroutine radiation_at(omega, ne)
  real(8), intent(in) :: omega
  integer, intent(in) :: ne

  ISOLV = 1
  call SetWaveFromOmega(omega)
  TNELEM = ne
  if (allocated(CGRN)) deallocate (CGRN)
  if (allocated(RKBN)) deallocate (RKBN)
  if (allocated(AMAT)) deallocate (AMAT)
  if (allocated(BRMAT)) deallocate (BRMAT)
  if (allocated(IPIV)) deallocate (IPIV)
  allocate (CGRN(ne, ne, NSYS, 4), RKBN(ne, ne, NSYS, 4), AMAT(ne, ne, NSYS), BRMAT(ne, 6, NSYS), IPIV(ne, NSYS))
  call CALGREEN
  call ASSB_LEFT(AMAT, IPIV, ne, NSYS)
  call ASSB_RBC(BRMAT, ne, NSYS)
  call RADIATION_SOLVER(AMAT, BRMAT, IPIV, MXPOT, ne, NSYS)
end subroutine radiation_at

!> psi_k of the radiation solve at omega (radiation_at first) at the inner points and its partial waves outside R, as
!> keep_fs_psi of SecondOrderDriver: one state per mode with psi_k as perturbation potential, dpsi_k/dn = n_k
subroutine psi_fs(mesh, pts, circle, centre, radius, depth, mmax, nev, omega, psi, psi_pts, pwp)
  type(panel_mesh), intent(in) :: mesh
  type(cs_points), intent(in) :: pts, circle
  real(8), intent(in) :: centre(2), radius, depth, omega
  integer, intent(in) :: mmax, nev
  complex(8), allocatable, intent(out) :: psi(:, :), psi_pts(:, :)
  type(pw_set), intent(out) :: pwp(6)
  type(first_order_state) :: ps(6)
  real(8) :: kev(nev), miss(6)
  integer :: k

  allocate (psi(mesh%n, 6), psi_pts(pts%n, 6))
  call radiation_psi(mesh, omega, MXPOT, psi)
  do k = 1, 6
    ps(k)%omega = omega
    ps(k)%k = WK
    ps(k)%nu = V
    ps(k)%h = H
    ps(k)%amp = 0d0
    ps(k)%beta = 0d0
    ps(k)%g = g_const
    ps(k)%rho = rho_w
    allocate (ps(k)%phi_p(mesh%n), ps(k)%dphi_p_dn(mesh%n))
    ps(k)%phi_p = psi(:, k)
    ps(k)%dphi_p_dn = mesh%nj(k, :)
  end do
  call fs_sample_phi(mesh, pts, ps, .false., psi_pts)
  call pw_evanescent_from_wvn(nev, kev)
  call fs_circle_modes(mesh, circle, ps, .false., centre, radius, depth, WK, kev, mmax, nev, pwp, miss)
end subroutine psi_fs

!> Rigid motions of three states (6, 3), different in every mode, so that every body-forcing term of T_pot is present
function rd_motions() result(x)
  complex(8) :: x(6, 3)
  x(:, 1) = [(0.3d0, 0d0), (-0.2d0, 0d0), (0.5d0, 0d0), (0.01d0, 0d0), (0.02d0, 0d0), (-0.015d0, 0d0)]
  x(:, 2) = [(0.1d0, 0.2d0), (0d0, 0.3d0), (-0.4d0, 0d0), (0d0, -0.02d0), (0.01d0, 0d0), (0.03d0, 0.01d0)]
  x(:, 3) = [(-0.2d0, 0.1d0), (0.15d0, -0.1d0), (0.3d0, 0.2d0), (0.02d0, 0.01d0), (-0.01d0, 0.02d0), (0d0, -0.02d0)]
end function rd_motions

!> A pair with Omega < 0 gives the conjugate of the swapped pair (header of QTFDirect: the pair is evaluated as the
!> swapped one, Q_D(2, 1) = conj Q_D(1, 2), and phi_D conjugated), on the 448-panel cylinder in depth 15, the states at
!> 1.1 rad/s (heading 0.6) and 0.8 rad/s (heading 0.2), Omega = 0.3: the free-surface part (parts slice 3) and the whole
!> load, to pair_order_limit of the largest. Without parts the free-surface right side is added to the body one and
!> solved with it: the load equals the sum of the three parts to the same limit. The free-surface part must be present
!> (above 1e-3 of the load), or the checks would hold for nothing. Measured: free-surface part 0 (the swapped pair is
!> solved from the same right side), load 2.8e-16, without parts 9.9e-16; free-surface share 0.17.
subroutine test_direct_fs_pair_order(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: om(2) = [1.1d0, 0.8d0], be(2) = [0.6d0, 0.2d0]
  type(panel_mesh) :: full, mesh, fm
  type(surface_geom) :: geom
  type(cs_points) :: pts, circle
  type(direct_fs_data) :: dfs
  type(first_order_state), allocatable :: store(:, :)
  complex(8), allocatable :: store_wl(:, :, :)
  type(fs_point_field), allocatable, target :: fa(:, :, :), ia(:, :, :)
  type(pw_set), allocatable, target :: pwa(:, :)
  real(8), allocatable :: delta(:)
  complex(8) :: t(6, 2), tn(6, 2), parts(6, 2, 3), xis(6, 2)
  integer :: pairs(4, 2), npc, mmax, nev
  real(8) :: centre(2), radius, omega, e_fs, e_all, e_sum, fs_share
  character(len=300) :: info

  call make_cylinder(10d0, 10d0, 32, 8, 6, full, .true.)
  call set_omega(fx_omega, fx_depth)
  call load_panels(full, 0, 0, mesh)
  call fs_common(mesh, fx_depth, fx_omax, geom, pts, fm, centre, radius, npc, delta, circle, mmax, nev, error)
  if (allocated(error)) return
  xis = rd_motions_two()
  call fs_states(mesh, geom, pts, delta, circle, centre, radius, fx_depth, mmax, nev, om, be, xis, 8, store, store_wl, &
                 fa, ia, pwa)
  omega = om(1) - om(2)
  ! Pair 1 has Omega < 0, pair 2 is it swapped
  pairs(:, 1) = [2, 1, 1, 1]
  pairs(:, 2) = [1, 1, 2, 1]
  call radiation_at(omega, mesh%nbase)
  NTHREAD = 8
  call direct_setup(mesh, DEFAULT_EDGE_DEG)
  call direct_fs_setup(mesh, pts, fm, centre, radius, fx_depth, omega, .false., dfs)
  dfs%f => fa
  dfs%inc => ia
  dfs%pw => pwa
  dfs%np = npc
  call direct_at_omega(mesh, geom, omega, store, store_wl, pairs, t, parts=parts, fs=dfs)
  call direct_at_omega(mesh, geom, omega, store, store_wl, pairs, tn, fs=dfs)
  NTHREAD = 1
  e_fs = maxval(abs(parts(:, 1, 3) - conjg(parts(:, 2, 3)))) / maxval(abs(parts(:, 2, 3)))
  e_all = maxval(abs(t(:, 1) - conjg(t(:, 2)))) / maxval(abs(t))
  e_sum = maxval(abs(tn - sum(parts, dim=3))) / maxval(abs(tn))
  fs_share = maxval(abs(parts(:, :, 3))) / maxval(abs(t))
  write (info, '(a,es10.2,a,es10.2,a,es10.2,a,es10.2)') "Omega < 0 against the swapped pair: free-surface part ", e_fs, &
    ", load ", e_all, "; load without parts against the sum of the parts ", e_sum, "; free-surface share ", fs_share
  print '(4x,a)', trim(info)
  call check(error, e_fs <= pair_order_limit .and. e_all <= pair_order_limit .and. e_sum <= pair_order_limit .and. &
             fs_share > 1d-3, trim(info))
end subroutine test_direct_fs_pair_order

!> With field points and the free-surface forcing (the 7b and 7c-1 merge, ledger ruling C15): direct_at_omega solves the
!> free-surface columns apart, so phi_field is the phi_D of the body forcing only. The setup of test_direct_fs_pair_order
!> (448-panel cylinder in depth 15, Omega = 0.3, a pair and its swap), three points off the hull. phi_field with fs must
!> equal phi_field without fs bit for bit, and the load with fs and field points must equal the load with fs and the
!> three parts (the same apart solve) bit for bit, with a free-surface share above 1e-3 so that the check is not empty.
subroutine test_direct_fs_field_body_only(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: om(2) = [1.1d0, 0.8d0], be(2) = [0.6d0, 0.2d0]
  real(8), parameter :: xp(3, 3) = reshape([14d0, 0d0, -2d0, 0d0, 13d0, -5d0, -12d0, 3d0, -0.5d0], [3, 3])
  type(panel_mesh) :: full, mesh, fm
  type(surface_geom) :: geom
  type(cs_points) :: pts, circle
  type(direct_fs_data) :: dfs
  type(first_order_state), allocatable :: store(:, :)
  complex(8), allocatable :: store_wl(:, :, :)
  type(fs_point_field), allocatable, target :: fa(:, :, :), ia(:, :, :)
  type(pw_set), allocatable, target :: pwa(:, :)
  real(8), allocatable :: delta(:)
  complex(8) :: t0(6, 2), t1(6, 2), t2(6, 2), parts(6, 2, 3), xis(6, 2), f0(3, 2), f1(3, 2)
  integer :: pairs(4, 2), npc, mmax, nev
  real(8) :: centre(2), radius, omega, fs_share
  character(len=300) :: info

  call make_cylinder(10d0, 10d0, 32, 8, 6, full, .true.)
  call set_omega(fx_omega, fx_depth)
  call load_panels(full, 0, 0, mesh)
  call fs_common(mesh, fx_depth, fx_omax, geom, pts, fm, centre, radius, npc, delta, circle, mmax, nev, error)
  if (allocated(error)) return
  xis = rd_motions_two()
  call fs_states(mesh, geom, pts, delta, circle, centre, radius, fx_depth, mmax, nev, om, be, xis, 8, store, store_wl, &
                 fa, ia, pwa)
  omega = om(1) - om(2)
  pairs(:, 1) = [2, 1, 1, 1]
  pairs(:, 2) = [1, 1, 2, 1]
  call radiation_at(omega, mesh%nbase)
  NTHREAD = 8
  call direct_setup(mesh, DEFAULT_EDGE_DEG)
  call direct_fs_setup(mesh, pts, fm, centre, radius, fx_depth, omega, .false., dfs)
  dfs%f => fa
  dfs%inc => ia
  dfs%pw => pwa
  dfs%np = npc
  call direct_at_omega(mesh, geom, omega, store, store_wl, pairs, t0, field_points=xp, phi_field=f0)
  call direct_at_omega(mesh, geom, omega, store, store_wl, pairs, t1, fs=dfs, field_points=xp, phi_field=f1)
  call direct_at_omega(mesh, geom, omega, store, store_wl, pairs, t2, parts=parts, fs=dfs)
  NTHREAD = 1
  fs_share = maxval(abs(parts(:, :, 3))) / maxval(abs(t2))
  write (info, '(a,es10.2,a,es10.2,a,es10.2)') "phi_field with fs minus without ", maxval(abs(f1 - f0)), &
    "; load with fs and field points minus with fs and parts ", maxval(abs(t1 - t2)), "; free-surface share ", fs_share
  print '(4x,a)', trim(info)
  call check(error, all(f1 == f0) .and. all(t1 == t2) .and. maxval(abs(f0)) > 0d0 .and. fs_share > 1d-3, trim(info))
end subroutine test_direct_fs_field_body_only

!> The difference kernel stays the default of the free-surface block, bit for bit (Phase 7c-2 Task 3): the setup of
!> test_direct_fs_pair_order (448-panel cylinder in depth 15, Omega = 0.3, the pair with Omega < 0 and its swap), the
!> free-surface part of the load (parts slice 3) of both columns against the values printed from the code before the
!> kernel argument (47a088c, 17 digits, round-trip). Printed again at the merge of main's single g (DISPERSION from
!> g = 9.807 to Const_mod's 9.80665): they moved by at most 1.2e-4 of the largest.
subroutine test_direct_fs_kernel_default(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: om(2) = [1.1d0, 0.8d0], be(2) = [0.6d0, 0.2d0]
  type(panel_mesh) :: full, mesh, fm
  type(surface_geom) :: geom
  type(cs_points) :: pts, circle
  type(direct_fs_data) :: dfs
  type(first_order_state), allocatable :: store(:, :)
  complex(8), allocatable :: store_wl(:, :, :)
  type(fs_point_field), allocatable, target :: fa(:, :, :), ia(:, :, :)
  type(pw_set), allocatable, target :: pwa(:, :)
  real(8), allocatable :: delta(:)
  complex(8), parameter :: ref(6, 2) = reshape([ &
    (-9.27995863068147810d+003, -1.58886198061486903d+003), (-5.53121799474604722d+003, -5.60510515086619216d+003), &
    (-9.64925858196820082d+003, 2.44823219856080686d+003), (-1.59832169016786684d+004, -1.54310962170309213d+004), &
    (2.72170912635856585d+004, 1.61555653903165421d+003), (-5.10002847458964345d-013, -2.20078491958482441d-012), &
    (-9.27995863068147810d+003, 1.58886198061486903d+003), (-5.53121799474604722d+003, 5.60510515086619216d+003), &
    (-9.64925858196820082d+003, -2.44823219856080686d+003), (-1.59832169016786684d+004, 1.54310962170309213d+004), &
    (2.72170912635856585d+004, -1.61555653903165421d+003), (-5.10002847458964345d-013, 2.20078491958482441d-012)], [6, 2])
  complex(8) :: t(6, 2), parts(6, 2, 3), xis(6, 2), parts_k(6, 2, 3), t_k(6, 2)
  integer :: pairs(4, 2), npc, mmax, nev, p, k
  real(8) :: centre(2), radius, omega

  call make_cylinder(10d0, 10d0, 32, 8, 6, full, .true.)
  call set_omega(fx_omega, fx_depth)
  call load_panels(full, 0, 0, mesh)
  call fs_common(mesh, fx_depth, fx_omax, geom, pts, fm, centre, radius, npc, delta, circle, mmax, nev, error)
  if (allocated(error)) return
  xis = rd_motions_two()
  call fs_states(mesh, geom, pts, delta, circle, centre, radius, fx_depth, mmax, nev, om, be, xis, 8, store, store_wl, &
                 fa, ia, pwa)
  omega = om(1) - om(2)
  pairs(:, 1) = [2, 1, 1, 1]
  pairs(:, 2) = [1, 1, 2, 1]
  call radiation_at(omega, mesh%nbase)
  NTHREAD = 8
  call direct_setup(mesh, DEFAULT_EDGE_DEG)
  call direct_fs_setup(mesh, pts, fm, centre, radius, fx_depth, omega, .false., dfs)
  dfs%f => fa
  dfs%inc => ia
  dfs%pw => pwa
  dfs%np = npc
  call direct_at_omega(mesh, geom, omega, store, store_wl, pairs, t, parts=parts, fs=dfs)
  call direct_at_omega(mesh, geom, omega, store, store_wl, pairs, t_k, parts=parts_k, fs=dfs, kernels=[KERNEL_DIFF, &
                       KERNEL_DIFF])
  NTHREAD = 1
  if (.not. all(parts(:, :, 3) == ref)) then
    do p = 1, 2
      do k = 1, 6
        print '(4x,a,i0,a,i0,a,2es25.17e3)', "pair ", p, " mode ", k, ": ", parts(k, p, 3)
      end do
    end do
  end if
  call check(error, all(parts(:, :, 3) == ref), "the free-surface part of the load is not bit-identical to 47a088c")
  if (allocated(error)) return
  call check(error, all(parts_k == parts) .and. all(t_k == t), "kernels = KERNEL_DIFF is not the default bit for bit")
end subroutine test_direct_fs_kernel_default

!> (1 / w_j) int_j G dS over inner panel p of fm by n x n Gauss-Legendre (as green_quadrature) of a reference G at the
!> solver globals: deep = .false., the eigenfunction series of finite depth summed directly (direct_sum_green, checked
!> against 30-digit mpmath to 3e-11) plus the 1/r + 1/r1 it leaves out; deep = .true., INFGREEN3D (flag 0)
complex(8) function green_quadrature_ref(fm, p, x, n, deep) result(s)
  type(panel_mesh), intent(in) :: fm
  integer, intent(in) :: p, n
  real(8), intent(in) :: x(3)
  logical, intent(in) :: deep
  real(8) :: gx(n), gw(n), xv(3, 4), xq(3), du(3), ds(3), c(3), u, t, rh
  complex(8) :: grn(4), gd(3)
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
      if (deep) then
        call INFGREEN3D(xq(1), x(1), xq(2), x(2), xq(3), x(3), V, grn, 0)
      else
        rh = hypot(xq(1) - x(1), xq(2) - x(2))
        call direct_sum_green(rh, xq(3), x(3), V, WVN(1), H, gd)
        grn(1) = gd(1) + 1d0 / norm2(xq - x) + 1d0 / sqrt(rh**2 + (xq(3) + x(3))**2)
      end if
      s = s + grn(1) * gw(iu) * gw(iv) * norm2(c)
    end do
  end do
end function green_quadrature_ref

!> int over inner panel p of fm of the rest G - 1/r - 1/r1 by n x n Gauss-Legendre (as green_quadrature): of the
!> eigenfunction series summed directly (series = .true., direct_sum_green), or of FINGREEN3D (flag 1, which leaves out
!> 1/r) less 1/r1, as panel_green of QTFDirectFS takes it
complex(8) function rest_quad(fm, p, x, n, series) result(s)
  type(panel_mesh), intent(in) :: fm
  integer, intent(in) :: p, n
  real(8), intent(in) :: x(3)
  logical, intent(in) :: series
  real(8) :: gx(n), gw(n), xv(3, 4), xq(3), du(3), ds(3), c(3), u, t, rh
  complex(8) :: grn(4), gd(3), val
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
      rh = hypot(xq(1) - x(1), xq(2) - x(2))
      if (series) then
        call direct_sum_green(rh, xq(3), x(3), V, WVN(1), H, gd)
        val = gd(1)
      else
        call FINGREEN3D(xq(1), x(1), xq(2), x(2), xq(3), x(3), V, WVN, NK, H, grn, 1)
        val = grn(1) - 1d0 / sqrt(rh**2 + (xq(3) + x(3))**2)
      end if
      s = s + val * gw(iu) * gw(iv) * norm2(c)
    end do
  end do
end function rest_quad

!> The inner Green matrix at large K h (Review Focus 3 of the 7c-2 plan: the sum frequencies reach K h far beyond those
!> of the difference): direct_fs_setup at Omega = 2.2 rad/s in h = 100 m (K h = 2.2^2 / 9.81 x 100 = 49, just below
!> the nu h = 50 limit of the series choice of FINGREEN3D) on the 448-panel cylinder and its inner mesh, gfs(i, j) at 5
!> (row, panel) pairs away from the waterline (as fs_green_matrix_matches) against (1 / w_j) int_j G dS by 20 x 20
!> Gauss-Legendre of the eigenfunction series summed directly (direct_sum_green). Measured 1.1e-7 to 7.8e-7 at three
!> pairs, but 2.3e-5 and 1.40e-4 at the two pairs of the smallest R/h (0.05 and 0.04). That is FINGREEN3D, not the
!> panel rule (fix round 1 of Task 3, the reviewer's question whether the Gauss order of fs_remainder_order, from d / a
!> alone, misses the oscillation of the rest on the scale 1 / K = 2 m, K a = 1.24 here):
!>   - the rule's own order (3 to 6 here) on the rest of the series is within 1.2e-7 of 40 x 40 points at all five
!>     pairs (1.3e-12 at the worst pair), so the order resolves the rest at this K a; gated at 1.6e-7 (1.3x, rounded up);
!>   - FINGREEN3D is what does not converge: on the worst pair its rest at the rule's order differs from 40 x 40 points by
!>     7.8e-5 and 20 x 20 from 40 x 40 by 4.0e-6 (it is not smooth on the panel), and at 40 x 40 it is 2.2e-4 off the
!>     series (1.9e-7 to 7.8e-7 at the three good pairs, 2.3e-5 at the other); the plain 20 x 20 Gauss of FINGREEN3D
!>     (flag 0) that the review proposed as the reference differs from gfs by 7.4e-5 there for the same reason;
!>   - at the centroids FINGREEN3D is off the series by 2.5e-5 and 1.3e-5 (1.6e-7 to 6.8e-7 at the other three), and at
!>     the R/h = 0.04 centroid its error grows with nu h: 3.4e-7, 8.5e-6 and 3.4e-5 at nu h = 1.6, 26 and 59.
!> All of these are printed per pair. The gate on gfs is 1.3x the worst, rounded up, 2e-4 (not the 1e-6 of
!> fs_green_matrix_matches; rulings T3-a and T3-c of the 7c-2 ledger).
!> INFGREEN3D, the deep-water Green function, is printed beside it but not held to 1e-6: it is
!> the global approximation of Wu et al. (2017, Eur. J. Mech. B 65, 54-64), about 7e-4 off the exact function (measured
!> 6e-4 to 9e-4 against FINGREEN3D at K h = 11 to 47, where the depth changes G by O(e^{-2 K h}) only [derived: the
!> bottom enters through cosh(k (z + h)) cosh(k (zeta + h)) / cosh^2(k h) = e^{k (z + zeta)} (1 + O(e^{-2 k (h + z)}))
!> per wave number, and nu = K tanh(K h) = K (1 + O(e^{-2 K h})) at the pole], so it is gated at 1e-2: the deep limit.
subroutine test_green_matrix_large_kh(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: omega = 2.2d0, depth = 100d0
  type(panel_mesh) :: mesh, fm
  type(cs_points) :: pts
  type(direct_fs_data) :: fs
  real(8) :: centre(2), radius, err, worst, err_deep, worst_deep, rh, pr, dd, gw_ij, e_rule, e_fin, e_fin_ser, worst_rule
  complex(8) :: ref, ref_deep, gp(4), gd(3), gc
  integer :: rows(5), cols(5), k, i, j, kk, nr
  character(len=250) :: info

  call fixture(32, 8, 6, 0, 0d0, mesh, pts, fm, centre, radius, error)
  if (allocated(error)) return
  call set_omega(omega, depth)
  print '(4x,a,f8.2)', "K h ", WVN(1) * depth
  call direct_fs_setup(mesh, pts, fm, centre, radius, depth, omega, .false., fs)
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
  worst = 0d0
  worst_deep = 0d0
  worst_rule = 0d0
  do k = 1, 5
    i = rows(k)
    j = cols(k)
    ref = green_quadrature_ref(fm, j, mesh%xc(:, i), 20, .false.) / pts%w(j)
    ref_deep = green_quadrature_ref(fm, j, mesh%xc(:, i), 20, .true.) / pts%w(j)
    err = abs(fs%gfs(i, j) - ref) / abs(ref)
    err_deep = abs(fs%gfs(i, j) - ref_deep) / abs(ref_deep)
    print '(4x,a,i0,a,i0,a,2es14.6,a,es10.2,a,es10.2)', "row ", i, ", inner panel ", j, ": gfs ", fs%gfs(i, j), &
      ", against the direct series ", err, ", against INFGREEN3D ", err_deep
    worst = max(worst, err)
    worst_deep = max(worst_deep, err_deep)
    ! The panel rule: its Gauss order (fs_remainder_order of d and the panel radius) on the rest G - 1/r - 1/r1 of the
    ! series, against 40 x 40 points; and the same with the rest of FINGREEN3D, as panel_green takes it
    pr = 0d0
    do kk = 1, fm%ncn(j)
      pr = max(pr, norm2(fm%xnode(:, fm%ncon(kk, j)) - pts%x(:, j)))
    end do
    dd = norm2(mesh%xc(:, i) - pts%x(:, j))
    nr = fs_remainder_order(dd, pr)
    gw_ij = abs(fs%gfs(i, j) * pts%w(j))
    e_rule = abs(rest_quad(fm, j, mesh%xc(:, i), nr, .true.) - rest_quad(fm, j, mesh%xc(:, i), 40, .true.)) / gw_ij
    e_fin = abs(rest_quad(fm, j, mesh%xc(:, i), nr, .false.) - rest_quad(fm, j, mesh%xc(:, i), 40, .false.)) / gw_ij
    e_fin_ser = abs(rest_quad(fm, j, mesh%xc(:, i), 40, .false.) - rest_quad(fm, j, mesh%xc(:, i), 40, .true.)) / gw_ij
    print '(6x,a,f6.2,a,f5.2,a,i0,a,es10.2,a,es10.2,a,es10.2)', "d / a ", dd / pr, ", K a ", WVN(1) * pr, ", order ", nr, &
      ": rule on the series ", e_rule, ", rule on FINGREEN3D ", e_fin, "; FINGREEN3D against the series (40 x 40) ", e_fin_ser
    worst_rule = max(worst_rule, e_rule)
    ! FINGREEN3D itself at the centroids against the series: the error of the Green function, not of the panel rule
    call FINGREEN3D(fm%xc(1, j), mesh%xc(1, i), fm%xc(2, j), mesh%xc(2, i), fm%xc(3, j), mesh%xc(3, i), V, WVN, NK, H, &
                    gp, 0)
    rh = hypot(fm%xc(1, j) - mesh%xc(1, i), fm%xc(2, j) - mesh%xc(2, i))
    call direct_sum_green(rh, fm%xc(3, j), mesh%xc(3, i), V, WVN(1), H, gd)
    gc = gd(1) + 1d0 / norm2(fm%xc(:, j) - mesh%xc(:, i)) + 1d0 / sqrt(rh**2 + (fm%xc(3, j) + mesh%xc(3, i))**2)
    print '(6x,a,f7.4,a,es10.2)', "R/h ", rh / depth, ": FINGREEN3D at the centroid against the series ", &
      abs(gp(1) - gc) / abs(gc)
  end do
  write (info, '(a,es10.3,a,es10.3,a,es10.3)') "gfs at K h = 49 against 20 x 20 Gauss of the direct eigenfunction "// &
    "series, worst relative ", worst, "; against INFGREEN3D ", worst_deep, "; the panel rule on the series ", worst_rule
  print '(4x,a)', trim(info)
  call check(error, worst <= 2d-4 .and. worst_deep <= 1d-2 .and. worst_rule <= 1.6d-7, trim(info))
end subroutine test_green_matrix_large_kh

!> A block that mixes sum and difference columns with the free-surface forcing (the stop of Task 2 removed, Phase 7c-2
!> Task 3): on the 448-panel cylinder in depth 15 with the states at 0.3, 0.5 and 1.1 rad/s (headings 0.6, 0.2, -0.4),
!> Omega = 0.8: the sum pair (0.3, 0.5), given with w1 < w2, the difference pair (1.1, 0.3) and its swap (0.3, 1.1),
!> and the sum pair (0.5, 0.3). Each column of the mixed block equals the same pair solved alone (a block of its own
!> kernel) to 1e-13 of the largest load (tdir and the three parts); the two orders of the sum pair agree to 1e-13
!> (direct_fs_block orders a sum pair by frequency, no conjugate), the two orders of the difference pair are conjugate;
!> and the free-surface part of the sum columns is present (above 1e-3 of the load), so the checks do not hold for
!> nothing.
subroutine test_direct_fs_mixed_block(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: om(3) = [0.3d0, 0.5d0, 1.1d0], be(3) = [0.6d0, 0.2d0, -0.4d0], omega = 0.8d0
  type(panel_mesh) :: full, mesh, fm
  type(surface_geom) :: geom
  type(cs_points) :: pts, circle
  type(direct_fs_data) :: dfs
  type(first_order_state), allocatable :: store(:, :)
  complex(8), allocatable :: store_wl(:, :, :)
  type(fs_point_field), allocatable, target :: fa(:, :, :), ia(:, :, :)
  type(pw_set), allocatable, target :: pwa(:, :)
  real(8), allocatable :: delta(:)
  complex(8) :: t(6, 4), parts(6, 4, 3), ts(6, 2), ps(6, 2, 3), td(6, 2), pd(6, 2, 3), xis(6, 3)
  integer :: pairs(4, 4), npc, mmax, nev
  real(8) :: centre(2), radius, big, e_alone, e_sum_swap, e_diff_swap, fs_share
  character(len=400) :: info

  call make_cylinder(10d0, 10d0, 32, 8, 6, full, .true.)
  call set_omega(fx_omega, fx_depth)
  call load_panels(full, 0, 0, mesh)
  call fs_common(mesh, fx_depth, fx_omax, geom, pts, fm, centre, radius, npc, delta, circle, mmax, nev, error)
  if (allocated(error)) return
  xis = rd_motions()
  call fs_states(mesh, geom, pts, delta, circle, centre, radius, fx_depth, mmax, nev, om, be, xis, 8, store, store_wl, &
                 fa, ia, pwa)
  pairs(:, 1) = [1, 1, 2, 1]
  pairs(:, 2) = [3, 1, 1, 1]
  pairs(:, 3) = [1, 1, 3, 1]
  pairs(:, 4) = [2, 1, 1, 1]
  call radiation_at(omega, mesh%nbase)
  NTHREAD = 8
  call direct_setup(mesh, DEFAULT_EDGE_DEG)
  call direct_fs_setup(mesh, pts, fm, centre, radius, fx_depth, omega, .false., dfs)
  dfs%f => fa
  dfs%inc => ia
  dfs%pw => pwa
  dfs%np = npc
  call direct_at_omega(mesh, geom, omega, store, store_wl, pairs, t, parts=parts, fs=dfs, &
                       kernels=[KERNEL_SUM, KERNEL_DIFF, KERNEL_DIFF, KERNEL_SUM])
  call direct_at_omega(mesh, geom, omega, store, store_wl, pairs(:, [1, 4]), ts, parts=ps, fs=dfs, &
                       kernels=[KERNEL_SUM, KERNEL_SUM])
  call direct_at_omega(mesh, geom, omega, store, store_wl, pairs(:, [2, 3]), td, parts=pd, fs=dfs)
  NTHREAD = 1
  big = maxval(abs(t))
  e_alone = max(maxval(abs(t(:, [1, 4]) - ts)), maxval(abs(t(:, [2, 3]) - td)), maxval(abs(parts(:, [1, 4], :) - ps)), &
                maxval(abs(parts(:, [2, 3], :) - pd))) / big
  e_sum_swap = max(maxval(abs(t(:, 1) - t(:, 4))), maxval(abs(parts(:, 1, :) - parts(:, 4, :)))) / big
  e_diff_swap = maxval(abs(t(:, 2) - conjg(t(:, 3)))) / big
  fs_share = maxval(abs(parts(:, [1, 4], 3))) / maxval(abs(t(:, [1, 4])))
  write (info, '(a,es10.2,a,es10.2,a,es10.2,a,es10.2)') "mixed block against each kernel alone ", e_alone, &
    "; sum (1, 2) against (2, 1) ", e_sum_swap, "; difference against conj of its swap ", e_diff_swap, &
    "; free-surface share of the sum columns ", fs_share
  print '(4x,a)', trim(info)
  call check(error, e_alone <= 1d-13 .and. e_sum_swap <= 1d-13 .and. e_diff_swap <= 1d-12 .and. fs_share > 1d-3, &
             trim(info))
end subroutine test_direct_fs_mixed_block

!> The free-surface part of the direct load of the sum against that of the indirect method (Phase 7c-2 Task 3): the
!> fixture of direct_fs_equals_indirect_full (the 480-panel rounded cylinder in depth 100,
!> but the inner mesh of build_fs_inner for Omega+ = 1.40), the states (0.70, heading 0) and (0.70, heading 30 deg), the
!> sum pair at Omega = 1.40 in both orders. By Green's
!> theorem with psi_k at Omega (QTFPotential header, the sum),
!>    i Omega rho int_S phi_D,FS n_k dS = s_FS (i Omega rho / g) int_F psi_k Q_S dS = FS_inner + FS_outer,
!> so parts slice 3 of direct_at_omega (KERNEL_SUM) equals FS_inner + FS_outer of fs_pair_terms (KERNEL_SUM) to the mesh
!> error, measured over max |T+_pot| (the body terms of potential_pair with KERNEL_SUM plus the free-surface terms).
!> With the opposite s_D the slice changes sign: the gap must then be of order 1 (above 10x the limit), so s_D carries
!> over to the sum.
subroutine test_direct_fs_sum_equals_indirect(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: om(2) = [0.70d0, 0.70d0], be(2) = [0d0, 0.5235987755982988d0], omega = 1.40d0
  type(panel_mesh) :: full, mesh, fm
  type(surface_geom) :: geom
  type(cs_points) :: pts, circle
  type(direct_fs_data) :: dfs
  type(first_order_state), allocatable :: store(:, :)
  complex(8), allocatable :: store_wl(:, :, :), psi(:, :), psi_pts(:, :), gpsi(:, :, :), pwl(:, :)
  type(fs_point_field), allocatable, target :: fa(:, :, :), ia(:, :, :)
  type(pw_set), allocatable, target :: pwa(:, :)
  type(pw_set) :: pwp(6)
  real(8), allocatable :: delta(:)
  complex(8) :: t(6, 2), parts(6, 2, 3), terms(6, NPOT_TERMS), tin(6), tout(6), tfs(6, 2), tpot(6, 2)
  integer :: pairs(4, 2), npc, mmax, nev, p
  real(8) :: centre(2), radius, gap, opp, rel_fs, body_gap
  integer(8) :: c0, c1, rate
  character(len=400) :: info

  call system_clock(c0, rate)
  call read_hull(round_mesh, [0d0, 0d0, -5d0], full, error)
  if (allocated(error)) return
  call set_omega(rd_omegas(3), rd_depth)
  call load_panels(full, 0, 0, mesh)
  call check(error, mesh%n == 480, "the fixture is not the 480-panel rounded cylinder")
  if (allocated(error)) return
  ! The inner mesh for Omega+ (panels of a tenth of the wavelength of psi at Omega+), not for the first-order 1.1 rad/s
  call fs_common(mesh, rd_depth, omega, geom, pts, fm, centre, radius, npc, delta, circle, mmax, nev, error)
  if (allocated(error)) return
  print '(4x,a,i0,a,i0,a,i0)', "inner points ", pts%n, ", M ", mmax, ", N_ev ", nev
  call fs_states(mesh, geom, pts, delta, circle, centre, radius, rd_depth, mmax, nev, om, be, rd_motions_two(), 8, &
                 store, store_wl, fa, ia, pwa)
  allocate (gpsi(3, mesh%n, 6), pwl(geom%nwl, 6))
  pairs(:, 1) = [1, 1, 2, 1]
  pairs(:, 2) = [2, 1, 1, 1]
  call radiation_at(omega, mesh%nbase)
  NTHREAD = 8
  call psi_fs(mesh, pts, circle, centre, radius, rd_depth, mmax, nev, omega, psi, psi_pts, pwp)
  call psi_surface_gradient(mesh, geom, psi, gpsi)
  call psi_waterline(mesh, geom, psi, gpsi, pwl)
  do p = 1, 2
    associate (s1 => store(pairs(1, p), 1), s2 => store(pairs(3, p), 1), i1 => pairs(1, p), i2 => pairs(3, p))
      call potential_pair(mesh, geom, s1, s2, store_wl(:, i1, 1), store_wl(:, i2, 1), psi, gpsi, pwl, omega, s1%rho, &
                          wave_gravity(s1%omega, s1%k, s1%h), terms, kernel=KERNEL_SUM)
      call fs_pair_terms(s1%omega, fa(:, i1, 1), ia(:, i1, 1), pwa(i1, 1), s1%beta, s1%k, &
                         s2%omega, fa(:, i2, 1), ia(:, i2, 1), pwa(i2, 1), s2%beta, s2%k, &
                         psi_pts, pwp, pts%w, npc, s1%g, s1%rho, tin, tout, kernel=KERNEL_SUM)
    end associate
    tfs(:, p) = tin + tout
    tpot(:, p) = sum(terms, dim=2) + tfs(:, p)
  end do
  call direct_setup(mesh, DEFAULT_EDGE_DEG)
  call direct_fs_setup(mesh, pts, fm, centre, radius, rd_depth, omega, .false., dfs)
  dfs%f => fa
  dfs%inc => ia
  dfs%pw => pwa
  dfs%np = npc
  call direct_at_omega(mesh, geom, omega, store, store_wl, pairs, t, parts=parts, fs=dfs, kernels=[KERNEL_SUM, KERNEL_SUM])
  NTHREAD = 1
  gap = maxval(abs(parts(:, :, 3) - tfs)) / maxval(abs(tpot))
  opp = maxval(abs(-parts(:, :, 3) - tfs)) / maxval(abs(tpot))
  rel_fs = maxval(abs(tfs)) / maxval(abs(tpot))
  ! The body-forced part alone (direct slices 1 + 2 against the body terms), the mesh error of Task 2 at this Omega
  body_gap = maxval(abs(parts(:, :, 1) + parts(:, :, 2) - (tpot - tfs))) / maxval(abs(tpot))
  do p = 1, 2
    print '(6x,a,i0,a,6es11.3)', "pair ", p, " slice 3 ", abs(parts(:, p, 3))
    print '(6x,a,i0,a,6es11.3)', "pair ", p, " T_FS    ", abs(tfs(:, p))
  end do
  call system_clock(c1)
  write (info, '(a,es10.3,a,es9.2,a,es10.3,a,es10.3,a,es10.3,a,i0,a,i0,a,i0,a,f7.1,a)') &
    "sum at Omega 1.40: |slice 3 - (FS_inner + FS_outer)| / max |T+_pot| ", gap, " (limit ", fs_sum_indirect_limit, &
    "); opposite s_D ", opp, "; max |T_FS| / max |T+_pot| ", rel_fs, "; body part gap ", body_gap, "; mg ", dfs%mg, &
    ", nev ", dfs%nev, ", mcut ", dfs%mcut, "; ", real(c1 - c0, 8) / real(rate, 8), " s"
  print '(4x,a)', trim(info)
  call check(error, gap <= fs_sum_indirect_limit .and. gap <= 1d-2 .and. opp > 10d0 * fs_sum_indirect_limit, trim(info))
end subroutine test_direct_fs_sum_equals_indirect

!> Rigid motions of two states (6, 2), the first two of rd_motions
function rd_motions_two() result(x)
  complex(8) :: x(6, 2), x3(6, 3)
  x3 = rd_motions()
  x = x3(:, 1:2)
end function rd_motions_two

!> The free-surface part of the direct load against that of the indirect method (the sign check of s_D, ledger Task 4):
!> on the 480-panel rounded cylinder of the 7a fixture qtf_direct_round (radius 10, draft 10, depth 100) with the inner
!> free surface of build_fs_inner at 1.1 rad/s, the states at 0.70, 0.95, 1.10 rad/s (heading 30 deg) and Omega = 0.25
!> (0.95, 0.70) and 0.15 (1.10, 0.95), both orders of each pair. By Green's theorem with psi_k (the radiation potential at
!> Omega, which solves the same equation as phi_D with n_k on the hull),
!>    i Omega rho int_S phi_D,FS n_k dS = s_FS (i Omega rho / g) int_F psi_k Q_D dS = FS_inner + FS_outer,
!> so parts slice 3 of direct_at_omega equals the terms of fs_pair_terms to the mesh error. Measured over max |T_pot| of
!> the pairs of each Omega (body terms of potential_pair plus the free-surface terms). With the opposite s_D the slice
!> changes sign, and the difference is twice the free-surface part: it must fail by order 1.
subroutine test_direct_fs_equals_indirect(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: hi(2) = [2, 3], lo(2) = [1, 2]
  type(panel_mesh) :: full, mesh, fm
  type(surface_geom) :: geom
  type(cs_points) :: pts, circle
  type(direct_fs_data) :: dfs
  type(first_order_state), allocatable :: store(:, :)
  complex(8), allocatable :: store_wl(:, :, :), psi(:, :), psi_pts(:, :), gpsi(:, :, :), pwl(:, :)
  type(fs_point_field), allocatable, target :: fa(:, :, :), ia(:, :, :)
  type(pw_set), allocatable, target :: pwa(:, :)
  type(pw_set) :: pwp(6)
  real(8), allocatable :: delta(:)
  complex(8) :: t(6, 2), parts(6, 2, 3), terms(6, NPOT_TERMS), tin(6), tout(6), tfs(6, 2), tpot(6, 2)
  integer :: pairs(4, 2), npc, mmax, nev, iom, p
  real(8) :: centre(2), radius, omega, gap(2), opp(2), rel_fs(2), betas(3)
  integer(8) :: c0, c1, rate
  character(len=400) :: info

  call system_clock(c0, rate)
  call read_hull(round_mesh, [0d0, 0d0, -5d0], full, error)
  if (allocated(error)) return
  call set_omega(rd_omegas(3), rd_depth)
  call load_panels(full, 0, 0, mesh)
  call check(error, mesh%n == 480, "the fixture is not the 480-panel rounded cylinder")
  if (allocated(error)) return
  call fs_common(mesh, rd_depth, rd_omegas(3), geom, pts, fm, centre, radius, npc, delta, circle, mmax, nev, error)
  if (allocated(error)) return
  print '(4x,a,i0,a,f6.2,a,i0,a,i0,a,i0)', "inner points ", pts%n, ", R ", radius, ", polygon ", npc, ", M ", mmax, &
    ", N_ev ", nev
  betas = rd_beta
  call fs_states(mesh, geom, pts, delta, circle, centre, radius, rd_depth, mmax, nev, rd_omegas, betas, rd_motions(), 8, &
                 store, store_wl, fa, ia, pwa)
  allocate (gpsi(3, mesh%n, 6), pwl(geom%nwl, 6))
  do iom = 1, 2
    omega = rd_omegas(hi(iom)) - rd_omegas(lo(iom))
    pairs(:, 1) = [hi(iom), 1, lo(iom), 1]
    pairs(:, 2) = [lo(iom), 1, hi(iom), 1]
    call radiation_at(omega, mesh%nbase)
    NTHREAD = 8
    call psi_fs(mesh, pts, circle, centre, radius, rd_depth, mmax, nev, omega, psi, psi_pts, pwp)
    call psi_surface_gradient(mesh, geom, psi, gpsi)
    call psi_waterline(mesh, geom, psi, gpsi, pwl)
    ! The indirect method: the body terms and FS_inner + FS_outer of each pair
    do p = 1, 2
      associate (s1 => store(pairs(1, p), 1), s2 => store(pairs(3, p), 1), i1 => pairs(1, p), i2 => pairs(3, p))
        call potential_pair(mesh, geom, s1, s2, store_wl(:, i1, 1), store_wl(:, i2, 1), psi, gpsi, pwl, omega, s1%rho, &
                            wave_gravity(s1%omega, s1%k, s1%h), terms)
        call fs_pair_terms(s1%omega, fa(:, i1, 1), ia(:, i1, 1), pwa(i1, 1), s1%beta, s1%k, &
                           s2%omega, fa(:, i2, 1), ia(:, i2, 1), pwa(i2, 1), s2%beta, s2%k, &
                           psi_pts, pwp, pts%w, npc, s1%g, s1%rho, tin, tout)
      end associate
      tfs(:, p) = tin + tout
      tpot(:, p) = sum(terms, dim=2) + tfs(:, p)
    end do
    ! The direct method with the free-surface forcing
    call direct_setup(mesh, DEFAULT_EDGE_DEG)
    call direct_fs_setup(mesh, pts, fm, centre, radius, rd_depth, omega, .false., dfs)
    dfs%f => fa
    dfs%inc => ia
    dfs%pw => pwa
    dfs%np = npc
    call direct_at_omega(mesh, geom, omega, store, store_wl, pairs, t, parts=parts, fs=dfs)
    NTHREAD = 1
    gap(iom) = maxval(abs(parts(:, :, 3) - tfs)) / maxval(abs(tpot))
    opp(iom) = maxval(abs(-parts(:, :, 3) - tfs)) / maxval(abs(tpot))
    rel_fs(iom) = maxval(abs(tfs)) / maxval(abs(tpot))
    write (info, '(a,f5.2,a,es10.3,a,es10.3,a,es10.3,a,i0,a,i0,a,i0,a,es9.2)') "Omega ", omega, &
      ": |slice 3 - (FS_inner + FS_outer)| / max |T_pot| ", gap(iom), " (opposite s_D ", opp(iom), &
      "; max |T_FS| / max |T_pot| ", rel_fs(iom), "); mg ", dfs%mg, ", nev ", dfs%nev, ", mcut ", dfs%mcut, &
      ", cut error ", dfs%cut_error
    print '(4x,a)', trim(info)
    do p = 1, 2
      print '(6x,a,i0,a,6es11.3)', "pair ", p, " slice 3 ", abs(parts(:, p, 3))
      print '(6x,a,i0,a,6es11.3)', "pair ", p, " T_FS    ", abs(tfs(:, p))
    end do
  end do
  call system_clock(c1)
  write (info, '(a,2es10.3,a,es9.2,a,2es10.3,a,f7.1,a)') "gap at Omega 0.25, 0.15: ", gap, " (limit ", fs_indirect_limit, &
    "); opposite s_D: ", opp, "; ", real(c1 - c0, 8) / real(rate, 8), " s"
  print '(4x,a)', trim(info)
  call check(error, DIRECT_FS_SIGN == 1 .and. all(gap <= fs_indirect_limit) .and. all(opp > 10d0 * fs_indirect_limit), &
             trim(info))
end subroutine test_direct_fs_equals_indirect

end module test_qtf_direct_fs
