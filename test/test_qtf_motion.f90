!> Unit tests of the second-order motion of free bodies (Phase 6c, QTFMotion): the added mass and damping at Omega from
!> the stored radiation potentials psi (radiation_coefficients) against RFORCE and RFORCE_MULTI of the first-order solve,
!> on a half mesh, for reciprocity, at two sum frequencies against the semi-analytical sharp cylinder, and the matrix of
!> the equation of motion (motion_matrix); the motion solve per Omega (motion_solve_omega) against SolveMotion at a
!> main-grid frequency, a closed form, the static limit, a singular system and column by column.
!> The radiation solves set module globals (the mesh arrays, the stored Green function, the wave), so the suite runs
!> serially (see unit_tester.f90).
module test_qtf_motion
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use HAMS_mod, only : NTHREAD, ISOLV
  use WaveDyn_mod, only : WK, W1, TP, BETA, AMP, OUFR, BLNR, BQDR
  use Body_mod, only : MATX, CRS, KSTF, REFL
  use Const_mod, only : RHO
  use PanelMesh_mod, only : ISX, ISY, NSYS, IPIV
  use Potentials_mod, only : MXPOT, MXPOT_MULTI_COMB
  use LinearMatrix_mod, only : AMAT
  use AssbMatx, only : ASSB_DBC, DIFFRACTION_SOLVER
  use PotentWavForce, only : RFORCE, EFORCE
  use PotentWavForceMulti, only : RFORCE_MULTI
  use SecondOrderTypes, only : panel_mesh, so_options, surface_geom, first_order_state, QTF_TERMS_BODY
  use FirstOrderState, only : radiation_psi
  use QTFMotion, only : radiation_coefficients, motion_matrix, motion_solve_omega
  use QTFDriver, only : qtf_init, qtf_keep, qtf_keep_psi, qtf_radiation_omegas, qtf_finish, potential_part, &
                        qtf_motion_part, qtf_motion_load
  use NearFieldDrift, only : KERNEL_DIFF
  use SurfaceGeometry, only : build_surface_geom, DEFAULT_SHARP_DEG
  use IO, only : CreateDirectory
  use test_qtf_potential, only : make_state
  use second_order_fixtures, only : read_pnl_hull
  use test_waterline_field, only : set_wave, load_panels, load_multi
  use test_qtf_direct, only : radiation_at
  implicit none
  private

  public :: collect_qtf_motion

  ! The single-body fixture of qtf_round_body (480 panels, rounded cylinder of radius 10 and draft 10, depth 100) and two
  ! of its main-grid frequencies
  character(len=*), parameter :: round_mesh = "test/data/qtf_round_body/input/HullMesh.pnl"
  real(8), parameter :: round_depth = 100d0, round_omegas(2) = [0.7d0, 0.95d0]
  ! The same panel sums as RFORCE in the same order: rounding only, relative to max |A| and max |B|
  real(8), parameter :: rforce_limit = 1d-12
  ! Half mesh (ISY = 1) against the whole: the same sums over the expanded images, rounding only
  real(8), parameter :: half_limit = 1d-10
  ! Half-mesh solve against the whole-mesh solve: the single-precision Green function at the mirrored points (6e-8 per
  ! value, see stored_limit of test_qtf_direct); 3.8e-7 measured on qtf_round_body at 0.7 rad/s
  real(8), parameter :: half_solve_limit = 1d-6
  ! Reciprocity, A_lk = A_kl and B_lk = B_kl (Newman 1977, section 6.15), relative to max |A| and max |B|. It holds for the
  ! exact potentials; the constant-panel collocation of HAMS (whose AMAS and BDMP these equal to rounding) meets it only to
  ! its mesh error. Measured: 6.0e-4 for qtf_round_body (480 panels, B at roll-sway (4, 2)) and 4.2e-2 for the
  ! two cylinders of load_multi (48 panels each, B at (9, 5)); the limits are 1.3 x these
  real(8), parameter :: symmetric_limit_one = 7.8d-4, symmetric_limit_two = 5.5d-2
  ! motion_matrix with and without B_ext: the difference is -i Omega B_ext, rounding only
  real(8), parameter :: blnr_limit = 1d-14
  ! The smallest sharp-cylinder fixture (qtf_direct_sharp_m1: 448 panels, radius 10, draft 10, depth 100, the reference
  ! point at the origin)
  character(len=*), parameter :: sharp_mesh = "test/data/qtf_direct_sharp_m1/input/HullMesh.pnl"
  real(8), parameter :: sharp_depth = 100d0
  ! Two sum frequencies off the main grid, where the 7c-2 free-surface error rises
  real(8), parameter :: sum_omegas(2) = [1.8d0, 2.8d0]
  ! The reference: added_mass_damping(cyl, omega, rads) of pre_postprocessing_python/SecondOrder/
  ! truncated_cylinder_reference.py at DEFAULT_RESOLUTION (n_modes = 64000, p_max = (20, 4, 4), orders = 15), with
  ! cyl = Cylinder(10.0, 10.0, 100.0) and rads = {j: solve_radiation(cyl, omega, j, 64000, (20, 4, 4)) for j in (1, 3, 5)},
  ! run as python -I ref_am.py <path of truncated_cylinder_reference.py> default 0.7 1.8 2.8 (a 20-line driver that imports
  ! the script by its path). Its A is A/rho and its B is B/(rho omega), moments about the origin: the AmssDamp.1 form.
  ! The entries (l, k) = (1, 1), (1, 5), (3, 3), (5, 1), (5, 5); the reference has the others zero (by symmetry).
  integer, parameter :: ref_l(5) = [1, 1, 3, 5, 5], ref_k(5) = [1, 5, 3, 1, 5]
  real(8), parameter :: ref_a(5, 2) = reshape([ &
    4.397910074760794d+02, -2.435306257555518d+03, 1.788230232032700d+03, -2.435306570088258d+03, 3.493824493134971d+04, &
    6.904053255396620d+02, -3.477723757647814d+03, 1.822455102718240d+03, -3.477724075227684d+03, 3.862359208671882d+04], &
    [5, 2])
  real(8), parameter :: ref_b(5, 2) = reshape([ &
    5.313649599917072d+02, -1.311205787299898d+03, 8.191339661078715d-01, -1.311205789505384d+03, 3.235555124988715d+03, &
    9.838944040298490d+01, -1.224452426278691d+02, 2.901765536094465d-05, -1.224452426347644d+02, 1.523825868064106d+02], &
    [5, 2])
  ! The gap of each entry is relative to the larger of |reference| and 1% of the largest reference entry of A or B at
  ! that frequency (heave B at 2.8 rad/s is 3e-5, against 152 for pitch). Limit: 1.3 x the gap measured at the first
  ! passing run: A 4.03e-2 (pitch-surge at 2.8 rad/s), B 1.30e-1 (pitch at 2.8 rad/s). The same check at the main-grid
  ! frequency 0.7 rad/s of the fixture gave A 3.1e-2 and B 5.0e-2 (reference (1, 1), (1, 5), (3, 3), (5, 5): 2492.18,
  ! -8852.26, 1748.00, 54844.27 for A and 544.31, -1706.77, 419.78, 5351.79 for B): the normalisation, rho and the mode
  ! order agree, and the gap is the mesh error of the 448 panels, not larger at the sum frequencies than on the grid
  ! except for B at 2.8 rad/s.
  real(8), parameter :: sum_limit_a = 5.24d-2, sum_limit_b = 1.68d-1
  ! The first-order limit (R4): motion_solve_omega at omega = 0.7 rad/s against conj(DSPL) of SolveMotion, relative to
  ! max |DSPL|. The two sides differ only by a, b of radiation_coefficients against AMAS, BDMP of RFORCE (rounding, see
  ! rforce_limit) and the two LU solves
  real(8), parameter :: first_order_limit = 1d-10
  ! The closed form of one degree of freedom and the static diagonal: rounding of one complex division per entry
  real(8), parameter :: closed_limit = 1d-14
  ! The free body of the first-order limit: M and C of qtf_round_body/input/Hydrostatic.in (the hull of round_mesh), an
  ! external spring in surge, sway, heave and yaw so that C + K_ext is regular, and an external linear damper in surge
  ! and heave so that the sign of -i omega B_ext is pinned as well. BQDR = 0 (SolveMotion's linear branch).
  real(8), parameter :: round_mass(6, 6) = reshape([ &
    3.00895d+06, 0d0, 0d0, 0d0, -1.45599d+07, 0d0, &
    0d0, 3.00895d+06, 0d0, 1.45599d+07, 0d0, 0d0, &
    0d0, 0d0, 3.00895d+06, 0d0, 0d0, 0d0, &
    0d0, 1.45599d+07, 0d0, 1.66881d+08, 0d0, 0d0, &
    -1.45599d+07, 0d0, 0d0, 0d0, 1.66881d+08, 0d0, &
    0d0, 0d0, 0d0, 0d0, 0d0, 1.45239d+08], [6, 6])
  real(8), parameter :: round_c_diag(6) = [0d0, 0d0, 3.06109d+06, 7.60372d+07, 7.60372d+07, 0d0]
  real(8), parameter :: round_k_diag(6) = [1d6, 1d6, 5d5, 0d0, 0d0, 1d8]
  real(8), parameter :: round_blnr_diag(6) = [2d5, 0d0, 1d5, 0d0, 0d0, 0d0]

contains

  !> Collect all exported unit tests
  subroutine collect_qtf_motion(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("radiation_coefficients_match_first_order", test_match_first_order), &
      new_unittest("radiation_coefficients_half", test_half), &
      new_unittest("radiation_coefficients_two_bodies", test_two_bodies), &
      new_unittest("radiation_coefficients_symmetric", test_symmetric), &
      new_unittest("motion_damping_includes_blnr", test_blnr), &
      new_unittest("radiation_coefficients_sum_frequency", test_sum_frequency), &
      new_unittest("motion_first_order_limit", test_first_order_limit), &
      new_unittest("motion_one_dof", test_one_dof), &
      new_unittest("motion_static_diagonal", test_static_diagonal), &
      new_unittest("motion_singular_flagged", test_singular_flagged), &
      new_unittest("motion_columns_independent", test_columns_independent), &
      new_unittest("motion_two_bodies_coupled", test_two_bodies_coupled) &
    ]
  end subroutine collect_qtf_motion

  !> The hull of a fixture, whole, into the single-body arrays (load_panels) with the symmetry isy_in
  subroutine load_hull(path, depth, isy_in, mesh, error)
    character(len=*), intent(in) :: path
    real(8), intent(in) :: depth
    integer, intent(in) :: isy_in
    type(panel_mesh), intent(out) :: mesh
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full
    logical :: ok

    call read_pnl_hull(path, [0d0, 0d0, -5d0], full, ok)
    call check(error, ok, "cannot open "//path)
    if (allocated(error)) return
    call set_wave(depth)
    NTHREAD = 4
    ISOLV = 1
    call load_panels(full, 0, isy_in, mesh)
  end subroutine load_hull

  !> radiation_coefficients of the radiation solve at omega on mesh (single-body or, with multi, the combined arrays), and
  !> the psi of radiation_psi it takes
  subroutine coefficients_at(mesh, omega, multi, a, b, psi)
    type(panel_mesh), intent(in) :: mesh
    real(8), intent(in) :: omega
    logical, intent(in) :: multi
    real(8), allocatable, intent(out) :: a(:, :), b(:, :)
    complex(8), allocatable, intent(out) :: psi(:, :)

    call radiation_at(omega, multi, mesh%nbase)
    allocate(psi(mesh%n, 6*mesh%nbody), a(6*mesh%nbody, 6*mesh%nbody), b(6*mesh%nbody, 6*mesh%nbody))
    if (multi) then
      call radiation_psi(mesh, omega, MXPOT_MULTI_COMB, psi)
    else
      call radiation_psi(mesh, omega, MXPOT, psi)
    end if
    call radiation_coefficients(mesh, psi, omega, RHO, a, b)
  end subroutine coefficients_at

  !> RFORCE (nbody = 0) or the RFORCE_MULTI blocks as HAMS_Prog assembles AMAS_MULTI_COMB (nbody bodies of ne_body base
  !> panels each) on the solver state of the last radiation solve. Their output files (units 31-36, 41-46 and 61, opened
  !> by InputFiles in a run) go to scratch files here.
  subroutine rforce_reference(nbody, ne_body, amas, bdmp)
    integer, intent(in) :: nbody, ne_body
    real(8), allocatable, intent(out) :: amas(:, :), bdmp(:, :)
    integer, parameter :: units(13) = [31, 32, 33, 34, 35, 36, 41, 42, 43, 44, 45, 46, 61]
    logical :: was_open(13)
    integer :: u, m, n

    do u = 1, size(units)
      inquire(unit=units(u), opened=was_open(u))
      if (.not. was_open(u)) open(unit=units(u), status="scratch")
    end do
    if (nbody == 0) then
      allocate(amas(6, 6), bdmp(6, 6))
      call RFORCE(WK, W1, TP, amas, bdmp)
    else
      allocate(amas(6*nbody, 6*nbody), bdmp(6*nbody, 6*nbody))
      do m = 1, nbody
        do n = 1, nbody
          call RFORCE_MULTI(WK, W1, TP, amas(6*(m - 1) + 1:6*m, 6*(n - 1) + 1:6*n), bdmp(6*(m - 1) + 1:6*m, &
                            6*(n - 1) + 1:6*n), (m - 1)*ne_body + 1, m*ne_body, n, m, n)
        end do
      end do
    end if
    do u = 1, size(units)
      if (.not. was_open(u)) close(units(u))
    end do
  end subroutine rforce_reference

  !> The larger of the differences of A and of B, each relative to the largest entry of the reference
  real(8) function rel_gap(a, b, a_ref, b_ref)
    real(8), intent(in) :: a(:, :), b(:, :), a_ref(:, :), b_ref(:, :)
    rel_gap = max(maxval(abs(a - a_ref))/maxval(abs(a_ref)), maxval(abs(b - b_ref))/maxval(abs(b_ref)))
  end function rel_gap

  subroutine restore_symmetry()
    ISOLV = 1
    ISX = 0; ISY = 0; NSYS = 1
  end subroutine restore_symmetry

  !> qtf_round_body at two main-grid frequencies: a, b of the stored psi equal AMAS, BDMP of RFORCE on the same solve
  subroutine test_match_first_order(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    real(8), allocatable :: a(:, :), b(:, :), amas(:, :), bdmp(:, :)
    complex(8), allocatable :: psi(:, :)
    real(8) :: worst, w
    integer :: i
    character(len=200) :: line

    call load_hull(round_mesh, round_depth, 0, mesh, error)
    if (allocated(error)) return
    worst = 0d0
    do i = 1, size(round_omegas)
      call coefficients_at(mesh, round_omegas(i), .false., a, b, psi)
      call rforce_reference(0, mesh%nbase, amas, bdmp)
      w = rel_gap(a, b, amas, bdmp)
      write(line, '(a,f5.2,a,i0,a,es10.2,a,es10.2,a,es10.2)') "omega ", round_omegas(i), ", ", mesh%n, &
        " panels: max |A| ", maxval(abs(amas)), ", max |B| ", maxval(abs(bdmp)), ", against RFORCE relative ", w
      print '(4x,a)', trim(line)
      worst = max(worst, w)
    end do
    call restore_symmetry()
    write(line, '(a,es10.2)') "radiation_coefficients against RFORCE, worst relative difference ", worst
    call check(error, worst <= rforce_limit, line)
  end subroutine test_match_first_order

  !> qtf_round_body as a half mesh (ISY = 1, the 240 panels of qtf_half_body) against its whole mesh: psi of the half
  !> solve, expanded over the images by radiation_psi, put on the matching panels of the whole mesh, gives the same A and
  !> B there, so the image signs that mesh%nj of the image panels carries are those of the whole hull; and against RFORCE
  !> on the half mesh, which applies them as the RY factors. The whole-mesh solve is printed besides: it differs by the
  !> single-precision storage of the Green function (CGRN, GRN_KIND = 4) at the mirrored points, not by rounding.
  subroutine test_half(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: whole, mesh
    real(8), allocatable :: a(:, :), b(:, :), ah(:, :), bh(:, :), aw(:, :), bw(:, :), amas(:, :), bdmp(:, :)
    complex(8), allocatable :: psi(:, :), psi_half(:, :), psi_whole(:, :)
    integer, allocatable :: match(:)
    real(8) :: w_whole, w_solve, w_rforce
    real(8), parameter :: om = 0.7d0
    integer :: i, j
    character(len=250) :: line

    call load_hull(round_mesh, round_depth, 0, whole, error)
    if (allocated(error)) return
    call coefficients_at(whole, om, .false., a, b, psi)
    call load_hull(round_mesh, round_depth, 1, mesh, error)
    if (allocated(error)) return
    call check(error, mesh%nbase == 240 .and. mesh%nsys == 2 .and. mesh%n == whole%n, &
               "the half mesh does not have 240 base panels and 2 images")
    if (allocated(error)) then
      call restore_symmetry()
      return
    end if
    call coefficients_at(mesh, om, .false., ah, bh, psi_half)
    call rforce_reference(0, mesh%nbase, amas, bdmp)
    call restore_symmetry()
    ! The panel of the whole mesh at the centroid of each panel of the expanded half mesh
    allocate(match(mesh%n), psi_whole(whole%n, 6), aw(6, 6), bw(6, 6))
    match = 0
    do i = 1, mesh%n
      do j = 1, whole%n
        if (norm2(mesh%xc(:, i) - whole%xc(:, j)) < 1d-6) match(i) = j
      end do
    end do
    call check(error, all(match > 0) .and. count([(any(match == j), j = 1, whole%n)]) == whole%n, &
               "the expanded half mesh does not cover the whole mesh panel by panel")
    if (allocated(error)) return
    psi_whole(match, :) = psi_half
    call radiation_coefficients(whole, psi_whole, om, RHO, aw, bw)
    w_whole = rel_gap(ah, bh, aw, bw)
    w_rforce = rel_gap(ah, bh, amas, bdmp)
    w_solve = rel_gap(ah, bh, a, b)
    write(line, '(a,es10.2,a,es10.2,a,es10.2)') "half mesh against the whole with the same psi, relative ", w_whole, &
      "; against RFORCE on the half ", w_rforce, "; against the whole-mesh solve ", w_solve
    print '(4x,a)', trim(line)
    call check(error, w_whole <= half_limit .and. w_rforce <= rforce_limit .and. w_solve <= half_solve_limit, line)
  end subroutine test_half

  !> Two cylinders in the combined arrays (load_multi, depth 50): a, b equal the AMAS_MULTI_COMB, BDMP_MULTI_COMB blocks of
  !> RFORCE_MULTI, the off-diagonal body blocks included
  subroutine test_two_bodies(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    real(8), allocatable :: a(:, :), b(:, :), amas(:, :), bdmp(:, :)
    complex(8), allocatable :: psi(:, :)
    real(8) :: w, off
    real(8), parameter :: om = 0.8d0
    character(len=200) :: line

    call set_wave(50d0)
    NTHREAD = 4
    ISOLV = 1
    call load_multi(0, mesh)
    call coefficients_at(mesh, om, .true., a, b, psi)
    call rforce_reference(2, mesh%nbase/2, amas, bdmp)
    call restore_symmetry()
    w = rel_gap(a, b, amas, bdmp)
    off = maxval(abs(amas(1:6, 7:12)))/maxval(abs(amas))
    write(line, '(a,i0,a,es10.2,a,es10.2)') "two bodies, ", mesh%n, " panels: against RFORCE_MULTI relative ", w, &
      "; off-diagonal block of A relative ", off
    print '(4x,a)', trim(line)
    call check(error, w <= rforce_limit .and. off > 1d-6, line)
  end subroutine test_two_bodies

  !> Reciprocity: A and B symmetric, qtf_round_body at a main-grid frequency and the two cylinders (both body blocks)
  subroutine test_symmetric(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    real(8), allocatable :: a(:, :), b(:, :)
    complex(8), allocatable :: psi(:, :)
    real(8) :: w1s, w2s
    character(len=200) :: line

    call load_hull(round_mesh, round_depth, 0, mesh, error)
    if (allocated(error)) return
    call coefficients_at(mesh, round_omegas(1), .false., a, b, psi)
    w1s = rel_gap(a, b, transpose(a), transpose(b))
    call set_wave(50d0)
    NTHREAD = 4
    call load_multi(0, mesh)
    call coefficients_at(mesh, 0.8d0, .true., a, b, psi)
    call restore_symmetry()
    w2s = rel_gap(a, b, transpose(a), transpose(b))
    write(line, '(a,es10.2,a,es10.2)') "|A - A^T|, |B - B^T| relative: one body ", w1s, ", two bodies ", w2s
    print '(4x,a)', trim(line)
    call check(error, w1s <= symmetric_limit_one .and. w2s <= symmetric_limit_two, line)
  end subroutine test_symmetric

  !> motion_matrix with B_ext /= 0 differs from B_ext = 0 by -i Omega B_ext alone, and its real part is
  !> -Omega^2 (M + A) + C + K_ext
  subroutine test_blnr(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: n = 12
    real(8) :: mass(n, n), a(n, n), b(n, n), b_ext(n, n), c(n, n), k_ext(n, n), zero(n, n)
    complex(8) :: lhs0(n, n), lhs1(n, n)
    real(8), parameter :: om = 1.3d0
    real(8) :: w_damp, w_real
    integer :: i, j
    character(len=200) :: line

    do j = 1, n
      do i = 1, n
        mass(i, j) = sin(0.3d0*i + 0.7d0*j)
        a(i, j) = cos(1.1d0*i - 0.4d0*j)
        b(i, j) = sin(0.9d0*i*j)
        b_ext(i, j) = 1d3*cos(0.2d0*i + 1.3d0*j)
        c(i, j) = 5d0*sin(1.7d0*i - 0.6d0*j)
        k_ext(i, j) = 2d0*cos(0.5d0*i*j)
      end do
    end do
    zero = 0d0
    call motion_matrix(om, mass, a, b, zero, c, k_ext, lhs0)
    call motion_matrix(om, mass, a, b, b_ext, c, k_ext, lhs1)
    w_damp = maxval(abs(lhs1 - lhs0 - cmplx(0d0, -om*b_ext, kind=8)))/maxval(abs(om*b_ext))
    w_real = maxval(abs(real(lhs0, 8) - (-om**2*(mass + a) + c + k_ext)))/maxval(abs(lhs0)) + &
             maxval(abs(aimag(lhs0) + om*b))/maxval(abs(om*b))
    write(line, '(a,es10.2,a,es10.2)') "lhs(B_ext) - lhs(0) + i Omega B_ext relative ", w_damp, &
      "; lhs(0) against -Omega^2 (M + A) - i Omega B + C + K_ext ", w_real
    print '(4x,a)', trim(line)
    call check(error, w_damp <= blnr_limit .and. w_real <= blnr_limit, line)
  end subroutine test_blnr

  !> The sharp cylinder (qtf_direct_sharp_m1) at the sum frequencies 1.8 and 2.8 rad/s: A/rho and B/(rho Omega) against
  !> the semi-analytical reference. Guards the motion at the sum frequencies against an error of psi at high Omega.
  subroutine test_sum_frequency(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    real(8), allocatable :: a(:, :), b(:, :)
    complex(8), allocatable :: psi(:, :)
    real(8) :: gap_a, gap_b, worst_a, worst_b, ha, hb, floor_a, floor_b
    integer :: i, e
    character(len=250) :: line

    call load_hull(sharp_mesh, sharp_depth, 0, mesh, error)
    if (allocated(error)) return
    worst_a = 0d0
    worst_b = 0d0
    do i = 1, size(sum_omegas)
      call coefficients_at(mesh, sum_omegas(i), .false., a, b, psi)
      floor_a = 1d-2*maxval(abs(ref_a(:, i)))
      floor_b = 1d-2*maxval(abs(ref_b(:, i)))
      gap_a = 0d0
      gap_b = 0d0
      do e = 1, size(ref_l)
        ha = a(ref_l(e), ref_k(e))/RHO
        hb = b(ref_l(e), ref_k(e))/(RHO*sum_omegas(i))
        write(line, '(a,f4.1,a,2i2,a,2es16.8,a,2es16.8)') "SUMFREQ ", sum_omegas(i), " (l, k)", ref_l(e), ref_k(e), &
          ": A/rho HAMS, ref ", ha, ref_a(e, i), "; B/(rho Omega) HAMS, ref ", hb, ref_b(e, i)
        print '(4x,a)', trim(line)
        gap_a = max(gap_a, abs(ha - ref_a(e, i))/max(abs(ref_a(e, i)), floor_a))
        gap_b = max(gap_b, abs(hb - ref_b(e, i))/max(abs(ref_b(e, i)), floor_b))
      end do
      write(line, '(a,f4.1,a,es10.3,a,es10.3)') "Omega ", sum_omegas(i), ": gap of A ", gap_a, ", of B ", gap_b
      print '(4x,a)', trim(line)
      worst_a = max(worst_a, gap_a)
      worst_b = max(worst_b, gap_b)
    end do
    call restore_symmetry()
    write(line, '(a,es10.3,a,es10.3,a,es10.3,a,es10.3)') "sharp cylinder against the reference: gap of A ", worst_a, &
      " (limit ", sum_limit_a, "), of B ", worst_b, " (limit ", sum_limit_b
    call check(error, worst_a <= sum_limit_a .and. worst_b <= sum_limit_b, trim(line)//")")
  end subroutine test_sum_frequency

  !> The first-order limit (R4). At the main-grid frequency 0.7 rad/s of qtf_round_body, HAMS_Prog's sequence in-process
  !> (the radiation solve, RFORCE, ASSB_DBC and DIFFRACTION_SOLVER on the factored AMAT, EFORCE, SolveMotion) gives DSPL
  !> in e^{+i omega t}. SolveMotion solves LEFT DSPL = RIGHT with
  !>    LEFT = -omega^2 (MATX + AMAS) + i omega (BDMP + BLNR) + CRS + KSTF,   RIGHT = -Im EXFC - i Re EXFC = -i conj(EXFC).
  !> The conjugate of that system is conj(LEFT) conj(DSPL) = i EXFC, and conj(LEFT) is the e^{-i omega t} matrix (2) of
  !> QTFMotion. So the internal exciting force in e^{-i Omega t} is F = i EXFC, and motion_solve_omega of motion_matrix
  !> (on a, b of radiation_coefficients of the same solve) with F gives xi = conj(DSPL). The e^{-i Omega t} side is built
  !> by motion_matrix alone, so this pins its -i Omega (B + B_ext) against SolveMotion's +i omega.
  subroutine test_first_order_limit(error)
    type(error_type), allocatable, intent(out) :: error
    interface
      subroutine SolveMotion(w, t, oufr_in, beta_in, amp_in, amas_in, bdmp_in, blnr_in, bqdr_in, exfc_in, dspl_out)
        real(8), intent(in) :: w, t, oufr_in, beta_in, amp_in
        real(8), intent(in) :: amas_in(6, 6), bdmp_in(6, 6), blnr_in(6, 6), bqdr_in(6, 6)
        complex(8), intent(in) :: exfc_in(6)
        complex(8), intent(out) :: dspl_out(6)
      end subroutine SolveMotion
    end interface
    ! Units written by EFORCE (21-26, 62) and SolveMotion (63), opened by InputFiles in a run
    integer, parameter :: units(8) = [21, 22, 23, 24, 25, 26, 62, 63]
    real(8), parameter :: om = 0.7d0
    type(panel_mesh) :: mesh
    real(8), allocatable :: a(:, :), b(:, :), amas(:, :), bdmp(:, :)
    complex(8), allocatable :: psi(:, :), bdmat(:, :)
    real(8) :: saved(6, 6, 5), refl_saved, rcond, worst, scale
    complex(8) :: exfc(6), dspl(6), lhs(6, 6), force(6, 1), xi(6, 1)
    logical :: was_open(size(units)), ok
    integer :: md, u
    character(len=200) :: line

    call load_hull(round_mesh, round_depth, 0, mesh, error)
    if (allocated(error)) return
    call coefficients_at(mesh, om, .false., a, b, psi)
    call rforce_reference(0, mesh%nbase, amas, bdmp)
    allocate(bdmat(mesh%nbase, NSYS))
    call ASSB_DBC(bdmat, mesh%nbase, NSYS)
    call DIFFRACTION_SOLVER(AMAT, bdmat, IPIV, MXPOT, mesh%nbase, NSYS)
    ! The free body in HAMS's globals, as ReadHydroStatic sets them; the old values are put back at the end
    saved(:, :, 1) = MATX
    saved(:, :, 2) = CRS
    saved(:, :, 3) = KSTF
    saved(:, :, 4) = BLNR
    saved(:, :, 5) = BQDR
    refl_saved = REFL
    MATX = round_mass
    CRS = 0d0
    KSTF = 0d0
    BLNR = 0d0
    do md = 1, 6
      CRS(md, md) = round_c_diag(md)
      KSTF(md, md) = round_k_diag(md)
      BLNR(md, md) = round_blnr_diag(md)
    end do
    BQDR = 0d0
    REFL = 1d0
    do u = 1, size(units)
      inquire(unit=units(u), opened=was_open(u))
      if (.not. was_open(u)) open(unit=units(u), status="scratch")
    end do
    call EFORCE(WK, W1, TP, BETA, AMP, exfc)
    call SolveMotion(W1, TP, OUFR, BETA, AMP, amas, bdmp, BLNR, BQDR, exfc, dspl)
    do u = 1, size(units)
      if (.not. was_open(u)) close(units(u))
    end do
    ! The e^{-i Omega t} side: Task 2's matrix on the stored psi, and F = i EXFC
    call motion_matrix(om, MATX, a, b, BLNR, CRS, KSTF, lhs)
    force(:, 1) = (0d0, 1d0)*exfc
    call motion_solve_omega(lhs, force, xi, ok, rcond)
    MATX = saved(:, :, 1)
    CRS = saved(:, :, 2)
    KSTF = saved(:, :, 3)
    BLNR = saved(:, :, 4)
    BQDR = saved(:, :, 5)
    REFL = refl_saved
    call restore_symmetry()
    scale = maxval(abs(dspl))
    worst = maxval(abs(xi(:, 1) - conjg(dspl)))/scale
    do md = 1, 6
      write(line, '(a,i0,a,es24.16,a,es24.16,a,es10.2)') "FIRSTORDER mode ", md, ", |conj(DSPL)| ", abs(dspl(md)), &
        ", |xi| ", abs(xi(md, 1)), ", relative ", abs(xi(md, 1) - conjg(dspl(md)))/max(abs(dspl(md)), tiny(1d0))
      print '(4x,a)', trim(line)
    end do
    write(line, '(a,l1,a,es10.2,a,es10.2)') "omega 0.7: ok ", ok, ", rcond ", rcond, &
      ", |xi - conj(DSPL)| relative to max |DSPL| ", worst
    print '(4x,a)', trim(line)
    call check(error, ok .and. worst <= first_order_limit, line)
  end subroutine test_first_order_limit

  !> One degree of freedom, m = 2, a = 0.5, b = 0.3, c = 4 at Omega = 1.1 with F = 1 + 2i, against the closed form
  !> xi = F/(c - Omega^2 (m + a) - i Omega b)
  subroutine test_one_dof(error)
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: m = 2d0, am = 0.5d0, bm = 0.3d0, cm = 4d0, om = 1.1d0
    complex(8), parameter :: f = (1d0, 2d0)
    real(8) :: zero(1, 1), rcond, w
    complex(8) :: lhs(1, 1), force(1, 1), xi(1, 1), expected
    logical :: ok
    character(len=200) :: line

    zero = 0d0
    call motion_matrix(om, reshape([m], [1, 1]), reshape([am], [1, 1]), reshape([bm], [1, 1]), zero, &
                       reshape([cm], [1, 1]), zero, lhs)
    force = f
    call motion_solve_omega(lhs, force, xi, ok, rcond)
    expected = f/cmplx(cm - om**2*(m + am), -om*bm, kind=8)
    w = abs(xi(1, 1) - expected)/abs(expected)
    write(line, '(a,l1,a,es10.2,a,es10.2)') "one dof: ok ", ok, ", rcond ", rcond, ", relative to the closed form ", w
    print '(4x,a)', trim(line)
    call check(error, ok .and. w <= closed_limit, line)
  end subroutine test_one_dof

  !> Omega = 0 with C + K_ext regular and diagonal: xi = (C + K_ext)^-1 F, whatever M, A and B are
  subroutine test_static_diagonal(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: n = 6, npair = 2
    real(8) :: mass(n, n), a(n, n), b(n, n), b_ext(n, n), c(n, n), k_ext(n, n), rcond, w
    complex(8) :: lhs(n, n), force(n, npair), xi(n, npair), expected(n, npair)
    logical :: ok
    integer :: i, j
    character(len=200) :: line

    c = 0d0
    k_ext = 0d0
    do j = 1, n
      do i = 1, n
        mass(i, j) = 1d3*sin(0.3d0*i + 0.7d0*j)
        a(i, j) = 1d2*cos(1.1d0*i - 0.4d0*j)
        b(i, j) = 1d2*sin(0.9d0*i*j)
        b_ext(i, j) = 1d1*cos(0.2d0*i + 1.3d0*j)
      end do
      c(j, j) = 1d5*j
      k_ext(j, j) = 3d4 + 1d3*j*j
      do i = 1, npair
        force(j, i) = cmplx(cos(1.7d0*i*j), sin(0.4d0*i + j), kind=8)
      end do
    end do
    call motion_matrix(0d0, mass, a, b, b_ext, c, k_ext, lhs)
    call motion_solve_omega(lhs, force, xi, ok, rcond)
    do j = 1, n
      expected(j, :) = force(j, :)/(c(j, j) + k_ext(j, j))
    end do
    w = maxval(abs(xi - expected)/abs(expected))
    write(line, '(a,l1,a,es10.2,a,es10.2)') "static diagonal: ok ", ok, ", rcond ", rcond, &
      ", relative to (C + K_ext)^-1 F ", w
    print '(4x,a)', trim(line)
    call check(error, ok .and. w <= closed_limit, line)
  end subroutine test_static_diagonal

  !> Omega = 0 with no stiffness in surge: flagged (ok = .false.) with xi = 0. Once with C(1, 1) = 0 (exactly singular,
  !> ZGETRF info > 0) and once with C(1, 1) = 1e-8 against 1e6 to 6e6 of the others (factored, rcond about 2e-15, below
  !> the bound 1e-12)
  subroutine test_singular_flagged(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: n = 6
    real(8) :: zero(n, n), c(n, n), rcond_exact, rcond_near
    complex(8) :: lhs(n, n), force(n, 2), xi_exact(n, 2), xi_near(n, 2)
    logical :: ok_exact, ok_near
    integer :: j
    character(len=200) :: line

    zero = 0d0
    c = 0d0
    do j = 1, n
      c(j, j) = 1d6*j
      force(j, :) = [(1d0, 1d0), (2d0, -1d0)]
    end do
    c(1, 1) = 0d0
    call motion_matrix(0d0, zero, zero, zero, zero, c, zero, lhs)
    xi_exact = (9d0, 9d0)
    call motion_solve_omega(lhs, force, xi_exact, ok_exact, rcond_exact)
    c(1, 1) = 1d-8
    call motion_matrix(0d0, zero, zero, zero, zero, c, zero, lhs)
    xi_near = (9d0, 9d0)
    call motion_solve_omega(lhs, force, xi_near, ok_near, rcond_near)
    write(line, '(a,l1,a,es10.2,a,l1,a,es10.2)') "C(1, 1) = 0: ok ", ok_exact, ", rcond ", rcond_exact, &
      "; C(1, 1) = 1e-8: ok ", ok_near, ", rcond ", rcond_near
    print '(4x,a)', trim(line)
    call check(error, .not. ok_exact .and. rcond_exact == 0d0 .and. all(xi_exact == (0d0, 0d0)) .and. &
               .not. ok_near .and. rcond_near < 1d-12 .and. all(xi_near == (0d0, 0d0)), line)
  end subroutine test_singular_flagged

  !> Each force column is solved on its own. Bit for bit: ok and rcond of a three-column call equal those of the single
  !> calls, lhs is left as it was, and changing column 2 changes no bit of columns 1 and 3 (no mixing between columns).
  !> Against the three single-column calls, xi agrees to rounding only (column_limit): MKL's ZGETRS takes other kernels
  !> for one right-hand side than for three (2.3e-15 measured), so the columns are equal to the bit only within one call.
  subroutine test_columns_independent(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: n = 12, npair = 3
    ! Rounding of the two triangular solves, relative to max |xi| of the column: 1e-14 (2.3e-15 measured, 1/rcond = 360)
    real(8), parameter :: column_limit = 1d-14
    real(8) :: mass(n, n), a(n, n), b(n, n), b_ext(n, n), c(n, n), k_ext(n, n), rcond, rcond1, rcond2, w
    complex(8) :: lhs(n, n), lhs_in(n, n), force(n, npair), xi(n, npair), xi1(n, 1), xi2(n, npair)
    real(8), parameter :: om = 0.9d0
    logical :: ok, ok1, ok2, same
    integer :: i, j
    character(len=200) :: line

    do j = 1, n
      do i = 1, n
        mass(i, j) = sin(0.3d0*i + 0.7d0*j)
        a(i, j) = cos(1.1d0*i - 0.4d0*j)
        b(i, j) = sin(0.9d0*i*j)
        b_ext(i, j) = cos(0.2d0*i + 1.3d0*j)
        c(i, j) = 5d0*sin(1.7d0*i - 0.6d0*j)
        k_ext(i, j) = 2d0*cos(0.5d0*i*j)
      end do
      mass(j, j) = mass(j, j) + 10d0
      do i = 1, npair
        force(j, i) = cmplx(cos(1.7d0*i*j), sin(0.4d0*i + j), kind=8)
      end do
    end do
    call motion_matrix(om, mass, a, b, b_ext, c, k_ext, lhs)
    lhs_in = lhs
    call motion_solve_omega(lhs, force, xi, ok, rcond)
    same = ok .and. all(lhs == lhs_in)
    w = 0d0
    do i = 1, npair
      call motion_solve_omega(lhs, force(:, i:i), xi1, ok1, rcond1)
      same = same .and. ok1 .and. rcond1 == rcond
      w = max(w, maxval(abs(xi1(:, 1) - xi(:, i)))/maxval(abs(xi(:, i))))
    end do
    force(:, 2) = 3d0*force(:, 2) + (0d0, 1d0)
    call motion_solve_omega(lhs, force, xi2, ok2, rcond2)
    same = same .and. ok2 .and. rcond2 == rcond .and. all(xi2(:, 1) == xi(:, 1)) .and. all(xi2(:, 3) == xi(:, 3)) .and. &
           any(xi2(:, 2) /= xi(:, 2))
    write(line, '(a,l1,a,es10.2,a,l1,a,es10.2)') "three columns: ok ", ok, ", rcond ", rcond, &
      ", bitwise checks ", same, "; against single calls, relative ", w
    print '(4x,a)', trim(line)
    call check(error, same .and. w <= column_limit, line)
  end subroutine test_columns_independent

  !> Two bodies in one QTF store (Task 4): the two cylinders of load_multi (radius 5, draft 5, 30 m apart, depth 50),
  !> qtf = 1, qtf_terms = body, grid 0.6, 0.8 rad/s, one heading, and qtf_init with the 12 x 12 block-diagonal body
  !> matrices that the multi-body adapter of SecondOrderDriver builds. xi2 of the pair (0.8, 0.6) solves the coupled
  !> system at Omega = 0.2 (backward error |lhs xi - F| / (|lhs| |xi|), infinity norms, 1e-13), lhs from the radiation
  !> solve of the test. Then only the force on body 2 changes (F(7:12) doubled): xi2 of body 1 changes by more than 1e-3
  !> of its size, through the off-diagonal body blocks of A and B; with A = B = 0 the same change leaves it as it was
  !> (1e-14), since M, C, K_ext and B_ext are block-diagonal.
  subroutine test_two_bodies_coupled(error)
    type(error_type), allocatable, intent(out) :: error
    character(len=*), parameter :: root = "hams_test_qtf_motion"
    real(8), parameter :: w(2) = [0.6d0, 0.8d0], om = 0.2d0, h50 = 50d0
    real(8), parameter :: m_diag(6) = [4.025d5, 4.025d5, 4.025d5, 2.0d6, 2.0d6, 5.0d6]
    real(8), parameter :: c_diag(6) = [0d0, 0d0, 7.9d5, 1.5d6, 1.5d6, 0d0]
    real(8), parameter :: k_diag(6) = [4d4, 4d4, 0d0, 0d0, 0d0, 1d6]
    real(8), parameter :: b_diag(6) = [2d4, 2d4, 0d0, 0d0, 0d0, 0d0]
    type(so_options) :: opts
    type(panel_mesh) :: mesh
    type(surface_geom) :: geom
    type(first_order_state) :: s
    complex(8), allocatable :: wl(:), psi(:, :)
    real(8), allocatable :: rad(:)
    real(8) :: mass(12, 12), c_hyd(12, 12), k_ext(12, 12), b_ext(12, 12), a(12, 12), b(12, 12), zero(12, 12)
    real(8) :: berr, change, change0, rcond, rho_state
    complex(8) :: xi(12), f(12), lhs(12, 12), fcol(12, 2), xcol(12, 2), xm(6, 2)
    character(len=:), allocatable :: msg
    character(len=250) :: line
    character(len=32) :: os_name
    logical :: ok
    integer :: i, md, bd, estat

    call set_wave(h50)
    NTHREAD = 4
    ISOLV = 1
    call load_multi(0, mesh)
    call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
    call check(error, ok .and. mesh%nbody == 2, "build_surface_geom stopped, or not two bodies")
    if (allocated(error)) return
    call CreateDirectory(root//"/Wamit_format", ok)
    call check(error, ok, "could not create the test folder")
    if (allocated(error)) return
    mass = 0d0
    c_hyd = 0d0
    k_ext = 0d0
    b_ext = 0d0
    zero = 0d0
    do bd = 0, 1
      do md = 1, 6
        mass(6*bd + md, 6*bd + md) = m_diag(md)
        c_hyd(6*bd + md, 6*bd + md) = c_diag(md)
        k_ext(6*bd + md, 6*bd + md) = k_diag(md)
        b_ext(6*bd + md, 6*bd + md) = b_diag(md)
      end do
    end do
    opts%mean_drift = 1; opts%qtf = 1; opts%qtf_terms = QTF_TERMS_BODY; opts%qtf_motion = 1
    call qtf_init(opts, 2, [0d0], geom%nwl, mesh%n, mesh%nbody, root, omegas_main=w, mass=mass, c_hyd=c_hyd, &
                  k_ext=k_ext, b_ext=b_ext)
    do i = 1, 2
      ! Rigid motions that differ per body and frequency
      xm(:, 1) = [(0.3d0, 0.1d0), (-0.2d0, 0.4d0), (0.5d0, -0.1d0), (0.02d0, 0.05d0), (-0.04d0, 0.03d0), (0.06d0, -0.02d0)]
      xm(:, 2) = [(-0.1d0, 0.2d0), (0.3d0, 0d0), (0.2d0, 0.2d0), (-0.03d0, 0.01d0), (0.05d0, 0d0), (0.01d0, 0.04d0)]
      xm = xm*exp((0d0, 0.7d0)*i)
      call make_state(mesh, geom, w(i), 0d0, h50, xm, s, wl)
      call qtf_keep(1, s, wl, w(i))
    end do
    rho_state = s%rho
    rad = qtf_radiation_omegas()
    call check(error, size(rad) == 1 .and. abs(rad(1) - om) <= 1d-12, "the radiation-only list is not [0.2]")
    if (allocated(error)) return
    allocate (psi(mesh%n, 12))
    call radiation_at(om, .true., mesh%nbase)
    call radiation_psi(mesh, om, MXPOT_MULTI_COMB, psi)
    call qtf_keep_psi(om, psi)
    call potential_part(mesh, geom, 1d0)
    call qtf_motion_part(mesh, geom)
    call qtf_motion_load(mesh, geom, KERNEL_DIFF, 1, 2, 1, 1, 1, xi, f, ok)
    call qtf_finish(mesh, geom, 1d0)
    call restore_symmetry()
    call check(error, ok, "Omega = 0.2 is flagged")
    if (allocated(error)) return
    call radiation_coefficients(mesh, psi, om, rho_state, a, b)
    call motion_matrix(om, mass, a, b, b_ext, c_hyd, k_ext, lhs)
    berr = maxval(abs(matmul(lhs, xi) - f))/(maxval(sum(abs(lhs), dim=2))*maxval(abs(xi)))
    ! Only body 2's force changes
    fcol(:, 1) = f
    fcol(:, 2) = f
    fcol(7:12, 2) = 2d0*f(7:12)
    call motion_solve_omega(lhs, fcol, xcol, ok, rcond)
    change = maxval(abs(xcol(1:6, 2) - xcol(1:6, 1)))/maxval(abs(xcol(1:6, 1)))
    call motion_matrix(om, mass, zero, zero, b_ext, c_hyd, k_ext, lhs)
    call motion_solve_omega(lhs, fcol, xcol, ok, rcond)
    change0 = maxval(abs(xcol(1:6, 2) - xcol(1:6, 1)))/maxval(abs(xcol(1:6, 1)))
    write (line, '(a,es10.2,a,es10.2,a,es10.2,a,es10.2)') "two bodies: backward error ", berr, &
      "; body 1 changes by ", change, " with A and B, by ", change0, " without; off-diagonal |A| ", maxval(abs(a(1:6, 7:12)))
    print '(4x,a)', trim(line)
    call get_environment_variable("OS", os_name)
    if (trim(os_name) == "Windows_NT") then
      call execute_command_line("if exist "//root//" rmdir /s /q "//root, exitstat=estat)
    else
      call execute_command_line("rm -rf "//root, exitstat=estat)
    end if
    call check(error, berr <= 1d-13 .and. change > 1d-3 .and. change0 <= 1d-14, line)
  end subroutine test_two_bodies_coupled

end module test_qtf_motion
