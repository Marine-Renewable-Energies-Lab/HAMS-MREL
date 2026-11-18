!> Unit tests of the second-order motion of free bodies (Phase 6c, QTFMotion): the added mass and damping at Omega from
!> the stored radiation potentials psi (radiation_coefficients) against RFORCE and RFORCE_MULTI of the first-order solve,
!> on a half mesh, for reciprocity, at two sum frequencies against the semi-analytical sharp cylinder, and the matrix of
!> the equation of motion (motion_matrix).
!> The radiation solves set module globals (the mesh arrays, the stored Green function, the wave), so the suite runs
!> serially (see unit_tester.f90).
module test_qtf_motion
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use HAMS_mod, only : NTHREAD, ISOLV
  use WaveDyn_mod, only : WK, W1, TP
  use Const_mod, only : RHO
  use PanelMesh_mod, only : ISX, ISY, NSYS
  use Potentials_mod, only : MXPOT, MXPOT_MULTI_COMB
  use PotentWavForce, only : RFORCE
  use PotentWavForceMulti, only : RFORCE_MULTI
  use SecondOrderTypes, only : panel_mesh
  use FirstOrderState, only : radiation_psi
  use QTFMotion, only : radiation_coefficients, motion_matrix
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
      new_unittest("radiation_coefficients_sum_frequency", test_sum_frequency) &
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

end module test_qtf_motion
