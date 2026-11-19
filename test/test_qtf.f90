!> Unit tests of QTFDriver: the stored pairs against the kernel, their symmetries, and the rows of QTF_Quadratic.10d/.10s.
!> The store is module state, so unit_tester runs this suite serially.
module test_qtf
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use Const_mod, only : G_SOLVER => G, RHO_SOLVER => RHO
  use IO, only : CreateDirectory
  use SecondOrderTypes, only : so_options, panel_mesh, first_order_state, surface_geom, cs_points, QTF_TERMS_BODY, &
                               QTF_TERMS_FULL, CS_FREE
  use SurfaceGeometry, only : build_surface_geom, DEFAULT_SHARP_DEG
  use NearFieldDrift, only : near_field_pair, KERNEL_DIFF, KERNEL_SUM
  use QTFDriver, only : qtf_init, qtf_keep, qtf_pair_load, qtf_finish, qtf_kept, qtf_radiation_list, qtf_radiation_omegas, &
                        qtf_direct_omegas, qtf_keep_psi, qtf_nonfinite_psi, qtf_nonfinite_pairs, qtf_psi_tolerance, &
                        qtf_psi_index, qtf_direct_missing, qtf_omega_list, irregular_frequency_estimate, sum_pair_index, &
                        qtf_keep_fs_state, qtf_keep_fs_psi, qtf_fs_state, qtf_fs_psi, qtf_fs_active, qtf_stored_mb, &
                        qtf_direct_sum_pairs, qtf_direct_at_omega, qtf_direct_sum_load, qtf_direct_pairs, &
                        qtf_direct_diff_load, qtf_sum_omega, qtf_diff_omega, qtf_psi_count, qtf_indirect_sum_load, &
                        potential_part_sum, psi_bad_files, potential_part, qtf_motion_part, qtf_motion_load
  use QTFMotion, only : radiation_coefficients, motion_matrix, motion_solve_omega
  use QTFFreeSurface, only : fs_point_field
  use HAMS_mod, only : NTHREAD, ISOLV
  use Potentials_mod, only : MXPOT
  use FirstOrderState, only : radiation_psi
  use QTFDirect, only : direct_setup
  use SharpEdge, only : DEFAULT_EDGE_DEG
  use test_waterline_field, only : set_wave, load_panels
  use test_qtf_potential, only : make_state, motions
  use test_qtf_direct, only : radiation_at, radiation_at_irr, append_lid, set_lid, clear_lid
  use PartialWaves, only : pw_set
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use SecondOrderOutput, only : write_qtf_row, write_motion_row
  use second_order_fixtures, only : make_cylinder, maccamy_fuchs_state, maccamy_fuchs_phi
  implicit none
  private

  public :: collect_qtf

  real(8), parameter :: radius = 5d0, depth = 10d0, g = 9.81d0, rho = 1025d0
  real(8), parameter :: omegas(3) = [0.8d0, 1.0d0, 1.3d0], betas(2) = [0d0, 30d0]
  ! The free-surface store of fs_store_order: 5 inner points, partial waves of M = 2 and N_ev = 1 about (1, 2), R = 7
  integer, parameter :: fs_n = 5, fs_mmax = 2, fs_nev = 1
  real(8), parameter :: fs_centre(2) = [1d0, 2d0], fs_radius = 7d0
  complex(8), parameter :: ci = (0d0, 1d0)
  character(len=*), parameter :: root = "hams_test_qtf_files"
  ! direct_sum_terms_waterline: D_I2 + D_D against P_BODY over max |P_BODY| on the sharp cylinder of radius 5, draft 5,
  ! depth 10, frequencies 0.4 and 0.5 rad/s. Measured at 112 / 448 / 1008 panels: all pairs (the diagonal at 0.8 and
  ! 1.0 rad/s is the largest) 8.76e-2 / 5.89e-2 / 4.54e-2, the pair (0.4, 0.5) 2.43e-2 / 1.62e-2 / 1.24e-2, an order of
  ! about 0.6 in the panel size. The difference kernel on the same runs (qtf = 3, .12d against .11d over max |.11d - .10d|)
  ! falls at the same order: 7.77e-3 / 4.93e-3 / 3.73e-3. So the gap is the mesh error of the waterline body, which both
  ! kernels share, not an error of the sum path (which would not fall). Limits at 448 panels 1.3x measured, rounded up;
  ! the ratio 112 to 448 (measured 1.49 for both) must be at least terms_waterline_ratio
  real(8), parameter :: terms_waterline_limit_all = 7.7d-2, terms_waterline_limit_pair = 2.1d-2
  real(8), parameter :: terms_waterline_ratio = 1.3d0

contains

subroutine collect_qtf(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("qtf_store_matches_kernel", test_store_matches_kernel), &
    new_unittest("qtf_store_symmetries", test_store_symmetries), &
    new_unittest("qtf_row_format", test_row_format), &
    new_unittest("radiation_list_uniform", test_radiation_list_uniform), &
    new_unittest("radiation_list_nonuniform", test_radiation_list_nonuniform), &
    new_unittest("radiation_list_period_grid", test_radiation_list_period_grid), &
    new_unittest("omega_list_union", test_omega_list_union), &
    new_unittest("omega_list_diff_unchanged", test_omega_list_diff_unchanged), &
    new_unittest("omega_list_sum_top", test_omega_list_sum_top), &
    new_unittest("irregular_estimate_cylinder", test_irregular_estimate_cylinder), &
    new_unittest("irregular_estimate_two_loops", test_irregular_estimate_two_loops), &
    new_unittest("irregular_cutoff_applies", test_irregular_cutoff_applies), &
    new_unittest("potential_part_period_grid", test_potential_part_period_grid), &
    new_unittest("psi_store_order", test_psi_store_order), &
    new_unittest("qtf_direct_omegas_all", test_qtf_direct_omegas_all), &
    new_unittest("qtf_direct_missing", test_qtf_direct_missing), &
    new_unittest("psi_nonfinite_counted", test_psi_nonfinite_counted), &
    new_unittest("psi_bad_files_text", test_psi_bad_files_text), &
    new_unittest("fs_store_order", test_fs_store_order), &
    new_unittest("qtf_sum_pairs_cover", test_sum_pairs_cover), &
    new_unittest("qtf_stored_mb_sum", test_stored_mb_sum), &
    new_unittest("radiation_list_sum_indirect", test_radiation_list_sum_indirect), &
    new_unittest("radiation_list_diff_unchanged", test_radiation_list_diff_unchanged), &
    new_unittest("qtf_sum_indirect_list_body", test_sum_indirect_list_body), &
    new_unittest("qtf_stored_mb_indirect_sum", test_stored_mb_indirect_sum), &
    new_unittest("qtf_sum_store_symmetric", test_sum_store_symmetric), &
    new_unittest("direct12s_rows", test_direct12s_rows), &
    new_unittest("direct_sum_terms_waterline", test_direct_sum_terms_waterline), &
    new_unittest("indirect11s_rows", test_indirect11s_rows), &
    new_unittest("indirect11s_equals_tpot", test_indirect11s_equals_tpot), &
    new_unittest("indirect11s_routes_agree", test_indirect11s_routes_agree), &
    new_unittest("qtf3_diff_unchanged", test_qtf3_diff_unchanged), &
    new_unittest("motion_rows_follow_force", test_motion_rows_follow_force), &
    new_unittest("motion_equals_h_times_f", test_motion_equals_h_times_f), &
    new_unittest("motion_scaling", test_motion_scaling), &
    new_unittest("motion_sum_symmetric", test_motion_sum_symmetric), &
    new_unittest("motion_sum_cutoff_rows", test_motion_sum_cutoff_rows), &
    new_unittest("qtf_stored_mb_motion", test_stored_mb_motion), &
    new_unittest("motion_off_unchanged", test_motion_off_unchanged), &
    new_unittest("motion_fixed_no_file", test_motion_fixed_no_file), &
    new_unittest("qtf_files", test_files) &
    ]
end subroutine collect_qtf

real(8) function wave_number(omega, h)
  real(8), intent(in) :: omega, h
  real(8) :: k
  integer :: it
  k = omega**2 / g
  do it = 1, 100
    k = k - (k * tanh(k * h) - omega**2 / g) / (tanh(k * h) + k * h / cosh(k * h)**2)
  end do
  wave_number = k
end function wave_number

!> Cylinder with the exact fixed-body fields at each frequency and heading, and made-up motions, kept in the store.
subroutine fill_store(opts, mesh, geom, states, wls, error, body, fs)
  type(so_options), intent(in) :: opts
  type(panel_mesh), intent(out) :: mesh
  type(surface_geom), intent(out) :: geom
  type(first_order_state), intent(out) :: states(size(omegas), size(betas))
  complex(8), allocatable, intent(out) :: wls(:,:,:)
  type(error_type), allocatable, intent(out) :: error
  logical, intent(in), optional :: body   ! qtf_terms = body: the store also takes psi at the differences of omegas
  ! qtf_terms = full: the store also takes the free-surface fields at these points (fs_code, pw_code), with body
  type(cs_points), intent(in), optional :: fs
  logical :: ok
  character(len=:), allocatable :: msg
  real(8) :: k, b
  integer :: i, ib, j

  call CreateDirectory(root//"/Wamit_format", ok)
  call check(error, ok, "could not create the test folder")
  if (allocated(error)) return
  call make_cylinder(radius, depth, 48, 10, 3, mesh, .true.)
  mesh%xr(:, 1) = [0.5d0, -0.2d0, 0.3d0]
  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "build_surface_geom stopped: "//msg)
  if (allocated(error)) return
  if (present(fs)) then
    call qtf_init(opts, size(omegas), betas, geom%nwl, mesh%n, mesh%nbody, root, omegas_main=omegas, fs_pts=fs, &
                  fs_centre=fs_centre, fs_radius=fs_radius, fs_depth=depth, fs_mmax=fs_mmax, fs_nev=fs_nev)
  else if (present(body)) then
    call qtf_init(opts, size(omegas), betas, geom%nwl, mesh%n, mesh%nbody, root, omegas_main=omegas)
  else
    call qtf_init(opts, size(omegas), betas, geom%nwl, mesh%n, mesh%nbody, root)
  end if
  allocate(wls(geom%nwl, size(omegas), size(betas)))
  do i = 1, size(omegas)
    k = wave_number(omegas(i), depth)
    do ib = 1, size(betas)
      b = betas(ib) * 3.141592653589793d0 / 180d0
      call maccamy_fuchs_state(mesh, radius, depth, k, omegas(i), 1d0, b, g, rho, states(i, ib))
      states(i, ib)%xi(:, 1) = exp(ci * (omegas(i) + b)) * [(0.3d0, 0.1d0), (-0.2d0, 0.4d0), (0.5d0, -0.1d0), &
                                                            (0.02d0, 0.05d0), (-0.04d0, 0.03d0), (0.06d0, -0.02d0)]
      do j = 1, geom%nwl
        wls(j, i, ib) = maccamy_fuchs_phi(geom%wl_xm(:, j), radius, depth, k, omegas(i), 1d0, b, g)
      end do
      call qtf_keep(ib, states(i, ib), wls(:, i, ib), omegas(i))
      if (present(fs)) call qtf_keep_fs_state(ib, fs_codes(fs%n, i, ib), pw_code(i, ib))
    end do
  end do
  call check(error, qtf_kept(), size(omegas), "frequencies kept")
end subroutine fill_store

subroutine test_store_matches_kernel(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: opts
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: states(size(omegas), size(betas))
  complex(8), allocatable :: wls(:,:,:)
  complex(8) :: total(6, 1), terms(6, 6, 1)
  integer :: kn

  opts%mean_drift = 2; opts%qtf = 3
  call fill_store(opts, mesh, geom, states, wls, error)
  if (allocated(error)) return
  do kn = KERNEL_DIFF, KERNEL_SUM
    call qtf_pair_load(mesh, geom, kn, 1, 2, 3, 1, total)
    call near_field_pair(mesh, geom, states(1, 2), states(3, 1), wls(:, 1, 2), wls(:, 3, 1), terms, kn)
    call check(error, maxval(abs(total(:, 1) - sum(terms(:, :, 1), dim=2))) <= 1d-12 * maxval(abs(total)), &
               "stored pair differs from the kernel on the original states")
    if (allocated(error)) return
  end do
  ! The diagonal of the difference kernel is the mean drift
  call qtf_pair_load(mesh, geom, KERNEL_DIFF, 2, 1, 2, 1, total)
  call near_field_pair(mesh, geom, states(2, 1), states(2, 1), wls(:, 2, 1), wls(:, 2, 1), terms)
  call check(error, all(total(:, 1) == sum(terms(:, :, 1), dim=2)), "diagonal is not the mean drift")
  if (allocated(error)) return
  call check(error, maxval(abs(aimag(total))) <= 1d-12 * maxval(abs(total)), "diagonal is not real")
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)   ! frees the store
  call remove_files()
end subroutine test_store_matches_kernel

subroutine test_store_symmetries(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: opts
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: states(size(omegas), size(betas))
  complex(8), allocatable :: wls(:,:,:)
  complex(8) :: t12(6, 1), t21(6, 1)

  opts%mean_drift = 2; opts%qtf = 3
  call fill_store(opts, mesh, geom, states, wls, error)
  if (allocated(error)) return
  call qtf_pair_load(mesh, geom, KERNEL_DIFF, 1, 1, 3, 2, t12)
  call qtf_pair_load(mesh, geom, KERNEL_DIFF, 3, 2, 1, 1, t21)
  call check(error, maxval(abs(t12 - conjg(t21))) <= 1d-12 * maxval(abs(t12)), "T- is not Hermitian")
  if (allocated(error)) return
  call qtf_pair_load(mesh, geom, KERNEL_SUM, 1, 1, 3, 2, t12)
  call qtf_pair_load(mesh, geom, KERNEL_SUM, 3, 2, 1, 1, t21)
  call check(error, maxval(abs(t12 - t21)) <= 1d-12 * maxval(abs(t12)), "T+ is not symmetric")
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)
  call remove_files()
end subroutine test_store_symmetries

subroutine test_row_format(error)
  type(error_type), allocatable, intent(out) :: error
  integer :: unit, idx
  real(8) :: w1, w2, b1, b2, md, ph, re, im, unit_load
  complex(8), parameter :: value = (3d4, -4d4)

  open(newunit=unit, status="scratch")
  call write_qtf_row(unit, 0.8d0, 1.2d0, 0d0, 30d0, 5, value, 2d0, .true.)
  rewind(unit)
  read(unit, *) w1, w2, b1, b2, idx, md, ph, re, im
  close(unit)
  call check(error, w1 == 0.8d0 .and. w2 == 1.2d0 .and. b1 == 0d0 .and. b2 == 30d0 .and. idx == 5, "row keys")
  if (allocated(error)) return
  unit_load = RHO_SOLVER * G_SOLVER * 4d0
  call check(error, abs(cmplx(re, im, 8) - conjg(value) / unit_load) <= 1d-6 * abs(value) / unit_load, &
             "row value is not conj(value)/(rho g refl^2)")
  if (allocated(error)) return
  call check(error, abs(md - abs(value) / unit_load) <= 1d-6 * md .and. &
             abs(ph - atan2(im, re) * 180d0 / 3.141592653589793d0) <= 1d-4, "modulus or phase")
end subroutine test_row_format

!> The 0.01 rad/s grid from 0.40 to 1.40: every difference 0.01 to 1.00 once, those from 0.40 up on the main grid. With
!> qtf_max_dw = 0.3 only 0.01 to 0.30, none of them on the main grid.
subroutine test_radiation_list_uniform(error)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: grid(101), tol
  real(8), allocatable :: need(:)
  integer, allocatable :: reuse(:)
  integer :: i, m

  do i = 1, 101
    grid(i) = 0.40d0 + 0.01d0 * real(i - 1, 8)
  end do
  tol = 1d-6 * maxval(grid)
  call qtf_radiation_list(grid, 0d0, tol, need, reuse)
  call check(error, size(need), 100, "distinct differences, max_dw = 0")
  if (allocated(error)) return
  call check(error, size(reuse), 100, "reuse indices")
  if (allocated(error)) return
  do m = 1, 100
    call check(error, abs(need(m) - 0.01d0 * real(m, 8)) <= 1d-9, "difference values")
    if (allocated(error)) return
    if (m >= 40) then
      call check(error, reuse(m) == m - 39, "a difference on the main grid is not reused")
      if (allocated(error)) return
      call check(error, need(m) == grid(reuse(m)), "a reused difference is not the main-grid value")
    else
      call check(error, reuse(m) == 0, "a difference off the main grid is reused")
    end if
    if (allocated(error)) return
  end do
  call check(error, count(reuse > 0), 61, "reused")
  if (allocated(error)) return
  call check(error, count(reuse == 0), 39, "radiation-only solves")
  if (allocated(error)) return
  call qtf_radiation_list(grid, 0.3d0, tol, need, reuse)
  call check(error, size(need), 30, "distinct differences, max_dw = 0.3")
  if (allocated(error)) return
  call check(error, abs(need(30) - 0.30d0) <= 1d-9 .and. abs(need(1) - 0.01d0) <= 1d-9, "range with max_dw = 0.3")
  if (allocated(error)) return
  call check(error, count(reuse > 0), 0, "reused with max_dw = 0.3")
end subroutine test_radiation_list_uniform

!> A non-uniform grid has its own differences, none on the grid.
subroutine test_radiation_list_nonuniform(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), allocatable :: need(:)
  integer, allocatable :: reuse(:)

  call qtf_radiation_list([0.5d0, 0.7d0, 1.3d0], 0d0, 1.3d-6, need, reuse)
  call check(error, size(need), 3, "distinct differences")
  if (allocated(error)) return
  call check(error, maxval(abs(need - [0.2d0, 0.6d0, 0.8d0])) <= 1d-12, "difference values")
  if (allocated(error)) return
  call check(error, all(reuse == 0), "a difference is reused")
end subroutine test_radiation_list_nonuniform

!> Phase 7c-2: the union of the difference and the sum frequencies of the grid 0.3, 0.6, 0.9. 0.6 is a difference (0.9 - 0.3,
!> 0.6 - 0.3) and a sum (0.3 + 0.3) and is solved once; 0.9 is a grid frequency, a sum (0.3 + 0.6) and is reused.
subroutine test_omega_list_union(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: grid(3) = [0.3d0, 0.6d0, 0.9d0]
  real(8), allocatable :: need(:)
  integer, allocatable :: reuse(:)
  logical, allocatable :: is_diff(:), is_sum(:)

  call qtf_omega_list(grid, 0d0, 0d0, qtf_psi_tolerance(grid), .true., .true., need, reuse, is_diff, is_sum)
  call check(error, size(need), 6, "frequencies")
  if (allocated(error)) return
  call check(error, maxval(abs(need - [0.3d0, 0.6d0, 0.9d0, 1.2d0, 1.5d0, 1.8d0])) <= 1d-12, "values")
  if (allocated(error)) return
  call check(error, all(is_diff .eqv. [.true., .true., .false., .false., .false., .false.]), "is_diff")
  if (allocated(error)) return
  call check(error, all(is_sum .eqv. [.false., .true., .true., .true., .true., .true.]), "is_sum")
  if (allocated(error)) return
  call check(error, all(reuse == [1, 2, 3, 0, 0, 0]), "reuse")
  if (allocated(error)) return
  ! The sum alone: the differences are not wanted
  call qtf_omega_list(grid, 0d0, 0d0, qtf_psi_tolerance(grid), .false., .true., need, reuse, is_diff, is_sum)
  call check(error, size(need), 5, "sums only")
  if (allocated(error)) return
  call check(error, .not. any(is_diff) .and. all(is_sum), "flags of the sums only")
end subroutine test_omega_list_union

!> Without the sum, qtf_omega_list is qtf_radiation_list to the bit, on the three grids of the tests above (uniform,
!> non-uniform, period grid), with qtf_max_dw 0 and 0.3.
subroutine test_omega_list_diff_unchanged(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: two_pi = 6.283185307179586d0
  real(8) :: g1(101), g2(3), g3(105), max_dw(2)
  real(8), allocatable :: need(:), need0(:)
  integer, allocatable :: reuse(:), reuse0(:)
  logical, allocatable :: is_diff(:), is_sum(:)
  integer :: i, ig, im

  do i = 1, 101
    g1(i) = 0.40d0 + 0.01d0 * real(i - 1, 8)
  end do
  g2 = [0.5d0, 0.7d0, 1.3d0]
  do i = 1, 105
    g3(i) = two_pi / (4d0 + 0.25d0 * real(i - 1, 8))
  end do
  max_dw = [0d0, 0.3d0]
  do ig = 1, 3
    do im = 1, 2
      select case (ig)
      case (1)
        call qtf_radiation_list(g1, max_dw(im), qtf_psi_tolerance(g1), need0, reuse0)
        call qtf_omega_list(g1, max_dw(im), 7d0, qtf_psi_tolerance(g1), .true., .false., need, reuse, is_diff, is_sum)
      case (2)
        call qtf_radiation_list(g2, max_dw(im), qtf_psi_tolerance(g2), need0, reuse0)
        call qtf_omega_list(g2, max_dw(im), 7d0, qtf_psi_tolerance(g2), .true., .false., need, reuse, is_diff, is_sum)
      case default
        call qtf_radiation_list(g3, max_dw(im), qtf_psi_tolerance(g3), need0, reuse0)
        call qtf_omega_list(g3, max_dw(im), 7d0, qtf_psi_tolerance(g3), .true., .false., need, reuse, is_diff, is_sum)
      end select
      call check(error, size(need) == size(need0), "size")
      if (allocated(error)) return
      call check(error, all(need == need0) .and. all(reuse == reuse0), "values or reuse differ")
      if (allocated(error)) return
      call check(error, all(is_diff) .and. .not. any(is_sum), "flags")
      if (allocated(error)) return
    end do
  end do
end subroutine test_omega_list_diff_unchanged

!> sum_top cuts the sums: the grid 0.3, 0.6, 0.9 with sum_top = 1.2 keeps the sums up to 1.2 (rounding allowed)
subroutine test_omega_list_sum_top(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: grid(3) = [0.3d0, 0.6d0, 0.9d0]
  real(8), allocatable :: need(:)
  integer, allocatable :: reuse(:)
  logical, allocatable :: is_diff(:), is_sum(:)

  call qtf_omega_list(grid, 0d0, 1.2d0, qtf_psi_tolerance(grid), .true., .true., need, reuse, is_diff, is_sum)
  call check(error, size(need), 4, "frequencies")
  if (allocated(error)) return
  call check(error, maxval(abs(need - [0.3d0, 0.6d0, 0.9d0, 1.2d0])) <= 1d-12, "values")
end subroutine test_omega_list_sum_top

!> The estimate of the lowest irregular frequency of the cylinder of radius 10 m with 32 waterline segments: sqrt(g j01 /
!> a_eq) with a_eq from the area of the polygon, and 1.536 rad/s within 0.5 %.
subroutine test_irregular_estimate_cylinder(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: pi = 3.141592653589793d0, j01 = 2.404826d0, r = 10d0
  integer, parameter :: n = 32
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8) :: area, expected, est
  logical :: ok
  character(len=:), allocatable :: msg

  call make_cylinder(r, 10d0, n, 4, 2, mesh, .true.)
  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "build_surface_geom stopped")
  if (allocated(error)) return
  call check(error, geom%nwl, n, "waterline segments")
  if (allocated(error)) return
  area = 0.5d0 * real(n, 8) * r**2 * sin(2d0 * pi / real(n, 8))
  expected = sqrt(g * j01 / sqrt(area / pi))
  est = irregular_frequency_estimate(geom, g)
  call check(error, abs(est - expected) <= 1d-12 * expected, "value against the polygon area")
  if (allocated(error)) return
  call check(error, abs(est - 1.536d0) <= 0.005d0 * 1.536d0, "within 0.5 % of 1.536")
end subroutine test_irregular_estimate_cylinder

!> One body with two closed loops, radii 6 and 3.25 m (the second one clockwise): the minimum is the 6 m loop, 1.983 rad/s
!> within 0.5 %. A third, open chain of two segments is skipped.
subroutine test_irregular_estimate_two_loops(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: pi = 3.141592653589793d0, j01 = 2.404826d0
  integer, parameter :: n = 64
  type(surface_geom) :: geom
  real(8) :: radius(2), a, b, area, expected
  integer :: l, k, iw

  geom%nwl = 2 * n + 2
  allocate(geom%wl_x1(3, geom%nwl), geom%wl_x2(3, geom%nwl), geom%wl_next(geom%nwl), geom%wl_body(geom%nwl))
  geom%wl_x1 = 0d0
  geom%wl_x2 = 0d0
  geom%wl_body = 1
  radius = [6d0, 3.25d0]
  do l = 1, 2
    do k = 1, n
      iw = (l - 1) * n + k
      a = 2d0 * pi * real(k - 1, 8) / real(n, 8)
      b = 2d0 * pi * real(k, 8) / real(n, 8)
      if (l == 2) then
        a = -a
        b = -b
      end if
      geom%wl_x1(1:2, iw) = radius(l) * [cos(a), sin(a)]
      geom%wl_x2(1:2, iw) = radius(l) * [cos(b), sin(b)]
      geom%wl_next(iw) = (l - 1) * n + mod(k, n) + 1
    end do
  end do
  ! The open chain: a 1 m square corner of two segments, whose "area" would give a far higher value anyway
  geom%wl_x1(:, 2*n + 1) = [0d0, 0d0, 0d0]
  geom%wl_x2(:, 2*n + 1) = [1d0, 0d0, 0d0]
  geom%wl_x1(:, 2*n + 2) = [1d0, 0d0, 0d0]
  geom%wl_x2(:, 2*n + 2) = [1d0, 1d0, 0d0]
  geom%wl_next(2*n + 1) = 2*n + 2
  geom%wl_next(2*n + 2) = 0
  area = 0.5d0 * real(n, 8) * radius(1)**2 * sin(2d0 * pi / real(n, 8))
  expected = sqrt(g * j01 / sqrt(area / pi))
  call check(error, abs(irregular_frequency_estimate(geom, g) - expected) <= 1d-12 * expected, "the 6 m loop")
  if (allocated(error)) return
  call check(error, abs(irregular_frequency_estimate(geom, g) - 1.983d0) <= 0.005d0 * 1.983d0, "within 0.5 % of 1.983")
  if (allocated(error)) return
  ! No closed loop: no cutoff
  geom%wl_next = 0
  call check(error, irregular_frequency_estimate(geom, g) == huge(1d0), "no closed loop")
end subroutine test_irregular_estimate_two_loops

!> IRSP = 0 with the estimate 1.536 rad/s of the cylinder of radius 10 m: the sum 0.8 + 0.8 = 1.6 is left out (0) and
!> 0.7 + 0.8 = 1.5 is kept (positive index); with IRSP = 1 both are kept. qtf_max_sw cuts the same way, inclusively.
subroutine test_irregular_cutoff_applies(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w(2) = [0.7d0, 0.8d0]
  type(so_options) :: opts
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: state
  complex(8), allocatable :: wl(:)
  real(8) :: k, c(2)
  logical :: ok, made
  character(len=:), allocatable :: msg
  integer :: irsp, i

  call CreateDirectory(root//"/Wamit_format", made)
  call check(error, made, "could not create the test folder")
  if (allocated(error)) return
  call make_cylinder(10d0, depth, 32, 4, 2, mesh, .true.)
  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "build_surface_geom stopped")
  if (allocated(error)) return
  opts%mean_drift = 1; opts%qtf = 2; opts%qtf_terms = QTF_TERMS_BODY; opts%qtf_direct = 1
  allocate(wl(geom%nwl))
  wl = (0d0, 0d0)
  do irsp = 0, 1
    call qtf_init(opts, 2, [0d0], geom%nwl, mesh%n, mesh%nbody, root, omegas_main=w, geom=geom, g=g, irsp=irsp)
    do i = 1, 2
      k = wave_number(w(i), depth)
      call maccamy_fuchs_state(mesh, 10d0, depth, k, w(i), 1d0, 0d0, g, rho, state)
      call qtf_keep(1, state, wl, w(i))
    end do
    if (irsp == 0) then
      call check(error, sum_pair_index(2, 2), 0, "0.8 + 0.8 above the estimate is not left out")
      if (allocated(error)) return
      call check(error, sum_pair_index(1, 2) > 0, "0.7 + 0.8 below the estimate is not kept")
      if (allocated(error)) return
      call check(error, size(qtf_direct_omegas()), 2, "the list stops below the estimate: 1.4 and 1.5, not 1.6")
    else
      call check(error, sum_pair_index(2, 2) > 0 .and. sum_pair_index(1, 2) > 0, "IRSP = 1 leaves a sum out")
      if (allocated(error)) return
      call check(error, size(qtf_direct_omegas()), 3, "IRSP = 1: 1.4, 1.5 and 1.6")
    end if
    if (allocated(error)) return
  end do
  ! qtf_max_sw = 1.5 keeps 0.7 + 0.8 (inclusive) and leaves 0.8 + 0.8 out, with IRSP = 1 too
  opts%qtf_max_sw = 1.5d0
  call qtf_init(opts, 2, [0d0], geom%nwl, mesh%n, mesh%nbody, root, omegas_main=w, geom=geom, g=g, irsp=1)
  do i = 1, 2
    k = wave_number(w(i), depth)
    call maccamy_fuchs_state(mesh, 10d0, depth, k, w(i), 1d0, 0d0, g, rho, state)
    call qtf_keep(1, state, wl, w(i))
  end do
  c = [real(sum_pair_index(1, 2), 8), real(sum_pair_index(2, 2), 8)]
  call check(error, c(1) > 0d0 .and. c(2) == 0d0, "qtf_max_sw = 1.5")
  call remove_files()
end subroutine test_irregular_cutoff_applies

!> The period grid of the final review (T = 4 to 30 s in 0.25 s steps, qtf_max_dw = 0.3): every pair finds a difference
!> frequency within qtf_psi_tolerance, the tolerance potential_pair is given, and within the 1e-6 max(omega_i, omega_j) of
!> its default check (the 1e-6 merge of the first version broke this at 52 pairs). Only rounding is merged: the list is
!> strictly ascending, and no two of its values are within the tolerance.
subroutine test_radiation_list_period_grid(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: two_pi = 6.283185307179586d0, max_dw = 0.3d0
  real(8) :: grid(105), tol, dw, best
  real(8), allocatable :: need(:)
  integer, allocatable :: reuse(:)
  integer :: i, j, m, far, pairs

  do i = 1, 105
    grid(i) = two_pi / (4d0 + 0.25d0 * real(i - 1, 8))
  end do
  tol = qtf_psi_tolerance(grid)
  call check(error, abs(tol - 1d-9 * maxval(grid)) <= 1d-24, "tolerance")
  if (allocated(error)) return
  call qtf_radiation_list(grid, max_dw, tol, need, reuse)
  call check(error, all(need(2:) - need(:size(need) - 1) > tol), "differences within the tolerance kept apart")
  if (allocated(error)) return
  far = 0
  pairs = 0
  do i = 1, 105
    do j = 1, 105
      dw = grid(i) - grid(j)
      if (dw <= 0d0 .or. dw > max_dw * (1d0 + 1d-9)) cycle
      pairs = pairs + 1
      best = huge(1d0)
      do m = 1, size(need)
        best = min(best, abs(need(m) - dw))
      end do
      if (best > tol .or. best > 1d-6 * grid(i)) far = far + 1
    end do
  end do
  call check(error, pairs > 3000, "pairs of the grid within qtf_max_dw")
  if (allocated(error)) return
  call check(error, far, 0, "pairs whose difference frequency is not within the tolerance")
end subroutine test_radiation_list_period_grid

!> The potential part on the periods 4.0, 4.6, 7.002 and 9.0738 s (the integration fixture qtf_period_body): two of their
!> differences are 1.1e-6 rad/s apart. The 1e-6 merge of the first version gave both pairs one psi, which potential_pair
!> then refused with a stop in qtf_finish. Now each Omega has its psi, every pair finds it (qtf_psi_index), and qtf_finish
!> writes every pair of QTF_Indirect.11d.
subroutine test_potential_part_period_grid(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: two_pi = 6.283185307179586d0
  real(8), parameter :: periods(4) = [4.0d0, 4.6d0, 7.002d0, 9.0738d0]
  type(so_options) :: opts
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: state
  real(8) :: w(4), k
  real(8), allocatable :: rad(:)
  complex(8), allocatable :: wl(:), psi(:,:)
  logical :: ok
  character(len=:), allocatable :: msg
  integer :: i, j, m, kept

  w = two_pi / periods
  call CreateDirectory(root//"/Wamit_format", ok)
  call check(error, ok, "could not create the test folder")
  if (allocated(error)) return
  call make_cylinder(radius, depth, 48, 10, 3, mesh, .true.)
  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "build_surface_geom stopped: "//msg)
  if (allocated(error)) return
  opts%mean_drift = 1; opts%qtf = 1; opts%qtf_terms = QTF_TERMS_BODY
  call qtf_init(opts, 4, [0d0], geom%nwl, mesh%n, mesh%nbody, root, omegas_main=w)
  allocate(wl(geom%nwl), psi(mesh%n, 6))
  do i = 1, 4
    k = wave_number(w(i), depth)
    call maccamy_fuchs_state(mesh, radius, depth, k, w(i), 1d0, 0d0, g, rho, state)
    state%xi(:, 1) = [(0.3d0, 0.1d0), (-0.2d0, 0.4d0), (0.5d0, -0.1d0), (0.02d0, 0.05d0), (-0.04d0, 0.03d0), &
                      (0.06d0, -0.02d0)]
    do j = 1, geom%nwl
      wl(j) = maccamy_fuchs_phi(geom%wl_xm(:, j), radius, depth, k, w(i), 1d0, 0d0, g)
    end do
    call qtf_keep(1, state, wl, w(i))
  end do
  ! Every pair finds its Omega
  do i = 1, 4
    do j = 1, 4
      if (i == j) cycle
      call check(error, qtf_psi_index(abs(w(i) - w(j))) > 0, "a pair has no difference frequency")
      if (allocated(error)) return
    end do
  end do
  ! psi of every Omega: a made-up value that depends on Omega, from the main grid where it is one, else radiation-only
  kept = 0
  do i = 1, 4
    psi = cmplx(0.3d0 + w(i), -0.1d0, kind=8)
    call qtf_keep_psi(w(i), psi, ok)
    if (ok) kept = kept + 1
  end do
  rad = qtf_radiation_omegas()
  do m = 1, size(rad)
    psi = cmplx(0.3d0 + rad(m), -0.1d0, kind=8)
    call qtf_keep_psi(rad(m), psi)
    kept = kept + 1
  end do
  call check(error, kept, 6, "difference frequencies: all 6 differences of the grid are distinct")
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)
  call check(error, qtf_nonfinite_pairs(), 0, "non-finite pairs")
  if (allocated(error)) return
  call check(error, count_rows(root//"/Wamit_format/QTF_Indirect.11d"), 16 * 6, "rows of .11d")
  call remove_files()
end subroutine test_potential_part_period_grid

!> qtf_direct_omegas gives every difference frequency, those on the main grid too; qtf_radiation_omegas keeps the others.
subroutine test_qtf_direct_omegas_all(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: opts
  real(8), allocatable :: w(:)
  logical :: made

  call CreateDirectory(root//"/Wamit_format", made)
  call check(error, made, "could not create the test folder")
  if (allocated(error)) return
  opts%mean_drift = 1; opts%qtf = 1; opts%qtf_terms = QTF_TERMS_BODY
  call qtf_init(opts, 3, [0d0], 4, 10, 1, root, omegas_main=[0.5d0, 0.7d0, 1.3d0])
  w = qtf_direct_omegas()
  call check(error, size(w), 3, "every difference frequency")
  if (allocated(error)) return
  call check(error, maxval(abs(w - [0.2d0, 0.6d0, 0.8d0])) <= 1d-12, "values")
  if (allocated(error)) return
  ! The grid 0.2, 0.4, 0.6: its differences 0.2 and 0.4 are on the grid
  call qtf_init(opts, 3, [0d0], 4, 10, 1, root, omegas_main=[0.2d0, 0.4d0, 0.6d0])
  w = qtf_direct_omegas()
  call check(error, size(w), 2, "on-grid differences: both")
  if (allocated(error)) return
  call check(error, size(qtf_radiation_omegas()), 0, "radiation-only list is empty on a uniform grid")
  if (allocated(error)) return
  ! qtf_max_dw limits both
  opts%qtf_max_dw = 0.3d0
  call qtf_init(opts, 3, [0d0], 4, 10, 1, root, omegas_main=[0.2d0, 0.4d0, 0.6d0])
  w = qtf_direct_omegas()
  call check(error, size(w), 1, "qtf_max_dw = 0.3")
  call remove_files()
end subroutine test_qtf_direct_omegas_all

!> qtf_direct_missing, the check of qtf_finish before QTF_Direct.12d is written (Phase 7a Task 4; its stop cannot be reached
!> from input, since RadiationOnlySolves solves every Omega of the list): before any direct solve it names the smallest
!> difference frequency of the grid, and without qtf_direct = 1 it is 0.
subroutine test_qtf_direct_missing(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: opts
  logical :: made

  call CreateDirectory(root//"/Wamit_format", made)
  call check(error, made, "could not create the test folder")
  if (allocated(error)) return
  opts%mean_drift = 1; opts%qtf = 1; opts%qtf_terms = QTF_TERMS_BODY; opts%qtf_direct = 1
  call qtf_init(opts, 3, [0d0], 4, 10, 1, root, omegas_main=[0.5d0, 0.7d0, 1.3d0])
  call check(error, abs(qtf_direct_missing() - 0.2d0) <= 1d-12, "no direct solve yet: the first Omega is missing")
  if (allocated(error)) return
  opts%qtf_direct = 0
  call qtf_init(opts, 3, [0d0], 4, 10, 1, root, omegas_main=[0.5d0, 0.7d0, 1.3d0])
  call check(error, qtf_direct_missing() == 0d0, "without qtf_direct nothing is missing")
  call remove_files()
end subroutine test_qtf_direct_missing

!> The psi store takes the differences of the grid only, at the size of the mesh.
subroutine test_psi_store_order(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: opts
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  complex(8) :: psi(10, 6), wrong(9, 6)
  real(8), allocatable :: w(:)
  logical :: ok, made

  call CreateDirectory(root//"/Wamit_format", made)
  call check(error, made, "could not create the test folder")
  if (allocated(error)) return
  opts%mean_drift = 1; opts%qtf = 1; opts%qtf_terms = QTF_TERMS_BODY
  call qtf_init(opts, 3, [0d0], 4, 10, 1, root, omegas_main=[0.5d0, 0.7d0, 1.3d0])
  w = qtf_radiation_omegas()
  call check(error, size(w), 3, "radiation-only frequencies")
  if (allocated(error)) return
  call check(error, maxval(abs(w - [0.2d0, 0.6d0, 0.8d0])) <= 1d-12, "radiation-only values")
  if (allocated(error)) return
  psi = (1d0, 2d0)
  call qtf_keep_psi(0.33d0, psi, ok)
  call check(error, .not. ok, "an omega that is not a difference was kept")
  if (allocated(error)) return
  wrong = (0d0, 0d0)
  call qtf_keep_psi(0.6d0, wrong, ok)
  call check(error, .not. ok, "psi of the wrong size was kept")
  if (allocated(error)) return
  call qtf_keep_psi(0.6d0, psi, ok)
  call check(error, ok, "psi at a difference was not kept")
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)   ! frees the store
  call remove_files()
end subroutine test_psi_store_order

!> The files named by the non-finite psi message: the terms files only when the run writes them (write_terms = 0 names
!> neither QTF_PotentialTerms.txt nor QTF_PotentialTermsSum.txt), the sum files for a sum frequency, both sets for a frequency
!> that is both
subroutine test_psi_bad_files_text(error)
  type(error_type), allocatable, intent(out) :: error

  call check(error, trim(psi_bad_files(.false., .false., .false., .false.)), "QTF_Indirect.11d", "diff, no terms")
  if (allocated(error)) return
  call check(error, trim(psi_bad_files(.false., .false., .true., .true.)), "QTF_Indirect.11d and QTF_PotentialTerms.txt", &
             "diff with terms")
  if (allocated(error)) return
  call check(error, trim(psi_bad_files(.true., .false., .false., .false.)), "QTF_Indirect.11s", "sum, no terms")
  if (allocated(error)) return
  call check(error, trim(psi_bad_files(.true., .false., .false., .true.)), "QTF_Indirect.11s and QTF_PotentialTermsSum.txt", &
             "sum with terms")
  if (allocated(error)) return
  call check(error, trim(psi_bad_files(.true., .true., .false., .false.)), "QTF_Indirect.11d and QTF_Indirect.11s", &
             "both, no terms")
  if (allocated(error)) return
  call check(error, trim(psi_bad_files(.true., .true., .true., .true.)), "QTF_Indirect.11d, QTF_PotentialTerms.txt, "// &
             "QTF_Indirect.11s and QTF_PotentialTermsSum.txt", "both with terms")
  if (allocated(error)) return
  call check(error, index(psi_bad_files(.true., .true., .false., .false.), "Terms") == 0, "terms file named at write_terms = 0")
end subroutine test_psi_bad_files_text

!> A psi with a NaN is kept but counted (qtf_nonfinite_psi), and the pairs of its Omega are counted by the final guard
!> (qtf_nonfinite_pairs): 0.5 rad/s = 1.3 - 0.8 has the pairs (1, 3) and (3, 1) at each of the 2 headings. psi at 0.2 and
!> 0.3 rad/s is finite and adds nothing. The main program turns either count into exit code 1.
subroutine test_psi_nonfinite_counted(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: opts
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: states(size(omegas), size(betas))
  complex(8), allocatable :: wls(:,:,:), psi(:,:)
  real(8), allocatable :: w(:)
  integer :: m

  opts%mean_drift = 1; opts%qtf = 1; opts%qtf_terms = QTF_TERMS_BODY; opts%write_terms = 1
  call fill_store(opts, mesh, geom, states, wls, error, body=.true.)
  if (allocated(error)) return
  w = qtf_radiation_omegas()
  call check(error, size(w), 3, "difference frequencies")
  if (allocated(error)) return
  allocate(psi(mesh%n, 6))
  do m = 1, 3
    psi = (0.4d0, -0.1d0)
    if (abs(w(m) - 0.5d0) < 1d-9) psi(7, 3) = cmplx(ieee_value(1d0, ieee_quiet_nan), 0d0, kind=8)
    call qtf_keep_psi(w(m), psi)
  end do
  call check(error, qtf_nonfinite_psi(), 1, "difference frequencies with a non-finite psi")
  if (allocated(error)) return
  call check(error, qtf_nonfinite_pairs(), 0, "non-finite pairs before qtf_finish")
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)
  call check(error, qtf_nonfinite_pairs(), 2 * 2, "frequency and heading pairs with a non-finite T_pot")
  if (allocated(error)) return
  call check(error, qtf_nonfinite_psi(), 1, "the psi count is kept after qtf_finish")
  call remove_files()
end subroutine test_psi_nonfinite_counted

!> The fields kept for inner point p of frequency i and heading ib: distinct values per point, frequency and heading.
function fs_codes(n, i, ib) result(f)
  integer, intent(in) :: n, i, ib
  type(fs_point_field) :: f(n)
  integer :: p, c
  do p = 1, n
    f(p)%phi = cmplx(100 * i + 10 * ib + p, -p, 8)
    do c = 1, 3
      f(p)%grad(c) = cmplx(1000 * c + 100 * i + 10 * ib + p, 0.5d0 * c, 8)
    end do
    f(p)%phizz = cmplx(-(100 * i + 10 * ib + p), 7d0 * i, 8)
  end do
end function fs_codes

!> The partial waves kept for frequency i and heading ib (code = 10 i + ib), or for psi_k at a difference frequency
!> (code = 1000 m + k): the sizes of the store, distinct coefficients.
function pw_code(i, ib) result(pw)
  integer, intent(in) :: i, ib
  type(pw_set) :: pw
  integer :: m, n
  pw%mmax = fs_mmax
  pw%nev = fs_nev
  ! k0 and kev are those that qtf_finish needs for the free-surface integral (Task 8): the wave number of the state for
  ! a state (i <= size(omegas)), and the same k0 and kev for every mode of a psi (i > 100), so that the outer integral
  ! accepts the sets; the amplitudes tell the frequencies, headings and modes apart
  if (i <= size(omegas)) then
    pw%k0 = wave_number(omegas(i), depth)
  else
    pw%k0 = 0.001d0 * i
  end if
  pw%h = depth
  pw%radius = fs_radius
  pw%centre = fs_centre
  allocate (pw%kev(fs_nev), pw%a(-fs_mmax:fs_mmax, 0:fs_nev))
  pw%kev = 0.3d0 + 0.001d0 * i
  do n = 0, fs_nev
    do m = -fs_mmax, fs_mmax
      pw%a(m, n) = cmplx(10 * i + ib + 0.1d0 * m, n - 0.01d0 * m, 8)
    end do
  end do
end function pw_code

!> True if two pw_sets hold the same modes, to the bit.
logical function same_pw(a, b)
  type(pw_set), intent(in) :: a, b
  same_pw = .false.
  if (.not. (allocated(a%a) .and. allocated(b%a) .and. allocated(a%kev) .and. allocated(b%kev))) return
  if (a%mmax /= b%mmax .or. a%nev /= b%nev .or. a%k0 /= b%k0 .or. a%h /= b%h .or. a%radius /= b%radius) return
  if (any(a%centre /= b%centre) .or. any(a%kev /= b%kev)) return
  if (any(shape(a%a) /= shape(b%a))) return
  same_pw = all(a%a == b%a)
end function same_pw

!> qtf_terms = full: the free-surface store keeps the fields at the inner points (5 complex values each) and the partial
!> waves of every kept frequency and heading, and psi_k at the inner points with the partial waves of each psi_k at every
!> difference frequency, and gives them back by frequency, heading and Omega, as qtf_store_matches_kernel for the states.
!> An Omega that is not a difference frequency is refused, and so are partial waves on another circle or depth. The size formula of qtf_stored_mb: 16 bytes x (5 n_fs + n_pw)
!> per frequency and heading, plus 16 x 6 nbody x (n_fs + n_pw) per Omega, n_pw = (2M + 1)(N_ev + 1). With
!> qtf_terms = body the store has no free-surface part.
subroutine test_fs_store_order(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: opts
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: states(size(omegas), size(betas))
  complex(8), allocatable :: wls(:,:,:), psi(:,:), psi_pts(:,:), psi_back(:,:)
  type(cs_points) :: fs
  type(fs_point_field), allocatable :: f_back(:), f_ref(:)
  type(pw_set) :: pw_back
  type(pw_set), allocatable :: pw_psi(:), pw_psi_back(:)
  real(8), allocatable :: w(:)
  real(8) :: mb0, mb1, expected
  logical :: ok, same
  integer :: i, ib, m, k, p, npw

  fs%n = fs_n
  allocate (fs%x(3, fs_n), fs%nrm(3, fs_n), fs%w(fs_n), fs%psize(fs_n), fs%kind(fs_n), fs%body(fs_n))
  do p = 1, fs_n
    fs%x(:, p) = [6d0 + p, 0.5d0 * p, 0d0]
  end do
  fs%nrm = 0d0
  fs%nrm(3, :) = -1d0
  fs%w = 1d0
  fs%psize = 1d0
  fs%kind = CS_FREE
  fs%body = 0

  ! The size formula
  npw = (2 * fs_mmax + 1) * (fs_nev + 1)
  mb0 = qtf_stored_mb(10, 4, 1, 3, 2, npsi=3, npairs=4)
  mb1 = qtf_stored_mb(10, 4, 1, 3, 2, npsi=3, npairs=4, nfs=fs_n, npw=npw)
  expected = 16d0 * (5 * fs_n + npw) * 3 * 2 / 1d6 + 16d0 * 6 * (fs_n + npw) * 3 / 1d6
  call check(error, abs((mb1 - mb0) - expected) <= 1d-12 * expected, "free-surface share of qtf_stored_mb")
  if (allocated(error)) return

  ! qtf_terms = body: no free-surface store
  opts%mean_drift = 2; opts%qtf = 1; opts%qtf_terms = QTF_TERMS_BODY; opts%write_terms = 0
  call fill_store(opts, mesh, geom, states, wls, error, body=.true.)
  if (allocated(error)) return
  call check(error, .not. qtf_fs_active(), "qtf_terms = body has a free-surface store")
  if (allocated(error)) return
  call keep_all_psi()
  call qtf_finish(mesh, geom, 1d0)

  opts%qtf_terms = QTF_TERMS_FULL
  call fill_store(opts, mesh, geom, states, wls, error, body=.true., fs=fs)
  if (allocated(error)) return
  call check(error, qtf_fs_active(), "qtf_terms = full has no free-surface store")
  if (allocated(error)) return
  ! The fields by frequency and heading
  same = .true.
  do i = 1, size(omegas)
    do ib = 1, size(betas)
      call qtf_fs_state(i, ib, f_back, pw_back)
      f_ref = fs_codes(fs_n, i, ib)
      if (size(f_back) /= fs_n) then
        same = .false.
        cycle
      end if
      do p = 1, fs_n
        if (f_back(p)%phi /= f_ref(p)%phi .or. any(f_back(p)%grad /= f_ref(p)%grad) .or. &
            f_back(p)%phizz /= f_ref(p)%phizz) same = .false.
      end do
      if (.not. same_pw(pw_back, pw_code(i, ib))) same = .false.
    end do
  end do
  call check(error, same, "the free-surface fields do not come back by frequency and heading")
  if (allocated(error)) return

  ! psi_k at the difference frequencies, and a frequency that is not one
  w = qtf_radiation_omegas()
  allocate (psi_pts(fs_n, 6), pw_psi(6))
  call qtf_keep_fs_psi(0.33d0, psi_pts, pw_psi, ok)
  call check(error, .not. ok, "psi at an omega that is not a difference frequency was kept")
  if (allocated(error)) return
  ! Partial waves on another circle or depth than the store's (review M3)
  do k = 1, 6
    pw_psi(k) = pw_code(100, k)
  end do
  pw_psi(3)%radius = fs_radius * (1d0 + 1d-9)
  call qtf_keep_fs_psi(w(1), psi_pts, pw_psi, ok)
  call check(error, .not. ok, "psi with partial waves on another radius was kept")
  if (allocated(error)) return
  pw_psi(3) = pw_code(100, 3)
  pw_psi(4)%centre(2) = fs_centre(2) + 1d-6
  call qtf_keep_fs_psi(w(1), psi_pts, pw_psi, ok)
  call check(error, .not. ok, "psi with partial waves about another centre was kept")
  if (allocated(error)) return
  pw_psi(4) = pw_code(100, 4)
  pw_psi(5)%h = depth + 1d-3
  call qtf_keep_fs_psi(w(1), psi_pts, pw_psi, ok)
  call check(error, .not. ok, "psi with partial waves of another depth was kept")
  if (allocated(error)) return
  do m = 1, size(w)
    do k = 1, 6
      do p = 1, fs_n
        psi_pts(p, k) = cmplx(m, 10 * k + p, 8)
      end do
      pw_psi(k) = pw_code(100 + m, k)
    end do
    call qtf_keep_fs_psi(w(m), psi_pts, pw_psi, ok)
    call check(error, ok, "psi at a difference frequency was not kept")
    if (allocated(error)) return
  end do
  do m = size(w), 1, -1
    call qtf_fs_psi(w(m), psi_back, pw_psi_back, ok)
    call check(error, ok, "psi at a difference frequency does not come back")
    if (allocated(error)) return
    same = size(psi_back, 1) == fs_n .and. size(psi_back, 2) == 6 .and. size(pw_psi_back) == 6
    if (same) then
      do k = 1, 6
        do p = 1, fs_n
          if (psi_back(p, k) /= cmplx(m, 10 * k + p, 8)) same = .false.
        end do
        if (.not. same_pw(pw_psi_back(k), pw_code(100 + m, k))) same = .false.
      end do
    end if
    call check(error, same, "psi at the inner points or its partial waves do not come back by Omega")
    if (allocated(error)) return
  end do
  call keep_all_psi()
  call qtf_finish(mesh, geom, 1d0)
  call check(error, .not. qtf_fs_active(), "qtf_finish does not free the free-surface store")
  call remove_files()

contains

  subroutine keep_all_psi()
    ! psi on the hull at every difference frequency, so that qtf_finish can form the potential part
    real(8), allocatable :: wd(:)
    integer :: md
    wd = qtf_radiation_omegas()
    allocate (psi(mesh%n, 6))
    psi = (0.4d0, -0.1d0)
    do md = 1, size(wd)
      call qtf_keep_psi(wd(md), psi)
    end do
    deallocate (psi)
  end subroutine keep_all_psi

end subroutine test_fs_store_order

!> Phase 7c-2, R6: on the grid [0.3, 0.6, 0.9] with two headings and mean_drift = 2, the canonical sum pairs of every
!> frequency of the list (qtf_direct_sum_pairs) hold each unordered pair of states (i, b1, j, b2) exactly once: 6 states,
!> 6 x 7 / 2 = 21 pairs, the diagonal ones (i, b, i, b) included. With the difference kernel on (qtf = 3) the list also has
!> the difference frequencies, and 0.6 is both (0.9 - 0.3 and 0.3 + 0.3): the count is the same.
subroutine test_sum_pairs_cover(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w(3) = [0.3d0, 0.6d0, 0.9d0], hb(2) = [0d0, 30d0]
  type(so_options) :: opts
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: state
  complex(8), allocatable :: wl(:)
  integer, allocatable :: pairs(:, :)
  real(8), allocatable :: om(:)
  integer :: seen(3, 2, 3, 2), i, ib, m, p, a(2), b(2), kq, n
  logical :: ok
  character(len=:), allocatable :: msg

  call CreateDirectory(root//"/Wamit_format", ok)
  call check(error, ok, "could not create the test folder")
  if (allocated(error)) return
  call make_cylinder(10d0, depth, 32, 4, 2, mesh, .true.)
  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "build_surface_geom stopped")
  if (allocated(error)) return
  allocate (wl(geom%nwl))
  wl = (0d0, 0d0)
  do kq = 2, 3
    opts%mean_drift = 2; opts%qtf = kq; opts%qtf_terms = QTF_TERMS_BODY; opts%qtf_direct = 1
    call qtf_init(opts, 3, hb, geom%nwl, mesh%n, mesh%nbody, root, omegas_main=w)
    do i = 1, 3
      do ib = 1, 2
        call maccamy_fuchs_state(mesh, 10d0, depth, wave_number(w(i), depth), w(i), 1d0, hb(ib)*0.0174533d0, g, rho, state)
        call qtf_keep(ib, state, wl, w(i))
      end do
    end do
    seen = 0
    n = 0
    om = qtf_direct_omegas()
    do m = 1, size(om)
      call qtf_direct_sum_pairs(om(m), pairs)
      do p = 1, size(pairs, 2)
        ! The unordered pair, as its two states sorted
        a = pairs(1:2, p)
        b = pairs(3:4, p)
        if (b(1) < a(1) .or. (b(1) == a(1) .and. b(2) < a(2))) then
          a = pairs(3:4, p)
          b = pairs(1:2, p)
        end if
        seen(a(1), a(2), b(1), b(2)) = seen(a(1), a(2), b(1), b(2)) + 1
        n = n + 1
        call check(error, abs(w(pairs(1, p)) + w(pairs(3, p)) - om(m)) <= 1d-12, "a pair is listed at another Omega")
        if (allocated(error)) return
      end do
    end do
    call check(error, n, 21, "canonical sum pairs over all Omega")
    if (allocated(error)) return
    call check(error, maxval(seen), 1, "an unordered pair is listed twice")
    if (allocated(error)) return
    call check(error, count(seen == 1), 21, "unordered pairs listed once")
    if (allocated(error)) return
  end do
  ! No qtf_finish: psi was not kept, which the potential parts need; the next qtf_init resets the store
end subroutine test_sum_pairs_cover

!> qtf_stored_mb with the sum stores (Phase 7c-2): tdir_sum and tpot_sum add 16 bytes x 6 nbody x nfreq^2 per slice,
!> nsum_pairs = 2 x the heading pairs; and the inner mesh built for the largest sum frequency adds 16 (5 dn_fs) bytes per
!> frequency and heading and 16 x 6 nbody dn_fs per frequency of psi for its dn_fs more points.
subroutine test_stored_mb_sum(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: nfreq = 3, nb = 2, nhp = 4, nbody = 1, npsi = 7, npw = 20, nfs0 = 10, nfs1 = 130
  real(8) :: mb0, mb1, expected

  mb0 = qtf_stored_mb(10, 4, nbody, nfreq, nb, npsi=npsi, npairs=0)
  mb1 = qtf_stored_mb(10, 4, nbody, nfreq, nb, npsi=npsi, npairs=0, nsum_pairs=2*nhp)
  expected = 2d0 * 16d0 * 6d0 * nbody * nfreq**2 * nhp / 1d6
  call check(error, abs((mb1 - mb0) - expected) <= 1d-12 * expected, "the sum stores in qtf_stored_mb")
  if (allocated(error)) return
  mb0 = qtf_stored_mb(10, 4, nbody, nfreq, nb, npsi=npsi, nfs=nfs0, npw=npw, nsum_pairs=2*nhp)
  mb1 = qtf_stored_mb(10, 4, nbody, nfreq, nb, npsi=npsi, nfs=nfs1, npw=npw, nsum_pairs=2*nhp)
  expected = 16d0 * (5d0 * (nfs1 - nfs0) * nfreq * nb + 6d0 * nbody * (nfs1 - nfs0) * npsi) / 1d6
  call check(error, abs((mb1 - mb0) - expected) <= 1d-12 * expected, "the larger inner mesh in qtf_stored_mb")
end subroutine test_stored_mb_sum

!> Phase 6b: qtf = 2 with qtf_direct = 0 keeps the sums of the list. Grid 0.3, 0.6, 0.9: the sums are 0.6, 0.9, 1.2, 1.2, 1.5,
!> 1.8; 0.6 and 0.9 are main-grid solves, so the radiation-only list is 1.2, 1.5, 1.8. qtf_direct_omegas still stops at
!> the direct solve (qtf_direct = 0 has none), but the list of psi is the same with and without qtf_direct, qtf = 3 adds the
!> differences 0.3 and 0.6 (both on the grid, so the radiation-only list does not change).
subroutine test_radiation_list_sum_indirect(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: opts
  real(8), allocatable :: w(:)
  integer :: kq, kd
  logical :: made

  call CreateDirectory(root//"/Wamit_format", made)
  call check(error, made, "could not create the test folder")
  if (allocated(error)) return
  do kq = 2, 3
    do kd = 0, 1
      opts%mean_drift = 1; opts%qtf = kq; opts%qtf_terms = QTF_TERMS_BODY; opts%qtf_direct = kd
      call qtf_init(opts, 3, [0d0], 4, 10, 1, root, omegas_main=[0.3d0, 0.6d0, 0.9d0])
      w = qtf_radiation_omegas()
      call check(error, size(w), 3, "radiation-only list size")
      if (allocated(error)) return
      call check(error, maxval(abs(w - [1.2d0, 1.5d0, 1.8d0])) <= 1d-12, "radiation-only list values")
      if (allocated(error)) return
    end do
  end do
  call remove_files()
end subroutine test_radiation_list_sum_indirect

!> qtf = 1 has no sum: with qtf_direct = 0 and 1 the lists are the difference lists of the grids of the earlier tests
subroutine test_radiation_list_diff_unchanged(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: opts
  real(8), allocatable :: w(:)
  integer :: kd
  logical :: made

  call CreateDirectory(root//"/Wamit_format", made)
  call check(error, made, "could not create the test folder")
  if (allocated(error)) return
  do kd = 0, 1
    opts%mean_drift = 1; opts%qtf = 1; opts%qtf_terms = QTF_TERMS_BODY; opts%qtf_direct = kd
    call qtf_init(opts, 3, [0d0], 4, 10, 1, root, omegas_main=[0.5d0, 0.7d0, 1.3d0])
    w = qtf_radiation_omegas()
    call check(error, size(w), 3, "difference list size")
    if (allocated(error)) return
    call check(error, maxval(abs(w - [0.2d0, 0.6d0, 0.8d0])) <= 1d-12, "difference list values")
    if (allocated(error)) return
    ! Nothing of the sum: no sum frequency has a psi
    call check(error, .not. qtf_sum_omega(1.0d0), "qtf = 1 has a sum frequency")
    if (allocated(error)) return
  end do
  call remove_files()
end subroutine test_radiation_list_diff_unchanged

!> qtf = 2, qtf_terms = body, qtf_direct = 0 (the path with the free-surface warning of sum_warnings): the sum frequencies are
!> flagged as sums, and no difference frequency is in the list
subroutine test_sum_indirect_list_body(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: opts
  logical :: made

  call CreateDirectory(root//"/Wamit_format", made)
  call check(error, made, "could not create the test folder")
  if (allocated(error)) return
  opts%mean_drift = 1; opts%qtf = 2; opts%qtf_terms = QTF_TERMS_BODY; opts%qtf_direct = 0
  call qtf_init(opts, 3, [0d0], 4, 10, 1, root, omegas_main=[0.3d0, 0.6d0, 0.9d0])
  call check(error, qtf_sum_omega(1.2d0) .and. qtf_sum_omega(0.6d0) .and. qtf_sum_omega(1.8d0), "sum frequencies flagged")
  if (allocated(error)) return
  call check(error, .not. qtf_diff_omega(1.2d0) .and. .not. qtf_diff_omega(0.6d0), "a difference flagged with qtf = 2")
  if (allocated(error)) return
  call check(error, qtf_psi_count(), 5, "5 distinct sums: 0.6, 0.9, 1.2, 1.5, 1.8")
  if (allocated(error)) return
  call remove_files()
end subroutine test_sum_indirect_list_body

!> qtf_stored_mb counts the stores of the sum by nsum_pairs; without the direct solve tpot_sum is the only one, 1 x the
!> heading pairs: 16 bytes x 6 nbody x nfreq^2 per slice
subroutine test_stored_mb_indirect_sum(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: nfreq = 3, nb = 2, nhp = 4, nbody = 1, npsi = 7
  real(8) :: mb0, mb1, mb2, expected

  mb0 = qtf_stored_mb(10, 4, nbody, nfreq, nb, npsi=npsi, npairs=0)
  mb1 = qtf_stored_mb(10, 4, nbody, nfreq, nb, npsi=npsi, npairs=0, nsum_pairs=nhp)
  mb2 = qtf_stored_mb(10, 4, nbody, nfreq, nb, npsi=npsi, npairs=0, nsum_pairs=2*nhp)
  expected = 16d0 * 6d0 * nbody * nfreq**2 * nhp / 1d6
  call check(error, abs((mb1 - mb0) - expected) <= 1d-12 * expected, "tpot_sum alone in qtf_stored_mb")
  if (allocated(error)) return
  call check(error, abs((mb2 - mb1) - expected) <= 1d-12 * expected, "tdir_sum is the second store")
end subroutine test_stored_mb_indirect_sum

!> The whole direct sum path of the driver on a waterline body (Phase 7c-2): the 112-panel cylinder of
!> direct_zero_forcing (radius 5, draft 5, depth 10), frequencies 0.4 and 0.5 rad/s, headings 0 and 30 degrees,
!> mean_drift = 2, qtf = 2, qtf_terms = body, qtf_direct = 1, write_terms = 1. States: the incident wave and rigid motions
!> (make_state). At every frequency of the list the radiation solve (radiation_at), psi into the store, and
!> qtf_direct_at_omega, as RadiationOnlySolves and so_after_radiation_only do it. max_sw is qtf_max_sw.
subroutine sum_pipeline(error, max_sw, mesh, geom, refine, kind_qtf, grid, lid, indirect_only, motion, fixed, irsp0)
  type(error_type), allocatable, intent(out) :: error
  real(8), intent(in) :: max_sw
  type(panel_mesh), intent(out) :: mesh
  type(surface_geom), intent(out) :: geom
  integer, intent(in), optional :: refine     ! 2: the cylinder with twice the panels in each direction (448)
  integer, intent(in), optional :: kind_qtf   ! qtf of SecondOrder.in, 2 by default
  real(8), intent(in), optional :: grid(:)    ! the frequencies, [0.4, 0.5] by default
  logical, intent(in), optional :: lid        ! IRSP = 1 with the lid of the cylinder (clear_lid after the test)
  logical, intent(in), optional :: indirect_only   ! qtf_direct = 0: psi of the radiation-only list, no direct solve
  ! qtf_motion = 1 (Phase 6c) with the matrices of cylinder_body: 1 the free body (no spring), 2 moored; 0 (default) off
  integer, intent(in), optional :: motion
  logical, intent(in), optional :: fixed      ! motion_state = fixed (with motion, no matrices are passed)
  logical, intent(in), optional :: irsp0      ! IRSP = 0 with the geometry: the sums at or above the irregular estimate are cut
  real(8) :: mass(6, 6), c_hyd(6, 6), k_ext(6, 6), b_ext(6, 6)
  integer :: rf, kq, nl, kd, kmot
  real(8), parameter :: hb(2) = [0d0, 30d0], h10 = 10d0
  real(8), allocatable :: w(:), pts(:, :)
  type(panel_mesh) :: full
  type(so_options) :: opts
  type(first_order_state) :: s
  complex(8), allocatable :: wl(:), psi(:, :)
  complex(8) :: x1(6, 1), x2(6, 1), xi(6, 1)
  real(8), allocatable :: om(:)
  integer :: i, ib, m
  logical :: ok, use_lid, use_irsp0
  character(len=:), allocatable :: msg

  call CreateDirectory(root//"/Wamit_format", ok)
  call check(error, ok, "could not create the test folder")
  if (allocated(error)) return
  rf = 1
  if (present(refine)) rf = refine
  kq = 2
  if (present(kind_qtf)) kq = kind_qtf
  if (present(grid)) then
    w = grid
  else
    w = [0.4d0, 0.5d0]
  end if
  use_lid = .false.
  if (present(lid)) use_lid = lid
  kd = 1
  if (present(indirect_only)) kd = merge(0, 1, indirect_only)
  call make_cylinder(5d0, 5d0, 16 * rf, 4 * rf, 3 * rf, full, .true.)
  call set_wave(h10)
  NTHREAD = 4
  ISOLV = 1
  call load_panels(full, 0, 0, mesh)
  nl = 0
  if (use_lid) then
    call append_lid(5d0, 16, 3, 0d0, 0, pts)
    call set_lid(pts, .false.)
    nl = size(pts, 2)
  end if
  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok .and. geom%nwl > 0, "build_surface_geom stopped, or no waterline")
  if (allocated(error)) return
  opts%mean_drift = 2; opts%qtf = kq; opts%qtf_terms = QTF_TERMS_BODY; opts%qtf_direct = kd; opts%write_terms = 1
  opts%qtf_max_sw = max_sw
  kmot = 0
  if (present(motion)) kmot = motion
  opts%qtf_motion = merge(1, 0, kmot > 0)
  if (present(fixed)) opts%motion_fixed = fixed
  use_irsp0 = .false.
  if (present(irsp0)) use_irsp0 = irsp0
  if (kmot > 0 .and. .not. opts%motion_fixed) call cylinder_body(kmot == 1, mass, c_hyd, k_ext, b_ext)
  if (use_irsp0) then
    ! Only with the motion on (motion_sum_cutoff_rows)
    call qtf_init(opts, size(w), hb, geom%nwl, mesh%n, mesh%nbody, root, omegas_main=w, geom=geom, g=G_SOLVER, irsp=0, &
                  mass=mass, c_hyd=c_hyd, k_ext=k_ext, b_ext=b_ext)
  else if (kmot > 0 .and. .not. opts%motion_fixed) then
    call qtf_init(opts, size(w), hb, geom%nwl, mesh%n, mesh%nbody, root, omegas_main=w, mass=mass, c_hyd=c_hyd, &
                  k_ext=k_ext, b_ext=b_ext)
  else
    call qtf_init(opts, size(w), hb, geom%nwl, mesh%n, mesh%nbody, root, omegas_main=w)
  end if
  call motions(x1, x2)
  do i = 1, size(w)
    do ib = 1, 2
      ! Motions that differ per frequency and heading, so that a slot error shows
      xi = merge(x1, x2, mod(i, 2) == 1) * exp(ci * 0.4d0 * ib)
      call make_state(mesh, geom, w(i), hb(ib) * 3.141592653589793d0 / 180d0, h10, xi, s, wl)
      call qtf_keep(ib, s, wl, w(i))
    end do
  end do
  if (kd == 1) then
    call direct_setup(mesh, DEFAULT_EDGE_DEG)
    om = qtf_direct_omegas()
  else
    om = qtf_radiation_omegas()
  end if
  allocate (psi(mesh%n, 6))
  do m = 1, size(om)
    if (use_lid) then
      call radiation_at_irr(om(m), .false., mesh%nbase, nl)
    else
      call radiation_at(om(m), .false., mesh%nbase)
    end if
    call radiation_psi(mesh, om(m), MXPOT, psi)
    call qtf_keep_psi(om(m), psi)
    if (kd == 1) call qtf_direct_at_omega(mesh, geom, om(m))
  end do
end subroutine sum_pipeline

!> Plan global constraint, review of 5e4b8da: with qtf = 3 the sum columns of an Omega that is also a difference frequency
!> must not change the difference results. On the grid [0.3, 0.6, 0.9], Omega = 0.6 is both a difference (0.9 - 0.3) and
!> a sum (0.3 + 0.3), so its direct solve has both kinds of columns with qtf = 3. With two headings and mean_drift = 2,
!> tdir and its two parts of every difference pair are the same to the bit with qtf = 1 and qtf = 3, without the lid and
!> with it (IRSP = 1, the least-squares solve of HAMS with A^T r). The fixture has no field points (qtf_field = 0).
subroutine test_qtf3_diff_unchanged(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w3(3) = [0.3d0, 0.6d0, 0.9d0]
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  complex(8) :: t(6, 3, 3, 2, 2, 2), p1(6, 3, 3, 2, 2, 2), p2(6, 3, 3, 2, 2, 2)
  integer, allocatable :: pd(:, :), ps(:, :)
  integer :: ilid, ik, i, j, b1, b2, nmixed
  character(len=200) :: line

  do ilid = 1, 2
    nmixed = 0
    do ik = 1, 2
      call sum_pipeline(error, 0d0, mesh, geom, kind_qtf=merge(1, 3, ik == 1), grid=w3, lid=ilid == 2)
      if (allocated(error)) then
        call clear_lid()
        return
      end if
      if (ik == 2) then
        call qtf_direct_pairs(0.6d0, pd)
        call qtf_direct_sum_pairs(0.6d0, ps)
        nmixed = min(size(pd, 2), size(ps, 2))
      end if
      do i = 1, 3
        do j = 1, 3
          do b1 = 1, 2
            do b2 = 1, 2
              call qtf_direct_diff_load(i, b1, j, b2, t(:, i, j, b1, b2, ik), p1(:, i, j, b1, b2, ik), &
                                        p2(:, i, j, b1, b2, ik))
            end do
          end do
        end do
      end do
      call qtf_finish(mesh, geom, 1d0)
      call remove_files()
    end do
    call clear_lid()
    write (line, '(a,i0,a,i0,a,es10.2,a,es10.2)') "IRSP = ", ilid - 1, ": Omega = 0.6 has difference and sum pairs ", &
      nmixed, "; max |T_dir| ", maxval(abs(t(:, :, :, :, :, 1))), ", max |qtf 3 - qtf 1| ", &
      maxval(abs(t(:, :, :, :, :, 2) - t(:, :, :, :, :, 1)))
    print '(4x,a)', trim(line)
    call check(error, nmixed > 0 .and. maxval(abs(t(:, :, :, :, :, 1))) > 0d0, "the fixture has no mixed Omega: "//trim(line))
    if (allocated(error)) return
    call check(error, all(t(:, :, :, :, :, 2) == t(:, :, :, :, :, 1)) .and. all(p1(:, :, :, :, :, 2) == p1(:, :, :, :, :, 1)) &
               .and. all(p2(:, :, :, :, :, 2) == p2(:, :, :, :, :, 1)), "qtf = 3 changes the difference results: "//trim(line))
    if (allocated(error)) return
  end do
end subroutine test_qtf3_diff_unchanged

!> R6: after qtf_direct_at_omega, tdir_sum(:, i, j, h(b1, b2)) = tdir_sum(:, j, i, h(b2, b1)) to the bit for every pair, and
!> no slot is left empty (both halves are filled from the canonical one).
subroutine test_sum_store_symmetric(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  complex(8) :: t1(6), t2(6)
  real(8) :: smallest
  integer :: i, j, b1, b2

  call sum_pipeline(error, 0d0, mesh, geom)
  if (allocated(error)) return
  smallest = huge(1d0)
  do i = 1, 2
    do j = 1, 2
      do b1 = 1, 2
        do b2 = 1, 2
          call qtf_direct_sum_load(i, b1, j, b2, t1)
          call qtf_direct_sum_load(j, b2, i, b1, t2)
          call check(error, all(t1 == t2), "tdir_sum(i, b1, j, b2) is not tdir_sum(j, b2, i, b1) to the bit")
          if (allocated(error)) return
          smallest = min(smallest, maxval(abs(t1)))
        end do
      end do
    end do
  end do
  print '(4x,a,es10.2)', "qtf_sum_store_symmetric: smallest max |T+_dir| over the 16 slots ", smallest
  call check(error, smallest > 0d0, "a slot of tdir_sum is empty")
  call qtf_finish(mesh, geom, 1d0)
  call remove_files()
end subroutine test_sum_store_symmetric

!> QTF_Direct.12s with qtf_max_sw = 0.9: the rows of .10s with omega_i + omega_j <= 0.9, in the order of .10s (the pair
!> (0.5, 0.5) is cut: 3 of 4 frequency pairs, 4 heading pairs, 6 modes = 72 rows of the 96 of .10s). The diagonal
!> (0.4, 0.4) is a row of its own and carries T+_dir (it is not .10s, which the diagonal of .12d is).
subroutine test_direct12s_rows(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8), allocatable :: rs(:, :), rd(:, :)
  integer :: i, k
  real(8) :: diag_change

  call sum_pipeline(error, 0.9d0, mesh, geom)
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)
  call read_qtf_rows(root//"/Wamit_format/QTF_Quadratic.10s", rs)
  call read_qtf_rows(root//"/Wamit_format/QTF_Direct.12s", rd)
  call check(error, size(rs, 2), 96, "rows of .10s")
  if (allocated(error)) return
  call check(error, size(rd, 2), 72, "rows of .12s with qtf_max_sw = 0.9")
  if (allocated(error)) return
  k = 0
  diag_change = 0d0
  do i = 1, size(rs, 2)
    if (rs(1, i) + rs(2, i) > 0.9d0 + 1d-9) cycle
    k = k + 1
    call check(error, all(rd(1:5, k) == rs(1:5, i)), ".12s is not in the order of .10s")
    if (allocated(error)) return
    if (rs(1, i) == 0.4d0 .and. rs(2, i) == 0.4d0) diag_change = max(diag_change, maxval(abs(rd(8:9, k) - rs(8:9, i))))
  end do
  call check(error, k, 72, "rows of .10s kept by the cut")
  if (allocated(error)) return
  print '(4x,a,es10.2)', "direct12s_rows: largest |.12s - .10s| on the diagonal (0.4, 0.4) ", diag_change
  call check(error, diag_change > 0d0, "the diagonal of .12s is .10s: T+_dir left out")
  call qtf_finish(mesh, geom, 1d0)
  call remove_files()
end subroutine test_direct12s_rows

!> The in-run indirect check on a body with a waterline (review minors of Tasks 2 and 3: the only direct-against-indirect
!> sum test so far used a submerged sphere and the diagonal pair). QTF_DirectTermsSum.txt has 11 columns (keys, D_I2, D_D,
!> P_BODY) and one row per row of .12s, the non-diagonal pair (0.4, 0.5) and its mirror included. Its parts add up to
!> .12s - .10s (to the 7 digits of the files), and D_I2 + D_D against P_BODY, over the largest |P_BODY|, is the mesh error
!> of the two methods: below the limits at 448 panels, and falling from 112 to 448 panels (limits and measured values at
!> terms_waterline_limit_all).
subroutine test_direct_sum_terms_waterline(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8), allocatable :: rs(:, :), rd(:, :), rt(:, :)
  real(8) :: v(11), gap_all(2), gap_off(2), pmax, pmax_off, sumerr, scale, dpart(2), ppart(2)
  integer :: i, k, ncol, rfn
  character(len=300) :: line
  character(len=*), parameter :: terms_path = root//"/Wamit_format/QTF_DirectTermsSum.txt"

  do rfn = 1, 2
    call sum_pipeline(error, 0d0, mesh, geom, refine=rfn)
    if (allocated(error)) return
    call qtf_finish(mesh, geom, 1d0)
    call read_qtf_rows(root//"/Wamit_format/QTF_Quadratic.10s", rs)
    call read_qtf_rows(root//"/Wamit_format/QTF_Direct.12s", rd)
    call check(error, size(rd, 2), size(rs, 2), ".12s has every row of .10s without a cut")
    if (allocated(error)) return
    call read_terms_rows(terms_path, rt, ncol)
    call check(error, ncol, 11, "columns of QTF_DirectTermsSum.txt with qtf_terms = body")
    if (allocated(error)) return
    call check(error, size(rt, 2), size(rd, 2), "rows of QTF_DirectTermsSum.txt")
    if (allocated(error)) return
    gap_all(rfn) = 0d0
    gap_off(rfn) = 0d0
    pmax = 0d0
    pmax_off = 0d0
    sumerr = 0d0
    scale = maxval(abs(rd(8:9, :)))
    do k = 1, size(rt, 2)
      v = rt(:, k)
      ! The row of .12s and .10s with the same keys
      do i = 1, size(rd, 2)
        if (all(rd(1:5, i) == v(1:5))) exit
      end do
      call check(error, i <= size(rd, 2), "a row of QTF_DirectTermsSum.txt has no row in .12s")
      if (allocated(error)) return
      dpart = v(6:7) + v(8:9)
      ppart = v(10:11)
      sumerr = max(sumerr, maxval(abs((rd(8:9, i) - rs(8:9, i)) - dpart)))
      gap_all(rfn) = max(gap_all(rfn), maxval(abs(dpart - ppart)))
      pmax = max(pmax, maxval(abs(ppart)))
      if (v(1) /= v(2)) then
        gap_off(rfn) = max(gap_off(rfn), maxval(abs(dpart - ppart)))
        pmax_off = max(pmax_off, maxval(abs(ppart)))
      end if
    end do
    gap_all(rfn) = gap_all(rfn) / pmax
    gap_off(rfn) = gap_off(rfn) / pmax_off
    write (line, '(a,i0,a,es10.3,a,es10.3,a,es10.2,a,es10.2,a)') "mesh ", mesh%n, " panels: D_I2 + D_D against P_BODY"// &
      " over max |P_BODY|: all pairs ", gap_all(rfn), ", the pair (0.4, 0.5) ", gap_off(rfn), &
      "; parts against .12s - .10s ", sumerr, " (max |.12s| ", scale, ")"
    print '(4x,a)', trim(line)
    call check(error, sumerr <= 2d-5 * scale .and. pmax_off > 0d0, &
               "the parts of QTF_DirectTermsSum.txt do not add up to .12s - .10s: "//trim(line))
    if (allocated(error)) return
    call remove_files()
  end do
  write (line, '(a,2es10.3,a,es10.3,a,3es10.3)') "the gap falls with the mesh (all, pair; 112 to 448 panels, ratio) ", &
    gap_all(1) / gap_all(2), gap_off(1) / gap_off(2), "; at 448 panels ", maxval([gap_all(2), gap_off(2)]), &
    "; limits ", terms_waterline_limit_all, terms_waterline_limit_pair, terms_waterline_ratio
  print '(4x,a)', trim(line)
  call check(error, gap_all(2) <= terms_waterline_limit_all .and. gap_off(2) <= terms_waterline_limit_pair .and. &
             gap_all(1) / gap_all(2) >= terms_waterline_ratio .and. gap_off(1) / gap_off(2) >= terms_waterline_ratio, line)
end subroutine test_direct_sum_terms_waterline

!> QTF_Indirect.11s with qtf_max_sw = 0.9 and qtf_direct = 0 (Phase 6b): the rows of .10s whose sum is kept, in the order of
!> .10s, as .12s (72 of the 96: the pair (0.5, 0.5) is cut). The diagonal (0.4, 0.4) is a row of its own and carries T+_pot
!> (it is not .10s). The run has no direct solve and writes no .12s.
subroutine test_indirect11s_rows(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8), allocatable :: rs(:, :), ri(:, :)
  integer :: i, k
  real(8) :: diag_change
  logical :: has12s

  call sum_pipeline(error, 0.9d0, mesh, geom, indirect_only=.true.)
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)
  call read_qtf_rows(root//"/Wamit_format/QTF_Quadratic.10s", rs)
  call read_qtf_rows(root//"/Wamit_format/QTF_Indirect.11s", ri)
  inquire (file=root//"/Wamit_format/QTF_Direct.12s", exist=has12s)
  call check(error, .not. has12s, "qtf_direct = 0 wrote QTF_Direct.12s")
  if (allocated(error)) return
  call check(error, size(rs, 2), 96, "rows of .10s")
  if (allocated(error)) return
  call check(error, size(ri, 2), 72, "rows of .11s with qtf_max_sw = 0.9")
  if (allocated(error)) return
  k = 0
  diag_change = 0d0
  do i = 1, size(rs, 2)
    if (rs(1, i) + rs(2, i) > 0.9d0 + 1d-9) cycle
    k = k + 1
    call check(error, all(ri(1:5, k) == rs(1:5, i)), ".11s is not in the order of .10s")
    if (allocated(error)) return
    if (rs(1, i) == 0.4d0 .and. rs(2, i) == 0.4d0) diag_change = max(diag_change, maxval(abs(ri(8:9, k) - rs(8:9, i))))
  end do
  call check(error, k, 72, "rows of .10s kept by the cut")
  if (allocated(error)) return
  print '(4x,a,es10.2)', "indirect11s_rows: largest |.11s - .10s| on the diagonal (0.4, 0.4) ", diag_change
  call check(error, diag_change > 0d0, "the diagonal of .11s is .10s: T+_pot left out")
  call remove_files()
end subroutine test_indirect11s_rows

!> .11s - .10s is T+_pot of the row: tpot_sum as the driver holds it (qtf_indirect_sum_load), conjugated and scaled as the
!> writer does (refl = 1: 1/(rho g) for forces and moments), against the rows of the files to their 7 digits, the diagonal
!> included. The terms file (16 digits) has 15 columns with qtf_terms = body (5 keys and 5 terms of 2 values) and its terms
!> add up to the same T+_pot to 1e-14 of the largest term of the row (rowmax), for every pair and its mirror. (The brief said
!> 15 and 11 columns, which does not fit its own list of 5 keys and 5 terms; ruling in the report, as in Phase 7c-2.)
subroutine test_indirect11s_equals_tpot(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8), parameter :: w(2) = [0.4d0, 0.5d0], hb(2) = [0d0, 30d0]
  complex(8) :: tp(6, 2, 2, 2, 2), t
  real(8), allocatable :: rs(:, :), ri(:, :), rt(:, :)
  real(8) :: scale, err_file, err_terms, tmax, sig(2), rowmax
  integer :: i, j, b1, b2, n, ncol, ia, ja, ba, bb, ks, c

  call sum_pipeline(error, 0d0, mesh, geom, indirect_only=.true.)
  if (allocated(error)) return
  ! tpot_sum is computed at qtf_finish; computing it here too only fills the same slots again, to the bit
  call potential_part_sum(mesh, geom, 1d0)
  do i = 1, 2
    do j = 1, 2
      do b1 = 1, 2
        do b2 = 1, 2
          call qtf_indirect_sum_load(i, b1, j, b2, tp(:, i, j, b1, b2))
        end do
      end do
    end do
  end do
  tmax = maxval(abs(tp)) / (RHO_SOLVER * G_SOLVER)
  call check(error, tmax > 0d0, "tpot_sum is empty without the direct solve")
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)
  call read_qtf_rows(root//"/Wamit_format/QTF_Quadratic.10s", rs)
  call read_qtf_rows(root//"/Wamit_format/QTF_Indirect.11s", ri)
  call read_terms_rows(root//"/Wamit_format/QTF_PotentialTermsSum.txt", rt, ncol)
  call check(error, size(ri, 2), size(rs, 2), ".11s has every row of .10s without a cut")
  if (allocated(error)) return
  call check(error, ncol, 15, "columns of QTF_PotentialTermsSum.txt with qtf_terms = body")
  if (allocated(error)) return
  call check(error, size(rt, 2), size(ri, 2), "rows of QTF_PotentialTermsSum.txt")
  if (allocated(error)) return
  scale = maxval(abs(ri(8:9, :)))
  err_file = 0d0
  do n = 1, size(ri, 2)
    call locate(ri(:, n), ia, ja, ba, bb, ks)
    t = conjg(tp(ks, ia, ja, ba, bb)) / (RHO_SOLVER * G_SOLVER)
    ! the matching row of .10s
    do i = 1, size(rs, 2)
      if (all(rs(1:5, i) == ri(1:5, n))) exit
    end do
    call check(error, i <= size(rs, 2), "a row of .11s has no row in .10s")
    if (allocated(error)) return
    err_file = max(err_file, abs(ri(8, n) - rs(8, i) - real(t, 8)), abs(ri(9, n) - rs(9, i) - aimag(t)))
  end do
  err_terms = 0d0
  do n = 1, size(rt, 2)
    call locate(rt(:, n), ia, ja, ba, bb, ks)
    t = conjg(tp(ks, ia, ja, ba, bb)) / (RHO_SOLVER * G_SOLVER)
    sig = 0d0
    do c = 0, 4
      sig = sig + rt(6 + 2 * c:7 + 2 * c, n)
    end do
    ! per row: the error over the largest of the row's terms and of T+_pot, not over the global maximum
    rowmax = max(maxval(abs(rt(6:15, n))), abs(real(t, 8)), abs(aimag(t)))
    if (rowmax > 0d0) err_terms = max(err_terms, abs(sig(1) - real(t, 8)) / rowmax, abs(sig(2) - aimag(t)) / rowmax)
  end do
  print '(4x,a,es10.2,a,es10.2,a,es10.2)', "indirect11s_equals_tpot: max |.11s - .10s - T+_pot| ", err_file, &
    " (max |.11s| ", scale, "), terms file against T+_pot ", err_terms
  call check(error, err_file <= 3d-6 * scale, ".11s - .10s is not T+_pot of tpot_sum")
  if (allocated(error)) return
  call check(error, err_terms <= 1d-14, "the terms do not add up to T+_pot")
  call remove_files()

contains

  !> The frequency, heading and mode indices of a row whose first five columns are PER1 PER2 BETA1 BETA2 I
  subroutine locate(row, i1, i2, h1, h2, kk)
    real(8), intent(in) :: row(:)
    integer, intent(out) :: i1, i2, h1, h2, kk

    i1 = merge(1, 2, row(1) == w(1))
    i2 = merge(1, 2, row(2) == w(1))
    h1 = merge(1, 2, row(3) == hb(1))
    h2 = merge(1, 2, row(4) == hb(1))
    kk = nint(row(5))
  end subroutine locate

end subroutine test_indirect11s_equals_tpot

!> R2 of Phase 6b: the indirect sum files of the run without the direct solve (qtf_direct = 0) and of the run with it
!> (qtf_direct = 1) are the same to the byte: QTF_Quadratic.10s, QTF_Indirect.11s and QTF_PotentialTermsSum.txt, with
!> qtf_max_sw = 0.9. Both runs solve the radiation problems of the same list, so psi, and with it every T+_pot, must be
!> the same; the run with the direct solve also has QTF_Direct.12s (checked to exist).
subroutine test_indirect11s_routes_agree(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  character(len=:), allocatable :: a1, a2, a3, b1, b2, b3
  logical :: has12s

  call sum_pipeline(error, 0.9d0, mesh, geom, indirect_only=.true.)
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)
  a1 = slurp(root//"/Wamit_format/QTF_Quadratic.10s")
  a2 = slurp(root//"/Wamit_format/QTF_Indirect.11s")
  a3 = slurp(root//"/Wamit_format/QTF_PotentialTermsSum.txt")
  call remove_files()
  call sum_pipeline(error, 0.9d0, mesh, geom)
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)
  b1 = slurp(root//"/Wamit_format/QTF_Quadratic.10s")
  b2 = slurp(root//"/Wamit_format/QTF_Indirect.11s")
  b3 = slurp(root//"/Wamit_format/QTF_PotentialTermsSum.txt")
  inquire (file=root//"/Wamit_format/QTF_Direct.12s", exist=has12s)
  call remove_files()
  call check(error, has12s, "qtf_direct = 1 wrote no QTF_Direct.12s")
  if (allocated(error)) return
  call check(error, len(a2) > 0 .and. len(a3) > 0, "empty .11s or terms file")
  if (allocated(error)) return
  print '(4x,a,2i8)', "indirect11s_routes_agree: bytes of .11s and of QTF_PotentialTermsSum.txt ", len(a2), len(a3)
  call check(error, a1 == b1, ".10s differs between qtf_direct = 0 and 1")
  if (allocated(error)) return
  call check(error, a2 == b2, ".11s differs between qtf_direct = 0 and 1")
  if (allocated(error)) return
  call check(error, a3 == b3, "QTF_PotentialTermsSum.txt differs between qtf_direct = 0 and 1")
end subroutine test_indirect11s_routes_agree

!> Phase 6c: the body matrices of the cylinder of sum_pipeline (radius 5, draft 5): the mass rho pi r^2 T in surge, sway and
!> heave, inertias in roll, pitch and yaw, the hydrostatic restoring in heave, roll and pitch, and a linear damper in surge
!> and sway. Moored (free = .false.) adds springs in surge, sway and yaw, so that C + K_ext is regular at Omega = 0; the
!> free body has none there, as an unmoored hull.
subroutine cylinder_body(free, mass, c_hyd, k_ext, b_ext)
  logical, intent(in) :: free
  real(8), intent(out) :: mass(6, 6), c_hyd(6, 6), k_ext(6, 6), b_ext(6, 6)
  real(8), parameter :: m_diag(6) = [4.025d5, 4.025d5, 4.025d5, 2.0d6, 2.0d6, 5.0d6]
  real(8), parameter :: c_diag(6) = [0d0, 0d0, 7.9d5, 1.5d6, 1.5d6, 0d0]
  real(8), parameter :: k_diag(6) = [4d4, 4d4, 0d0, 0d0, 0d0, 1d6]
  real(8), parameter :: b_diag(6) = [2d4, 2d4, 0d0, 0d0, 0d0, 0d0]
  integer :: md

  mass = 0d0
  c_hyd = 0d0
  k_ext = 0d0
  b_ext = 0d0
  do md = 1, 6
    mass(md, md) = m_diag(md)
    c_hyd(md, md) = c_diag(md)
    b_ext(md, md) = b_diag(md)
    if (.not. free) k_ext(md, md) = k_diag(md)
  end do
end subroutine cylinder_body

!> The rows of a QTF_Motion file, r(9, n) as read_qtf_rows, and its header (the lines that start with #), joined
subroutine read_motion_rows(path, r, header)
  character(len=*), intent(in) :: path
  real(8), allocatable, intent(out) :: r(:, :)
  character(len=:), allocatable, intent(out) :: header
  character(len=1000) :: line
  real(8) :: v(8)
  integer :: unit, ios, n, idx
  logical :: exists

  header = ""
  allocate (r(9, 0))
  inquire (file=path, exist=exists)
  if (.not. exists) return
  n = 0
  open (newunit=unit, file=path, status="old", action="read")
  do
    read (unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0) cycle
    if (line(1:1) == "#") then
      header = header//trim(line)//" "
    else
      n = n + 1
    end if
  end do
  deallocate (r)
  allocate (r(9, n))
  rewind (unit)
  n = 0
  do
    read (unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0 .or. line(1:1) == "#") cycle
    n = n + 1
    read (line, *) v(1:4), idx, v(5:8)
    r(1:4, n) = v(1:4)
    r(5, n) = idx
    r(6:9, n) = v(5:8)
  end do
  close (unit)
end subroutine read_motion_rows

!> Phase 6c: QTF_Motion.11d (and with qtf_direct = 1 QTF_Motion.12d) of the free cylinder has the rows of QTF_Indirect.11d
!> (QTF_Direct.12d) in their order, minus those of Omega = 0: C + K_ext has no surge, sway and yaw, so rcond < 1e-12 there
!> (R6) and the header names Omega = 0 and the modes 1 2 6. Grid 0.4, 0.5 rad/s, two headings, mean_drift = 2: 96 rows
!> of the force file, 48 of them on the diagonal.
subroutine test_motion_rows_follow_force(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8), allocatable :: rf(:, :), rm(:, :), rm_ind(:, :)
  character(len=:), allocatable :: header
  character(len=*), parameter :: force(2) = [character(len=16) :: "QTF_Indirect.11d", "QTF_Direct.12d"]
  character(len=*), parameter :: motion(2) = [character(len=16) :: "QTF_Motion.11d", "QTF_Motion.12d"]
  integer :: kd, f, i, k
  real(8) :: route_gap

  do kd = 1, 2
    call sum_pipeline(error, 0d0, mesh, geom, kind_qtf=1, indirect_only=kd == 1, motion=1)
    if (allocated(error)) return
    call qtf_finish(mesh, geom, 1d0)
    do f = 1, kd
      call read_qtf_rows(root//"/Wamit_format/"//trim(force(f)), rf)
      call read_motion_rows(root//"/Wamit_format/"//trim(motion(f)), rm, header)
      call check(error, size(rf, 2), 96, "rows of "//trim(force(f)))
      if (allocated(error)) return
      call check(error, size(rm, 2), 48, "rows of "//trim(motion(f))//" (the diagonal left out)")
      if (allocated(error)) return
      k = 0
      do i = 1, size(rf, 2)
        if (rf(1, i) == rf(2, i)) cycle
        k = k + 1
        call check(error, all(rm(1:5, k) == rf(1:5, i)), trim(motion(f))//" is not in the order of "//trim(force(f)))
        if (allocated(error)) return
      end do
      call check(error, index(header, "Omega =   0.00000 rad/s") > 0 .and. index(header, "modes 1 2 6") > 0, &
                 trim(motion(f))//": the header does not name Omega = 0 and the modes 1 2 6: "//header)
      if (allocated(error)) return
      call check(error, index(header, trim(force(f))) > 0, trim(motion(f))//": the header does not name its force file")
      if (allocated(error)) return
    end do
    if (kd == 2) then
      ! The two routes have their own force, so their own motion
      call read_motion_rows(root//"/Wamit_format/QTF_Motion.11d", rm_ind, header)
      route_gap = maxval(abs(rm(8:9, :) - rm_ind(8:9, :)))
      print '(4x,a,es10.2,a,es10.2)', "motion_rows_follow_force: max |QTF_Motion.12d - .11d| ", route_gap, &
        ", max |QTF_Motion.11d| ", maxval(abs(rm_ind(8:9, :)))
      call check(error, route_gap > 0d0, "QTF_Motion.12d is QTF_Motion.11d: the direct force is not used")
      if (allocated(error)) return
    end if
    call remove_files()
  end do
end subroutine test_motion_rows_follow_force

!> Phase 6c: xi2 of a pair solves the equation of motion with the force of its row (R4, R6). Moored cylinder, qtf = 3,
!> grid 0.4, 0.5: the difference pairs at Omega = +0.1 and -0.1 (the conjugate system of R4) and at Omega = 0 (C + K_ext
!> only), and two sum pairs at 0.9 and 1.0 rad/s. Unscaled: the backward error |lhs xi - F| / (|lhs| |xi|) (infinity
!> norms) of the stored xi and force is rounding, 1e-13, with lhs from the radiation solve of the test at |Omega|.
!> Rescaled: the rows of QTF_Motion.11d/.11s and QTF_Indirect.11d/.11s, written with refl = 2 and taken back to N and m
!> (force times rho g L^k, motion over L^k, conjugated), solve it to the 7 digits of the files (1e-5). Also prints rcond
!> of each lhs, as is and equilibrated (rows, then columns, to unit largest entry), for parked item J.
subroutine test_motion_equals_h_times_f(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: refl = 2d0
  integer, parameter :: ncase = 5
  integer, parameter :: kernels(ncase) = [KERNEL_DIFF, KERNEL_DIFF, KERNEL_DIFF, KERNEL_SUM, KERNEL_SUM]
  integer, parameter :: pairs(4, ncase) = reshape([2, 1, 1, 2,  1, 2, 2, 1,  1, 1, 1, 2,  1, 1, 2, 2,  2, 2, 2, 1], [4, ncase])
  real(8), parameter :: w(2) = [0.4d0, 0.5d0], hb(2) = [0d0, 30d0]
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8) :: mass(6, 6), c_hyd(6, 6), k_ext(6, 6), b_ext(6, 6), zero(6, 6), a(6, 6), b(6, 6)
  complex(8) :: xi(6), f(6, 1), lhs(6, 6), xt(6, 1), lhs_eq(6, 6), xi_file(6), f_file(6)
  complex(8), allocatable :: psi(:, :), xis(:, :), fs(:, :), lhss(:, :, :)
  real(8), allocatable :: rm(:, :), rf(:, :)
  character(len=:), allocatable :: header
  real(8) :: om, berr, berr_file, worst, worst_file, rc, rc_eq, s
  logical :: ok, okt
  integer :: p, md, n, kk
  character(len=250) :: line

  call sum_pipeline(error, 0d0, mesh, geom, kind_qtf=3, indirect_only=.true., motion=2)
  if (allocated(error)) return
  call cylinder_body(.false., mass, c_hyd, k_ext, b_ext)
  zero = 0d0
  ! The parts qtf_finish computes, here before it, so that the stores can be read (qtf_finish computes them again)
  call potential_part(mesh, geom, refl)
  call potential_part_sum(mesh, geom, refl)
  call qtf_motion_part(mesh, geom)
  allocate (xis(6, ncase), fs(6, ncase), lhss(6, 6, ncase), psi(mesh%n, 6))
  worst = 0d0
  do p = 1, ncase
    call qtf_motion_load(mesh, geom, kernels(p), 1, pairs(1, p), pairs(2, p), pairs(3, p), pairs(4, p), xi, f(:, 1), ok)
    call check(error, ok, "a pair of a regular Omega is flagged")
    if (allocated(error)) return
    if (kernels(p) == KERNEL_DIFF) then
      om = w(pairs(1, p)) - w(pairs(3, p))
    else
      om = w(pairs(1, p)) + w(pairs(3, p))
    end if
    if (om == 0d0) then
      call motion_matrix(0d0, mass, zero, zero, b_ext, c_hyd, k_ext, lhs)
    else
      call radiation_at(abs(om), .false., mesh%nbase)
      call radiation_psi(mesh, abs(om), MXPOT, psi)
      call radiation_coefficients(mesh, psi, abs(om), rho, a, b)
      call motion_matrix(abs(om), mass, a, b, b_ext, c_hyd, k_ext, lhs)
      ! Omega < 0: the matrix of (2) at -|Omega| with A, B of |Omega| is the conjugate (R4)
      if (om < 0d0) lhs = conjg(lhs)
    end if
    berr = maxval(abs(matmul(lhs, xi) - f(:, 1)))/(maxval(sum(abs(lhs), dim=2))*maxval(abs(xi)))
    worst = max(worst, berr)
    call motion_solve_omega(lhs, f, xt, okt, rc)
    ! Equilibrated: rows to unit largest entry, then columns
    lhs_eq = lhs
    do md = 1, 6
      lhs_eq(md, :) = lhs_eq(md, :)/maxval(abs(lhs_eq(md, :)))
    end do
    do md = 1, 6
      lhs_eq(:, md) = lhs_eq(:, md)/maxval(abs(lhs_eq(:, md)))
    end do
    call motion_solve_omega(lhs_eq, f, xt, okt, rc_eq)
    write (line, '(a,i0,a,f6.2,a,es10.2,a,es10.3,a,es10.3)') "kernel ", kernels(p), ", Omega ", om, &
      ": backward error ", berr, ", rcond ", rc, ", equilibrated ", rc_eq
    print '(4x,a)', trim(line)
    xis(:, p) = xi
    fs(:, p) = f(:, 1)
    lhss(:, :, p) = lhs
  end do
  call qtf_finish(mesh, geom, refl)
  ! Rescaled: the rows of the files at the keys of each pair
  worst_file = 0d0
  do p = 1, ncase
    if (kernels(p) == KERNEL_DIFF) then
      call read_motion_rows(root//"/Wamit_format/QTF_Motion.11d", rm, header)
      call read_qtf_rows(root//"/Wamit_format/QTF_Indirect.11d", rf)
    else
      call read_motion_rows(root//"/Wamit_format/QTF_Motion.11s", rm, header)
      call read_qtf_rows(root//"/Wamit_format/QTF_Indirect.11s", rf)
    end if
    do md = 1, 6
      kk = 0
      do n = 1, size(rm, 2)
        if (rm(1, n) == w(pairs(1, p)) .and. rm(2, n) == w(pairs(3, p)) .and. rm(3, n) == hb(pairs(2, p)) .and. &
            rm(4, n) == hb(pairs(4, p)) .and. nint(rm(5, n)) == md) kk = n
      end do
      call check(error, kk > 0, "a pair has no row in the motion file")
      if (allocated(error)) return
      s = merge(refl**2, refl, md > 3)
      xi_file(md) = conjg(cmplx(rm(8, kk), rm(9, kk), kind=8))/s
      do n = 1, size(rf, 2)
        if (all(rf(1:5, n) == rm(1:5, kk))) exit
      end do
      call check(error, n <= size(rf, 2), "a motion row has no force row")
      if (allocated(error)) return
      f_file(md) = conjg(cmplx(rf(8, n), rf(9, n), kind=8))*RHO_SOLVER*G_SOLVER*s
    end do
    berr_file = maxval(abs(matmul(lhss(:, :, p), xi_file) - f_file))/(maxval(sum(abs(lhss(:, :, p)), dim=2))* &
                maxval(abs(xi_file)))
    worst_file = max(worst_file, berr_file)
  end do
  write (line, '(a,es10.2,a,es10.2)') "motion_equals_h_times_f: worst backward error unscaled ", worst, &
    ", of the file rows (refl = 2) ", worst_file
  print '(4x,a)', trim(line)
  call check(error, worst <= 1d-13 .and. worst_file <= 1d-5, line)
  call remove_files()
end subroutine test_motion_equals_h_times_f

!> Phase 6c, R4: write_motion_row scales a translation by refl (xi over A^2/L) and a rotation by refl^2 (over A^2/L^2), and
!> conjugates as write_qtf_row; checked against a row scaled by hand, with refl = 2
subroutine test_motion_scaling(error)
  type(error_type), allocatable, intent(out) :: error
  complex(8), parameter :: value = (0.3d0, -0.4d0)
  integer :: unit, idx, md
  real(8) :: w1, w2, b1, b2, mo, ph, re, im, s
  complex(8) :: hand

  open (newunit=unit, status="scratch")
  call write_motion_row(unit, 0.8d0, 1.2d0, 0d0, 30d0, 2, value, 2d0, .false.)
  call write_motion_row(unit, 0.8d0, 1.2d0, 0d0, 30d0, 5, value, 2d0, .true.)
  rewind (unit)
  do md = 1, 2
    read (unit, *) w1, w2, b1, b2, idx, mo, ph, re, im
    s = merge(4d0, 2d0, md == 2)
    hand = conjg(value)*s
    call check(error, w1 == 0.8d0 .and. w2 == 1.2d0 .and. b1 == 0d0 .and. b2 == 30d0 .and. idx == merge(5, 2, md == 2), &
               "keys of the motion row")
    if (allocated(error)) exit
    call check(error, abs(re - real(hand, 8)) <= 1d-6*abs(hand) .and. abs(im - aimag(hand)) <= 1d-6*abs(hand) .and. &
               abs(mo - abs(hand)) <= 1d-6*abs(hand) .and. &
               abs(ph - atan2(aimag(hand), real(hand, 8))*180d0/3.141592653589793d0) <= 1d-4, &
               "the motion row is not conj(xi) refl^k")
    if (allocated(error)) exit
  end do
  close (unit)
end subroutine test_motion_scaling

!> Phase 6c: QTF_Motion.11s of the moored cylinder (qtf = 2) has row (i, b1, j, b2) equal to row (j, b2, i, b1): the sum
!> force is symmetric (T+), and so is the matrix, which depends on Omega+ only. To the 7 digits of the file.
subroutine test_motion_sum_symmetric(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8), allocatable :: rm(:, :)
  character(len=:), allocatable :: header
  real(8) :: scale, gap
  integer :: i, j
  logical :: found

  call sum_pipeline(error, 0d0, mesh, geom, kind_qtf=2, indirect_only=.true., motion=2)
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)
  call read_motion_rows(root//"/Wamit_format/QTF_Motion.11s", rm, header)
  call check(error, size(rm, 2), 96, "rows of QTF_Motion.11s")
  if (allocated(error)) return
  scale = maxval(abs(rm(8:9, :)))
  gap = 0d0
  do i = 1, size(rm, 2)
    found = .false.
    do j = 1, size(rm, 2)
      if (rm(1, j) /= rm(2, i) .or. rm(2, j) /= rm(1, i) .or. rm(3, j) /= rm(4, i) .or. rm(4, j) /= rm(3, i)) cycle
      if (rm(5, j) /= rm(5, i)) cycle
      found = .true.
      gap = max(gap, maxval(abs(rm(8:9, j) - rm(8:9, i))))
      exit
    end do
    call check(error, found, "a row of QTF_Motion.11s has no mirror")
    if (allocated(error)) return
  end do
  print '(4x,a,es10.2,a,es10.2)', "motion_sum_symmetric: max |row - mirror| ", gap, ", max |QTF_Motion.11s| ", scale
  call check(error, scale > 0d0 .and. gap <= 1d-6*scale, "QTF_Motion.11s is not symmetric")
  call remove_files()
end subroutine test_motion_sum_symmetric

!> Phase 6c with IRSP = 0: the sum pairs at or above the estimate of the lowest irregular frequency (2.2 rad/s for the
!> cylinder of radius 5) are cut from QTF_Indirect.11s (sum_pair_index), and from QTF_Motion.11s too. Grid 0.7, 1.3 (no
!> sum on the main grid, whose psi the pipeline does not keep): the pair (1.3, 1.3) is cut, 72 of the 96 rows are kept;
!> the motion file has the same keys in the same order.
subroutine test_motion_sum_cutoff_rows(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8), allocatable :: rs(:, :), rm(:, :)
  character(len=:), allocatable :: header
  real(8) :: cut

  call sum_pipeline(error, 0d0, mesh, geom, kind_qtf=2, grid=[0.7d0, 1.3d0], indirect_only=.true., motion=2, irsp0=.true.)
  if (allocated(error)) return
  cut = irregular_frequency_estimate(geom, G_SOLVER)
  call qtf_finish(mesh, geom, 1d0)
  call read_qtf_rows(root//"/Wamit_format/QTF_Indirect.11s", rs)
  call read_motion_rows(root//"/Wamit_format/QTF_Motion.11s", rm, header)
  print '(4x,a,f8.4,a,i0,a,i0)', "motion_sum_cutoff_rows: cut at ", cut, " rad/s, rows of .11s ", size(rs, 2), &
    ", of QTF_Motion.11s ", size(rm, 2)
  call check(error, size(rs, 2), 72, "rows of .11s with the cut")
  if (allocated(error)) return
  call check(error, size(rm, 2), size(rs, 2), "rows of QTF_Motion.11s")
  if (allocated(error)) return
  call check(error, all(rm(1:5, :) == rs(1:5, :)), "QTF_Motion.11s does not have the keys of .11s")
  if (allocated(error)) return
  call check(error, all(rm(1, :) + rm(2, :) < cut), "QTF_Motion.11s has a row above the cut")
  call remove_files()
end subroutine test_motion_sum_cutoff_rows

!> Phase 6c: qtf_stored_mb with nmotion grows by xi2 and xi2_sum, 16 bytes x 6 nbody x nfreq^2 per slice
subroutine test_stored_mb_motion(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: nfreq = 3, nb = 2, nbody = 2, npsi = 7, nmotion = 8
  real(8) :: mb0, mb1, expected

  mb0 = qtf_stored_mb(10, 4, nbody, nfreq, nb, npsi=npsi, npairs=4, nsum_pairs=4)
  mb1 = qtf_stored_mb(10, 4, nbody, nfreq, nb, npsi=npsi, npairs=4, nsum_pairs=4, nmotion=nmotion)
  expected = 16d0 * 6d0 * nbody * nfreq**2 * nmotion / 1d6
  call check(error, abs((mb1 - mb0) - expected) <= 1d-12 * expected, "the motion stores in qtf_stored_mb")
end subroutine test_stored_mb_motion

!> Phase 6c: qtf_motion = 0 writes no QTF_Motion file, and the force files of qtf_motion = 1 are those of qtf_motion = 0
!> to the byte (qtf = 3, the moored cylinder, qtf_direct = 0)
subroutine test_motion_off_unchanged(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  character(len=*), parameter :: names(4) = [character(len=20) :: "QTF_Quadratic.10d", "QTF_Indirect.11d", &
                                             "QTF_Quadratic.10s", "QTF_Indirect.11s"]
  type :: text_t
    character(len=:), allocatable :: s
  end type text_t
  type(text_t) :: off(4)
  logical :: has_d, has_s
  integer :: f, km

  do km = 0, 2, 2
    call sum_pipeline(error, 0d0, mesh, geom, kind_qtf=3, indirect_only=.true., motion=km)
    if (allocated(error)) return
    call qtf_finish(mesh, geom, 1d0)
    inquire (file=root//"/Wamit_format/QTF_Motion.11d", exist=has_d)
    inquire (file=root//"/Wamit_format/QTF_Motion.11s", exist=has_s)
    if (km > 0) then
      call check(error, has_d .and. has_s, "no QTF_Motion.11d or .11s with qtf_motion = 1")
    else
      call check(error, .not. (has_d .or. has_s), "a QTF_Motion file with qtf_motion = 0")
    end if
    if (allocated(error)) return
    do f = 1, 4
      if (km == 0) then
        off(f)%s = slurp(root//"/Wamit_format/"//trim(names(f)))
        call check(error, len(off(f)%s) > 0, trim(names(f))//" is empty")
      else
        call check(error, slurp(root//"/Wamit_format/"//trim(names(f))) == off(f)%s, &
                   trim(names(f))//" differs between qtf_motion = 0 and 1")
      end if
      if (allocated(error)) return
    end do
    call remove_files()
  end do
end subroutine test_motion_off_unchanged

!> Phase 6c, plan review focus 2: qtf_motion = 1 with motion_state = fixed computes nothing and writes no QTF_Motion file;
!> qtf_init takes no body matrices then, and the force files are written as usual
subroutine test_motion_fixed_no_file(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  logical :: has_d, has_s, has_force

  call sum_pipeline(error, 0d0, mesh, geom, kind_qtf=3, indirect_only=.true., motion=2, fixed=.true.)
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)
  inquire (file=root//"/Wamit_format/QTF_Motion.11d", exist=has_d)
  inquire (file=root//"/Wamit_format/QTF_Motion.11s", exist=has_s)
  inquire (file=root//"/Wamit_format/QTF_Indirect.11d", exist=has_force)
  call remove_files()
  call check(error, has_force .and. .not. (has_d .or. has_s), "a fixed body wrote a QTF_Motion file, or no .11d")
end subroutine test_motion_fixed_no_file

!> The whole file as one string (stream access); empty if it does not exist
function slurp(path) result(text)
  character(len=*), intent(in) :: path
  character(len=:), allocatable :: text
  integer :: unit, n
  logical :: exists

  inquire (file=path, exist=exists)
  if (.not. exists) then
    text = ""
    return
  end if
  open (newunit=unit, file=path, access="stream", form="unformatted", status="old", action="read")
  inquire (unit=unit, size=n)
  allocate (character(len=n) :: text)
  read (unit) text
  close (unit)
end function slurp

!> The rows of a QTF file, r(9, n): PER1 PER2 BETA1 BETA2 I Mod Pha Re Im
subroutine read_qtf_rows(path, r)
  character(len=*), intent(in) :: path
  real(8), allocatable, intent(out) :: r(:, :)
  integer :: n, unit, i, idx
  real(8) :: v(8)

  n = count_rows(path)
  allocate (r(9, max(n, 0)))
  if (n <= 0) return
  open (newunit=unit, file=path, status="old", action="read")
  do i = 1, n
    read (unit, *) v(1:4), idx, v(5:8)
    r(1:4, i) = v(1:4)
    r(5, i) = idx
    r(6:9, i) = v(5:8)
  end do
  close (unit)
end subroutine read_qtf_rows

!> The data rows of a terms file (lines without #), r(ncol, n), and the number of columns of the first row
subroutine read_terms_rows(path, r, ncol)
  character(len=*), intent(in) :: path
  real(8), allocatable, intent(out) :: r(:, :)
  integer, intent(out) :: ncol
  character(len=1000) :: line
  integer :: unit, ios, n, i
  logical :: blank

  ncol = 0
  n = 0
  open (newunit=unit, file=path, status="old", action="read")
  do
    read (unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0 .or. line(1:1) == "#") cycle
    n = n + 1
    if (n == 1) then
      blank = .true.
      do i = 1, len_trim(line)
        if (line(i:i) /= " " .and. blank) ncol = ncol + 1
        blank = line(i:i) == " "
      end do
    end if
  end do
  allocate (r(ncol, n))
  rewind (unit)
  i = 0
  do
    read (unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) == 0 .or. line(1:1) == "#") cycle
    i = i + 1
    read (line, *) r(:, i)
  end do
  close (unit)
end subroutine read_terms_rows

!> Rows of both files: the full square of frequency pairs, heading pairs for mean_drift = 2, and qtf_max_dw on .10d only.
subroutine test_files(error)
  type(error_type), allocatable, intent(out) :: error
  type(so_options) :: opts
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: states(size(omegas), size(betas))
  complex(8), allocatable :: wls(:,:,:)

  opts%mean_drift = 2; opts%qtf = 3; opts%qtf_max_dw = 0.25d0
  call fill_store(opts, mesh, geom, states, wls, error)
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)
  ! .10s: 3 x 3 frequency pairs, 2 x 2 heading pairs, 6 rows each. .10d: |dw| <= 0.25 keeps the diagonal and (0.8, 1.0) twice
  call check(error, count_rows(root//"/Wamit_format/QTF_Quadratic.10s"), 9 * 4 * 6, "rows of .10s")
  if (allocated(error)) return
  call check(error, count_rows(root//"/Wamit_format/QTF_Quadratic.10d"), 5 * 4 * 6, "rows of .10d")
  if (allocated(error)) return
  call check_file_symmetry(error, root//"/Wamit_format/QTF_Quadratic.10d", .true.)
  if (allocated(error)) return
  call check_file_symmetry(error, root//"/Wamit_format/QTF_Quadratic.10s", .false.)
  if (allocated(error)) return

  ! mean_drift = 1: beta1 = beta2 only, and qtf = 1 writes no .10s
  call remove_files()
  opts%mean_drift = 1; opts%qtf = 1; opts%qtf_max_dw = 0d0
  call fill_store(opts, mesh, geom, states, wls, error)
  if (allocated(error)) return
  call qtf_finish(mesh, geom, 1d0)
  call check(error, count_rows(root//"/Wamit_format/QTF_Quadratic.10d"), 9 * 2 * 6, "rows of .10d, one heading pair each")
  if (allocated(error)) return
  call check(error, count_rows(root//"/Wamit_format/QTF_Quadratic.10s"), -1, ".10s written for qtf = 1")
  call remove_files()
  call remove_root()
end subroutine test_files

!> qtf_files runs last: it removes the test folder.
subroutine remove_root()
  character(len=32) :: os_name
  integer :: estat

  call get_environment_variable("OS", os_name)
  if (trim(os_name) == "Windows_NT") then
    call execute_command_line("if exist "//root//" rmdir /s /q "//root, exitstat=estat)
  else
    call execute_command_line("rm -rf "//root, exitstat=estat)
  end if
end subroutine remove_root

!> Number of rows of a file, or -1 if it does not exist.
integer function count_rows(path)
  character(len=*), intent(in) :: path
  character(len=256) :: line
  logical :: exists
  integer :: unit, ios

  inquire(file=path, exist=exists)
  count_rows = -1
  if (.not. exists) return
  count_rows = 0
  open(newunit=unit, file=path, status="old", action="read")
  do
    read(unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    if (len_trim(line) > 0) count_rows = count_rows + 1
  end do
  close(unit)
end function count_rows

!> Every row (w1, w2, b1, b2, I) has its partner (w2, w1, b2, b1, I): the conjugate for .10d, the same value for .10s.
subroutine check_file_symmetry(error, path, hermitian)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), intent(in) :: path
  logical, intent(in) :: hermitian
  real(8), allocatable :: r(:,:)
  integer :: n, unit, i, j, idx
  real(8) :: v(8), scale, diff
  logical :: found

  n = count_rows(path)
  allocate(r(9, n))
  open(newunit=unit, file=path, status="old", action="read")
  do i = 1, n
    read(unit, *) v(1:4), idx, v(5:8)
    r(1:4, i) = v(1:4)
    r(5, i) = idx
    r(6:9, i) = v(5:8)
  end do
  close(unit)
  scale = maxval(abs(r(8:9, :)))
  do i = 1, n
    found = .false.
    do j = 1, n
      if (r(1, j) /= r(2, i) .or. r(2, j) /= r(1, i) .or. r(3, j) /= r(4, i) .or. r(4, j) /= r(3, i)) cycle
      if (r(5, j) /= r(5, i)) cycle
      found = .true.
      if (hermitian) then
        diff = abs(r(8, j) - r(8, i)) + abs(r(9, j) + r(9, i))
      else
        diff = abs(r(8, j) - r(8, i)) + abs(r(9, j) - r(9, i))
      end if
      call check(error, diff <= 1d-5 * scale, "row and its partner do not have the symmetry of "//path)
      if (allocated(error)) return
      exit
    end do
    call check(error, found, "row has no partner in "//path)
    if (allocated(error)) return
  end do
end subroutine check_file_symmetry

subroutine remove_files()
  character(len=*), parameter :: names(13) = [character(len=25) :: "QTF_Quadratic.10d", "QTF_Quadratic.10s", &
                                                 "QTF_Direct.12s", "QTF_DirectTermsSum.txt", "QTF_Indirect.11s", &
                                                 "QTF_PotentialTermsSum.txt", "QTF_Indirect.11d", "QTF_Direct.12d", &
                                                 "QTF_Motion.11d", "QTF_Motion.11s", "QTF_Motion.12d", "QTF_Motion.12s", &
                                                 "QTF_PotentialTerms.txt"]
  integer :: unit, ios, i
  logical :: exists

  do i = 1, size(names)
    inquire(file=root//"/Wamit_format/"//trim(names(i)), exist=exists)
    if (.not. exists) cycle
    open(newunit=unit, file=root//"/Wamit_format/"//trim(names(i)), status="old", iostat=ios)
    if (ios == 0) close(unit, status="delete")
  end do
end subroutine remove_files

end module test_qtf
