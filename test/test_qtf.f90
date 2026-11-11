!> Unit tests of QTFDriver: the stored pairs against the kernel, their symmetries, and the rows of QTF_Quadratic.10d/.10s.
!> The store is module state, so unit_tester runs this suite serially.
module test_qtf
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use Const_mod, only : G_SOLVER => G, RHO_SOLVER => RHO
  use IO, only : CreateDirectory
  use SecondOrderTypes, only : so_options, panel_mesh, first_order_state, surface_geom, QTF_TERMS_BODY
  use SurfaceGeometry, only : build_surface_geom, DEFAULT_SHARP_DEG
  use NearFieldDrift, only : near_field_pair, KERNEL_DIFF, KERNEL_SUM
  use QTFDriver, only : qtf_init, qtf_keep, qtf_pair_load, qtf_finish, qtf_kept, qtf_radiation_list, qtf_radiation_omegas, &
                        qtf_keep_psi, qtf_nonfinite_psi, qtf_nonfinite_pairs, qtf_psi_tolerance, qtf_psi_index
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use SecondOrderOutput, only : write_qtf_row
  use second_order_fixtures, only : make_cylinder, maccamy_fuchs_state, maccamy_fuchs_phi
  implicit none
  private

  public :: collect_qtf

  real(8), parameter :: radius = 5d0, depth = 10d0, g = 9.81d0, rho = 1025d0
  real(8), parameter :: omegas(3) = [0.8d0, 1.0d0, 1.3d0], betas(2) = [0d0, 30d0]
  complex(8), parameter :: ci = (0d0, 1d0)
  character(len=*), parameter :: root = "hams_test_qtf_files"

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
    new_unittest("potential_part_period_grid", test_potential_part_period_grid), &
    new_unittest("psi_store_order", test_psi_store_order), &
    new_unittest("psi_nonfinite_counted", test_psi_nonfinite_counted), &
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
subroutine fill_store(opts, mesh, geom, states, wls, error, body)
  type(so_options), intent(in) :: opts
  type(panel_mesh), intent(out) :: mesh
  type(surface_geom), intent(out) :: geom
  type(first_order_state), intent(out) :: states(size(omegas), size(betas))
  complex(8), allocatable, intent(out) :: wls(:,:,:)
  type(error_type), allocatable, intent(out) :: error
  logical, intent(in), optional :: body   ! qtf_terms = body: the store also takes psi at the differences of omegas
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
  if (present(body)) then
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
  character(len=*), parameter :: names(2) = [character(len=17) :: "QTF_Quadratic.10d", "QTF_Quadratic.10s"]
  integer :: unit, ios, i
  logical :: exists

  do i = 1, 2
    inquire(file=root//"/Wamit_format/"//names(i), exist=exists)
    if (.not. exists) cycle
    open(newunit=unit, file=root//"/Wamit_format/"//names(i), status="old", iostat=ios)
    if (ios == 0) close(unit, status="delete")
  end do
end subroutine remove_files

end module test_qtf
