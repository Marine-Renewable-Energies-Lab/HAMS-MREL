!> Unit tests of the right side of the direct second-order problem (Phase 7a): ASSB_RHS_FIELD against ASSB_RBC, and the
!> B2 part (QTFDirect) in the Stokes form with G against the direct form int G n.[(X.grad)u] dS on analytic fields, and
!> the closed-form panel moments of 1/r (rankine_moments) against adaptive quadrature.
!> These tests set module globals (the mesh arrays, the stored Green function, the wave), so the suite runs serially
!> (see unit_tester.f90).
module test_qtf_direct
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use HAMS_mod, only : NTHREAD, ISOLV, IRSP
  use WaveDyn_mod, only : VWAVE => V
  use PanelMesh_mod, only : ISX, ISY, NSYS, NELEM, DS, PNSZ, XYZ_P, DXYZ_P, DS_MULTI_COMB, PNSZ_MULTI_COMB, &
                            XYZ_GLOBAL_MULTI_COMB_P, DXYZ_MULTI_COMB, IPIV, IPIV_MULTI_COMB, NELEM_TOTAL, &
                            iNELEM_TOTAL, iXYZ_GLOBAL_MULTI_COMB_P, iDXYZ_MULTI_COMB
  use Inerfs_mod, only : TNELEM, iNELEM, iXYZ_P, iDXYZ_P
  use Potentials_mod, only : CGRN, RKBN, CGRN_MULTI_COMB, RKBN_MULTI_COMB, MXPOT, MXPOT_MULTI_COMB, DGRN, PKBN, &
                             DGRN_MULTI_COMB, PKBN_MULTI_COMB
  use LinearMatrix_mod, only : AMAT, BRMAT, AMAT_MULTI, BRMAT_MULTI, CMAT, DRMAT, CMAT_MULTI, DRMAT_MULTI
  use Const_mod, only : RXY
  use AssbMatx, only : ASSB_RBC, ASSB_RHS_FIELD, ASSB_RHS_FIELD_IRR, ASSB_LEFT, RADIATION_SOLVER
  use AssbMatxMulti, only : ASSB_RBC_MULTI, ASSB_RHS_FIELD_MULTI, ASSB_RHS_FIELD_MULTI_IRR, ASSB_LEFT_MULTI, &
                            RADIATION_SOLVER_MULTI
  use AssbMatx_irr, only : ASSB_LEFT_IRR, ASSB_RBC_IRR, RADIATION_SOLVER_IRR
  use AssbMatx_irrMulti, only : ASSB_LEFT_IRR_MULTI, ASSB_RBC_IRR_MULTI
  use SecondOrderTypes, only : so_options, QTF_TERMS_BODY
  use SecondOrderInput, only : validate_second_order
  use ImplementSubs, only : SetWaveFromOmega
  use SecondOrderTypes, only : panel_mesh, surface_geom, first_order_state
  use FirstOrderState, only : image_signs, radiation_psi, wave_gravity
  use QTFPotential, only : potential_pair, psi_surface_gradient, psi_waterline, NPOT_TERMS
  use SurfaceVelocity, only : tangential_fit
  use NearFieldDrift, only : pair_product, rigid_displacement, KERNEL_DIFF
  use QTFDirect, only : direct_rhs_b2, b2_stokes_densities, b2_stokes_at_point, green_panel_integrals, rankine_moments, &
                        sharp_edge_panels, direct_rhs_b2_block, direct_setup, direct_edge_mask, direct_factor, &
                        direct_solve, direct_at_omega, direct_lid_points
  use SharpEdge, only : DEFAULT_EDGE_DEG
  use CalGreenFunc, only : CALGREEN, CALGREEN_IRR
  use CalGreenFuncMulti, only : CALGREEN_MULTI, CALGREEN_MULTI_IRR
  use second_order_fixtures, only : make_cylinder, set_panel
  use test_waterline_field, only : set_wave, load_panels, load_multi
  use test_qtf_potential, only : make_sphere, motions, velocity, normal_hessian, make_state, geometry
  implicit none
  private

  public :: collect_qtf_direct

  ! The first-order states of the pair, as test_qtf_potential: frequency and heading
  real(8), parameter :: w1 = 0.8d0, b1 = 0d0, w2 = 0.7d0, b2 = 0.3d0
  ! ASSB_RHS_FIELD sums the same terms as ASSB_RBC in another order
  real(8), parameter :: rbc_limit = 1d-12
  ! Stokes form with G against the direct form at collocation points, depth 50, relative: below b2_limit at 2048 panels and
  ! an observed order of at least b2_order from 128 to 2048 panels (controller ruling, Phase 7a Task 2)
  ! The level is limited by the fitted density gradient on the pole fan triangles (1.24e-3 measured); with the analytic
  ! gradient the error is 2.8e-4. See the Phase 7a ledger, Task 2.
  real(8), parameter :: b2_limit = 1.5d-3, b2_order = 1.5d0
  ! The same 10 m off the surface with a near-Rankine Green function: O(h^2) quadrature error
  real(8), parameter :: off_limit = 5d-4, off_order = 1.7d0
  ! Closed-form moments against adaptive quadrature
  real(8), parameter :: moment_limit = 1d-10
  ! Half mesh against whole: the same sums in another order, and the Green function at mirrored points
  real(8), parameter :: half_limit = 1d-8
  ! Stored integrals (CALGREEN, single precision, GRN_KIND = 4) against green_panel_integrals (double), relative to the
  ! largest entry of the right side. The flag rules differ only for mirrored pairs beyond 50 panel sizes, which these
  ! small fixtures do not have, so the difference is the single-precision storage (about 6e-8 per value).
  real(8), parameter :: stored_limit = 1d-5
  ! Sphere of radius 5 at centre depth 20: 128, 512 and 2048 panels
  integer, parameter :: sph_nlat(3) = [8, 16, 32]
  ! Phase 7a Task 3. B2 by the Stokes form against ASSB_RHS_FIELD of the analytic panel values (512-panel sphere): the
  ! mesh error, 7.9e-3 measured; a sign error gives about 2
  real(8), parameter :: b2_sign_limit = 2d-2
  ! b = n_k against psi_k: the same matrix and right side, rounding only
  real(8), parameter :: unit_normal_limit = 1d-10
  ! The guard effect on the production path against the same B2 parts through direct_solve: rounding only
  real(8), parameter :: guard_limit = 1d-8
  ! T_dir against T_pot of the indirect method (512-panel sphere): the mesh error, 9.5e-4 measured
  real(8), parameter :: indirect_limit = 3d-3
  ! Phase 7a Task 3b, irregular frequency removal. The cylinder of load_multi (test_waterline_field), for its lids
  real(8), parameter :: multi_radius = 5d0, multi_spacing = 30d0
  integer, parameter :: multi_ntheta = 12, multi_nr = 2
  ! direct_vs_indirect_irr: the rounded cylinder of radius 5 and draft 5 at depth 50 (448 panels, 192 lid points). Its
  ! discrete first irregular frequency, irr_omegas(3), is where the lid changes T_dir most: a scan from 1.90 to 2.70 rad/s
  ! by 0.02 and from 2.188 to 2.212 by 0.001 gave a peak of 15.2% at 2.195 (half width about 0.002 rad/s) on a background
  ! of 1-1.5%. irr_omegas(1) is a low Omega, irr_omegas(2) a regular one next to the irregular frequency. State 2 is at
  ! irr_low rad/s, state 1 at irr_low + Omega.
  real(8), parameter :: irr_depth = 50d0, irr_low = 0.6d0, irr_omegas(3) = [1.0d0, 2.25d0, 2.195d0]
  ! With the lid, T_dir against T_pot at the low Omega: the mesh error, 2.07e-3 measured (2.35e-3 without the lid). At
  ! 2.25 and 2.195 rad/s the mesh error is 3.1% and 3.0% (state 1 at 2.8 rad/s, 8 panels per wavelength on the side):
  ! with the lid the irregular frequency adds nothing to it (ratio 0.98 measured, limit irr_ratio_limit), while without
  ! the lid T_dir against T_pot is 14.7% there. Without the lid against with, T_dir: 1.75e-3 at the low Omega (limit
  ! irr_regular_limit), 15.1% at the irregular one (at least irr_spike_min).
  real(8), parameter :: irr_indirect_limit = 5d-3, irr_mesh_limit = 5d-2, irr_ratio_limit = 1.2d0
  real(8), parameter :: irr_regular_limit = 5d-3
  ! Phase 7a Task 5, the spike over a window: irr_window spans the peak by steps of 0.0025 rad/s, so that one Omega is
  ! within 0.00125 rad/s of the peak wherever it lies in the window; irr_ends are 0.05 rad/s either side. Measured: with
  ! the lid T_dir is 1.06e-2 to 1.20e-2 from the straight line through the ends (the curvature over 0.1 rad/s; 1.20e-2
  ! at 2.195), limit irr_smooth_limit (1.3x measured); without the lid it is 15.1% from the lid result at 2.195, 6.2% at
  ! 0.0025 rad/s from it and 1.6-1.9% at the window ends, and irr_spike_min is below the value 0.0025 rad/s off the peak
  real(8), parameter :: irr_window(9) = [2.185d0, 2.1875d0, 2.19d0, 2.1925d0, 2.195d0, 2.1975d0, 2.2d0, 2.2025d0, 2.205d0]
  real(8), parameter :: irr_ends(2) = [2.145d0, 2.245d0]
  real(8), parameter :: irr_smooth_limit = 1.6d-2, irr_spike_min = 5d-2
  ! rhs_b2_lid_rows: the lid rows of the half mesh against the whole mesh, both fresh: 1.05e-8 measured (the hull rows of
  ! rhs_b2_half_matches_whole have half_limit), 2x measured, rounded up
  real(8), parameter :: lid_half_limit = 2.1d-8

contains

  subroutine collect_qtf_direct(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("rhs_field_matches_rbc", test_rhs_field_matches_rbc), &
      new_unittest("rhs_field_matches_rbc_multi", test_rhs_field_matches_rbc_multi), &
      new_unittest("rankine_moments_quadrature", test_rankine_moments), &
      new_unittest("rhs_b2_stokes_vs_direct_smooth", test_b2_smooth), &
      new_unittest("rhs_b2_stokes_vs_direct_off_surface", test_b2_off_surface), &
      new_unittest("rhs_b2_half_matches_whole", test_b2_half), &
      new_unittest("rhs_b2_stored_matches_direct", test_b2_stored), &
      new_unittest("rhs_b2_sharp_edge_fallback", test_b2_sharp_edge), &
      new_unittest("rhs_b2_sign_through_rhs_field", test_b2_sign), &
      new_unittest("direct_zero_forcing", test_direct_zero_forcing), &
      new_unittest("direct_unit_normal_is_radiation", test_direct_unit_normal), &
      new_unittest("direct_edge_guard_production", test_direct_edge_guard), &
      new_unittest("direct_matches_indirect_pair", test_direct_vs_indirect), &
      new_unittest("direct_irr_accepted", test_direct_irr_accepted), &
      new_unittest("direct_irr_unit_normal", test_direct_irr_unit_normal), &
      new_unittest("rhs_b2_lid_rows", test_b2_lid_rows), &
      new_unittest("direct_vs_indirect_irr", test_direct_vs_indirect_irr) &
      ]
  end subroutine collect_qtf_direct

  ! ------------------------------------------------------------------------------------------------ ASSB_RHS_FIELD

  !> Factors of the generalized normals n_1..6 on symmetry image is: (fx, fy, 1, fy, fx, fx fy), as expand_images.
  function mode_signs(is) result(fj)
    integer, intent(in) :: is
    real(8) :: fj(6), f(3)
    f = image_signs(is, ISX, ISY)
    fj = [f(1), f(2), 1d0, f(2), f(1), f(1)*f(2)]
  end function mode_signs

  !> Synthetic stored Green function and panels; ne panels, near and far ones (the 50 panel-size switch).
  subroutine synthetic_single(ne)
    integer, intent(in) :: ne
    integer :: i, j, is, k
    if (allocated(CGRN)) deallocate(CGRN)
    if (allocated(RKBN)) deallocate(RKBN)
    if (allocated(XYZ_P)) deallocate(XYZ_P)
    if (allocated(DXYZ_P)) deallocate(DXYZ_P)
    if (allocated(DS)) deallocate(DS)
    if (allocated(PNSZ)) deallocate(PNSZ)
    allocate(CGRN(ne, ne, NSYS, 4), RKBN(ne, ne, NSYS, 4), XYZ_P(ne, 3), DXYZ_P(ne, 6), DS(ne), PNSZ(ne))
    do i = 1, ne
      XYZ_P(i, :) = [10d0*sin(1.7d0*i), 3d0*cos(0.9d0*i), -1d0 - real(i, 8)]
      DS(i) = 0.5d0 + 0.1d0*i
      PNSZ(i) = 0.02d0*i      ! panel i is near only within 0.02 i * 50 = i m
      do k = 1, 6
        DXYZ_P(i, k) = cos(0.3d0*i + 1.1d0*k)
      end do
    end do
    do k = 1, 4
      do is = 1, NSYS
        do j = 1, ne
          do i = 1, ne
            RKBN(i, j, is, k) = real(sin(1d0*i + 2d0*j + 3d0*k + 0.7d0*is), kind(RKBN))
            CGRN(i, j, is, k) = cmplx(cos(1d0*i - 1d0*j + k + is), sin(3d0*i + j - k - 2d0*is), kind(CGRN))
          end do
        end do
      end do
    end do
  end subroutine synthetic_single

  !> With V = n_k on every image, ASSB_RHS_FIELD equals ASSB_RBC column k, without symmetry and with each image.
  subroutine test_rhs_field_matches_rbc(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: ne = 9
    integer :: icase, is, k
    integer, parameter :: sx(3) = [0, 0, 1], sy(3) = [0, 1, 0]
    complex(8), allocatable :: brmat(:, :, :), vv(:, :, :), b(:, :, :)
    real(8) :: worst, scale
    character(len=120) :: line

    NTHREAD = 2
    worst = 0d0
    do icase = 1, 3
      ISX = sx(icase)
      ISY = sy(icase)
      NSYS = 1 + ISX + ISY
      NELEM = ne
      call synthetic_single(ne)
      allocate(brmat(ne, 6, NSYS), vv(ne, 6, NSYS), b(ne, 6, NSYS))
      call ASSB_RBC(brmat, ne, NSYS)
      do is = 1, NSYS
        do k = 1, 6
          vv(:, k, is) = mode_signs_k(is, k)*DXYZ_P(:, k)
        end do
      end do
      call ASSB_RHS_FIELD(vv, b, ne, NSYS, 6)
      scale = maxval(abs(brmat))
      worst = max(worst, maxval(abs(b - brmat))/scale)
      deallocate(brmat, vv, b)
    end do
    write(line, '(a,es10.2)') "ASSB_RHS_FIELD against ASSB_RBC, worst relative difference ", worst
    print '(4x,a)', trim(line)
    call check(error, worst <= rbc_limit, line)
    ISX = 0; ISY = 0; NSYS = 1
  end subroutine test_rhs_field_matches_rbc

  real(8) function mode_signs_k(is, k)
    integer, intent(in) :: is, k
    real(8) :: fj(6)
    fj = mode_signs(is)
    mode_signs_k = fj(k)
  end function mode_signs_k

  !> As rhs_field_matches_rbc in the combined multi-body arrays: two bodies, each with its six modes from ASSB_RBC_MULTI
  !> over its own panels, against one call of ASSB_RHS_FIELD_MULTI with n_k on that body and 0 on the other.
  subroutine test_rhs_field_matches_rbc_multi(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: nper = 5, ne = 2*nper
    integer :: icase, is, k, i, j, bd, first
    integer, parameter :: sy(2) = [0, 1]
    complex(8), allocatable :: brmat(:, :, :), vv(:, :, :), b(:, :, :)
    real(8) :: worst
    character(len=120) :: line

    NTHREAD = 2
    worst = 0d0
    do icase = 1, 2
      ISX = 0
      ISY = sy(icase)
      NSYS = 1 + ISY
      TNELEM = ne
      if (allocated(CGRN_MULTI_COMB)) deallocate(CGRN_MULTI_COMB)
      if (allocated(RKBN_MULTI_COMB)) deallocate(RKBN_MULTI_COMB)
      if (allocated(XYZ_GLOBAL_MULTI_COMB_P)) deallocate(XYZ_GLOBAL_MULTI_COMB_P)
      if (allocated(DXYZ_MULTI_COMB)) deallocate(DXYZ_MULTI_COMB)
      if (allocated(DS_MULTI_COMB)) deallocate(DS_MULTI_COMB)
      if (allocated(PNSZ_MULTI_COMB)) deallocate(PNSZ_MULTI_COMB)
      allocate(CGRN_MULTI_COMB(ne, ne, NSYS, 4), RKBN_MULTI_COMB(ne, ne, NSYS, 4), XYZ_GLOBAL_MULTI_COMB_P(ne, 3), &
               DXYZ_MULTI_COMB(ne, 6), DS_MULTI_COMB(ne), PNSZ_MULTI_COMB(ne))
      do i = 1, ne
        XYZ_GLOBAL_MULTI_COMB_P(i, :) = [10d0*sin(1.7d0*i), 3d0*cos(0.9d0*i), -1d0 - real(i, 8)]
        DS_MULTI_COMB(i) = 0.5d0 + 0.1d0*i
        PNSZ_MULTI_COMB(i) = 0.02d0*i
        do k = 1, 6
          DXYZ_MULTI_COMB(i, k) = cos(0.3d0*i + 1.1d0*k)
        end do
      end do
      do k = 1, 4
        do is = 1, NSYS
          do j = 1, ne
            do i = 1, ne
              RKBN_MULTI_COMB(i, j, is, k) = real(sin(1d0*i + 2d0*j + 3d0*k + 0.7d0*is), kind(RKBN_MULTI_COMB))
              CGRN_MULTI_COMB(i, j, is, k) = cmplx(cos(1d0*i - 1d0*j + k + is), sin(3d0*i + j - k - 2d0*is), &
                                                   kind(CGRN_MULTI_COMB))
            end do
          end do
        end do
      end do
      allocate(brmat(ne, 12, NSYS), vv(ne, 12, NSYS), b(ne, 12, NSYS))
      brmat = (0d0, 0d0)
      vv = (0d0, 0d0)
      do bd = 1, 2
        first = (bd - 1)*nper
        call ASSB_RBC_MULTI(brmat(:, 6*(bd - 1) + 1:6*bd, :), ne, NSYS, first + 1, first + nper)
        do is = 1, NSYS
          do k = 1, 6
            vv(first + 1:first + nper, 6*(bd - 1) + k, is) = mode_signs_k(is, k)*DXYZ_MULTI_COMB(first + 1:first + nper, k)
          end do
        end do
      end do
      call ASSB_RHS_FIELD_MULTI(vv, b, ne, NSYS, 12)
      worst = max(worst, maxval(abs(b - brmat))/maxval(abs(brmat)))
      deallocate(brmat, vv, b)
    end do
    write(line, '(a,es10.2)') "ASSB_RHS_FIELD_MULTI against ASSB_RBC_MULTI, worst relative difference ", worst
    print '(4x,a)', trim(line)
    call check(error, worst <= rbc_limit, line)
    ISX = 0; ISY = 0; NSYS = 1
  end subroutine test_rhs_field_matches_rbc_multi

  ! ------------------------------------------------------------------------------------------------ B2

  !> First-order state with total velocity u (in grad_i) and rigid motion xi, for b2_stokes_densities.
  subroutine velocity_state(mesh, omega, beta, xi, s)
    type(panel_mesh), intent(in) :: mesh
    real(8), intent(in) :: omega, beta
    complex(8), intent(in) :: xi(6, 1)
    type(first_order_state), intent(out) :: s
    complex(8), allocatable :: u(:, :)
    call velocity(mesh, omega, beta, u)
    s%omega = omega
    allocate(s%grad_i(3, mesh%n), s%grad_p(3, mesh%n), s%xi(6, 1))
    s%grad_i = u
    s%grad_p = (0d0, 0d0)
    s%xi = xi
  end subroutine velocity_state

  !> Density of the direct form on each panel: wd = n_p.[(X.grad)u] (kernel pair) with the analytic Hessian at the
  !> centroid, and its in-plane gradient wg, from tangential_fit of the vector (X.grad)u dotted with n_p, so that the
  !> reference is integrated as accurately as the Stokes form (linear density on the near panels).
  subroutine direct_density(mesh, geom, x1, x2, wd, wg)
    type(panel_mesh), intent(in) :: mesh
    type(surface_geom), intent(in) :: geom
    complex(8), intent(in) :: x1(6, 1), x2(6, 1)
    complex(8), allocatable, intent(out) :: wd(:), wg(:, :)
    complex(8), allocatable :: hv(:, :)
    complex(8) :: d1(3), d2(3), h1(3), h2(3), g1, g2, dv1(3), dv2(3)
    real(8) :: e(3), t1(3), t2(3)
    integer :: i, j, k
    logical :: ok
    allocate(hv(mesh%n, 3), wd(mesh%n), wg(3, mesh%n))
    do i = 1, mesh%n
      d1 = rigid_displacement(x1(:, 1), mesh%xc(:, i) - mesh%xr(:, 1))
      d2 = rigid_displacement(x2(:, 1), mesh%xc(:, i) - mesh%xr(:, 1))
      do k = 1, 3
        e = 0d0
        e(k) = 1d0
        h1 = normal_hessian(mesh%xc(:, i), e, w1, b1)     ! row k of the Hessian (it is symmetric)
        h2 = normal_hessian(mesh%xc(:, i), e, w2, b2)
        hv(i, k) = (0d0, 0d0)
        do j = 1, 3
          hv(i, k) = hv(i, k) + pair_product(KERNEL_DIFF, d1(j), h1(j), d2(j), h2(j))
        end do
      end do
      wd(i) = sum(mesh%nrm(:, i)*hv(i, :))
    end do
    do i = 1, mesh%n
      wg(:, i) = (0d0, 0d0)
      ok = .true.
      do k = 1, 3
        call tangential_fit(mesh, geom, i, hv(:, k), t1, t2, g1, g2, ok)
        if (.not. ok) exit
        dv1(k) = g1
        dv2(k) = g2
      end do
      if (ok) wg(:, i) = sum(mesh%nrm(:, i)*dv1)*t1 + sum(mesh%nrm(:, i)*dv2)*t2
    end do
  end subroutine direct_density

  !> The direct form at x: sum_p [ int_p G(x, .) dS wd_p + wg_p.int_p (1/r + 1/r1)(xi - xi_p) dS (near panels) ]. On the
  !> panels within 3 panel sizes of x the Rankine part is summed over sub-panels (quads 3 x 3, triangles 4 by their
  !> midpoints, so that a centroid lies inside one), each with the analytic density at its centroid and the fitted
  !> gradient, and the wave part keeps the panel centroid. This makes the reference finer than the Stokes form where
  !> both lose most: the reference with one panel had differences of order 1e-3 of the result between its own variants.
  complex(8) function direct_at_point(mesh, x, wd, wg, x1, x2) result(out)
    type(panel_mesh), intent(in) :: mesh
    real(8), intent(in) :: x(3)
    complex(8), intent(in) :: wd(:), wg(:, :), x1(6, 1), x2(6, 1)
    complex(8) :: g_int, grad_int(3)
    real(8) :: xv(3, 4), sv(3, 4, 9), m0, m0i, g0(3), m1(3, 3), mg(3), mgi(3), sc(3), xim(3)
    integer :: i, nv, ns, k
    out = (0d0, 0d0)
    xim = [x(1), x(2), -x(3)]
    do i = 1, mesh%n
      call green_panel_integrals(mesh, x, i, .false., g_int, grad_int)
      nv = mesh%ncn(i)
      xv(:, 1:nv) = mesh%xnode(:, mesh%ncon(1:nv, i))
      if (norm2(x - mesh%xc(:, i)) > 3d0*mesh%psize(i)) then
        out = out + g_int*wd(i)
        if (norm2(x - mesh%xc(:, i)) > 50d0*mesh%psize(i)) cycle
        call rankine_moments(xv, nv, mesh%nrm(:, i), x, mesh%xc(:, i), m0, g0, m1, mg)
        call rankine_moments(xv, nv, mesh%nrm(:, i), xim, mesh%xc(:, i), m0, g0, m1, mgi)
        out = out + sum(wg(:, i)*(mg + mgi))
        cycle
      end if
      ! wave part at the panel centroid, the Rankine part over the sub-panels
      call rankine_moments(xv, nv, mesh%nrm(:, i), x, mesh%xc(:, i), m0, g0, m1, mg)
      call rankine_moments(xv, nv, mesh%nrm(:, i), xim, mesh%xc(:, i), m0i, g0, m1, mgi)
      out = out + (g_int - m0 - m0i)*wd(i)
      call subpanels(xv, nv, sv, ns)
      do k = 1, ns
        sc = area_centroid3(sv(:, :, k), nv)
        call rankine_moments(sv(:, :, k), nv, mesh%nrm(:, i), x, sc, m0, g0, m1, mg)
        call rankine_moments(sv(:, :, k), nv, mesh%nrm(:, i), xim, sc, m0i, g0, m1, mgi)
        out = out + (m0 + m0i)*point_density(sc, mesh%nrm(:, i)) + sum(wg(:, i)*(mg + mgi))
      end do
    end do

  contains

    complex(8) function point_density(y, nrm) result(w)
      real(8), intent(in) :: y(3), nrm(3)
      complex(8) :: d1(3), d2(3), h1(3), h2(3)
      integer :: j
      d1 = rigid_displacement(x1(:, 1), y - mesh%xr(:, 1))
      d2 = rigid_displacement(x2(:, 1), y - mesh%xr(:, 1))
      h1 = normal_hessian(y, nrm, w1, b1)
      h2 = normal_hessian(y, nrm, w2, b2)
      w = (0d0, 0d0)
      do j = 1, 3
        w = w + pair_product(KERNEL_DIFF, d1(j), h1(j), d2(j), h2(j))
      end do
    end function point_density

  end function direct_at_point

  !> Sub-panels of a quad (3 x 3, bilinear in the vertices) or a triangle (4, by the edge midpoints).
  subroutine subpanels(xv, nv, sv, ns)
    real(8), intent(in) :: xv(3, 4)
    integer, intent(in) :: nv
    real(8), intent(out) :: sv(3, 4, 9)
    integer, intent(out) :: ns
    real(8) :: m12(3), m23(3), m31(3)
    integer :: iu, iv
    sv = 0d0
    if (nv == 4) then
      ns = 0
      do iv = 0, 2
        do iu = 0, 2
          ns = ns + 1
          sv(:, 1, ns) = bil(iu/3d0, iv/3d0)
          sv(:, 2, ns) = bil((iu + 1)/3d0, iv/3d0)
          sv(:, 3, ns) = bil((iu + 1)/3d0, (iv + 1)/3d0)
          sv(:, 4, ns) = bil(iu/3d0, (iv + 1)/3d0)
        end do
      end do
    else
      m12 = 0.5d0*(xv(:, 1) + xv(:, 2))
      m23 = 0.5d0*(xv(:, 2) + xv(:, 3))
      m31 = 0.5d0*(xv(:, 3) + xv(:, 1))
      ns = 4
      sv(:, 1:3, 1) = reshape([xv(:, 1), m12, m31], [3, 3])
      sv(:, 1:3, 2) = reshape([m12, xv(:, 2), m23], [3, 3])
      sv(:, 1:3, 3) = reshape([m31, m23, xv(:, 3)], [3, 3])
      sv(:, 1:3, 4) = reshape([m12, m23, m31], [3, 3])
    end if
  contains
    function bil(u, v) result(y)
      real(8), intent(in) :: u, v
      real(8) :: y(3)
      y = (1 - u)*(1 - v)*xv(:, 1) + u*(1 - v)*xv(:, 2) + u*v*xv(:, 3) + (1 - u)*v*xv(:, 4)
    end function bil
  end subroutine subpanels

  !> Area centroid of a flat polygon of 3 or 4 nodes.
  function area_centroid3(v, nc) result(xa)
    real(8), intent(in) :: v(3, 4)
    integer, intent(in) :: nc
    real(8) :: xa(3), cr(3), e1(3), e2(3), at, asum
    integer :: m
    xa = 0d0
    asum = 0d0
    do m = 2, nc - 1
      e1 = v(:, m) - v(:, 1)
      e2 = v(:, m + 1) - v(:, 1)
      cr = [e1(2)*e2(3) - e1(3)*e2(2), e1(3)*e2(1) - e1(1)*e2(3), e1(1)*e2(2) - e1(2)*e2(1)]
      at = 0.5d0*sqrt(sum(cr**2))
      xa = xa + at*(v(:, 1) + v(:, m) + v(:, m + 1))/3d0
      asum = asum + at
    end do
    xa = xa/asum
  end function area_centroid3

  !> Relative difference of the Stokes form with G and the direct form at five points of the sphere of 2 nlat^2 panels, on
  !> the surface (dist = 0) or moved by dist into the fluid along -n. The Green function is that of the solver globals.
  real(8) function sphere_difference(nlat, dist) result(rel)
    integer, intent(in) :: nlat
    real(8), intent(in) :: dist
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    type(first_order_state) :: s1, s2
    type(error_type), allocatable :: gerr
    complex(8) :: x1(6, 1), x2(6, 1), st(1), di
    complex(8), allocatable :: a(:, :), c(:, :, :), cl(:, :, :, :), e(:, :), wl(:), wd(:), wg(:, :)
    real(8) :: xo(3)
    integer :: ipt, i

    call motions(x1, x2)
    call make_sphere(5d0, [0d0, 0d0, -20d0], nlat, full)
    call load_panels(full, 0, 0, mesh)
    call geometry(mesh, geom, gerr)
    if (allocated(gerr)) error stop 'sphere_difference: build_surface_geom failed'
    if (geom%nwl /= 0) error stop 'sphere_difference: the submerged sphere has waterline segments'
    call velocity_state(mesh, w1, b1, x1, s1)
    call velocity_state(mesh, w2, b2, x2, s2)
    allocate(a(mesh%n, 1), c(3, mesh%n, 1), cl(3, 4, mesh%n, 1), e(0, 1), wl(0))
    call b2_stokes_densities(mesh, geom, s1, s2, wl, wl, a(:, 1), c(:, :, 1), e(:, 1), cl(:, :, :, 1))
    call direct_density(mesh, geom, x1, x2, wd, wg)
    rel = 0d0
    do ipt = 1, 5
      i = 1 + (ipt - 1)*(mesh%n - 1)/4
      if (ipt == 3) i = mesh%n/2 + nlat/2     ! a panel at the equator
      xo = mesh%xc(:, i) - dist*mesh%nrm(:, i)
      call b2_stokes_at_point(mesh, geom, xo, a, c, e, .false., st, cl)
      di = direct_at_point(mesh, xo, wd, wg, x1, x2)
      rel = max(rel, abs(st(1) - di)/abs(di))
    end do
  end function sphere_difference

  !> Submerged sphere (radius 5, centre depth 20, area centroids), the fields of stokes_identity_submerged, at five
  !> collocation points on the sphere: the Stokes form with G(x_i, .) against the direct form int G n.[(X.grad)u] dS.
  !> Depth 50 is gated: below b2_limit at 2048 panels and an observed order of at least b2_order from 128 to 2048 panels.
  !> Deep water is printed only: there the value and the gradient of the Green function differ by up to 0.46%
  !> (test_source_velocity), and the two forms use one each.
  subroutine test_b2_smooth(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: e(3), ed(3), order, order_deep
    integer :: m
    character(len=250) :: line

    call set_wave(50d0)
    NTHREAD = 1
    do m = 1, 3
      e(m) = sphere_difference(sph_nlat(m), 0d0)
    end do
    call set_wave(-1d0)
    NTHREAD = 1
    do m = 1, 3
      ed(m) = sphere_difference(sph_nlat(m), 0d0)
    end do
    order = log(e(1)/e(3))/log(4d0)
    order_deep = log(ed(1)/ed(3))/log(4d0)
    write(line, '(a,3es11.3,a,f6.2)') "depth 50, at 128, 512, 2048 panels: ", e, ", order ", order
    print '(4x,a)', trim(line)
    write(line, '(a,3es11.3,a,f6.2,a)') "deep water (not gated): ", ed, ", order ", order_deep
    print '(4x,a)', trim(line)
    write(line, '(a,3es11.3,a,f6.2)') "depth 50: ", e, ", order ", order
    ! Level limited by the fitted density gradient on the pole fan triangles; analytic gradient 2.8e-4 (7a ledger)
    call check(error, e(3) < b2_limit, "not below the limit at 2048 panels; " // line)
    if (allocated(error)) return
    call check(error, order >= b2_order, "observed order below the limit; " // line)
  end subroutine test_b2_smooth

  !> As rhs_b2_stokes_vs_direct_smooth with the points 10 m off the sphere, where G is smooth on S, and a Green function
  !> close to its Rankine part (deep water, omega = 0.05 rad/s): the two forms then differ by the quadrature error of
  !> smooth integrands, O(h^2), as stokes_identity_submerged shows for psi.
  subroutine test_b2_off_surface(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: e(3), order
    integer :: m
    character(len=200) :: line

    call set_wave(-1d0)
    VWAVE = 0.05d0**2/9.80665d0
    NTHREAD = 1
    do m = 1, 3
      e(m) = sphere_difference(sph_nlat(m), 10d0)
    end do
    order = log(e(1)/e(3))/log(4d0)
    write(line, '(a,3es11.3,a,f6.2)') "10 m off the surface, at 128, 512, 2048 panels: ", e, ", order ", order
    print '(4x,a)', trim(line)
    call check(error, e(3) < off_limit .and. order >= off_order, line)
  end subroutine test_b2_off_surface

  ! ------------------------------------------------------------------------------------------------ rankine_moments

  !> rankine_moments against adaptive quadrature in polar coordinates about the projected field point, on a flat quad
  !> (a trapezoid in a tilted plane) and a triangle, for field points on the panel (principal value) and at heights 0.01,
  !> 0.3 and 2 above it, one of them outside the panel outline.
  subroutine test_rankine_moments(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: quad(3, 4), tri(3, 4), nq(3), nt(3), worst, rel
    real(8) :: m0, g0(3), m1(3, 3), mg(3), b0, bg(3), bm(3, 3), bq(3), xc(3), x(3)
    real(8), parameter :: heights(4) = [0d0, 0.01d0, 0.3d0, 2d0]
    real(8) :: offsets(2, 3)
    integer :: ip, ih, io, nv
    character(len=200) :: line

    call tilted_panel([0d0, 0d0], [1.2d0, 0.1d0], [0.9d0, 1d0], [-0.1d0, 0.8d0], 4, quad, nq)
    call tilted_panel([0d0, 0d0], [1d0, 0.2d0], [0.3d0, 0.9d0], [0d0, 0d0], 3, tri, nt)
    worst = 0d0
    do ip = 1, 2
      if (ip == 1) then
        nv = 4
        xc = sum(quad, 2)/4d0
        offsets = reshape([0d0, 0d0, 0.2d0, -0.15d0, 1.5d0, 0.4d0], [2, 3])
      else
        nv = 3
        xc = sum(tri(:, 1:3), 2)/3d0
        offsets = reshape([0d0, 0d0, 0.1d0, 0.12d0, -0.7d0, 0.3d0], [2, 3])
      end if
      do io = 1, 3
        do ih = 1, 4
          if (io == 3 .and. ih == 1) then
            ! in the plane and outside the outline: the integrand is regular, plain 2D Gauss on the panel's triangles
            if (ip == 1) then
              x = xc + offsets(1, io)*in_plane(quad, nq, 1) + offsets(2, io)*in_plane(quad, nq, 2)
              call rankine_moments(quad, nv, nq, x, xc, m0, g0, m1, mg)
              call gauss_moments(quad, nv, nq, x, xc, b0, bg, bm, bq)
            else
              x = xc + offsets(1, io)*in_plane(tri, nt, 1) + offsets(2, io)*in_plane(tri, nt, 2)
              call rankine_moments(tri, nv, nt, x, xc, m0, g0, m1, mg)
              call gauss_moments(tri, nv, nt, x, xc, b0, bg, bm, bq)
            end if
            rel = max(abs(m0 - b0)/abs(b0), norm2(g0 - bg)/max(norm2(bg), abs(b0)), &
                      maxval(abs(m1 - bm))/max(maxval(abs(bm)), abs(b0)), norm2(mg - bq)/max(norm2(bq), abs(b0)))
            worst = max(worst, rel)
            if (rel > moment_limit) print '(4x,a,2i3,es10.2)', "rankine_moments: panel, case, relative difference ", &
              ip, 10*io + ih, rel
            cycle
          end if
          if (ip == 1) then
            x = xc + offsets(1, io)*in_plane(quad, nq, 1) + offsets(2, io)*in_plane(quad, nq, 2) + heights(ih)*nq
            call rankine_moments(quad, nv, nq, x, xc, m0, g0, m1, mg)
            call brute_moments(quad, nv, nq, x, xc, b0, bg, bm, bq)
          else
            x = xc + offsets(1, io)*in_plane(tri, nt, 1) + offsets(2, io)*in_plane(tri, nt, 2) + heights(ih)*nt
            call rankine_moments(tri, nv, nt, x, xc, m0, g0, m1, mg)
            call brute_moments(tri, nv, nt, x, xc, b0, bg, bm, bq)
          end if
          rel = max(abs(m0 - b0)/abs(b0), norm2(g0 - bg)/max(norm2(bg), abs(b0)), &
                    maxval(abs(m1 - bm))/max(maxval(abs(bm)), abs(b0)), norm2(mg - bq)/max(norm2(bq), abs(b0)))
          worst = max(worst, rel)
          if (rel > moment_limit) print '(4x,a,2i3,es10.2)', "rankine_moments: panel, case, relative difference ", ip, &
            10*io + ih, rel
        end do
      end do
    end do
    write(line, '(a,es10.2)') "closed form against adaptive quadrature, worst relative difference ", worst
    print '(4x,a)', trim(line)
    call check(error, worst <= moment_limit, line)
  end subroutine test_rankine_moments

  !> A flat panel from 2D vertices placed in a tilted plane, and its unit normal.
  subroutine tilted_panel(p1, p2, p3, p4, nv, xv, nrm)
    real(8), intent(in) :: p1(2), p2(2), p3(2), p4(2)
    integer, intent(in) :: nv
    real(8), intent(out) :: xv(3, 4), nrm(3)
    real(8) :: e1(3), e2(3), o(3), p(2, 4)
    integer :: k
    o = [1d0, -2d0, -3d0]
    e1 = [0.8d0, 0.6d0, 0d0]
    e2 = [-0.36d0, 0.48d0, 0.8d0]
    nrm = [e1(2)*e2(3) - e1(3)*e2(2), e1(3)*e2(1) - e1(1)*e2(3), e1(1)*e2(2) - e1(2)*e2(1)]
    p(:, 1) = p1; p(:, 2) = p2; p(:, 3) = p3; p(:, 4) = p4
    xv = 0d0
    do k = 1, nv
      xv(:, k) = o + p(1, k)*e1 + p(2, k)*e2
    end do
  end subroutine tilted_panel

  !> In-plane unit vector m (1 or 2) of the panel, from its first edge.
  function in_plane(xv, nrm, m) result(t)
    real(8), intent(in) :: xv(3, 4), nrm(3)
    integer, intent(in) :: m
    real(8) :: t(3)
    t = xv(:, 2) - xv(:, 1)
    t = t/norm2(t)
    if (m == 2) t = [nrm(2)*t(3) - nrm(3)*t(2), nrm(3)*t(1) - nrm(1)*t(3), nrm(1)*t(2) - nrm(2)*t(1)]
  end function in_plane

  !> The moments of rankine_moments by adaptive Gauss-Legendre quadrature over the triangles (1, k, k+1) of the panel, in
  !> the collapsed coordinates s in [0, 1], t in [0, 1 - s]: for x away from the panel, where the integrand is regular.
  subroutine gauss_moments(xv, nv, nrm, x, xc, b0, bg, bm, bq)
    real(8), intent(in) :: xv(3, 4), nrm(3), x(3), xc(3)
    integer, intent(in) :: nv
    real(8), intent(out) :: b0, bg(3), bm(3, 3), bq(3)
    real(8) :: e1(3), e2(3), pa(3), pb(3), pc(3), jac, ss, tot(9), res(9)
    integer :: k, a, b
    e1 = xv(:, 2) - xv(:, 1)
    e1 = e1/norm2(e1)
    e2 = [nrm(2)*e1(3) - nrm(3)*e1(2), nrm(3)*e1(1) - nrm(1)*e1(3), nrm(1)*e1(2) - nrm(2)*e1(1)]
    tot = 0d0
    do k = 2, nv - 1
      pa = xv(:, 1)
      pb = xv(:, k)
      pc = xv(:, k + 1)
      jac = norm2([(pb(2) - pa(2))*(pc(3) - pa(3)) - (pb(3) - pa(3))*(pc(2) - pa(2)), &
                   (pb(3) - pa(3))*(pc(1) - pa(1)) - (pb(1) - pa(1))*(pc(3) - pa(3)), &
                   (pb(1) - pa(1))*(pc(2) - pa(2)) - (pb(2) - pa(2))*(pc(1) - pa(1))])
      call adapt(s_fun, 0d0, 1d0, 1d-13, 9, res, 0)
      tot = tot + res
    end do
    b0 = tot(1)
    bg = tot(2)*e1 + tot(3)*e2
    bq = tot(8)*e1 + tot(9)*e2
    do a = 1, 3
      do b = 1, 3
        bm(a, b) = (tot(4)*e1(a) + tot(5)*e2(a))*e1(b) + (tot(6)*e1(a) + tot(7)*e2(a))*e2(b)
      end do
    end do
  contains
    recursive function s_fun(sv) result(v)
      real(8), intent(in) :: sv
      real(8) :: v(9)
      ss = sv
      call adapt(t_fun, 0d0, 1d0 - sv, 1d-14, 9, v, 0)
    end function s_fun
    recursive function t_fun(tv) result(v)
      real(8), intent(in) :: tv
      real(8) :: v(9), y(3), d(3), r, ga(2), q(2)
      y = pa + ss*(pb - pa) + tv*(pc - pa)
      d = y - x
      r = norm2(d)
      ga = -[dot_product(d, e1), dot_product(d, e2)]/r**3
      q = [dot_product(y - xc, e1), dot_product(y - xc, e2)]
      v(1) = 1d0/r
      v(2:3) = ga
      v(4) = ga(1)*q(1)
      v(5) = ga(2)*q(1)
      v(6) = ga(1)*q(2)
      v(7) = ga(2)*q(2)
      v(8:9) = q/r
      v = v*jac
    end function t_fun
  end subroutine gauss_moments

  !> The moments of rankine_moments by adaptive Gauss-Legendre quadrature: the panel as the triangles from the projected
  !> point X to each edge (signed by the side of X), each in polar coordinates (theta, rho) about X. In the plane of x
  !> (Z = 0) the rho integrals are done in closed form, with the principal value of int grad(1/r) taken as
  !> -e_a ln(rho_max) per direction (the ln of the inner radius cancels over the full circle, X inside the panel).
  subroutine brute_moments(xv, nv, nrm, x, xc, b0, bg, bm, bq)
    real(8), intent(in) :: xv(3, 4), nrm(3), x(3), xc(3)
    integer, intent(in) :: nv
    real(8), intent(out) :: b0, bg(3), bm(3, 3), bq(3)
    real(8) :: xp(3), z, q(3), u1(3), u2(3), nus(3), hk, delta, res(9), tot(9), p1(3), p2(3), ed(3), sg
    real(8) :: e1(3), e2(3), theta
    integer :: k, a, b

    z = dot_product(x - xv(:, 1), nrm)
    xp = x - z*nrm
    q = xp - (xc - dot_product(xc - xv(:, 1), nrm)*nrm)      ! X - xc in the plane
    e1 = xv(:, 2) - xv(:, 1)
    e1 = e1/norm2(e1)
    e2 = [nrm(2)*e1(3) - nrm(3)*e1(2), nrm(3)*e1(1) - nrm(1)*e1(3), nrm(1)*e1(2) - nrm(2)*e1(1)]
    tot = 0d0
    do k = 1, nv
      p1 = xv(:, k)
      p2 = xv(:, mod(k, nv) + 1)
      ed = (p2 - p1)/norm2(p2 - p1)
      ! normal of the edge line in the plane, pointing from X to the line
      nus = (p1 - xp) - dot_product(p1 - xp, ed)*ed
      hk = norm2(nus)
      if (hk < 1d-12) cycle
      nus = nus/hk
      ! side of X: + if X is on the inner side of this edge (the polygon is convex, vertices in order)
      sg = sign(1d0, dot_product(nus, sum(xv(:, 1:nv), 2)/nv - p1)*(-1d0))
      u1 = (p1 - xp)/norm2(p1 - xp)
      u2 = (p2 - xp) - dot_product(p2 - xp, u1)*u1
      u2 = u2/norm2(u2)
      delta = acos(max(-1d0, min(1d0, dot_product((p2 - xp)/norm2(p2 - xp), u1))))
      call adapt(theta_fun, 0d0, delta, 1d-13, 9, res, 0)
      tot = tot + sg*res
    end do
    b0 = tot(1)
    bg = tot(2)*e1 + tot(3)*e2
    bq = tot(8)*e1 + tot(9)*e2
    do a = 1, 3
      do b = 1, 3
        bm(a, b) = (tot(4)*e1(a) + tot(5)*e2(a))*e1(b) + (tot(6)*e1(a) + tot(7)*e2(a))*e2(b)
      end do
    end do

  contains

    !> The rho integral at direction theta: values m0, g0 (2), m1 (a, b) in the order (1,1), (2,1), (1,2), (2,2), and
    !> int (1/r)(xi - xc) (2).
    recursive function theta_fun(t) result(v)
      real(8), intent(in) :: t
      real(8) :: v(9), dir(3), rmax, ea(2), qb(2)
      dir = cos(t)*u1 + sin(t)*u2
      rmax = hk/dot_product(dir, nus)
      ea = [dot_product(dir, e1), dot_product(dir, e2)]
      qb = [dot_product(q, e1), dot_product(q, e2)]
      theta = t
      if (abs(z) <= 1d-12) then      ! in the plane, to rounding
        ! 1/r = 1/rho: int rho/rho = rmax; d_a(1/r) rho drho = -e_a/rho (PV: -e_a ln rmax); with (rho e_b + q_b)
        v(1) = rmax
        v(2:3) = -ea*log(rmax)
        v(4) = -ea(1)*(ea(1)*rmax + qb(1)*log(rmax))
        v(5) = -ea(2)*(ea(1)*rmax + qb(1)*log(rmax))
        v(6) = -ea(1)*(ea(2)*rmax + qb(2)*log(rmax))
        v(7) = -ea(2)*(ea(2)*rmax + qb(2)*log(rmax))
        v(8:9) = ea*0.5d0*rmax*rmax + qb*rmax
      else
        call adapt(rho_fun, 0d0, rmax, 1d-14, 9, v, 0)
      end if
    end function theta_fun

    recursive function rho_fun(r) result(v)
      real(8), intent(in) :: r
      real(8) :: v(9), dir(3), ea(2), qb(2), rr
      dir = cos(theta)*u1 + sin(theta)*u2
      ea = [dot_product(dir, e1), dot_product(dir, e2)]
      qb = [dot_product(q, e1), dot_product(q, e2)]
      rr = sqrt(r*r + z*z)
      ! grad_xi(1/r) in the plane = -(xi - x)_inplane/r^3 = -r e/rr^3; Jacobian r
      v(1) = r/rr
      v(2:3) = -r*r*ea/rr**3
      v(4) = -r*r*ea(1)*(r*ea(1) + qb(1))/rr**3
      v(5) = -r*r*ea(2)*(r*ea(1) + qb(1))/rr**3
      v(6) = -r*r*ea(1)*(r*ea(2) + qb(2))/rr**3
      v(7) = -r*r*ea(2)*(r*ea(2) + qb(2))/rr**3
      v(8:9) = r*(r*ea + qb)/rr
    end function rho_fun

  end subroutine brute_moments

  !> Adaptive Gauss-Legendre: 10 points on [a, b] against 10 on each half, bisecting until they agree to tol times the
  !> larger of 1 and the result on [a, b]. The tolerance is not halved on the halves: GL10 on an analytic integrand gains
  !> many digits per bisection, so the sum of the leaves stays near tol. At most 18 levels.
  recursive subroutine adapt(f, a, b, tol, n, res, depth)
    interface
      function f(t) result(v)
        real(8), intent(in) :: t
        real(8) :: v(9)
      end function f
    end interface
    real(8), intent(in) :: a, b, tol
    integer, intent(in) :: n, depth
    real(8), intent(out) :: res(n)
    real(8) :: whole(9), left(9), right(9), m, rl(9), rr(9)
    m = 0.5d0*(a + b)
    whole = gl10(f, a, b)
    left = gl10(f, a, m)
    right = gl10(f, m, b)
    if (maxval(abs(left + right - whole)) <= tol*max(1d0, maxval(abs(whole))) .or. depth >= 18) then
      res = left + right
    else
      call adapt(f, a, m, tol, n, rl, depth + 1)
      call adapt(f, m, b, tol, n, rr, depth + 1)
      res = rl + rr
    end if
  end subroutine adapt

  recursive function gl10(f, a, b) result(s)
    interface
      function f(t) result(v)
        real(8), intent(in) :: t
        real(8) :: v(9)
      end function f
    end interface
    real(8), intent(in) :: a, b
    real(8) :: s(9)
    real(8), parameter :: xg(5) = [0.1488743389816312d0, 0.4333953941292472d0, 0.6794095682990244d0, &
                                   0.8650633666889845d0, 0.9739065285171717d0]
    real(8), parameter :: wg(5) = [0.2955242247147529d0, 0.2692667193099963d0, 0.2190863625159820d0, &
                                   0.1494513491505806d0, 0.0666713443086881d0]
    real(8) :: c, h
    integer :: k
    c = 0.5d0*(a + b)
    h = 0.5d0*(b - a)
    s = 0d0
    do k = 1, 5
      s = s + wg(k)*(f(c + h*xg(k)) + f(c - h*xg(k)))
    end do
    s = s*h
  end function gl10

  !> direct_rhs_b2 on the half cylinder (symmetry in y, then in x) against the whole one, with the waterline: the blocks
  !> of the half mesh, put back on the images, are the right side of the whole mesh.
  subroutine test_b2_half(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, whole, half
    type(surface_geom) :: geom_w, geom_h
    type(first_order_state) :: s1w, s2w, s1h, s2h
    complex(8) :: x1(6, 1), x2(6, 1)
    complex(8), allocatable :: wl1w(:), wl2w(:), wl1h(:), wl2h(:), rw(:, :), rh(:, :), back(:)
    real(8) :: worst
    integer :: icase, p, i, im, ip, nb
    integer, parameter :: sx(2) = [0, 1], sy(2) = [1, 0]
    character(len=200) :: line

    call make_cylinder(5d0, 5d0, 16, 4, 3, full, .true.)
    call motions(x1, x2)
    call set_wave(-1d0)
    NTHREAD = 4
    call load_panels(full, 0, 0, whole)
    call geometry(whole, geom_w, error)
    if (allocated(error)) return
    call make_state(whole, geom_w, w1, b1, -1d0, x1, s1w, wl1w)
    call make_state(whole, geom_w, w2, b2, -1d0, x2, s2w, wl2w)
    allocate(rw(whole%nbase, 1))
    call direct_rhs_b2(whole, geom_w, s1w, s2w, wl1w, wl2w, rw, stored=.false.)
    worst = 0d0
    do icase = 1, 2
      call load_panels(full, sx(icase), sy(icase), half)
      NTHREAD = 4
      call geometry(half, geom_h, error)
      if (allocated(error)) return
      call make_state(half, geom_h, w1, b1, -1d0, x1, s1h, wl1h)
      call make_state(half, geom_h, w2, b2, -1d0, x2, s2h, wl2h)
      nb = half%nbase
      allocate(rh(nb, half%nsys), back(half%n))
      call direct_rhs_b2(half, geom_h, s1h, s2h, wl1h, wl2h, rh, stored=.false.)
      ! The full panel (im, i) of the half mesh: sum_ip RXY(im, ip) rh(i, ip) / nsys, as RADIATION_SOLVER
      do im = 1, half%nsys
        do i = 1, nb
          back((im - 1)*nb + i) = (0d0, 0d0)
          do ip = 1, half%nsys
            back((im - 1)*nb + i) = back((im - 1)*nb + i) + RXY(im, ip)*rh(i, ip)
          end do
          back((im - 1)*nb + i) = back((im - 1)*nb + i)/half%nsys
        end do
      end do
      ! Match each panel of the half mesh to the panel of the whole mesh at the same centroid
      do p = 1, half%n
        i = minloc(norm2(whole%xc - spread(half%xc(:, p), 2, whole%n), dim=1), 1)
        worst = max(worst, abs(back(p) - rw(i, 1))/maxval(abs(rw)))
      end do
      deallocate(rh, back)
    end do
    write(line, '(a,es10.2)') "half mesh against whole, worst difference relative to the largest entry ", worst
    print '(4x,a)', trim(line)
    call check(error, worst <= half_limit, line)
    ISX = 0; ISY = 0; NSYS = 1
  end subroutine test_b2_half

  !> direct_rhs_b2 from the stored integrals of CALGREEN against green_panel_integrals: the cylinder of the half test,
  !> whole and half (ISY = 1), depth 10 and deep water, and two cylinders in the multi-body arrays (CALGREEN_MULTI).
  subroutine test_b2_stored(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    type(first_order_state) :: s1, s2
    complex(8) :: x1(6, 1), x2(6, 1), xm1(6, 2), xm2(6, 2)
    complex(8), allocatable :: wl1(:), wl2(:), ra(:, :), rb(:, :)
    real(8) :: worst, d
    integer :: icase, isy_in, ne
    real(8), parameter :: depths(2) = [10d0, -1d0]
    character(len=200) :: line

    call make_cylinder(5d0, 5d0, 16, 4, 3, full, .true.)
    call motions(x1, x2)
    worst = 0d0
    do icase = 1, 6
      d = depths(1 + mod(icase - 1, 2))
      isy_in = merge(0, 1, icase <= 2 .or. icase == 5)
      call set_wave(d)
      NTHREAD = 4
      if (icase <= 4) then
        call load_panels(full, 0, isy_in, mesh)
        if (allocated(CGRN)) deallocate(CGRN)
        if (allocated(RKBN)) deallocate(RKBN)
        allocate(CGRN(NELEM, NELEM, NSYS, 4), RKBN(NELEM, NELEM, NSYS, 4))
        call CALGREEN
        call geometry(mesh, geom, error)
        if (allocated(error)) return
        call make_state(mesh, geom, w1, b1, d, x1, s1, wl1)
        call make_state(mesh, geom, w2, b2, d, x2, s2, wl2)
      else
        call load_multi(isy_in, mesh)
        ne = mesh%nbase
        TNELEM = ne
        if (allocated(CGRN_MULTI_COMB)) deallocate(CGRN_MULTI_COMB)
        if (allocated(RKBN_MULTI_COMB)) deallocate(RKBN_MULTI_COMB)
        allocate(CGRN_MULTI_COMB(ne, ne, NSYS, 4), RKBN_MULTI_COMB(ne, ne, NSYS, 4))
        call CALGREEN_MULTI
        call geometry(mesh, geom, error)
        if (allocated(error)) return
        xm1(:, 1) = x1(:, 1)
        xm1(:, 2) = 0.5d0*x2(:, 1)
        xm2(:, 1) = x2(:, 1)
        xm2(:, 2) = (0d0, 1d0)*x1(:, 1)
        call make_state(mesh, geom, w1, b1, d, xm1, s1, wl1)
        call make_state(mesh, geom, w2, b2, d, xm2, s2, wl2)
      end if
      allocate(ra(mesh%nbase, mesh%nsys), rb(mesh%nbase, mesh%nsys))
      call direct_rhs_b2(mesh, geom, s1, s2, wl1, wl2, ra)
      call direct_rhs_b2(mesh, geom, s1, s2, wl1, wl2, rb, stored=.false.)
      worst = max(worst, maxval(abs(ra - rb))/maxval(abs(rb)))
      deallocate(ra, rb)
    end do
    write(line, '(a,es10.2)') "stored against fresh panel integrals, worst difference relative to the largest entry ", worst
    print '(4x,a)', trim(line)
    call check(error, worst <= stored_limit, line)
    ISX = 0; ISY = 0; NSYS = 1
  end subroutine test_b2_stored

  !> The sharp-cylinder fixture (radius 10, draft 10, 448 panels): the panels at the bottom rim, the bottom row of the wall
  !> and the outer ring of the bottom (2 x 32), keep the one-point rule (cl = 0); every other panel gets a linear part.
  subroutine test_b2_sharp_edge(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    type(first_order_state) :: s1, s2
    complex(8) :: x1(6, 1), x2(6, 1)
    complex(8), allocatable :: wl1(:), wl2(:), a(:), c(:, :), e(:), cl(:, :, :)
    logical, allocatable :: mask(:)
    integer :: i, nmask, nlinear, wrong
    character(len=200) :: line

    call make_cylinder(10d0, 10d0, 32, 8, 6, full, .true.)
    call set_wave(-1d0)
    call load_panels(full, 0, 0, mesh)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    call motions(x1, x2)
    call make_state(mesh, geom, w1, b1, -1d0, x1, s1, wl1)
    call make_state(mesh, geom, w2, b2, -1d0, x2, s2, wl2)
    call sharp_edge_panels(mesh, DEFAULT_EDGE_DEG, mask)
    allocate(a(mesh%n), c(3, mesh%n), e(geom%nwl), cl(3, 4, mesh%n))
    call b2_stokes_densities(mesh, geom, s1, s2, wl1, wl2, a, c, e, cl, mask)
    nmask = count(mask)
    nlinear = 0
    wrong = 0
    do i = 1, mesh%n
      if (mask(i)) then
        if (any(cl(:, :, i) /= (0d0, 0d0))) wrong = wrong + 1
        ! a rim panel has a node on the rim, z = -10 and radius 10
        if (.not. any(abs(mesh%xnode(3, mesh%ncon(1:mesh%ncn(i), i)) + 10d0) < 1d-9 .and. &
                      abs(norm2(mesh%xnode(1:2, mesh%ncon(1:mesh%ncn(i), i)), dim=1) - 10d0) < 1d-9)) wrong = wrong + 1
      else if (any(cl(:, :, i) /= (0d0, 0d0))) then
        nlinear = nlinear + 1
      end if
    end do
    write(line, '(a,i4,a,i4,a,i4,a,i3)') "panels ", mesh%n, ", at the edge ", nmask, ", with a linear part ", nlinear, &
      ", wrong ", wrong
    print '(4x,a)', trim(line)
    call check(error, nmask == 64 .and. wrong == 0 .and. nlinear == mesh%n - nmask, line)
  end subroutine test_b2_sharp_edge

  ! ------------------------------------------------------------------------------------------------ the solve at Omega

  !> The radiation problem at omega on the ne base panels of load_panels (multi = .false.) or load_multi (two bodies), as
  !> SolveRadiation and SolveRadiationMulti of HAMS_Prog run it with IRSP = 0 and the ISOLV of the caller: the stored Green
  !> function, AMAT (factored when ISOLV = 1), BRMAT and MXPOT at omega. The wave of the solver globals is set to omega
  !> by SetWaveFromOmega, as RadiationOnlySolves does it; the depth is that of the last set_wave.
  subroutine radiation_at(omega, multi, ne)
    real(8), intent(in) :: omega
    logical, intent(in) :: multi
    integer, intent(in) :: ne
    call SetWaveFromOmega(omega)
    TNELEM = ne
    if (.not. multi) then
      if (allocated(CGRN)) deallocate(CGRN)
      if (allocated(RKBN)) deallocate(RKBN)
      if (allocated(AMAT)) deallocate(AMAT)
      if (allocated(BRMAT)) deallocate(BRMAT)
      if (allocated(IPIV)) deallocate(IPIV)
      allocate(CGRN(ne, ne, NSYS, 4), RKBN(ne, ne, NSYS, 4), AMAT(ne, ne, NSYS), BRMAT(ne, 6, NSYS), IPIV(ne, NSYS))
      call CALGREEN
      call ASSB_LEFT(AMAT, IPIV, ne, NSYS)
      call ASSB_RBC(BRMAT, ne, NSYS)
      call RADIATION_SOLVER(AMAT, BRMAT, IPIV, MXPOT, ne, NSYS)
    else
      if (allocated(CGRN_MULTI_COMB)) deallocate(CGRN_MULTI_COMB)
      if (allocated(RKBN_MULTI_COMB)) deallocate(RKBN_MULTI_COMB)
      if (allocated(AMAT_MULTI)) deallocate(AMAT_MULTI)
      if (allocated(BRMAT_MULTI)) deallocate(BRMAT_MULTI)
      if (allocated(IPIV_MULTI_COMB)) deallocate(IPIV_MULTI_COMB)
      allocate(CGRN_MULTI_COMB(ne, ne, NSYS, 4), RKBN_MULTI_COMB(ne, ne, NSYS, 4), AMAT_MULTI(ne, ne, NSYS), &
               BRMAT_MULTI(ne, 12, NSYS), IPIV_MULTI_COMB(ne, NSYS))
      call CALGREEN_MULTI
      call ASSB_LEFT_MULTI(AMAT_MULTI, IPIV_MULTI_COMB, ne, NSYS)
      BRMAT_MULTI = (0d0, 0d0)
      call ASSB_RBC_MULTI(BRMAT_MULTI(:, 1:6, :), ne, NSYS, 1, ne/2)
      call ASSB_RBC_MULTI(BRMAT_MULTI(:, 7:12, :), ne, NSYS, ne/2 + 1, ne)
      call RADIATION_SOLVER_MULTI(AMAT_MULTI, BRMAT_MULTI, IPIV_MULTI_COMB, MXPOT_MULTI_COMB, ne, NSYS)
    end if
  end subroutine radiation_at

  !> The QTF store of two states (frequency 1 and 2, one heading) and the pairs (1, 2) and (2, 1), Omega = +-(w1 - w2).
  subroutine two_state_store(mesh, geom, depth, x1, x2, store, store_wl, pairs, wa, wb)
    type(panel_mesh), intent(in) :: mesh
    type(surface_geom), intent(in) :: geom
    real(8), intent(in) :: depth
    complex(8), intent(in) :: x1(:, :), x2(:, :)
    type(first_order_state), allocatable, intent(out) :: store(:, :)
    complex(8), allocatable, intent(out) :: store_wl(:, :, :)
    integer, allocatable, intent(out) :: pairs(:, :)
    real(8), intent(in), optional :: wa, wb   ! the two frequencies, w1 and w2 by default
    complex(8), allocatable :: wl1(:), wl2(:)
    real(8) :: fa, fb
    fa = w1
    fb = w2
    if (present(wa)) fa = wa
    if (present(wb)) fb = wb
    allocate(store(2, 1), store_wl(geom%nwl, 2, 1))
    call make_state(mesh, geom, fa, b1, depth, x1, store(1, 1), wl1)
    call make_state(mesh, geom, fb, b2, depth, x2, store(2, 1), wl2)
    store_wl(:, 1, 1) = wl1
    store_wl(:, 2, 1) = wl2
    pairs = reshape([1, 1, 2, 1, 2, 1, 1, 1], [4, 2])
  end subroutine two_state_store

  !> The unit normals n_k of each mode as panel fields in the layout of ASSB_RHS_FIELD: V(i, k, a) on image a of base
  !> panel i, n_k (mesh%nj) on the panels of the body of mode k and 0 on the others.
  subroutine unit_normal_fields(mesh, vf)
    type(panel_mesh), intent(in) :: mesh
    complex(8), allocatable, intent(out) :: vf(:, :, :)
    integer :: f, k, kb
    allocate(vf(mesh%nbase, 6*mesh%nbody, mesh%nsys))
    vf = (0d0, 0d0)
    do f = 1, mesh%n
      do k = 1, 6*mesh%nbody
        kb = (k - 1)/6 + 1
        if (mesh%body(f) == kb) vf(mesh%base(f), k, mesh%image(f)) = mesh%nj(k - 6*(kb - 1), f)
      end do
    end do
  end subroutine unit_normal_fields

  !> M4: the B2 right side of direct_rhs_b2 (Stokes form with G, stored integrals) against ASSB_RHS_FIELD with the
  !> analytic panel values B2 = -n.[(X.grad)u], on the sphere of the smooth test (512 panels, depth 50). Both are
  !> int G B2 dS at the collocation points, so they agree to the mesh error, and a sign error in either would make
  !> them differ by twice the right side. This pins the sign of B2 on the path where it is added to the other parts of b.
  subroutine test_b2_sign(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    type(first_order_state) :: s1, s2
    complex(8) :: x1(6, 1), x2(6, 1)
    complex(8), allocatable :: wd(:), wg(:, :), vf(:, :, :), rf(:, :, :), r(:, :), wl(:)
    real(8) :: rel, rel_flipped
    character(len=200) :: line

    call set_wave(50d0)
    NTHREAD = 4
    call make_sphere(5d0, [0d0, 0d0, -20d0], sph_nlat(2), full)
    call load_panels(full, 0, 0, mesh)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    if (allocated(CGRN)) deallocate(CGRN)
    if (allocated(RKBN)) deallocate(RKBN)
    allocate(CGRN(NELEM, NELEM, NSYS, 4), RKBN(NELEM, NELEM, NSYS, 4))
    call CALGREEN
    call motions(x1, x2)
    call velocity_state(mesh, w1, b1, x1, s1)
    call velocity_state(mesh, w2, b2, x2, s2)
    call direct_density(mesh, geom, x1, x2, wd, wg)
    allocate(vf(mesh%nbase, 1, 1), rf(mesh%nbase, 1, 1), r(mesh%nbase, 1), wl(0))
    vf(:, 1, 1) = -wd
    call ASSB_RHS_FIELD(vf, rf, mesh%nbase, 1, 1)
    call direct_rhs_b2(mesh, geom, s1, s2, wl, wl, r)
    rel = maxval(abs(rf(:, 1, 1) - r(:, 1)))/maxval(abs(r))
    rel_flipped = maxval(abs(rf(:, 1, 1) + r(:, 1)))/maxval(abs(r))
    write(line, '(a,es10.2,a,es10.2)') "B2 by the Stokes form against ASSB_RHS_FIELD of the panel values, relative ", rel, &
      "; with the sign flipped ", rel_flipped
    print '(4x,a)', trim(line)
    call check(error, rel < b2_sign_limit .and. rel_flipped > 1d0, line)
  end subroutine test_b2_sign

  !> With b = 0 the load is zero: both bodies fixed (xi = 0, so B1 = B2 = B3 = 0) and phi_I2 = 0 (the test hook of
  !> potential_pair), on the production path direct_at_omega, for Omega > 0 and < 0. With the motions the load is not zero.
  subroutine test_direct_zero_forcing(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    type(first_order_state), allocatable :: store(:, :)
    complex(8), allocatable :: store_wl(:, :, :), t(:, :)
    integer, allocatable :: pairs(:, :)
    complex(8) :: zx(6, 1), x1(6, 1), x2(6, 1)
    character(len=200) :: line

    call make_cylinder(5d0, 5d0, 16, 4, 3, full, .true.)
    call set_wave(10d0)
    NTHREAD = 4
    ISOLV = 1
    call load_panels(full, 0, 0, mesh)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    call radiation_at(w1 - w2, .false., mesh%nbase)
    call direct_setup(mesh, DEFAULT_EDGE_DEG)
    zx = (0d0, 0d0)
    call two_state_store(mesh, geom, 10d0, zx, zx, store, store_wl, pairs)
    allocate(t(6, 2))
    call direct_at_omega(mesh, geom, w1 - w2, store, store_wl, pairs, t, phi2_constant=(0d0, 0d0))
    write(line, '(a,es10.2)') "b = 0: largest |T_dir| ", maxval(abs(t))
    print '(4x,a)', trim(line)
    call check(error, all(t == (0d0, 0d0)), line)
    if (allocated(error)) return
    call motions(x1, x2)
    call two_state_store(mesh, geom, 10d0, x1, x2, store, store_wl, pairs)
    call direct_at_omega(mesh, geom, w1 - w2, store, store_wl, pairs, t, phi2_constant=(0d0, 0d0))
    write(line, '(a,es10.2)') "with the motions: smallest pair max |T_dir| ", minval(maxval(abs(t), dim=1))
    print '(4x,a)', trim(line)
    call check(error, minval(maxval(abs(t), dim=1)) > 0d0, line)
  end subroutine test_direct_zero_forcing

  !> With b = n_k the solution of direct_solve (right side by ASSB_RHS_FIELD, the factored AMAT at Omega) is psi_k of the
  !> radiation solve at Omega (RADIATION_SOLVER, radiation_psi): the same matrix and the same right side, so to rounding.
  !> The cylinder of the half test, whole at depth 10, half (ISY = 1) in deep water, two cylinders in the multi-body arrays
  !> at depth 10, whole with ISOLV = 2 (AMAT unfactored by ASSB_LEFT, factored by direct_factor), and the two cylinders
  !> as half meshes (ISY = 1), whose images direct_solve expands with RXY_MULTI.
  subroutine test_direct_unit_normal(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    complex(8), allocatable :: vf(:, :, :), rf(:, :, :), phi(:, :), psi(:, :)
    real(8) :: worst, w, d
    real(8), parameter :: om = 0.4d0
    integer :: icase, nb
    logical :: multi
    character(len=200) :: line

    call make_cylinder(5d0, 5d0, 16, 4, 3, full, .true.)
    worst = 0d0
    do icase = 1, 5
      d = merge(-1d0, 10d0, icase == 2)
      multi = icase == 3 .or. icase == 5
      call set_wave(d)
      NTHREAD = 4
      ISOLV = 1
      if (multi) then
        call load_multi(merge(1, 0, icase == 5), mesh)
      else
        call load_panels(full, 0, merge(1, 0, icase == 2), mesh)
      end if
      nb = mesh%nbase
      call radiation_at(om, multi, nb)
      allocate(psi(mesh%n, 6*mesh%nbody), phi(mesh%n, 6*mesh%nbody), rf(nb, 6*mesh%nbody, mesh%nsys))
      if (multi) then
        call radiation_psi(mesh, om, MXPOT_MULTI_COMB, psi)
      else
        call radiation_psi(mesh, om, MXPOT, psi)
      end if
      if (icase == 4) then
        ISOLV = 2
        call ASSB_LEFT(AMAT, IPIV, nb, NSYS)
        call direct_factor(mesh)
      end if
      call unit_normal_fields(mesh, vf)
      if (multi) then
        call ASSB_RHS_FIELD_MULTI(vf, rf, nb, mesh%nsys, 6*mesh%nbody)
      else
        call ASSB_RHS_FIELD(vf, rf, nb, mesh%nsys, 6*mesh%nbody)
      end if
      call direct_solve(mesh, rf, phi, multi)
      w = maxval(abs(phi - psi))/maxval(abs(psi))
      write(line, '(a,i0,a,i0,a,i0,a,es10.2)') "case ", icase, ": ", mesh%n, " panels, nsys ", mesh%nsys, &
        ", |phi - psi_k| relative ", w
      print '(4x,a)', trim(line)
      worst = max(worst, w)
      deallocate(psi, phi, rf, vf)
    end do
    ISOLV = 1
    ISX = 0; ISY = 0; NSYS = 1
    write(line, '(a,es10.2)') "b = n_k against psi_k, worst relative difference ", worst
    call check(error, worst <= unit_normal_limit, line)
  end subroutine test_direct_unit_normal

  !> M9: the sharp-edge guard on the production path. On the sharp cylinder (448 panels), direct_setup at the default
  !> angle masks the 64 rim panels, and the load of direct_at_omega changes from the unmasked run (direct_setup at 180
  !> degrees, no panel masked) by exactly the B2 linear parts of the masked panels: their right side through direct_solve
  !> and the load integral.
  subroutine test_direct_edge_guard(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    type(first_order_state), allocatable :: store(:, :)
    complex(8), allocatable :: store_wl(:, :, :), tm(:, :), tn(:, :), a(:, :), c(:, :, :), e(:, :), clm(:, :, :, :), &
                               cln(:, :, :, :), rm(:, :, :), rn(:, :, :), dphi(:, :)
    integer, allocatable :: pairs(:, :)
    logical, allocatable :: mask(:)
    complex(8) :: x1(6, 1), x2(6, 1), delta(6)
    real(8) :: worst, om
    integer :: nmask, nnone, f, j
    character(len=250) :: line

    call make_cylinder(10d0, 10d0, 32, 8, 6, full, .true.)
    call set_wave(-1d0)
    NTHREAD = 4
    ISOLV = 1
    call load_panels(full, 0, 0, mesh)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    om = w1 - w2
    call radiation_at(om, .false., mesh%nbase)
    call motions(x1, x2)
    call two_state_store(mesh, geom, -1d0, x1, x2, store, store_wl, pairs)
    pairs = pairs(:, 1:1)
    allocate(tm(6, 1), tn(6, 1))
    call direct_setup(mesh, DEFAULT_EDGE_DEG, nmask)
    mask = direct_edge_mask()
    call direct_at_omega(mesh, geom, om, store, store_wl, pairs, tm)
    call direct_setup(mesh, 180d0, nnone)
    call direct_at_omega(mesh, geom, om, store, store_wl, pairs, tn)
    ! The B2 linear parts on the masked panels alone, through the same solve and load
    allocate(a(mesh%n, 1), c(3, mesh%n, 1), e(geom%nwl, 1), clm(3, 4, mesh%n, 1), cln(3, 4, mesh%n, 1), &
             rm(mesh%nbase, 1, 1), rn(mesh%nbase, 1, 1), dphi(mesh%n, 1))
    call b2_stokes_densities(mesh, geom, store(1, 1), store(2, 1), store_wl(:, 1, 1), store_wl(:, 2, 1), a(:, 1), &
                             c(:, :, 1), e(:, 1), clm(:, :, :, 1), mask)
    call b2_stokes_densities(mesh, geom, store(1, 1), store(2, 1), store_wl(:, 1, 1), store_wl(:, 2, 1), a(:, 1), &
                             c(:, :, 1), e(:, 1), cln(:, :, :, 1))
    call direct_rhs_b2_block(mesh, geom, a, c, e, rm, cl=clm)
    call direct_rhs_b2_block(mesh, geom, a, c, e, rn, cl=cln)
    rm = rm - rn
    call direct_solve(mesh, rm, dphi)
    delta = (0d0, 0d0)
    do f = 1, mesh%n
      do j = 1, 6
        delta(j) = delta(j) + dphi(f, 1)*mesh%nj(j, f)*mesh%area(f)
      end do
    end do
    delta = (0d0, 1d0)*om*store(1, 1)%rho*delta
    worst = maxval(abs(tm(:, 1) - tn(:, 1) - delta))/maxval(abs(delta))
    write(line, '(a,i0,a,i0,a,es10.2,a,es10.2)') "masked ", nmask, " (none at 180 deg: ", nnone, &
      "); guard effect relative to T_dir ", maxval(abs(delta))/maxval(abs(tm)), ", against the B2 parts alone ", worst
    print '(4x,a)', trim(line)
    call check(error, nmask == 64 .and. count(mask) == 64 .and. nnone == 0 .and. &
               maxval(abs(delta)) > 1d-6*maxval(abs(tm)) .and. worst <= guard_limit, line)
    call direct_setup(mesh, DEFAULT_EDGE_DEG)
  end subroutine test_direct_edge_guard

  !> End to end: T_dir of direct_at_omega against T_pot of potential_pair (the indirect method, qtf_b2 = body, psi_k from
  !> the radiation solve at Omega) on the submerged sphere (512 panels, depth 50), for the pairs (1, 2) and (2, 1). The
  !> two are equal by Green's theorem, so they differ by the mesh error only; a sign or slot error would be of order 1.
  subroutine test_direct_vs_indirect(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    type(first_order_state), allocatable :: store(:, :)
    complex(8), allocatable :: store_wl(:, :, :), t(:, :), psi(:, :), gpsi(:, :, :), pwl(:, :)
    integer, allocatable :: pairs(:, :)
    complex(8) :: x1(6, 1), x2(6, 1), terms(6, NPOT_TERMS), tp(6, 2)
    real(8) :: om, rel, seconds, seconds_rad
    integer :: p
    integer(8) :: c0, c1, rate
    character(len=250) :: line

    call set_wave(50d0)
    NTHREAD = 4
    ISOLV = 1
    call make_sphere(5d0, [0d0, 0d0, -20d0], sph_nlat(2), full)
    call load_panels(full, 0, 0, mesh)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    om = w1 - w2
    call system_clock(c0, rate)
    call radiation_at(om, .false., mesh%nbase)
    call system_clock(c1)
    seconds_rad = real(c1 - c0, 8)/real(rate, 8)
    allocate(psi(mesh%n, 6), gpsi(3, mesh%n, 6), pwl(geom%nwl, 6), t(6, 2))
    call radiation_psi(mesh, om, MXPOT, psi)
    call psi_surface_gradient(mesh, geom, psi, gpsi)
    call psi_waterline(mesh, geom, psi, gpsi, pwl)
    call motions(x1, x2)
    call two_state_store(mesh, geom, 50d0, x1, x2, store, store_wl, pairs)
    call direct_setup(mesh, DEFAULT_EDGE_DEG)
    call system_clock(c0)
    call direct_at_omega(mesh, geom, om, store, store_wl, pairs, t)
    call system_clock(c1)
    seconds = real(c1 - c0, 8)/real(rate, 8)
    do p = 1, 2
      associate (s1 => store(pairs(1, p), 1), s2 => store(pairs(3, p), 1))
        call potential_pair(mesh, geom, s1, s2, store_wl(:, pairs(1, p), 1), store_wl(:, pairs(3, p), 1), psi, gpsi, &
                            pwl, om, s1%rho, wave_gravity(s1%omega, s1%k, s1%h), terms)
      end associate
      tp(:, p) = sum(terms, dim=2)
    end do
    rel = maxval(abs(t - tp))/maxval(abs(tp))
    write(line, '(a,es10.2,a,2es10.2,a,f7.3,a,f7.3,a)') "T_dir against T_pot, relative to max |T_pot| ", rel, &
      " (|T_pot| of the pairs", maxval(abs(tp), dim=1), "); direct_at_omega ", seconds, " s, radiation solve ", &
      seconds_rad, " s"
    print '(4x,a)', trim(line)
    call check(error, rel <= indirect_limit, line)
  end subroutine test_direct_vs_indirect

  ! ------------------------------------------------------------------------------------- irregular frequency removal

  !> Lid points of irregular frequency removal on the waterplane of a cylinder of make_cylinder (radius, ntheta, nr), centred
  !> at (x0, 0): the centroids of its bottom panels at z = 0, on the side the symmetry plane keeps (isy_in = 1: y > 0),
  !> appended to pts (3, n). The IRR routines of HAMS take the lid as field points only (iXYZ_P; iDXYZ_P is read, not used).
  subroutine append_lid(radius, ntheta, nr, x0, isy_in, pts)
    real(8), intent(in) :: radius, x0
    integer, intent(in) :: ntheta, nr, isy_in
    real(8), allocatable, intent(inout) :: pts(:, :)
    type(panel_mesh) :: c
    real(8), allocatable :: tmp(:, :)
    integer :: p, n
    if (.not. allocated(pts)) allocate(pts(3, 0))
    call make_cylinder(radius, 1d0, ntheta, 1, nr, c, .true.)
    do p = 1, c%n
      if (abs(c%nrm(3, p)) < 0.99d0) cycle
      if (isy_in == 1 .and. c%xc(2, p) <= 0d0) cycle
      n = size(pts, 2)
      allocate(tmp(3, n + 1))
      tmp(:, 1:n) = pts
      tmp(:, n + 1) = [c%xc(1, p) + x0, c%xc(2, p), 0d0]
      call move_alloc(tmp, pts)
    end do
  end subroutine append_lid

  !> IRSP = 1 with the lid points pts (3, nl): in the single-body arrays (iNELEM, iXYZ_P) or the combined multi-body ones
  subroutine set_lid(pts, multi)
    real(8), intent(in) :: pts(:, :)
    logical, intent(in) :: multi
    integer :: nl
    nl = size(pts, 2)
    IRSP = 1
    if (multi) then
      if (.not. allocated(iNELEM_TOTAL)) allocate(iNELEM_TOTAL)
      iNELEM_TOTAL = nl
      if (allocated(iXYZ_GLOBAL_MULTI_COMB_P)) deallocate(iXYZ_GLOBAL_MULTI_COMB_P)
      if (allocated(iDXYZ_MULTI_COMB)) deallocate(iDXYZ_MULTI_COMB)
      allocate(iXYZ_GLOBAL_MULTI_COMB_P(nl, 3), iDXYZ_MULTI_COMB(nl, 6))
      iXYZ_GLOBAL_MULTI_COMB_P = transpose(pts)
      iDXYZ_MULTI_COMB = 0d0
    else
      iNELEM = nl
      if (allocated(iXYZ_P)) deallocate(iXYZ_P)
      if (allocated(iDXYZ_P)) deallocate(iDXYZ_P)
      allocate(iXYZ_P(nl, 3), iDXYZ_P(nl, 6))
      iXYZ_P = transpose(pts)
      iDXYZ_P = 0d0
    end if
  end subroutine set_lid

  !> Back to IRSP = 0 after a test with a lid
  subroutine clear_lid()
    IRSP = 0
    iNELEM = 0
    if (allocated(iNELEM_TOTAL)) iNELEM_TOTAL = 0
  end subroutine clear_lid

  !> As radiation_at with IRSP = 1 and the nl lid points of set_lid, as SolveRadiation and SolveRadiationMulti run it:
  !> CALGREEN_IRR, ASSB_LEFT_IRR (AMAT with the lid rows, CMAT = AMAT^T AMAT, factored when ISOLV = 1), ASSB_RBC_IRR and
  !> RADIATION_SOLVER_IRR, or their multi-body forms (RADIATION_SOLVER_MULTI on CMAT_MULTI, as HAMS_Prog).
  subroutine radiation_at_irr(omega, multi, ne, nl)
    real(8), intent(in) :: omega
    logical, intent(in) :: multi
    integer, intent(in) :: ne, nl
    call SetWaveFromOmega(omega)
    TNELEM = ne + nl
    if (.not. multi) then
      if (allocated(CGRN)) deallocate(CGRN)
      if (allocated(RKBN)) deallocate(RKBN)
      if (allocated(DGRN)) deallocate(DGRN)
      if (allocated(PKBN)) deallocate(PKBN)
      if (allocated(AMAT)) deallocate(AMAT)
      if (allocated(BRMAT)) deallocate(BRMAT)
      if (allocated(CMAT)) deallocate(CMAT)
      if (allocated(DRMAT)) deallocate(DRMAT)
      if (allocated(IPIV)) deallocate(IPIV)
      allocate(CGRN(ne, ne, NSYS, 4), RKBN(ne, ne, NSYS, 4), DGRN(nl, ne, NSYS, 4), PKBN(nl, ne, NSYS, 4), &
               AMAT(TNELEM, TNELEM, NSYS), BRMAT(TNELEM, 6, NSYS), CMAT(ne, ne, NSYS), DRMAT(ne, 6, NSYS), IPIV(ne, NSYS))
      call CALGREEN_IRR
      call ASSB_LEFT_IRR(AMAT, CMAT, IPIV, ne, TNELEM, NSYS)
      call ASSB_RBC_IRR(BRMAT, DRMAT, AMAT, ne, TNELEM, NSYS)
      call RADIATION_SOLVER_IRR(CMAT, DRMAT, IPIV, MXPOT, ne, NSYS)
    else
      if (.not. allocated(NELEM_TOTAL)) allocate(NELEM_TOTAL)
      NELEM_TOTAL = ne
      if (allocated(CGRN_MULTI_COMB)) deallocate(CGRN_MULTI_COMB)
      if (allocated(RKBN_MULTI_COMB)) deallocate(RKBN_MULTI_COMB)
      if (allocated(DGRN_MULTI_COMB)) deallocate(DGRN_MULTI_COMB)
      if (allocated(PKBN_MULTI_COMB)) deallocate(PKBN_MULTI_COMB)
      if (allocated(AMAT_MULTI)) deallocate(AMAT_MULTI)
      if (allocated(BRMAT_MULTI)) deallocate(BRMAT_MULTI)
      if (allocated(CMAT_MULTI)) deallocate(CMAT_MULTI)
      if (allocated(DRMAT_MULTI)) deallocate(DRMAT_MULTI)
      if (allocated(IPIV_MULTI_COMB)) deallocate(IPIV_MULTI_COMB)
      allocate(CGRN_MULTI_COMB(ne, ne, NSYS, 4), RKBN_MULTI_COMB(ne, ne, NSYS, 4), DGRN_MULTI_COMB(nl, ne, NSYS, 4), &
               PKBN_MULTI_COMB(nl, ne, NSYS, 4), AMAT_MULTI(TNELEM, TNELEM, NSYS), BRMAT_MULTI(TNELEM, 12, NSYS), &
               CMAT_MULTI(ne, ne, NSYS), DRMAT_MULTI(ne, 12, NSYS), IPIV_MULTI_COMB(ne, NSYS))
      call CALGREEN_MULTI_IRR
      call ASSB_LEFT_IRR_MULTI(AMAT_MULTI, CMAT_MULTI, IPIV_MULTI_COMB, ne, nl, TNELEM, NSYS)
      BRMAT_MULTI = (0d0, 0d0)
      DRMAT_MULTI = (0d0, 0d0)
      call ASSB_RBC_IRR_MULTI(BRMAT_MULTI(:, 1:6, :), DRMAT_MULTI(:, 1:6, :), AMAT_MULTI, ne, nl, TNELEM, NSYS, 1, ne/2)
      call ASSB_RBC_IRR_MULTI(BRMAT_MULTI(:, 7:12, :), DRMAT_MULTI(:, 7:12, :), AMAT_MULTI, ne, nl, TNELEM, NSYS, &
                              ne/2 + 1, ne)
      call RADIATION_SOLVER_MULTI(CMAT_MULTI, DRMAT_MULTI, IPIV_MULTI_COMB, MXPOT_MULTI_COMB, ne, NSYS)
    end if
  end subroutine radiation_at_irr

  !> qtf_direct = 1 with irregular frequency removal (IRSP = 1) passes the input check (Task 3 refused it), and so does
  !> qtf_direct = 1 without it
  subroutine test_direct_irr_accepted(error)
    type(error_type), allocatable, intent(out) :: error
    type(so_options) :: opts
    logical :: ok
    character(len=200) :: msg
    integer :: irsp_in

    opts%active = .true.; opts%mean_drift = 1; opts%do_far = .true.; opts%do_near = .true.
    opts%qtf = 1; opts%qtf_terms = QTF_TERMS_BODY; opts%qtf_direct = 1
    do irsp_in = 0, 1
      call validate_second_order(opts, 1, -1d0, 1, irsp_in, .false., ok, msg)
      call check(error, ok, "qtf_direct = 1 refused with IRSP = "//achar(48 + irsp_in)//": "//trim(msg))
      if (allocated(error)) return
    end do
  end subroutine test_direct_irr_accepted

  !> With IRSP = 1 and b = n_k, the direct solve with the lid rows (ASSB_RHS_FIELD_IRR, the least-squares combination
  !> AMAT^T r and the factored CMAT) is psi_k of the radiation solve with the lid at Omega (RADIATION_SOLVER_IRR): the same
  !> matrix and the same right side, so to rounding. The cylinder of direct_unit_normal_is_radiation with its bottom
  !> panels lifted to z = 0 as the lid: whole at depth 10, half (ISY = 1) in deep water, whole with ISOLV = 2 (CMAT
  !> unfactored by ASSB_LEFT_IRR, factored by direct_factor), and two cylinders with their lids in the multi-body arrays,
  !> whole and as half meshes. In case 1 psi_k of the solve without the lid is printed too: the lid rows change psi_k.
  subroutine test_direct_irr_unit_normal(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    complex(8), allocatable :: vf(:, :, :), rf(:, :, :), phi(:, :), psi(:, :)
    real(8), allocatable :: pts(:, :)
    real(8) :: worst, w, d, lid_effect
    real(8), parameter :: om = 0.4d0
    integer :: icase, nb, nl, isy_in, b
    logical :: multi
    character(len=200) :: line

    call make_cylinder(5d0, 5d0, 16, 4, 3, full, .true.)
    worst = 0d0
    lid_effect = 0d0
    do icase = 1, 5
      d = merge(-1d0, 10d0, icase == 2)
      multi = icase >= 4
      isy_in = merge(1, 0, icase == 2 .or. icase == 5)
      call set_wave(d)
      NTHREAD = 4
      ISOLV = 1
      if (allocated(pts)) deallocate(pts)
      if (multi) then
        call load_multi(isy_in, mesh)
        do b = 1, 2
          call append_lid(multi_radius, multi_ntheta, multi_nr, multi_spacing*(b - 1), isy_in, pts)
        end do
      else
        call load_panels(full, 0, isy_in, mesh)
        call append_lid(5d0, 16, 3, 0d0, isy_in, pts)
      end if
      call set_lid(pts, multi)
      nb = mesh%nbase
      nl = size(pts, 2)
      call radiation_at_irr(om, multi, nb, nl)
      allocate(psi(mesh%n, 6*mesh%nbody), phi(mesh%n, 6*mesh%nbody), rf(nb + nl, 6*mesh%nbody, mesh%nsys))
      if (multi) then
        call radiation_psi(mesh, om, MXPOT_MULTI_COMB, psi)
      else
        call radiation_psi(mesh, om, MXPOT, psi)
      end if
      if (icase == 3) then
        ISOLV = 2
        call ASSB_LEFT_IRR(AMAT, CMAT, IPIV, nb, TNELEM, NSYS)
        call direct_factor(mesh)
      end if
      call unit_normal_fields(mesh, vf)
      if (multi) then
        call ASSB_RHS_FIELD_MULTI_IRR(vf, rf, nb, nl, mesh%nsys, 6*mesh%nbody)
      else
        call ASSB_RHS_FIELD_IRR(vf, rf, nb, nl, mesh%nsys, 6*mesh%nbody)
      end if
      call check(error, direct_lid_points(multi) == nl, "direct_lid_points is not the number of lid points")
      if (allocated(error)) then
        call clear_lid()
        return
      end if
      call direct_solve(mesh, rf, phi, multi)
      w = maxval(abs(phi - psi))/maxval(abs(psi))
      write(line, '(a,i0,a,i0,a,i0,a,i0,a,es10.2)') "case ", icase, ": ", mesh%n, " panels, ", nl, " lid points, nsys ", &
        mesh%nsys, ", |phi - psi_k| relative ", w
      print '(4x,a)', trim(line)
      worst = max(worst, w)
      if (icase == 1) then
        call clear_lid()
        call radiation_at(om, .false., nb)
        call radiation_psi(mesh, om, MXPOT, phi)
        lid_effect = maxval(abs(phi - psi))/maxval(abs(psi))
      end if
      deallocate(psi, phi, rf, vf)
    end do
    ISOLV = 1
    call clear_lid()
    ISX = 0; ISY = 0; NSYS = 1
    write(line, '(a,es10.2,a,es10.2)') "IRSP = 1, b = n_k against psi_k, worst relative difference ", worst, &
      "; psi_k with and without the lid differ by ", lid_effect
    print '(4x,a)', trim(line)
    call check(error, worst <= unit_normal_limit .and. lid_effect > 1d3*unit_normal_limit, line)
  end subroutine test_direct_irr_unit_normal

  !> The lid rows of the B2 right side (direct_rhs_b2_block with nbase + nl rows, IRSP = 1), Phase 7a Task 5: from the
  !> stored lid integrals (DGRN and PKBN of CALGREEN_IRR, CALGREEN_MULTI_IRR) against fresh ones (green_panel_integrals at
  !> the lid point), and the half mesh (ISY = 1, its lid points y > 0, the rows expanded by RXY) against the whole mesh, both
  !> fresh, at the full lid points. One cylinder (the one of the half test, 48 lid points whole) and two cylinders in the
  !> multi-body arrays (2 x 24 lid points whole), depth 10, the pair of test_b2_stored with the linear densities and the
  !> sharp-edge mask of the production path. The hull rows are compared stored against fresh too.
  subroutine test_b2_lid_rows(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    type(first_order_state) :: s1, s2
    complex(8) :: x1(6, 1), x2(6, 1), xm1(6, 2), xm2(6, 2)
    complex(8), allocatable :: wl1(:), wl2(:), a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :), rs(:, :, :), rf(:, :, :)
    complex(8), allocatable :: whole_lid(:)
    real(8), allocatable :: pts(:, :), whole_pts(:, :)
    logical, allocatable :: mask(:)
    complex(8) :: v
    real(8) :: stored_hull, stored_lid, half_lid, lid_size, xl(3)
    real(8), parameter :: d = 10d0
    integer :: icase, isy_in, nb, nl, l, im, ip, q
    logical :: multi
    character(len=250) :: line

    call make_cylinder(5d0, 5d0, 16, 4, 3, full, .true.)
    call motions(x1, x2)
    xm1(:, 1) = x1(:, 1)
    xm1(:, 2) = 0.5d0*x2(:, 1)
    xm2(:, 1) = x2(:, 1)
    xm2(:, 2) = (0d0, 1d0)*x1(:, 1)
    stored_hull = 0d0
    stored_lid = 0d0
    half_lid = 0d0
    lid_size = huge(1d0)
    do icase = 1, 4
      ! 1, 2: one cylinder whole and half; 3, 4: two cylinders whole and half
      multi = icase >= 3
      isy_in = merge(1, 0, mod(icase, 2) == 0)
      call set_wave(d)
      NTHREAD = 4
      ISOLV = 1
      if (allocated(pts)) deallocate(pts)
      if (multi) then
        call load_multi(isy_in, mesh)
        call append_lid(multi_radius, multi_ntheta, multi_nr, 0d0, isy_in, pts)
        call append_lid(multi_radius, multi_ntheta, multi_nr, multi_spacing, isy_in, pts)
      else
        call load_panels(full, 0, isy_in, mesh)
        call append_lid(5d0, 16, 3, 0d0, isy_in, pts)
      end if
      call geometry(mesh, geom, error)
      if (allocated(error)) exit
      call set_lid(pts, multi)
      nb = mesh%nbase
      nl = size(pts, 2)
      if (multi) then
        call make_state(mesh, geom, w1, b1, d, xm1, s1, wl1)
        call make_state(mesh, geom, w2, b2, d, xm2, s2, wl2)
      else
        call make_state(mesh, geom, w1, b1, d, x1, s1, wl1)
        call make_state(mesh, geom, w2, b2, d, x2, s2, wl2)
      end if
      call radiation_at_irr(w1 - w2, multi, nb, nl)
      call sharp_edge_panels(mesh, DEFAULT_EDGE_DEG, mask)
      allocate(a(mesh%n, 1), c(3, mesh%n, 1), e(geom%nwl, 1), cl(3, 4, mesh%n, 1), rs(nb + nl, 1, mesh%nsys), &
               rf(nb + nl, 1, mesh%nsys))
      call b2_stokes_densities(mesh, geom, s1, s2, wl1, wl2, a(:, 1), c(:, :, 1), e(:, 1), cl(:, :, :, 1), mask)
      call direct_rhs_b2_block(mesh, geom, a, c, e, rs, multi, cl)
      call direct_rhs_b2_block(mesh, geom, a, c, e, rf, multi, cl, stored=.false.)
      stored_hull = max(stored_hull, maxval(abs(rs(1:nb, :, :) - rf(1:nb, :, :)))/maxval(abs(rf(1:nb, :, :))))
      stored_lid = max(stored_lid, maxval(abs(rs(nb + 1:, :, :) - rf(nb + 1:, :, :)))/maxval(abs(rf(nb + 1:, :, :))))
      lid_size = min(lid_size, maxval(abs(rf(nb + 1:, :, :)))/maxval(abs(rf(1:nb, :, :))))
      if (isy_in == 0) then
        ! The whole mesh: its lid rows are the full lid points
        whole_pts = pts
        whole_lid = rf(nb + 1:, 1, 1)
      else
        ! The half mesh: full lid point (image im of lid point l) = sum_ip RXY(im, ip) rf(nb + l, ip) / nsys, at the lid
        ! point mirrored in y for im = 2 (as ASSB_RBC_IRR and RADIATION_SOLVER_IRR expand the blocks)
        do im = 1, mesh%nsys
          do l = 1, nl
            v = (0d0, 0d0)
            do ip = 1, mesh%nsys
              v = v + RXY(im, ip)*rf(nb + l, 1, ip)
            end do
            v = v/mesh%nsys
            xl = pts(:, l)
            if (im == 2) xl(2) = -xl(2)
            q = minloc(norm2(whole_pts - spread(xl, 2, size(whole_pts, 2)), dim=1), 1)
            if (norm2(whole_pts(:, q) - xl) > 1d-9) then
              call check(error, .false., "a lid point of the half mesh is not a lid point of the whole mesh")
              exit
            end if
            half_lid = max(half_lid, abs(v - whole_lid(q))/maxval(abs(whole_lid)))
          end do
          if (allocated(error)) exit
        end do
      end if
      write(line, '(a,i0,a,i0,a,i0,a,i0,a,es10.2,a,es10.2)') "case ", icase, ": ", mesh%n, " panels, ", nl, &
        " lid points, nsys ", mesh%nsys, "; stored against fresh, hull rows ", &
        maxval(abs(rs(1:nb, :, :) - rf(1:nb, :, :)))/maxval(abs(rf(1:nb, :, :))), ", lid rows ", &
        maxval(abs(rs(nb + 1:, :, :) - rf(nb + 1:, :, :)))/maxval(abs(rf(nb + 1:, :, :)))
      print '(4x,a)', trim(line)
      deallocate(a, c, e, cl, rs, rf)
      if (allocated(error)) exit
    end do
    call clear_lid()
    ISX = 0; ISY = 0; NSYS = 1
    if (allocated(error)) return
    write(line, '(a,es10.2,a,es10.2,a,es10.2,a,es10.2)') "B2 lid rows: stored against fresh ", stored_lid, &
      " (hull rows ", stored_hull, "); half against whole ", half_lid, "; smallest lid rows over hull rows ", lid_size
    print '(4x,a)', trim(line)
    call check(error, stored_lid <= stored_limit .and. stored_hull <= stored_limit .and. half_lid <= lid_half_limit .and. &
               lid_size > 1d-3, line)
  end subroutine test_b2_lid_rows

  !> Direct against indirect with the lid, through the irregular frequency, on a cylinder with a rounded bilge (radius 5,
  !> draft 5, bilge radius 1.5, 448 panels, depth 50; 192 lid points). Its first irregular frequency is near that of the
  !> sharp cylinder, J0(k a) = 0 for the interior problem with nu = k coth(k d): 2.19 rad/s; the discrete one of this mesh
  !> is irr_omegas(3) (measured, see there). At a low Omega, at a regular one next to the irregular frequency and at it:
  !>   with the lid, T_dir of direct_at_omega equals T_pot of potential_pair with psi of the lid solve to the mesh error,
  !>   and at the irregular frequency that error is the one of its regular neighbour;
  !>   without the lid, T_dir agrees with the lid result at the low Omega. T_pot without the lid is printed as well.
  !> The spike that the lid removes is checked over a window of Omega about the irregular frequency, not at one Omega
  !> inside a peak 0.002 rad/s wide (Phase 7a Task 5): with the lid, T_dir at every Omega of irr_window is close to the
  !> straight line through its values at irr_ends, 0.05 rad/s either side; without the lid, T_dir is far from the lid
  !> result at some Omega of the window. A shift of the discrete irregular frequency within the window changes neither.
  subroutine test_direct_vs_indirect_irr(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    type(first_order_state), allocatable :: store(:, :)
    complex(8), allocatable :: store_wl(:, :, :)
    real(8), allocatable :: pts(:, :)
    integer, allocatable :: pairs(:, :)
    complex(8) :: x1(6, 1), x2(6, 1), t(6, 2, 2), tp(6, 2, 2), tend(6, 2, 2), tw(6, 2, 2), line_t(6, 2)
    real(8) :: om, rel(3, 2), spike(3), spike_pot(3), smooth, spike_window, f, smooth_nolid
    integer :: io, lid, nb, nl, k
    character(len=300) :: line

    call irr_hull(full, mesh, geom, pts, error)
    if (allocated(error)) return
    nb = mesh%nbase
    nl = size(pts, 2)
    call motions(x1, x2)
    call direct_setup(mesh, DEFAULT_EDGE_DEG)
    do io = 1, 3
      om = irr_omegas(io)
      call two_state_store(mesh, geom, irr_depth, x1, x2, store, store_wl, pairs, irr_low + om, irr_low)
      do lid = 1, 2
        ! lid = 1 without the lid, 2 with it
        if (lid == 1) then
          call clear_lid()
          call radiation_at(om, .false., nb)
        else
          call set_lid(pts, .false.)
          call radiation_at_irr(om, .false., nb, nl)
        end if
        call both_methods(mesh, geom, om, store, store_wl, pairs, t(:, :, lid), tp(:, :, lid))
        rel(io, lid) = maxval(abs(t(:, :, lid) - tp(:, :, lid)))/maxval(abs(tp(:, :, lid)))
      end do
      spike(io) = maxval(abs(t(:, :, 1) - t(:, :, 2)))/maxval(abs(t(:, :, 2)))
      spike_pot(io) = maxval(abs(tp(:, :, 1) - tp(:, :, 2)))/maxval(abs(tp(:, :, 2)))
      write(line, '(a,f7.4,a,es10.2,a,es10.2,a,es10.2,a,es10.2,a,es10.2)') "Omega ", om, &
        ": T_dir against T_pot with the lid ", rel(io, 2), ", without ", rel(io, 1), &
        "; without the lid against with: T_dir ", spike(io), ", T_pot ", spike_pot(io), "; max |T_dir| ", &
        maxval(abs(t(:, :, 2)))
      print '(4x,a)', trim(line)
    end do
    ! The window: T_dir with the lid at its two ends, then at each Omega of it with and without the lid
    do k = 1, 2
      call direct_only(irr_ends(k), .true., tend(:, :, k))
    end do
    smooth = 0d0
    smooth_nolid = 0d0
    spike_window = 0d0
    do k = 1, size(irr_window)
      om = irr_window(k)
      call direct_only(om, .false., tw(:, :, 1))
      call direct_only(om, .true., tw(:, :, 2))
      f = (om - irr_ends(1))/(irr_ends(2) - irr_ends(1))
      line_t = (1d0 - f)*tend(:, :, 1) + f*tend(:, :, 2)
      smooth = max(smooth, maxval(abs(tw(:, :, 2) - line_t))/maxval(abs(tw(:, :, 2))))
      smooth_nolid = max(smooth_nolid, maxval(abs(tw(:, :, 1) - line_t))/maxval(abs(tw(:, :, 2))))
      spike_window = max(spike_window, maxval(abs(tw(:, :, 1) - tw(:, :, 2)))/maxval(abs(tw(:, :, 2))))
      write(line, '(a,f7.4,a,es10.2,a,es10.2,a,es10.2)') "window Omega ", om, ": with the lid from the line ", &
        maxval(abs(tw(:, :, 2) - line_t))/maxval(abs(tw(:, :, 2))), ", without the lid from the line ", &
        maxval(abs(tw(:, :, 1) - line_t))/maxval(abs(tw(:, :, 2))), ", without against with ", &
        maxval(abs(tw(:, :, 1) - tw(:, :, 2)))/maxval(abs(tw(:, :, 2)))
      print '(4x,a)', trim(line)
    end do
    call clear_lid()
    write(line, '(a,3es10.2,a,f6.3,a,es10.2)') "T_dir against T_pot with the lid ", rel(:, 2), &
      " (irregular over its neighbour ", rel(3, 2)/rel(2, 2), "); T_dir without the lid against with, low ", spike(1)
    print '(4x,a)', trim(line)
    write(line, '(a,es10.2,a,es10.2,a,es10.2)') "window: with the lid at most ", smooth, &
      " from the line through the ends (without the lid ", smooth_nolid, "); without the lid against with, largest ", &
      spike_window
    print '(4x,a)', trim(line)
    call check(error, rel(1, 2) <= irr_indirect_limit .and. maxval(rel(2:3, 2)) <= irr_mesh_limit .and. &
               rel(3, 2) <= irr_ratio_limit*rel(2, 2) .and. spike(1) <= irr_regular_limit .and. &
               smooth <= irr_smooth_limit .and. spike_window >= irr_spike_min, line)

  contains

    !> T_dir of the two pairs at Omega om, with or without the lid
    subroutine direct_only(om_in, with_lid, td)
      real(8), intent(in) :: om_in
      logical, intent(in) :: with_lid
      complex(8), intent(out) :: td(:, :)
      call two_state_store(mesh, geom, irr_depth, x1, x2, store, store_wl, pairs, irr_low + om_in, irr_low)
      if (with_lid) then
        call set_lid(pts, .false.)
        call radiation_at_irr(om_in, .false., nb, nl)
      else
        call clear_lid()
        call radiation_at(om_in, .false., nb)
      end if
      call direct_at_omega(mesh, geom, om_in, store, store_wl, pairs, td)
    end subroutine direct_only

  end subroutine test_direct_vs_indirect_irr

  !> The hull, geometry and lid points of direct_vs_indirect_irr, with the solver globals at depth irr_depth
  subroutine irr_hull(full, mesh, geom, pts, error)
    type(panel_mesh), intent(out) :: full, mesh
    type(surface_geom), intent(out) :: geom
    real(8), allocatable, intent(out) :: pts(:, :)
    type(error_type), allocatable, intent(out) :: error
    call make_rounded_cylinder(5d0, 5d0, 1.5d0, 32, 6, 4, 4, full)
    call set_wave(irr_depth)
    NTHREAD = 4
    ISOLV = 1
    call load_panels(full, 0, 0, mesh)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    call append_lid(5d0, 32, 6, 0d0, 0, pts)
  end subroutine irr_hull

  !> T_dir (direct_at_omega) and T_pot (potential_pair, the sum of its terms, psi of the radiation solve just run at om)
  !> of the two pairs of two_state_store
  subroutine both_methods(mesh, geom, om, store, store_wl, pairs, t, tp)
    type(panel_mesh), intent(in) :: mesh
    type(surface_geom), intent(in) :: geom
    real(8), intent(in) :: om
    type(first_order_state), intent(in) :: store(:, :)
    complex(8), intent(in) :: store_wl(:, :, :)
    integer, intent(in) :: pairs(:, :)
    complex(8), intent(out) :: t(:, :), tp(:, :)
    complex(8), allocatable :: psi(:, :), gpsi(:, :, :), pwl(:, :)
    complex(8) :: terms(6, NPOT_TERMS)
    integer :: p
    allocate(psi(mesh%n, 6), gpsi(3, mesh%n, 6), pwl(geom%nwl, 6))
    call direct_at_omega(mesh, geom, om, store, store_wl, pairs, t)
    call radiation_psi(mesh, om, MXPOT, psi)
    call psi_surface_gradient(mesh, geom, psi, gpsi)
    call psi_waterline(mesh, geom, psi, gpsi, pwl)
    do p = 1, 2
      associate (s1 => store(pairs(1, p), 1), s2 => store(pairs(3, p), 1))
        call potential_pair(mesh, geom, s1, s2, store_wl(:, pairs(1, p), 1), store_wl(:, pairs(3, p), 1), psi, gpsi, &
                            pwl, om, s1%rho, wave_gravity(s1%omega, s1%k, s1%h), terms)
      end associate
      tp(:, p) = sum(terms, dim=2)
    end do
  end subroutine both_methods


  !> A vertical cylinder with a rounded bilge (no edge), about the z axis, waterline at z = 0, normals into the body, as
  !> make_cylinder: the profile from (radius, 0) down the side (nz rows) to z = -(draft - bilge), a quarter circle of
  !> radius bilge (narc rows) to (radius - bilge, -draft), and the flat bottom (nr rings, triangles at the axis), revolved
  !> with ntheta panels per ring.
  subroutine make_rounded_cylinder(radius, draft, bilge, ntheta, nz, narc, nr, mesh)
    real(8), intent(in) :: radius, draft, bilge
    integer, intent(in) :: ntheta, nz, narc, nr
    type(panel_mesh), intent(out) :: mesh
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), allocatable :: pr(:), pz(:)
    real(8) :: th, ang, v(3, 4), d1(3), d2(3), nv(3), inside(3)
    integer :: nk, k, it, itn, ip, nnode, npan, m, i

    nk = nz + narc + nr
    allocate(pr(0:nk), pz(0:nk))
    do k = 0, nz
      pr(k) = radius
      pz(k) = -(draft - bilge)*real(k, 8)/real(nz, 8)
    end do
    do k = 1, narc
      ang = 0.5d0*pi*real(k, 8)/real(narc, 8)
      pr(nz + k) = radius - bilge + bilge*cos(ang)
      pz(nz + k) = -(draft - bilge) - bilge*sin(ang)
    end do
    do k = 1, nr
      pr(nz + narc + k) = (radius - bilge)*real(nr - k, 8)/real(nr, 8)
      pz(nz + narc + k) = -draft
    end do
    nnode = ntheta*nk + 1
    npan = ntheta*nk
    mesh%n = npan; mesh%nbase = npan; mesh%nnode = nnode; mesh%nbody = 1; mesh%nsys = 1
    allocate(mesh%xc(3, npan), mesh%nrm(3, npan), mesh%area(npan), mesh%psize(npan), mesh%nj(6, npan))
    allocate(mesh%xnode(3, nnode), mesh%ncn(npan), mesh%ncon(4, npan), mesh%body(npan), mesh%base(npan), mesh%image(npan))
    allocate(mesh%xr(3, 1))
    mesh%xr = 0d0
    mesh%body = 1
    mesh%image = 1
    mesh%base = [(i, i = 1, npan)]
    mesh%ncon = 0
    do k = 0, nk - 1
      do it = 1, ntheta
        th = 2d0*pi*real(it - 1, 8)/real(ntheta, 8)
        mesh%xnode(:, k*ntheta + it) = [pr(k)*cos(th), pr(k)*sin(th), pz(k)]
      end do
    end do
    mesh%xnode(:, nnode) = [0d0, 0d0, -draft]
    inside = [0d0, 0d0, -0.5d0*draft]
    ip = 0
    do k = 1, nk
      do it = 1, ntheta
        itn = mod(it, ntheta) + 1
        ip = ip + 1
        if (k < nk) then
          mesh%ncn(ip) = 4
          mesh%ncon(1:4, ip) = [(k - 1)*ntheta + it, k*ntheta + it, k*ntheta + itn, (k - 1)*ntheta + itn]
        else
          mesh%ncn(ip) = 3
          mesh%ncon(1:3, ip) = [nnode, (k - 1)*ntheta + itn, (k - 1)*ntheta + it]
        end if
        do m = 1, mesh%ncn(ip)
          v(:, m) = mesh%xnode(:, mesh%ncon(m, ip))
        end do
        if (mesh%ncn(ip) == 4) then
          d1 = v(:, 3) - v(:, 1)
          d2 = v(:, 4) - v(:, 2)
        else
          d1 = v(:, 2) - v(:, 1)
          d2 = v(:, 3) - v(:, 1)
        end if
        nv = [d1(2)*d2(3) - d1(3)*d2(2), d1(3)*d2(1) - d1(1)*d2(3), d1(1)*d2(2) - d1(2)*d2(1)]
        nv = nv/norm2(nv)
        if (sum(nv*(inside - sum(v(:, 1:mesh%ncn(ip)), dim=2)/mesh%ncn(ip))) < 0d0) nv = -nv
        call set_panel(mesh, ip, nv)
      end do
    end do
  end subroutine make_rounded_cylinder

end module test_qtf_direct
