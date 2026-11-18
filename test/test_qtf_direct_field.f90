!> Unit tests of QTFDirectField (Phase 7c-1): phi_D of the direct solution at field points off the hull by Green's formula,
!> 4 pi phi = sum [ int G b - (int dG/dn) phi_D ], with the B2 part of b in the Stokes form.
!> These tests set module globals (the mesh arrays, the stored Green function, the wave), so the suite runs serially (see
!> unit_tester.f90).
module test_qtf_direct_field
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use HAMS_mod, only : NTHREAD, ISOLV
  use WaveDyn_mod, only : VWAVE => V, WK
  use ImplementSubs, only : SetWaveFromOmega
  use PanelMesh_mod, only : ISX, ISY, NSYS
  use Potentials_mod, only : MXPOT, MXPOT_MULTI_COMB
  use AssbMatx, only : ASSB_RHS_FIELD, ASSB_RHS_FIELD_IRR
  use SecondOrderTypes, only : panel_mesh, surface_geom, first_order_state, cs_points, cs_fields, CS_LINE
  use FirstOrderState, only : radiation_psi
  use FieldVelocity, only : fields_green
  use NearFieldDrift, only : pair_product, pair_cross, rigid_displacement, KERNEL_DIFF, KERNEL_SUM
  use QTFDirect, only : b2_stokes_densities, b2_stokes_at_point, direct_rhs_b2_block, direct_solve, rankine_moments, &
                        sharp_edge_panels, direct_setup, direct_edge_mask, direct_at_omega, &
                        direct_field_at_points
  use SharpEdge, only : DEFAULT_EDGE_DEG
  use INFG3D_Open, only : INFGREEN3D
  use second_order_fixtures, only : make_cylinder
  use test_waterline_field, only : set_wave, load_panels, load_multi
  use test_qtf_potential, only : motions, make_state, geometry
  use test_qtf_direct, only : radiation_at, append_lid, set_lid, clear_lid, radiation_at_irr
  use test_sharp_edge, only : make_closed_box
  implicit none
  private

  public :: collect_qtf_direct_field

  real(8), parameter :: pi = 3.141592653589793d0
  ! The difference frequency of the solves, and the two first-order states of the B2 densities (frequency, heading)
  real(8), parameter :: om = 0.4d0, w1 = 0.8d0, b1 = 0d0, w2 = 0.7d0, b2 = 0.3d0
  ! The truncated cylinder of test_qtf_direct: radius 5, draft 5, 112 panels, depth 10
  real(8), parameter :: cyl_radius = 5d0, cyl_draft = 5d0, cyl_depth = 10d0
  ! Against fields_green on the same rows (both in field mode): the same sums in another order; 1.8e-16 measured
  real(8), parameter :: radiation_limit = 1d-10
  ! Stokes form against the plain form, on fields for which both are integrated exactly (test_b2_identity); 9.8e-8
  ! measured
  real(8), parameter :: b2_limit = 1d-6
  ! The B2 part of direct_field_at_points against b2_stokes_at_point in field mode: the same sums, so to rounding; and the
  ! default route must differ from it by more than route_min_diff at these points, so that the test sees the route
  ! (ruling C7 of the 7c-1 ledger). Measured: the default route differs by 2.1e-7; route_min_diff is a tenth of it
  real(8), parameter :: route_limit = 1d-12, route_min_diff = 2d-8
  ! The field 0.25 panel sizes off a panel centroid against the panel value (spec, Phase 7 item 2); 1.5e-3 measured
  real(8), parameter :: approach_limit = 2d-2
  ! Half mesh against whole: the same sums in another order, and the Green function at mirrored points; 1.7e-14 measured
  real(8), parameter :: symmetry_limit = 1d-9
  ! The pair (2, 1) against the conjugate of the pair (1, 2): the same sums conjugated, so to rounding (Review Focus 1)
  real(8), parameter :: conjugate_limit = 1d-12
  ! Two cylinders in the multi-body arrays (load_multi): phi_field against fields_green, as radiation_limit
  ! (3.6e-16 measured; the limit is 2.8e3 times it, and with is_multi = .false. the difference is 0.75)
  real(8), parameter :: two_body_limit = 1d-12
  ! The hook with the lid (IRSP = 1) against without it, at Omega = w1 - w2, far from an irregular frequency: the lid rows
  ! change phi_D by the mesh-level amount, so the field differs by that order; it must differ (lid_min_diff: the lid solve
  ! is not the plain one) and at most lid_limit. Measured: 8.75e-3 in phi_field (1.05e-3 in T_dir, the mesh-level change of
  ! irr_regular_limit of test_qtf_direct); lid_limit is 2.3 times it, lid_min_diff 1e-4 times it
  real(8), parameter :: lid_limit = 2d-2, lid_min_diff = 1d-6

contains

  subroutine collect_qtf_direct_field(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("direct_field_radiation_identity", test_radiation_identity), &
      new_unittest("direct_field_b2_identity", test_b2_identity), &
      new_unittest("direct_field_b2_field_route", test_b2_field_route), &
      new_unittest("direct_field_approaches_panel", test_approaches_panel), &
      new_unittest("direct_field_symmetry", test_symmetry), &
      new_unittest("direct_field_hook_matches", test_hook_matches), &
      new_unittest("direct_field_sum_hook_matches", test_hook_sum_matches), &
      new_unittest("direct_field_negative_omega", test_negative_omega), &
      new_unittest("direct_field_hook_absent", test_hook_absent), &
      new_unittest("direct_field_two_bodies", test_two_bodies), &
      new_unittest("direct_field_lid", test_lid) &
      ]
  end subroutine collect_qtf_direct_field

  !> The panel forcing in the layout of pair_forcing for one column: vf(base(f), 1, image(f)) = field(f).
  subroutine panel_forcing(mesh, field, vf)
    type(panel_mesh), intent(in) :: mesh
    complex(8), intent(in) :: field(:)
    complex(8), allocatable, intent(out) :: vf(:, :, :)
    integer :: f
    allocate(vf(mesh%nbase, 1, mesh%nsys))
    vf = (0d0, 0d0)
    do f = 1, mesh%n
      vf(mesh%base(f), 1, mesh%image(f)) = field(f)
    end do
  end subroutine panel_forcing

  !> B2 densities of one column, all zero.
  subroutine zero_densities(mesh, geom, a, c, e, cl)
    type(panel_mesh), intent(in) :: mesh
    type(surface_geom), intent(in) :: geom
    complex(8), allocatable, intent(out) :: a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :)
    allocate(a(mesh%n, 1), c(3, mesh%n, 1), e(geom%nwl, 1), cl(3, 4, mesh%n, 1))
    a = (0d0, 0d0)
    c = (0d0, 0d0)
    e = (0d0, 0d0)
    cl = (0d0, 0d0)
  end subroutine zero_densities

  !> The B2 densities of one column for the pair of unit incident waves (w1, b1) and (w2, b2) with the rigid motions of
  !> test_qtf_potential, at the given depth, with the sharp-edge mask at the default angle, as direct_at_omega has them.
  subroutine pair_densities(mesh, geom, depth, a, c, e, cl)
    type(panel_mesh), intent(in) :: mesh
    type(surface_geom), intent(in) :: geom
    real(8), intent(in) :: depth
    complex(8), allocatable, intent(out) :: a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :)
    type(first_order_state) :: s1, s2
    complex(8) :: x1(6, 1), x2(6, 1)
    complex(8), allocatable :: wl1(:), wl2(:)
    logical, allocatable :: mask(:)
    call motions(x1, x2)
    call make_state(mesh, geom, w1, b1, depth, x1, s1, wl1)
    call make_state(mesh, geom, w2, b2, depth, x2, s2, wl2)
    call sharp_edge_panels(mesh, DEFAULT_EDGE_DEG, mask)
    allocate(a(mesh%n, 1), c(3, mesh%n, 1), e(geom%nwl, 1), cl(3, 4, mesh%n, 1))
    call b2_stokes_densities(mesh, geom, s1, s2, wl1, wl2, a(:, 1), c(:, :, 1), e(:, 1), cl(:, :, :, 1), mask)
  end subroutine pair_densities

  !> The index of a panel on the side wall of the cylinder (horizontal normal) in the row nearest z = z0.
  integer function side_panel(mesh, z0) result(ip)
    type(panel_mesh), intent(in) :: mesh
    real(8), intent(in) :: z0
    real(8) :: best
    integer :: i
    ip = 0
    best = huge(1d0)
    do i = 1, mesh%n
      if (abs(mesh%nrm(3, i)) > 0.1d0) cycle
      if (abs(mesh%xc(3, i) - z0) < best) then
        best = abs(mesh%xc(3, i) - z0)
        ip = i
      end if
    end do
  end function side_panel

  ! ------------------------------------------------------------------------------------------------ radiation identity

  !> With vf = n_3, no B2 and phi_D = psi_3 of the heave radiation solve at Omega, Green's formula of
  !> direct_field_at_points is the formula of fields_green for the state with phi_p = psi_3, dphi_p/dn = n_3 and no
  !> incident wave (amp = 0): the same rows, so to rounding. The truncated cylinder at depth 10 and in deep water, at 6
  !> points: 2 below z = 0 away from the body, 2 on z = 0, 2 half a panel size off the hull along -n (into the fluid).
  !> This pins the 4 pi and the sign of the dG/dn term against the first-order field.
  subroutine test_radiation_identity(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    type(first_order_state) :: states(1)
    type(cs_points) :: pts
    type(cs_fields) :: fields(1)
    complex(8), allocatable :: psi(:, :), vf(:, :, :), a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :), phi_field(:, :)
    real(8) :: xp(3, 6), worst, rel
    real(8), parameter :: depths(2) = [cyl_depth, -1d0]
    integer :: icase, ia, ib
    character(len=200) :: line

    call make_cylinder(cyl_radius, cyl_draft, 16, 4, 3, full, .true.)
    worst = 0d0
    do icase = 1, 2
      call set_wave(depths(icase))
      NTHREAD = 4
      ISOLV = 1
      call load_panels(full, 0, 0, mesh)
      call geometry(mesh, geom, error)
      if (allocated(error)) return
      call radiation_at(om, .false., mesh%nbase)
      allocate(psi(mesh%n, 6))
      call radiation_psi(mesh, om, MXPOT, psi)
      call panel_forcing(mesh, cmplx(mesh%nj(3, :), 0d0, kind=8), vf)
      call zero_densities(mesh, geom, a, c, e, cl)
      ia = side_panel(mesh, -0.5d0*cyl_draft)
      ib = minloc(mesh%xc(3, :), 1)                 ! a bottom panel
      xp(:, 1) = [12d0, 3d0, -2d0]
      xp(:, 2) = [-4d0, -9d0, -7d0]
      xp(:, 3) = [8d0, -2d0, 0d0]
      xp(:, 4) = [-3d0, 11d0, 0d0]
      xp(:, 5) = mesh%xc(:, ia) - 0.5d0*mesh%psize(ia)*mesh%nrm(:, ia)
      xp(:, 6) = mesh%xc(:, ib) - 0.5d0*mesh%psize(ib)*mesh%nrm(:, ib)
      allocate(phi_field(6, 1))
      call direct_field_at_points(mesh, geom, xp, psi(:, 3:3), vf, a, c, e, cl, .false., phi_field)

      states(1)%omega = om
      states(1)%k = WK
      states(1)%h = depths(icase)
      states(1)%nu = VWAVE
      states(1)%amp = 0d0
      states(1)%phi_p = psi(:, 3)
      states(1)%dphi_p_dn = cmplx(mesh%nj(3, :), 0d0, kind=8)
      pts%n = 6
      pts%x = xp
      allocate(pts%nrm(3, 6), pts%w(6), pts%psize(6), pts%kind(6), pts%body(6))
      pts%nrm = 0d0
      pts%w = 1d0
      pts%psize = 1d0
      pts%kind = CS_LINE        ! phi at the centre at every depth
      pts%body = 1
      call fields_green(mesh, pts, states, 0.01d0, .false., fields)
      rel = maxval(abs(phi_field(:, 1) - fields(1)%phi))/maxval(abs(fields(1)%phi))
      write(line, '(a,f6.1,a,es10.2,a,es10.2)') "depth ", depths(icase), ": phi_field against fields_green, relative ", &
        rel, "; largest |phi| ", maxval(abs(fields(1)%phi))
      print '(4x,a)', trim(line)
      worst = max(worst, rel)
      deallocate(psi, vf, a, c, e, cl, phi_field, pts%x, pts%nrm, pts%w, pts%psize, pts%kind, pts%body)
    end do
    write(line, '(a,es10.2)') "radiation identity, worst relative difference ", worst
    call check(error, worst <= radiation_limit, line)
  end subroutine test_radiation_identity

  ! ------------------------------------------------------------------------------------------------ B2 identity

  !> The B2 part of direct_field_at_points (vf = 0, phi_D = 0, so 4 pi phi = -r_B2) against the plain form
  !> -(1/4 pi) int_S G n.[(X.grad)u] dS by 8 x 8 sub-panels, at 3 points off a submerged closed box (8 x 6 x 4 m, top at
  !> z = -1, 1 m panels, no waterline). The fields are chosen so that both forms are integrated exactly on the flat panels,
  !> and the comparison is of the Stokes form and its implementation, not of the panel discretisation (which the 7a tests
  !> rhs_b2_stokes_vs_direct_* measure): state 1 has a uniform u1 and a translation X1, state 2 the linear harmonic
  !> u2 = A2 x (A2 symmetric, trace 0, so curl u2 = div u2 = 0) and a translation with a rotation X2. Under the kernel
  !> pair_minus the products are X1 x conj(u2) and conj(X2) x u1, both linear on a panel, so c = (X x u) x n and
  !> a = n.(alpha x u) are linear, which the fit over the coplanar neighbours of a face and the moments of the near
  !> panels integrate exactly; the plain density n.[(X.grad)u] = (1/4) X1.conj(A2 n) is constant per panel. The Green
  !> function is close to its Rankine part (deep water, 0.005 rad/s), whose panel and sub-panel integrals are in closed
  !> form; the rest, of order nu, is taken at one point per panel by the Stokes form (its gradient) and at the sub-panel
  !> centroids by the reference (its value). That rest sets the difference: 1.0e-5 at 0.05 rad/s (the frequency of
  !> rhs_b2_stokes_vs_direct_off_surface), 9.8e-8 at 0.005 rad/s, so nu / 100 gives 1 / 100 (measured); with the sign of
  !> B2 flipped the difference is 2.
  subroutine test_b2_identity(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    type(first_order_state) :: s1, s2
    complex(8), allocatable :: vf(:, :, :), a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :), phi_d(:, :), &
                               phi_field(:, :), wl(:)
    complex(8) :: u1(3), x1(6, 1), x2(6, 1), a2(3, 3), u2(3), plain(3)
    real(8) :: xp(3, 3), rel
    integer :: i, ipt
    character(len=250) :: line

    call set_wave(-1d0)
    VWAVE = 0.005d0**2/9.80665d0
    NTHREAD = 4
    call make_closed_box(8d0, 6d0, 4d0, 8, 6, 4, full)
    call load_panels(full, 0, 0, mesh)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    if (geom%nwl /= 0) error stop 'test_b2_identity: the submerged box has waterline segments'

    u1 = [(0.3d0, -0.2d0), (0.1d0, 0.4d0), (-0.5d0, 0.1d0)]
    a2 = reshape([(0.20d0, 0.10d0), (0.05d0, -0.10d0), (-0.07d0, 0.02d0), &
                  (0.05d0, -0.10d0), (-0.12d0, 0.03d0), (0.09d0, 0.04d0), &
                  (-0.07d0, 0.02d0), (0.09d0, 0.04d0), (-0.08d0, -0.13d0)], [3, 3])
    x1(:, 1) = [(0.3d0, 0.1d0), (-0.2d0, 0d0), (0.5d0, -0.3d0), (0d0, 0d0), (0d0, 0d0), (0d0, 0d0)]
    x2(:, 1) = [(0.1d0, 0.2d0), (0d0, 0.3d0), (-0.4d0, 0d0), (0d0, -0.02d0), (0.01d0, 0d0), (0.03d0, 0.01d0)]
    call box_state(w1, x1, s1, .true.)
    call box_state(w2, x2, s2, .false.)
    allocate(a(mesh%n, 1), c(3, mesh%n, 1), e(0, 1), cl(3, 4, mesh%n, 1), wl(0))
    call b2_stokes_densities(mesh, geom, s1, s2, wl, wl, a(:, 1), c(:, :, 1), e(:, 1), cl(:, :, :, 1))
    allocate(vf(mesh%nbase, 1, mesh%nsys), phi_d(mesh%n, 1), phi_field(3, 1))
    vf = (0d0, 0d0)
    phi_d = (0d0, 0d0)
    xp(:, 1) = [6.5d0, 1d0, -3d0]
    xp(:, 2) = [-1d0, 5d0, -2d0]
    xp(:, 3) = [1d0, -0.5d0, -8d0]
    call direct_field_at_points(mesh, geom, xp, phi_d, vf, a, c, e, cl, .false., phi_field)
    do ipt = 1, 3
      plain(ipt) = plain_form(xp(:, ipt))
    end do
    rel = maxval(abs(phi_field(:, 1) + plain/(4d0*pi)))/maxval(abs(plain/(4d0*pi)))
    write(line, '(a,es10.2,a,es10.2,a,es10.2)') "Stokes form against the plain form, relative ", rel, &
      "; with the sign of B2 flipped ", maxval(abs(phi_field(:, 1) - plain/(4d0*pi)))/maxval(abs(plain/(4d0*pi))), &
      "; largest |phi| ", maxval(abs(plain/(4d0*pi)))
    print '(4x,a)', trim(line)
    call check(error, rel <= b2_limit, line)

  contains

    !> A state with velocity u at the panel centroids (in grad_i) and rigid motion xi: u1 if uniform, else A2 x.
    subroutine box_state(omega, xi, s, uniform)
      real(8), intent(in) :: omega
      complex(8), intent(in) :: xi(6, 1)
      type(first_order_state), intent(out) :: s
      logical, intent(in) :: uniform
      s%omega = omega
      allocate(s%grad_i(3, mesh%n), s%grad_p(3, mesh%n), s%xi(6, 1))
      s%grad_p = (0d0, 0d0)
      s%xi = xi
      do i = 1, mesh%n
        if (uniform) then
          s%grad_i(:, i) = u1
        else
          s%grad_i(:, i) = matmul(a2, cmplx(mesh%xc(:, i), 0d0, kind=8))
        end if
      end do
    end subroutine box_state

    !> int_S G(x, .) n.[(X.grad)u] dS: on each panel the constant density sum_j pair_minus(X1_j, 0, X2_j, (A2 n)_j)
    !> times int G over 8 x 8 sub-panels, 1/r and 1/r1 in closed form (rankine_moments) and the rest at the sub-panel
    !> centroid.
    complex(8) function plain_form(x) result(out)
      real(8), intent(in) :: x(3)
      real(8) :: xv(3, 4), sv(3, 4), sc(3), xim(3), m0, m0i, g0(3), m1(3, 3), mg(3), d(3), r1, us(2), vs(2)
      complex(8) :: grn(4), dens, h2(3), d1(3), d2(3)
      integer :: p, j, iu, iv, k
      integer, parameter :: nsub = 8
      out = (0d0, 0d0)
      xim = [x(1), x(2), -x(3)]
      do p = 1, mesh%n
        d1 = rigid_displacement(x1(:, 1), mesh%xc(:, p))
        d2 = rigid_displacement(x2(:, 1), mesh%xc(:, p))
        h2 = matmul(cmplx(mesh%nrm(:, p), 0d0, kind=8), a2)
        dens = (0d0, 0d0)
        do j = 1, 3
          dens = dens + pair_product(KERNEL_DIFF, d1(j), (0d0, 0d0), d2(j), h2(j))
        end do
        do k = 1, 4
          xv(:, k) = mesh%xnode(:, mesh%ncon(k, p))
        end do
        do iv = 0, nsub - 1
          do iu = 0, nsub - 1
            us = [real(iu, 8), real(iu + 1, 8)]/nsub
            vs = [real(iv, 8), real(iv + 1, 8)]/nsub
            sv(:, 1) = bilinear(xv, us(1), vs(1))
            sv(:, 2) = bilinear(xv, us(2), vs(1))
            sv(:, 3) = bilinear(xv, us(2), vs(2))
            sv(:, 4) = bilinear(xv, us(1), vs(2))
            sc = 0.25d0*sum(sv, 2)
            call rankine_moments(sv, 4, mesh%nrm(:, p), x, sc, m0, g0, m1, mg)
            call rankine_moments(sv, 4, mesh%nrm(:, p), xim, sc, m0i, g0, m1, mg)
            call INFGREEN3D(sc(1), x(1), sc(2), x(2), sc(3), x(3), VWAVE, grn, 1)
            d = [sc(1) - x(1), sc(2) - x(2), sc(3) + x(3)]
            r1 = norm2(d)
            ! flag = 1: GRN(1) holds 1/r1 and the rest, as in green_panel_integrals
            out = out + dens*(m0 + m0i + (grn(1) - 1d0/r1)*mesh%area(p)/nsub**2)
          end do
        end do
      end do
    end function plain_form

  end subroutine test_b2_identity

  !> The point (u, v) in [0, 1]^2 of the quad xv, bilinear in its vertices
  function bilinear(xv, uu, vv) result(y)
    real(8), intent(in) :: xv(3, 4), uu, vv
    real(8) :: y(3)
    y = (1 - uu)*(1 - vv)*xv(:, 1) + uu*(1 - vv)*xv(:, 2) + uu*vv*xv(:, 3) + (1 - uu)*vv*xv(:, 4)
  end function bilinear

  ! ------------------------------------------------------------------------------------------------ B2 route

  !> The B2 part of direct_field_at_points (vf = 0, phi_D = 0, so phi = -r_B2/(4 pi)) takes the field-mode route of
  !> FINGREEN3D, as the rows do: at depth 10, where the routes differ, on the truncated cylinder with its waterline (so
  !> the line integral is in), with the densities of pair_densities, at the 6 points of the radiation identity. It equals
  !> -b2_stokes_at_point(field_mode = .true.)/(4 pi) to rounding, and differs from the default route (field_mode absent,
  !> as Phase 7a calls it) by more than route_min_diff, so a B2 part on the default route fails the first check.
  subroutine test_b2_field_route(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    complex(8), allocatable :: vf(:, :, :), a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :), phi_d(:, :), phi_field(:, :)
    complex(8) :: r_field(1), r_default(1), ref_field(6), ref_default(6)
    real(8) :: xp(3, 6), rel_field, rel_default
    integer :: ia, ib, ipt
    character(len=250) :: line

    call make_cylinder(cyl_radius, cyl_draft, 16, 4, 3, full, .true.)
    call set_wave(cyl_depth)
    call SetWaveFromOmega(om)
    NTHREAD = 4
    call load_panels(full, 0, 0, mesh)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    if (geom%nwl == 0) error stop 'test_b2_field_route: the cylinder has no waterline segments'
    call pair_densities(mesh, geom, cyl_depth, a, c, e, cl)
    allocate(vf(mesh%nbase, 1, mesh%nsys), phi_d(mesh%n, 1), phi_field(6, 1))
    vf = (0d0, 0d0)
    phi_d = (0d0, 0d0)
    ia = side_panel(mesh, -0.5d0*cyl_draft)
    ib = minloc(mesh%xc(3, :), 1)
    xp(:, 1) = [12d0, 3d0, -2d0]
    xp(:, 2) = [-4d0, -9d0, -7d0]
    xp(:, 3) = [8d0, -2d0, 0d0]
    xp(:, 4) = [-3d0, 11d0, 0d0]
    xp(:, 5) = mesh%xc(:, ia) - 0.5d0*mesh%psize(ia)*mesh%nrm(:, ia)
    xp(:, 6) = mesh%xc(:, ib) - 0.5d0*mesh%psize(ib)*mesh%nrm(:, ib)
    call direct_field_at_points(mesh, geom, xp, phi_d, vf, a, c, e, cl, .false., phi_field)
    do ipt = 1, 6
      call b2_stokes_at_point(mesh, geom, xp(:, ipt), a, c, e, .false., r_field, cl, field_mode=.true.)
      call b2_stokes_at_point(mesh, geom, xp(:, ipt), a, c, e, .false., r_default, cl)
      ref_field(ipt) = -r_field(1)/(4d0*pi)
      ref_default(ipt) = -r_default(1)/(4d0*pi)
    end do
    rel_field = maxval(abs(phi_field(:, 1) - ref_field))/maxval(abs(ref_field))
    rel_default = maxval(abs(phi_field(:, 1) - ref_default))/maxval(abs(ref_field))
    write(line, '(a,es10.2,a,es10.2,a,es10.2)') "B2 part against the field route ", rel_field, &
      ", against the default route ", rel_default, "; largest |phi| ", maxval(abs(ref_field))
    print '(4x,a)', trim(line)
    call check(error, rel_field <= route_limit .and. rel_default >= route_min_diff, line)
  end subroutine test_b2_field_route

  ! ------------------------------------------------------------------------------------------------ approach

  !> The direct solution with b = n_3 + B2 (the densities of pair_densities) on the truncated cylinder at depth 10, solved
  !> as direct_at_omega does (ASSB_RHS_FIELD, direct_rhs_b2_block from the stored integrals, direct_solve). The field at
  !> 2, 1, 0.5 and 0.25 panel sizes from the centroid of a side panel at mid draft, along -n into the fluid, tends to the
  !> panel value: the difference decreases and is below approach_limit at 0.25 (spec: phi(2) at field points tends to the
  !> direct solution on the body).
  subroutine test_approaches_panel(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    complex(8), allocatable :: vf(:, :, :), rf(:, :, :), r2(:, :, :), a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :), &
                               phi(:, :), phi_field(:, :)
    real(8), parameter :: dist(4) = [2d0, 1d0, 0.5d0, 0.25d0]
    real(8) :: xp(3, 4), diff(4)
    integer :: ip, k, nb
    character(len=250) :: line

    call make_cylinder(cyl_radius, cyl_draft, 16, 4, 3, full, .true.)
    call set_wave(cyl_depth)
    NTHREAD = 4
    ISOLV = 1
    call load_panels(full, 0, 0, mesh)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    nb = mesh%nbase
    call radiation_at(om, .false., nb)
    call panel_forcing(mesh, cmplx(mesh%nj(3, :), 0d0, kind=8), vf)
    call pair_densities(mesh, geom, cyl_depth, a, c, e, cl)
    allocate(rf(nb, 1, mesh%nsys), r2(nb, 1, mesh%nsys), phi(mesh%n, 1), phi_field(4, 1))
    call ASSB_RHS_FIELD(vf, rf, nb, mesh%nsys, 1)
    call direct_rhs_b2_block(mesh, geom, a, c, e, r2, .false., cl)
    rf = rf + r2
    call direct_solve(mesh, rf, phi, .false.)
    ip = side_panel(mesh, -0.5d0*cyl_draft)
    do k = 1, 4
      xp(:, k) = mesh%xc(:, ip) - dist(k)*mesh%psize(ip)*mesh%nrm(:, ip)
    end do
    call direct_field_at_points(mesh, geom, xp, phi, vf, a, c, e, cl, .false., phi_field)
    diff = abs(phi_field(:, 1) - phi(ip, 1))/abs(phi(ip, 1))
    write(line, '(a,i0,a,4es10.2)') "panel ", ip, ": |phi_field - phi_D| relative at 2, 1, 0.5, 0.25 panel sizes ", diff
    print '(4x,a)', trim(line)
    call check(error, diff(2) < diff(1) .and. diff(3) < diff(2) .and. diff(4) < diff(3) .and. diff(4) < approach_limit, &
               line)
  end subroutine test_approaches_panel

  ! ------------------------------------------------------------------------------------------------ symmetry

  !> The truncated cylinder whole and as half meshes (ISY = 1, then ISX = 1; HAMS takes one plane at a time, nsys <= 2) at
  !> depth 10 and Omega: the same panel data on the full panels (phi_D and vf smooth functions of the centroid, the B2
  !> densities of pair_densities) give the same phi_field at 4 points off the symmetry planes, one of them on z = 0.
  subroutine test_symmetry(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: full, mesh
    type(surface_geom) :: geom
    complex(8), allocatable :: vf(:, :, :), a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :), phi_d(:, :), phi_field(:, :), &
                               whole(:), pv(:)
    real(8) :: xp(3, 4), worst, rel
    integer :: icase, f
    integer, parameter :: sx(3) = [0, 0, 1], sy(3) = [0, 1, 0]
    character(len=200) :: line

    call make_cylinder(cyl_radius, cyl_draft, 16, 4, 3, full, .true.)
    xp(:, 1) = [8d0, 3d0, -2d0]
    xp(:, 2) = [-7d0, -4d0, -1d0]
    xp(:, 3) = [3d0, -9d0, -6d0]
    xp(:, 4) = [-6d0, 6d0, 0d0]
    worst = 0d0
    do icase = 1, 3
      call set_wave(cyl_depth)
      call SetWaveFromOmega(om)
      NTHREAD = 4
      call load_panels(full, sx(icase), sy(icase), mesh)
      call geometry(mesh, geom, error)
      if (allocated(error)) return
      allocate(phi_d(mesh%n, 1), pv(mesh%n), phi_field(4, 1))
      do f = 1, mesh%n
        phi_d(f, 1) = smooth(mesh%xc(:, f), 0.3d0)
        pv(f) = smooth(mesh%xc(:, f), -0.7d0)
      end do
      call panel_forcing(mesh, pv, vf)
      call pair_densities(mesh, geom, cyl_depth, a, c, e, cl)
      call direct_field_at_points(mesh, geom, xp, phi_d, vf, a, c, e, cl, .false., phi_field)
      if (icase == 1) then
        whole = phi_field(:, 1)
      else
        rel = maxval(abs(phi_field(:, 1) - whole))/maxval(abs(whole))
        write(line, '(a,i0,a,i0,a,es10.2)') "ISX = ", sx(icase), ", ISY = ", sy(icase), &
          ": phi_field against the whole mesh, relative ", rel
        print '(4x,a)', trim(line)
        worst = max(worst, rel)
      end if
      deallocate(phi_d, pv, phi_field, vf, a, c, e, cl)
    end do
    ISX = 0; ISY = 0; NSYS = 1
    write(line, '(a,es10.2)') "half meshes against the whole, worst relative difference ", worst
    call check(error, worst <= symmetry_limit, line)

  contains

    !> A smooth complex function of the position with no symmetry in x or y
    complex(8) function smooth(x, s)
      real(8), intent(in) :: x(3), s
      smooth = cmplx(cos(0.31d0*x(1) + 0.17d0*x(2) + s), sin(0.23d0*x(2) - 0.11d0*x(1) + 0.2d0*x(3)), kind=8)
    end function smooth

  end subroutine test_symmetry

  ! ------------------------------------------------------------------------------------------------ the hook

  !> The fixture of the hook tests: the truncated cylinder at depth 10 (test_qtf_direct), the radiation solve at
  !> Omega = w1 - w2 (or omega, for the sum), the two states with the motions of test_qtf_potential and the pairs (1, 2) with
  !> Omega > 0 and (2, 1)
  !> with Omega < 0; xp holds 4 points off the hull.
  subroutine hook_fixture(mesh, geom, store, store_wl, pairs, xp, error, omega)
    type(panel_mesh), intent(out) :: mesh
    type(surface_geom), intent(out) :: geom
    type(first_order_state), allocatable, intent(out) :: store(:, :)
    complex(8), allocatable, intent(out) :: store_wl(:, :, :)
    integer, allocatable, intent(out) :: pairs(:, :)
    real(8), intent(out) :: xp(3, 4)
    type(error_type), allocatable, intent(out) :: error
    real(8), intent(in), optional :: omega
    type(panel_mesh) :: full
    complex(8) :: x1(6, 1), x2(6, 1)
    complex(8), allocatable :: wl1(:), wl2(:)

    call make_cylinder(cyl_radius, cyl_draft, 16, 4, 3, full, .true.)
    call set_wave(cyl_depth)
    NTHREAD = 4
    ISOLV = 1
    call load_panels(full, 0, 0, mesh)
    call geometry(mesh, geom, error)
    if (allocated(error)) return
    call radiation_at(merge(omega, w1 - w2, present(omega)), .false., mesh%nbase)
    call direct_setup(mesh, DEFAULT_EDGE_DEG)
    call motions(x1, x2)
    allocate(store(2, 1), store_wl(geom%nwl, 2, 1))
    call make_state(mesh, geom, w1, b1, cyl_depth, x1, store(1, 1), wl1)
    call make_state(mesh, geom, w2, b2, cyl_depth, x2, store(2, 1), wl2)
    store_wl(:, 1, 1) = wl1
    store_wl(:, 2, 1) = wl2
    pairs = reshape([1, 1, 2, 1, 2, 1, 1, 1], [4, 2])
    xp(:, 1) = [8d0, 3d0, -2d0]
    xp(:, 2) = [-7d0, -4d0, -1d0]
    xp(:, 3) = [3d0, -9d0, -6d0]
    xp(:, 4) = [-6d0, 6d0, 0d0]
  end subroutine hook_fixture

  !> phi_field of the hook of direct_at_omega, for both pairs (one block) with phi_I2 replaced by zero (the test hook
  !> phi2_constant, so that the panel forcing needs no incident field), against the same columns built here: the forcing
  !> B1 + B3 of pair_forcing, the B2 densities, ASSB_RHS_FIELD_IRR and direct_rhs_b2_block, one direct_solve of the block,
  !> direct_field_at_points, and the conjugation of the column of Omega < 0. The same operations on the same data, so equal
  !> to the bit, for the pair with Omega > 0 and the pair with Omega < 0.
  subroutine test_hook_matches(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    type(surface_geom) :: geom
    type(first_order_state), allocatable :: store(:, :)
    complex(8), allocatable :: store_wl(:, :, :), t(:, :), phi_hook(:, :), phi_ref(:, :), vf(:, :, :), rf(:, :, :), &
                               r2(:, :, :), a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :), phi(:, :)
    integer, allocatable :: pairs(:, :)
    logical, allocatable :: mask(:)
    real(8) :: xp(3, 4), dw, nv(3), x(3)
    complex(8) :: d1(3), d2(3), v1(3), v2(3), q1(3), q2(3), hx(3), b1v, b3v, field, zero3(3)
    integer :: p, f, b, n, nb, ns
    character(len=200) :: line

    call hook_fixture(mesh, geom, store, store_wl, pairs, xp, error)
    if (allocated(error)) return
    n = mesh%n
    nb = mesh%nbase
    ns = mesh%nsys
    allocate(t(6, 2), phi_hook(4, 2), phi_ref(4, 2))
    call direct_at_omega(mesh, geom, w1 - w2, store, store_wl, pairs, t, phi2_constant=(0d0, 0d0), &
                         field_points=xp, phi_field=phi_hook)

    allocate(vf(nb, 2, ns), rf(nb, 2, ns), r2(nb, 2, ns), a(n, 2), c(3, n, 2), e(geom%nwl, 2), cl(3, 4, n, 2), phi(n, 2))
    mask = direct_edge_mask()
    zero3 = (0d0, 0d0)
    do p = 1, 2
      associate (s1 => store(pairs(1, p), pairs(2, p)), s2 => store(pairs(3, p), pairs(4, p)))
        dw = s1%omega - s2%omega
        do f = 1, n
          b = mesh%body(f)
          nv = mesh%nrm(:, f)
          x = mesh%xc(:, f) - mesh%xr(:, b)
          d1 = rigid_displacement(s1%xi(:, b), x)
          d2 = rigid_displacement(s2%xi(:, b), x)
          v1 = -cmplx(0d0, 1d0, kind=8)*s1%omega*d1 - (s1%grad_i(:, f) + s1%grad_p(:, f))
          v2 = -cmplx(0d0, 1d0, kind=8)*s2%omega*d2 - (s2%grad_i(:, f) + s2%grad_p(:, f))
          b1v = sum(nv*pair_cross(KERNEL_DIFF, v1, s1%xi(4:6, b), v2, s2%xi(4:6, b)))
          q1 = rigid_displacement([zero3, s1%xi(4:6, b)], x)
          q2 = rigid_displacement([zero3, s2%xi(4:6, b)], x)
          hx = 0.5d0*pair_cross(KERNEL_DIFF, s1%xi(4:6, b), q1, s2%xi(4:6, b), q2)
          b3v = -cmplx(0d0, 1d0, kind=8)*dw*sum(nv*hx)
          field = b1v + b3v - (0d0, 0d0)
          if (dw < 0d0) field = conjg(field)
          vf(mesh%base(f), p, mesh%image(f)) = field
        end do
        call b2_stokes_densities(mesh, geom, s1, s2, store_wl(:, pairs(1, p), pairs(2, p)), &
                                 store_wl(:, pairs(3, p), pairs(4, p)), a(:, p), c(:, :, p), e(:, p), cl(:, :, :, p), mask)
        if (dw < 0d0) then
          a(:, p) = conjg(a(:, p))
          c(:, :, p) = conjg(c(:, :, p))
          e(:, p) = conjg(e(:, p))
          cl(:, :, :, p) = conjg(cl(:, :, :, p))
        end if
      end associate
    end do
    call ASSB_RHS_FIELD_IRR(vf, rf, nb, 0, ns, 2)
    call direct_rhs_b2_block(mesh, geom, a, c, e, r2, .false., cl)
    rf = rf + r2
    call direct_solve(mesh, rf, phi, .false.)
    call direct_field_at_points(mesh, geom, xp, phi, vf, a, c, e, cl, .false., phi_ref)
    phi_ref(:, 2) = conjg(phi_ref(:, 2))
    write(line, '(a,es10.2,a,es10.2,a,es10.2)') "hook against the separate call, largest difference ", &
      maxval(abs(phi_hook - phi_ref)), "; largest |phi_field| ", maxval(abs(phi_ref(:, 1))), " and ", &
      maxval(abs(phi_ref(:, 2)))
    print '(4x,a)', trim(line)
    call check(error, all(phi_hook == phi_ref) .and. maxval(abs(phi_ref(:, 1))) > 0d0 .and. &
               maxval(abs(phi_ref(:, 2))) > 0d0, line)
  end subroutine test_hook_matches

  !> The hook for sum columns (Phase 7c-2): phi_field of two sum pairs (1, 2) and (2, 1) at Omega+ = w1 + w2 against the same
  !> columns built here as in test_hook_matches, with the sum kernel (pair_cross and b2_stokes_densities with KERNEL_SUM,
  !> B3 with dw = Omega+) and no conjugation (Omega+ > 0): the same operations on the same data, so equal to the bit.
  subroutine test_hook_sum_matches(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    type(surface_geom) :: geom
    type(first_order_state), allocatable :: store(:, :)
    complex(8), allocatable :: store_wl(:, :, :), t(:, :), phi_hook(:, :), phi_ref(:, :), vf(:, :, :), rf(:, :, :), &
                               r2(:, :, :), a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :), phi(:, :)
    integer, allocatable :: pairs(:, :)
    logical, allocatable :: mask(:)
    real(8) :: xp(3, 4), dw, nv(3), x(3)
    complex(8) :: d1(3), d2(3), v1(3), v2(3), q1(3), q2(3), hx(3), b1v, b3v, field, zero3(3)
    integer :: p, f, b, n, nb, ns
    character(len=200) :: line

    call hook_fixture(mesh, geom, store, store_wl, pairs, xp, error, w1 + w2)
    if (allocated(error)) return
    n = mesh%n
    nb = mesh%nbase
    ns = mesh%nsys
    allocate(t(6, 2), phi_hook(4, 2), phi_ref(4, 2))
    call direct_at_omega(mesh, geom, w1 + w2, store, store_wl, pairs, t, phi2_constant=(0d0, 0d0), &
                         field_points=xp, phi_field=phi_hook, kernels=[KERNEL_SUM, KERNEL_SUM])

    allocate(vf(nb, 2, ns), rf(nb, 2, ns), r2(nb, 2, ns), a(n, 2), c(3, n, 2), e(geom%nwl, 2), cl(3, 4, n, 2), phi(n, 2))
    mask = direct_edge_mask()
    zero3 = (0d0, 0d0)
    do p = 1, 2
      associate (s1 => store(pairs(1, p), pairs(2, p)), s2 => store(pairs(3, p), pairs(4, p)))
        dw = s1%omega + s2%omega
        do f = 1, n
          b = mesh%body(f)
          nv = mesh%nrm(:, f)
          x = mesh%xc(:, f) - mesh%xr(:, b)
          d1 = rigid_displacement(s1%xi(:, b), x)
          d2 = rigid_displacement(s2%xi(:, b), x)
          v1 = -cmplx(0d0, 1d0, kind=8)*s1%omega*d1 - (s1%grad_i(:, f) + s1%grad_p(:, f))
          v2 = -cmplx(0d0, 1d0, kind=8)*s2%omega*d2 - (s2%grad_i(:, f) + s2%grad_p(:, f))
          b1v = sum(nv*pair_cross(KERNEL_SUM, v1, s1%xi(4:6, b), v2, s2%xi(4:6, b)))
          q1 = rigid_displacement([zero3, s1%xi(4:6, b)], x)
          q2 = rigid_displacement([zero3, s2%xi(4:6, b)], x)
          hx = 0.5d0*pair_cross(KERNEL_SUM, s1%xi(4:6, b), q1, s2%xi(4:6, b), q2)
          b3v = -cmplx(0d0, 1d0, kind=8)*dw*sum(nv*hx)
          field = b1v + b3v - (0d0, 0d0)
          if (dw < 0d0) field = conjg(field)
          vf(mesh%base(f), p, mesh%image(f)) = field
        end do
        call b2_stokes_densities(mesh, geom, s1, s2, store_wl(:, pairs(1, p), pairs(2, p)), &
                                 store_wl(:, pairs(3, p), pairs(4, p)), a(:, p), c(:, :, p), e(:, p), cl(:, :, :, p), mask, &
                                 KERNEL_SUM)
        if (dw < 0d0) then
          a(:, p) = conjg(a(:, p))
          c(:, :, p) = conjg(c(:, :, p))
          e(:, p) = conjg(e(:, p))
          cl(:, :, :, p) = conjg(cl(:, :, :, p))
        end if
      end associate
    end do
    call ASSB_RHS_FIELD_IRR(vf, rf, nb, 0, ns, 2)
    call direct_rhs_b2_block(mesh, geom, a, c, e, r2, .false., cl)
    rf = rf + r2
    call direct_solve(mesh, rf, phi, .false.)
    call direct_field_at_points(mesh, geom, xp, phi, vf, a, c, e, cl, .false., phi_ref)
    write(line, '(a,es10.2,a,es10.2,a,es10.2)') "sum hook against the separate call, largest difference ", &
      maxval(abs(phi_hook - phi_ref)), "; largest |phi_field| ", maxval(abs(phi_ref(:, 1))), " and ", &
      maxval(abs(phi_ref(:, 2)))
    print '(4x,a)', trim(line)
    call check(error, all(phi_hook == phi_ref) .and. maxval(abs(phi_ref(:, 1))) > 0d0 .and. &
               maxval(abs(phi_ref(:, 2))) > 0d0, line)
  end subroutine test_hook_sum_matches

  !> The pair (2, 1) is the pair (1, 2) with the waves exchanged, so its phi_D is the complex conjugate (T_dir of the
  !> negative frequency is the conjugate of that of the positive one). phi_field of the hook, with the real phi_I2, for the
  !> two pairs of one block: (2, 1) against conj of (1, 2).
  subroutine test_negative_omega(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    type(surface_geom) :: geom
    type(first_order_state), allocatable :: store(:, :)
    complex(8), allocatable :: store_wl(:, :, :), t(:, :), phi_field(:, :)
    integer, allocatable :: pairs(:, :)
    real(8) :: xp(3, 4), rel
    character(len=200) :: line

    call hook_fixture(mesh, geom, store, store_wl, pairs, xp, error)
    if (allocated(error)) return
    allocate(t(6, 2), phi_field(4, 2))
    call direct_at_omega(mesh, geom, w1 - w2, store, store_wl, pairs, t, field_points=xp, phi_field=phi_field)
    rel = maxval(abs(phi_field(:, 2) - conjg(phi_field(:, 1))))/maxval(abs(phi_field(:, 1)))
    write(line, '(a,es10.2,a,es10.2)') "phi_field of (2, 1) against conj of (1, 2), relative ", rel, &
      "; largest |phi_field| ", maxval(abs(phi_field(:, 1)))
    print '(4x,a)', trim(line)
    call check(error, rel <= conjugate_limit .and. maxval(abs(phi_field(:, 1))) > 0d0, line)
  end subroutine test_negative_omega

  !> Without field_points and phi_field, tdir and parts are those of the call with them, to the bit: the hook adds
  !> nothing to the solve. (The values themselves are those pinned by the tests of test_qtf_direct, which run the routine
  !> without the arguments.)
  subroutine test_hook_absent(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    type(surface_geom) :: geom
    type(first_order_state), allocatable :: store(:, :)
    complex(8), allocatable :: store_wl(:, :, :), t0(:, :), t1(:, :), parts0(:, :, :), parts1(:, :, :), phi_field(:, :)
    integer, allocatable :: pairs(:, :)
    real(8) :: xp(3, 4)
    character(len=200) :: line

    call hook_fixture(mesh, geom, store, store_wl, pairs, xp, error)
    if (allocated(error)) return
    allocate(t0(6, 2), t1(6, 2), parts0(6, 2, 2), parts1(6, 2, 2), phi_field(4, 2))
    call direct_at_omega(mesh, geom, w1 - w2, store, store_wl, pairs, t0, parts=parts0)
    call direct_at_omega(mesh, geom, w1 - w2, store, store_wl, pairs, t1, parts=parts1, field_points=xp, phi_field=phi_field)
    write(line, '(a,es10.2,a,es10.2)') "without the hook against with it, largest difference in tdir ", &
      maxval(abs(t0 - t1)), " and in parts ", maxval(abs(parts0 - parts1))
    print '(4x,a)', trim(line)
    call check(error, all(t0 == t1) .and. all(parts0 == parts1) .and. maxval(abs(t0)) > 0d0, line)
  end subroutine test_hook_absent

  ! ------------------------------------------------------------------------------------------------ two bodies

  !> The radiation identity of test_radiation_identity for two cylinders in the combined multi-body arrays (load_multi,
  !> is_multi = .true.): with vf = n_3 on the panels of body ib (0 on the other) and phi_D = psi of the heave mode of body
  !> ib from the multi-body radiation solve, Green's formula of direct_field_at_points is that of fields_green with
  !> multi = .true. (phi_p = psi, dphi_p/dn = vf, no incident wave). Heave of both bodies, depth 10 and deep water, 6 points
  !> off the hulls: between the bodies, beside each, one on z = 0, and half a panel size off a side panel of each body.
  subroutine test_two_bodies(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    type(surface_geom) :: geom
    type(first_order_state) :: states(1)
    type(cs_points) :: pts
    type(cs_fields) :: fields(1)
    complex(8), allocatable :: psi(:, :), vf(:, :, :), a(:, :), c(:, :, :), e(:, :), cl(:, :, :, :), phi_field(:, :), vn(:)
    real(8) :: xp(3, 6), worst, rel
    real(8), parameter :: depths(2) = [cyl_depth, -1d0]
    integer :: icase, ib, k, ia(2), f
    character(len=200) :: line

    worst = 0d0
    do icase = 1, 2
      call set_wave(depths(icase))
      NTHREAD = 4
      ISOLV = 1
      call load_multi(0, mesh)
      call geometry(mesh, geom, error)
      if (allocated(error)) return
      call radiation_at(om, .true., mesh%nbase)
      allocate(psi(mesh%n, 6*mesh%nbody), vn(mesh%n))
      call radiation_psi(mesh, om, MXPOT_MULTI_COMB, psi)
      do ib = 1, 2
        ia(ib) = side_panel_of_body(mesh, ib, -0.5d0*cyl_draft)
      end do
      xp(:, 1) = [15d0, 3d0, -2d0]                 ! between the bodies (centres at x = 0 and 30)
      xp(:, 2) = [-9d0, -4d0, -7d0]                ! beside body 1
      xp(:, 3) = [40d0, 2d0, -3d0]                 ! beside body 2
      xp(:, 4) = [17d0, -6d0, 0d0]                 ! on z = 0
      xp(:, 5) = mesh%xc(:, ia(1)) - 0.5d0*mesh%psize(ia(1))*mesh%nrm(:, ia(1))
      xp(:, 6) = mesh%xc(:, ia(2)) - 0.5d0*mesh%psize(ia(2))*mesh%nrm(:, ia(2))
      do ib = 1, 2
        k = 6*(ib - 1) + 3
        ! b = n_3 on the panels of body ib, 0 on the other (the layout of the multi-body unit normal fields)
        vn = (0d0, 0d0)
        do f = 1, mesh%n
          if (mesh%body(f) == ib) vn(f) = cmplx(mesh%nj(3, f), 0d0, kind=8)
        end do
        call panel_forcing(mesh, vn, vf)
        call zero_densities(mesh, geom, a, c, e, cl)
        allocate(phi_field(6, 1))
        call direct_field_at_points(mesh, geom, xp, psi(:, k:k), vf, a, c, e, cl, .true., phi_field)
        states(1)%omega = om
        states(1)%k = WK
        states(1)%h = depths(icase)
        states(1)%nu = VWAVE
        states(1)%amp = 0d0
        states(1)%phi_p = psi(:, k)
        states(1)%dphi_p_dn = vn
        pts%n = 6
        pts%x = xp
        allocate(pts%nrm(3, 6), pts%w(6), pts%psize(6), pts%kind(6), pts%body(6))
        pts%nrm = 0d0
        pts%w = 1d0
        pts%psize = 1d0
        pts%kind = CS_LINE
        pts%body = 1
        call fields_green(mesh, pts, states, 0.01d0, .true., fields)
        rel = maxval(abs(phi_field(:, 1) - fields(1)%phi))/maxval(abs(fields(1)%phi))
        write(line, '(a,f6.1,a,i0,a,es10.2,a,es10.2)') "depth ", depths(icase), ", heave of body ", ib, &
          ": phi_field against fields_green, relative ", rel, "; largest |phi| ", maxval(abs(fields(1)%phi))
        print '(4x,a)', trim(line)
        worst = max(worst, rel)
        deallocate(phi_field, vf, a, c, e, cl, pts%x, pts%nrm, pts%w, pts%psize, pts%kind, pts%body)
      end do
      deallocate(psi, vn)
    end do
    write(line, '(a,es10.2)') "two bodies, radiation identity, worst relative difference ", worst
    call check(error, worst <= two_body_limit, line)
  end subroutine test_two_bodies

  !> The index of a panel of body ib on the side wall (horizontal normal) in the row nearest z = z0.
  integer function side_panel_of_body(mesh, ib, z0) result(ip)
    type(panel_mesh), intent(in) :: mesh
    integer, intent(in) :: ib
    real(8), intent(in) :: z0
    real(8) :: best
    integer :: i
    ip = 0
    best = huge(1d0)
    do i = 1, mesh%n
      if (mesh%body(i) /= ib .or. abs(mesh%nrm(3, i)) > 0.1d0) cycle
      if (abs(mesh%xc(3, i) - z0) < best) then
        best = abs(mesh%xc(3, i) - z0)
        ip = i
      end if
    end do
  end function side_panel_of_body

  ! ------------------------------------------------------------------------------------------------ the lid

  !> The hook with irregular frequency removal (IRSP = 1, the lid points of the cylinder waterplane) against the same call
  !> without it (IRSP = 0): the lid is inside the body, so it takes no part in Green's formula over S, and only the solved
  !> phi_D differs (the least-squares solve with the lid rows against the plain one; mesh-level, far from an irregular
  !> frequency). The pairs (1, 2) and (2, 1) at Omega = w1 - w2, the real phi_I2, 4 points off the hull. The difference
  !> must be non-zero (the lid solve really ran) and at most lid_limit.
  subroutine test_lid(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    type(surface_geom) :: geom
    type(first_order_state), allocatable :: store(:, :)
    complex(8), allocatable :: store_wl(:, :, :), t0(:, :), t1(:, :), phi0(:, :), phi1(:, :)
    integer, allocatable :: pairs(:, :)
    real(8), allocatable :: lid(:, :)
    real(8) :: xp(3, 4), rel, tdiff
    character(len=200) :: line

    call hook_fixture(mesh, geom, store, store_wl, pairs, xp, error)
    if (allocated(error)) return
    allocate(t0(6, 2), t1(6, 2), phi0(4, 2), phi1(4, 2))
    call direct_at_omega(mesh, geom, w1 - w2, store, store_wl, pairs, t0, field_points=xp, phi_field=phi0)
    call append_lid(cyl_radius, 16, 3, 0d0, 0, lid)
    call set_lid(lid, .false.)
    call radiation_at_irr(w1 - w2, .false., mesh%nbase, size(lid, 2))
    call direct_at_omega(mesh, geom, w1 - w2, store, store_wl, pairs, t1, field_points=xp, phi_field=phi1)
    call clear_lid()
    rel = maxval(abs(phi1 - phi0))/maxval(abs(phi0))
    tdiff = maxval(abs(t1 - t0))/maxval(abs(t0))
    write(line, '(a,i0,a,es10.2,a,es10.2,a,es10.2)') "IRSP = 1 (", size(lid, 2), " lid points) against IRSP = 0: phi_field ", &
      rel, ", T_dir ", tdiff, "; largest |phi_field| ", maxval(abs(phi0))
    print '(4x,a)', trim(line)
    call check(error, rel >= lid_min_diff .and. rel <= lid_limit .and. maxval(abs(phi0)) > 0d0, line)
  end subroutine test_lid

end module test_qtf_direct_field
