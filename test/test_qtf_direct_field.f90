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
  use Potentials_mod, only : MXPOT
  use AssbMatx, only : ASSB_RHS_FIELD
  use SecondOrderTypes, only : panel_mesh, surface_geom, first_order_state, cs_points, cs_fields, CS_LINE
  use FirstOrderState, only : radiation_psi
  use FieldVelocity, only : fields_green
  use NearFieldDrift, only : pair_product, rigid_displacement, KERNEL_DIFF
  use QTFDirect, only : b2_stokes_densities, direct_rhs_b2_block, direct_solve, rankine_moments, sharp_edge_panels
  use QTFDirectField, only : direct_field_at_points
  use SharpEdge, only : DEFAULT_EDGE_DEG
  use INFG3D_Open, only : INFGREEN3D
  use second_order_fixtures, only : make_cylinder
  use test_waterline_field, only : set_wave, load_panels
  use test_qtf_potential, only : motions, make_state, geometry
  use test_qtf_direct, only : radiation_at
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
  ! The field 0.25 panel sizes off a panel centroid against the panel value (spec, Phase 7 item 2); 1.5e-3 measured
  real(8), parameter :: approach_limit = 2d-2
  ! Half mesh against whole: the same sums in another order, and the Green function at mirrored points; 1.7e-14 measured
  real(8), parameter :: symmetry_limit = 1d-9

contains

  subroutine collect_qtf_direct_field(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("direct_field_radiation_identity", test_radiation_identity), &
      new_unittest("direct_field_b2_identity", test_b2_identity), &
      new_unittest("direct_field_approaches_panel", test_approaches_panel), &
      new_unittest("direct_field_symmetry", test_symmetry) &
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

end module test_qtf_direct_field
