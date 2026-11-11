!> Unit tests of the right side of the direct second-order problem (Phase 7a): ASSB_RHS_FIELD against ASSB_RBC, and the
!> B2 part (QTFDirect) in the Stokes form with G against the direct form int G n.[(X.grad)u] dS on analytic fields, and
!> the closed-form panel moments of 1/r (rankine_moments) against adaptive quadrature.
!> These tests set module globals (the mesh arrays, the stored Green function, the wave), so the suite runs serially
!> (see unit_tester.f90).
module test_qtf_direct
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use HAMS_mod, only : NTHREAD
  use WaveDyn_mod, only : VWAVE => V
  use PanelMesh_mod, only : ISX, ISY, NSYS, NELEM, DS, PNSZ, XYZ_P, DXYZ_P, DS_MULTI_COMB, PNSZ_MULTI_COMB, &
                            XYZ_GLOBAL_MULTI_COMB_P, DXYZ_MULTI_COMB
  use Inerfs_mod, only : TNELEM
  use Potentials_mod, only : CGRN, RKBN, CGRN_MULTI_COMB, RKBN_MULTI_COMB
  use Const_mod, only : RXY
  use AssbMatx, only : ASSB_RBC, ASSB_RHS_FIELD
  use AssbMatxMulti, only : ASSB_RBC_MULTI, ASSB_RHS_FIELD_MULTI
  use SecondOrderTypes, only : panel_mesh, surface_geom, first_order_state
  use FirstOrderState, only : image_signs
  use SurfaceVelocity, only : tangential_fit
  use NearFieldDrift, only : pair_product, rigid_displacement, KERNEL_DIFF
  use QTFDirect, only : direct_rhs_b2, b2_stokes_densities, b2_stokes_at_point, green_panel_integrals, rankine_moments
  use second_order_fixtures, only : make_cylinder
  use test_waterline_field, only : set_wave, load_panels
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
  ! Sphere of radius 5 at centre depth 20: 128, 512 and 2048 panels
  integer, parameter :: sph_nlat(3) = [8, 16, 32]

contains

  subroutine collect_qtf_direct(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("rhs_field_matches_rbc", test_rhs_field_matches_rbc), &
      new_unittest("rhs_field_matches_rbc_multi", test_rhs_field_matches_rbc_multi), &
      new_unittest("rankine_moments_quadrature", test_rankine_moments), &
      new_unittest("rhs_b2_stokes_vs_direct_smooth", test_b2_smooth), &
      new_unittest("rhs_b2_stokes_vs_direct_off_surface", test_b2_off_surface), &
      new_unittest("rhs_b2_half_matches_whole", test_b2_half) &
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
          if (io == 3 .and. ih == 1) cycle      ! outside the outline only above the plane
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
    call direct_rhs_b2(whole, geom_w, s1w, s2w, wl1w, wl2w, rw)
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
      call direct_rhs_b2(half, geom_h, s1h, s2h, wl1h, wl2h, rh)
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

end module test_qtf_direct
