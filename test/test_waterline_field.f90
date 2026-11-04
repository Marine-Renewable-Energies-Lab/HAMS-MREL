!> The perturbation potential by field evaluation against CalPotential and CalPotentialMulti, on the same panels, the same
!> potentials and the same wave. The panel potentials are synthetic, since the routines only apply a linear formula to them.
!> These tests set module globals, so the suite runs serially (see unit_tester.f90). It runs last, because it leaves the mesh
!> arrays at its own sizes.
module test_waterline_field
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use HAMS_mod, only : NTHREAD, ISOL, IRSP
  use WaveDyn_mod, only : AMP, BETA, H, NBODY, NK, TP, V, W1, WK, WVN
  use Body_mod, only : XW
  use Const_mod, only : G, RHO
  use PanelMesh_mod, only : DS, DXYZ_MULTI_COMB, DXYZ_P, DS_MULTI_COMB, GEONRM_MULTI_COMB, GEONRM_P, ISX, ISY, NCN, NCN_MULTI_COMB, &
                            NCON, NCON_MULTI_COMB, NELEM, NSYS, NTND, PNSZ, PNSZ_MULTI_COMB, XYZ, XYZ_GLOBAL_MULTI_COMB, &
                            XYZ_GLOBAL_MULTI_COMB_P, XYZ_P
  use Inerfs_mod, only : NELEM_PE
  use Potentials_mod, only : MXPOT, MXPOT_MULTI_COMB
  use PatcVelct, only : VINP
  use PressureElevation, only : CalPotential
  use PressureElevationMulti, only : CalPotentialMulti
  use SecondOrderTypes, only : panel_mesh, first_order_state, surface_geom
  use FirstOrderState, only : expand_images, build_state, incident_phi, wave_gravity
  use SurfaceGeometry, only : build_surface_geom, DEFAULT_SHARP_DEG
  use NearFieldDrift, only : near_field_pair, T_WATERLINE
  use WaterlineField, only : field_influence_row, field_potential, waterline_rows, waterline_field, field_point_inside_body
  use second_order_fixtures, only : make_cylinder, maccamy_fuchs_phi, maccamy_fuchs_state, maccamy_fuchs_waterline
  implicit none
  private

  public :: collect_waterline_field
  public :: set_wave, load_single, load_multi, load_panels

  complex(8), parameter :: ci = (0d0, 1d0)
  ! The two sides differ only in the order of the sums and the panels' symmetry bookkeeping
  real(8), parameter :: tol = 1d-8
  real(8), parameter :: radius = 5d0, draft = 5d0
  integer, parameter :: ntheta = 12, nz = 2, nr = 2
  real(8), parameter :: omega = 1.1d0, heading = 0.6d0
  real(8), parameter :: offsets(2) = [0.01d0, 0.1d0]
  ! Further than 50 panel sizes from every panel, so the routines use the Green function with no panel integrals
  real(8), parameter :: far_point(3) = [400d0, 300d0, -2d0]
  real(8), parameter :: body_spacing = 30d0   ! along x, so that a mirror image in y stays on its own body

contains

  subroutine collect_waterline_field(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("field_matches_existing", test_finite_depth), &
      new_unittest("field_matches_existing_deep", test_deep), &
      new_unittest("field_matches_existing_y_symmetry", test_y_symmetry), &
      new_unittest("field_matches_existing_x_symmetry", test_x_symmetry), &
      new_unittest("field_matches_existing_multi", test_multi), &
      new_unittest("field_matches_existing_multi_y_symmetry", test_multi_y_symmetry), &
      new_unittest("waterline_reference", test_waterline_reference), &
      new_unittest("field_point_inside_body", test_inside_body) &
      ]
  end subroutine collect_waterline_field

  !> The limit is a wall that widens downward by 45 degrees: n_z = -sqrt(1/2) for a normal into the body.
  subroutine test_inside_body(error)
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: limit = sqrt(0.5d0)

    call check(error, field_point_inside_body(0d0), .false., "vertical wall")
    if (allocated(error)) return
    call check(error, field_point_inside_body(0.5d0), .false., "wall narrowing downward")
    if (allocated(error)) return
    call check(error, field_point_inside_body(-0.9d0*limit), .false., "wall widening by less than 45 degrees")
    if (allocated(error)) return
    call check(error, field_point_inside_body(-1.1d0*limit), .true., "wall widening by more than 45 degrees")
    if (allocated(error)) return
    call check(error, field_point_inside_body(-1d0), .true., "horizontal floor-like panel")
    if (allocated(error)) return
    call check(error, all(field_point_inside_body([0d0, -0.5d0, -0.8d0]) .eqv. [.false., .false., .true.]), "array form")
  end subroutine test_inside_body

  subroutine test_finite_depth(error)
    type(error_type), allocatable, intent(out) :: error
    call compare_single(10d0, 0, 0, error)
  end subroutine test_finite_depth

  subroutine test_deep(error)
    type(error_type), allocatable, intent(out) :: error
    call compare_single(-1d0, 0, 0, error)
  end subroutine test_deep

  subroutine test_y_symmetry(error)
    type(error_type), allocatable, intent(out) :: error
    call compare_single(10d0, 0, 1, error)
  end subroutine test_y_symmetry

  subroutine test_x_symmetry(error)
    type(error_type), allocatable, intent(out) :: error
    call compare_single(-1d0, 1, 0, error)
  end subroutine test_x_symmetry

  subroutine test_multi(error)
    type(error_type), allocatable, intent(out) :: error
    call compare_multi(10d0, 0, error)
  end subroutine test_multi

  subroutine test_multi_y_symmetry(error)
    type(error_type), allocatable, intent(out) :: error
    call compare_multi(-1d0, 1, error)
  end subroutine test_multi_y_symmetry

  !> The analytic value of the integration test is the near-field term with the exact potential, on a fine waterline.
  subroutine test_waterline_reference(error)
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: a = 5d0, depth = 10d0, om = 1.6d0, rho_local = 1000d0, g_local = 9.81d0
    type(panel_mesh) :: mesh
    type(surface_geom) :: geom
    type(first_order_state) :: s
    complex(8), allocatable :: phi_wl(:), terms(:,:,:)
    real(8) :: k
    logical :: ok
    character(len=:), allocatable :: msg
    integer :: j, it

    k = om**2/g_local
    do it = 1, 100
      k = k - (k*tanh(k*depth) - om**2/g_local)/(tanh(k*depth) + k*depth/cosh(k*depth)**2)
    end do
    call make_cylinder(a, depth, 360, 2, 1, mesh, .false.)
    call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
    call check(error, ok, "build_surface_geom stopped")
    if (allocated(error)) return
    call maccamy_fuchs_state(mesh, a, depth, k, om, 1d0, 0d0, g_local, rho_local, s)
    allocate (phi_wl(geom%nwl), terms(6, 6, 1))
    do j = 1, geom%nwl
      phi_wl(j) = maccamy_fuchs_phi(geom%wl_xm(:, j), a, depth, k, om, 1d0, 0d0, g_local)
    end do
    call near_field_pair(mesh, geom, s, s, phi_wl, phi_wl, terms)
    call check(error, abs(real(terms(1, T_WATERLINE, 1), 8)/(rho_local*g_local) &
                          /maccamy_fuchs_waterline(a, depth, k, om, g_local, 720) - 1d0) < 1d-3, &
               "maccamy_fuchs_waterline differs from the near-field term with the exact potential")
  end subroutine test_waterline_reference

  ! ---------------------------------------------------------------------------------------------- the comparisons

  subroutine compare_single(depth, isx_in, isy_in, error)
    real(8), intent(in) :: depth
    integer, intent(in) :: isx_in, isy_in
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    type(surface_geom) :: geom
    type(first_order_state) :: s
    complex(8) :: dspl(6)

    call set_wave(depth)
    call load_single(isx_in, isy_in, mesh)
    call motions(dspl)
    call build_state(mesh, W1, WK, V, H, 1d0, BETA, G, RHO, MXPOT, dspl, .false., s)
    call surface_geometry(mesh, geom, error)
    if (allocated(error)) return
    call compare_all(mesh, geom, s, .false., error)
  end subroutine compare_single

  subroutine compare_multi(depth, isy_in, error)
    real(8), intent(in) :: depth
    integer, intent(in) :: isy_in
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    type(surface_geom) :: geom
    type(first_order_state) :: s
    complex(8) :: dspl(12)

    call set_wave(depth)
    call load_multi(isy_in, mesh)
    call motions(dspl)
    call build_state(mesh, W1, WK, V, H, 1d0, BETA, G, RHO, MXPOT_MULTI_COMB, dspl, .false., s)
    call surface_geometry(mesh, geom, error)
    if (allocated(error)) return
    call compare_all(mesh, geom, s, .true., error)
  end subroutine compare_multi

  !> At every waterline point at each offset, close to the panels, and at one point far from them.
  subroutine compare_all(mesh, geom, s, is_multi, error)
    type(panel_mesh), intent(in) :: mesh
    type(surface_geom), intent(in) :: geom
    type(first_order_state), intent(in) :: s
    logical, intent(in) :: is_multi
    type(error_type), allocatable, intent(out) :: error
    complex(8), allocatable :: rows_g(:,:,:), rows_dg(:,:,:), row_g(:,:), row_dg(:,:), phi_wl(:)
    complex(8) :: existing, field
    real(8) :: x(3), xm(3), length
    integer :: io, iw
    character(len=120) :: info

    allocate (rows_g(mesh%nbase, mesh%nsys, geom%nwl), rows_dg(mesh%nbase, mesh%nsys, geom%nwl), phi_wl(geom%nwl))
    do io = 1, size(offsets)
      call waterline_rows(mesh, geom, offsets(io), is_multi, rows_g, rows_dg)
      call waterline_field(mesh, geom, s, rows_g, rows_dg, phi_wl)
      do iw = 1, geom%nwl
        length = geom%wl_len(iw)
        x = geom%wl_xm(:, iw) - offsets(io)*length*geom%wl_nh(:, iw)
        x(3) = -offsets(io)*length
        xm = [geom%wl_xm(1, iw), geom%wl_xm(2, iw), 0d0]
        existing = existing_potential(x, s, mesh, is_multi)
        field = phi_wl(iw) - incident_phi(xm, s%omega, s%k, s%h, s%amp, s%beta, wave_gravity(s%omega, s%k, s%h))
        write (info, '(a,f6.3,a,i0,a,2es11.3,a,es10.3)') "offset ", offsets(io), " segment ", iw, " existing ", &
          existing, " difference ", abs(field - existing)
        call check(error, abs(field - existing) <= tol*max(abs(existing), 1d0), "Waterline: "//trim(info))
        if (allocated(error)) return
      end do
    end do

    allocate (row_g(mesh%nbase, mesh%nsys), row_dg(mesh%nbase, mesh%nsys))
    call field_influence_row(mesh, far_point, is_multi, row_g, row_dg)
    existing = existing_potential(far_point, s, mesh, is_multi)
    field = field_potential(mesh, s, row_g, row_dg)
    write (info, '(a,2es11.3,a,es10.3)') "existing ", existing, " difference ", abs(field - existing)
    call check(error, abs(field - existing) <= tol*max(abs(existing), 1d0), "Far point: "//trim(info))
  end subroutine compare_all

  !> Perturbation potential from the existing routines: the scattered part with F0 subtracted, as a physical potential (i times
  !> the HAMS value), and the radiated part of each mode with velocity -i omega xi_j.
  complex(8) function existing_potential(x, s, mesh, is_multi)
    real(8), intent(in) :: x(3)
    type(first_order_state), intent(in) :: s
    type(panel_mesh), intent(in) :: mesh
    logical, intent(in) :: is_multi
    complex(8) :: pot
    integer, allocatable :: factor_normal(:)
    integer :: jb, j

    if (is_multi) then
      allocate (factor_normal(NELEM_PE))
      factor_normal = 0
      call CalPotentialMulti(x, 'Diffraction', 6*NBODY + 1, pot, NBODY, 0, factor_normal)
    else
      call CalPotential(x, 'Diffraction', 7, pot)
    end if
    existing_potential = ci*(pot - VINP(x(1), x(2), x(3), XW(1), XW(2), BETA))
    do jb = 1, mesh%nbody
      do j = 1, 6
        if (is_multi) then
          factor_normal = merge(1, 0, mesh%body(1:mesh%nbase) == jb)
          call CalPotentialMulti(x, 'Radiation', j, pot, NBODY, jb, factor_normal)
        else
          call CalPotential(x, 'Radiation', j, pot)
        end if
        existing_potential = existing_potential - ci*s%omega*s%xi(j, jb)*pot
      end do
    end do
  end function existing_potential

  ! ---------------------------------------------------------------------------------------------- set-up

  subroutine set_wave(depth)
    real(8), intent(in) :: depth

    NTHREAD = 1
    ISOL = 1
    IRSP = 0
    AMP = 1d0
    BETA = heading
    XW = 0d0
    W1 = omega
    TP = 2d0*acos(-1d0)/W1
    H = depth
    V = W1**2/G
    if (H < 0d0) then
      WK = V
    else
      call DISPERSION(WVN, NK, W1, H)
      WK = WVN(1)
    end if
  end subroutine set_wave

  subroutine surface_geometry(mesh, geom, error)
    type(panel_mesh), intent(in) :: mesh
    type(surface_geom), intent(out) :: geom
    type(error_type), allocatable, intent(out) :: error
    logical :: ok
    character(len=:), allocatable :: msg

    call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
    call check(error, ok, "build_surface_geom stopped")
    if (allocated(error)) return
    call check(error, geom%nwl > 0, "No waterline")
  end subroutine surface_geometry

  !> The panels of the cylinder on the side of the symmetry plane that is kept.
  subroutine kept_panels(full, isx_in, isy_in, keep)
    type(panel_mesh), intent(in) :: full
    integer, intent(in) :: isx_in, isy_in
    logical, allocatable, intent(out) :: keep(:)

    allocate (keep(full%n))
    keep = .true.
    if (isy_in == 1) keep = full%xc(2, :) > 0d0
    if (isx_in == 1) keep = full%xc(1, :) > 0d0
  end subroutine kept_panels

  subroutine set_symmetry(isx_in, isy_in)
    integer, intent(in) :: isx_in, isy_in

    ISX = isx_in
    ISY = isy_in
    NSYS = 1 + isx_in + isy_in
  end subroutine set_symmetry

  subroutine load_single(isx_in, isy_in, mesh)
    integer, intent(in) :: isx_in, isy_in
    type(panel_mesh), intent(out) :: mesh
    type(panel_mesh) :: full

    call make_cylinder(radius, draft, ntheta, nz, nr, full, .true.)
    call load_panels(full, isx_in, isy_in, mesh)
  end subroutine load_single

  !> Puts a whole hull in the single-body arrays, or the side of it that a symmetry plane keeps.
  subroutine load_panels(full, isx_in, isy_in, mesh)
    type(panel_mesh), intent(in) :: full
    integer, intent(in) :: isx_in, isy_in
    type(panel_mesh), intent(out) :: mesh
    logical, allocatable :: keep(:)
    integer, allocatable :: body(:)
    integer :: i, k, ne

    call set_symmetry(isx_in, isy_in)
    NBODY = 1
    call kept_panels(full, isx_in, isy_in, keep)
    ne = count(keep)
    NELEM = ne
    NELEM_PE = ne
    NTND = full%nnode
    if (allocated(XYZ)) deallocate (XYZ)
    if (allocated(NCN)) deallocate (NCN)
    if (allocated(NCON)) deallocate (NCON)
    if (allocated(PNSZ)) deallocate (PNSZ)
    if (allocated(XYZ_P)) deallocate (XYZ_P)
    if (allocated(DXYZ_P)) deallocate (DXYZ_P)
    if (allocated(GEONRM_P)) deallocate (GEONRM_P)
    if (allocated(DS)) deallocate (DS)
    if (allocated(MXPOT)) deallocate (MXPOT)
    allocate (XYZ(NTND, 3), NCN(ne), NCON(ne, 4), PNSZ(ne), XYZ_P(ne, 3), DXYZ_P(ne, 6), GEONRM_P(ne, 3), DS(ne), MXPOT(ne, 7, NSYS))
    XYZ = transpose(full%xnode)
    NCN = pack(full%ncn, keep)
    PNSZ = pack(full%psize, keep)
    DS = pack(full%area, keep)
    do k = 1, 4
      NCON(:, k) = pack(full%ncon(k, :), keep)
    end do
    do k = 1, 3
      XYZ_P(:, k) = pack(full%xc(k, :), keep)
      GEONRM_P(:, k) = pack(full%nrm(k, :), keep)
    end do
    do k = 1, 6
      DXYZ_P(:, k) = pack(full%nj(k, :), keep)
    end do
    call fill_potentials(MXPOT)
    allocate (body(ne))
    body = 1
    call expand_images(XYZ_P, GEONRM_P, DXYZ_P, DS, PNSZ, XYZ, NCN, NCON, body, reshape([(0d0, i = 1, 3)], [3, 1]), ISX, ISY, NSYS, mesh)
  end subroutine load_panels

  !> Two of the cylinders, in the combined arrays. The panel moment arms are about the centre of each cylinder.
  subroutine load_multi(isy_in, mesh)
    integer, intent(in) :: isy_in
    type(panel_mesh), intent(out) :: mesh
    type(panel_mesh) :: full
    logical, allocatable :: keep(:)
    integer, allocatable :: body(:)
    real(8) :: xr(3, 2)
    integer :: k, b, ne, nn, first, node_first

    call set_symmetry(0, isy_in)
    NBODY = 2
    call make_cylinder(radius, draft, ntheta, nz, nr, full, .true.)
    call kept_panels(full, 0, isy_in, keep)
    ne = count(keep)
    nn = full%nnode
    NELEM_PE = 2*ne
    if (allocated(XYZ_GLOBAL_MULTI_COMB)) deallocate (XYZ_GLOBAL_MULTI_COMB)
    if (allocated(XYZ_GLOBAL_MULTI_COMB_P)) deallocate (XYZ_GLOBAL_MULTI_COMB_P)
    if (allocated(NCN_MULTI_COMB)) deallocate (NCN_MULTI_COMB)
    if (allocated(NCON_MULTI_COMB)) deallocate (NCON_MULTI_COMB)
    if (allocated(PNSZ_MULTI_COMB)) deallocate (PNSZ_MULTI_COMB)
    if (allocated(DS_MULTI_COMB)) deallocate (DS_MULTI_COMB)
    if (allocated(DXYZ_MULTI_COMB)) deallocate (DXYZ_MULTI_COMB)
    if (allocated(GEONRM_MULTI_COMB)) deallocate (GEONRM_MULTI_COMB)
    if (allocated(MXPOT_MULTI_COMB)) deallocate (MXPOT_MULTI_COMB)
    allocate (XYZ_GLOBAL_MULTI_COMB(2*nn, 3), XYZ_GLOBAL_MULTI_COMB_P(2*ne, 3), NCN_MULTI_COMB(2*ne), NCON_MULTI_COMB(2*ne, 4), &
              PNSZ_MULTI_COMB(2*ne), DS_MULTI_COMB(2*ne), DXYZ_MULTI_COMB(2*ne, 6), GEONRM_MULTI_COMB(2*ne, 3), &
              MXPOT_MULTI_COMB(2*ne, 13, NSYS), body(2*ne))
    do b = 1, 2
      first = (b - 1)*ne
      node_first = (b - 1)*nn
      xr(:, b) = [body_spacing*(b - 1), 0d0, 0d0]
      do k = 1, 3
        XYZ_GLOBAL_MULTI_COMB(node_first + 1:node_first + nn, k) = full%xnode(k, :) + xr(k, b)
        XYZ_GLOBAL_MULTI_COMB_P(first + 1:first + ne, k) = pack(full%xc(k, :), keep) + xr(k, b)
        GEONRM_MULTI_COMB(first + 1:first + ne, k) = pack(full%nrm(k, :), keep)
      end do
      do k = 1, 4
        NCON_MULTI_COMB(first + 1:first + ne, k) = pack(full%ncon(k, :), keep) + node_first
      end do
      do k = 1, 6
        DXYZ_MULTI_COMB(first + 1:first + ne, k) = pack(full%nj(k, :), keep)
      end do
      NCN_MULTI_COMB(first + 1:first + ne) = pack(full%ncn, keep)
      PNSZ_MULTI_COMB(first + 1:first + ne) = pack(full%psize, keep)
      DS_MULTI_COMB(first + 1:first + ne) = pack(full%area, keep)
      body(first + 1:first + ne) = b
    end do
    call fill_potentials(MXPOT_MULTI_COMB)
    call expand_images(XYZ_GLOBAL_MULTI_COMB_P, GEONRM_MULTI_COMB, DXYZ_MULTI_COMB, DS_MULTI_COMB, PNSZ_MULTI_COMB, &
                       XYZ_GLOBAL_MULTI_COMB, NCN_MULTI_COMB, NCON_MULTI_COMB, body, xr, ISX, ISY, NSYS, mesh)
  end subroutine load_multi

  !> Reproducible values of the size of unit, different for each panel, mode and image.
  subroutine fill_potentials(mx)
    complex(8), intent(out) :: mx(:,:,:)
    integer :: i, md, is

    do is = 1, size(mx, 3)
      do md = 1, size(mx, 2)
        do i = 1, size(mx, 1)
          mx(i, md, is) = cmplx(sin(1.3d0*i + 0.7d0*md + 0.4d0*is), cos(0.9d0*i - 1.1d0*md + 2.1d0*is), kind=8)
        end do
      end do
    end do
  end subroutine fill_potentials

  subroutine motions(dspl)
    complex(8), intent(out) :: dspl(:)
    integer :: j

    do j = 1, size(dspl)
      dspl(j) = cmplx(cos(0.8d0*j), sin(1.7d0*j), kind=8)
    end do
  end subroutine motions

end module test_waterline_field
