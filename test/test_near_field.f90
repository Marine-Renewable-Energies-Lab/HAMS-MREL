!> Unit tests for the near-field mean drift terms, with exact or synthetic fields only (no BEM solve).
module test_near_field
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use SecondOrderTypes, only : panel_mesh, first_order_state, surface_geom
  use FirstOrderState, only : incident_phi
  use SurfaceGeometry, only : build_surface_geom, DEFAULT_SHARP_DEG
  use FarFieldDrift, only : group_ratio
  use NearFieldDrift, only : pair_minus, waterline_extrapolate, near_field_pair
  use second_order_fixtures, only : make_cylinder, set_panel, maccamy_fuchs_state, maccamy_fuchs_phi, maccamy_fuchs_drift
  implicit none
  private

  public :: collect_near_field

  real(8), parameter :: pi = 3.141592653589793d0
  real(8), parameter :: g = 9.81d0, rho = 1025d0
  real(8), parameter :: radius = 5d0, depth = 10d0, omega = 1d0
  complex(8), parameter :: ci = (0d0, 1d0)

  ! Positions in the term index of near_field_pair
  integer, parameter :: T_I = 1, T_IA = 2, T_II = 3, T_III = 4, T_IV = 5, T_V = 6

  ! Box of hydrostatic_exact_box: footprint [bx1,bx2] x [by1,by2], z down to -bdraft, reference point bxr
  real(8), parameter :: bx1 = -1d0, bx2 = 5d0, by1 = -3d0, by2 = 1d0, bdraft = 3d0
  real(8), parameter :: bxr(3) = [0.5d0, 0.25d0, 0.7d0]

contains

subroutine collect_near_field(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("pair_mean_identity", test_pair_mean_identity), &
    new_unittest("pair_time_average", test_pair_time_average), &
    new_unittest("waterline_extrapolate_linear", test_waterline_extrapolate_linear), &
    new_unittest("fixed_cylinder_analytic", test_fixed_cylinder_analytic), &
    new_unittest("term3_uniform", test_term3_uniform), &
    new_unittest("term4_cross", test_term4_cross), &
    new_unittest("hydrostatic_exact_box", test_hydrostatic_exact_box), &
    new_unittest("flare_only_vertical", test_flare_only_vertical), &
    new_unittest("form_equivalence", test_form_equivalence), &
    new_unittest("pair_symmetry", test_pair_symmetry), &
    new_unittest("pair_superposition", test_pair_superposition), &
    new_unittest("reference_point", test_reference_point), &
    new_unittest("two_bodies_independent", test_two_bodies_independent) &
    ]
end subroutine collect_near_field

! ---------------------------------------------------------------------------------------------- helpers

!> Wave number from omega^2/g = k tanh(kh), by Newton iteration.
real(8) function wave_number(h)
  real(8), intent(in) :: h
  real(8) :: nu, k, fk
  integer :: it
  nu = omega**2 / g
  k = max(nu, sqrt(nu / h))
  do it = 1, 100
    fk = k * tanh(k * h) - nu
    k = k - fk / (tanh(k * h) + k * h / cosh(k * h)**2)
    if (abs(fk) < 1d-15) exit
  end do
  wave_number = k
end function wave_number

!> Reproducible values in [-1, 1], so the tests do not touch the shared random number state under OpenMP.
real(8) function pseudo(i)
  integer, intent(in) :: i
  pseudo = sin(12.9898d0 * real(i, 8) + 0.5d0) * cos(78.233d0 * real(i, 8))
end function pseudo

complex(8) function pseudo_c(i)
  integer, intent(in) :: i
  pseudo_c = cmplx(pseudo(2 * i), pseudo(2 * i + 1), kind=8)
end function pseudo_c

pure function cross_r(a, b) result(c)
  real(8), intent(in) :: a(3), b(3)
  real(8) :: c(3)
  c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
end function cross_r

pure function cross_c(a, b) result(c)
  complex(8), intent(in) :: a(3), b(3)
  complex(8) :: c(3)
  c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
end function cross_c

subroutine geometry(mesh, geom, error)
  type(panel_mesh), intent(in) :: mesh
  type(surface_geom), intent(out) :: geom
  type(error_type), allocatable, intent(out) :: error
  logical :: ok
  character(len=:), allocatable :: msg

  call build_surface_geom(mesh, 1d-6, DEFAULT_SHARP_DEG, geom, ok, msg)
  call check(error, ok, "build_surface_geom stopped: " // msg)
end subroutine geometry

!> State with all fields zero, on the panels of mesh.
subroutine zero_state(mesh, s)
  type(panel_mesh), intent(in) :: mesh
  type(first_order_state), intent(out) :: s

  s%omega = omega; s%k = omega**2 / g; s%nu = s%k; s%h = -1d0; s%amp = 1d0; s%beta = 0d0; s%g = g; s%rho = rho
  allocate(s%phi_i(mesh%n), s%phi_p(mesh%n), s%dphi_p_dn(mesh%n), s%dphi_dn(mesh%n))
  allocate(s%grad_i(3, mesh%n), s%grad_p(3, mesh%n), s%xi(6, mesh%nbody))
  s%phi_i = (0d0, 0d0); s%phi_p = (0d0, 0d0); s%dphi_p_dn = (0d0, 0d0); s%dphi_dn = (0d0, 0d0)
  s%grad_i = (0d0, 0d0); s%grad_p = (0d0, 0d0); s%xi = (0d0, 0d0)
end subroutine zero_state

!> Closed cylinder with analytic fields and made-up motions about xr: the fields are not a solution for this motion,
!> but every term below is an algebraic statement about given fields.
subroutine free_cylinder(xr, beta, mesh, geom, s, phi_wl, error)
  real(8), intent(in) :: xr(3), beta
  type(panel_mesh), intent(out) :: mesh
  type(surface_geom), intent(out) :: geom
  type(first_order_state), intent(out) :: s
  complex(8), allocatable, intent(out) :: phi_wl(:)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: k
  integer :: j

  call make_cylinder(radius, depth, 120, 20, 4, mesh, .true.)
  mesh%xr(:, 1) = xr
  call geometry(mesh, geom, error)
  if (allocated(error)) return
  k = wave_number(depth)
  call maccamy_fuchs_state(mesh, radius, depth, k, omega, 1d0, beta, g, rho, s)
  call free_motion(beta, s%xi(:, 1))
  allocate(phi_wl(geom%nwl))
  do j = 1, geom%nwl
    phi_wl(j) = maccamy_fuchs_phi(geom%wl_xm(:, j), radius, depth, k, omega, 1d0, beta, g)
  end do
end subroutine free_cylinder

subroutine free_motion(beta, xi)
  real(8), intent(in) :: beta
  complex(8), intent(out) :: xi(6)
  xi = exp(ci * beta) * [(0.3d0, 0.1d0), (-0.2d0, 0.4d0), (0.5d0, -0.1d0), (0.02d0, 0.05d0), (-0.04d0, 0.03d0), (0.06d0, -0.02d0)]
end subroutine free_motion

!> Sets panel ip of mesh from four nodes, stored as nodes 4*ip-3 .. 4*ip.
subroutine add_quad(mesh, ip, x1, x2, x3, x4, nrm)
  type(panel_mesh), intent(inout) :: mesh
  integer, intent(in) :: ip
  real(8), intent(in) :: x1(3), x2(3), x3(3), x4(3), nrm(3)
  mesh%xnode(:, 4*ip - 3) = x1
  mesh%xnode(:, 4*ip - 2) = x2
  mesh%xnode(:, 4*ip - 1) = x3
  mesh%xnode(:, 4*ip) = x4
  mesh%ncn(ip) = 4
  mesh%ncon(:, ip) = [4*ip - 3, 4*ip - 2, 4*ip - 1, 4*ip]
  call set_panel(mesh, ip, nrm)
end subroutine add_quad

subroutine alloc_mesh(npan, xr, mesh)
  integer, intent(in) :: npan
  real(8), intent(in) :: xr(3)
  type(panel_mesh), intent(out) :: mesh
  integer :: i

  mesh%n = npan; mesh%nbase = npan; mesh%nnode = 4 * npan; mesh%nbody = 1; mesh%nsys = 1
  allocate(mesh%xc(3,npan), mesh%nrm(3,npan), mesh%area(npan), mesh%psize(npan), mesh%nj(6,npan))
  allocate(mesh%xnode(3,4*npan), mesh%ncn(npan), mesh%ncon(4,npan), mesh%body(npan), mesh%base(npan), mesh%image(npan))
  allocate(mesh%xr(3,1))
  mesh%xr(:, 1) = xr
  mesh%body = 1
  mesh%image = 1
  mesh%base = [(i, i = 1, npan)]
end subroutine alloc_mesh

!> Wall-sided box of the hydrostatic test: open at z = 0, closed below, with normals into the body.
subroutine make_box(xr, mesh)
  real(8), intent(in) :: xr(3)
  type(panel_mesh), intent(out) :: mesh
  integer, parameter :: nx = 4, ny = 3, nz = 2
  integer :: ip, i, j
  real(8) :: xa, xb, ya, yb, za, zb

  call alloc_mesh(2 * nx * nz + 2 * ny * nz + nx * ny, xr, mesh)
  ip = 0
  do j = 1, nz
    za = -bdraft * real(j - 1, 8) / nz
    zb = -bdraft * real(j, 8) / nz
    do i = 1, nx
      xa = bx1 + (bx2 - bx1) * real(i - 1, 8) / nx
      xb = bx1 + (bx2 - bx1) * real(i, 8) / nx
      ip = ip + 1
      call add_quad(mesh, ip, [xa, by1, za], [xb, by1, za], [xb, by1, zb], [xa, by1, zb], [0d0, 1d0, 0d0])
      ip = ip + 1
      call add_quad(mesh, ip, [xa, by2, za], [xb, by2, za], [xb, by2, zb], [xa, by2, zb], [0d0, -1d0, 0d0])
    end do
    do i = 1, ny
      ya = by1 + (by2 - by1) * real(i - 1, 8) / ny
      yb = by1 + (by2 - by1) * real(i, 8) / ny
      ip = ip + 1
      call add_quad(mesh, ip, [bx1, ya, za], [bx1, yb, za], [bx1, yb, zb], [bx1, ya, zb], [1d0, 0d0, 0d0])
      ip = ip + 1
      call add_quad(mesh, ip, [bx2, ya, za], [bx2, yb, za], [bx2, yb, zb], [bx2, ya, zb], [-1d0, 0d0, 0d0])
    end do
  end do
  do j = 1, ny
    ya = by1 + (by2 - by1) * real(j - 1, 8) / ny
    yb = by1 + (by2 - by1) * real(j, 8) / ny
    do i = 1, nx
      xa = bx1 + (bx2 - bx1) * real(i - 1, 8) / nx
      xb = bx1 + (bx2 - bx1) * real(i, 8) / nx
      ip = ip + 1
      call add_quad(mesh, ip, [xa, ya, -bdraft], [xb, ya, -bdraft], [xb, yb, -bdraft], [xa, yb, -bdraft], [0d0, 0d0, 1d0])
    end do
  end do
end subroutine make_box

!> Frustum side wall with r(z) = r0 + slope z, open at both ends, with normals into the body.
subroutine make_frustum(r0, slope, draft, nth, nz, mesh)
  real(8), intent(in) :: r0, slope, draft
  integer, intent(in) :: nth, nz
  type(panel_mesh), intent(out) :: mesh
  integer :: i, j, ip
  real(8) :: t1, t2, tm, z1, z2, r1, r2, den

  call alloc_mesh(nth * nz, [0d0, 0d0, 0d0], mesh)
  den = sqrt(1d0 + slope**2)
  ip = 0
  do j = 1, nz
    z1 = -draft * real(j - 1, 8) / nz
    z2 = -draft * real(j, 8) / nz
    r1 = r0 + slope * z1
    r2 = r0 + slope * z2
    do i = 1, nth
      t1 = 2d0 * pi * real(i - 1, 8) / nth
      t2 = 2d0 * pi * real(i, 8) / nth
      tm = 0.5d0 * (t1 + t2)
      ip = ip + 1
      call add_quad(mesh, ip, [r1*cos(t1), r1*sin(t1), z1], [r1*cos(t2), r1*sin(t2), z1], &
                    [r2*cos(t2), r2*sin(t2), z2], [r2*cos(t1), r2*sin(t1), z2], [-cos(tm)/den, -sin(tm)/den, slope/den])
    end do
  end do
end subroutine make_frustum

!> Rotation of eq. 15.49 as implied by F_S and M_S: R = Rx(a1) Ry(a2) Rz(a3), whose second-order part is lower
!> triangular, with row 3 = (a1 a3, a2 a3, -(a1^2 + a2^2)/2).
pure function rotation(a) result(r)
  real(8), intent(in) :: a(3)
  real(8) :: r(3,3), c(3), s(3)
  c = cos(a)
  s = sin(a)
  r(1,:) = [c(2)*c(3), -c(2)*s(3), s(2)]
  r(2,:) = [c(1)*s(3) + s(1)*s(2)*c(3), c(1)*c(3) - s(1)*s(2)*s(3), -s(1)*c(2)]
  r(3,:) = [s(1)*s(3) - c(1)*s(2)*c(3), s(1)*c(3) + c(1)*s(2)*s(3), c(1)*c(2)]
end function rotation

!> Exact buoyancy load (F_z, M_x, M_y about O) on the box moved to X = R(alpha) x + xi + O, from columns of the box
!> in body axes (x relative to O). Each column runs from the box bottom up to the plane X_3 = 0, and is linear in z,
!> so it integrates exactly. Weight is not part of the load.
pure function exact_box_load(q) result(load)
  real(8), intent(in) :: q(6)
  real(8) :: load(3)
  integer, parameter :: ngrid = 400
  real(8) :: r(3,3), dx, dy, da, x, y, ztop, zbot, h, zm, xi_, yi_
  integer :: i, j

  r = rotation(q(4:6))
  dx = (bx2 - bx1) / ngrid
  dy = (by2 - by1) / ngrid
  da = dx * dy
  zbot = -bdraft - bxr(3)
  load = 0d0
  do j = 1, ngrid
    y = by1 + (real(j, 8) - 0.5d0) * dy - bxr(2)
    do i = 1, ngrid
      x = bx1 + (real(i, 8) - 0.5d0) * dx - bxr(1)
      ztop = -(r(3,1) * x + r(3,2) * y + q(3) + bxr(3)) / r(3,3)
      h = ztop - zbot
      zm = 0.5d0 * (ztop + zbot)
      xi_ = r(1,1) * x + r(1,2) * y + r(1,3) * zm + q(1)
      yi_ = r(2,1) * x + r(2,2) * y + r(2,3) * zm + q(2)
      load = load + [h, yi_ * h, -xi_ * h]
    end do
  end do
  load = load * rho * g * da
end function exact_box_load

!> Mean second-order part of the exact load for in-phase motion of amplitude q: the part of F(q cos t) in q^2 has mean
!> (1/2) [ (F(q) + F(-q))/2 - F(0) ].
pure function exact_box_mean(q, f0) result(mean)
  real(8), intent(in) :: q(6), f0(3)
  real(8) :: mean(3)
  mean = 0.5d0 * (0.5d0 * (exact_box_load(q) + exact_box_load(-q)) - f0)
end function exact_box_mean

!> (F_z, M_x, M_y) of terms V + Ia for real motions q in still water.
subroutine formula_box_mean(mesh, geom, q, mean)
  type(panel_mesh), intent(in) :: mesh
  type(surface_geom), intent(in) :: geom
  real(8), intent(in) :: q(6)
  real(8), intent(out) :: mean(3)
  type(first_order_state) :: s
  complex(8) :: phi_wl(geom%nwl), terms(6, 6, 1)

  call zero_state(mesh, s)
  s%xi(:, 1) = cmplx(q, 0d0, kind=8)
  phi_wl = (0d0, 0d0)
  call near_field_pair(mesh, geom, s, s, phi_wl, phi_wl, terms)
  mean = real(terms(3:5, T_V, 1) + terms(3:5, T_IA, 1))
end subroutine formula_box_mean

! ---------------------------------------------------------------------------------------------- tests

subroutine test_pair_mean_identity(error)
  type(error_type), allocatable, intent(out) :: error
  complex(8) :: f, gg, p
  integer :: i

  do i = 1, 20
    f = pseudo_c(i)
    gg = pseudo_c(i + 50)
    p = pair_minus(f, gg, f, gg)
    call check(error, abs(real(p) - 0.5d0 * real(f * conjg(gg))) <= 1d-15, "pair_minus(f,g,f,g) is not Re(f g*)/2")
    if (allocated(error)) return
    call check(error, abs(aimag(p)) <= 1d-15, "pair_minus(f,g,f,g) is not real")
    if (allocated(error)) return
  end do
end subroutine test_pair_mean_identity

subroutine test_pair_time_average(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: nt = 64
  complex(8) :: f1, f2, g1, g2, a1, a2, ph, ft, gt
  real(8) :: avg, expected, wt
  integer :: icase, it

  do icase = 1, 20
    f1 = pseudo_c(10 * icase + 1)
    f2 = pseudo_c(10 * icase + 2)
    g1 = pseudo_c(10 * icase + 3)
    g2 = pseudo_c(10 * icase + 4)
    a1 = pseudo_c(10 * icase + 5)
    a2 = pseudo_c(10 * icase + 6)
    avg = 0d0
    do it = 1, nt
      wt = 2d0 * pi * real(it - 1, 8) / nt
      ph = exp(-ci * wt)
      ft = (a1 * f1 + a2 * f2) * ph
      gt = (a1 * g1 + a2 * g2) * ph
      avg = avg + real(ft) * real(gt) / nt
    end do
    expected = abs(a1)**2 * real(pair_minus(f1, g1, f1, g1)) + abs(a2)**2 * real(pair_minus(f2, g2, f2, g2)) &
               + 2d0 * real(a1 * conjg(a2) * pair_minus(f1, g1, f2, g2))
    call check(error, abs(avg - expected) <= 1d-10, "time average differs from the pair expansion")
    if (allocated(error)) return
  end do
end subroutine test_pair_time_average

subroutine test_waterline_extrapolate_linear(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s
  complex(8), allocatable :: phi_wl(:)
  complex(8) :: a(3), b, ref
  real(8) :: h, k, g_wave, xz(3)
  integer :: i, j

  call make_cylinder(radius, depth, 24, 4, 1, mesh, .false.)
  call geometry(mesh, geom, error)
  if (allocated(error)) return
  ! k from a slightly different g, as HAMS does: the incident wave must use the gravity of the solver's wave
  h = 8d0
  k = 0.2d0
  g_wave = omega**2 / (k * tanh(k * h))
  call zero_state(mesh, s)
  s%k = k
  s%h = h
  a = [(0.3d0, -0.1d0), (0.2d0, 0.4d0), (-0.5d0, 0.2d0)]
  b = (1.5d0, -0.7d0)
  do i = 1, mesh%n
    s%phi_p(i) = b + sum(a * mesh%xc(:, i))
    s%grad_p(:, i) = a
    s%phi_i(i) = (1d3, 1d3)   ! the incident part must come from z = 0, not from the panel
  end do
  allocate(phi_wl(geom%nwl))
  call waterline_extrapolate(mesh, geom, s, phi_wl)
  do j = 1, geom%nwl
    xz = [geom%wl_xm(1:2, j), 0d0]
    ref = incident_phi(xz, omega, k, h, 1d0, 0d0, g_wave) + b + sum(a * xz)
    call check(error, abs(phi_wl(j) - ref) <= 1d-12 * abs(ref), "extrapolated waterline potential is wrong")
    if (allocated(error)) return
  end do
end subroutine test_waterline_extrapolate_linear

subroutine test_fixed_cylinder_analytic(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s
  complex(8), allocatable :: phi_wl(:)
  complex(8), allocatable :: terms(:,:,:)
  real(8) :: k, fref, total
  integer :: j

  call make_cylinder(radius, depth, 120, 20, 1, mesh, .false.)
  call geometry(mesh, geom, error)
  if (allocated(error)) return
  k = wave_number(depth)
  call maccamy_fuchs_state(mesh, radius, depth, k, omega, 1d0, 0d0, g, rho, s)
  allocate(phi_wl(geom%nwl), terms(6, 6, 1))
  do j = 1, geom%nwl
    phi_wl(j) = maccamy_fuchs_phi(geom%wl_xm(:, j), radius, depth, k, omega, 1d0, 0d0, g)
  end do
  call near_field_pair(mesh, geom, s, s, phi_wl, phi_wl, terms)

  fref = maccamy_fuchs_drift(radius, k, 0.5d0 * rho * g * group_ratio(k, depth), 720)
  total = real(sum(terms(1, :, 1)))
  call check(error, abs(total - fref) < 2d-2 * fref, "near-field F_x differs from the analytic drift by more than 2%")
  if (allocated(error)) return
  call check(error, maxval(abs(terms(:, [T_III, T_IV, T_V], 1))) <= 1d-12 * fref, &
             "terms III, IV and V are not zero on a fixed body")
  if (allocated(error)) return
  call check(error, abs(sum(terms(2, :, 1))) <= 1d-6 * fref, "lateral force is not zero on a symmetric body")
end subroutine test_fixed_cylinder_analytic

subroutine test_term3_uniform(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s
  complex(8), allocatable :: phi_wl(:), terms(:,:,:)
  complex(8) :: v(3), x(3)
  real(8) :: ref
  integer :: i

  call make_cylinder(radius, depth, 60, 4, 4, mesh, .true.)
  call geometry(mesh, geom, error)
  if (allocated(error)) return
  call zero_state(mesh, s)
  v = [(0.4d0, 0.3d0), (-0.2d0, 0.1d0), (0.7d0, -0.5d0)]
  x = [(0.1d0, -0.2d0), (0.3d0, 0.2d0), (-0.6d0, 0.4d0)]
  do i = 1, mesh%n
    s%grad_i(:, i) = 0.25d0 * v
    s%grad_p(:, i) = 0.75d0 * v
  end do
  s%xi(1:3, 1) = x
  allocate(phi_wl(geom%nwl), terms(6, 6, 1))
  phi_wl = (0d0, 0d0)
  call near_field_pair(mesh, geom, s, s, phi_wl, phi_wl, terms)
  ref = -rho * 0.5d0 * real(sum(conjg(x) * (-ci * omega * v))) * geom%awp(1)
  call check(error, abs(real(terms(3, T_III, 1)) - ref) <= 1d-3 * abs(ref), "term III differs from the uniform-field value")
end subroutine test_term3_uniform

subroutine test_term4_cross(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s
  complex(8), allocatable :: phi_wl(:), terms(:,:,:)
  complex(8) :: c(3), gg(3), hh(3), alpha(3), xi(3)
  real(8) :: fref(3), mref(3), n(3), area(3), arm(3)
  integer :: i, m

  call make_box(bxr, mesh)
  call geometry(mesh, geom, error)
  if (allocated(error)) return
  call zero_state(mesh, s)
  ! Constant potential c_m on the one flat face with normal +e_m, so G and H have closed forms
  c = [(0.4d0, 0.3d0), (-0.2d0, 0.1d0), (0.7d0, -0.5d0)]
  area = [(by2 - by1) * bdraft, (bx2 - bx1) * bdraft, (bx2 - bx1) * (by2 - by1)]
  gg = (0d0, 0d0)
  hh = (0d0, 0d0)
  do m = 1, 3
    n = 0d0
    n(m) = 1d0
    arm = 0.5d0 * [bx1 + bx2, by1 + by2, -bdraft] - bxr
    if (m == 1) arm(1) = bx1 - bxr(1)
    if (m == 2) arm(2) = by1 - bxr(2)
    if (m == 3) arm(3) = -bdraft - bxr(3)
    gg(m) = -rho * (-ci * omega * c(m)) * area(m)
    hh = hh - rho * (-ci * omega * c(m)) * area(m) * cmplx(cross_r(arm, n), 0d0, kind=8)
  end do
  do i = 1, mesh%n
    do m = 1, 3
      if (mesh%nrm(m, i) > 0.5d0) s%phi_p(i) = 0.5d0 * c(m)
      if (mesh%nrm(m, i) > 0.5d0) s%phi_i(i) = 0.5d0 * c(m)
    end do
  end do
  alpha = [(0.02d0, 0.05d0), (-0.04d0, 0.03d0), (0.06d0, -0.02d0)]
  xi = [(0.3d0, 0.1d0), (-0.2d0, 0.4d0), (0.5d0, -0.1d0)]
  s%xi(1:3, 1) = xi
  s%xi(4:6, 1) = alpha
  allocate(phi_wl(geom%nwl), terms(6, 6, 1))
  phi_wl = (0d0, 0d0)
  call near_field_pair(mesh, geom, s, s, phi_wl, phi_wl, terms)
  fref = 0.5d0 * real(cross_c(conjg(alpha), gg))
  mref = 0.5d0 * real(cross_c(conjg(alpha), hh)) + 0.5d0 * real(cross_c(conjg(xi), gg))
  call check(error, maxval(abs(terms(1:3, T_IV, 1) - fref)) <= 1d-12 * maxval(abs(fref)), "term IV force is wrong")
  if (allocated(error)) return
  call check(error, maxval(abs(terms(4:6, T_IV, 1) - mref)) <= 1d-12 * maxval(abs(mref)), "term IV moment is wrong")
end subroutine test_term4_cross

subroutine test_hydrostatic_exact_box(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: eps = 1d-3
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  real(8) :: f0(3), ex(3, 6), fo(3, 6), q(6), ex2(3), fo2(3), cross_ex(3), cross_fo(3), floor_f, floor_m, tol(3)
  integer :: a, b

  call make_box(bxr, mesh)
  call geometry(mesh, geom, error)
  if (allocated(error)) return
  f0 = exact_box_load([0d0, 0d0, 0d0, 0d0, 0d0, 0d0])
  ! A small floor, so that a component that is exactly zero is compared with its round-off
  floor_f = 1d-6 * rho * g * eps**2 * geom%awp(1)
  floor_m = floor_f * (bx2 - bx1)

  do a = 1, 6
    q = 0d0
    q(a) = eps
    ex(:, a) = exact_box_mean(q, f0)
    call formula_box_mean(mesh, geom, q, fo(:, a))
    tol = 1d-3 * abs(ex(:, a)) + [floor_f, floor_m, floor_m]
    call check(error, all(abs(fo(:, a) - ex(:, a)) <= tol), "V + Ia differs from the exact box for a single mode")
    if (allocated(error)) return
  end do

  ! Cross terms: move both modes together, in phase, and remove the single-mode parts
  do a = 1, 5
    do b = a + 1, 6
      q = 0d0
      q(a) = eps
      q(b) = eps
      ex2 = exact_box_mean(q, f0)
      call formula_box_mean(mesh, geom, q, fo2)
      cross_ex = ex2 - ex(:, a) - ex(:, b)
      cross_fo = fo2 - fo(:, a) - fo(:, b)
      tol = 1d-3 * abs(cross_ex) + [floor_f, floor_m, floor_m]
      call check(error, all(abs(cross_fo - cross_ex) <= tol), "V + Ia differs from the exact box for a pair of modes")
      if (allocated(error)) return
    end do
  end do
end subroutine test_hydrostatic_exact_box

subroutine test_flare_only_vertical(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: r0 = 3d0, slope = 0.25d0, eps = 0.1d0
  type(panel_mesh) :: mesh, wall
  type(surface_geom) :: geom, gwall
  type(first_order_state) :: s
  complex(8), allocatable :: phi_wl(:), terms(:,:,:)
  real(8) :: ref
  integer :: j

  call make_frustum(r0, slope, 2d0, 120, 3, mesh)
  call geometry(mesh, geom, error)
  if (allocated(error)) return
  call zero_state(mesh, s)
  allocate(phi_wl(geom%nwl), terms(6, 6, 1))

  ! Pure heave in still water: the extra volume is (1/2) dA/dz h^2, so the mean extra force is
  ! rho g (1/2) (2 pi r0 slope) <h^2> = rho g pi r0 slope eps^2 / 2, positive on a hull widening upward.
  s%xi(3, 1) = eps
  phi_wl = (0d0, 0d0)
  call near_field_pair(mesh, geom, s, s, phi_wl, phi_wl, terms)
  ref = rho * g * pi * r0 * slope * eps**2 / 2d0
  call check(error, abs(real(terms(3, T_IA, 1)) - ref) <= 1d-3 * ref, "flare force differs from the heave volume change")
  if (allocated(error)) return

  ! A varying elevation and rotations still give no horizontal force
  s%xi(:, 1) = [(0.1d0, 0d0), (0d0, 0.2d0), (0.1d0, 0.1d0), (0.03d0, 0d0), (0d0, 0.02d0), (0.01d0, 0d0)]
  do j = 1, geom%nwl
    phi_wl(j) = cmplx(cos(0.7d0 * geom%wl_xm(1, j)), sin(0.4d0 * geom%wl_xm(2, j)), kind=8)
  end do
  call near_field_pair(mesh, geom, s, s, phi_wl, phi_wl, terms)
  call check(error, maxval(abs(terms(1:2, T_IA, 1))) == 0d0, "flare term has a horizontal force")
  if (allocated(error)) return
  call check(error, abs(terms(3, T_IA, 1)) > 0d0, "flare term is zero on a flared hull")
  if (allocated(error)) return

  ! Wall-sided hull: no flare term at all
  call make_cylinder(r0, 2d0, 60, 3, 1, wall, .false.)
  call geometry(wall, gwall, error)
  if (allocated(error)) return
  deallocate(phi_wl)
  allocate(phi_wl(gwall%nwl))
  phi_wl = (1d0, 0.5d0)
  call zero_state(wall, s)
  s%xi(3, 1) = eps
  call near_field_pair(wall, gwall, s, s, phi_wl, phi_wl, terms)
  call check(error, maxval(abs(terms(:, T_IA, 1))) == 0d0, "flare term is not zero on a wall-sided hull")
end subroutine test_flare_only_vertical

subroutine test_form_equivalence(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s
  complex(8), allocatable :: phi_wl(:), terms(:,:,:)
  complex(8) :: alt(2), zr, gg(3), f_hs, q(6)
  real(8) :: x(3), fx
  integer :: j, i

  call free_cylinder([1d0, -0.5d0, 0d0], 0d0, mesh, geom, s, phi_wl, error)
  if (allocated(error)) return
  allocate(terms(6, 6, 1))
  call near_field_pair(mesh, geom, s, s, phi_wl, phi_wl, terms)
  q = s%xi(:, 1)

  ! Alternative form: (zeta - X3)^2 on the waterline, and the rotation of dynamic plus first-order hydrostatic force
  alt = (0d0, 0d0)
  do j = 1, geom%nwl
    x = geom%wl_xm(:, j) - mesh%xr(:, 1)
    zr = ci * omega / g * phi_wl(j) - (q(3) + q(4) * x(2) - q(5) * x(1))
    alt = alt + 0.5d0 * rho * g * geom%wl_nh(1:2, j) * geom%wl_len(j) * pair_minus(zr, zr, zr, zr)
  end do
  gg = (0d0, 0d0)
  do i = 1, mesh%n
    gg = gg - rho * mesh%nrm(:, i) * mesh%area(i) * (-ci * omega) * (s%phi_i(i) + s%phi_p(i))
  end do
  f_hs = -rho * g * geom%awp(1) * (q(3) + geom%yf(1) * q(4) - geom%xf(1) * q(5))
  gg(3) = gg(3) + f_hs
  alt(1) = alt(1) + pair_minus(q(5), gg(3), q(5), gg(3)) - pair_minus(q(6), gg(2), q(6), gg(2))
  alt(2) = alt(2) + pair_minus(q(6), gg(1), q(6), gg(1)) - pair_minus(q(4), gg(3), q(4), gg(3))

  fx = abs(terms(1, T_I, 1) + terms(1, T_IV, 1))
  call check(error, fx > 0d0, "degenerate test: no horizontal force")
  if (allocated(error)) return
  call check(error, abs(terms(1, T_I, 1) + terms(1, T_IV, 1) - alt(1)) <= 1d-3 * fx, "the two forms differ in F_x")
  if (allocated(error)) return
  call check(error, abs(terms(2, T_I, 1) + terms(2, T_IV, 1) - alt(2)) <= 1d-3 * fx, "the two forms differ in F_y")
end subroutine test_form_equivalence

subroutine test_pair_symmetry(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s1, s2
  complex(8), allocatable :: w1(:), w2(:), t12(:,:,:), t21(:,:,:), t11(:,:,:)

  call free_cylinder([1d0, -0.5d0, 0.3d0], 0d0, mesh, geom, s1, w1, error)
  if (allocated(error)) return
  call free_cylinder([1d0, -0.5d0, 0.3d0], 0.7d0, mesh, geom, s2, w2, error)
  if (allocated(error)) return
  allocate(t12(6, 6, 1), t21(6, 6, 1), t11(6, 6, 1))
  call near_field_pair(mesh, geom, s1, s2, w1, w2, t12)
  call near_field_pair(mesh, geom, s2, s1, w2, w1, t21)
  call check(error, maxval(abs(t12)) > 0d0 .and. maxval(abs(t12 - conjg(t21))) <= 1d-12 * maxval(abs(t12)), &
             "terms(s1,s2) is not conj(terms(s2,s1))")
  if (allocated(error)) return
  call near_field_pair(mesh, geom, s1, s1, w1, w1, t11)
  call check(error, maxval(abs(aimag(t11))) <= 1d-12 * maxval(abs(t11)), "diagonal terms are not real")
end subroutine test_pair_symmetry

!> Every term is a sum of products of first-order quantities, so the total for a superposed state must expand exactly
!> into the diagonal and pair terms, including IV, whose factors are integrals.
subroutine test_pair_superposition(error)
  type(error_type), allocatable, intent(out) :: error
  complex(8), parameter :: a1 = (0.8d0, -0.3d0), a2 = (0.4d0, 0.9d0)
  type(panel_mesh) :: mesh
  type(surface_geom) :: geom
  type(first_order_state) :: s1, s2, st
  complex(8), allocatable :: w1(:), w2(:), wt(:), t11(:,:,:), t22(:,:,:), t12(:,:,:), tt(:,:,:), ref(:,:,:)

  call free_cylinder([1d0, -0.5d0, 0.3d0], 0d0, mesh, geom, s1, w1, error)
  if (allocated(error)) return
  call free_cylinder([1d0, -0.5d0, 0.3d0], 0.7d0, mesh, geom, s2, w2, error)
  if (allocated(error)) return
  st = s1
  st%phi_i = a1 * s1%phi_i + a2 * s2%phi_i
  st%phi_p = a1 * s1%phi_p + a2 * s2%phi_p
  st%grad_i = a1 * s1%grad_i + a2 * s2%grad_i
  st%grad_p = a1 * s1%grad_p + a2 * s2%grad_p
  st%xi = a1 * s1%xi + a2 * s2%xi
  wt = a1 * w1 + a2 * w2
  allocate(t11(6, 6, 1), t22(6, 6, 1), t12(6, 6, 1), tt(6, 6, 1))
  call near_field_pair(mesh, geom, s1, s1, w1, w1, t11)
  call near_field_pair(mesh, geom, s2, s2, w2, w2, t22)
  call near_field_pair(mesh, geom, s1, s2, w1, w2, t12)
  call near_field_pair(mesh, geom, st, st, wt, wt, tt)
  ref = abs(a1)**2 * t11 + abs(a2)**2 * t22 + 2d0 * real(a1 * conjg(a2) * t12)
  call check(error, maxval(abs(tt - ref)) <= 1d-10 * maxval(abs(tt)), "superposed state differs from the pair expansion")
end subroutine test_pair_superposition

subroutine test_reference_point(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: d(3) = [2d0, -1d0, 0.5d0]
  type(panel_mesh) :: mesh, moved
  type(surface_geom) :: geom, gmoved
  type(first_order_state) :: s, sm
  complex(8), allocatable :: w(:), t(:,:,:), tm(:,:,:)
  complex(8) :: f(3), mexpected(3)
  real(8) :: scale_m
  integer :: it

  call free_cylinder([1d0, -0.5d0, 0.3d0], 0.4d0, mesh, geom, s, w, error)
  if (allocated(error)) return
  allocate(t(6, 6, 1), tm(6, 6, 1))
  call near_field_pair(mesh, geom, s, s, w, w, t)

  ! The same rigid motion about O' = O + d: xi' = xi + alpha x d, alpha' = alpha.
  moved = mesh
  moved%xr(:, 1) = mesh%xr(:, 1) + d
  call geometry(moved, gmoved, error)
  if (allocated(error)) return
  sm = s
  sm%xi(1:3, 1) = s%xi(1:3, 1) + cross_c(s%xi(4:6, 1), cmplx(d, 0d0, kind=8))
  call near_field_pair(moved, gmoved, sm, sm, w, w, tm)

  f = sum(t(1:3, :, 1), dim=2)
  scale_m = maxval(abs(f)) * norm2(d)
  call check(error, abs(f(1)) + abs(f(2)) > 0d0 .and. &
             abs(sum(tm(6, :, 1)) - (sum(t(6, :, 1)) - (d(1) * f(2) - d(2) * f(1)))) <= 1d-6 * scale_m, &
             "M_z did not change by -(d x F)_z under a shift of the reference point")
  if (allocated(error)) return

  ! Terms I to IV keep all components exactly; V depends on O at second order, so only its M_z is covered above.
  do it = T_I, T_IV
    f = t(1:3, it, 1)
    mexpected = t(4:6, it, 1) - cross_c(cmplx(d, 0d0, kind=8), f)
    call check(error, maxval(abs(tm(1:3, it, 1) - f)) <= 1d-6 * maxval(abs(t(1:3, :, 1))), &
               "a force term changed under a shift of the reference point")
    if (allocated(error)) return
    call check(error, maxval(abs(tm(4:6, it, 1) - mexpected)) <= 1d-6 * maxval(abs(t(4:6, :, 1))), &
               "a moment term did not change by -d x F under a shift of the reference point")
    if (allocated(error)) return
  end do
end subroutine test_reference_point

subroutine test_two_bodies_independent(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: shift(3) = [40d0, 7d0, 0d0]
  type(panel_mesh) :: mesh, both
  type(surface_geom) :: geom, gboth
  type(first_order_state) :: s1, s2, st
  complex(8), allocatable :: w1(:), w2(:), wt(:), t1(:,:,:), t2(:,:,:), tt(:,:,:)
  integer :: n, nn

  ! Each single-body load, with its own motion
  call free_cylinder([1d0, -0.5d0, 0.3d0], 0d0, mesh, geom, s1, w1, error)
  if (allocated(error)) return
  s2 = s1
  w2 = w1
  s2%xi(:, 1) = (0.5d0, -0.25d0) * s1%xi(:, 1)
  allocate(t1(6, 6, 1), t2(6, 6, 1), tt(6, 6, 2))
  call near_field_pair(mesh, geom, s1, s1, w1, w1, t1)
  call near_field_pair(mesh, geom, s2, s2, w2, w2, t2)

  ! The same two bodies in one mesh, the second a translated copy with the same fields panel by panel
  n = mesh%n
  nn = mesh%nnode
  both = mesh
  both%n = 2 * n; both%nbase = 2 * n; both%nnode = 2 * nn; both%nbody = 2
  both%xc = reshape([mesh%xc, mesh%xc + spread(shift, 2, n)], [3, 2 * n])
  both%nrm = reshape([mesh%nrm, mesh%nrm], [3, 2 * n])
  both%nj = reshape([mesh%nj, mesh%nj], [6, 2 * n])
  both%area = [mesh%area, mesh%area]
  both%psize = [mesh%psize, mesh%psize]
  both%xnode = reshape([mesh%xnode, mesh%xnode + spread(shift, 2, nn)], [3, 2 * nn])
  both%ncn = [mesh%ncn, mesh%ncn]
  both%ncon = reshape([mesh%ncon, merge(mesh%ncon + nn, 0, mesh%ncon > 0)], [4, 2 * n])
  both%body = [spread(1, 1, n), spread(2, 1, n)]
  both%base = [mesh%base, mesh%base + n]
  both%image = [mesh%image, mesh%image]
  both%xr = reshape([mesh%xr(:, 1), mesh%xr(:, 1) + shift], [3, 2])
  call geometry(both, gboth, error)
  if (allocated(error)) return

  st = s1
  st%phi_i = [s1%phi_i, s1%phi_i]
  st%phi_p = [s1%phi_p, s1%phi_p]
  st%grad_i = reshape([s1%grad_i, s1%grad_i], [3, 2 * n])
  st%grad_p = reshape([s1%grad_p, s1%grad_p], [3, 2 * n])
  st%xi = reshape([s1%xi(:, 1), s2%xi(:, 1)], [6, 2])
  wt = [w1, w2]
  call near_field_pair(both, gboth, st, st, wt, wt, tt)
  call check(error, maxval(abs(tt(:, :, 1) - t1(:, :, 1))) <= 1d-8 * maxval(abs(t1)), "body 1 differs from its single-body load")
  if (allocated(error)) return
  call check(error, maxval(abs(tt(:, :, 2) - t2(:, :, 1))) <= 1d-8 * maxval(abs(t2)), "body 2 differs from its single-body load")
end subroutine test_two_bodies_independent

end module test_near_field
