!> Unit tests of QTFField: second-order pressure and elevation of a pair of waves at a point (quadratic part plus the
!> second-order incident potential).
module test_qtf_field
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use QTFField, only : field_pair, build_field_points, field_direct, qfield_init, qfield_keep, qfield_keep_direct, &
                       qfield_direct_active, qfield_finish
  use Const_mod, only : G_SOLVER => G, RHO_SOLVER => RHO
  use IO, only : CreateDirectory
  use SecondOrderTypes, only : panel_mesh, cs_points, so_options, first_order_state, cs_fields
  use second_order_fixtures, only : make_cylinder, make_two_bodies, set_panel
  use SecondOrderIncident, only : incident2, incident2_pair, incident2_phi
  use FirstOrderState, only : incident_phi, incident_grad
  use NearFieldDrift, only : KERNEL_DIFF, KERNEL_SUM
  implicit none
  private

  public :: collect_qtf_field, collect_qtf_field_serial

  real(8), parameter :: g = 9.81d0, rho = 1025d0, pi = 3.14159265358979323846d0, deg = pi / 180d0
  complex(8), parameter :: ci = (0d0, 1d0)

contains

subroutine collect_qtf_field(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("field_stokes_elevation", test_stokes_elevation), &
    new_unittest("field_stokes_pressure", test_stokes_pressure), &
    new_unittest("field_mean_pressure", test_mean_pressure), &
    new_unittest("field_swap_symmetry", test_swap_symmetry), &
    new_unittest("field_bichromatic_time", test_bichromatic_time), &
    new_unittest("field_sum_continuity", test_sum_continuity), &
    new_unittest("points_classified", test_points_classified), &
    new_unittest("points_near_zero_is_surface", test_points_near_zero), &
    new_unittest("points_above_rejected", test_points_above), &
    new_unittest("points_on_hull_rejected", test_points_on_hull), &
    new_unittest("points_inside_rejected", test_points_inside), &
    new_unittest("points_inside_waterplane_rejected", test_points_waterplane), &
    new_unittest("points_two_bodies", test_points_two_bodies), &
    new_unittest("points_mixed_panel_sizes", test_points_mixed_sizes), &
    new_unittest("field_direct_pressure_elevation", test_direct_pressure_elevation) &
    ]
end subroutine collect_qtf_field

!> The tests that run the store of QTFField (module state) and write test folders: serial, one suite
subroutine collect_qtf_field_serial(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("field_direct_file", test_direct_file), &
    new_unittest("field_direct_sum_rows", test_direct_sum_rows), &
    new_unittest("field_direct_sum_symmetric", test_direct_sum_symmetric), &
    new_unittest("field_direct_sum_cutoff", test_direct_sum_cutoff), &
    new_unittest("field_direct_sum_only", test_direct_sum_only), &
    new_unittest("field_direct_both_kernels", test_direct_both_kernels), &
    new_unittest("field_direct_off", test_direct_off), &
    new_unittest("field_header_qtf1_c716ed7", test_header_qtf1) &
    ]
end subroutine collect_qtf_field_serial

!> Wave number of omega, from omega^2 = g k tanh(k h) (h < 0: deep water).
real(8) function wave_number(omega, h)
  real(8), intent(in) :: omega, h
  real(8) :: k, nu
  integer :: it
  nu = omega**2 / g
  k = nu
  if (h < 0d0) then
    wave_number = k
    return
  end if
  k = max(nu, sqrt(nu / h))
  do it = 1, 100
    k = k - (k * tanh(k * h) - nu) / (tanh(k * h) + k * h / cosh(k * h)**2)
  end do
  wave_number = k
end function wave_number

!> Second-order incident potential of the pair (omega, heading in degrees) of unit waves in depth h, under the kernel.
subroutine make_pair(w1, b1, w2, b2, h, inc, kernel)
  real(8), intent(in) :: w1, b1, w2, b2, h
  type(incident2), intent(out) :: inc
  integer, intent(in) :: kernel
  call incident2_pair(w1, wave_number(w1, h), b1 * deg, w2, wave_number(w2, h), b2 * deg, h, g, inc, kernel)
end subroutine make_pair

!> First-order potential and gradient of the unit wave (omega, heading in degrees) at x.
subroutine wave_at(w, b, h, x, phi, u)
  real(8), intent(in) :: w, b, h, x(3)
  complex(8), intent(out) :: phi, u(3)
  real(8) :: k
  k = wave_number(w, h)
  phi = incident_phi(x, w, k, h, 1d0, b * deg, g)
  call incident_grad(x, w, k, h, 1d0, b * deg, g, u)
end subroutine wave_at

!> Quadratic and incident parts of the pair (w1, b1), (w2, b2) of unit waves at x.
subroutine pair_total(kernel, w1, b1, w2, b2, h, x, surface, quad, pot)
  integer, intent(in) :: kernel
  real(8), intent(in) :: w1, b1, w2, b2, h, x(3)
  logical, intent(in) :: surface
  complex(8), intent(out) :: quad, pot
  type(incident2) :: inc
  complex(8) :: phi1, u1(3), phi2, u2(3)
  call wave_at(w1, b1, h, x, phi1, u1)
  call wave_at(w2, b2, h, x, phi2, u2)
  call make_pair(w1, b1, w2, b2, h, inc, kernel)
  call field_pair(kernel, w1, phi1, u1, w2, phi2, u2, inc, x, surface, rho, g, quad, pot)
end subroutine pair_total

!> Second-order Stokes elevation (Dean and Dalrymple 1991), finite depth and deep water.
subroutine test_stokes_elevation(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w = 0.9d0, b = 20d0
  real(8), parameter :: x(3) = [0.7d0, 0.3d0, 0d0]
  complex(8) :: quad, pot, expected, ph
  real(8) :: k, hs
  integer :: ih

  do ih = 1, 2
    hs = merge(12d0, -1d0, ih == 1)
    k = wave_number(w, hs)
    call pair_total(KERNEL_SUM, w, b, w, b, hs, x, .true., quad, pot)
    ph = exp(2d0 * ci * k * (x(1) * cos(b * deg) + x(2) * sin(b * deg)))
    if (ih == 1) then
      expected = 0.25d0 * k * cosh(k * hs) * (2d0 + cosh(2d0 * k * hs)) / sinh(k * hs)**3 * ph
    else
      expected = 0.5d0 * k * ph
    end if
    call check(error, abs(quad + pot - expected) <= 1d-10 * abs(expected), "second-order Stokes elevation")
    if (allocated(error)) return
  end do
end subroutine test_stokes_elevation

!> Second-order pressure of a Stokes wave from Bernoulli with the Stokes potential (Dean and Dalrymple 1991).
subroutine test_stokes_pressure(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w = 0.9d0, b = 20d0, h = 12d0
  real(8), parameter :: x(3) = [0.7d0, 0.3d0, -3d0]
  complex(8) :: quad, pot, expected
  real(8) :: k

  k = wave_number(w, h)
  call pair_total(KERNEL_SUM, w, b, w, b, h, x, .false., quad, pot)
  expected = rho * (0.75d0 * w**2 * cosh(2d0 * k * (x(3) + h)) / sinh(k * h)**4 - g * k / (2d0 * sinh(2d0 * k * h))) &
             * exp(2d0 * ci * k * (x(1) * cos(b * deg) + x(2) * sin(b * deg)))
  call check(error, abs(quad + pot - expected) <= 1d-10 * abs(expected), "second-order Stokes pressure")
end subroutine test_stokes_pressure

subroutine test_mean_pressure(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w = 0.9d0, h = 12d0
  real(8), parameter :: x(3) = [0.7d0, 0.3d0, -3d0]
  type(incident2) :: inc
  complex(8) :: phi1, u1(3), quad, pot, expected
  real(8) :: k

  ! random complex first-order data, the same wave twice
  phi1 = (0.3d0, -1.1d0)
  u1 = [(0.2d0, 0.9d0), (-1.3d0, 0.1d0), (0.5d0, -0.6d0)]
  call make_pair(0.8d0, 0d0, 0.8d0, 0d0, h, inc, KERNEL_DIFF)
  call field_pair(KERNEL_DIFF, 0.8d0, phi1, u1, 0.8d0, phi1, u1, inc, x, .false., rho, g, quad, pot)
  expected = -0.25d0 * rho * sum(abs(u1)**2)
  call check(error, abs(quad - expected) <= 1d-14 * abs(expected), "mean pressure of one wave")
  if (allocated(error)) return
  ! Exact by construction: the difference-kernel diagonal sets coef = 0 and Omega = 0
  call check(error, abs(pot) == 0d0, "potential part on the diagonal is not zero")
  if (allocated(error)) return

  k = wave_number(w, h)
  call pair_total(KERNEL_DIFF, w, 0d0, w, 0d0, h, x, .false., quad, pot)
  expected = -rho * g * k * cosh(2d0 * k * (x(3) + h)) / (2d0 * sinh(2d0 * k * h))
  call check(error, abs(quad - expected) <= 1d-12 * abs(expected), "mean pressure of the incident wave")
  if (allocated(error)) return
  ! Exact by construction: the difference-kernel diagonal sets coef = 0 and Omega = 0
  call check(error, abs(pot) == 0d0, "incident mean pressure has a potential part")
end subroutine test_mean_pressure

subroutine test_swap_symmetry(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w1 = 0.8d0, w2 = 1.1d0, h = 12d0
  real(8), parameter :: x(3) = [0.4d0, -0.9d0, -2d0]
  type(incident2) :: a, b
  complex(8) :: phi1, u1(3), phi2, u2(3), q12, p12, q21, p21
  integer :: kernel, is
  logical :: surf

  phi1 = (0.3d0, -1.1d0); phi2 = (-0.7d0, 0.4d0)
  u1 = [(0.2d0, 0.9d0), (-1.3d0, 0.1d0), (0.5d0, -0.6d0)]
  u2 = [(0.8d0, -0.2d0), (0.1d0, 0.7d0), (-0.9d0, -0.4d0)]
  do kernel = KERNEL_DIFF, KERNEL_SUM
    call make_pair(w1, 0d0, w2, 30d0, h, a, kernel)
    call make_pair(w2, 30d0, w1, 0d0, h, b, kernel)
    do is = 1, 2
      surf = is == 1
      call field_pair(kernel, w1, phi1, u1, w2, phi2, u2, a, x, surf, rho, g, q12, p12)
      call field_pair(kernel, w2, phi2, u2, w1, phi1, u1, b, x, surf, rho, g, q21, p21)
      if (kernel == KERNEL_DIFF) then
        q12 = conjg(q12)
        p12 = conjg(p12)
      end if
      call check(error, abs(q21 - q12) <= 1d-13 * abs(q12), "swap symmetry of the quadratic part")
      if (allocated(error)) return
      call check(error, abs(p21 - p12) <= 1d-13 * abs(p12), "swap symmetry of the potential part")
      if (allocated(error)) return
    end do
  end do
end subroutine test_swap_symmetry

!> Time-domain Bernoulli expressions from the first-order waves and the second-order potential (all four ordered pairs,
!> both kernels), sampled over the common period. The e^{-i Omega t} amplitudes are 2 a1 a2 (quad + pot) of the kernels.
subroutine test_bichromatic_time(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: nt = 64
  real(8), parameter :: amp(2) = [1d0, 0.7d0], ws(2) = [1d0, 0.8d0], bs(2) = [0d0, 30d0], h = 12d0
  real(8), parameter :: xs(3, 2) = reshape([0.7d0, 0.3d0, 0d0, 0.7d0, 0.3d0, -2d0], [3, 2])
  type(incident2) :: inc
  complex(8) :: phi(2), u(3, 2), e(2), c(2, 2, 2), comp(2), quad, pot, expected
  real(8) :: x(3), t, uu(3), phi2t, eta1, wzt, val, om(2), period
  integer :: ip, it, i, j, m, kn, ik
  logical :: surf

  om = [ws(1) - ws(2), ws(1) + ws(2)]
  period = 2d0 * pi / (ws(1) - ws(2))
  do ip = 1, 2
    x = xs(:, ip)
    surf = ip == 1
    do i = 1, 2
      call wave_at(ws(i), bs(i), h, x, phi(i), u(:, i))
    end do
    do kn = KERNEL_DIFF, KERNEL_SUM
      do i = 1, 2
        do j = 1, 2
          call make_pair(ws(i), bs(i), ws(j), bs(j), h, inc, kn)
          c(i, j, kn) = incident2_phi(inc, x)
        end do
      end do
    end do
    comp = (0d0, 0d0)
    do it = 0, nt - 1
      t = period * real(it, 8) / real(nt, 8)
      do i = 1, 2
        e(i) = amp(i) * exp(-ci * ws(i) * t)
      end do
      do m = 1, 3
        uu(m) = real(sum(e * u(m, :)))
      end do
      phi2t = 0d0
      do i = 1, 2
        do j = 1, 2
          phi2t = phi2t + real(-ci * (ws(i) - ws(j)) * amp(i) * amp(j) * c(i, j, KERNEL_DIFF) &
                               * exp(-ci * (ws(i) - ws(j)) * t))
          phi2t = phi2t + real(-ci * (ws(i) + ws(j)) * amp(i) * amp(j) * c(i, j, KERNEL_SUM) &
                               * exp(-ci * (ws(i) + ws(j)) * t))
        end do
      end do
      if (surf) then
        eta1 = real(sum(ci * ws / g * e * phi))
        wzt = real(sum(-ci * ws * e * u(3, :)))
        val = -(phi2t + 0.5d0 * sum(uu**2) + eta1 * wzt) / g
      else
        val = -rho * (phi2t + 0.5d0 * sum(uu**2))
      end if
      do ik = 1, 2
        comp(ik) = comp(ik) + val * exp(ci * om(ik) * t)
      end do
    end do
    ! the field holds Re[C e^{-i Omega t}], so the projection on e^{+i Omega t} gives C / 2
    comp = 2d0 * comp / real(nt, 8)
    do kn = KERNEL_DIFF, KERNEL_SUM
      call pair_total(kn, ws(1), bs(1), ws(2), bs(2), h, x, surf, quad, pot)
      expected = 2d0 * amp(1) * amp(2) * (quad + pot)
      call check(error, abs(comp(kn) - expected) <= 1d-10 * abs(expected), &
                 "time-sampled second-order field differs from the pair kernel")
      if (allocated(error)) return
    end do
  end do
end subroutine test_bichromatic_time

subroutine test_sum_continuity(error)
  type(error_type), allocatable, intent(out) :: error
  real(8), parameter :: w = 0.9d0, b = 20d0, h = 12d0
  real(8), parameter :: x(3) = [0.7d0, 0.3d0, 0d0]
  complex(8) :: quad, pot, expected
  real(8) :: k

  k = wave_number(w, h)
  call pair_total(KERNEL_SUM, w, b, w + 1d-5, b, h, x, .true., quad, pot)
  expected = 0.25d0 * k * cosh(k * h) * (2d0 + cosh(2d0 * k * h)) / sinh(k * h)**3 &
             * exp(2d0 * ci * k * (x(1) * cos(b * deg) + x(2) * sin(b * deg)))
  call check(error, abs(quad + pot - expected) <= 1d-4 * abs(expected), "sum-frequency pair does not tend to the Stokes value")
end subroutine test_sum_continuity

!> Cylinder of radius 5, draft 5, with bottom.
subroutine one_cylinder(mesh)
  type(panel_mesh), intent(out) :: mesh
  call make_cylinder(5d0, 5d0, 24, 10, 5, mesh, .true.)
end subroutine one_cylinder

subroutine test_points_classified(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  call one_cylinder(mesh)
  call build_field_points(mesh, reshape([15d0, 0d0, 0d0, 0d0, -20d0, -3d0], [3, 2]), pts, surface, ok, msg)
  call check(error, ok, "classified points rejected: " // msg)
  if (allocated(error)) return
  call check(error, surface(1) .and. .not. surface(2), "surface flags")
  if (allocated(error)) return
  call check(error, all(pts%psize > 0d0) .and. all(pts%body == 1), "psize or body")
end subroutine test_points_classified

subroutine test_points_near_zero(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  call one_cylinder(mesh)
  call build_field_points(mesh, reshape([15d0, 0d0, -1d-7], [3, 1]), pts, surface, ok, msg)
  call check(error, ok .and. surface(1) .and. pts%x(3, 1) == 0d0, "z = -1e-7 is not a surface point at z = 0")
end subroutine test_points_near_zero

subroutine test_points_above(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  call one_cylinder(mesh)
  call build_field_points(mesh, reshape([15d0, 0d0, 0.1d0], [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "above the free surface") > 0, "above the free surface not rejected")
end subroutine test_points_above

subroutine test_points_on_hull(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  integer :: i
  call one_cylinder(mesh)
  do i = 1, mesh%n
    if (abs(mesh%nrm(3, i)) < 0.1d0) exit
  end do
  call build_field_points(mesh, reshape(mesh%xc(:, i), [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "hull") > 0, "point on the hull not rejected")
end subroutine test_points_on_hull

subroutine test_points_inside(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  call one_cylinder(mesh)
  call build_field_points(mesh, reshape([0d0, 0d0, -2d0], [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "inside body 1") > 0, "point inside not rejected")
end subroutine test_points_inside

subroutine test_points_waterplane(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  call one_cylinder(mesh)
  call build_field_points(mesh, reshape([0d0, 0d0, 0d0], [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "inside") > 0, "centre of the waterplane not rejected")
  if (allocated(error)) return
  call build_field_points(mesh, reshape([2d0, 1d0, 0d0], [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "inside") > 0, "waterplane point not rejected")
end subroutine test_points_waterplane

subroutine test_points_two_bodies(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh, m1
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  call one_cylinder(m1)
  call make_two_bodies(m1, 40d0, mesh)
  call build_field_points(mesh, reshape([40d0, 0d0, -2d0], [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "body 2") > 0, "point inside body 2 not rejected")
  if (allocated(error)) return
  call build_field_points(mesh, reshape([20d0, 0d0, 0d0], [3, 1]), pts, surface, ok, msg)
  call check(error, ok .and. pts%body(1) >= 1, "point between the bodies rejected")
end subroutine test_points_two_bodies

!> Box 20 x 10 x 5 below z = 0 (open top): side walls of 0.5 m panels, bottom of 5 x 3 panels of 4 m x 3.33 m.
subroutine mixed_box(mesh)
  type(panel_mesh), intent(out) :: mesh
  integer :: ip, inode
  mesh%n = 2*200 + 2*400 + 15; mesh%nbase = mesh%n; mesh%nnode = 4*mesh%n; mesh%nbody = 1
  allocate(mesh%xc(3, mesh%n), mesh%nrm(3, mesh%n), mesh%area(mesh%n), mesh%psize(mesh%n), mesh%nj(6, mesh%n))
  allocate(mesh%xnode(3, mesh%nnode), mesh%ncn(mesh%n), mesh%ncon(4, mesh%n), mesh%body(mesh%n), mesh%base(mesh%n), &
           mesh%image(mesh%n), mesh%xr(3, 1))
  mesh%xr = 0d0; mesh%body = 1; mesh%image = 1; mesh%ncn = 4
  ip = 0; inode = 0
  call add_face([0d0, 0d0, -5d0], [0d0, 10d0, 0d0], [0d0, 0d0, 5d0], 20, 10, [1d0, 0d0, 0d0])
  call add_face([20d0, 0d0, -5d0], [0d0, 10d0, 0d0], [0d0, 0d0, 5d0], 20, 10, [-1d0, 0d0, 0d0])
  call add_face([0d0, 0d0, -5d0], [20d0, 0d0, 0d0], [0d0, 0d0, 5d0], 40, 10, [0d0, 1d0, 0d0])
  call add_face([0d0, 10d0, -5d0], [20d0, 0d0, 0d0], [0d0, 0d0, 5d0], 40, 10, [0d0, -1d0, 0d0])
  call add_face([0d0, 0d0, -5d0], [20d0, 0d0, 0d0], [0d0, 10d0, 0d0], 5, 3, [0d0, 0d0, 1d0])
  mesh%base = [(ip, ip = 1, mesh%n)]
contains
  subroutine add_face(o, u, v, nu, nv, nin)
    real(8), intent(in) :: o(3), u(3), v(3), nin(3)
    integer, intent(in) :: nu, nv
    integer :: i, j
    do i = 0, nu - 1
      do j = 0, nv - 1
        ip = ip + 1
        mesh%xnode(:, inode + 1) = o + u*i/nu + v*j/nv
        mesh%xnode(:, inode + 2) = o + u*(i + 1)/nu + v*j/nv
        mesh%xnode(:, inode + 3) = o + u*(i + 1)/nu + v*(j + 1)/nv
        mesh%xnode(:, inode + 4) = o + u*i/nu + v*(j + 1)/nv
        mesh%ncon(:, ip) = inode + [1, 2, 3, 4]
        inode = inode + 4
        call set_panel(mesh, ip, nin)
      end do
    end do
  end subroutine add_face
end subroutine mixed_box

!> Below the box, 3.5 m under the bottom and 0.3 m inside the plane of a side wall: the nearest centroid is on the fine
!> wall, whose plane puts the point inside. The point is outside and must be accepted; points in the box and in its
!> waterplane must still be rejected.
subroutine test_points_mixed_sizes(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: mesh
  type(cs_points) :: pts
  logical, allocatable :: surface(:)
  logical :: ok
  character(len=:), allocatable :: msg
  call mixed_box(mesh)
  call build_field_points(mesh, reshape([0.3d0, 5d0, -8.5d0, 25d0, 5d0, 0d0], [3, 2]), pts, surface, ok, msg)
  call check(error, ok, "point below the box rejected: " // msg)
  if (allocated(error)) return
  call build_field_points(mesh, reshape([10d0, 5d0, -0.5d0], [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "inside body 1") > 0, "point in the box not rejected")
  if (allocated(error)) return
  call build_field_points(mesh, reshape([10d0, 5d0, 0d0], [3, 1]), pts, surface, ok, msg)
  call check(error, .not. ok .and. index(msg, "inside body 1") > 0, "waterplane point not rejected")
end subroutine test_points_mixed_sizes

!> The direct part of the pressure and of the elevation of a pair, from phi_D of the pair, is i rho Omega phi_D and
!> (i Omega / g) phi_D, with Omega = omega_1 - omega_2 signed: the same forms as the second-order incident part of field_pair
!> (pot = i rho Omega phi_2, i Omega/g phi_2), so with phi_D replaced by phi_2 the two are equal, for i /= j in both orders
!> (Omega of either sign) and for i = j (Omega = 0). Per unit amplitudes; the amplitude 2 a_i a_j* (i /= j) or |a_i|^2
!> (i = j) multiplies the row afterwards, as for quad and pot.
subroutine test_direct_pressure_elevation(error)
  type(error_type), allocatable, intent(out) :: error
  complex(8), parameter :: phi = (0.7d0, -0.2d0)
  complex(8), parameter :: unit_grad(3) = [(1d0, 0d0), (0d0, 0d0), (0d0, 0d0)]
  real(8), parameter :: x(3) = [2d0, -1d0, -3d0], xs(3) = [2d0, -1d0, 0d0]
  type(incident2) :: inc
  complex(8) :: p2, quad, pot, d
  real(8) :: w(2), om
  integer :: i, j, ks

  ! The values by hand: Omega = 0.3, phi = 0.7 - 0.2 i, i phi = 0.2 + 0.7 i
  d = field_direct(0.3d0, phi, rho, g, .false.)
  call check(error, abs(d - rho*0.3d0*(0.2d0, 0.7d0)) < 1d-9*rho, "pressure part is not i rho Omega phi_D")
  if (allocated(error)) return
  d = field_direct(0.3d0, phi, rho, g, .true.)
  call check(error, abs(d - 0.3d0/g*(0.2d0, 0.7d0)) < 1d-12, "elevation part is not (i Omega/g) phi_D")
  if (allocated(error)) return
  d = field_direct(-0.3d0, conjg(phi), rho, g, .false.)
  call check(error, abs(d - conjg(field_direct(0.3d0, phi, rho, g, .false.))) < 1d-9*rho, "pair (j,i) is not the conjugate")
  if (allocated(error)) return
  call check(error, abs(field_direct(0d0, phi, rho, g, .false.)) == 0d0, "Omega = 0 does not give 0")
  if (allocated(error)) return

  ! The forms of field_pair, with phi_D = phi_2, for the pairs of (0.5, 0.7) in both orders and for i = j
  w = [0.5d0, 0.7d0]
  do ks = 1, 2
    do i = 1, 2
      do j = 1, 2
        call incident2_pair(w(i), wave_number(w(i), -1d0), 0d0, w(j), wave_number(w(j), -1d0), 20d0*deg, -1d0, g, inc, &
                            KERNEL_DIFF)
        p2 = incident2_phi(inc, merge(xs, x, ks == 2))
        call field_pair(KERNEL_DIFF, w(i), (1d0, 0d0), unit_grad, w(j), (1d0, 0d0), unit_grad, inc, &
                        merge(xs, x, ks == 2), ks == 2, rho, g, quad, pot)
        om = w(i) - w(j)
        d = field_direct(om, p2, rho, g, ks == 2)
        call check(error, abs(d - pot) <= 1d-12*abs(pot), "direct form differs from the incident part with phi_D = phi_2")
        if (allocated(error)) return
        call check(error, abs(inc%omega - om) < 1d-14, "Omega of the pair is not omega_i - omega_j")
        if (allocated(error)) return
      end do
    end do
  end do
end subroutine test_direct_pressure_elevation

!> Runs the field store for two frequencies (0.5, 0.7), one heading and two points (pressure at z = -1, elevation at z = 0),
!> qtf = 3, with the first-order fields set by hand. With direct, phi_D of the pair (1, 2) is phi12(ip) and that of (2, 1) its
!> conjugate, as direct_at_omega gives them (Omega < 0 conjugated); the pairs (1, 1), (2, 2) are not given (Omega = 0).
subroutine run_field_store(root, direct, phi12, phisum, sumkept, qtf_kind, both)
  character(len=*), intent(in) :: root
  logical, intent(in) :: direct
  complex(8), intent(in) :: phi12(2)
  complex(8), intent(in), optional :: phisum(2, 3)   ! sum phi_D at the 2 points of the pairs (1, 1), (1, 2), (2, 2)
  logical, intent(in), optional :: sumkept(3)         ! which of the three sum pairs the direct solve gives
  integer, intent(in), optional :: qtf_kind           ! opts%qtf (default 3)
  logical, intent(in), optional :: both               ! with phisum: also keep the difference pairs
  type(so_options) :: opts
  type(cs_points) :: pts
  type(first_order_state) :: st(1)
  type(cs_fields) :: fl(1)
  logical :: ok, surface(2)
  integer :: i
  real(8), parameter :: w(2) = [0.5d0, 0.7d0]

  call remove_dir(root)
  call CreateDirectory(root // "/Wamit_format", ok)
  opts%qtf = 3
  if (present(qtf_kind)) opts%qtf = qtf_kind
  opts%qtf_field = 1
  opts%qtf_direct = merge(1, 0, direct)
  opts%mean_drift = 1
  pts%n = 2
  allocate (pts%x(3, 2), pts%nrm(3, 2), pts%w(2), pts%psize(2), pts%kind(2), pts%body(2))
  pts%x = reshape([2d0, -1d0, -1d0, 2d0, -1d0, 0d0], [3, 2])
  surface = [.false., .true.]
  call qfield_init(opts, 2, [0d0], pts, surface, root)
  do i = 1, 2
    st(1)%omega = w(i)
    st(1)%k = w(i)**2 / G_SOLVER
    st(1)%h = -1d0
    st(1)%beta = 0d0
    st(1)%rho = RHO_SOLVER
    allocate (fl(1)%phi(2), fl(1)%grad(3, 2))
    fl(1)%phi = [cmplx(0.3d0, 0.1d0*i, 8), cmplx(-0.2d0, 0.4d0, 8)]
    fl(1)%grad = reshape([cmplx(0.1d0*i, 0d0, 8), cmplx(0d0, 0.2d0, 8), cmplx(0.3d0, 0.1d0, 8), &
                          cmplx(0.2d0, 0d0, 8), cmplx(0d0, 0.1d0*i, 8), cmplx(0.1d0, 0.2d0, 8)], [3, 2])
    call qfield_keep(st, fl, 0.1d0*i)
    deallocate (fl(1)%phi, fl(1)%grad)
  end do
  if (direct .and. (.not. present(phisum) .or. present(both))) then
    call qfield_keep_direct(1, 1, 2, 1, phi12)
    call qfield_keep_direct(2, 1, 1, 1, conjg(phi12))
  end if
  if (present(phisum)) then
    ! One call per symmetric pair (w1 >= w2 order of the driver): the store fills both slots
    if (sumkept(1)) call qfield_keep_direct(1, 1, 1, 1, phisum(:, 1), KERNEL_SUM)
    if (sumkept(2)) call qfield_keep_direct(2, 1, 1, 1, phisum(:, 2), KERNEL_SUM)
    if (sumkept(3)) call qfield_keep_direct(2, 1, 2, 1, phisum(:, 3), KERNEL_SUM)
  end if
  call qfield_finish(1d0)
end subroutine run_field_store

!> Removes a test folder (Windows, as test_qtf), if it is there
subroutine remove_dir(root)
  character(len=*), intent(in) :: root
  integer :: estat
  call execute_command_line("if exist " // root // " rmdir /s /q " // root, exitstat=estat)
end subroutine remove_dir

!> Reads the data rows of a field file into a character array
subroutine read_rows(path, rows, n)
  character(len=*), intent(in) :: path
  character(len=600), intent(out) :: rows(:)
  integer, intent(out) :: n
  character(len=600) :: line
  integer :: unit, ios

  n = 0
  open (newunit=unit, file=path, status='old', action='read')
  do
    read (unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    if (line(1:1) == '#' .or. len_trim(line) == 0) cycle
    n = n + 1
    rows(n) = line
  end do
  close (unit)
end subroutine read_rows

!> QTF_Field.txt of the store with direct parts: the columns of QTF_FieldPartial.txt first (unchanged), then Re_dir, Im_dir and
!> has_direct. A solved pair (1, 2) or (2, 1) has the part of field_direct (i rho Omega phi_D, (i Omega/g) phi_D, Omega signed,
!> conjugated and scaled as the other columns) and has_direct = 1, the pair (2, 1) with the conjugate phi_D of the pair (1, 2);
!> the mean rows (i = j) and the sum rows have 0 and has_direct = 0.
subroutine test_direct_file(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "qtf_field_direct_tmp"
  complex(8), parameter :: phi12(2) = [(0.4d0, -0.3d0), (-0.1d0, 0.6d0)]
  real(8), parameter :: w(2) = [0.5d0, 0.7d0]
  character(len=600) :: rows(32), prows(32)
  real(8) :: dv(2)
  complex(8) :: e
  integer :: n, np, r, ih, kk, ii, jj, pp

  call run_field_store(root, .true., phi12)
  call check(error, .not. qfield_direct_active(), "The store is still active after qfield_finish")
  if (allocated(error)) return
  call read_rows(root // "/Wamit_format/QTF_Field.txt", rows, n)
  call read_rows(root // "/Wamit_format/QTF_FieldPartial.txt", prows, np)
  call check(error, n == 16 .and. np == 16, "2 kernels, 4 frequency pairs, 2 points are 16 rows in both files")
  if (allocated(error)) return
  ! Order: kernel (1 difference, 2 sum), heading pair, frequency pair (i, j), point
  r = 0
  do kk = 1, 2
    do ii = 1, 2
      do jj = 1, 2
        do pp = 1, 2
          r = r + 1
          call check(error, rows(r)(1:len_trim(prows(r))) == trim(prows(r)), &
                     "The 5b columns of QTF_Field.txt are not those of QTF_FieldPartial.txt")
          if (allocated(error)) return
          read (rows(r)(len_trim(prows(r)) + 1:), '(2ES14.6,I6)') dv, ih
          if (kk == 1 .and. ii /= jj) then
            ! Unit amplitudes, e^{-i Omega t}: i rho Omega phi_D or (i Omega/g) phi_D, then conjugated, scaled as the others
            e = merge(phi12(pp), conjg(phi12(pp)), ii == 1)
            if (pp == 1) then
              e = (0d0, 1d0) * RHO_SOLVER * (w(ii) - w(jj)) * e
              e = conjg(e) / (RHO_SOLVER * G_SOLVER)
            else
              e = (0d0, 1d0) * (w(ii) - w(jj)) / G_SOLVER * e
              e = conjg(e)
            end if
            call check(error, ih == 1, "has_direct of a solved pair is not 1")
            if (allocated(error)) return
            call check(error, abs(cmplx(dv(1), dv(2), 8) - e) < 1d-5 * abs(e), "direct part of a solved pair")
            if (allocated(error)) return
          else
            call check(error, ih == 0 .and. dv(1) == 0d0 .and. dv(2) == 0d0, "mean and sum rows have 0 and has_direct = 0")
            if (allocated(error)) return
          end if
        end do
      end do
    end do
  end do
  call remove_dir(root)
end subroutine test_direct_file

!> The sum rows (K = 2) of QTF_Field.txt with stored sum phi_D (Phase 7c-2): pressure and elevation of the pairs (1, 1), (1, 2),
!> (2, 1), (2, 2) are i rho Omega phi_D and (i Omega/g) phi_D with Omega = w_i + w_j, conjugated and scaled as the other columns,
!> has_direct = 1; the difference rows (nothing kept) have 0 and has_direct = 0. i = j and i /= j. The file has 7
!> digits, so the rows are compared to 1e-5; the forms themselves, at the sum Omega, to 1e-14 below.
subroutine test_direct_sum_rows(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "qtf_field_dsum_tmp"
  complex(8), parameter :: phi12(2) = [(0.4d0, -0.3d0), (-0.1d0, 0.6d0)]
  complex(8), parameter :: ps(2, 3) = reshape([(0.2d0, 0.5d0), (-0.7d0, 0.1d0), (0.9d0, -0.4d0), (0.3d0, 0.3d0), &
                                               (-0.5d0, -0.2d0), (0.1d0, 0.8d0)], [2, 3])
  real(8), parameter :: w(2) = [0.5d0, 0.7d0]
  character(len=600) :: rows(32), prows(32)
  real(8) :: dv(2), om
  complex(8) :: e
  integer :: n, np, r, ih, kk, ii, jj, pp, idx

  ! The forms at the sum Omega, by hand: i rho Omega phi_D and (i Omega/g) phi_D, i /= j (0.5 + 0.7) and i = j (2 x 0.5)
  e = field_direct(w(1) + w(2), (0.2d0, 0.5d0), RHO_SOLVER, G_SOLVER, .false.)
  call check(error, abs(e - RHO_SOLVER*1.2d0*(-0.5d0, 0.2d0)) < 1d-14*RHO_SOLVER, "sum pressure form, i /= j")
  if (allocated(error)) return
  e = field_direct(2d0*w(1), (0.2d0, 0.5d0), RHO_SOLVER, G_SOLVER, .true.)
  call check(error, abs(e - 1d0/G_SOLVER*(-0.5d0, 0.2d0)) < 1d-14, "sum elevation form, i = j")
  if (allocated(error)) return

  call run_field_store(root, .true., phi12, ps, [.true., .true., .true.])
  call read_rows(root // "/Wamit_format/QTF_Field.txt", rows, n)
  call read_rows(root // "/Wamit_format/QTF_FieldPartial.txt", prows, np)
  call check(error, n == 16 .and. np == 16, "16 rows in both files")
  if (allocated(error)) return
  r = 0
  do kk = 1, 2
    do ii = 1, 2
      do jj = 1, 2
        do pp = 1, 2
          r = r + 1
          call check(error, rows(r)(1:len_trim(prows(r))) == trim(prows(r)), "5b columns changed by the sum direct part")
          if (allocated(error)) return
          read (rows(r)(len_trim(prows(r)) + 1:), '(2ES14.6,I6)') dv, ih
          if (kk == 2) then
            idx = ii + jj - 1          ! pairs (1, 1), (1, 2) = (2, 1), (2, 2) are columns 1, 2, 3
            om = w(ii) + w(jj)
            if (pp == 1) then
              e = (0d0, 1d0) * RHO_SOLVER * om * ps(pp, idx)
              e = conjg(e) / (RHO_SOLVER * G_SOLVER)
            else
              e = (0d0, 1d0) * om / G_SOLVER * ps(pp, idx)
              e = conjg(e)
            end if
            call check(error, ih == 1, "has_direct of a kept sum pair is not 1")
            if (allocated(error)) return
            call check(error, abs(cmplx(dv(1), dv(2), 8) - e) < 1d-5 * abs(e), "direct part of a sum row")
            if (allocated(error)) return
          else
            call check(error, ih == 0 .and. dv(1) == 0d0 .and. dv(2) == 0d0, "difference rows have no direct part here")
            if (allocated(error)) return
          end if
        end do
      end do
    end do
  end do
  call remove_dir(root)
end subroutine test_direct_sum_rows

!> The pair (i, j) and the pair (j, i) of the sum carry the same direct columns, to the character (R6: one call fills both
!> slots): rows (1, 2) and (2, 1) at both points.
subroutine test_direct_sum_symmetric(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "qtf_field_dsym_tmp"
  complex(8), parameter :: phi12(2) = [(0.4d0, -0.3d0), (-0.1d0, 0.6d0)]
  complex(8), parameter :: ps(2, 3) = reshape([(0.2d0, 0.5d0), (-0.7d0, 0.1d0), (0.9d0, -0.4d0), (0.3d0, 0.3d0), &
                                               (-0.5d0, -0.2d0), (0.1d0, 0.8d0)], [2, 3])
  character(len=600) :: rows(32), prows(32), a, b
  integer :: n, np, pp

  call run_field_store(root, .true., phi12, ps, [.true., .true., .true.])
  call read_rows(root // "/Wamit_format/QTF_Field.txt", rows, n)
  call read_rows(root // "/Wamit_format/QTF_FieldPartial.txt", prows, np)
  ! Kernel 2 starts at row 9; the order is (1, 1), (1, 2), (2, 1), (2, 2) with 2 points each: (1, 2) at rows 11-12, (2, 1) at 13-14
  do pp = 1, 2
    a = rows(10 + pp)(len_trim(prows(10 + pp)) + 1:)
    b = rows(12 + pp)(len_trim(prows(12 + pp)) + 1:)
    call check(error, trim(a) == trim(b) .and. len_trim(a) > 0, "sum rows (1, 2) and (2, 1) differ in the direct columns")
    if (allocated(error)) return
    call check(error, rows(10 + pp)(1:14) /= rows(12 + pp)(1:14), "rows (1, 2) and (2, 1) are the same row")
    if (allocated(error)) return
  end do
  call remove_dir(root)
end subroutine test_direct_sum_symmetric

!> A sum row the direct solve did not keep (cut by qtf_max_sw or the irregular-frequency cutoff, so never given to
!> qfield_keep_direct) has has_direct = 0 and zero direct columns, as have the rows of the pairs that were kept elsewhere: only
!> (1, 1) and (1, 2) = (2, 1) kept, so (2, 2) is cut.
subroutine test_direct_sum_cutoff(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "qtf_field_dcut_tmp"
  complex(8), parameter :: phi12(2) = [(0.4d0, -0.3d0), (-0.1d0, 0.6d0)]
  complex(8), parameter :: ps(2, 3) = reshape([(0.2d0, 0.5d0), (-0.7d0, 0.1d0), (0.9d0, -0.4d0), (0.3d0, 0.3d0), &
                                               (-0.5d0, -0.2d0), (0.1d0, 0.8d0)], [2, 3])
  character(len=600) :: rows(32), prows(32)
  real(8) :: dv(2)
  integer :: n, np, r, ih, k, expect

  call run_field_store(root, .true., phi12, ps, [.true., .true., .false.])
  call read_rows(root // "/Wamit_format/QTF_Field.txt", rows, n)
  call read_rows(root // "/Wamit_format/QTF_FieldPartial.txt", prows, np)
  do r = 9, 16
    read (rows(r)(len_trim(prows(r)) + 1:), '(2ES14.6,I6)') dv, ih
    k = (r - 9)/2 + 1               ! the pair of the row: 1 = (1, 1), 2 = (1, 2), 3 = (2, 1), 4 = (2, 2)
    expect = merge(0, 1, k == 4)
    call check(error, ih == expect, "has_direct of a sum row is not 1 for a kept pair and 0 for the cut one")
    if (allocated(error)) return
    if (expect == 0) then
      call check(error, dv(1) == 0d0 .and. dv(2) == 0d0, "a cut sum row has direct columns")
      if (allocated(error)) return
    else
      call check(error, abs(dv(1)) + abs(dv(2)) > 0d0, "a kept sum row has zero direct columns")
      if (allocated(error)) return
    end if
  end do
  call remove_dir(root)
end subroutine test_direct_sum_cutoff

!> qtf = 2 (sum only, no difference store): qfield_keep_direct with the sum kernel must work and the file has the 8 sum rows
!> with has_direct = 1 for the kept pairs. Fix round 1: the bounds check once read the unallocated difference store here.
subroutine test_direct_sum_only(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "qtf_field_dsonly_tmp"
  complex(8), parameter :: phi12(2) = [(0.4d0, -0.3d0), (-0.1d0, 0.6d0)]
  complex(8), parameter :: ps(2, 3) = reshape([(0.2d0, 0.5d0), (-0.7d0, 0.1d0), (0.9d0, -0.4d0), (0.3d0, 0.3d0), &
                                               (-0.5d0, -0.2d0), (0.1d0, 0.8d0)], [2, 3])
  character(len=600) :: rows(32), prows(32)
  real(8) :: dv(2)
  integer :: n, np, r, ih

  call run_field_store(root, .true., phi12, ps, [.true., .true., .true.], qtf_kind=2)
  call read_rows(root // "/Wamit_format/QTF_Field.txt", rows, n)
  call read_rows(root // "/Wamit_format/QTF_FieldPartial.txt", prows, np)
  call check(error, n == 8 .and. np == 8, "qtf = 2 has 8 sum rows")
  if (allocated(error)) return
  do r = 1, 8
    read (rows(r)(len_trim(prows(r)) + 1:), '(2ES14.6,I6)') dv, ih
    call check(error, ih == 1 .and. abs(dv(1)) + abs(dv(2)) > 0d0, "sum-only row without direct part")
    if (allocated(error)) return
  end do
  call remove_dir(root)
end subroutine test_direct_sum_only

!> qtf = 3 with both kernels kept in one run: difference rows (i /= j) and all sum rows have has_direct = 1, mean rows of the
!> difference 0, and the two stores do not mix (the difference value is not the sum value).
subroutine test_direct_both_kernels(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "qtf_field_dboth_tmp"
  complex(8), parameter :: phi12(2) = [(0.4d0, -0.3d0), (-0.1d0, 0.6d0)]
  complex(8), parameter :: ps(2, 3) = reshape([(0.2d0, 0.5d0), (-0.7d0, 0.1d0), (0.9d0, -0.4d0), (0.3d0, 0.3d0), &
                                               (-0.5d0, -0.2d0), (0.1d0, 0.8d0)], [2, 3])
  character(len=600) :: rows(32), prows(32)
  real(8) :: dv(2)
  integer :: n, np, r, ih, k, ii, jj
  logical :: expect

  call run_field_store(root, .true., phi12, ps, [.true., .true., .true.], both=.true.)
  call read_rows(root // "/Wamit_format/QTF_Field.txt", rows, n)
  call read_rows(root // "/Wamit_format/QTF_FieldPartial.txt", prows, np)
  call check(error, n == 16, "16 rows")
  if (allocated(error)) return
  do r = 1, 16
    read (rows(r)(len_trim(prows(r)) + 1:), '(2ES14.6,I6)') dv, ih
    k = (r - 1)/8 + 1
    ii = mod((r - 1)/4, 2) + 1
    jj = mod((r - 1)/2, 2) + 1
    expect = k == 2 .or. ii /= jj
    call check(error, (ih == 1) .eqv. expect, "has_direct with both kernels kept")
    if (allocated(error)) return
  end do
  call remove_dir(root)
end subroutine test_direct_both_kernels

!> qtf_field = 1 with qtf_direct = 0: the store is the one of 5b. QTF_Field.txt is not written, and QTF_FieldPartial.txt is
!> the file of the run with direct parts, to the byte (the direct parts are in QTF_Field.txt only), whose header names no
!> direct part. The integration tests qtf_field_* pin that file against the Phase 5b rows.
subroutine test_direct_off(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "qtf_field_off_tmp", root2 = "qtf_field_direct_tmp"
  complex(8), parameter :: phi12(2) = [(0.4d0, -0.3d0), (-0.1d0, 0.6d0)]
  character(len=600) :: a(32), b(32), line
  logical :: exists
  integer :: n, m, r, unit, ios

  call run_field_store(root, .false., phi12)
  call run_field_store(root2, .true., phi12)
  inquire (file=root // "/Wamit_format/QTF_Field.txt", exist=exists)
  call check(error, .not. exists, "qtf_direct = 0 wrote QTF_Field.txt")
  if (allocated(error)) return
  call read_rows(root // "/Wamit_format/QTF_FieldPartial.txt", a, n)
  call read_rows(root2 // "/Wamit_format/QTF_FieldPartial.txt", b, m)
  call check(error, n == m .and. n == 16, "rows of QTF_FieldPartial.txt")
  if (allocated(error)) return
  do r = 1, n
    call check(error, a(r) == b(r), "QTF_FieldPartial.txt differs with direct parts on")
    if (allocated(error)) return
  end do
  ! The header lines too: none names a direct part
  open (newunit=unit, file=root // "/Wamit_format/QTF_FieldPartial.txt", status='old', action='read')
  do
    read (unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    if (index(line, "direct") > 0 .or. index(line, "Re_dir") > 0) then
      call check(error, .false., "The header of QTF_FieldPartial.txt names a direct part")
      exit
    end if
  end do
  close (unit)
  call remove_dir(root)
  call remove_dir(root2)
end subroutine test_direct_off

!> The three header lines of QTF_Field.txt that Phase 7c-2 reworded: with qtf = 1 (no sum column) they are the text of
!> c716ed7, byte for byte; with the sum on (qtf = 3) the sum wording.
subroutine test_header_qtf1(error)
  type(error_type), allocatable, intent(out) :: error
  character(len=*), parameter :: root = "qtf_field_hdr1_tmp", root3 = "qtf_field_hdr3_tmp"
  complex(8), parameter :: phi12(2) = [(0.4d0, -0.3d0), (-0.1d0, 0.6d0)]
  character(len=700) :: l1, l2, l3, m1, m2, m3

  call run_field_store(root, .true., phi12, qtf_kind=1)
  call header_lines(root // "/Wamit_format/QTF_Field.txt", l1, l2, l3)
  call check(error, trim(l1) == "# Re_dir Im_dir: the direct part, i rho Omega phi_D (pressure) and (i Omega/g) phi_D (elevation),"// &
             " Omega = omega_1 - omega_2, per unit amplitude pair, conjugated and scaled as Re Im. Add it to Re Im.", &
             "qtf = 1: Re_dir line is not the one of c716ed7")
  if (allocated(error)) return
  call check(error, trim(l2) == "# phi_D of the direct solve with body forcing; free-surface forcing not included (it is in .12d with"// &
             " qtf_terms = full); radiated wave of the second-order motion not included (Phase 6c)", &
             "qtf = 1: phi_D line is not the one of c716ed7")
  if (allocated(error)) return
  call check(error, trim(l3) == "# has_direct: 1 where the direct solve covers the row (difference frequency, omega_1 /= omega_2),"// &
             " else 0 with Re_dir = Im_dir = 0 (mean rows and sum rows).", "qtf = 1: has_direct line is not the one of c716ed7")
  if (allocated(error)) return
  call remove_dir(root)
  call run_field_store(root3, .true., phi12, qtf_kind=3)
  call header_lines(root3 // "/Wamit_format/QTF_Field.txt", m1, m2, m3)
  call check(error, index(m1, "(K = 1) or omega_1 + omega_2 (K = 2)") > 0 .and. index(m2, "for both kernels") > 0 .and. &
             index(m3, "K = 2: the sum pairs of QTF_Direct.12s") > 0, "qtf = 3: the header does not name the sum")
  call remove_dir(root3)
end subroutine test_header_qtf1

!> The header lines of a field file that start with "# Re_dir Im_dir:", "# phi_D" and "# has_direct:"
subroutine header_lines(path, l1, l2, l3)
  character(len=*), intent(in) :: path
  character(len=*), intent(out) :: l1, l2, l3
  character(len=700) :: line
  integer :: unit, ios

  l1 = ""
  l2 = ""
  l3 = ""
  open (newunit=unit, file=path, status='old', action='read')
  do
    read (unit, '(a)', iostat=ios) line
    if (ios /= 0) exit
    if (line(1:15) == "# Re_dir Im_dir") l1 = line
    if (line(1:8) == "# phi_D ") l2 = line
    if (line(1:13) == "# has_direct:") l3 = line
  end do
  close (unit)
end subroutine header_lines

end module test_qtf_field
