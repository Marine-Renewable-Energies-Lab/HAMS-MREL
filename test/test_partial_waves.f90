!> Unit tests of PartialWaves (Phase 6a, Task 5): the partial-wave decomposition of a first-order field on a vertical
!> cylinder r = R, and its reconstruction on z = 0 outside R. pw_hams_field sets module globals (the mesh arrays and the wave
!> of WaveDyn_mod), so the suite runs serially (see unit_tester.f90).
module test_partial_waves
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use HAMS_mod, only : NTHREAD
  use WaveDyn_mod, only : BETA, H, NK, W1, WK, WVN
  use SecondOrderTypes, only : panel_mesh, first_order_state, cs_points, cs_fields, CS_WALL, CS_LINE, CS_FREE
  use FirstOrderState, only : wave_gravity
  use FieldVelocity, only : fields_green
  use PartialWaves, only : pw_set, pw_sample_points, pw_decompose, pw_eval_surface, pw_bessel_k, pw_evanescent_from_wvn
  use second_order_fixtures, only : make_cylinder, maccamy_fuchs_state
  use test_waterline_field, only : set_wave, load_panels
  implicit none
  private

  public :: collect_partial_waves

  real(8), parameter :: pi = 3.14159265358979323846d0
  complex(8), parameter :: ci = (0d0, 1d0)
  ! The g of DISPERSION (WavDynSubs.f90), which fills WVN
  real(8), parameter :: g_dispersion = 9.807d0

  ! K_m(x), m = 0..10 (rows), x = 0.1, 1, 10 (columns). Generated with mpmath 1.2.1 at 30 digits:
  !   python -I -c "import mpmath as mp; mp.mp.dps=30; [print([mp.nstr(mp.besselk(m, mp.mpf(x)), 17) for m in range(11)])
  !                 for x in ['0.1', '1', '10']]"
  real(8), parameter :: k_table(0:10, 3) = reshape([ &
    2.4270690247020166d0, 9.8538447808706061d0, 1.9950396464211414d+2, 7.9900124304654362d+3, 4.7960024979256828d+5, &
    3.8376009995835928d+7, 3.8380805998333854d+9, 4.6060804799000208d+11, 6.4488964799200125d+13, 1.031869497592001d+16, &
    1.857429584630401d+18, &
    4.2102443824070833d-1, 6.0190723019723457d-1, 1.6248388986351775d0, 7.1012628247379445d0, 4.4232415847062845d+1, &
    3.609605896012407d+2, 3.6538383118594699d+3, 4.4207020331914879d+4, 6.2255212295866777d+5, 1.0005040987670599d+7, &
    1.8071328990102945d+8, &
    1.7780062316167652d-5, 1.8648773453825585d-5, 2.1509817006932769d-5, 2.7252700256598692d-5, 3.7861437160891984d-5, &
    5.7541849985312279d-5, 9.5403287146204263d-5, 1.720257945607574d-4, 3.3623939953126462d-4, 7.1000883381078078d-4, &
    1.61425530039067d-3], [11, 3])
  real(8), parameter :: k_args(3) = [0.1d0, 1d0, 10d0]

  ! pw_hams_field: the 448-panel sharp cylinder (radius 10, draft 10, sharp bottom edge) of the QTF fixtures, in depth 15
  ! (finite depth, so that the field has evanescent modes; 5 m of water under the keel). Modes chosen so the gate passes.
  real(8), parameter :: hf_depth = 15d0, hf_radius = 15d0
  integer, parameter :: hf_mmax = 16, hf_nev = 8
  ! The gate of the brief, and a limit for the horizontal gradient fixed before the measurement
  real(8), parameter :: hf_phi_limit = 1d-4, hf_grad_limit = 1d-3
  ! With the near-panel integrals of fields_green (see test_hams_field): measured 4.9e-4, the Task 6 warning level
  real(8), parameter :: hf_stock_limit = 1d-3

contains

subroutine collect_partial_waves(testsuite)
  type(unittest_type), allocatable, intent(out) :: testsuite(:)

  testsuite = [ &
    new_unittest("pw_bessel_k_values", test_bessel_k_values), &
    new_unittest("pw_round_trip", test_round_trip), &
    new_unittest("pw_hams_field", test_hams_field) &
    ]
end subroutine collect_partial_waves

! ---------------------------------------------------------------------------------------------- K_m

!> K_m(x) for m = 0..10 at x = 0.1, 1 and 10 against mpmath, 1e-12 relative.
subroutine test_bessel_k_values(error)
  type(error_type), allocatable, intent(out) :: error
  real(8) :: kv(0:10), worst
  integer :: ix, m
  character(len=120) :: info

  worst = 0d0
  do ix = 1, 3
    call pw_bessel_k(10, k_args(ix), kv)
    do m = 0, 10
      worst = max(worst, abs(kv(m) - k_table(m, ix)) / k_table(m, ix))
    end do
  end do
  print '(4x,a,es10.2)', "K_m against mpmath, worst relative error ", worst
  write (info, '(a,es10.3)') "K_m against mpmath, worst relative error ", worst
  call check(error, worst <= 1d-12, trim(info))
end subroutine test_bessel_k_values

! ---------------------------------------------------------------------------------------------- round trip

!> H_m^(1)(x) for m = 0..mmax + 1, by the intrinsics.
subroutine hankel_all(mmax, x, hk)
  integer, intent(in) :: mmax
  real(8), intent(in) :: x
  complex(8), intent(out) :: hk(0:mmax + 1)
  integer :: m
  do m = 0, mmax + 1
    hk(m) = cmplx(bessel_jn(m, x), bessel_yn(m, x), 8)
  end do
end subroutine hankel_all

!> A field of known modes, sampled on the cylinder and decomposed, gives back its coefficients to 1e-12; and pw_eval_surface
!> reproduces phi, grad phi, phi_z and phi_zz of the same modes at r = 1.5 R on z = 0 to 1e-12. The reference gradient uses
!> other derivative identities (H_m' = H_{m-1} - m H_m / x, K_m' = -K_{m-1} - m K_m / x) than the module.
subroutine test_round_trip(error)
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: mmax = 6, nev = 4, npt = 7
  real(8), parameter :: h = 20d0, radius = 12d0, omega = 0.9d0, centre(2) = [3d0, -2d0]
  real(8) :: wvn_local(nev + 1), k0, kev(nev), nu, zf, r, th, xk, worst_root, worst_a, x(3)
  real(8) :: kvals(0:mmax + 1), kr(0:mmax + 1, nev)
  complex(8) :: b(-mmax:mmax, 0:nev), a_true(-mmax:mmax, 0:nev), hk(0:mmax + 1), hr(0:mmax + 1)
  complex(8), allocatable :: samples(:)
  complex(8) :: phi, grad(3), phizz, phi_ref, gr_ref(3), phizz_ref, dr, dth, e, term, dterm, sgn
  real(8) :: err_phi, err_grad, err_zz, max_phi, max_grad, max_zz
  type(cs_points) :: pts
  type(pw_set) :: pw
  integer :: m, n, i, ip, am
  character(len=160) :: info
  external :: DISPERSION

  call DISPERSION(wvn_local, nev + 1, omega, h)
  k0 = wvn_local(1)
  kev = wvn_local(2:)
  nu = omega**2 / g_dispersion
  ! The roots of WVN satisfy the dispersion relations, so they can be used as the vertical modes
  worst_root = abs(k0 * tanh(k0 * h) - nu)
  do n = 1, nev
    worst_root = max(worst_root, abs(kev(n) * tan(kev(n) * h) + nu))
  end do
  write (info, '(a,es10.3)') "WVN roots off the dispersion relation, relative to omega^2/g ", worst_root / nu
  call check(error, worst_root <= 1d-12 * nu, trim(info))
  if (allocated(error)) return

  ! Coefficients at r = R, O(1), and the stored coefficients a that multiply H_m(k0 r) and K_m(kn r)
  call hankel_all(mmax, k0 * radius, hk)
  do n = 0, nev
    if (n > 0) call pw_bessel_k(mmax + 1, kev(n) * radius, kr(:, n))
    do m = -mmax, mmax
      b(m, n) = cmplx(sin(1.3d0 * (m + 7 * n) + 0.2d0), cos(0.7d0 * m - 1.1d0 * n), 8)
      am = abs(m)
      if (n == 0) then
        sgn = merge(-1d0, 1d0, m < 0 .and. mod(am, 2) == 1)
        a_true(m, n) = b(m, n) / (sgn * hk(am))
      else
        a_true(m, n) = b(m, n) / kr(am, n)
      end if
    end do
  end do

  ! Samples: phi(R, theta, z) = sum b_mn e^{i m theta} Z_n(z)
  call pw_sample_points(centre, radius, h, mmax, nev, pts)
  call check(error, pts%n == (4 * mmax + 4) * 2 * (nev + 1), "wrong number of sample points")
  if (allocated(error)) return
  call check(error, all(pts%kind == CS_WALL), "sample points are not CS_WALL")
  if (allocated(error)) return
  allocate (samples(pts%n))
  do i = 1, pts%n
    x = pts%x(:, i)
    call check(error, abs(norm2(x(1:2) - centre) - radius) <= 1d-12 * radius .and. x(3) < 0d0 .and. x(3) > -h, &
               "sample point not on the cylinder")
    if (allocated(error)) return
    call check(error, abs(sum(pts%nrm(1:2, i) * (x(1:2) - centre)) / radius - 1d0) <= 1d-12, "normal not outward")
    if (allocated(error)) return
    th = atan2(x(2) - centre(2), x(1) - centre(1))
    samples(i) = (0d0, 0d0)
    do n = 0, nev
      if (n == 0) then
        zf = (exp(k0 * x(3)) + exp(-k0 * (x(3) + 2d0 * h))) / (1d0 + exp(-2d0 * k0 * h))
      else
        zf = cos(kev(n) * (x(3) + h)) / cos(kev(n) * h)
      end if
      do m = -mmax, mmax
        samples(i) = samples(i) + b(m, n) * exp(ci * m * th) * zf
      end do
    end do
  end do

  call pw_decompose(samples, centre, radius, h, k0, kev, mmax, nev, pw)
  worst_a = maxval(abs(pw%a - a_true) / abs(a_true))
  print '(4x,a,es10.2)', "round trip: worst relative error of the coefficients ", worst_a
  write (info, '(a,es10.3)') "coefficients not recovered, worst relative error ", worst_a
  call check(error, worst_a <= 1d-12, trim(info))
  if (allocated(error)) return

  ! Reconstruction at r = 1.5 R on z = 0
  err_phi = 0d0; err_grad = 0d0; err_zz = 0d0; max_phi = 0d0; max_grad = 0d0; max_zz = 0d0
  r = 1.5d0 * radius
  do ip = 1, npt
    th = 2d0 * pi * (ip - 0.37d0) / npt
    x = [centre(1) + r * cos(th), centre(2) + r * sin(th), 0d0]
    call pw_eval_surface(pw, x, phi, grad, phizz)
    call hankel_all(mmax, k0 * r, hr)
    phi_ref = (0d0, 0d0); dr = (0d0, 0d0); dth = (0d0, 0d0); phizz_ref = (0d0, 0d0); gr_ref(3) = (0d0, 0d0)
    do n = 0, nev
      if (n > 0) then
        xk = kev(n) * r
        call pw_bessel_k(mmax + 1, xk, kvals)
      end if
      do m = -mmax, mmax
        am = abs(m)
        e = exp(ci * m * th)
        if (n == 0) then
          sgn = merge(-1d0, 1d0, m < 0 .and. mod(am, 2) == 1)
          term = sgn * hr(am)
          if (am == 0) then
            dterm = -sgn * hr(1) * k0
          else
            dterm = sgn * (hr(am - 1) - am * hr(am) / (k0 * r)) * k0
          end if
          gr_ref(3) = gr_ref(3) + a_true(m, n) * e * term * k0 * tanh(k0 * h)
          phizz_ref = phizz_ref + a_true(m, n) * e * term * k0**2
        else
          term = kvals(am)
          if (am == 0) then
            dterm = -kvals(1) * kev(n)
          else
            dterm = (-kvals(am - 1) - am * kvals(am) / xk) * kev(n)
          end if
          gr_ref(3) = gr_ref(3) - a_true(m, n) * e * term * kev(n) * tan(kev(n) * h)
          phizz_ref = phizz_ref - a_true(m, n) * e * term * kev(n)**2
        end if
        phi_ref = phi_ref + a_true(m, n) * e * term
        dr = dr + a_true(m, n) * e * dterm
        dth = dth + a_true(m, n) * e * term * ci * m
      end do
    end do
    gr_ref(1) = cos(th) * dr - sin(th) * dth / r
    gr_ref(2) = sin(th) * dr + cos(th) * dth / r
    max_phi = max(max_phi, abs(phi_ref))
    max_grad = max(max_grad, norm2(abs(gr_ref)))
    max_zz = max(max_zz, abs(phizz_ref))
    err_phi = max(err_phi, abs(phi - phi_ref))
    err_grad = max(err_grad, norm2(abs(grad - gr_ref)))
    err_zz = max(err_zz, abs(phizz - phizz_ref))
  end do
  err_phi = err_phi / max_phi
  err_grad = err_grad / max_grad
  err_zz = err_zz / max_zz
  print '(4x,a,3es10.2)', "round trip at r = 1.5 R: phi, grad phi, phi_zz errors relative to the largest ", &
    err_phi, err_grad, err_zz
  write (info, '(a,3es10.3)') "reconstruction at r = 1.5 R off (phi, grad, phi_zz) ", err_phi, err_grad, err_zz
  call check(error, err_phi <= 1d-12 .and. err_grad <= 1d-12 .and. err_zz <= 1d-12, trim(info))
end subroutine test_round_trip

! ---------------------------------------------------------------------------------------------- HAMS field

!> The scattered field of one state on the 448-panel sharp cylinder, by fields_green: sampled on r = R = 15 m, decomposed,
!> and evaluated at r = 20 m and 30 m on z = 0, against fields_green there. Hull data: the MacCamy-Fuchs scattered wave of
!> depth 15 m on the hull panels (maccamy_fuchs_state, with amp = 0 so that only the panel integrals enter, as in
!> green_perturbation_only). The unit tests have no BEM solve; whatever the hull data, a sum of Green functions over the
!> panels satisfies Laplace, the free-surface and bottom conditions and the radiation condition outside the hull, which is
!> all the decomposition assumes.
!> The gate (1e-4) is on that Green-function sum: the mesh is given to fields_green with psize = 0.09 m, so that no panel is
!> within 50 psize (4.5 m) of a point (the nearest centroid is 5.05 m from r = R) and every panel takes its centroid value.
!> With the mesh as it is, panels within 50 psize take the panel integrals of 1/r and of its image in z = 0 in place of
!> their centroid values (field_influence_row). That correction satisfies phi_z = 0 on z = 0, not the free-surface
!> condition, so it is not a sum of the modes; at r = R it is extrapolated by the propagating modes. Measured with the
!> mesh as it is: 4.9e-4 (phi) and 3.8e-4 (gradient), the same for mmax = 16..32 and nev = 8..16; the correction alone is
!> 2.5e-4 at the rings. That route is checked against 1e-3, the level at which Task 6 warns about the mode check.
subroutine test_hams_field(error)
  type(error_type), allocatable, intent(out) :: error
  type(panel_mesh) :: full, mesh, centroid_mesh
  type(first_order_state) :: states(1)
  real(8) :: err_phi, err_grad
  integer :: nthread_saved
  character(len=160) :: info

  call set_wave(hf_depth)
  nthread_saved = NTHREAD
  NTHREAD = 4
  call make_cylinder(10d0, 10d0, 32, 8, 6, full, .true.)
  call load_panels(full, 0, 0, mesh)
  call check(error, mesh%n == 448, "the fixture is not the 448-panel cylinder")
  if (allocated(error)) return
  call maccamy_fuchs_state(mesh, 10d0, hf_depth, WK, W1, 1d0, BETA, wave_gravity(W1, WK, hf_depth), 1000d0, states(1))
  states(1)%amp = 0d0
  centroid_mesh = mesh
  centroid_mesh%psize = 0.09d0

  call hams_field_errors(centroid_mesh, states, err_phi, err_grad, error)
  if (allocated(error)) then
    NTHREAD = nthread_saved
    return
  end if
  print '(4x,a,i0,a,i0,a,2es10.2)', "partial waves (mmax = ", hf_mmax, ", nev = ", hf_nev, &
    ") against the Green-function sum at r = 20, 30 m: phi and horizontal gradient errors, relative to the largest ", &
    err_phi, err_grad
  write (info, '(a,es10.3)') "phi against fields_green (centroid values), relative to max |phi_P| ", err_phi
  call check(error, err_phi <= hf_phi_limit, trim(info))
  if (allocated(error)) then
    NTHREAD = nthread_saved
    return
  end if
  write (info, '(a,es10.3)') "horizontal gradient against fields_green (centroid values), relative to the largest ", err_grad
  call check(error, err_grad <= hf_grad_limit, trim(info))
  if (allocated(error)) then
    NTHREAD = nthread_saved
    return
  end if

  call hams_field_errors(mesh, states, err_phi, err_grad, error)
  NTHREAD = nthread_saved
  if (allocated(error)) return
  print '(4x,a,2es10.2)', "the same with the near-panel integrals of fields_green (mesh as it is): phi and gradient errors ", &
    err_phi, err_grad
  write (info, '(a,2es10.3)') "partial waves against fields_green with near-panel integrals (phi, gradient) ", &
    err_phi, err_grad
  call check(error, err_phi <= hf_stock_limit .and. err_grad <= hf_stock_limit, trim(info))
end subroutine test_hams_field

!> pw_hams_field for one mesh: sample by fields_green on r = R, decompose, and compare pw_eval_surface with fields_green at
!> 16 points on each of r = 20 m and 30 m, z = 0. Errors relative to the largest |phi_P| and |grad_h phi_P| there.
subroutine hams_field_errors(mesh, states, err_phi, err_grad, error)
  type(panel_mesh), intent(in) :: mesh
  type(first_order_state), intent(in) :: states(1)
  real(8), intent(out) :: err_phi, err_grad
  type(error_type), allocatable, intent(out) :: error
  integer, parameter :: nring = 16
  real(8), parameter :: rings(2) = [20d0, 30d0], step = 0.02d0
  type(cs_points) :: pts, ring
  type(cs_fields) :: fields(1), ref(1)
  type(pw_set) :: pw
  real(8) :: kev(hf_nev), th, max_phi, max_grad
  complex(8) :: phi, grad(3), phizz
  integer :: i, l, j

  ! Sample on the cylinder; fields_green leaves phi zero at CS_WALL points below the surface, so it is given the points as
  ! CS_FREE (the kind only decides whether phi is computed)
  call pw_sample_points([0d0, 0d0], hf_radius, hf_depth, hf_mmax, hf_nev, pts)
  pts%kind = CS_FREE
  call fields_green(mesh, pts, states, step, .false., fields)
  call pw_evanescent_from_wvn(hf_nev, kev)
  call check(error, all(kev == WVN(2:hf_nev + 1)), "pw_evanescent_from_wvn does not return WVN(2:nev+1)")
  if (allocated(error)) return
  call pw_decompose(fields(1)%phi, [0d0, 0d0], hf_radius, hf_depth, WK, kev, hf_mmax, hf_nev, pw)

  ! The reference at r = 20 m and 30 m on z = 0
  ring%n = 2 * nring
  allocate (ring%x(3, ring%n), ring%nrm(3, ring%n), ring%w(ring%n), ring%psize(ring%n), ring%kind(ring%n), ring%body(ring%n))
  ring%w = 1d0
  ring%psize = 0.5d0
  ring%kind = CS_LINE
  ring%body = 0
  i = 0
  do l = 1, 2
    do j = 1, nring
      i = i + 1
      th = 2d0 * pi * (j - 0.3d0) / nring
      ring%x(:, i) = [rings(l) * cos(th), rings(l) * sin(th), 0d0]
      ring%nrm(:, i) = [cos(th), sin(th), 0d0]
    end do
  end do
  call fields_green(mesh, ring, states, step, .false., ref)

  max_phi = 0d0; max_grad = 0d0; err_phi = 0d0; err_grad = 0d0
  do i = 1, ring%n
    call pw_eval_surface(pw, ring%x(:, i), phi, grad, phizz)
    max_phi = max(max_phi, abs(ref(1)%phi(i)))
    max_grad = max(max_grad, norm2(abs(ref(1)%grad(1:2, i))))
    err_phi = max(err_phi, abs(phi - ref(1)%phi(i)))
    err_grad = max(err_grad, norm2(abs(grad(1:2) - ref(1)%grad(1:2, i))))
  end do
  err_phi = err_phi / max_phi
  err_grad = err_grad / max_grad
end subroutine hams_field_errors

end module test_partial_waves
