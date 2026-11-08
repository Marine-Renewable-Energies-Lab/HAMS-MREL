!> The body velocity by source solve: the gradient of the Green function panel integrals against a finite difference of the same
!> integrals in the field point, then the solve against MacCamy-Fuchs and half meshes against whole ones.
!> These tests set module globals, so the suite runs serially (see unit_tester.f90). The limits follow what the Green function
!> routines themselves hold, so that a limit fails on a wrong sign and not on their noise. See the comments at each.
module test_source_velocity
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use WaveDyn_mod, only : BETA, H, W1, WK
  use PanelMesh_mod, only : ISX, ISY
  use SecondOrderTypes, only : panel_mesh, first_order_state
  use FirstOrderState, only : image_signs, incident_grad, wave_gravity
  use SurfaceVelocity, only : source_operator, source_panel_integrals, build_source_operator, velocity_source
  use second_order_fixtures, only : make_cylinder, maccamy_fuchs_state
  use test_waterline_field, only : set_wave, load_single, load_multi, load_panels
  implicit none
  private

  public :: collect_source_velocity

  ! Field points in the fluid around the cylinder of load_single, at all the distances the Green function routines separate
  real(8), parameter :: points(3, 7) = reshape([ &
    5.3d0, 0.4d0, -0.15d0, &    ! next to the top row, close to the free surface
    5.2d0, 1.5d0, -2.1d0, &     ! next to the wall
    2.d0, 1.d0, -5.4d0, &       ! under the bottom
    6.d0, 4.d0, -1.d0, &
    12.d0, 8.d0, -3.d0, &
    60.d0, 20.d0, -2.d0, &
    400.d0, 300.d0, -2.d0], [3, 7])

  ! Depth 10: EIGENE at most of these distances, with no noise that a small step shows. The largest value measured is 5.7e-4.
  real(8), parameter :: step_clean = 1d-3, limit_clean = 1d-3
  ! Deep water: the value and the gradient differ by up to 0.46% of the gradient over the near panels. The difference converges
  ! with the step, so it is not noise. Inference, not shown: they come from separate approximations in HavelockGF. The worst
  ! component is z (4.6e-3), and GRN(4) comes from GF, like the value, so the inference does not explain that component.
  real(8), parameter :: limit_deep = 6d-3
  ! Depth 100: PIDCOCK near the body. A move of 1e-9 changes the value by up to 1.8e-4 of it, so the step has to be wide.
  real(8), parameter :: step_wide = 5d-2, limit_wide = 1d-2
  ! The same, against deep water, where only the sea bed differs. The largest value measured at the points of the loop is 1.73e-3.
  real(8), parameter :: limit_against_deep = 5d-3
  real(8), parameter :: solve_limit = 1d-9
  ! sigma_full against grad_p: the same sums in another order, in double precision
  real(8), parameter :: sigma_full_limit = 1d-10
  real(8), parameter :: pi = 3.141592653589793d0
  ! Whole mesh against half mesh: the two sum the same terms in another order
  real(8), parameter :: symmetry_limit = 1d-8
  ! MacCamy-Fuchs: rms over the wall of the velocity difference, relative to the rms total velocity. Measured 1.68% on 32 x 8
  ! panels and 0.77% on 64 x 16, which is first order in panel size.
  real(8), parameter :: analytic_limit = 3d-2, analytic_fine_limit = 1.5d-2

contains

  subroutine collect_source_velocity(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("gradient_shallow", test_gradient_shallow), &
      new_unittest("gradient_deep", test_gradient_deep), &
      new_unittest("gradient_depth_100", test_gradient_depth_100), &
      new_unittest("gradient_multi_arrays", test_gradient_multi), &
      new_unittest("gradient_depth_100_against_deep", test_against_deep), &
      new_unittest("gradient_on_own_panel", test_gradient_on_own_panel), &
      new_unittest("solve_recovers_normal_derivative", test_solve_normal), &
      new_unittest("half_mesh_y_matches_whole", test_half_y), &
      new_unittest("half_mesh_x_matches_whole", test_half_x), &
      new_unittest("maccamy_fuchs_velocity", test_maccamy_fuchs), &
      new_unittest("sigma_full_reproduces_velocity", test_sigma_full) &
      ]
  end subroutine collect_source_velocity

  subroutine test_gradient_shallow(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    call set_wave(10d0)
    call load_single(1, 0, mesh)
    call compare_gradient(mesh, .false., "depth 10", step_clean, limit_clean, error)
  end subroutine test_gradient_shallow

  subroutine test_gradient_deep(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    call set_wave(-1d0)
    call load_single(0, 1, mesh)
    call compare_gradient(mesh, .false., "deep", step_clean, limit_deep, error)
  end subroutine test_gradient_deep

  subroutine test_gradient_depth_100(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    call set_wave(100d0)
    call load_single(0, 1, mesh)
    call compare_gradient(mesh, .false., "depth 100", step_wide, limit_wide, error)
  end subroutine test_gradient_depth_100

  subroutine test_gradient_multi(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    call set_wave(10d0)
    call load_multi(1, mesh)
    call compare_gradient(mesh, .true., "multi-body arrays", step_clean, limit_clean, error)
  end subroutine test_gradient_multi

  !> Against the difference of the same integrals, for each panel, field point, symmetry image and flag. Flag 0 means one point on
  !> the panel, which the solver uses far away. Here it also runs close, where the value and the gradient still have to agree.
  subroutine compare_gradient(mesh, is_multi, label, step, limit, error)
    type(panel_mesh), intent(in) :: mesh
    logical, intent(in) :: is_multi
    character(len=*), intent(in) :: label
    real(8), intent(in) :: step, limit
    type(error_type), allocatable, intent(out) :: error
    real(8) :: xf(3), offset(3), worst(0:1, 2), relative
    complex(8) :: g, g_plus, g_minus, grad(3), unused(3), difference
    integer :: ip, jel, is, flag, c
    character(len=150) :: info

    worst = 0d0
    do ip = 1, size(points, 2)
      do is = 1, mesh%nsys
        xf = image_signs(is, ISX, ISY)*points(:, ip)
        do jel = 1, mesh%nbase
          do flag = 0, 1
            call source_panel_integrals(mesh, xf, jel, flag, is_multi, g, grad)
            do c = 1, 3
              offset = 0d0
              offset(c) = step
              call source_panel_integrals(mesh, xf + offset, jel, flag, is_multi, g_plus, unused)
              call source_panel_integrals(mesh, xf - offset, jel, flag, is_multi, g_minus, unused)
              difference = (g_plus - g_minus)/(2d0*step)
              relative = abs(grad(c) - difference)/sqrt(sum(abs(grad)**2))
              worst(flag, is) = max(worst(flag, is), relative)
              write (info, '(a,a,4(a,i0),a,i0,a,es10.3)') label, ":", " point ", ip, " image ", is, " panel ", jel, " flag ", &
                flag, " component ", c, " relative difference ", relative
              call check(error, relative <= limit, trim(info))
              if (allocated(error)) return
            end do
          end do
        end do
      end do
    end do
    print '(4x,a,a,4es10.2)', label, ": worst relative difference, flags 0 and 1 of image 1, then of image 2 ", worst
  end subroutine compare_gradient

  !> At depth 100, the gradient with the noise that a difference would need a wide step for, against the deep water routine.
  subroutine test_against_deep(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    complex(8) :: g, grad(3), grad_deep(3)
    real(8) :: worst
    integer :: ip, jel

    call set_wave(100d0)
    call load_single(0, 0, mesh)
    worst = 0d0
    do ip = 1, size(points, 2) - 1
      do jel = 1, mesh%nbase
        call set_wave(100d0)
        call source_panel_integrals(mesh, points(:, ip), jel, 1, .false., g, grad)
        call set_wave(-1d0)
        call source_panel_integrals(mesh, points(:, ip), jel, 1, .false., g, grad_deep)
        worst = max(worst, sqrt(sum(abs(grad - grad_deep)**2))/sqrt(sum(abs(grad)**2)))
      end do
    end do
    print '(4x,a,es10.2)', "depth 100 against deep: worst relative difference ", worst
    call check(error, worst <= limit_against_deep, "Gradient at depth 100 against deep water")
  end subroutine test_against_deep

  !> At the centroid of a panel, within its plane: the 1/r part is continuous there.
  !> Only the in-plane part is checked, not the normal part.
  subroutine test_gradient_on_own_panel(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    real(8) :: t(3), a(3), nv(3), worst
    complex(8) :: g, g_plus, g_minus, grad(3), unused(3), difference
    integer :: jel

    call set_wave(10d0)
    call load_single(0, 0, mesh)
    worst = 0d0
    do jel = 1, mesh%nbase
      nv = mesh%nrm(:, jel)
      a = 0d0
      a(minloc(abs(nv), 1)) = 1d0
      t = a - dot_product(a, nv)*nv
      t = t/sqrt(dot_product(t, t))
      call source_panel_integrals(mesh, mesh%xc(:, jel), jel, 1, .false., g, grad)
      call source_panel_integrals(mesh, mesh%xc(:, jel) + step_clean*t, jel, 1, .false., g_plus, unused)
      call source_panel_integrals(mesh, mesh%xc(:, jel) - step_clean*t, jel, 1, .false., g_minus, unused)
      difference = (g_plus - g_minus)/(2d0*step_clean)
      worst = max(worst, abs(sum(grad*t) - difference)/sqrt(sum(abs(grad)**2)))
    end do
    print '(4x,a,es10.2)', "own panel: worst relative difference ", worst
    call check(error, worst <= limit_clean, "In-plane gradient at the centroid of the own panel")
  end subroutine test_gradient_on_own_panel

  !> The normal derivative of the perturbation on a fixed body in the wave of set_wave.
  subroutine fixed_body_state(mesh, s)
    type(panel_mesh), intent(in) :: mesh
    type(first_order_state), intent(out) :: s
    real(8) :: g_wave
    integer :: i

    allocate (s%grad_p(3, mesh%n), s%dphi_p_dn(mesh%n), s%grad_i(3, mesh%n))
    g_wave = wave_gravity(W1, WK, H)
    do i = 1, mesh%n
      call incident_grad(mesh%xc(:, i), W1, WK, H, 1d0, BETA, g_wave, s%grad_i(:, i))
      s%dphi_p_dn(i) = -sum(s%grad_i(:, i)*mesh%nrm(:, i))
    end do
  end subroutine fixed_body_state

  subroutine test_solve_normal(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    type(first_order_state) :: s
    type(source_operator) :: op
    integer :: i
    real(8) :: worst

    call set_wave(10d0)
    call load_single(0, 1, mesh)
    call fixed_body_state(mesh, s)
    call build_source_operator(mesh, .false., .false., op)
    call velocity_source(mesh, op, s)
    worst = 0d0
    do i = 1, mesh%n
      worst = max(worst, abs(sum(s%grad_p(:, i)*mesh%nrm(:, i)) - s%dphi_p_dn(i)))
    end do
    worst = worst/maxval(abs(s%dphi_p_dn))
    print '(4x,a,es10.2)', "normal derivative: worst difference, relative to the largest ", worst
    call check(error, worst <= solve_limit, "The gradient does not recover the normal derivative")
  end subroutine test_solve_normal

  subroutine test_half_y(error)
    type(error_type), allocatable, intent(out) :: error
    call compare_half(10d0, 0, 1, error)
  end subroutine test_half_y

  subroutine test_half_x(error)
    type(error_type), allocatable, intent(out) :: error
    call compare_half(-1d0, 1, 0, error)
  end subroutine test_half_x

  !> Every panel of the whole mesh against the same panel of a half mesh or of its image, in gradient and potential.
  subroutine compare_half(depth, isx_in, isy_in, error)
    real(8), intent(in) :: depth
    integer, intent(in) :: isx_in, isy_in
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: whole, half
    type(first_order_state) :: s_whole, s_half
    type(source_operator) :: op_whole, op_half
    complex(8), allocatable :: phi_whole(:), phi_half(:)
    real(8) :: worst, scale
    integer :: i, j, matched

    call set_wave(depth)
    call load_single(0, 0, whole)
    call fixed_body_state(whole, s_whole)
    allocate (phi_whole(whole%n))
    call build_source_operator(whole, .false., .true., op_whole)
    call velocity_source(whole, op_whole, s_whole, phi_whole)

    call load_single(isx_in, isy_in, half)
    call fixed_body_state(half, s_half)
    allocate (phi_half(half%n))
    call build_source_operator(half, .false., .true., op_half)
    call velocity_source(half, op_half, s_half, phi_half)

    call check(error, half%n == whole%n, "The half mesh and its images do not make the whole mesh")
    if (allocated(error)) return
    scale = maxval(abs(s_whole%grad_p))
    worst = 0d0
    matched = 0
    do i = 1, whole%n
      do j = 1, half%n
        if (sum((whole%xc(:, i) - half%xc(:, j))**2) > 1d-16) cycle
        matched = matched + 1
        worst = max(worst, maxval(abs(s_whole%grad_p(:, i) - s_half%grad_p(:, j)))/scale, &
                    abs(phi_whole(i) - phi_half(j))/maxval(abs(phi_whole)))
      end do
    end do
    print '(4x,a,es10.2)', "half against whole mesh: worst relative difference ", worst
    call check(error, matched == whole%n, "Not every panel of the whole mesh has its match")
    if (allocated(error)) return
    call check(error, worst <= symmetry_limit, "Half mesh against whole mesh")
  end subroutine compare_half

  !> With a symmetry plane (ISY = 1), the gradient at every panel centroid, base and image, from sigma_full: the jump term
  !> 2 pi sigma n, plus the sum over every panel (image is of base panel j, full index (is-1) nbase + j) of sigma M grad_int at
  !> the mirrored centroid, over 4 pi. It must equal grad_p of velocity_source.
  subroutine test_sigma_full(error)
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    type(first_order_state) :: s
    type(source_operator) :: op
    complex(8), allocatable :: sigma_full(:)
    complex(8) :: g_int, grad_int(3), grad(3)
    real(8) :: xf(3), f(3), worst, scale
    integer :: i, is, j, flag

    call set_wave(10d0)
    call load_single(0, 1, mesh)
    call fixed_body_state(mesh, s)
    allocate (sigma_full(mesh%n))
    call build_source_operator(mesh, .false., .false., op)
    call velocity_source(mesh, op, s, sigma_full=sigma_full)
    call check(error, mesh%nsys == 2, "The fixture should have one symmetry plane")
    if (allocated(error)) return
    scale = maxval(abs(s%grad_p))
    worst = 0d0
    do i = 1, mesh%n
      grad = 2d0*pi*sigma_full(i)*mesh%nrm(:, i)
      do is = 1, mesh%nsys
        f = image_signs(is, ISX, ISY)
        xf = f*mesh%xc(:, i)
        do j = 1, mesh%nbase
          flag = merge(1, 0, norm2(xf - mesh%xc(:, j)) <= 50d0*mesh%psize(j))
          call source_panel_integrals(mesh, xf, j, flag, .false., g_int, grad_int)
          grad = grad + sigma_full((is - 1)*mesh%nbase + j)*f*grad_int
        end do
      end do
      worst = max(worst, maxval(abs(grad/(4d0*pi) - s%grad_p(:, i)))/scale)
    end do
    print '(4x,a,es10.2)', "sigma_full: worst difference to grad_p, relative to the largest ", worst
    call check(error, worst <= sigma_full_limit, "sigma_full does not reproduce grad_p")
  end subroutine test_sigma_full

  subroutine test_maccamy_fuchs(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: relative

    call velocity_against_exact(32, 8, relative)
    print '(4x,a,es10.2)', "MacCamy-Fuchs, 32 x 8: rms velocity difference over rms velocity ", relative
    call check(error, relative <= analytic_limit, "MacCamy-Fuchs velocity on 32 x 8 panels")
    if (allocated(error)) return
    call velocity_against_exact(64, 16, relative)
    print '(4x,a,es10.2)', "MacCamy-Fuchs, 64 x 16: rms velocity difference over rms velocity ", relative
    call check(error, relative <= analytic_fine_limit, "MacCamy-Fuchs velocity on 64 x 16 panels")
  end subroutine test_maccamy_fuchs

  !> A fixed bottom-mounted cylinder has an exact velocity, which the solve never sees: it only takes the normal derivative.
  subroutine velocity_against_exact(ntheta, nz, relative)
    integer, intent(in) :: ntheta, nz
    real(8), intent(out) :: relative
    real(8), parameter :: radius = 5d0, depth = 10d0
    type(panel_mesh) :: full, mesh
    type(first_order_state) :: s
    type(source_operator) :: op
    complex(8), allocatable :: exact(:, :)
    real(8) :: rms_diff, rms_total
    integer :: i

    call set_wave(depth)
    call make_cylinder(radius, depth, ntheta, nz, 1, full, .false.)
    call load_panels(full, 0, 0, mesh)
    call maccamy_fuchs_state(mesh, radius, depth, WK, W1, 1d0, BETA, wave_gravity(W1, WK, depth), 1000d0, s)
    exact = s%grad_p
    call build_source_operator(mesh, .false., .false., op)
    call velocity_source(mesh, op, s)
    rms_diff = 0d0
    rms_total = 0d0
    do i = 1, mesh%n
      rms_diff = rms_diff + mesh%area(i)*sum(abs(s%grad_p(:, i) - exact(:, i))**2)
      rms_total = rms_total + mesh%area(i)*sum(abs(s%grad_i(:, i) + exact(:, i))**2)
    end do
    relative = sqrt(rms_diff/rms_total)
  end subroutine velocity_against_exact

end module test_source_velocity
