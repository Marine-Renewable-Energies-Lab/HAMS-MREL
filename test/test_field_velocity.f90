!> The total potential and its gradient at points in the fluid, by central differences of the Green-formula field potential
!> (FieldVelocity, Chen 2007 eq. 39 to 41). The rows are checked against the analytic panel gradient, the stencil against the
!> incident wave, and the whole route against the exact MacCamy-Fuchs field.
!> These tests set module globals, so the suite runs serially (see unit_tester.f90).
module test_field_velocity
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use WaveDyn_mod, only : BETA, H, V, W1, WK
  use PanelMesh_mod, only : ISX, ISY
  use SecondOrderTypes, only : panel_mesh, first_order_state, cs_points, cs_fields, CS_WALL, CS_LINE
  use FirstOrderState, only : image_signs, incident_phi, incident_grad, wave_gravity
  use WaterlineField, only : field_influence_row
  use FinGreen3D_Open, only : FINGREEN3D
  use SurfaceVelocity, only : source_panel_integrals, source_operator, build_source_operator, velocity_source
  use FieldVelocity, only : fields_green, fields_source
  use second_order_fixtures, only : make_cylinder, maccamy_fuchs_state, maccamy_fuchs_phi, maccamy_fuchs_field
  use test_waterline_field, only : set_wave, load_single, load_panels
  implicit none
  private

  public :: collect_field_velocity

  real(8), parameter :: pi = 3.141592653589793d0
  ! The step of the stencil, as a fraction of the panel size, as the default cs_fd_step
  real(8), parameter :: step = 0.02d0
  ! Distances of the points off a panel, in panel sizes
  real(8), parameter :: distances(3) = [1d0, 3d0, 8d0]
  ! Differenced rows against the analytic panel gradient. Finite depth, 8 sizes off: G is smooth there, and the error falls as
  ! delta^2 (measured 6.4e-6). Finite depth, 1 and 3 sizes off: the value of FINGREEN3D jumps by about 1.5e-5 in z (measured,
  ! probe at R/h = 0.58), and the difference turns a jump into jump/delta (measured worst 7.7e-3). Inferred cause: EIGEN stops
  ! its series on an absolute rule (FinGreen3D.f90:373). Deep water: the value (GF) and the gradient (GFh) of INFGREEN3D are
  ! separate approximations, and they differ by up to 3% in the real part, the same for every step (measured worst 3.2e-2).
  ! test_source_velocity already allows 6e-3 nearer the hull.
  real(8), parameter :: row_limit = 1d-4, row_jump_limit = 2d-2, row_deep_limit = 5d-2
  real(8), parameter :: incident_grad_limit = 1d-5, incident_phi_limit = 1d-12
  real(8), parameter :: exact_phi_limit = 1d-12, exact_fd_limit = 1d-6
  real(8), parameter :: green_phi_limit = 1d-2, green_grad_limit = 2d-2
  real(8), parameter :: source_phi_limit = 2d-2, source_grad_limit = 3d-2
  ! Half mesh with its images against the whole mesh: the same sums in another order (half_mesh_y_matches_whole measured 7.35e-10)
  real(8), parameter :: source_half_limit = 1d-8
  ! The size given to the test points, a typical control-surface panel
  real(8), parameter :: point_size = 0.5d0

contains

  subroutine collect_field_velocity(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("x_decide_sets_the_flag", test_x_decide), &
      new_unittest("stencil_takes_one_green_series", test_one_series), &
      new_unittest("field_green_is_smooth_in_z", test_green_smooth_in_z), &
      new_unittest("field_rows_match_analytic_gradient", test_rows_gradient), &
      new_unittest("incident_stencil", test_incident_stencil), &
      new_unittest("maccamy_fuchs_field_consistent", test_maccamy_fuchs_field), &
      new_unittest("green_field_matches_exact", test_green_field), &
      new_unittest("source_field_matches_exact", test_source_field), &
      new_unittest("source_field_half_matches_full", test_source_field_half) &
      ]
  end subroutine collect_field_velocity

  ! ---------------------------------------------------------------------------------------------- rows

  !> x_decide = xp gives the row of xp. x_decide far from every panel gives the row of a mesh whose panel sizes make every
  !> panel far, which takes the one-point value everywhere; the point itself is 1 to 3 panel sizes from the hull, where every
  !> panel is near. Both comparisons are exact. x_decide also sets the series of FINGREEN3D, so the reference row passes the
  !> same far x_decide: on that mesh it changes only the series.
  subroutine test_x_decide(error)
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: xp(3) = [6d0, 1d0, -2d0], far(3) = [400d0, 300d0, -2d0]
    type(panel_mesh) :: mesh, all_far
    complex(8), allocatable :: row_g(:, :), row_dg(:, :), same_g(:, :), same_dg(:, :), far_g(:, :), far_dg(:, :)
    complex(8), allocatable :: ref_g(:, :), ref_dg(:, :)
    real(8) :: separation

    call set_wave(10d0)
    call load_single(0, 1, mesh)
    allocate (row_g(mesh%nbase, mesh%nsys), row_dg(mesh%nbase, mesh%nsys), same_g(mesh%nbase, mesh%nsys), &
              same_dg(mesh%nbase, mesh%nsys), far_g(mesh%nbase, mesh%nsys), far_dg(mesh%nbase, mesh%nsys), &
              ref_g(mesh%nbase, mesh%nsys), ref_dg(mesh%nbase, mesh%nsys))
    call field_influence_row(mesh, xp, .false., row_g, row_dg)
    call field_influence_row(mesh, xp, .false., same_g, same_dg, xp)
    call check(error, all(same_g == row_g) .and. all(same_dg == row_dg), "x_decide = xp changes the row")
    if (allocated(error)) return
    call field_influence_row(mesh, xp, .false., far_g, far_dg, far)
    all_far = mesh
    all_far%psize = 1d-9
    call field_influence_row(all_far, xp, .false., ref_g, ref_dg, far)
    call check(error, all(far_g == ref_g) .and. all(far_dg == ref_dg), "A far x_decide does not give the one-point row")
    if (allocated(error)) return
    ! The two rows must differ, or the test would not see the decision
    separation = maxval(abs(far_g - row_g)/abs(row_g))
    print '(4x,a,es10.2)', "x_decide: largest relative change of an entry between near and far ", separation
    call check(error, separation > 1d-3, "The near and far rows are too close to tell apart")
  end subroutine test_x_decide

  !> FINGREEN3D picks its series by R/h (FINGREEN3D_R), and the series do not agree exactly at a switch. A central
  !> difference whose two points straddle a switch then takes the mismatch over 2 delta. With RR_DECIDE at the centre both
  !> points take the centre's series, and the difference matches dG/dx at the centre. Shown at R/h = 0.5 and 0.05,
  !> h = 50, omega = 0.4: the case of the box at m2 that failed the drift gate. The calls take FIELD, as the stencil points
  !> of field_influence_row do. Without FIELD the 0.05 case fails: the Limes noise of EIGENE just above R/h = 0.05 gives
  !> 2.8e-2 with one series (measured). With FIELD the jumps are smaller (2.6e-8 and 5e-9 of G, measured), so the check
  !> that the test sees them asks for the error with two series to be 10 times that with one, not above the limit.
  subroutine test_one_series(error)
    type(error_type), allocatable, intent(out) :: error

    call series_case(0.5d0, error)
    if (allocated(error)) return
    call series_case(0.05d0, error)
  end subroutine test_one_series

  subroutine series_case(r_over_h, error)
    real(8), intent(in) :: r_over_h
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: nk_local = 200
    real(8), parameter :: depth = 50d0, om = 0.4d0, g_disp = 9.807d0, delta = 1d-3, zq = -3d0, zs = -8d0
    real(8), parameter :: series_limit = 1d-3
    real(8) :: wvn_local(nk_local), nu, r0, rel_one, rel_free
    complex(8) :: g0(4), gp(4), gm(4), one, free
    character(len=160) :: info
    external :: DISPERSION

    nu = om**2/g_disp
    call DISPERSION(wvn_local, nk_local, om, depth)
    r0 = r_over_h*depth
    call FINGREEN3D(r0, 0d0, 0d0, 0d0, zq, zs, nu, wvn_local, nk_local, depth, g0, 1, FIELD=.true.)
    ! The series of the centre at both points
    call FINGREEN3D(r0 + delta, 0d0, 0d0, 0d0, zq, zs, nu, wvn_local, nk_local, depth, gp, 1, RR_DECIDE=r0, FIELD=.true.)
    call FINGREEN3D(r0 - delta, 0d0, 0d0, 0d0, zq, zs, nu, wvn_local, nk_local, depth, gm, 1, RR_DECIDE=r0, FIELD=.true.)
    one = (gp(1) - gm(1))/(2d0*delta)
    rel_one = abs(one - g0(2))/abs(g0(2))
    ! Each point its own series
    call FINGREEN3D(r0 + delta, 0d0, 0d0, 0d0, zq, zs, nu, wvn_local, nk_local, depth, gp, 1, FIELD=.true.)
    call FINGREEN3D(r0 - delta, 0d0, 0d0, 0d0, zq, zs, nu, wvn_local, nk_local, depth, gm, 1, FIELD=.true.)
    free = (gp(1) - gm(1))/(2d0*delta)
    rel_free = abs(free - g0(2))/abs(g0(2))
    print '(4x,a,f5.2,a,es10.2,a,2es10.2)', "series switch at R/h = ", r_over_h, ": jump in G ", &
      abs(free - one)*2d0*delta, ", relative error of the difference with one series and with two ", rel_one, rel_free
    write (info, '(a,f5.2,a,2es10.3)') "R/h = ", r_over_h, ": relative error with one series and with two ", &
      rel_one, rel_free
    call check(error, rel_one <= series_limit, trim(info))
    if (allocated(error)) return
    ! Without the forced choice the difference must be much worse, or the test would not show the fix matters
    call check(error, rel_free > 10d0*rel_one, "No jump seen: "//trim(info))
  end subroutine series_case

  !> The epsilon algorithm (Limes) in EIGENE returns isolated outliers: at h = 50, omega = 0.4, dG/dz_q of a bottom panel of
  !> the m1 box (centroid 296, R/h = 0.061 from control-surface point 867) jumps from 6.12e-3 to 7.37e-3 over 5e-6 m in z
  !> (measured). A difference stencil that lands on one takes it over 2 delta. With FIELD the series is summed directly. GRN(4)
  !> is scanned over +-0.01 m around the point in 2001 steps: its second difference must stay below 1e-8 of |GRN(4)|, while
  !> the default route shows the spike. The default route must also stay as it was, bit for bit, at five points that cover
  !> LINTON, PIDCOCK, EIGENE (twice, one on the spike) and EIGEN (values from the code before FIELD).
  subroutine test_green_smooth_in_z(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: nk_local = 200, nscan = 2000
    real(8), parameter :: depth = 50d0, om = 0.4d0, g_disp = 9.807d0, half_width = 0.01d0, smooth_limit = 1d-8
    real(8), parameter :: q(3) = [-4.5d0, 0.5d0, -5d0], p0(3) = [-1.702359704562689d0, -0.7631552924796104d0, -10d0]
    real(8) :: wvn_local(nk_local), nu, rr, z, worst_field, worst_default
    complex(8) :: g(4), gz_field(0:nscan), gz_default(0:nscan)
    integer :: k
    character(len=160) :: info
    external :: DISPERSION

    nu = om**2/g_disp
    call DISPERSION(wvn_local, nk_local, om, depth)
    rr = sqrt((q(1) - p0(1))**2 + (q(2) - p0(2))**2)
    do k = 0, nscan
      z = p0(3) - half_width + k*(2d0*half_width/nscan)
      call FINGREEN3D(q(1), p0(1), q(2), p0(2), q(3), z, nu, wvn_local, nk_local, depth, g, 1, RR_DECIDE=rr, FIELD=.true.)
      gz_field(k) = g(4)
      call FINGREEN3D(q(1), p0(1), q(2), p0(2), q(3), z, nu, wvn_local, nk_local, depth, g, 1, RR_DECIDE=rr)
      gz_default(k) = g(4)
    end do
    worst_field = 0d0
    worst_default = 0d0
    do k = 1, nscan - 1
      worst_field = max(worst_field, abs(gz_field(k + 1) - 2d0*gz_field(k) + gz_field(k - 1))/abs(gz_field(k)))
      worst_default = max(worst_default, abs(gz_default(k + 1) - 2d0*gz_default(k) + gz_default(k - 1))/abs(gz_default(k)))
    end do
    print '(4x,a,2es10.2)', "largest second difference of dG/dz_q over |dG/dz_q|, field and default route: ", &
      worst_field, worst_default
    write (info, '(a,es10.3)') "dG/dz_q is not smooth in z with FIELD: largest relative second difference ", worst_field
    call check(error, worst_field <= smooth_limit, trim(info))
    if (allocated(error)) return
    ! Without FIELD the spike must show, or the scan would not test the fix
    call check(error, worst_default > smooth_limit, "No spike seen on the default route")
    if (allocated(error)) return
    call check_default_route(wvn_local, nu, error)
  end subroutine test_green_smooth_in_z

  subroutine check_default_route(wvn_local, nu, error)
    real(8), intent(in) :: wvn_local(:), nu
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: depth = 50d0
    ! Panel q and field point p of each sample, and Re, Im of GRN(1:4) from the code before FIELD, built with the fpm flags.
    ! Samples 1, 4 and 5 were regenerated after the merge of main (LINTON in units of the depth, 18a3e9c, and the large-kh
    ! fix, d2e1066), which changed them by round-off only (at most 1.4e-17, 1.5e-16 relative).
    ! The values come from ifort on Windows. Another compiler rounds differently (ifx on Linux moved sample 1 by 1.4e-17,
    ! one ulp), so each sample is compared to 1e-12 of its largest value, which still catches any change of the route.
    ! Sample 3 sits on a Limes spike, which amplifies a 1-ulp change of the input, so it gets 1e-8.
    real(8), parameter :: q(3, 5) = reshape([0.5d0, 0d0, -3d0, 1.5d0, 0d0, -3d0, -4.5d0, 0.5d0, -5d0, 15d0, 0d0, -3d0, &
      40d0, 0d0, -3d0], [3, 5])
    real(8), parameter :: p(3, 5) = reshape([0d0, 0d0, -8d0, 0d0, 0d0, -8d0, &
      -1.702359704562689d0, -0.7631552924796104d0, -10.00422d0, 0d0, 0d0, -8d0, 0d0, 0d0, -8d0], [3, 5])
    real(8), parameter :: before(8, 5) = reshape([ &
      9.25266101066841012d-002, 9.29525227078842192d-002, -4.56357227742748713d-004, -1.01637852617775684d-005, &
      0.00000000000000000d+000, 0.00000000000000000d+000, 1.12007766434191365d-002, 1.46638253745496705d-003, &
      9.16220409919376588d-002, 9.29321962493177972d-002, -1.33919685829752091d-003, -3.04880218954885826d-005, &
      0.00000000000000000d+000, 0.00000000000000000d+000, 1.09740218534911289d-002, 1.46606187521771065d-003, &
      5.69473667994984528d-002, 8.74261263947343076d-002, 1.03443110685737553d-003, 5.35149698859913403d-005, &
      -4.67053298261200503d-004, -2.41624048483368242d-005, 7.36706485289481683d-003, 1.34518738883298498d-003, &
      4.30195128038764821d-002, 9.06822078938584741d-002, -3.54496613176616838d-003, -3.01182276470973818d-004, &
      0.00000000000000000d+000, 0.00000000000000000d+000, 3.26219296111452273d-003, 1.43056694148373112d-003, &
      -1.51596673435841420d-002, 7.73903715617725030d-002, -1.55335636902876207d-003, -7.44031633373910860d-004, &
      0.00000000000000000d+000, 0.00000000000000000d+000, 2.90773069836218977d-004, 1.22088014525407779d-003 &
      ], [8, 5])
    complex(8) :: g(4)
    real(8) :: now(8), tol
    integer :: k, m
    character(len=200) :: info

    do k = 1, 5
      call FINGREEN3D(q(1, k), p(1, k), q(2, k), p(2, k), q(3, k), p(3, k), nu, wvn_local, size(wvn_local), depth, g, 1)
      now = [(real(g(m), 8), aimag(g(m)), m = 1, 4)]
      write (info, '(a,i0,a,8es9.1)') "default route changed at sample ", k, ": ", now - before(:, k)
      tol = merge(1d-8, 1d-12, k == 3)*maxval(abs(before(:, k)))
      call check(error, all(abs(now - before(:, k)) <= tol), trim(info))
      if (allocated(error)) return
    end do
  end subroutine check_default_route

  !> The central difference of row_g, with the near-or-far flag decided at the centre, against grad_int of
  !> source_panel_integrals for the same panel, image and flag. The gradient of G(M x, q) in x is M grad_p G.
  subroutine test_rows_gradient(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: worst_smooth, worst_jump, worst_deep
    character(len=160) :: info_smooth, info_jump, info_deep

    ! All three run before the checks, so that a failure reports each
    call rows_against_gradient(10d0, 0, 1, [3], worst_smooth, info_smooth, error)
    if (allocated(error)) return
    call rows_against_gradient(10d0, 0, 1, [1, 2], worst_jump, info_jump, error)
    if (allocated(error)) return
    call rows_against_gradient(-1d0, 1, 0, [1, 2, 3], worst_deep, info_deep, error)
    if (allocated(error)) return
    call check(error, worst_smooth <= row_limit, "8 sizes off, "//trim(info_smooth))
    if (allocated(error)) return
    call check(error, worst_jump <= row_jump_limit, "1 and 3 sizes off, "//trim(info_jump))
    if (allocated(error)) return
    call check(error, worst_deep <= row_deep_limit, "deep, "//trim(info_deep))
  end subroutine test_rows_gradient

  subroutine rows_against_gradient(depth, isx_in, isy_in, ids, worst, info, error)
    real(8), intent(in) :: depth
    integer, intent(in) :: isx_in, isy_in
    integer, intent(in) :: ids(:)   ! the entries of distances to use
    real(8), intent(out) :: worst
    character(len=*), intent(out) :: info
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    complex(8), allocatable :: row_plus(:, :), row_minus(:, :), unused(:, :), difference(:, :, :)
    complex(8) :: g_int, grad_int(3)
    real(8) :: x0(3), xs(3), xf(3), delta, relative
    integer :: jp, id, k, c, is, jel, flag, npoint

    call set_wave(depth)
    call load_single(isx_in, isy_in, mesh)
    allocate (row_plus(mesh%nbase, mesh%nsys), row_minus(mesh%nbase, mesh%nsys), unused(mesh%nbase, mesh%nsys))
    allocate (difference(mesh%nbase, mesh%nsys, 3))
    worst = 0d0
    npoint = 0
    info = "no entry"
    do jp = 1, mesh%nbase
      do k = 1, size(ids)
        id = ids(k)
        ! Off the panel into the fluid: the normals point into the body
        x0 = mesh%xc(:, jp) - distances(id)*mesh%psize(jp)*mesh%nrm(:, jp)
        delta = step*mesh%psize(jp)
        if (depth > 0d0 .and. x0(3) - delta < -depth) cycle   ! below the seabed
        npoint = npoint + 1
        do c = 1, 3
          xs = x0
          xs(c) = x0(c) + delta
          call field_influence_row(mesh, xs, .false., row_plus, unused, x0, field_mode=.true.)
          xs(c) = x0(c) - delta
          call field_influence_row(mesh, xs, .false., row_minus, unused, x0, field_mode=.true.)
          difference(:, :, c) = (row_plus - row_minus)/(2d0*delta)
        end do
        do is = 1, mesh%nsys
          xf = image_signs(is, ISX, ISY)*x0
          do jel = 1, mesh%nbase
            flag = merge(1, 0, norm2(xf - mesh%xc(:, jel)) <= 50d0*mesh%psize(jel))
            call source_panel_integrals(mesh, xf, jel, flag, .false., g_int, grad_int)
            grad_int = image_signs(is, ISX, ISY)*grad_int
            relative = maxval(abs(difference(jel, is, :) - grad_int))/sqrt(sum(abs(grad_int)**2))
            if (relative > worst) then
              worst = relative
              write (info, '(a,es9.2,a,i0,a,f4.1,a,i0,a,i0,a,i0,a,es10.3)') "depth ", depth, ": off panel ", jp, " at ", &
                distances(id), " sizes, image ", is, " panel ", jel, " flag ", flag, ": worst relative difference ", relative
            end if
          end do
        end do
      end do
    end do
    print '(4x,a)', trim(info)
    call check(error, npoint >= mesh%nbase*size(ids)/2, "Too few points in the fluid")
  end subroutine rows_against_gradient



  ! ---------------------------------------------------------------------------------------------- incident

  !> With no perturbation the stencil returns the incident wave. The states carry nu = V, as build_state does, which differs
  !> from k tanh kh of the incident wave by 3.6e-5 in finite depth: dphi/dz at z = 0 takes the incident part exactly.
  subroutine test_incident_stencil(error)
    type(error_type), allocatable, intent(out) :: error

    call incident_case(10d0, error)
    if (allocated(error)) return
    call incident_case(-1d0, error)
  end subroutine test_incident_stencil

  subroutine incident_case(depth, error)
    real(8), intent(in) :: depth
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: headings(2) = [0.6d0, 2.d0]
    type(panel_mesh) :: mesh
    type(cs_points) :: pts
    type(first_order_state) :: states(2)
    type(cs_fields) :: fields(2)
    complex(8) :: phi_ref, grad_ref(3)
    real(8) :: g_w, worst_phi, worst_grad, rel
    integer :: i, ist
    character(len=120) :: info

    call set_wave(depth)
    call load_single(0, 0, mesh)
    call ring_points(8, 12, [9d0, 6.5d0, 11d0, 15d0], [0d0, -0.6d0, -2.5d0, -7d0], pts)
    g_w = wave_gravity(W1, WK, H)
    do ist = 1, 2
      call zero_state(mesh, headings(ist), states(ist))
    end do
    call fields_green(mesh, pts, states, step, .false., fields)
    worst_phi = 0d0
    worst_grad = 0d0
    do ist = 1, 2
      do i = 1, pts%n
        phi_ref = incident_phi(pts%x(:, i), W1, WK, H, 1d0, headings(ist), g_w)
        call incident_grad(pts%x(:, i), W1, WK, H, 1d0, headings(ist), g_w, grad_ref)
        rel = sqrt(sum(abs(fields(ist)%grad(:, i) - grad_ref)**2))/sqrt(sum(abs(grad_ref)**2))
        worst_grad = max(worst_grad, rel)
        write (info, '(a,es9.2,a,i0,a,i0,a,es10.3)') "depth ", depth, " heading ", ist, " point ", i, " gradient ", rel
        call check(error, rel <= incident_grad_limit, trim(info))
        if (allocated(error)) return
        if (pts%kind(i) == CS_WALL .and. pts%x(3, i) < -step*pts%psize(i)) then
          ! phi is not computed at a wall point below the surface
          call check(error, fields(ist)%phi(i) == (0d0, 0d0), "Wall point below the surface: phi should stay zero")
        else
          rel = abs(fields(ist)%phi(i) - phi_ref)/abs(phi_ref)
          worst_phi = max(worst_phi, rel)
          write (info, '(a,es9.2,a,i0,a,i0,a,es10.3)') "depth ", depth, " heading ", ist, " point ", i, " phi ", rel
          call check(error, rel <= incident_phi_limit, trim(info))
        end if
        if (allocated(error)) return
      end do
    end do
    print '(4x,a,es9.2,a,2es10.2)', "incident at depth ", depth, ": worst relative phi and gradient ", worst_phi, worst_grad
  end subroutine incident_case

  !> A state of the mesh with no perturbation, in the wave of set_wave at the given heading.
  subroutine zero_state(mesh, heading, s)
    type(panel_mesh), intent(in) :: mesh
    real(8), intent(in) :: heading
    type(first_order_state), intent(out) :: s

    s%omega = W1
    s%k = WK
    s%h = H
    s%amp = 1d0
    s%beta = heading
    s%g = 9.81d0
    s%rho = 1000d0
    s%nu = V
    allocate (s%phi_p(mesh%n), s%dphi_p_dn(mesh%n))
    s%phi_p = (0d0, 0d0)
    s%dphi_p_dn = (0d0, 0d0)
  end subroutine zero_state

  !> Points on circles about the z axis: for each level z(l), nper points at radius r(l). Points at z = 0 are CS_LINE, the
  !> others CS_WALL. All have the size point_size.
  subroutine ring_points(nsurface, nbelow, r, z, pts)
    integer, intent(in) :: nsurface, nbelow
    real(8), intent(in) :: r(:), z(:)
    type(cs_points), intent(out) :: pts
    integer :: l, j, i, nper

    pts%n = nsurface + nbelow
    allocate (pts%x(3, pts%n), pts%nrm(3, pts%n), pts%w(pts%n), pts%psize(pts%n), pts%kind(pts%n), pts%body(pts%n))
    pts%w = 1d0
    pts%psize = point_size
    pts%body = 1
    i = 0
    do l = 1, size(z)
      if (l == 1) then
        nper = nsurface
      else
        nper = nbelow/(size(z) - 1)
        if (l == size(z)) nper = nbelow - (size(z) - 2)*(nbelow/(size(z) - 1))
      end if
      do j = 1, nper
        i = i + 1
        pts%x(:, i) = [r(l)*cos(2d0*pi*(j - 0.3d0)/nper), r(l)*sin(2d0*pi*(j - 0.3d0)/nper), z(l)]
        pts%nrm(:, i) = [-cos(2d0*pi*(j - 0.3d0)/nper), -sin(2d0*pi*(j - 0.3d0)/nper), 0d0]
        pts%kind(i) = merge(CS_LINE, CS_WALL, z(l) == 0d0)
      end do
    end do
  end subroutine ring_points

  ! ---------------------------------------------------------------------------------------------- exact field

  !> maccamy_fuchs_field reduces to maccamy_fuchs_phi on r = a, its radial velocity vanishes there (the body condition), and
  !> its gradient matches a central difference of its potential.
  subroutine test_maccamy_fuchs_field(error)
    type(error_type), allocatable, intent(out) :: error

    call exact_case(10d0, error)
    if (allocated(error)) return
    call exact_case(-1d0, error)
  end subroutine test_maccamy_fuchs_field

  subroutine exact_case(depth, error)
    real(8), intent(in) :: depth
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: a = 5d0, om = 1.1d0, g_local = 9.81d0, beta_local = 0.6d0
    real(8), parameter :: radii(4) = [5d0, 6d0, 7.5d0, 12d0], zs(3) = [0d0, -1.5d0, -6d0]
    complex(8) :: phi, grad(3), phi_plus, phi_minus, unused(3), fd
    real(8) :: k, x(3), xs(3), th, worst_on, worst_radial, worst_fd, hstep
    integer :: ir, iz, it, c
    character(len=120) :: info

    k = wave_number(om, depth, g_local)
    hstep = 1d-5*a
    worst_on = 0d0
    worst_radial = 0d0
    worst_fd = 0d0
    do ir = 1, size(radii)
      do iz = 1, size(zs)
        do it = 1, 7
          th = 2d0*pi*(it - 0.6d0)/7d0
          x = [radii(ir)*cos(th), radii(ir)*sin(th), zs(iz)]
          call maccamy_fuchs_field(x, a, depth, k, om, 1d0, beta_local, g_local, phi, grad)
          if (ir == 1) then
            worst_on = max(worst_on, abs(phi - maccamy_fuchs_phi(x, a, depth, k, om, 1d0, beta_local, g_local))/abs(phi))
            worst_radial = max(worst_radial, abs(grad(1)*cos(th) + grad(2)*sin(th))/sqrt(sum(abs(grad)**2)))
          end if
          do c = 1, 3
            xs = x
            xs(c) = x(c) + hstep
            call maccamy_fuchs_field(xs, a, depth, k, om, 1d0, beta_local, g_local, phi_plus, unused)
            xs(c) = x(c) - hstep
            call maccamy_fuchs_field(xs, a, depth, k, om, 1d0, beta_local, g_local, phi_minus, unused)
            fd = (phi_plus - phi_minus)/(2d0*hstep)
            worst_fd = max(worst_fd, abs(fd - grad(c))/sqrt(sum(abs(grad)**2)))
          end do
        end do
      end do
    end do
    print '(4x,a,es9.2,a,3es10.2)', "exact field at depth ", depth, ": on r = a, radial velocity, difference ", &
      worst_on, worst_radial, worst_fd
    write (info, '(a,es9.2,a,es10.3)') "depth ", depth, ": against maccamy_fuchs_phi on r = a ", worst_on
    call check(error, worst_on <= exact_phi_limit, trim(info))
    if (allocated(error)) return
    write (info, '(a,es9.2,a,es10.3)') "depth ", depth, ": radial velocity on r = a ", worst_radial
    call check(error, worst_radial <= exact_fd_limit, trim(info))
    if (allocated(error)) return
    write (info, '(a,es9.2,a,es10.3)') "depth ", depth, ": gradient against a central difference ", worst_fd
    call check(error, worst_fd <= exact_fd_limit, trim(info))
  end subroutine exact_case

  !> The root of omega^2 = g k tanh kh by Newton's method, or omega^2/g in deep water (h < 0).
  real(8) function wave_number(om, depth, g_local)
    real(8), intent(in) :: om, depth, g_local
    integer :: it

    wave_number = om**2/g_local
    if (depth < 0d0) return
    do it = 1, 100
      wave_number = wave_number - (wave_number*tanh(wave_number*depth) - om**2/g_local) &
                    /(tanh(wave_number*depth) + wave_number*depth/cosh(wave_number*depth)**2)
    end do
  end function wave_number

  ! ---------------------------------------------------------------------------------------------- against the exact field

  !> The field of the exact hull potentials of maccamy_fuchs_state, on the 64 x 16 bottom-mounted mesh of drift_field_* and
  !> velocity_against_exact (a = 5, h = 10), at 40 points on r = 7.5, 10 of them at z = 0.
  !> Measured: 6.99e-4 for phi and 1.43e-3 for the gradient; limits 1e-2 and 2e-2, fixed before the measurement.
  subroutine test_green_field(error)
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: radius = 5d0, depth = 10d0
    type(panel_mesh) :: full, mesh
    type(first_order_state) :: states(1)
    type(cs_points) :: pts
    type(cs_fields) :: fields(1)
    complex(8) :: phi_ref, grad_ref(3)
    real(8) :: g_w, max_phi, max_grad, err_phi, err_grad
    integer :: i
    character(len=120) :: info

    call set_wave(depth)
    call make_cylinder(radius, depth, 64, 16, 1, full, .false.)
    call load_panels(full, 0, 0, mesh)
    g_w = wave_gravity(W1, WK, depth)
    call maccamy_fuchs_state(mesh, radius, depth, WK, W1, 1d0, BETA, g_w, 1000d0, states(1))
    call ring_points(10, 30, [7.5d0, 7.5d0, 7.5d0, 7.5d0], [0d0, -1d0, -4d0, -8d0], pts)
    call fields_green(mesh, pts, states, step, .false., fields)
    max_phi = 0d0
    max_grad = 0d0
    err_phi = 0d0
    err_grad = 0d0
    do i = 1, pts%n
      call maccamy_fuchs_field(pts%x(:, i), radius, depth, WK, W1, 1d0, BETA, g_w, phi_ref, grad_ref)
      max_grad = max(max_grad, sqrt(sum(abs(grad_ref)**2)))
      err_grad = max(err_grad, sqrt(sum(abs(fields(1)%grad(:, i) - grad_ref)**2)))
      if (pts%kind(i) /= CS_WALL) then
        max_phi = max(max_phi, abs(phi_ref))
        err_phi = max(err_phi, abs(fields(1)%phi(i) - phi_ref))
      end if
    end do
    err_phi = err_phi/max_phi
    err_grad = err_grad/max_grad
    print '(4x,a,2es10.2)', "Green field against MacCamy-Fuchs: phi and gradient errors, relative to the largest ", &
      err_phi, err_grad
    write (info, '(a,es10.3)') "phi against MacCamy-Fuchs, relative to max |phi| ", err_phi
    call check(error, err_phi <= green_phi_limit, trim(info))
    if (allocated(error)) return
    write (info, '(a,es10.3)') "gradient against MacCamy-Fuchs, relative to max |gradient| ", err_grad
    call check(error, err_grad <= green_grad_limit, trim(info))
  end subroutine test_green_field

  !> As test_green_field, by the source route: sigma from the solve of the exact normal derivative, then fields_source.
  !> Measured: 6.41e-3 for phi and 8.28e-3 for the gradient; limits 2e-2 and 3e-2, fixed before the measurement.
  subroutine test_source_field(error)
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: radius = 5d0, depth = 10d0
    type(panel_mesh) :: full, mesh
    type(first_order_state) :: states(1)
    type(source_operator) :: op
    type(cs_points) :: pts
    type(cs_fields) :: fields(1)
    complex(8), allocatable :: sigmas(:, :)
    complex(8) :: phi_ref, grad_ref(3)
    real(8) :: g_w, max_phi, max_grad, err_phi, err_grad
    integer :: i
    character(len=120) :: info

    call set_wave(depth)
    call make_cylinder(radius, depth, 64, 16, 1, full, .false.)
    call load_panels(full, 0, 0, mesh)
    g_w = wave_gravity(W1, WK, depth)
    call maccamy_fuchs_state(mesh, radius, depth, WK, W1, 1d0, BETA, g_w, 1000d0, states(1))
    allocate (sigmas(mesh%n, 1))
    call build_source_operator(mesh, .false., .false., op)
    call velocity_source(mesh, op, states(1), sigma_full=sigmas(:, 1))
    call ring_points(10, 30, [7.5d0, 7.5d0, 7.5d0, 7.5d0], [0d0, -1d0, -4d0, -8d0], pts)
    call fields_source(mesh, pts, states, sigmas, .false., fields)
    max_phi = 0d0
    max_grad = 0d0
    err_phi = 0d0
    err_grad = 0d0
    do i = 1, pts%n
      call maccamy_fuchs_field(pts%x(:, i), radius, depth, WK, W1, 1d0, BETA, g_w, phi_ref, grad_ref)
      max_grad = max(max_grad, sqrt(sum(abs(grad_ref)**2)))
      err_grad = max(err_grad, sqrt(sum(abs(fields(1)%grad(:, i) - grad_ref)**2)))
      ! phi is computed at every point by this route, wall points too
      max_phi = max(max_phi, abs(phi_ref))
      err_phi = max(err_phi, abs(fields(1)%phi(i) - phi_ref))
    end do
    err_phi = err_phi/max_phi
    err_grad = err_grad/max_grad
    print '(4x,a,2es10.2)', "Source field against MacCamy-Fuchs: phi and gradient errors, relative to the largest ", &
      err_phi, err_grad
    write (info, '(a,es10.3)') "phi against MacCamy-Fuchs, relative to max |phi| ", err_phi
    call check(error, err_phi <= source_phi_limit, trim(info))
    if (allocated(error)) return
    write (info, '(a,es10.3)') "gradient against MacCamy-Fuchs, relative to max |gradient| ", err_grad
    call check(error, err_grad <= source_grad_limit, trim(info))
  end subroutine test_source_field

  !> The symmetry-image path of fields_source: a hull with a plane of symmetry (nsys = 2, the images of sigma_full and the mirrored
  !> field point) against the same hull whole (nsys = 1), same state and points, 20 points around the hull with y < 0 and x < 0
  !> among them, 8 at z = 0. Both planes: ISY = 1 at depth 10 and ISX = 1 in deep water.
  subroutine test_source_field_half(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: rel_phi, rel_grad

    call half_against_whole(10d0, 0, 1, rel_phi, rel_grad, error)
    if (allocated(error)) return
    print '(4x,a,2es10.2)', "Source field, ISY = 1 half against whole: phi and gradient ", rel_phi, rel_grad
    call check(error, max(rel_phi, rel_grad) <= source_half_limit, "ISY = 1 half mesh against whole mesh")
    if (allocated(error)) return
    call half_against_whole(-1d0, 1, 0, rel_phi, rel_grad, error)
    if (allocated(error)) return
    print '(4x,a,2es10.2)', "Source field, ISX = 1 half against whole: phi and gradient ", rel_phi, rel_grad
    call check(error, max(rel_phi, rel_grad) <= source_half_limit, "ISX = 1 half mesh against whole mesh")
  end subroutine test_source_field_half

  subroutine half_against_whole(depth, isx_in, isy_in, rel_phi, rel_grad, error)
    real(8), intent(in) :: depth
    integer, intent(in) :: isx_in, isy_in
    real(8), intent(out) :: rel_phi, rel_grad
    type(error_type), allocatable, intent(out) :: error
    type(panel_mesh) :: mesh
    type(first_order_state) :: states(1)
    type(source_operator) :: op
    type(cs_points) :: pts
    type(cs_fields) :: fields_w(1), fields_h(1)
    complex(8), allocatable :: sigmas(:, :)
    real(8) :: g_w
    integer :: pass, i

    call ring_points(8, 12, [9d0, 6.5d0, 11d0, 15d0], [0d0, -0.6d0, -2.5d0, -7d0], pts)
    rel_phi = 0d0
    rel_grad = 0d0
    do pass = 1, 2
      call set_wave(depth)
      if (pass == 1) then
        call load_single(0, 0, mesh)
      else
        call load_single(isx_in, isy_in, mesh)
      end if
      g_w = wave_gravity(W1, WK, H)
      call zero_state(mesh, 0.6d0, states(1))
      allocate (states(1)%grad_p(3, mesh%n), states(1)%grad_i(3, mesh%n))
      do i = 1, mesh%n
        call incident_grad(mesh%xc(:, i), W1, WK, H, 1d0, 0.6d0, g_w, states(1)%grad_i(:, i))
        states(1)%dphi_p_dn(i) = -sum(states(1)%grad_i(:, i)*mesh%nrm(:, i))
      end do
      allocate (sigmas(mesh%n, 1))
      call build_source_operator(mesh, .false., .false., op)
      call velocity_source(mesh, op, states(1), sigma_full=sigmas(:, 1))
      if (pass == 1) then
        call check(error, mesh%nsys == 1, "The whole mesh should have no symmetry")
        call fields_source(mesh, pts, states, sigmas, .false., fields_w)
      else
        call check(error, mesh%nsys == 2, "The half mesh should have one plane")
        call fields_source(mesh, pts, states, sigmas, .false., fields_h)
      end if
      if (allocated(error)) return
      deallocate (sigmas, op%lu, op%ipiv, op%gradient)
    end do
    rel_phi = maxval(abs(fields_h(1)%phi - fields_w(1)%phi))/maxval(abs(fields_w(1)%phi))
    rel_grad = maxval(abs(fields_h(1)%grad - fields_w(1)%grad))/maxval(abs(fields_w(1)%grad))
  end subroutine half_against_whole

end module test_field_velocity
