!> Rankine integrals of a constant normal-dipole panel (DipoleIntgr): the solid angle Omega and the finite-part integral
!> of d2(1/r)/dn_F dn_Q, which is n_F . grad Omega (docs/developer_documentation.md).
!> The references here are independent of the module: Gauss quadrature of the analytic kernels, the closed form of the
!> finite part on a square from polar coordinates, central differences of Omega, and the solid angle of a closed cube.
module test_dipole_kernel

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use DipoleIntgr, only : DIPOLE_SOLID_ANGLE, DIPOLE_GRAD_OMEGA, DIPOLE_RANKINE_NN

  implicit none
  private
  public :: collect_dipole_kernel

  real(8), parameter :: pi = 3.141592653589793d0

contains

  subroutine collect_dipole_kernel(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("dipole_nn_far_matches_quadrature", test_far_quadrature), &
      new_unittest("dipole_omega_far_matches_quadrature", test_omega_quadrature), &
      new_unittest("dipole_nn_self_square", test_self_square), &
      new_unittest("dipole_nn_matches_omega_differences", test_grad_fd), &
      new_unittest("dipole_nn_near_edge_finite", test_near_edge), &
      new_unittest("dipole_image_term", test_image), &
      new_unittest("dipole_cube_solid_angle", test_cube), &
      new_unittest("dipole_normal_flip", test_flip) &
    ]
  end subroutine collect_dipole_kernel

  !> Nodes and weights of the n-point Gauss-Legendre rule on (-1, 1), by Newton on the Legendre recurrence.
  subroutine gauss_legendre(n, x, w)
    integer, intent(in) :: n
    real(8), intent(out) :: x(n), w(n)
    integer :: i, j, it
    real(8) :: z, p0, p1, p2, dp
    do i = 1, n
      z = cos(pi*(i - 0.25d0)/(n + 0.5d0))
      do it = 1, 100
        p0 = 1d0
        p1 = z
        do j = 2, n
          p2 = ((2*j - 1)*z*p1 - (j - 1)*p0)/j
          p0 = p1
          p1 = p2
        end do
        dp = n*(z*p1 - p0)/(z*z - 1d0)
        z = z - p1/dp
        if (abs(p1/dp) < 1d-16) exit
      end do
      x(i) = z
      w(i) = 2d0/((1d0 - z*z)*dp*dp)
    end do
  end subroutine gauss_legendre

  !> A skewed flat quadrilateral in the plane z = z0 tilted about x, with its unit normal along the node order.
  subroutine skew_quad(a, nq)
    real(8), intent(out) :: a(4,3), nq(3)
    real(8) :: c, s, u(3), w(3)
    integer :: k
    real(8) :: p(4,2)
    p(1,:) = [0.0d0, 0.0d0]
    p(2,:) = [1.2d0, 0.1d0]
    p(3,:) = [1.0d0, 0.9d0]
    p(4,:) = [-0.1d0, 0.7d0]
    c = cos(0.4d0)
    s = sin(0.4d0)
    do k = 1, 4
      a(k,:) = [p(k,1), c*p(k,2), -3d0 + s*p(k,2)]
    end do
    u = a(3,:) - a(1,:)
    w = a(4,:) - a(2,:)
    nq = [u(2)*w(3) - u(3)*w(2), u(3)*w(1) - u(1)*w(3), u(1)*w(2) - u(2)*w(1)]
    nq = nq/norm2(nq)
  end subroutine skew_quad

  !> Gauss quadrature over the bilinear map of the panel of f(x, xi), with n x n points.
  function quad_kernel(a, xf, nf, nq, n, kind) result(s)
    real(8), intent(in) :: a(4,3), xf(3), nf(3), nq(3)
    integer, intent(in) :: n, kind
    real(8) :: s
    real(8) :: gx(n), gw(n), u, v, q(3), du(3), dv(3), cn(3), r(3), rr
    integer :: i, j
    call gauss_legendre(n, gx, gw)
    s = 0d0
    do i = 1, n
      do j = 1, n
        u = 0.5d0*(1d0 + gx(i))
        v = 0.5d0*(1d0 + gx(j))
        q = (1-u)*(1-v)*a(1,:) + u*(1-v)*a(2,:) + u*v*a(3,:) + (1-u)*v*a(4,:)
        du = (1-v)*(a(2,:) - a(1,:)) + v*(a(3,:) - a(4,:))
        dv = (1-u)*(a(4,:) - a(1,:)) + u*(a(3,:) - a(2,:))
        cn = [du(2)*dv(3) - du(3)*dv(2), du(3)*dv(1) - du(1)*dv(3), du(1)*dv(2) - du(2)*dv(1)]
        r = xf - q
        rr = norm2(r)
        if (kind == 1) then        ! d(1/r)/dn_Q = (x - xi).n_Q / r^3
          s = s + 0.25d0*gw(i)*gw(j)*norm2(cn)*dot_product(r, nq)/rr**3
        else                       ! d2(1/r)/dn_F dn_Q = n_F.n_Q / r^3 - 3 (n_F.r)(n_Q.r) / r^5
          s = s + 0.25d0*gw(i)*gw(j)*norm2(cn)*(dot_product(nf, nq)/rr**3 &
                  - 3d0*dot_product(nf, r)*dot_product(nq, r)/rr**5)
        end if
      end do
    end do
  end function quad_kernel

  !> Five panel sizes away, the closed form against 20 x 20 Gauss points of the analytic kernel.
  subroutine test_far_quadrature(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: a(4,3), nq(3), xf(3), nf(3), enn, ref
    call skew_quad(a, nq)
    xf = [0.5d0, 0.4d0, -3d0] + 5.5d0*[0.3d0, -0.5d0, 0.81d0]/norm2([0.3d0, -0.5d0, 0.81d0])
    nf = [0.2d0, 0.6d0, -0.77d0]
    nf = nf/norm2(nf)
    call DIPOLE_RANKINE_NN(xf, nf, a, 4, nq, 0d0, enn)
    ref = quad_kernel(a, xf, nf, nq, 20, 2)
    call check(error, abs(enn - ref) <= 1d-10*abs(ref), "far: closed form vs Gauss quadrature")
  end subroutine test_far_quadrature

  !> The solid angle against Gauss quadrature of (x - xi).n_Q / r^3, on both sides of the panel.
  subroutine test_omega_quadrature(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: a(4,3), nq(3), xf(3), om, ref
    integer :: k
    call skew_quad(a, nq)
    do k = -1, 1, 2
      xf = [0.5d0, 0.4d0, -3d0] + k*4d0*[0.3d0, -0.5d0, 0.81d0]
      om = DIPOLE_SOLID_ANGLE(xf, a, 4, nq)
      ref = quad_kernel(a, xf, nq, nq, 20, 1)
      call check(error, abs(om - ref) <= 1d-10*abs(ref), "solid angle vs Gauss quadrature")
      if (allocated(error)) return
    end do
  end subroutine test_omega_quadrature

  !> Unit square at its own centroid: -8 sqrt(2) by hand (design note Sec. 3.2), and the finite part from polar
  !> coordinates: f.p. over the inscribed disc of radius 1/2 is -2 pi/(1/2), plus the regular rest by Gauss quadrature.
  subroutine test_self_square(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: a(4,3), nq(3), xf(3), enn, rest, th, rmax
    real(8) :: gx(40), gw(40)
    integer :: i
    a(1,:) = [-0.5d0, -0.5d0, -2d0]
    a(2,:) = [0.5d0, -0.5d0, -2d0]
    a(3,:) = [0.5d0, 0.5d0, -2d0]
    a(4,:) = [-0.5d0, 0.5d0, -2d0]
    nq = [0d0, 0d0, 1d0]
    xf = [0d0, 0d0, -2d0]
    call DIPOLE_RANKINE_NN(xf, nq, a, 4, nq, 0d0, enn)
    call check(error, abs(enn + 8d0*sqrt(2d0)) <= 1d-12, "self term of the unit square is -8 sqrt(2)")
    if (allocated(error)) return
    ! In plane, n_F.n_Q / r^3 = 1/rho^3. The rest, square minus disc: 8 x int_0^(pi/4) int_(1/2)^(rmax) rho^-2 drho dth
    call gauss_legendre(40, gx, gw)
    rest = 0d0
    do i = 1, 40
      th = 0.125d0*pi*(1d0 + gx(i))
      rmax = 0.5d0/cos(th)
      rest = rest + 0.125d0*pi*gw(i)*(1d0/0.5d0 - 1d0/rmax)
    end do
    rest = 8d0*rest
    call check(error, abs(enn - (-2d0*pi/0.5d0 + rest)) <= 1d-6*abs(enn), "self term vs polar finite part")
  end subroutine test_self_square

  !> n_F . grad Omega against central differences of Omega, at points around the skewed panel (also checks orientation).
  subroutine test_grad_fd(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: a(4,3), nq(3), xf(3), g(3), fd, h
    integer :: k, d
    real(8) :: pts(3,4)
    call skew_quad(a, nq)
    pts(:,1) = [0.5d0, 0.3d0, -2.2d0]
    pts(:,2) = [0.5d0, 0.3d0, -3.9d0]
    pts(:,3) = [1.6d0, 0.2d0, -3.0d0]
    pts(:,4) = [0.2d0, 0.5d0, -2.75d0]
    h = 1d-5
    do k = 1, 4
      xf = pts(:,k)
      call DIPOLE_GRAD_OMEGA(xf, a, 4, nq, g)
      do d = 1, 3
        fd = (DIPOLE_SOLID_ANGLE(xf + h*unit(d), a, 4, nq) - DIPOLE_SOLID_ANGLE(xf - h*unit(d), a, 4, nq))/(2d0*h)
        call check(error, abs(g(d) - fd) <= 1d-6*max(1d0, abs(fd)), "grad Omega vs central differences")
        if (allocated(error)) return
      end do
    end do
  contains
    function unit(d) result(e)
      integer, intent(in) :: d
      real(8) :: e(3)
      e = 0d0
      e(d) = 1d0
    end function unit
  end subroutine test_grad_fd

  !> Points 1e-6, 1e-3 and 0.1 sizes off an edge of the unit square, in its plane (outside) and above the edge: the result
  !> is finite and equals central differences of Omega (so continuous with it); exactly on the edge it is finite.
  subroutine test_near_edge(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: a(4,3), nq(3), xf(3), enn, fd, h, d, e0
    real(8), parameter :: ds(3) = [1d-6, 1d-3, 1d-1]
    integer :: k, m
    a(1,:) = [0d0, 0d0, -1d0]
    a(2,:) = [1d0, 0d0, -1d0]
    a(3,:) = [1d0, 1d0, -1d0]
    a(4,:) = [0d0, 1d0, -1d0]
    nq = [0d0, 0d0, 1d0]
    do m = 1, 2
      do k = 1, 3
        d = ds(k)
        if (m == 1) then
          xf = [1d0 + d, 0.37d0, -1d0]            ! in the plane, outside the edge x = 1
        else
          xf = [1d0, 0.37d0, -1d0 + d]            ! above the edge
        end if
        call DIPOLE_RANKINE_NN(xf, nq, a, 4, nq, 0d0, enn)
        h = 1d-3*d
        fd = (DIPOLE_SOLID_ANGLE(xf + [0d0, 0d0, h], a, 4, nq) - DIPOLE_SOLID_ANGLE(xf - [0d0, 0d0, h], a, 4, nq))/(2d0*h)
        call check(error, abs(enn) < huge(1d0) .and. enn == enn, "finite near an edge")
        if (allocated(error)) return
        print '(2x,a,2i3,3es14.6)', 'near edge', m, k, d, enn, fd
        if (m == 2 .and. k == 1) then
          ! 1e-6 above the edge, differences of Omega with a 1e-9 step are round-off; check continuity with the edge
          ! value (the edge skipped: its own field above it is parallel to the plane, so it adds nothing to n . grad)
          call DIPOLE_RANKINE_NN([1d0, 0.37d0, -1d0], nq, a, 4, nq, 0d0, e0)
          call check(error, abs(enn - e0) <= 1d-4*abs(e0), "1e-6 above an edge: continuous with the edge value")
        else
          call check(error, abs(enn - fd) <= 1d-4*abs(fd), "near an edge: vs differences of Omega")
        end if
        if (allocated(error)) return
      end do
    end do
    xf = [1d0, 0.37d0, -1d0]
    call DIPOLE_RANKINE_NN(xf, nq, a, 4, nq, 0d0, enn)
    call check(error, abs(enn) < huge(1d0) .and. enn == enn, "finite exactly on an edge")
  end subroutine test_near_edge

  !> The image part (SV) is the field normal derivative of Omega at the mirrored point: differences of Omega(xbar(x)).
  subroutine test_image(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: a(4,3), nq(3), xf(3), nf(3), e0, e1, em, fd, h
    call skew_quad(a, nq)
    xf = [0.4d0, 0.9d0, -2.1d0]
    nf = [0.3d0, -0.2d0, 0.93d0]
    nf = nf/norm2(nf)
    call DIPOLE_RANKINE_NN(xf, nf, a, 4, nq, 0d0, e0)
    call DIPOLE_RANKINE_NN(xf, nf, a, 4, nq, 1d0, e1)
    call DIPOLE_RANKINE_NN(xf, nf, a, 4, nq, -1d0, em)
    h = 1d-5
    fd = (DIPOLE_SOLID_ANGLE(mirror(xf + h*nf), a, 4, nq) - DIPOLE_SOLID_ANGLE(mirror(xf - h*nf), a, 4, nq))/(2d0*h)
    call check(error, abs((e1 - e0) - fd) <= 1d-6*max(1d0, abs(fd)), "image term vs differences")
    if (allocated(error)) return
    call check(error, abs((e0 - em) - fd) <= 1d-6*max(1d0, abs(fd)), "image term with SV = -1")
  contains
    function mirror(x) result(y)
      real(8), intent(in) :: x(3)
      real(8) :: y(3)
      y = [x(1), x(2), -x(3)]
    end function mirror
  end subroutine test_image

  !> A closed unit cube with outward normals, faces as quads and one face as two triangles: the solid angles sum to -4 pi
  !> inside (every face seen from its -n side) and to 0 outside; the gradient sums to 0 everywhere off the faces.
  subroutine test_cube(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: f(4,3,6), nf(3,6), om, g(3), gs(3), xf(3), t(4,3)
    integer :: k, p
    real(8) :: pts(3,2), want(2)
    call cube(f, nf)
    pts(:,1) = [0.3d0, 0.6d0, -0.45d0]
    pts(:,2) = [1.7d0, 0.2d0, -0.3d0]
    want = [-4d0*pi, 0d0]
    do p = 1, 2
      xf = pts(:,p)
      om = 0d0
      gs = 0d0
      do k = 1, 5
        om = om + DIPOLE_SOLID_ANGLE(xf, f(:,:,k), 4, nf(:,k))
        call DIPOLE_GRAD_OMEGA(xf, f(:,:,k), 4, nf(:,k), g)
        gs = gs + g
      end do
      t = 0d0
      t(1:3,:) = f([1,2,3],:,6)
      om = om + DIPOLE_SOLID_ANGLE(xf, t, 3, nf(:,6))
      call DIPOLE_GRAD_OMEGA(xf, t, 3, nf(:,6), g)
      gs = gs + g
      t(1:3,:) = f([1,3,4],:,6)
      om = om + DIPOLE_SOLID_ANGLE(xf, t, 3, nf(:,6))
      call DIPOLE_GRAD_OMEGA(xf, t, 3, nf(:,6), g)
      gs = gs + g
      call check(error, abs(om - want(p)) <= 1d-12, "solid angle of the closed cube")
      if (allocated(error)) return
      call check(error, norm2(gs) <= 1d-11, "gradient of the closed cube's solid angle")
      if (allocated(error)) return
    end do
  end subroutine test_cube

  subroutine cube(f, nf)
    real(8), intent(out) :: f(4,3,6), nf(3,6)
    real(8), parameter :: z0 = -1d0
    ! bottom, top, x = 0, x = 1, y = 0, y = 1 (node orders mixed on purpose: the routines orient by the normal)
    f(:,:,1) = reshape([0d0,1d0,1d0,0d0, 0d0,0d0,1d0,1d0, z0,z0,z0,z0], [4,3])
    nf(:,1) = [0d0, 0d0, -1d0]
    f(:,:,2) = reshape([0d0,0d0,1d0,1d0, 0d0,1d0,1d0,0d0, 0d0,0d0,0d0,0d0], [4,3])
    nf(:,2) = [0d0, 0d0, 1d0]
    f(:,:,3) = reshape([0d0,0d0,0d0,0d0, 0d0,1d0,1d0,0d0, z0,z0,0d0,0d0], [4,3])
    nf(:,3) = [-1d0, 0d0, 0d0]
    f(:,:,4) = reshape([1d0,1d0,1d0,1d0, 0d0,0d0,1d0,1d0, z0,0d0,0d0,z0], [4,3])
    nf(:,4) = [1d0, 0d0, 0d0]
    f(:,:,5) = reshape([0d0,1d0,1d0,0d0, 0d0,0d0,0d0,0d0, z0,z0,0d0,0d0], [4,3])
    nf(:,5) = [0d0, -1d0, 0d0]
    f(:,:,6) = reshape([0d0,1d0,1d0,0d0, 1d0,1d0,1d0,1d0, z0,z0,0d0,0d0], [4,3])
    nf(:,6) = [0d0, 1d0, 0d0]
  end subroutine cube

  !> Reversing the panel normal reverses Omega and its gradient, whatever the node order.
  subroutine test_flip(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: a(4,3), nq(3), xf(3), e1, e2, o1, o2
    call skew_quad(a, nq)
    xf = [0.4d0, 0.9d0, -2.1d0]
    call DIPOLE_RANKINE_NN(xf, nq, a, 4, nq, 1d0, e1)
    call DIPOLE_RANKINE_NN(xf, nq, a, 4, -nq, 1d0, e2)
    o1 = DIPOLE_SOLID_ANGLE(xf, a, 4, nq)
    o2 = DIPOLE_SOLID_ANGLE(xf, a, 4, -nq)
    call check(error, e1 == -e2 .and. o1 == -o2, "normal flip changes the sign only")
  end subroutine test_flip

end module test_dipole_kernel
