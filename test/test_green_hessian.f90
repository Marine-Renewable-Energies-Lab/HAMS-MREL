!> Mixed second derivative of the wave part of the Green function, GNN = d2 G_w / dn_F dn_Q (DipoleIntgr, GREEN_WAVE_NN,
!> central differences), and the field gradient by symmetry (GREEN_WAVE_FIELD_GRAD), against independent references:
!> deep water, the integral form of Wehausen and Laitone (1960, Sec. 13) by quadrature; finite depth, the eigenfunction
!> series summed directly, as test/direct_sum_green.f90 (John 1950). G_w is G minus 1/r and 1/r1 (design note Sec. 1.2).
!> Each reference is first checked on G_w itself, which fixes the sign of the imaginary part to that of HAMS.
!> The tolerances are those of the Green function routines themselves, not of the differences. Against these references
!> (measured 2026-11-24) the source gradient of INFGREEN3D is off by 1e-3 to 7e-3 and its value by up to 2.2e-3 (the
!> deep-water approximation); dG/dzeta of FINGREEN3D is off by 1.2e-4 at R/h = 0.5 and 2e-6 at R/h = 1, its value by
!> 7e-7. The error of the differences alone is checked against Richardson extrapolation (hessian_fd_step_error).
module test_green_hessian

  use testdrive, only : new_unittest, unittest_type, error_type, check
  use DipoleIntgr, only : GREEN_WAVE_GRAD, GREEN_WAVE_FIELD_GRAD, GREEN_WAVE_NN
  use FinGrnExtSubs, only : IK01A, JY01B
  use Const_mod, only : GRAV => G   ! the g of DISPERSION

  implicit none
  private
  public :: collect_green_hessian

  integer, parameter :: NK = 200
  real(8), parameter :: PI = 3.141592653589793d0
  complex(8), parameter :: CI = (0d0, 1d0)

contains

  subroutine collect_green_hessian(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("hessian_deep_fd_vs_series", test_deep), &
      new_unittest("hessian_finite_eigen", test_finite), &
      new_unittest("hessian_symmetric", test_symmetric), &
      new_unittest("hessian_fd_step_error", test_step), &
      new_unittest("field_gradient_by_symmetry", test_field_grad) &
    ]
  end subroutine collect_green_hessian

  !> GNN from the derivatives of G_w(R, z, zeta) (z field, zeta source, R the horizontal distance), with e the unit
  !> horizontal vector from source to field [derived: chain rule; dR/dx_h = e, dR/dxi_h = -e, de/dx_h = (I - e e^T)/R].
  function gnn_from(gr, grr, grz, grq, gzq, e, rr, nf, nq) result(g)
    complex(8), intent(in) :: gr, grr, grz, grq, gzq
    real(8), intent(in) :: e(2), rr, nf(3), nq(3)
    complex(8) :: g
    real(8) :: ef, eq
    ef = e(1)*nf(1) + e(2)*nf(2)
    eq = e(1)*nq(1) + e(2)*nq(2)
    g = -(grr*ef*eq + gr/rr*(nf(1)*nq(1) + nf(2)*nq(2) - ef*eq)) - grz*nf(3)*eq + grq*ef*nq(3) + gzq*nf(3)*nq(3)
  end function gnn_from

  !> The same derivatives for 1/r (sgn = -1, r with z - zeta) and 1/r1 (sgn = +1, z + zeta), to remove them.
  subroutine rankine_derivs(rr, z, q, sgn, g, gr, grr, grz, grq, gzq)
    real(8), intent(in) :: rr, z, q, sgn
    real(8), intent(out) :: g, gr, grr, grz, grq, gzq
    real(8) :: w, r
    w = z + sgn*q
    r = sqrt(rr*rr + w*w)
    g = 1d0/r
    gr = -rr/r**3
    grr = -1d0/r**3 + 3d0*rr*rr/r**5
    grz = 3d0*rr*w/r**5
    grq = sgn*grz
    gzq = sgn*(-1d0/r**3 + 3d0*w*w/r**5)
  end subroutine rankine_derivs

  !> Nodes and weights of the n-point Gauss-Legendre rule on (-1, 1).
  subroutine gauss_legendre(n, x, w)
    integer, intent(in) :: n
    real(8), intent(out) :: x(n), w(n)
    integer :: i, j, it
    real(8) :: z, p0, p1, p2, dp
    do i = 1, n
      z = cos(PI*(i - 0.25d0)/(n + 0.5d0))
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

  !> Deep water: G_w = 2 nu PV int_0^inf e^{kZ} J0(kR)/(k - nu) dk + s 2 pi i nu e^{nu Z} J0(nu R), Z = z + zeta < 0,
  !> and its derivatives (Wehausen and Laitone 1960, Sec. 13, with (k + nu)/(k - nu) = 1 + 2 nu/(k - nu) and
  !> int e^{kZ} J0 dk = 1/r1). d/dZ gives k, d/dR J0(kR) = -k J1, d2/dR2 J0(kR) = -k^2 (J0 - J1/(kR)).
  !> The PV integral: (f(k) - f(nu))/(k - nu) on (0, 2 nu), then f/(k - nu) to k Z = -45; composite 20-point Gauss.
  !> d = [G_w, G_R, G_RR, G_RZ]; G_ZZ follows from Laplace in the caller.
  subroutine deep_reference(nu, rr, zz, s, d)
    real(8), intent(in) :: nu, rr, zz, s
    complex(8), intent(out) :: d(4)
    real(8) :: gx(20), gw(20), k, a, b, kmax, f(4), f0(4), acc(4), h
    integer :: np, ip, i
    call gauss_legendre(20, gx, gw)
    call kernel(nu, f0)
    acc = 0d0
    np = max(40, ceiling(2d0*nu*max(rr, abs(zz))*2d0))
    h = 2d0*nu/np
    do ip = 1, np
      a = (ip - 1)*h
      b = ip*h
      do i = 1, 20
        k = 0.5d0*(a + b) + 0.5d0*(b - a)*gx(i)
        call kernel(k, f)
        acc = acc + 0.5d0*(b - a)*gw(i)*(f - f0)/(k - nu)
      end do
    end do
    kmax = 45d0/abs(zz)
    np = max(40, ceiling((kmax - 2d0*nu)*max(rr, abs(zz))*2d0))
    h = (kmax - 2d0*nu)/np
    do ip = 1, np
      a = 2d0*nu + (ip - 1)*h
      b = a + h
      do i = 1, 20
        k = 0.5d0*(a + b) + 0.5d0*(b - a)*gx(i)
        call kernel(k, f)
        acc = acc + 0.5d0*(b - a)*gw(i)*f/(k - nu)
      end do
    end do
    call kernel(nu, f)
    d = 2d0*nu*acc + s*2d0*PI*CI*nu*f
  contains
    subroutine kernel(kk, fk)
      real(8), intent(in) :: kk
      real(8), intent(out) :: fk(4)
      real(8) :: ez, j0, j1, x
      ez = exp(kk*zz)
      x = kk*rr
      j0 = bessel_j0(x)
      j1 = bessel_j1(x)
      fk(1) = ez*j0
      fk(2) = -kk*ez*j1
      if (x > 1d-12) then
        fk(3) = -kk*kk*ez*(j0 - j1/x)
      else
        fk(3) = -0.5d0*kk*kk*ez
      end if
      fk(4) = -kk*kk*ez*j1
    end subroutine kernel
  end subroutine deep_reference

  subroutine unit_normal(n, a, b, c)
    real(8), intent(out) :: n(3)
    real(8), intent(in) :: a, b, c
    n = [a, b, c]
    n = n/norm2(n)
  end subroutine unit_normal

  !> Deep water, nu R from 0.1 to 30, both points at least 1/nu deep: GREEN_WAVE_NN against the integral form, to the
  !> accuracy of INFGREEN3D's own gradient (5e-3).
  subroutine test_deep(error)
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: nus(2) = [0.4d0, 1.5d0], krs(4) = [0.1d0, 1d0, 5d0, 30d0]
    real(8) :: nu, rr, xf(3), xq(3), nf(3), nq(3), e(2), s, wvn(NK), err
    complex(8) :: d(4), gw(4), gnn, ref
    character(len=160) :: msg
    integer :: iv, ir
    wvn = 0d0
    call unit_normal(nf, 0.3d0, -0.4d0, 0.87d0)
    call unit_normal(nq, -0.5d0, 0.2d0, 0.84d0)
    do iv = 1, 2
      nu = nus(iv)
      do ir = 1, 4
        rr = krs(ir)/nu
        e = [0.6d0, 0.8d0]
        xq = [1d0, 2d0, -1.3d0/nu]
        xf = [xq(1) + rr*e(1), xq(2) + rr*e(2), -1.1d0/nu]
        ! G_w itself, which fixes s, the sign of the imaginary part
        call GREEN_WAVE_GRAD(xf, xq, nu, wvn, NK, -1d0, gw)
        s = sign(1d0, aimag(gw(1))*bessel_j0(nu*rr))
        call deep_reference(nu, rr, xf(3) + xq(3), s, d)
        err = abs(gw(1) - d(1))/abs(d(1))
        write(msg, '(a,f5.2,a,f6.2,a,es9.2)') 'deep G_w, nu = ', nu, ', nu R = ', krs(ir), ': rel err ', err
        print '(2x,a)', trim(msg)
        call check(error, err < 3d-3, trim(msg))
        if (allocated(error)) return
        ! G_ZZ = -(G_RR + G_R/R) (Laplace); G_Rz = G_Rzeta = G_RZ; G_zzeta = G_ZZ
        ref = gnn_from(d(2), d(3), d(4), d(4), -(d(3) + d(2)/rr), e, rr, nf, nq)
        call GREEN_WAVE_NN(xf, nf, xq, nq, nu, wvn, NK, -1d0, gnn)
        err = abs(gnn - ref)/abs(ref)
        write(msg, '(a,f5.2,a,f6.2,a,es9.2)') 'deep GNN, nu = ', nu, ', nu R = ', krs(ir), ': rel err ', err
        print '(2x,a)', trim(msg)
        call check(error, err < 5d-3, trim(msg))
        if (allocated(error)) return
      end do
    end do
  end subroutine test_deep

  !> Finite depth, R/h from 0.5 to 3: the eigenfunction series (as direct_sum_green) differentiated term by term, minus the
  !> Rankine terms. Mode 0: pi/N0 cosh k0(z+h) cosh k0(zeta+h) (i J0 - Y0)(k0 R), N0 = h/2 (1 + sinh(2 k0 h)/(2 k0 h));
  !> modes m >= 1: 4/(h (1 + sin(2 u h)/(2 u h))) cos u(z+h) cos u(zeta+h) K0(u R), u tan(u h) = -nu.
  !> Tolerance: 8 times the error of FINGREEN3D's own dG/dzeta at each R/h.
  subroutine test_finite(error)
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: hh = 30d0, nuhs(2) = [0.5d0, 4d0], rhs(3) = [0.5d0, 1d0, 3d0]
    real(8), parameter :: tol_fin(3) = [1d-3, 1d-4, 1d-6]
    real(8) :: nu, wvn(NK), rr, z, q, e(2), xf(3), xq(3), nf(3), nq(3), err
    real(8) :: g, gr, grr, grz, grq, gzq
    complex(8) :: c(6), gnn, ref, gw(4)
    character(len=160) :: msg
    integer :: iv, ir
    call unit_normal(nf, 0.3d0, -0.4d0, 0.87d0)
    call unit_normal(nq, -0.5d0, 0.2d0, 0.84d0)
    do iv = 1, 2
      nu = nuhs(iv)/hh
      call DISPERSION_LOCAL(wvn, nu, hh)
      do ir = 1, 3
        rr = rhs(ir)*hh
        e = [0.8d0, -0.6d0]
        z = -8d0
        q = -21d0
        xq = [0d0, 0d0, q]
        xf = [rr*e(1), rr*e(2), z]
        call finite_reference(rr, z, q, nu, wvn(1), hh, c)
        call GREEN_WAVE_GRAD(xf, xq, nu, wvn, NK, hh, gw)
        err = abs(gw(1) - c(1))/abs(c(1))
        write(msg, '(a,f5.2,a,f5.2,a,es9.2)') 'finite G_w, nu h = ', nuhs(iv), ', R/h = ', rhs(ir), ': rel err ', err
        print '(2x,a)', trim(msg)
        call check(error, err < 1d-6, trim(msg))
        if (allocated(error)) return
        ref = gnn_from(c(2), c(3), c(4), c(5), c(6), e, rr, nf, nq)
        call GREEN_WAVE_NN(xf, nf, xq, nq, nu, wvn, NK, hh, gnn)
        err = abs(gnn - ref)/abs(ref)
        write(msg, '(a,f5.2,a,f5.2,a,es9.2)') 'finite GNN, nu h = ', nuhs(iv), ', R/h = ', rhs(ir), ': rel err ', err
        print '(2x,a)', trim(msg)
        call check(error, err < tol_fin(ir), trim(msg))
        if (allocated(error)) return
      end do
    end do
  end subroutine test_finite

  !> The wave part from the series: eigen_derivs minus 1/r and 1/r1, c = [G_w, G_R, G_RR, G_Rz, G_Rzeta, G_zzeta].
  subroutine finite_reference(rr, z, q, nu, k0, hh, c)
    real(8), intent(in) :: rr, z, q, nu, k0, hh
    complex(8), intent(out) :: c(6)
    real(8) :: g, gr, grr, grz, grq, gzq
    call eigen_derivs(rr, z, q, nu, k0, hh, c)
    call rankine_derivs(rr, z, q, -1d0, g, gr, grr, grz, grq, gzq)
    c = c - [g, gr, grr, grz, grq, gzq]
    call rankine_derivs(rr, z, q, 1d0, g, gr, grr, grz, grq, gzq)
    c = c - [g, gr, grr, grz, grq, gzq]
  end subroutine finite_reference

  !> The roots of the dispersion relation from HAMS's DISPERSION (Const_mod's G, as in test_fingreen3d).
  subroutine DISPERSION_LOCAL(wvn, nu, hh)
    real(8), intent(out) :: wvn(NK)
    real(8), intent(in) :: nu, hh
    external :: DISPERSION
    call DISPERSION(wvn, NK, sqrt(nu*GRAV), hh)
  end subroutine DISPERSION_LOCAL

  !> Full G (eigen series) and its derivatives in R, z (field) and zeta (source): c = [G, G_R, G_RR, G_Rz, G_Rzeta, G_zzeta].
  !> d/dx (i J0 - Y0) = -(i J1 - Y1), d/dx (i J1 - Y1) = (i J0 - Y0) - (i J1 - Y1)/x; K0' = -K1, K0'' = K0 + K1/x.
  subroutine eigen_derivs(rr, z, q, nu, k0, hh, c)
    real(8), intent(in) :: rr, z, q, nu, k0, hh
    complex(8), intent(out) :: c(6)
    real(8) :: j0, j1, y0, y1, bk0, bk1, u, tw, pm, x
    complex(8) :: nm0, a, h0, h1, h1p
    integer :: m
    nm0 = hh/2d0*(1d0 + sinh(2d0*k0*hh)/(2d0*k0*hh))
    x = k0*rr
    call JY01B(x, j0, j1, y0, y1)
    h0 = CI*j0 - y0
    h1 = CI*j1 - y1
    h1p = h0 - h1/x
    a = PI/nm0
    c(1) = a*cosh(k0*(z + hh))*cosh(k0*(q + hh))*h0
    c(2) = -a*k0*cosh(k0*(z + hh))*cosh(k0*(q + hh))*h1
    c(3) = -a*k0*k0*cosh(k0*(z + hh))*cosh(k0*(q + hh))*h1p
    c(4) = -a*k0*k0*sinh(k0*(z + hh))*cosh(k0*(q + hh))*h1
    c(5) = -a*k0*k0*cosh(k0*(z + hh))*sinh(k0*(q + hh))*h1
    c(6) = a*k0*k0*sinh(k0*(z + hh))*sinh(k0*(q + hh))*h0
    m = 0
    do
      m = m + 1
      u = um_root(m, nu*hh, hh)
      if (u*rr > 40d0) exit
      tw = 2d0*u*hh
      pm = 4d0/(1d0 + sin(tw)/tw)/hh
      call IK01A(u*rr, bk0, bk1)
      c(1) = c(1) + pm*bk0*cos(u*(z + hh))*cos(u*(q + hh))
      c(2) = c(2) - pm*u*bk1*cos(u*(z + hh))*cos(u*(q + hh))
      c(3) = c(3) + pm*u*u*(bk0 + bk1/(u*rr))*cos(u*(z + hh))*cos(u*(q + hh))
      c(4) = c(4) + pm*u*u*bk1*sin(u*(z + hh))*cos(u*(q + hh))
      c(5) = c(5) + pm*u*u*bk1*cos(u*(z + hh))*sin(u*(q + hh))
      c(6) = c(6) + pm*u*u*bk0*sin(u*(z + hh))*sin(u*(q + hh))
    end do
  end subroutine eigen_derivs

  !> m-th positive root u of u tan(u h) = -nu (as in direct_sum_green).
  real(8) function um_root(m, vh, h)
    integer, intent(in) :: m
    real(8), intent(in) :: vh, h
    real(8) :: d, f, y
    integer :: it
    d = atan2(vh, m*PI)
    do it = 1, 60
      y = m*PI - d
      f = (d - atan2(vh, y))/(1d0 - vh/(y*y + vh*vh))
      d = d - f
      if (abs(f) < 1d-15) exit
    end do
    um_root = (m*PI - d)/h
  end function um_root

  !> GNN(F, Q) = GNN(Q, F) with the normals swapped: deep water to 1e-3 (INFGREEN3D's gradient is not the exact
  !> derivative of its value, see the header), finite depth to 1e-6.
  subroutine test_symmetric(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: xf(3), xq(3), nf(3), nq(3), wvn(NK), nu, hh
    complex(8) :: a, b
    integer :: ih
    call unit_normal(nf, 0.3d0, -0.4d0, 0.87d0)
    call unit_normal(nq, -0.5d0, 0.2d0, 0.84d0)
    xf = [3d0, -1d0, -4d0]
    xq = [-2d0, 1.5d0, -9d0]
    nu = 0.3d0
    do ih = 1, 2
      wvn = 0d0
      hh = -1d0
      if (ih == 2) then
        hh = 40d0
        call DISPERSION_LOCAL(wvn, nu, hh)
      end if
      call GREEN_WAVE_NN(xf, nf, xq, nq, nu, wvn, NK, hh, a)
      call GREEN_WAVE_NN(xq, nq, xf, nf, nu, wvn, NK, hh, b)
      print '(2x,a,i2,a,es9.2)', 'symmetry, depth case', ih, ': rel diff ', abs(a - b)/abs(a)
      call check(error, abs(a - b) <= merge(1d-3, 1d-6, ih == 1)*abs(a), "GNN symmetric in the two points")
      if (allocated(error)) return
    end do
  end subroutine test_symmetric

  !> The error of the central differences alone: Richardson extrapolation of the same differences with steps delta and
  !> 2 delta, (4 D(delta) - D(2 delta))/3, removes the delta^2 term; GNN must agree with it to 1e-7 (plan Task 2: the
  !> differences are kept while their own error is below 1e-6). Deep and finite depth, a near pair and a far pair.
  subroutine test_step(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: xf(3), xq(3), nf(3), nq(3), wvn(NK), nu, hh, r1, k0, dl
    complex(8) :: gnn, d2, ext, gp(4), gm(4)
    integer :: ih, ip
    character(len=120) :: msg
    call unit_normal(nf, 0.3d0, -0.4d0, 0.87d0)
    call unit_normal(nq, -0.5d0, 0.2d0, 0.84d0)
    do ih = 1, 2
      nu = 0.3d0
      wvn = 0d0
      hh = -1d0
      if (ih == 2) then
        hh = 40d0
        call DISPERSION_LOCAL(wvn, nu, hh)
      end if
      do ip = 1, 2
        xq = [0d0, 0d0, -6d0]
        if (ip == 1) then
          xf = [0.4d0, -0.3d0, -5.5d0]
        else
          xf = [25d0, 10d0, -3d0]
        end if
        call GREEN_WAVE_NN(xf, nf, xq, nq, nu, wvn, NK, hh, gnn)
        r1 = sqrt(sum((xq(1:2) - xf(1:2))**2) + (xq(3) + xf(3))**2)
        k0 = merge(nu, wvn(1), hh < 0d0)
        dl = max(1d-4*min(r1, 1d0/k0), 1d-5)
        call grad_at(xf + 2d0*dl*nf, gp)
        call grad_at(xf - 2d0*dl*nf, gm)
        d2 = sum((gp(2:4) - gm(2:4))*nq)/(4d0*dl)
        ext = (4d0*gnn - d2)/3d0
        write(msg, '(a,i2,a,i2,a,es9.2)') 'step error, depth case', ih, ', pair', ip, ': ', abs(gnn - ext)/abs(ext)
        print '(2x,a)', trim(msg)
        call check(error, abs(gnn - ext) <= 1d-7*abs(ext), trim(msg))
        if (allocated(error)) return
      end do
    end do
  contains
    subroutine grad_at(x, g4)
      real(8), intent(in) :: x(3)
      complex(8), intent(out) :: g4(4)
      if (hh < 0d0) then
        call GREEN_WAVE_GRAD(x, xq, nu, wvn, NK, hh, g4)
      else
        call GREEN_WAVE_GRAD(x, xq, nu, wvn, NK, hh, g4, RRD=sqrt(sum((xq(1:2) - xf(1:2))**2)), FIELD=.true.)
      end if
    end subroutine grad_at
  end subroutine test_step

  !> The field gradient by symmetry, against the references: d/dx_h = G_R e (deep water nu R = 1 to 5e-3, finite depth
  !> R/h = 1 to 1e-6); and it is exactly the source gradient with the points exchanged.
  subroutine test_field_grad(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: xf(3), xq(3), wvn(NK), nu, hh, rr, e(2), s
    complex(8) :: gf(4), d(4), c(6), gw(4)
    e = [0.6d0, 0.8d0]
    nu = 0.4d0
    rr = 1d0/nu
    wvn = 0d0
    xq = [1d0, 2d0, -1.3d0/nu]
    xf = [xq(1) + rr*e(1), xq(2) + rr*e(2), -1.1d0/nu]
    call GREEN_WAVE_GRAD(xf, xq, nu, wvn, NK, -1d0, gw)
    s = sign(1d0, aimag(gw(1))*bessel_j0(nu*rr))
    call deep_reference(nu, rr, xf(3) + xq(3), s, d)
    call GREEN_WAVE_FIELD_GRAD(xf, xq, nu, wvn, NK, -1d0, gf)
    call check(error, abs(gf(2) - d(2)*e(1)) <= 5d-3*abs(d(2)) .and. abs(gf(3) - d(2)*e(2)) <= 5d-3*abs(d(2)), &
               "deep field gradient, horizontal")
    if (allocated(error)) return
    hh = 30d0
    nu = 0.5d0/hh
    call DISPERSION_LOCAL(wvn, nu, hh)
    rr = hh
    e = [0.8d0, -0.6d0]
    xq = [0d0, 0d0, -21d0]
    xf = [rr*e(1), rr*e(2), -8d0]
    call finite_reference(rr, xf(3), xq(3), nu, wvn(1), hh, c)
    call GREEN_WAVE_FIELD_GRAD(xf, xq, nu, wvn, NK, hh, gf)
    call check(error, abs(gf(2) - c(2)*e(1)) <= 1d-6*abs(c(2)) .and. abs(gf(3) - c(2)*e(2)) <= 1d-6*abs(c(2)), &
               "finite-depth field gradient, horizontal")
    if (allocated(error)) return
    call GREEN_WAVE_GRAD(xq, xf, nu, wvn, NK, hh, gw)
    call check(error, all(gf == gw), "field gradient is the source gradient with the points exchanged")
  end subroutine test_field_grad

end module test_green_hessian
