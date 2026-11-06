!> Reference for the finite-depth Green function: the eigenfunction series summed directly,
!> with no acceleration. Same return convention as FINGREEN3D_R with TAG = 1 (the wave part,
!> with 1/SR and 1/SR1 removed, as in EIGEN). The sum stops where um*R > XCUT, where K0 is
!> about 1e-18. Checked against 30-digit mpmath: G to 3e-11, dG/dR to 1.5e-8.
!> Slow (up to 20000 terms at R/h = 0.0006), so for tests only.
module direct_sum_green_mod
  use FinGrnExtSubs, only : IK01A, JY01B
  implicit none
  private
  public :: direct_sum_green

  real(8), parameter :: XCUT = 40d0, PI = 3.141592653589793d0

contains

  !> m-th positive root um of um*tan(um*h) = -V, by Newton iteration on um*h = m*pi - d.
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

  !> Neumaier compensated addition: s + c holds the sum.
  pure subroutine kadd(s, c, x)
    real(8), intent(inout) :: s, c
    real(8), intent(in) :: x
    real(8) :: t
    t = s + x
    if (abs(s) >= abs(x)) then
      c = c + ((s - t) + x)
    else
      c = c + ((x - t) + s)
    end if
    s = t
  end subroutine kadd

  subroutine direct_sum_green(R, ZF, ZP, V, WK, H, GRN)
    real(8), intent(in) :: R, ZF, ZP, V, WK, H
    complex(8), intent(out) :: GRN(3)
    real(8) :: J0, J1, Y0, Y1, K0, K1, um, tw, pm, cf, cp, sf, sr, sr1
    real(8) :: sg, sgr, sgz, cg, cgr, cgz
    complex(8) :: p0, dp0, nm0
    complex(8), parameter :: CI = (0d0, 1d0)
    integer :: m

    nm0 = H/2d0*(1d0 + sinh(2d0*WK*H)/(2d0*WK*H))
    p0 = PI/nm0*cosh(WK*(ZF + H))*cosh(WK*(ZP + H))
    dp0 = PI/nm0*sinh(WK*(ZF + H))*cosh(WK*(ZP + H))
    call JY01B(WK*R, J0, J1, Y0, Y1)
    GRN(1) = p0*(CI*J0 - Y0)
    GRN(2) = -WK*p0*(CI*J1 - Y1)
    GRN(3) = WK*dp0*(CI*J0 - Y0)

    sg = 0d0; sgr = 0d0; sgz = 0d0; cg = 0d0; cgr = 0d0; cgz = 0d0
    m = 0
    do
      m = m + 1
      um = um_root(m, V*H, H)
      if (um*R > XCUT) exit
      tw = 2d0*um*H
      pm = 4d0/(1d0 + sin(tw)/tw)/H
      call IK01A(um*R, K0, K1)
      cf = cos(um*(ZF + H)); cp = cos(um*(ZP + H)); sf = sin(um*(ZF + H))
      call kadd(sg, cg, pm*K0*cf*cp)
      call kadd(sgr, cgr, -pm*um*K1*cf*cp)
      call kadd(sgz, cgz, -pm*um*K0*sf*cp)
    end do

    sr = sqrt(R**2 + (ZF - ZP)**2)
    sr1 = sqrt(R**2 + (ZF + ZP)**2)
    GRN(1) = GRN(1) + (sg + cg) - (1d0/sr + 1d0/sr1)
    GRN(2) = GRN(2) + (sgr + cgr) + (R/sr**3 + R/sr1**3)
    GRN(3) = GRN(3) + (sgz + cgz) + ((ZF - ZP)/sr**3 + (ZF + ZP)/sr1**3)
  end subroutine direct_sum_green

end module direct_sum_green_mod
