!> Unit tests for the finite-depth Green function (FinGreen3D.f90).
!>  1. G at the self point (R=0, z=z') must equal G at R=1e-4 m (LINTON SR->0 limit).
!>  2. G must be continuous across R/h = 0.02, where LINTON hands over to PIDCOCK.
!>  3. Below R/h = 0.02 and up to nu*h = 50, G and its derivatives must not scatter when R
!>     moves by rounding size. The epsilon algorithm in PIDCOCK did, by up to 1e5 in dG/dR*h^2.
!>  4. There they must match a direct sum of the eigenfunction series at any depth. LINTON
!>     run in metres had errors growing with h (G*h up to 1 at h = 1000), from its absolute
!>     tolerances; it is now run in units of h.
!>  5. The routine choice follows the limits R/h < 0.02 and nu*h <= 50 exactly.
module test_fingreen3d
  use testdrive, only : new_unittest, unittest_type, error_type, check
  use FinGreen3D_Open, only : FINGREEN3D_R, LINTON_SCALED, PIDCOCK
  use direct_sum_green_mod, only : direct_sum_green
  implicit none
  private
  public :: collect_fingreen3d

  integer, parameter :: NK = 200
  real(8), parameter :: GRAV = 9.807d0
  real(8), parameter :: OMEGAS(3) = [0.3d0, 0.8d0, 1.5d0]
  real(8), parameter :: DEPTHS(2) = [70d0, 150d0]

  ! Grid for the small-R tests, in units of h. Depth pairs (zf/h, zp/h) near the surface,
  ! at the same depth, and deeper. nu*h = 0.018 is example 5 at 0.05 rad/s.
  real(8), parameter :: SR_DEPTHS(3) = [10d0, 100d0, 1000d0]
  real(8), parameter :: SR_NUH(4) = [0.018d0, 1d0, 20d0, 50d0]
  real(8), parameter :: SR_RH(5) = [0.0006d0, 0.002d0, 0.005d0, 0.01d0, 0.015d0]
  real(8), parameter :: SR_ZF(6) = [-0.003d0, -0.01d0, -0.03d0, -0.1d0, -0.3d0, -0.06d0]
  real(8), parameter :: SR_ZP(6) = [-0.03d0, -0.003d0, -0.03d0, -0.2d0, -0.5d0, -0.06d0]
  ! Irregular-frequency lid points lie on z = 0.
  real(8), parameter :: LID_ZF(2) = [0d0, 0d0]
  real(8), parameter :: LID_ZP(2) = [-0.03d0, 0d0]

contains

  subroutine collect_fingreen3d(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [ &
      new_unittest("selfterm_continuity", test_selfterm_continuity), &
      new_unittest("linton_pidcock_continuity", test_linton_pidcock_continuity), &
      new_unittest("small_r_no_scatter", test_small_r_no_scatter), &
      new_unittest("small_r_matches_direct_sum", test_small_r_matches_direct_sum), &
      new_unittest("linton_limit", test_linton_limit) ]
  end subroutine collect_fingreen3d

  subroutine test_selfterm_continuity(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: WVN(NK), V, W, H
    complex(8) :: g0(3), g1(3)
    integer :: i, j
    external :: DISPERSION
    do j = 1, size(DEPTHS)
      do i = 1, size(OMEGAS)
        H = DEPTHS(j); W = OMEGAS(i); V = W*W/GRAV
        call DISPERSION(WVN, NK, W, H)
        call FINGREEN3D_R(0d0,  -20d0, -20d0, V, WVN, NK, H, g0, 1)
        call FINGREEN3D_R(1d-4, -20d0, -20d0, V, WVN, NK, H, g1, 1)
        call check(error, abs(g0(1) - g1(1)) < 1d-6*abs(g1(1)), &
                   "LINTON self-term (R=0) not continuous with R=1e-4")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_selfterm_continuity

  subroutine test_linton_pidcock_continuity(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: WVN(NK), V, W, H
    complex(8) :: gl(3), gp(3)
    integer :: i, j
    external :: DISPERSION
    do j = 1, size(DEPTHS)
      do i = 1, size(OMEGAS)
        H = DEPTHS(j); W = OMEGAS(i); V = W*W/GRAV
        call DISPERSION(WVN, NK, W, H)
        call FINGREEN3D_R(0.99999d0*2d-2*H, -20d0, -24d0, V, WVN, NK, H, gl, 1)
        call FINGREEN3D_R(1.00001d0*2d-2*H, -20d0, -24d0, V, WVN, NK, H, gp, 1)
        call check(error, abs(gl(1) - gp(1)) < 1d-2*abs(gp(1)), &
                   "G jumps across R/h = 0.02 (LINTON -> PIDCOCK)")
        if (allocated(error)) return
      end do
    end do
  end subroutine test_linton_pidcock_continuity

  !> Set V and the wave numbers for a given nu*h and depth.
  subroutine set_frequency(nuh, H, V, WVN)
    real(8), intent(in) :: nuh, H
    real(8), intent(out) :: V, WVN(NK)
    external :: DISPERSION
    V = nuh/H
    call DISPERSION(WVN, NK, sqrt(V*GRAV), H)
  end subroutine set_frequency

  subroutine test_small_r_no_scatter(error)
    type(error_type), allocatable, intent(out) :: error
    integer, parameter :: NS = 20
    real(8) :: WVN(NK), V, H, R, s(4, NS), spread, zf(8), zp(8)
    complex(8) :: g(3)
    integer :: ih, iv, ir, iz, k
    character(len=160) :: msg
    zf = [SR_ZF, LID_ZF]; zp = [SR_ZP, LID_ZP]
    do ih = 1, size(SR_DEPTHS)
      H = SR_DEPTHS(ih)
      do iv = 1, size(SR_NUH)
        call set_frequency(SR_NUH(iv), H, V, WVN)
        do ir = 1, size(SR_RH)
          do iz = 1, size(zf)
            do k = 1, NS
              R = SR_RH(ir)*H*(1d0 + k*2d-11)
              call FINGREEN3D_R(R, zf(iz)*H, zp(iz)*H, V, WVN, NK, H, g, 1)
              s(:, k) = [real(g(2)), aimag(g(2)), real(g(3)), aimag(g(3))]*H**2
            end do
            spread = maxval(maxval(s, dim=2) - minval(s, dim=2))
            write(msg, '(A,ES9.2,A,F6.1,A,F7.4,A,2F7.3)') "dG/dR or dG/dz *h^2 scatters by ", &
              spread, " at h = ", H, ", R/h = ", SR_RH(ir), ", z/h = ", zf(iz), zp(iz)
            call check(error, spread < 1d0, trim(msg))
            if (allocated(error)) return
          end do
        end do
      end do
    end do
  end subroutine test_small_r_no_scatter

  subroutine test_small_r_matches_direct_sum(error)
    type(error_type), allocatable, intent(out) :: error
    real(8) :: WVN(NK), V, H, R, zf(8), zp(8), e(3)
    complex(8) :: g(3), ref(3)
    integer :: ih, iv, ir, iz
    character(len=200) :: msg
    zf = [SR_ZF, LID_ZF]; zp = [SR_ZP, LID_ZP]
    do ih = 1, size(SR_DEPTHS)
      H = SR_DEPTHS(ih)
      do iv = 1, size(SR_NUH)
        call set_frequency(SR_NUH(iv), H, V, WVN)
        do ir = 1, size(SR_RH)
          R = SR_RH(ir)*H
          do iz = 1, size(zf)
            call FINGREEN3D_R(R, zf(iz)*H, zp(iz)*H, V, WVN, NK, H, g, 1)
            call direct_sum_green(R, zf(iz)*H, zp(iz)*H, V, WVN(1), H, ref)
            e = [abs(g(1) - ref(1))*H, abs(g(2) - ref(2))*H**2, abs(g(3) - ref(3))*H**2]
            write(msg, '(A,3ES9.2,A,F6.1,A,F6.3,A,F7.4,A,2F7.3)') "error against direct sum (G*h, dG/dR*h^2, dG/dz*h^2) ", &
              e, " at h = ", H, ", nu*h = ", SR_NUH(iv), ", R/h = ", SR_RH(ir), ", z/h = ", zf(iz), zp(iz)
            ! Largest errors measured: 2.6e-5, 4e-4 and 1.1e-3 (at nu*h = 50), so the limits
            ! leave a margin of about 10. PIDCOCK misses them by up to 1e3.
            call check(error, e(1) < 1d-3 .and. e(2) < 5d-3 .and. e(3) < 1d-2, trim(msg))
            if (allocated(error)) return
          end do
        end do
      end do
    end do
  end subroutine test_small_r_matches_direct_sum

  subroutine test_linton_limit(error)
    type(error_type), allocatable, intent(out) :: error
    real(8), parameter :: H = 100d0, ZF = -1d0, ZP = -3d0
    real(8), parameter :: NUH(6) = [49.9d0, 50.1d0, 60d0, 60d0, 20d0, 20d0]
    real(8), parameter :: RH(6) = [0.01d0, 0.01d0, 0.01d0, 0.0004d0, 0.0199d0, 0.0201d0]
    logical, parameter :: WANT_LINTON(6) = [.true., .false., .false., .true., .true., .false.]
    real(8) :: WVN(NK), V, dl, dp
    complex(8) :: g(3), wl(3), wp(3)
    integer :: i
    character(len=100) :: msg
    ! The result must match the expected routine to rounding (inlining can change the last
    ! bits) and differ from the other one by much more than that.
    do i = 1, size(NUH)
      call set_frequency(NUH(i), H, V, WVN)
      call FINGREEN3D_R(RH(i)*H, ZF, ZP, V, WVN, NK, H, g, 1)
      call LINTON_SCALED(RH(i)*H, ZF, ZP, V, WVN, NK, H, wl)
      call PIDCOCK(RH(i)*H, ZF, ZP, V, WVN(1), WVN, NK, H, wp)
      dl = maxval(abs(g - wl)/abs(wl)); dp = maxval(abs(g - wp)/abs(wp))
      write(msg, '(A,F5.1,A,F7.4,A,L1,A,2ES9.2)') "wrong routine at nu*h = ", NUH(i), ", R/h = ", RH(i), &
        ", LINTON expected: ", WANT_LINTON(i), ", differences ", dl, dp
      if (WANT_LINTON(i)) then
        call check(error, dl < 1d-10 .and. dp > 1d-6, trim(msg))
      else
        call check(error, dp < 1d-10 .and. dl > 1d-6, trim(msg))
      end if
      if (allocated(error)) return
    end do
  end subroutine test_linton_limit

end module test_fingreen3d
