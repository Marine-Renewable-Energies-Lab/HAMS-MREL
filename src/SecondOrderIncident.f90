!  ------------------------------------------------------------------------------------------------------
!
!    Program HAMS-MREL for the diffraction and radiation of waves
!    by 3D structures.
!
!  License:
!
!    This routine is part of HAMS-MREL.
!
!    HAMS-MREL is a free software framework: you can redistribute it and/or modify it
!    under the terms of the Apache License, Version 2.0 (the "License"); you may
!    not use these subroutines except in compliance with the License. The software
!    is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND.
!
!    You should have received a copy of the Apache License, Version 2.0, along with
!    HAMS-MREL. If not, see <http://www.apache.org/licenses/LICENSE-2.0>.
!
!  Code Original Author:
!
!    Vaibhav Raghavan
!
!  ------------------------------------------------------------------------------------------------------
MODULE SecondOrderIncident
   ! Second-order incident potential phi_I2 of a pair of unit-amplitude first-order waves (Molin 1979; Sharma and Dean 1981),
   ! time factor e^{-i omega t}, at the difference frequency Omega = omega_i - omega_j (KERNEL_DIFF) or the sum frequency
   ! Omega = omega_i + omega_j (KERNEL_SUM). It is the solution of
   !   -Omega^2 phi + g phi_z = Q  on z = 0,   Q = -d/dt |grad phi1|^2 + (1/g) phi1_t d/dz (phi1_tt + g phi1_z)
   ! with the dependence e^{i K x} Z_K(z), K = k_i e_i - k_j e_j (difference) or k_i e_i + k_j e_j (sum). Q is formed with the
   ! pair kernel of NearFieldDrift (pair_minus, pair_plus), so the factor of the pair is not derived here. In the sum case
   ! the time derivative of the product carries the factor i Omega, and D = -Omega^2 + g K tanh(K h) does not vanish: there
   ! is no diagonal short-cut. For one wave twice (omega_i = omega_j, same heading) the sum potential is the second-order
   ! Stokes wave, -i (3/8) omega cosh(2kh) / sinh^4(kh) for Z normalised to 1 at z = 0 (Dean and Dalrymple 1991), which
   ! vanishes in deep water. Per unit amplitudes the potential is coef e^{i K x} Z_K(z); the pair of waves a_i, a_j gives
   ! a_i conj(a_j) times it for the difference and a_i a_j times it for the sum, as the other quantities of that kernel.

   USE FirstOrderState, ONLY: incident_phi, incident_grad, depth_factor
   USE NearFieldDrift, ONLY: pair_product, KERNEL_DIFF, KERNEL_SUM
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: incident2, incident2_pair, incident2_phi, incident2_grad

   TYPE :: incident2
      COMPLEX(8) :: coef = (0D0, 0D0)   ! Q / D, per unit amplitudes
      REAL(8) :: kvec(2) = 0D0          ! K, horizontal wave vector of the difference or sum
      REAL(8) :: kmag = 0D0             ! |K|
      REAL(8) :: omega = 0D0            ! Omega = omega_i -/+ omega_j (difference/sum)
      REAL(8) :: h = -1D0               ! depth, h < 0: deep water
   END TYPE incident2

CONTAINS

   SUBROUTINE incident2_pair(w1, k1, b1, w2, k2, b2, h, g, inc, kernel)
      ! Waves (omega, wave number, heading in radians) 1 and 2 in depth h (h < 0: deep) with gravity g of the waves. The
      ! exact diagonal (same wave twice) has no second-order potential at the difference frequency: coef = 0, K = 0. The
      ! kernel (KERNEL_DIFF, the default, or KERNEL_SUM) selects the difference or sum frequency; the sum has no diagonal short-cut.
      REAL(8), INTENT(IN) :: w1, k1, b1, w2, k2, b2, h, g
      TYPE(incident2), INTENT(OUT) :: inc
      INTEGER, INTENT(IN), OPTIONAL :: kernel
      REAL(8), PARAMETER :: ORIGIN(3) = 0D0
      COMPLEX(8) :: phi(2), u(3, 2), s(2), q, d
      REAL(8) :: omega(2), k(2), beta(2), dd
      INTEGER :: i, m, kn

      kn = KERNEL_DIFF
      IF (PRESENT(kernel)) kn = kernel
      omega = [w1, w2]
      k = [k1, k2]
      beta = [b1, b2]
      inc%h = h
      IF (kn == KERNEL_SUM) THEN
         inc%omega = w1 + w2
         inc%kvec = [k1*COS(b1) + k2*COS(b2), k1*SIN(b1) + k2*SIN(b2)]
      ELSE
         inc%omega = w1 - w2
         inc%kvec = [k1*COS(b1) - k2*COS(b2), k1*SIN(b1) - k2*SIN(b2)]
      END IF
      inc%kmag = SQRT(inc%kvec(1)**2 + inc%kvec(2)**2)
      IF (kn /= KERNEL_SUM .AND. w1 == w2 .AND. k1 == k2 .AND. b1 == b2) THEN
         inc%kvec = 0D0
         inc%kmag = 0D0
         inc%coef = (0D0, 0D0)
         RETURN
      END IF

      DO i = 1, 2
         phi(i) = incident_phi(ORIGIN, omega(i), k(i), h, 1D0, beta(i), g)
         CALL incident_grad(ORIGIN, omega(i), k(i), h, 1D0, beta(i), g, u(:, i))
         ! phi_zz = k^2 phi for a first-order wave, at any depth
         s(i) = -omega(i)**2*u(3, i) + g*k(i)**2*phi(i)
      END DO
      q = (0D0, 0D0)
      DO m = 1, 3
         q = q + pair_product(kn, u(m, 1), u(m, 1), u(m, 2), u(m, 2))
      END DO
      q = (0D0, 1D0)*inc%omega*q &
          + pair_product(kn, -(0D0, 1D0)*omega(1)*phi(1), s(1), -(0D0, 1D0)*omega(2)*phi(2), s(2))/g
      IF (h < 0D0) THEN
         dd = -inc%omega**2 + g*inc%kmag
      ELSE
         dd = -inc%omega**2 + g*inc%kmag*TANH(inc%kmag*h)
      END IF
      d = CMPLX(dd, 0D0, KIND=8)
      inc%coef = q/d
   END SUBROUTINE incident2_pair

   PURE COMPLEX(8) FUNCTION incident2_phi(inc, x)
      ! Second-order incident potential of the pair, per unit amplitudes, at x(3).
      TYPE(incident2), INTENT(IN) :: inc
      REAL(8), INTENT(IN) :: x(3)
      REAL(8) :: zf, dzf

      CALL depth_factor(inc%kmag, x(3), inc%h, zf, dzf)
      incident2_phi = inc%coef*EXP(CMPLX(0D0, inc%kvec(1)*x(1) + inc%kvec(2)*x(2), KIND=8))*zf
   END FUNCTION incident2_phi

   PURE SUBROUTINE incident2_grad(inc, x, grad)
      ! Gradient of incident2_phi at x(3).
      TYPE(incident2), INTENT(IN) :: inc
      REAL(8), INTENT(IN) :: x(3)
      COMPLEX(8), INTENT(OUT) :: grad(3)
      REAL(8) :: zf, dzf
      COMPLEX(8) :: pref

      CALL depth_factor(inc%kmag, x(3), inc%h, zf, dzf)
      pref = inc%coef*EXP(CMPLX(0D0, inc%kvec(1)*x(1) + inc%kvec(2)*x(2), KIND=8))
      grad(1) = CMPLX(0D0, inc%kvec(1), KIND=8)*pref*zf
      grad(2) = CMPLX(0D0, inc%kvec(2), KIND=8)*pref*zf
      grad(3) = pref*dzf
   END SUBROUTINE incident2_grad

END MODULE SecondOrderIncident
