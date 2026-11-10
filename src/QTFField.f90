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
MODULE QTFField
   ! Second-order pressure and free-surface elevation of a pair of first-order waves at a field point: the quadratic part
   ! (products of first-order quantities) plus the second-order incident potential. The second-order diffracted and radiated
   ! potentials are omitted, so near a body the result is partial.
   !
   ! Bernoulli's equation expanded about z = 0 (Dean and Dalrymple 1991), time factor e^{-i omega t}:
   !   pressure at a fixed point below z = 0:  p2 = -rho [phi2_t + 1/2 |grad phi1|^2]
   !   elevation on z = 0:  eta2 = -(1/g) [phi2_t + 1/2 |grad phi1|^2 + eta1 phi1_zt],
   !   with eta1 = (i omega/g) phi1 and phi1_zt = -i omega phi1_z.
   ! For the ordered pair (i, j) and the kernel P of NearFieldDrift (pair_minus: difference, pair_plus: sum), per unit
   ! amplitudes:
   !   quad_p   = -(rho/2) sum_m P(u_im, u_im, u_jm, u_jm)
   !   quad_eta = -(1/g) [ 1/2 sum_m P(u_im, u_im, u_jm, u_jm) + P(eta_i, -i omega_i u_i3, eta_j, -i omega_j u_j3) ]
   !   pot_p    = i rho Omega phi2,   pot_eta = (i Omega/g) phi2,
   ! with Omega = omega_i -/+ omega_j and phi2 the second-order incident potential of the pair (SecondOrderIncident). The
   ! physical e^{-i Omega t} amplitude of two waves a_i /= a_j is 2 a_i a_j (quad + pot); for the sum with i = j it is
   ! a_i^2 (quad + pot).

   USE SecondOrderIncident, ONLY: incident2, incident2_phi
   USE NearFieldDrift, ONLY: pair_product
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: field_pair

CONTAINS

   PURE SUBROUTINE field_pair(kernel, w1, phi1, grad1, w2, phi2, grad2, inc, x, surface, rho, g, quad, pot)
      ! Quadratic part quad and second-order incident part pot of the pressure (surface = .FALSE.) or the elevation
      ! (surface = .TRUE., x(3) = 0) of the pair of unit waves with frequencies w1, w2, first-order potentials phi and
      ! gradients grad at x. inc is the second-order incident potential of the pair, made with the same kernel.
      INTEGER, INTENT(IN) :: kernel
      REAL(8), INTENT(IN) :: w1, w2, x(3), rho, g
      COMPLEX(8), INTENT(IN) :: phi1, grad1(3), phi2, grad2(3)
      TYPE(incident2), INTENT(IN) :: inc
      LOGICAL, INTENT(IN) :: surface
      COMPLEX(8), INTENT(OUT) :: quad, pot
      COMPLEX(8) :: speed, p2, eta1, eta2
      COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
      INTEGER :: m

      speed = (0D0, 0D0)
      DO m = 1, 3
         speed = speed + pair_product(kernel, grad1(m), grad1(m), grad2(m), grad2(m))
      END DO
      p2 = incident2_phi(inc, x)
      IF (surface) THEN
         eta1 = CI*w1/g*phi1
         eta2 = CI*w2/g*phi2
         quad = -(0.5D0*speed + pair_product(kernel, eta1, -CI*w1*grad1(3), eta2, -CI*w2*grad2(3)))/g
         pot = CI*inc%omega/g*p2
      ELSE
         quad = -0.5D0*rho*speed
         pot = CI*rho*inc%omega*p2
      END IF
   END SUBROUTINE field_pair

END MODULE QTFField
