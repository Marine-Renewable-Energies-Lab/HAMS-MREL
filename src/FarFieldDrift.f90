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
MODULE FarFieldDrift
   ! Mean wave drift force and yaw moment from the far-field (Kochin) amplitude of the perturbation potential.
   ! Formulae: docs/superpowers/specs/2026-11-04-second-order-waves-design.md, section 5.3. Time factor e^{-i omega t}.
   ! The amplitudes scale with s%amp; the callers pass states with amp = 1 to get loads per unit amplitude.

   USE SecondOrderTypes, ONLY: panel_mesh, first_order_state, so_options
   USE FirstOrderState, ONLY: depth_factor
   IMPLICIT NONE
   PRIVATE

   PUBLIC :: far_amplitude, group_ratio, far_field_pair

   REAL(8), PARAMETER :: PI = 3.141592653589793D0
   COMPLEX(8), PARAMETER :: CI = (0D0, 1D0)
   INTEGER, PARAMETER :: MAX_AZIMUTH = 2**14

CONTAINS

   PURE REAL(8) FUNCTION group_ratio(k, h)
      ! C_g/C_p: 1/2 in deep water (h < 0), otherwise 1/2 (1 + 2kh/sinh 2kh), written to be safe for large kh.
      REAL(8), INTENT(IN) :: k, h
      REAL(8) :: x

      IF (h < 0D0) THEN
         group_ratio = 0.5D0
      ELSE
         x = 2D0*k*h
         IF (x <= 0D0) THEN
            group_ratio = 1D0
         ELSE
            group_ratio = 0.5D0*(1D0 + 2D0*x*EXP(-x)/(1D0 - EXP(-2D0*x)))
         END IF
      END IF
   END FUNCTION group_ratio

   PURE REAL(8) FUNCTION kochin_c(k, nu, h)
      ! Wave-part strength of the Green function: k in deep water, k^2/((k^2-nu^2)h + nu) otherwise.
      REAL(8), INTENT(IN) :: k, nu, h

      IF (h < 0D0) THEN
         kochin_c = k
      ELSE
         kochin_c = k*k/((k*k - nu*nu)*h + nu)
      END IF
   END FUNCTION kochin_c

   SUBROUTINE far_amplitude(mesh, s, theta, a, da)
      ! Far-field amplitude a(theta) of the outgoing perturbation wave, eta_P ~ a e^{ikr}/sqrt(r), and da = da/dtheta.
      ! Computed from phi_p and dphi_p_dn only.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(first_order_state), INTENT(IN) :: s
      REAL(8), INTENT(IN) :: theta
      COMPLEX(8), INTENT(OUT) :: a, da
      REAL(8) :: c, sn, zf, dzf, xx, yy, zz, nd, dnd
      COMPLEX(8) :: ph, e, dn_e, q, de, ddn_e, hh, dhh, pref
      INTEGER :: i

      c = COS(theta)
      sn = SIN(theta)
      hh = (0D0, 0D0)
      dhh = (0D0, 0D0)
      DO i = 1, mesh%n
         xx = mesh%xc(1, i)
         yy = mesh%xc(2, i)
         zz = mesh%xc(3, i)
         CALL depth_factor(s%k, zz, s%h, zf, dzf)
         ph = EXP(CMPLX(0D0, -s%k*(xx*c + yy*sn), KIND=8))
         nd = mesh%nrm(1, i)*c + mesh%nrm(2, i)*sn
         dnd = -mesh%nrm(1, i)*sn + mesh%nrm(2, i)*c
         e = ph*zf
         dn_e = ph*(-CI*s%k*zf*nd + mesh%nrm(3, i)*dzf)           ! n . grad E
         q = -CI*s%k*(-xx*sn + yy*c)                               ! d(phase)/dtheta
         de = e*q
         ddn_e = dn_e*q + ph*zf*(-CI*s%k)*dnd
         hh = hh + (s%dphi_p_dn(i)*e - s%phi_p(i)*dn_e)*mesh%area(i)
         dhh = dhh + (s%dphi_p_dn(i)*de - s%phi_p(i)*ddn_e)*mesh%area(i)
      END DO
      pref = -(s%omega*kochin_c(s%k, s%nu, s%h)/(2D0*s%g))*SQRT(2D0/(PI*s%k))*EXP(-CI*PI/4D0)
      a = pref*hh
      da = pref*dhh
   END SUBROUTINE far_amplitude

   SUBROUTINE far_field_pair(mesh, s1, s2, opts, f)
      ! Mean load F_12 = (F_x, F_y, M_z) from two states of the same frequency, per unit amplitude pair, with M_z about
      ! the global origin. For s1 = s2 this is the mean drift of one heading. The azimuth integrals use the trapezoid
      ! rule, doubling the number of points from opts%n_azimuth until the change is below opts%far_tolerance
      ! (absolute after scaling by rho g, or relative) or 2^14 points are reached. a(theta) is evaluated once per point.
      TYPE(panel_mesh), INTENT(IN) :: mesh
      TYPE(first_order_state), INTENT(IN) :: s1, s2
      TYPE(so_options), INTENT(IN) :: opts
      COMPLEX(8), INTENT(OUT) :: f(3)
      COMPLEX(8), ALLOCATABLE :: a1(:), da1(:), a2(:), da2(:), t1(:), dt1(:), t2(:), dt2(:)
      COMPLEX(8) :: pt(3), iv(3), fold(3), ep, a1b2, da1b2, a2b1, da2b1
      REAL(8) :: kappa, sq, k, scale
      INTEGER :: n, j

      k = s1%k
      kappa = 0.5D0*s1%rho*s1%g*group_ratio(k, s1%h)
      sq = SQRT(2D0*PI/k)
      ep = EXP(CI*PI/4D0)

      ! Terms evaluated at the incident headings
      CALL far_amplitude(mesh, s2, s1%beta, a2b1, da2b1)
      CALL far_amplitude(mesh, s1, s2%beta, a1b2, da1b2)
      pt(1:2) = -sq*([COS(s1%beta), SIN(s1%beta)]*CONJG(a2b1*ep) + [COS(s2%beta), SIN(s2%beta)]*a1b2*ep)
      pt(3) = sq*(CI*CONJG(da2b1*ep) - CI*da1b2*ep)

      n = MAX(opts%n_azimuth, 4)
      ALLOCATE (a1(n), da1(n), a2(n), da2(n))
      DO j = 1, n
         CALL far_amplitude(mesh, s1, 2D0*PI*REAL(j - 1, 8)/REAL(n, 8), a1(j), da1(j))
         CALL far_amplitude(mesh, s2, 2D0*PI*REAL(j - 1, 8)/REAL(n, 8), a2(j), da2(j))
      END DO
      iv = azimuth_integrals(a1, da1, a2, da2)
      f = assemble(iv)
      scale = s1%rho*s1%g

      DO
         IF (2*n > MAX_AZIMUTH) THEN
            WRITE (*, '(a,i0,a)') 'far_field_pair: azimuth integral not converged with ', n, ' points'
            EXIT
         END IF
         ! Double the grid: old points sit at the odd positions, only the even ones are new
         ALLOCATE (t1(2*n), dt1(2*n), t2(2*n), dt2(2*n))
         DO j = 1, n
            t1(2*j - 1) = a1(j)
            dt1(2*j - 1) = da1(j)
            t2(2*j - 1) = a2(j)
            dt2(2*j - 1) = da2(j)
            CALL far_amplitude(mesh, s1, 2D0*PI*REAL(2*j - 1, 8)/REAL(2*n, 8), t1(2*j), dt1(2*j))
            CALL far_amplitude(mesh, s2, 2D0*PI*REAL(2*j - 1, 8)/REAL(2*n, 8), t2(2*j), dt2(2*j))
         END DO
         CALL MOVE_ALLOC(t1, a1)
         CALL MOVE_ALLOC(dt1, da1)
         CALL MOVE_ALLOC(t2, a2)
         CALL MOVE_ALLOC(dt2, da2)
         n = 2*n
         fold = f
         iv = azimuth_integrals(a1, da1, a2, da2)
         f = assemble(iv)
         IF (ALL(ABS(f - fold) <= opts%far_tolerance*MAX(scale, ABS(f)))) EXIT
      END DO

   CONTAINS

      PURE FUNCTION azimuth_integrals(x1, dx1, x2, dx2) RESULT(r)
         ! Trapezoid rule on [0, 2pi): int x1 x2* e_theta and int (x2* dx1 - x1 dx2*)/(2i)
         COMPLEX(8), INTENT(IN) :: x1(:), dx1(:), x2(:), dx2(:)
         COMPLEX(8) :: r(3)
         INTEGER :: m, jj
         REAL(8) :: th

         m = SIZE(x1)
         r = (0D0, 0D0)
         DO jj = 1, m
            th = 2D0*PI*REAL(jj - 1, 8)/REAL(m, 8)
            r(1) = r(1) + x1(jj)*CONJG(x2(jj))*COS(th)
            r(2) = r(2) + x1(jj)*CONJG(x2(jj))*SIN(th)
            r(3) = r(3) + (CONJG(x2(jj))*dx1(jj) - x1(jj)*CONJG(dx2(jj)))/(2D0*CI)
         END DO
         r = r*(2D0*PI/REAL(m, 8))
      END FUNCTION azimuth_integrals

      PURE FUNCTION assemble(ints) RESULT(r)
         COMPLEX(8), INTENT(IN) :: ints(3)
         COMPLEX(8) :: r(3)

         r(1:2) = kappa*(pt(1:2) - ints(1:2))
         r(3) = -(kappa/k)*(ints(3) + pt(3))
      END FUNCTION assemble

   END SUBROUTINE far_field_pair

END MODULE FarFieldDrift
